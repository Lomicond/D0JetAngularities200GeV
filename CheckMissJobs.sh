#!/usr/bin/env bash
# CheckMissJobs v5.4 -- SUMS i obnova bez session, sipkove menu; Python 3.6+ a ROOT.
# Pouziti: bash CheckMissJobs.sh [SESSION.session.xml] [--scan-only] [--help]
# Bez session: bash CheckMissJobs.sh --recover SUMS_HASH [--request-memory 12000]
# --scan-only vytvori reporty a navrh .condor, nikdy nemeni frontu.
# --jobs 12,19,30-35 omezi pouze kandidaty k novemu odeslani (ne volbu HELD).
# --resume-report ADRESAR obnovi vyber; znovu kontroluje jen vybrane joby.
# --add-jobs 27583 prida rucni duvod opakovani; nemeni klasifikaci picoDst.
# Recovery pouziva puvodni kod/ZIP a kopirovani picoDst. Nic v nich neprepisuje.
# ROOT soubory se jen ctou. Bez potvrzeni se nic neodesila.
# Stare TXNetFile chyby samy o sobe NEJSOU dukazem nezpracovaneho vstupu.
command -v python3 >/dev/null 2>&1 || {
    echo 'Chybi python3 (potreba Python 3.6+).' >&2
    exit 1
}
exec python3 - "$@" <<'CHECKMISS_PYTHON'
import argparse
import collections
import csv
import datetime
import getpass
import hashlib
import json
import os
from pathlib import Path
import re
import select
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET


def say(text):
    print(text, flush=True)


def ask(prompt, default=''):
    # Python source arrives through stdin; read interactive replies from /dev/tty.
    try:
        with open('/dev/tty', 'r+') as tty:
            tty.write(prompt)
            tty.flush()
            value = tty.readline().strip()
            return value or default
    except OSError:
        return default


def menu_available():
    if os.environ.get('TERM', '') in ('', 'dumb'):
        return False
    try:
        fd = os.open('/dev/tty', os.O_RDWR)
        try:
            return os.isatty(fd)
        finally:
            os.close(fd)
    except OSError:
        return False


def menu(title, options, multiple=False, checked=None, default=0, minimum=0):
    """Arrow/space/Enter menu on the controlling terminal, independent of stdin.
    Alternate screen preserves the scan output. Viewport scrolls with selection.
    Raw mode and cursor are always restored, including after Ctrl-C or exceptions.
    """
    import termios
    import tty
    fd = os.open('/dev/tty', os.O_RDWR)
    previous = termios.tcgetattr(fd)
    selected = min(max(default, 0), len(options)-1)
    marked = set(checked or ())
    total = len(options) + (1 if multiple else 0)
    top = 0
    notice = ''

    def emit(value):
        data = value.encode('utf-8')
        while data:
            data = data[os.write(fd, data):]

    def safe(value):
        return re.sub(r'[\x00-\x1f\x7f]', ' ', str(value))

    def read_key():
        key = os.read(fd, 1)
        if key != b'\x1b':
            return key
        sequence = key
        # SSH can deliver escape sequences in multiple packets.
        while select.select([fd], [], [], 0.15)[0]:
            sequence += os.read(fd, 1)
            if len(sequence) >= 3 and (sequence[-1:].isalpha() or sequence.endswith(b'~')):
                break
            if len(sequence) >= 16:
                break
        return sequence

    try:
        tty.setraw(fd)
        emit('\x1b[?1049h\x1b[?25l')
        while True:
            size = os.get_terminal_size(fd)
            width = max(10, (size.columns or 80)-2)
            height = max(1, (size.lines or 24)-8)
            if selected < top:
                top = selected
            if selected >= top+height:
                top = selected-height+1
            lines = [safe(title)[:width], '',
                     'Nahoru/dolu: vyber | ' + ('Mezernik: oznacit | ' if multiple else '') + 'Enter: potvrdit',
                     'q / Esc: zrusit', '']
            for i in range(top, min(total, top+height)):
                all_row = multiple and i == len(options)
                label = 'Vse' if all_row else options[i]
                checked_now = len(marked) == len(options) if all_row else i in marked
                prefix = ('[x] ' if checked_now else '[ ] ') if multiple else ''
                row = ('> ' if i == selected else '  ') + prefix + safe(label)
                row = row[:width]
                lines.append(('\x1b[7m' + row + '\x1b[0m') if i == selected else row)
            lines.append('')
            lines.append(safe(notice or 'Polozka {}/{}{}'.format(selected+1, total,
                ' | Oznaceno: {}'.format(len(marked)) if multiple else ''))[:width])
            emit('\x1b[H\x1b[2J' + '\r\n'.join(lines))
            key = read_key()
            if key in (b'\x03', b'q', b'Q', b'\x1b', b'\x04', b''):
                raise KeyboardInterrupt
            if key in (b'\x1b[A', b'\x1bOA', b'k'):
                selected = (selected-1) % total
            elif key in (b'\x1b[B', b'\x1bOB', b'j'):
                selected = (selected+1) % total
            elif key == b'\x1b[5~':
                selected = max(0, selected-height)
            elif key == b'\x1b[6~':
                selected = min(total-1, selected+height)
            elif multiple and key == b' ':
                notice = ''
                if selected == len(options):
                    marked = set() if len(marked) == len(options) else set(range(len(options)))
                elif selected in marked:
                    marked.remove(selected)
                else:
                    marked.add(selected)
            elif key in (b'\r', b'\n'):
                if multiple and len(marked) < minimum:
                    notice = 'Vyber alespon jednu moznost (mezernik).'
                else:
                    return sorted(marked) if multiple else selected
    finally:
        try:
            emit('\x1b[0m\x1b[?25h\x1b[?1049l')
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, previous)
            os.close(fd)


def confirm(title, use_menu=False):
    if use_menu:
        return menu(title, ['Ne', 'Ano'], default=0) == 1
    return ask(title + ' [a/N]: ', 'n').lower() in ('a', 'ano', 'y', 'yes')


def write_tsv(path, header, rows):
    with path.open('w', newline='', encoding='utf-8') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(header)
        writer.writerows(rows)


def fingerprint(path):
    try:
        s = path.stat()
        return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns)
    except FileNotFoundError:
        return None


def run(argv, **kw):
    return subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True, errors='replace', **kw)


def session_info(path):
    # Parse XML as data, never execute Java XMLDecoder or shell content.
    tree = ET.parse(str(path))
    indices = set()
    for element in tree.iter():
        if element.attrib.get('property') == 'jobIndex':
            for child in element:
                if child.tag in ('int', 'long') and child.text:
                    n = int(child.text.strip())
                    if n < 0:
                        raise ValueError('Zaporny jobIndex v session.')
                    indices.add(n)
    strings = [e.text or '' for e in tree.iter('string')]
    text = '\n'.join(strings)
    # Version commonly occurs in the trailing XML comment, not a string node.
    with path.open('rb') as stream:
        stream.seek(max(0, path.stat().st_size - 8192))
        tail = stream.read().decode('utf-8', 'replace')
    versions = re.findall(r'Meta Scheduler\s+(\d+\.\d+\.\d+)', tail)
    return indices, text, versions[-1] if versions else None


def find_pattern(text, suffix, override):
    if override:
        return override
    choices = set()
    for token in re.findall(r'''[^\s"'<>\\]+\.''' + suffix + r'\b', text):
        if re.search(r'\$(?:\{JOBINDEX\}|JOBINDEX\b)', token):
            # The directory is separately configurable; names must be unambiguous.
            choices.add(token.rsplit('/', 1)[-1])
    if len(choices) != 1:
        raise ValueError('Nelze jednoznacne urcit vzor .' + suffix +
                         ': ' + repr(sorted(choices)) +
                         '. Pouzij --' + ('output' if suffix == 'root' else 'log') + '-pattern.')
    return choices.pop()


def expand_pattern(pattern, index, directory):
    result = pattern.replace('${JOBINDEX}', str(index))
    result = re.sub(r'\$JOBINDEX\b', str(index), result)
    if '$' in result:
        raise ValueError('Nevyresena promenna ve vzoru: ' + result)
    path = Path(result)
    return path if path.is_absolute() else directory / path


def index_files(base, token, indices):
    """Only session-specific list/script names can supplement missing default index 0."""
    result = collections.defaultdict(list)
    pattern = re.compile(re.escape(token) + r'_(\d+)\.(?:list|lis|csh)$', re.I)
    for directory in (base / 'list', base / 'csh', base):
        if not directory.is_dir():
            continue
        for path in directory.iterdir():
            match = pattern.search(path.name)
            if match and path.is_file():
                n = int(match.group(1))
                indices.add(n)
                if path.suffix in ('.list', '.lis'):
                    result[n].append(path)
    return result


def select_indices(value, available):
    if value == 'all':
        return set(available)
    result = set()
    for item in value.split(','):
        match = re.fullmatch(r'\s*(\d+)(?:-(\d+))?\s*', item)
        if not match:
            raise ValueError('Indexy: pouzij napr. 12,19,30-35 nebo all.')
        first, last = int(match.group(1)), int(match.group(2) or match.group(1))
        if first > last or last-first > 1000000:
            raise ValueError('Neplatny rozsah indexu: ' + item)
        result.update(range(first, last+1))
    if result - set(available):
        raise ValueError('Nezname indexy: ' + str(sorted(result-set(available))[:30]))
    return result


def discover_recovery(base):
    names = {}
    for directory in (base / 'csh', base):
        if directory.is_dir():
            for path in directory.iterdir():
                match = re.fullmatch(r'(.*?([A-Fa-f0-9]{32}))_\d+\.csh', path.name)
                if match:
                    names[match.group(2)] = match.group(1)
    return names


def condor_blocks(path, token):
    """Parse the concrete SUMS format, never evaluate shell/submit macros.
    A block is independently complete; inherited fields, includes and advanced
    Queue forms are rejected rather than silently changing their meaning.
    """
    allowed = {'universe', 'notification', 'executable', 'arguments', 'output',
               'error', 'request_memory', 'log', 'initialdir', 'kill_sig',
               'periodicremove', 'periodic_remove', 'priority'}
    required = {'universe','executable','arguments','output','error','request_memory','initialdir'}
    block = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.lower() == 'queue':
            if not required <= set(block):
                raise ValueError('Neuplny samostatny Queue blok: ' + str(path))
            found = re.findall(re.escape(token) + r'_(\d+)\.csh\b', block['arguments'], re.I)
            if len(found) != 1:
                raise ValueError('Nejednoznacny index v Arguments: ' + str(path))
            if '$(' in '\n'.join(block.values()):
                raise ValueError('Submit makra nejsou v recovery podporovana: ' + str(path))
            yield int(found[0]), dict(block)
            block = {}
            continue
        match = re.fullmatch(r'(\w+)\s*=\s*(.+)', line)
        if not match or match.group(1).lower() not in allowed:
            raise ValueError('Nepodporovany prikaz v recovery .condor: ' + line[:160])
        key, value = match.group(1).lower(), match.group(2)
        if key in block:
            raise ValueError('Opakovane pole v Queue bloku: ' + key)
        block[key] = value
    if block:
        raise ValueError('Nedokonceny Queue blok: ' + str(path))


def recovery_info(base, token, args, selected=None):
    names = discover_recovery(base)
    matches = [t for t in names if t.lower() == token.lower()]
    if len(matches) != 1:
        raise ValueError('Nenalezeny jednoznacne puvodni .csh pro --recover ' + token)
    token = matches[0]
    name = names[token]
    scripts, blocks, lists, issues = {}, {}, {}, collections.defaultdict(list)
    indices, counts, nodes = set(), set(), set()
    root_patterns, log_patterns = set(), set()
    guards = {}
    duplicates = []
    def guard(path):
        guards[path] = fingerprint(path)
    for directory in (base / 'csh', base):
        if not directory.is_dir():
            continue
        for path in directory.glob(name + '_*.condor'):
            guard(path)
            for n, fields in condor_blocks(path, token):
                if selected is not None and n not in selected:
                    continue
                indices.add(n)
                if n in blocks:
                    same = fields == blocks[n][1]
                    duplicates.append((n,'CONDOR',str(blocks[n][0]),str(path),'IDENTICAL' if same else 'CONFLICT'))
                    if not same:
                        issues[n].append('DUPLICATE_CONDOR_BLOCK')
                else:
                    blocks[n] = (path, fields)
        for path in directory.glob(name + '_*.csh'):
            match = re.fullmatch(re.escape(name) + r'_(\d+)\.csh', path.name)
            if not match:
                continue
            n = int(match.group(1))
            if selected is not None and n not in selected:
                continue
            indices.add(n)
            if n in scripts:
                guard(path)
                same = hashlib.sha256(path.read_bytes()).hexdigest() == scripts[n]['sha256']
                duplicates.append((n,'CSH',str(scripts[n]['path']),str(path),'IDENTICAL' if same else 'CONFLICT'))
                if not same:
                    issues[n].append('DUPLICATE_CSH')
                continue
            guard(path)
            raw = path.read_bytes()
            text = raw.decode('utf-8',errors='replace')
            env = dict(re.findall(r'^setenv\s+(\w+)\s+"([^"\n]*)"\s*$', text, re.M))
            scripts[n] = {'path':path, 'env':env, 'dependencies':[], 'sha256':hashlib.sha256(raw).hexdigest()}
            if env.get('JOBINDEX') != str(n) or env.get('REQUESTID','').lower() != token.lower() or env.get('JOBID','').lower() != (token+'_'+str(n)).lower():
                issues[n].append('CSH_IDENTITY_MISMATCH')
            count = env.get('SUMS_nProcesses','')
            if count.isdigit() and 0 < int(count) <= 1000000:
                counts.add(int(count))
            if env.get('SUBMITTINGNODE'):
                nodes.add(env['SUBMITTINGNODE'])
            submit_dir = env.get('SUBMITTINGDIRECTORY','')
            if not submit_dir or Path(submit_dir).resolve() != base:
                issues[n].append('SUBMITTINGDIRECTORY_MISMATCH')
            list_value = env.get('FILELIST','')
            lp = Path(list_value)
            if not list_value or not lp.is_absolute() or '$' in list_value:
                issues[n].append('INVALID_FILELIST_PATH')
            else:
                lists[n] = [lp]
                scripts[n]['dependencies'].append(lp)
                guard(lp)
                if not lp.is_file() or lp.stat().st_size == 0:
                    issues[n].append('MISSING_INPUT_LIST')
                else:
                    actual,state = load_inputs([lp])
                    expected = {canonical(v) for k,v in env.items() if re.fullmatch(r'INPUTFILE\d+',k)}
                    count_inputs = env.get('INPUTFILECOUNT','')
                    if state != 'INPUT_LIST_OK' or set(actual) != expected or not count_inputs.isdigit() or len(expected) != int(count_inputs):
                        issues[n].append('INPUT_LIST_CSH_MISMATCH')
            packages = re.findall(r'^/usr/bin/unzip\s+(\S+)\s+-d\s+\$SCRATCH\s*$', text, re.M)
            if not packages:
                issues[n].append('UNVERIFIED_SANDBOX')
            for package in packages:
                pp = Path(package)
                scripts[n]['dependencies'].append(pp)
                guard(pp)
                if not pp.is_absolute() or not pp.is_file():
                    issues[n].append('MISSING_SANDBOX')
            usercode = text.split('# User command BEGIN',1)[-1].split('# <------------------------------User command END',1)[0]
            for suffix, choices, override in [('root', root_patterns, args.output_pattern), ('log', log_patterns, args.log_pattern)]:
                try:
                    choices.add(find_pattern(usercode,suffix,override))
                except ValueError:
                    issues[n].append('AMBIGUOUS_' + suffix.upper() + '_PATTERN')
            if len(scripts) % 2000 == 0:
                say('Recovery: nacteno {} puvodnich .csh.'.format(len(scripts)))
    if len(counts) > 1 or len(nodes) > 1:
        raise ValueError('Rozpor v SUMS_nProcesses nebo SUBMITTINGNODE mezi puvodnimi skripty.')
    if counts:
        expected = next(iter(counts))
        if (indices and max(indices) >= expected) or (selected and max(selected) >= expected):
            raise ValueError('Index presahuje puvodni SUMS_nProcesses.')
        indices.update(range(expected) if selected is None else selected)
    if selected is not None:
        indices.update(selected)
    if not indices:
        raise ValueError('Nenalezeny joby pro recovery.')
    for n in indices:
        if n not in scripts:
            issues[n].append('MISSING_CSH')
        if n not in blocks:
            issues[n].append('MISSING_CONDOR_BLOCK')
        if n not in scripts or n not in blocks:
            continue
        fields = blocks[n][1]
        if Path(fields['initialdir'].strip('"')).resolve() != base:
            issues[n].append('CONDOR_INITIALDIR_MISMATCH')
        try:
            words = shlex.split(fields['arguments'].strip('"'))
        except ValueError:
            words = []
        if not words or Path(words[-1]).resolve() != scripts[n]['path'].resolve():
            issues[n].append('CONDOR_CSH_MISMATCH')
        if not re.fullmatch(r'\d+', fields['request_memory']):
            issues[n].append('UNSUPPORTED_MEMORY_EXPRESSION')
    if len(root_patterns) != 1 or len(log_patterns) != 1:
        raise ValueError('Vzory vystupu/logu nejsou jednoznacne. Zadej --output-pattern a --log-pattern.')
    return {'name':name, 'token':token, 'indices':indices, 'scripts':scripts,
            'blocks':blocks, 'lists':lists, 'issues':issues, 'guards':guards,'duplicates':duplicates,
            'root_pattern':next(iter(root_patterns)), 'log_pattern':next(iter(log_patterns)),
            'node':next(iter(nodes)) if nodes else None}


def repaired_periodic_remove(expression):
    """Only repair the two exact SUMS remove-after-start clauses.
    Preserve CPU watchdog and held-expiry, and any unrelated policy.
    """
    original = (expression or '').strip()
    if original.lower().startswith('/expr(') and original.endswith(')/'):
        original = original[6:-2]
    result = re.sub(r'\(\s*(?:NumJobStarts|JobRunCount)\s*>=\s*1\s*&&\s*JobStatus\s*==\s*1\s*\)',
                    'false',original,flags=re.I)
    return result or 'false', result != original


def policy_key(expression):
    """Ignore ClassAd pretty-print whitespace/case outside quoted strings."""
    value = str(expression).strip()
    if value.lower().startswith('/expr(') and value.endswith(')/'):
        value = value[6:-2]
    pieces = re.split(r'("(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\')',value)
    return ''.join(p if i % 2 else re.sub(r'\s+','',p).lower() for i,p in enumerate(pieces))


def memory_menu(title, current=None):
    labels = ['Ponechat puvodni RequestMemory', 'Alespon 12000 MiB', 'Alespon 16000 MiB', 'Zadat vlastni minimum v MiB']
    choice = menu(title, labels)
    if choice == 0:
        return current
    if choice in (1,2):
        return (12000,16000)[choice-1]
    value = ask('Minimalni RequestMemory v MiB: ')
    if not value.isdigit() or int(value) <= 0:
        raise ValueError('Pamet musi byt kladne cele cislo v MiB.')
    return int(value)


def recovery_unchanged(info, n):
    paths = [info['scripts'][n]['path'], info['blocks'][n][0]] + info['scripts'][n]['dependencies']
    paths += [Path(row[3]) for row in info['duplicates'] if row[0] == n]
    return all(fingerprint(p) == info['guards'][p] for p in paths)


def prepare_condor(info, ids, report, memory):
    prepared = {}
    directory = report / 'condor'
    directory.mkdir(exist_ok=True)
    for n in ids:
        if info['issues'].get(n):
            continue
        fields = dict(info['blocks'][n][1])
        if 'periodicremove' in fields:
            if 'periodic_remove' in fields:
                raise ValueError('Duplicitni zapis PeriodicRemove pro job {}.'.format(n))
            fields['periodic_remove'] = fields.pop('periodicremove')
        if memory:
            fields['request_memory'] = str(max(memory, int(fields['request_memory'])))
        for key in ('periodicremove','periodic_remove'):
            if key in fields:
                fixed, changed = repaired_periodic_remove(fields[key])
                if changed:
                    fields[key] = fixed
        # Keep previous scheduler stdout/stderr and put new diagnostics in report.
        for key, suffix in [('output','out'),('error','err'),('log','events.log')]:
            fields[key] = str(directory / ('job_{}.{}'.format(n,suffix)))
        target = directory / ('job_{}.condor'.format(n))
        target.write_text(''.join('{} = {}\n'.format(k,v) for k,v in fields.items()) + 'Queue\n')
        prepared[n] = target
    return prepared


def read_report_rows(path, required):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError('Neuplna hlavicka reportu: ' + str(path))
        rows = list(reader)
    if any(None in row or any(row.get(key) is None for key in required) for row in rows):
        raise ValueError('Neuplny radek reportu: ' + str(path))
    return rows


def load_resume(args):
    """Legacy reports have no file fingerprints. Recheck selected jobs rather
    than treating their old classification or generated submit files as current.
    """
    source = args.resume_report.resolve(strict=True)
    config = json.loads((source / 'scan_config.json').read_text())
    if not isinstance(config, dict):
        raise ValueError('Neplatny scan_config.json.')
    token = config.get('recover')
    if not isinstance(token, str) or not re.fullmatch(r'[A-Fa-f0-9]{32}', token):
        raise ValueError('--resume-report podporuje reporty recovery bez session.')
    rows = read_report_rows(source / 'jobs.tsv', ['jobIndex', 'resubmit_reasons', 'root_path', 'log_path'])
    by_id = {}
    for row in rows:
        n = int(row['jobIndex'])
        if n < 0 or n in by_id:
            raise ValueError('Duplicitni/neplatny jobIndex v jobs.tsv.')
        by_id[n] = row
    if not by_id:
        raise ValueError('Prazdny jobs.tsv.')
    exported = set(config.get('exported_jobs', []))
    manual = select_indices(args.add_jobs, by_id) if args.add_jobs else set()
    if manual & exported:
        raise ValueError('Rucni vyber obsahuje joby exportovane pro picoDst doplneni. Nejdrive vyres pripadne dvojite zapocteni.')
    candidates = {n for n, row in by_id.items() if row['resubmit_reasons'] and n not in exported}
    candidates.update(manual)
    if args.jobs:
        candidates &= select_indices(args.jobs, by_id)
    if not candidates:
        raise ValueError('Ve vyberu nejsou kandidati; --jobs pouze filtruje, --add-jobs pridava rucni opakovani.')
    # A resumed manual decision must be explicit on each invocation.
    # Ordinary candidates are reclassified below against current logs/outputs.
    for key in ('output_pattern', 'log_pattern'):
        if not isinstance(config.get(key), str) or not config[key]:
            raise ValueError('V reportu chybi ' + key)
        setattr(args, key, config[key])
    args.recover = token
    args.schedd = args.schedd or config.get('schedd')
    base = (args.project_dir or source.parent.parent.parent).resolve()
    if config.get('project_dir') and Path(config['project_dir']).resolve() != base:
        raise ValueError('Projekt neodpovida puvodnimu reportu.')
    args.project_dir = base
    # Recover custom directories from the actual recorded paths, also for v5.3.
    for field, option in [('root_path', 'output_dir'), ('log_path', 'log_dir')]:
        parents = {p.resolve() for p in {Path(row[field]).parent for row in by_id.values() if row[field]}}
        if len(parents) != 1:
            raise ValueError('Nejednoznacne adresare v jobs.tsv: ' + field)
        setattr(args, option, parents.pop())
    if config.get('rules'):
        raise ValueError('Report pouzival vlastni --log-rules; bez ulozeneho obsahu pravidel je nutna nova kontrola.')
    keyfile = source / 'root_required_keys.txt'
    args.require_key = config.get('require_keys', keyfile.read_text().splitlines() if keyfile.is_file() else [])
    if not isinstance(args.require_key, list) or any(not isinstance(key, str) for key in args.require_key):
        raise ValueError('Neplatny seznam pozadovanych ROOT objektu.')
    args.require_key = [key for key in args.require_key if key]
    args.resume_since = config.get('resume_baseline_ns') or (source / 'jobs.tsv').stat().st_mtime_ns
    if not isinstance(args.resume_since, int) or args.resume_since <= 0:
        raise ValueError('Neplatny cas puvodni kontroly v reportu.')
    args.resume_report = source
    return {'rows': by_id, 'selected': candidates, 'manual': manual & candidates,
            'root_pattern':args.output_pattern, 'log_pattern':args.log_pattern}


def previous_resume_attempts(report, since):
    """Called under the production lock held until all submissions finish.
    Even a job that already left the queue must not be resubmitted from the
    same old report. An unresolved submit intent is treated conservatively.
    """
    found = {}
    for path in report.parent.glob('*/resubmit_status.tsv'):
        if path.parent == report or path.stat().st_mtime_ns < since:
            continue
        for row in read_report_rows(path, ['jobIndex', 'state']):
            if row['state'] in ('SUBMITTING', 'DISPATCH_REPORTED', 'SUBMIT_UNCONFIRMED'):
                found[int(row['jobIndex'])] = str(path)
    return found


def submit_condor(ids, info, report, base, jobs, observed, roots, args):
    # Same session lock as the SUMS path; a second process cannot race this one.
    import fcntl
    with (report.parent / 'resubmit.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('Pro tuto produkci uz bezi jina akce CheckMissJobs.')
        if getattr(args, 'resume_report', None):
            prior = previous_resume_attempts(report, args.resume_since)
            excluded = [(n, 'PREVIOUS_SUBMIT_ATTEMPT', prior[n]) for n in ids if n in prior]
            write_tsv(report / 'resume_excluded.tsv', ['jobIndex', 'reason', 'previous_status'], excluded)
            ids = [n for n in ids if n not in prior]
            (report / 'resubmit_jobs.list').write_text(''.join(str(n)+'\n' for n in ids))
            if excluded:
                say('Resume: {} jobu uz ma novy pokus o odeslani; vynechavam. Viz resume_excluded.tsv.'.format(len(excluded)))
        prepared = prepare_condor(info, ids, report, args.request_memory)
        blocked = [(n,','.join(info['issues'].get(n,[]))) for n in ids if n not in prepared]
        write_tsv(report / 'recovery_blocked.tsv', ['jobIndex','reason'], blocked)
        if blocked:
            say('Recovery: {} jobu nelze pripravit; viz recovery_blocked.tsv.'.format(len(blocked)))
        if args.scan_only:
            say('Pripraveno {} jednotlivych .condor v {}. Nic neodeslano.'.format(len(prepared),report/'condor'))
            return 0
        if not prepared:
            say('Zadny job nelze odeslat; viz recovery_blocked.tsv.')
            return 2 if blocked else 0
        prefix = ['condor_submit'] + (['-name',args.schedd] if args.schedd else [])
        say('Overuji pripraveny submit pomoci condor_submit -dry-run (bez odeslani).')
        for position, (n, path) in enumerate(sorted(prepared.items()), 1):
            if not recovery_unchanged(info,n):
                raise ValueError('Puvodni soubory jobu {} se od kontroly zmenily.'.format(n))
            result = run(prefix + ['-dry-run',str(path.with_suffix('.classad')), '-maxjobs','1',str(path)], cwd=str(base), timeout=60)
            path.with_suffix('.dryrun.log').write_text(result.stdout)
            if result.returncode:
                raise ValueError('Condor dry-run selhal pro job {}; viz {}.'.format(n,path.with_suffix('.dryrun.log')))
            if position % 100 == 0 or position == len(prepared):
                say('Condor dry-run: {}/{}.'.format(position, len(prepared)))
        prepared_stamps = {n:fingerprint(p) for n,p in prepared.items()}
        say('Pripraveno {} jobu. Minimum pameti: {}. Puvodni vystupy vybranych jobu budou nahrazeny.'.format(len(prepared),args.request_memory or 'puvodni'))
        if not args.yes and not confirm('Odeslat {} vybranych jobu pres Condor?'.format(len(prepared)),args.use_menu):
            say('Nic nebylo odeslano; pripraveny navrh zustava v reportu.')
            return 0
        attempts = []
        failed = bool(blocked)
        for n, path in sorted(prepared.items()):
            fresh, reliable, note = queue_snapshot(info['token'],base,args.schedd)
            if not reliable:
                say('Zastavuji: fronta neni overena: ' + note)
                failed = True; break
            if n in fresh:
                attempts.append((n,'SKIPPED_IN_QUEUE','','',''))
                continue
            lp = Path(jobs[n]['log_path']) if jobs[n]['log_path'] else None
            if fingerprint(path) != prepared_stamps[n] or not recovery_unchanged(info,n) or observed[n]['root'] != fingerprint(roots[n]) or (lp and observed[n]['log'] != fingerprint(lp)):
                attempts.append((n,'SKIPPED_CHANGED','','',''))
                failed = True; continue
            if lp and lp.is_file():
                archive = report/'previous_logs'; archive.mkdir(exist_ok=True)
                shutil.copy2(str(lp),str(archive/('{}_'.format(n)+lp.name)))
            output = report/('resubmit_{}.log'.format(n))
            command = prefix + ['-terse','-maxjobs','1',str(path)]
            say('Condor resubmit SUMS {} (1 job)'.format(n))
            # Record intent before invoking: an interrupted/uncertain acceptance
            # must be reconciled with queue/history before retrying this index.
            attempts.append((n,'SUBMITTING','','',str(output)))
            write_tsv(report/'resubmit_status.tsv',['jobIndex','state','exit_code','condor_id','log'],attempts)
            with output.open('w') as stream:
                process = subprocess.run(command,cwd=str(base),stdout=stream,stderr=subprocess.STDOUT)
            response = output.read_text(errors='replace')
            accepted = re.findall(r'^\s*(\d+\.\d+)\s*-\s*(\d+\.\d+)\s*$',response,re.M)
            success = process.returncode == 0 and len(accepted) == 1 and accepted[0][0] == accepted[0][1]
            state = 'DISPATCH_REPORTED' if success else 'SUBMIT_UNCONFIRMED'
            attempts[-1] = (n,state,process.returncode,accepted[0][0] if success else '',str(output))
            write_tsv(report/'resubmit_status.tsv',['jobIndex','state','exit_code','condor_id','log'],attempts)
            if not success:
                say('Nejasny vysledek odeslani; zastavuji. Over frontu/historii a log pred opakovanim.')
                failed = True; break
        write_tsv(report/'resubmit_status.tsv',['jobIndex','state','exit_code','condor_id','log'],attempts)
        say('Hotovo. Novy vystup nahrazuje puvodni blok; nemerguj obe verze.')
        return 2 if failed else 0


URL = re.compile(r'''(?:root|file)://[^\s"'<>]+?\.picoDst\.root|/[^\s"'<>]+?\.picoDst\.root''', re.I)


def canonical(url):
    # Redirector host is not file identity; preserve the full logical path.
    value = url.strip().strip('"\'')
    value = re.sub(r'^(?:root|file)://[^/]*', '', value, flags=re.I)
    return re.sub(r'/+', '/', value)


def load_inputs(paths):
    variants = []
    for path in paths:
        urls = []
        for line in path.read_text(errors='replace').splitlines():
            if line.lstrip().startswith('#'):
                continue
            urls.extend(URL.findall(line))
        variants.append({canonical(u): u for u in urls})
    if not variants:
        return {}, 'NO_INPUT_LIST'
    if any(set(v) != set(variants[0]) for v in variants[1:]):
        return {}, 'AMBIGUOUS_INPUT_LIST'
    return variants[0], 'INPUT_LIST_OK' if variants[0] else 'EMPTY_INPUT_LIST'


def load_rules(path):
    """Optional regex JSON. Each regex MUST contain a named (?P<file>...) group.
    Keys: processed (fully read, not merely opened), skipped (zero contribution),
    partial (some events may already be included). No guessed success patterns.
    Only rules checked against real reader logs should be supplied.
    """
    if not path:
        return {}
    data = json.loads(path.read_text())
    unknown = set(data) - {'processed', 'skipped', 'partial'}
    if unknown:
        raise ValueError('Nezname klice log-rules: ' + repr(unknown))
    rules = {}
    for key, values in data.items():
        if not isinstance(values, list):
            raise ValueError('Pravidla musi byt seznam regexu: ' + key)
        rules[key] = []
        for value in values:
            regex = re.compile(value)
            if 'file' not in regex.groupindex:
                raise ValueError('Regex vyzaduje pojmenovanou skupinu (?P<file>...): ' + value)
            rules[key].append(regex)
    return rules


def scan_log(path, inputs, rules, job_token=None, job_index=None):
    """Recognize the SL22c input-initialization block, not per-file event completion.
    Read-in + Total is an admission ledger. Error-only files excluded from that
    completed ledger are SKIPPED; admitted files are NEVER called PROCESSED.
    Basenames are mapped only for this job's SUMS INPUTFILES staging directory
    and only if that basename identifies exactly one original input.
    """
    meta = {'issues': [], 'complete': False, 'init_complete': False,
            'mapped_all': True, 'identity_ok': True}
    if not path.is_file():
        return [], 'LOG_MISSING', meta
    before = fingerprint(path)
    records = {}
    basenames = collections.defaultdict(list)
    for key in inputs:
        basenames[Path(key).name].append(key)
    staged_job = re.compile(r'/([A-Fa-f0-9]{32})_(\d+)/INPUTFILES/')
    starts, totals, accepted = [], [], []
    last_file = None
    saw_close, saw_goodbye = False, False
    expect_codes = False

    def resolve(url):
        key = canonical(url)
        if key in inputs:
            return key
        match = staged_job.search(key)
        local = match and (not job_token or match.group(1).lower() == job_token.lower())
        local = local and (job_index is None or int(match.group(2)) == job_index)
        candidates = basenames[Path(key).name]
        if local and len(candidates) == 1:
            return candidates[0]
        return key

    def event(url, kind, lineno):
        key = resolve(url)
        record = records.setdefault(key, {'url': inputs.get(key, url), 'errors': 0,
            'terminal': [], 'line': 0, 'last_error': 0, 'last_terminal': 0,
            'admitted': [], 'failures': [], 'damaged': False, 'zombie': False})
        record['line'] = lineno
        if kind in ('access_error', 'damaged', 'zombie'):
            record['last_error'] = lineno
            record['failures'].append(lineno)
            if kind == 'access_error': record['errors'] += 1
            if kind == 'damaged': record['damaged'] = True
            if kind == 'zombie': record['zombie'] = True
        elif kind == 'admitted':
            record['admitted'].append(lineno)
            accepted.append((key, lineno))
        else:
            record['terminal'].append(kind)
            record['last_terminal'] = lineno
        return key

    with path.open(errors='replace') as stream:
        for lineno, line in enumerate(stream, 1):
            if re.match(r'^Processing\s+runD0JetAna\.C\(', line):
                starts.append(lineno)
                ident = staged_job.search(line)
                if ident and ((job_token and ident.group(1).lower() != job_token.lower()) or
                              (job_index is not None and int(ident.group(2)) != job_index)):
                    meta['identity_ok'] = False
            if 'StPicoD0JetAnaMaker - writing data and closing output file' in line:
                saw_close = True
            if 'This is the end of ROOT -- Goodbye' in line:
                saw_goodbye = True
            match = re.search(r'StInfo:\s+Read in picoDst file\s+(\S+\.picoDst\.root)', line)
            if match:
                last_file = (match.group(1), lineno)
                event(match.group(1), 'admitted', lineno)
            match = re.search(r'StInfo:\s+Total\s+(\d+) files have been read in\.', line)
            if match:
                totals.append((int(match.group(1)), lineno))
            urls = URL.findall(line)
            if ('Error in <TXNetFile::CreateXClient>' in line or
                    'Error in <TXNetSystem::Connect>' in line):
                for url in urls:
                    event(url, 'access_error', lineno)
            # ROOT may reject an incomplete local copy before it can read the
            # header. These diagnostics were observed in SL22c logs 20792/21101.
            # Record a failed access, NOT an unconditional terminal skip: the
            # completed initialization ledger below must exclude this input.
            # A later successful admission still counts as ADMITTED_AFTER_RETRY.
            header_error = ('Error in <TFile::Init>' in line and
                            ('not a ROOT file' in line or
                             'failed to read the file type data' in line))
            buffer_error = ('Error in <TFile::ReadBuffer>' in line and
                            'error reading all requested bytes from file' in line)
            if header_error or buffer_error:
                for url in urls:
                    event(url, 'access_error', lineno)
            damaged = 'is truncated at' in line and '<TFile::Init>' in line
            if damaged:
                for url in urls:
                    event(url, 'damaged', lineno)
                    last_file = (url, lineno)
            if 'no keys recovered, file has been made a Zombie' in line:
                # This diagnostic omits the filename; bind only to the immediately
                # preceding truncated-file line, never an arbitrarily old file.
                if last_file and last_file[1] == lineno - 1:
                    event(last_file[0], 'zombie', lineno)
            marker = re.match(r'^CHECKMISS\s+(FILE_DONE|FILE_SKIPPED|FILE_PARTIAL)\s+(\S+)\s*$', line)
            if marker:
                kind = {'FILE_DONE': 'processed', 'FILE_SKIPPED': 'skipped',
                        'FILE_PARTIAL': 'partial'}[marker.group(1)]
                event(marker.group(2), kind, lineno)
            for kind, patterns in rules.items():
                for pattern in patterns:
                    match = pattern.search(line)
                    if match:
                        event(match.group('file'), kind, lineno)
            if expect_codes:
                code_row = re.match(r'StInfo:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$', line)
                if code_row and (int(code_row.group(4)) or int(code_row.group(5))):
                    meta['issues'].append((lineno, 'CHAIN_ERROR', line.strip()))
                expect_codes = False
            if 'nStOK' in line and 'nStFatal' in line:
                expect_codes = True
            if re.search(r'(?:\*\*\* Break \*\*\*.*segmentation|Segmentation (?:fault|violation)|'
                         r'terminate called after throwing|std::bad_alloc|^Killed\s*$)', line, re.I):
                meta['issues'].append((lineno, 'CRASH', line.strip()))
            elif 'Error in <' in line or 'fatal error:' in line.lower() or 'no keys recovered' in line:
                category = 'INPUT_DIAGNOSTIC' if (urls or 'TXNet' in line or 'no keys recovered' in line) else 'REVIEW'
                meta['issues'].append((lineno, category, line.strip()))
            elif re.search(r'Input/output error|No space left on device|Disk quota exceeded|'
                           r'Permission denied|cannot open shared object file|undefined symbol', line, re.I):
                meta['issues'].append((lineno, 'REVIEW', line.strip()))

    changed = before != fingerprint(path)
    one_run = len(starts) == 1
    total_line = totals[0][1] if len(totals) == 1 else 0
    admitted_before = [k for k, n in accepted if n < total_line]
    meta['init_complete'] = (one_run and len(totals) == 1 and starts[0] < total_line and
        totals[0][0] == len(admitted_before) == len(set(admitted_before)) and
        all(starts[0] < n < total_line for _, n in accepted) and meta['identity_ok'])
    meta['complete'] = one_run and saw_close and saw_goodbye and meta['identity_ok']
    if len(starts) > 1:
        meta['issues'].append((0, 'REVIEW', 'Multiple analysis attempts in one log'))
    if starts and not meta['complete']:
        meta['issues'].append((0, 'INCOMPLETE', 'Analysis start without both output-close and ROOT goodbye'))
    rows = []
    for key, record in records.items():
        terminals = record['terminal']
        admitted = record['admitted']
        failure_in_init = (one_run and bool(record['failures']) and
                           all(starts[0] < n < total_line for n in record['failures']))
        skipped_at_init = meta['init_complete'] and failure_in_init and not admitted
        if changed:
            state = 'LOG_CHANGING'
        elif not meta['identity_ok']:
            state = 'WRONG_JOB_LOG'
        elif len(starts) > 1 or len(terminals) > 1 or (terminals and record['last_error'] > record['last_terminal']):
            state = 'AMBIGUOUS_ATTEMPTS'
        elif key not in inputs:
            # Whole-block retry is possible even when the staging-to-source map
            # is absent. Export must never contain temporary INPUTFILES paths.
            state = 'SKIPPED_UNMAPPED' if skipped_at_init else 'UNMAPPED'
            meta['mapped_all'] = False
        elif terminals:
            state = {'processed': 'PROCESSED', 'skipped': 'SKIPPED', 'partial': 'PARTIAL'}[terminals[0]]
        elif skipped_at_init:
            state = 'SKIPPED'
        elif record['damaged']:
            state = 'DAMAGED_INPUT'
        elif admitted:
            if record['last_error'] > max(admitted):
                state = 'ACCESS_UNCERTAIN'
            else:
                state = 'ADMITTED_AFTER_RETRY' if record['errors'] else 'ADMITTED'
        else:
            state = 'ACCESS_UNCERTAIN'
        rows.append((record['url'], state, record['errors'], record['line']))
    if inputs and meta['init_complete']:
        known = {k for k,r in records.items() if r['admitted'] or (r['failures'] and all(n < total_line for n in r['failures']))}
        for key in sorted(set(inputs) - known):
            rows.append((inputs[key], 'NOT_ACCOUNTED', 0, total_line))
            meta['mapped_all'] = False
    if changed:
        state = 'LOG_CHANGING'
    elif not meta['identity_ok']:
        state = 'WRONG_JOB_LOG'
    elif len(starts) > 1:
        state = 'MULTIPLE_ATTEMPTS'
    elif any(item[1] in ('CRASH','CHAIN_ERROR') for item in meta['issues']):
        state = 'LOG_FAILURE'
    elif starts and not meta['complete']:
        state = 'LOG_INCOMPLETE'
    else:
        state = 'LOG_FINISHED' if meta['complete'] else 'LOG_READ'
    return rows, state, meta

STATUS = {1: 'IDLE', 2: 'RUNNING', 3: 'REMOVING', 4: 'COMPLETED_IN_QUEUE',
          5: 'HELD', 6: 'TRANSFERRING', 7: 'SUSPENDED'}
UNCERTAIN_PICO = ('ACCESS_UNCERTAIN', 'UNMAPPED', 'AMBIGUOUS_ATTEMPTS',
                  'LOG_CHANGING', 'WRONG_JOB_LOG', 'NOT_ACCOUNTED', 'SKIPPED_UNMAPPED')


def queue_snapshot(token, base, schedd):
    command = ['condor_q']
    if schedd:
        command += ['-name', schedd]
    command += [getpass.getuser(), '-json', '-attributes',
                'ClusterId,ProcId,JobStatus,Cmd,Args,Arguments,Iwd,JobBatchName,HoldReason,'
                'HoldReasonCode,HoldReasonSubCode,RequestMemory,MemoryUsage,NumJobStarts,PeriodicRemove']
    try:
        process = subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                 universal_newlines=True,errors='replace',timeout=45)
        if process.returncode:
            return {}, False, 'condor_q exit {}: {}'.format(process.returncode,(process.stdout+'\n'+process.stderr)[-2000:])
        raw = process.stdout.strip()
        diagnostic = process.stderr.strip()
        # Observed on starsub04: a successful query with zero jobs emits no
        # JSON at all. Only accept that form with rc=0 and clean stderr.
        if not raw and diagnostic:
            return {}, False, 'condor_q: prazdny stdout s diagnostikou: ' + diagnostic[-2000:]
        try:
            ads = json.loads(raw) if raw else []
        except ValueError as exc:
            return {}, False, 'Neplatny JSON condor_q: {}; stdout={!r}; stderr={!r}'.format(str(exc),raw[:500],diagnostic[:500])
        if not isinstance(ads, list) or any(not isinstance(ad,dict) for ad in ads):
            raise ValueError('Vystup condor_q neni seznam job ClassAds.')
    except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
        return {}, False, str(exc)
    result = collections.defaultdict(list)
    unassigned = []
    pattern = re.compile(re.escape(token) + r'_(\d+)(?=\D|$)', re.I)
    for ad in ads:
        text = ' '.join(str(ad.get(field, '')) for field in
                        ('Cmd', 'Args', 'Arguments', 'JobBatchName'))
        found = set(int(n) for n in pattern.findall(text))
        if not found:
            # A job from this project/session that cannot be mapped blocks automation.
            iwd = ad.get('Iwd')
            same_dir = bool(iwd) and os.path.realpath(iwd) == str(base)
            if token.lower() in text.lower() or same_dir:
                unassigned.append(str(ad.get('ClusterId')) + '.' + str(ad.get('ProcId')))
        for n in found:
            result[n].append({'id': '{}.{}'.format(ad.get('ClusterId'), ad.get('ProcId')),
                              'state': STATUS.get(ad.get('JobStatus'), 'UNKNOWN'),
                              'reason': ad.get('HoldReason', ''),
                              'code': ad.get('HoldReasonCode', ''),
                              'subcode': ad.get('HoldReasonSubCode', ''),
                              'request_memory': ad.get('RequestMemory', ''),
                              'memory_usage': ad.get('MemoryUsage', ''),
                              'starts': ad.get('NumJobStarts', ''),
                              'periodic_remove': ad.get('PeriodicRemove', False)})
    note = ('Nelze priradit joby tohoto projektu: ' + ', '.join(unassigned)) if unassigned else ''
    if diagnostic:
        note += ('; ' if note else '') + 'condor_q stderr: ' + diagnostic[-2000:]
    elif not raw:
        note = 'Prazdna fronta: condor_q vratil exit 0, prazdny stdout a prazdny stderr.'
    return result, not unassigned, note


def held_advice(ad):
    reason = ad['reason'].lower()
    if re.search(r'memory|out of mem|\boom\b', reason):
        return 'Over pamet a limit; samotny release nemeni RequestMemory.'
    if 'no space' in reason or 'quota' in reason:
        return 'Nejdrive uvolnit misto/kvotu v uvedenem umisteni.'
    if 'permission' in reason or 'denied' in reason:
        return 'Nejdrive opravit pristup k uvedene ceste.'
    if 'no such file' in reason or 'not found' in reason:
        return 'Over existenci vstupu/executable a dostupnost mountu.'
    if re.search(r'timeout|timed out|connection|transfer', reason):
        return 'Over dostupnost serveru a smer prenosu; po naprave lze release.'
    return 'Procti HoldReason; po odstraneni priciny lze release.'


def manage_held(queue, reliable, token, base, report, jobs, args):
    held = [(n, ad) for n in sorted(jobs) for ad in queue.get(n, []) if ad['state'] == 'HELD']
    write_tsv(report / 'held.tsv', ['jobIndex','condor_id','HoldReasonCode','HoldReasonSubCode',
        'HoldReason','RequestMemory_MiB','MemoryUsage_MiB','NumJobStarts','PeriodicRemove','advice'],
        [(n,a['id'],a['code'],a['subcode'],a['reason'],a['request_memory'],a['memory_usage'],
          a['starts'],a['periodic_remove'],held_advice(a)) for n,a in held])
    if not held:
        if args.release_held:
            say('Nenalezen zadny HELD job teto session k uvolneni.')
        return False
    say('\nHELD: {} jobu; podrobnosti v held.tsv.'.format(len(held)))
    for n, ad in held[:15]:
        say('  SUMS {} -> Condor {} [{} / {}]: {}'.format(n, ad['id'],ad['code'],ad['subcode'],ad['reason']))
        say('    ' + held_advice(ad))
    if args.scan_only:
        return False
    selection = args.release_held
    if not selection and not args.yes:
        if args.use_menu:
            choices = sorted({n for n,_ in held})
            labels = []
            for n in choices:
                details = [ad for index,ad in held if index == n]
                labels.append('SUMS {} | Condor {} | {}'.format(n,
                    ','.join(ad['id'] for ad in details), details[0]['reason']))
            marked = menu('HELD: oznac joby k uvolneni po naprave priciny; prazdny vyber = nic',
                          labels, multiple=True)
            selection = ','.join(str(choices[i]) for i in marked)
        else:
            selection = ask('Uvolnit vybrane HELD po naprave priciny? SUMS indexy oddelene carkou, Enter = nic: ')
    if not selection:
        return False
    if not reliable:
        say('Release zablokovan: fronta neni spolehlive zmapovana.')
        return True
    available = {n for n,_ in held}
    if selection == 'all':
        selected = available
    else:
        try:
            selected = {int(x.strip()) for x in selection.split(',')}
        except ValueError:
            raise ValueError('--release-held: zadej SUMS indexy oddelene carkou nebo all.')
    if selected - available:
        raise ValueError('Tyto indexy nebyly HELD v teto session: ' + str(sorted(selected-available)))
    memory = args.request_memory
    if memory is None and args.use_menu:
        memory = memory_menu('Pamet pro vybrane HELD joby pred release')
    say('Release vrati existujici joby do fronty. Disk a prava musis opravit predem.')
    plans = {}
    for n, ad in held:
        if n not in selected:
            continue
        expression = str(ad['periodic_remove'])
        fixed, changed = repaired_periodic_remove(expression)
        # Unknown idle-removal policies must not be guessed at. Known SUMS
        # clauses are removed; unrelated CPU watchdog/held-expiry stay intact.
        if re.search(r'jobstatus(?:==|=\?=)1\b|1==jobstatus',re.sub(r'\s+','',fixed).lower()):
            raise ValueError('HELD {} ma nerozpoznane pravidlo odstraneni v idle; nejdrive uprav PeriodicRemove rucne.'.format(ad['id']))
        edits = []
        if changed:
            edits.append(('PeriodicRemove',fixed))
            say('  {}: odstranim z PeriodicRemove dve zname SUMS podminky idle po predchozim startu.'.format(ad['id']))
        if memory:
            try:
                desired = max(memory,int(ad['request_memory']))
            except (ValueError,TypeError):
                raise ValueError('HELD {} nema ciselny RequestMemory; nelze bezpecne menit.'.format(ad['id']))
            if desired != ad['request_memory']:
                edits.append(('RequestMemory',str(desired)))
            say('  {}: RequestMemory {} -> {} MiB'.format(ad['id'],ad['request_memory'],desired))
        plans[ad['id']] = edits
    (report/'held_edit_plan.json').write_text(json.dumps(plans,indent=2))
    if not args.yes and not confirm('Uvolnit {}?'.format(sorted(selected)), args.use_menu):
        return False
    import fcntl
    action_lock = (report.parent / 'resubmit.lock').open('a')
    try:
        fcntl.flock(action_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise ValueError('Pro tuto session uz bezi jina akce CheckMissJobs.')
    results = []
    failed = False
    for n, old in held:
        if n not in selected:
            continue
        fresh, ok, note = queue_snapshot(token, base, args.schedd)
        current = fresh.get(n, [])
        if not ok:
            say('Release zastaven: ' + note)
            failed = True
            break
        if len(queue.get(n, [])) != 1 or len(current) != 1:
            results.append((n,old['id'],'SKIPPED_MULTIPLE_OR_GONE',''))
            continue
        ad = current[0]
        if any(ad[k] != old[k] for k in ('id','state','reason','code','subcode','request_memory','periodic_remove')):
            results.append((n,old['id'],'SKIPPED_CHANGED',''))
            continue
        if not re.fullmatch(r'\d+\.\d+', ad['id']):
            raise ValueError('Neplatne Condor ID: ' + ad['id'])
        edits = plans[ad['id']]
        for attribute,value in edits:
            command = ['condor_qedit'] + (['-name',args.schedd] if args.schedd else []) + [ad['id'],attribute,value]
            try:
                changed = run(command,cwd=str(base),timeout=45)
                code,output = changed.returncode,changed.stdout
            except (OSError,subprocess.TimeoutExpired) as exc:
                code,output = -1,str(exc)
            editlog = report/('qedit_{}_{}.log'.format(ad['id'],attribute))
            editlog.write_text(output)
            if code:
                results.append((n,ad['id'],'QEDIT_UNCONFIRMED',str(editlog)))
                failed = True; break
        if failed:
            break
        # A zero qedit exit code alone is not sufficient: verify values and
        # held identity before releasing the same existing Cluster.Proc.
        if edits:
            updated,ok,note = queue_snapshot(token,base,args.schedd)
            rows = updated.get(n,[])
            good = ok and len(rows)==1 and all(rows[0][k]==ad[k] for k in ('id','state','reason','code','subcode'))
            for attr,value in edits:
                if good and attr == 'RequestMemory':
                    good = str(rows[0]['request_memory']) == value
                if good and attr == 'PeriodicRemove':
                    actual,_ = repaired_periodic_remove(str(rows[0]['periodic_remove']))
                    _,still_needs_fix = repaired_periodic_remove(str(rows[0]['periodic_remove']))
                    good = not still_needs_fix and policy_key(actual) == policy_key(value)
            if not good:
                results.append((n,ad['id'],'QEDIT_VERIFY_FAILED',''))
                failed=True; break
        logpath = Path(jobs[n]['log_path']) if jobs[n]['log_path'] else None
        if logpath and logpath.is_file():
            oldlogs = report / 'previous_logs'
            oldlogs.mkdir(exist_ok=True)
            shutil.copy2(str(logpath), str(oldlogs / ('held_{}_'.format(n) + logpath.name)))
        command = ['condor_release'] + (['-name', args.schedd] if args.schedd else []) + [ad['id']]
        say('Release SUMS {} / Condor {}'.format(n,ad['id']))
        try:
            proc = run(command, cwd=str(base), timeout=45)
            output, code = proc.stdout, proc.returncode
        except (OSError, subprocess.TimeoutExpired) as exc:
            output, code = str(exc), -1
        log = report / ('release_{}.log'.format(ad['id']))
        log.write_text(output)
        results.append((n,ad['id'],'RELEASE_REPORTED' if code == 0 else 'RELEASE_UNCONFIRMED',str(log)))
        write_tsv(report / 'release_status.tsv', ['jobIndex','condor_id','state','log'],results)
        if code:
            say('Release nepotvrzen; over frontu a log. Dalsi uvolnovani zastavuji.')
            failed = True
            break
    write_tsv(report / 'release_status.tsv', ['jobIndex','condor_id','state','log'],results)
    action_lock.close()
    return failed


ROOT_MACRO = r'''
#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>
#include <TSystem.h>
#include <fstream>
#include <string>
#include <vector>
void CheckMissRoot() {
  std::ifstream input(gSystem->Getenv("CHECKMISS_ROOT_INPUT"));
  std::ifstream required(gSystem->Getenv("CHECKMISS_ROOT_KEYS"));
  std::ofstream out(gSystem->Getenv("CHECKMISS_ROOT_RESULT"));
  std::vector<std::string> keys; std::string row;
  while (std::getline(required,row)) if (!row.empty()) keys.push_back(row);
  while (std::getline(input,row)) {
    std::string::size_type tab = row.find('\t');
    if (tab == std::string::npos) continue;
    std::string id = row.substr(0,tab), path = row.substr(tab+1);
    TFile *file = TFile::Open(path.c_str(), "READ");
    std::string state = "READABLE", detail = "structural check only";
    if (!file || file->IsZombie()) {state="BROKEN"; detail="open failed or zombie";}
    else if (file->TestBit(TFile::kRecovered)) {state="RECOVERED"; detail="file required recovery";}
    else if (!file->GetListOfKeys() || file->GetListOfKeys()->GetSize()==0) {
      state="BROKEN"; detail="no saved keys";
    } else {
      for (unsigned int i=0; i<keys.size(); ++i) {
        if (!file->Get(keys[i].c_str())) {state="BROKEN"; detail="missing/unreadable key: "+keys[i]; break;}
      }
    }
    if (file) {file->Close(); delete file;}
    out << id << '\t' << state << '\t' << detail << std::endl;
  }
}
'''


def root_check(paths, executable, keys, report, timeout):
    results = {}
    if not paths:
        return results
    if not shutil.which(executable):
        return {n: ('UNCHECKED', 'ROOT executable not found: ' + executable) for n in paths}
    (report / 'CheckMissRoot.C').write_text(ROOT_MACRO)
    (report / 'root_required_keys.txt').write_text('\n'.join(keys) + '\n')
    items = sorted(paths.items())
    # Batches amortize ROOT startup. A crash/timeout never marks unvisited files OK.
    for start in range(0, len(items), 100):
        chunk = items[start:start + 100]
        inp = report / 'root_inputs.tsv'
        out = report / ('root_result_{:06d}.tsv'.format(start))
        inp.write_text(''.join('{}\t{}\n'.format(n, p) for n, p in chunk))
        env = os.environ.copy()
        env.update(CHECKMISS_ROOT_INPUT=str(inp), CHECKMISS_ROOT_RESULT=str(out),
                   CHECKMISS_ROOT_KEYS=str(report / 'root_required_keys.txt'))
        before = {n: fingerprint(p) for n, p in chunk}
        log = report / ('root_check_{:06d}.log'.format(start))
        try:
            with log.open('w') as stream:
                proc = subprocess.run([executable, '-b', '-l', '-q', str(report / 'CheckMissRoot.C')],
                                      env=env, stdout=stream, stderr=subprocess.STDOUT, timeout=timeout)
            message = 'ROOT checker ended with status ' + str(proc.returncode)
        except (OSError, subprocess.TimeoutExpired) as exc:
            message = str(exc)
        parsed = {}
        if out.exists():
            for line in out.read_text(errors='replace').splitlines():
                fields = line.split('\t', 2)
                if len(fields) == 3 and fields[0].isdigit() and fields[1] in ('READABLE','BROKEN','RECOVERED'):
                    parsed[int(fields[0])] = (fields[1], fields[2])
        for n, path in chunk:
            results[n] = parsed.get(n, ('UNCHECKED', message + '; see ' + log.name))
            if before[n] != fingerprint(path):
                results[n] = ('CHANGING', 'file changed during ROOT check')
        say('ROOT: {}/{} souboru zkontrolovano.'.format(min(start + 100, len(items)), len(items)))
    return results


def main():
    parser = argparse.ArgumentParser(prog='CheckMissJobs.sh', description='Kontrola SUMS vystupu a opatrny resubmit celych jobu.',
        epilog='Prvni beh: bash CheckMissJobs.sh --recover HASH --scan-only (bez session), '
               'nebo bash CheckMissJobs.sh SESSION.session.xml --scan-only. '
               'Spoustej na puvodnim submit uzlu, nebo zadej --schedd. '
               'Nejasna TXNetFile chyba sama nikdy nespousti resubmit. '
               'Log pravidla je nutne overit na skutecnych logach readeru.')
    parser.add_argument('session', nargs='?', type=Path)
    parser.add_argument('--recover', metavar='SUMS_HASH', help='bez session; puvodni .csh/.condor/.list v projektu')
    parser.add_argument('--resume-report', type=Path, metavar='ADRESAR', help='pokracovat z recovery reportu; znovu proverit pouze vybrane joby')
    parser.add_argument('--add-jobs', help='rucne pridat indexy k opakovani (napr. 27583); klasifikace picoDst se nemeni')
    parser.add_argument('--request-memory', type=int, metavar='MiB', help='minimalni pamet novych Condor jobu a vybranych HELD; napr. 12000')
    parser.add_argument('--jobs', help='omezit opakovani na indexy, napr. 12,19,30-35; pouze zjisteni kandidati')
    parser.add_argument('--scan-only', action='store_true', help='jen kontrola a reporty, bez odesilani')
    parser.add_argument('--project-dir', type=Path, help='default: adresar session, pri --recover aktualni adresar')
    parser.add_argument('--output-dir', type=Path, help='default: PROJECT/production')
    parser.add_argument('--log-dir', type=Path, help='default: PROJECT/log')
    parser.add_argument('--output-pattern', help='napr. D0_Output_2014_lowMid_${JOBINDEX}.root (v jednoduchych uvozovkach)')
    parser.add_argument('--log-pattern', help='napr. Log_2014_lowMid_${JOBINDEX}.log')
    parser.add_argument('--checks', help='missing,root,pico nebo podmnozina; bez argumentu vyber v menu')
    parser.add_argument('--no-menu', action='store_true', help='vypnout sipkove menu; vhodne pro logovani/skripty')
    parser.add_argument('--root-command', default='root' if shutil.which('root') else 'root4star')
    parser.add_argument('--root-timeout', type=int, default=180, help='sekundy na davku max. 100 ROOT souboru')
    parser.add_argument('--require-key', action='append', default=[], help='pozadovany ROOT objekt; lze opakovat')
    parser.add_argument('--log-rules', type=Path, help='JSON regexy processed/skipped/partial s (?P<file>...)')
    parser.add_argument('--export-picodst', action='store_true', help='export pouze potvrzenych SKIPPED; neresubmituje tyto pico-only joby')
    parser.add_argument('--schedd', help='jmeno puvodniho Condor schedd, pokud nejsi na puvodnim uzlu')
    parser.add_argument('--submit-command', default='star-submit-beta')
    parser.add_argument('--release-held', help='po naprave priciny uvolnit uvedene SUMS indexy (napr. 3321,3379), nebo all pro tuto session')
    parser.add_argument('--yes', action='store_true', help='odeslat jednoznacne kandidaty bez interaktivni otazky')
    args = parser.parse_args()
    if args.resume_report and (args.session or args.recover):
        parser.error('--resume-report se nepouziva soucasne se session nebo --recover.')
    if args.resume_report and (args.export_picodst or args.release_held or args.checks or args.log_rules or args.require_key or
                               args.output_dir or args.log_dir or args.output_pattern or args.log_pattern):
        parser.error('--resume-report obnovuje cele joby s nastavenim reportu; nelze menit kontroly/cesty, exportovat picoDst ani uvolnovat HELD.')
    if args.session and args.recover:
        parser.error('Zadej bud session, nebo --recover SUMS_HASH.')
    if args.recover and not re.fullmatch(r'[A-Fa-f0-9]{32}',args.recover):
        parser.error('--recover vyzaduje 32znakovy SUMS hash.')
    if args.request_memory is not None and args.request_memory <= 0:
        parser.error('--request-memory musi byt kladne cislo v MiB.')
    args.use_menu = not args.no_menu and not args.yes and menu_available()
    resume = load_resume(args) if args.resume_report else None
    checks = set(args.checks.split(',')) if args.checks is not None else {'missing','root','pico'}
    if not checks or checks - {'missing', 'root', 'pico'}:
        parser.error('--checks: povoleno missing,root,pico')
    if args.root_timeout <= 0:
        parser.error('--root-timeout musi byt kladny')
    if args.export_picodst and 'pico' not in checks:
        parser.error('--export-picodst vyzaduje kontrolu pico')
    if not args.session and not args.recover:
        search_base = (args.project_dir or Path.cwd()).resolve()
        sessions = sorted(search_base.glob('D0*.session.xml'), key=lambda p: p.stat().st_mtime, reverse=True)
        labels = []
        for i, path in enumerate(sessions, 1):
            stamp = datetime.datetime.fromtimestamp(path.stat().st_mtime).strftime('%Y-%m-%d %H:%M')
            labels.append('[{}] {}'.format(stamp, path.name))
        recovery_choices = sorted(discover_recovery(search_base).items())
        choices = [('session',p) for p in sessions] + [('recover',t) for t,_ in recovery_choices]
        labels += ['BEZ SESSION | '+name for _,name in recovery_choices]
        if not choices:
            raise ValueError('Nenalezena session ani puvodni SUMS .csh; zadej --project-dir.')
        if args.use_menu:
            selected = menu('Vyber produkci / zpusob obnovy',labels)
        else:
            say('Vyber session (cas posledni zmeny, nikoli nutne vytvoreni):')
            for i,label in enumerate(labels, 1):
                say('  {}. {}'.format(i,label))
            answer = ask('Cislo produkce: ', '1' if len(choices) == 1 else '')
            if not answer.isdigit() or not 1 <= int(answer) <= len(choices):
                raise ValueError('Zadej session nebo --recover HASH, nebo vyber platne cislo.')
            selected = int(answer)-1
        mode,value = choices[selected]
        if mode == 'session':
            args.session = value
        else:
            args.recover = value
    if args.use_menu and args.checks is None and not resume:
        check_names = ['missing','root','pico']
        check_labels = ['Missing ROOT files - chybejici vystupy',
                        'Broken ROOT files - otevreni, Zombie, recovery',
                        'Missing picoDst - kontrola vstupu a logu']
        marked = menu('Vyber kontroly (mezernik = zaskrtnout)', check_labels,
                      multiple=True, minimum=1,
                      checked=[2] if args.export_picodst else [])
        checks = {check_names[i] for i in marked}
        if args.export_picodst and 'pico' not in checks:
            parser.error('--export-picodst vyzaduje vybranou kontrolu Missing picoDst')
    session = args.session.resolve(strict=True) if args.session else None
    session_stamp = fingerprint(session) if session else None
    base = (args.project_dir or (session.parent if session else Path.cwd())).resolve()
    output_dir = (args.output_dir or base / 'production').resolve()
    log_dir = (args.log_dir or base / 'log').resolve()
    recovery = None
    if args.recover:
        if resume:
            say('Resume: {} vybranych jobu z reportu {}. Proverim jejich soucasne vstupy, logy a ROOT vystupy.'.format(
                len(resume['selected']), args.resume_report))
            # Infer from today's scripts; an old report must not override a
            # changed output target and thereby validate/backup the wrong file.
            args.output_pattern = args.log_pattern = None
        recovery = recovery_info(base,args.recover,args, resume['selected'] if resume else None)
        name,token,indices,lists = (recovery[k] for k in ('name','token','indices','lists'))
        version = None
        output_pattern,log_pattern = recovery['root_pattern'],recovery['log_pattern']
        if resume and (output_pattern != resume['root_pattern'] or log_pattern != resume['log_pattern']):
            raise ValueError('Vzory vystupu/logu v puvodnich skriptech se lisi od reportu; je nutna nova kontrola.')
        if not args.schedd:
            args.schedd = recovery['node']
        if not args.schedd:
            raise ValueError('Nelze urcit puvodni schedd; zadej --schedd.')
    else:
        name = session.name[:-12]
        match = re.search(r'([A-Fa-f0-9]{32})$', name)
        if not match:
            raise ValueError('V nazvu session neni SUMS hash (32 hex znaku). Nelze bezpecne mapovat frontu.')
        token = match.group(1)
        indices, xml_text, version = session_info(session)
        lists = index_files(base, token, indices)
        output_pattern = find_pattern(xml_text, 'root', args.output_pattern)
        log_pattern = find_pattern(xml_text, 'log', args.log_pattern) if 'pico' in checks else None
    if not indices:
        raise ValueError('Nenalezeny skutecne jobIndex ani vstupni .list/.csh teto session.')
    manual_ids = resume['manual'] if resume else (select_indices(args.add_jobs, indices) if args.add_jobs else set())
    if resume:
        for n in indices:
            row = resume['rows'][n]
            for key, pattern, directory in [('root_path',output_pattern,output_dir), ('log_path',log_pattern,log_dir)]:
                if not row[key] or Path(row[key]).resolve() != expand_pattern(pattern,n,directory).resolve():
                    raise ValueError('Cesta v reportu neodpovida puvodnimu jobu {}: {}'.format(n,key))
    rules = load_rules(args.log_rules)
    report_base = base / 'CheckMissReports' / name
    report_base.mkdir(parents=True, exist_ok=True)
    report = Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime('%Y%m%d_%H%M%S_'), dir=str(report_base)))
    say('{}: {} | {} indexu ({} az {})'.format('Recovery' if recovery else 'Session', name, len(indices), min(indices), max(indices)))
    if recovery:
        write_tsv(report/'duplicate_definitions.tsv',['jobIndex','type','primary_path','other_path','comparison'],recovery['duplicates'])
        same_count = sum(row[-1] == 'IDENTICAL' for row in recovery['duplicates'])
        if same_count:
            say('Recovery: {} shodnych duplicitnich definic slouceno; viz duplicate_definitions.tsv.'.format(same_count))
        write_tsv(report/'recovery_inventory.tsv',['jobIndex','csh','condor_source','issues'],
                  [(n,str(recovery['scripts'].get(n,{}).get('path','')),str(recovery['blocks'][n][0]) if n in recovery['blocks'] else '',
                    ','.join(recovery['issues'].get(n,[]))) for n in sorted(indices)])
        say('Puvodni soubory se nemeni. Recovery prebere jejich kod i chovani kopirovani picoDst.')
    if resume:
        shutil.copy2(str(args.resume_report/'jobs.tsv'), str(report/'previous_jobs.tsv'))
        (report/'resume_source.json').write_text(json.dumps({'report':str(args.resume_report),
            'jobs_mtime_ns':args.resume_since, 'selected':sorted(indices), 'manual':sorted(manual_ids)}, indent=2))
    write_tsv(report/'manual_jobs.tsv', ['jobIndex','reason'], [(n,'MANUAL_RETRY') for n in sorted(manual_ids)])
    say('ROOT vzor: {} / {}'.format(output_dir, output_pattern))
    say('Report: ' + str(report))
    say('Fronta: ' + (args.schedd or 'lokalni schedd; spoustej na puvodnim submit uzlu'))
    if 'pico' in checks:
        say('PICO: kontrola Read in/Total pri inicializaci; samostatna chyba pripojeni zustava nejista.')
    if 'root' in checks:
        say('ROOT: kontrola otevreni/keys/recovery; necita vsechny udalosti a nedokazuje uplnost produkce.')
    queue, queue_ok, queue_note = queue_snapshot(token, base, args.schedd)
    (report / 'queue_initial.json').write_text(json.dumps({'jobs': queue, 'reliable': queue_ok, 'note': queue_note}, indent=2))
    if not queue_ok:
        say('FRONTA NEOVERENA: ' + queue_note + '\nReport vznikne, automaticky resubmit bude zablokovan.')
    roots = {n: expand_pattern(output_pattern, n, output_dir) for n in indices}
    root_results = {}
    to_check = {}
    observed = {}
    for n in sorted(indices):
        path = roots[n]
        observed[n] = {'root': fingerprint(path)}
        if n in queue:
            root_results[n] = ('IN_QUEUE', 'job still in queue; output not inspected')
        elif not path.exists():
            root_results[n] = ('MISSING', 'output absent')
        elif not path.is_file():
            root_results[n] = ('UNCHECKED', 'output path is not a regular file')
        elif 'root' in checks:
            to_check[n] = path
        else:
            root_results[n] = ('EXISTS_UNCHECKED', 'ROOT check not selected')
    root_results.update(root_check(to_check, args.root_command, args.require_key, report, args.root_timeout))
    records = []
    log_issues = []
    jobs = {}
    input_owners = collections.defaultdict(set)
    progress_start = progress_last = time.monotonic()
    say('Logy a vstupni seznamy: zacinam kontrolu {} jobu.'.format(len(indices)))
    for position, n in enumerate(sorted(indices), 1):
        root_state, root_detail = root_results[n]
        input_paths = lists.get(n, [])
        if any(not p.is_file() for p in input_paths):
            inputs,input_state = {},'NO_INPUT_LIST'
        else:
            inputs, input_state = load_inputs(input_paths) if 'pico' in checks else ({}, 'NOT_CHECKED')
        for key in inputs:
            input_owners[key].add(n)
        logpath = expand_pattern(log_pattern, n, log_dir) if log_pattern else None
        observed[n]['log'] = fingerprint(logpath) if logpath else None
        if n in queue:
            pico_rows, log_state = [], 'IN_QUEUE'
            log_meta = {'issues': []}
        elif 'pico' in checks:
            pico_rows, log_state, log_meta = scan_log(logpath, inputs, rules, token, n)
        else:
            pico_rows, log_state = [], 'NOT_CHECKED'
            log_meta = {'issues': []}
        for line, category, message in log_meta['issues']:
            log_issues.append((n, line, category, message, str(logpath)))
        counts = collections.Counter(row[1] for row in pico_rows)
        reasons = []
        if root_state == 'MISSING' and 'missing' in checks:
            reasons.append('MISSING_ROOT')
        if root_state in ('BROKEN', 'RECOVERED') and 'root' in checks:
            reasons.append(root_state + '_ROOT')
        if counts['SKIPPED']:
            reasons.append('SKIPPED_PICO')
        if counts['PARTIAL']:
            reasons.append('PARTIAL_PICO')
        if counts['SKIPPED_UNMAPPED']:
            reasons.append('SKIPPED_PICO_UNMAPPED')
        if counts['DAMAGED_INPUT']:
            reasons.append('DAMAGED_INPUT_PICO')
        if log_state in ('LOG_FAILURE', 'LOG_INCOMPLETE'):
            reasons.append(log_state)
        if n in manual_ids:
            reasons.append('MANUAL_RETRY')
        for url, state, errors, line in pico_rows:
            records.append((n, url, state, errors, line, str(logpath)))
        jobs[n] = {'reasons': reasons, 'root_state': root_state, 'root_detail': root_detail,
                   'inputs': input_state, 'log': log_state, 'counts': counts,
                   'log_meta': log_meta,
                   'root_path': str(roots[n]), 'log_path': str(logpath) if logpath else '',
                   'queue': queue.get(n, []), 'pico_rows': pico_rows}
        now = time.monotonic()
        if position % 200 == 0 or now-progress_last >= 10 or position == len(indices):
            say('Logy a vstupy: {}/{} jobu; {:.0f} s.'.format(position,len(indices),now-progress_start))
            progress_last = now
    candidates = sorted(n for n, j in jobs.items() if j['reasons'] and n not in queue)
    report_rows = []
    for n, j in sorted(jobs.items()):
        report_rows.append((n, ','.join(x['state'] + ':' + x['id'] for x in j['queue']),
                           j['root_state'], j['root_detail'], j['inputs'], j['log'],
                           j['counts']['SKIPPED'], j['counts']['PARTIAL'],
                           sum(j['counts'][k] for k in UNCERTAIN_PICO),
                           ','.join(j['reasons']), j['root_path'], j['log_path']))
    write_tsv(report / 'jobs.tsv', ['jobIndex','queue','root_state','root_detail','input_list','log_state',
                                   'skipped_pico','partial_pico','uncertain_pico','resubmit_reasons','root_path','log_path'], report_rows)
    write_tsv(report / 'pico_status.tsv', ['jobIndex','url','state','access_errors','last_line','log_path'], records)
    write_tsv(report / 'log_issues.tsv', ['jobIndex','line','category','message','log_path'], log_issues)
    (report / 'candidate_jobs.list').write_text(''.join(str(n)+'\n' for n in candidates))
    summary = collections.Counter(j['root_state'] for j in jobs.values())
    say('\nROOT stavy: ' + ', '.join('{}={}'.format(k,v) for k,v in sorted(summary.items())))
    pico_summary = collections.Counter(row[2] for row in records)
    if 'pico' in checks:
        say('PICO stavy: ' + (', '.join('{}={}'.format(k,v) for k,v in sorted(pico_summary.items())) or 'zadne rozpoznane zaznamy'))
        say('Chybejici logy: {}'.format(sum(j['log'] == 'LOG_MISSING' for j in jobs.values())))
        say('Neoverene vstupni seznamy: {}'.format(sum(j['inputs'] not in ('INPUT_LIST_OK','NOT_CHECKED') for j in jobs.values())))
    say('Kandidati na opakovani celeho jobu: {} {}'.format(len(candidates),
        str(candidates[:30]) + (' ...' if len(candidates) > 30 else '')))
    say('Podrobnosti: jobs.tsv, pico_status.tsv a log_issues.tsv v adresari reportu.')
    say('READABLE znamena citelny ROOT, nikoli potvrzene zpracovani vsech vstupu.')
    if not resume and not args.scan_only and not args.export_picodst and any(j['counts']['SKIPPED'] for j in jobs.values()):
        if args.use_menu:
            args.export_picodst = menu('Jak vyresit potvrzene vynechane picoDst?',
                ['Opakovat cely puvodni job (vychozi)',
                 'Vytvorit seznam vynechanych picoDst pro samostatne doplneni']) == 1
        elif not args.yes:
            args.export_picodst = confirm('Misto opakovani pico-only jobu exportovat potvrzene vynechane vstupy?')
    export_rows = []
    export_jobs = set()
    if args.export_picodst:
        # Only cleanly SKIPPED inputs from structurally readable, stable, inactive jobs.
        # Other failures still require replacing the whole block, avoiding double counting.
        per_url_states = collections.defaultdict(set)
        for n, url, state, _, _, _ in records:
            per_url_states[canonical(url)].add(state)
        for n in candidates:
            j = jobs[n]
            uncertain = any(j['counts'][s] for s in UNCERTAIN_PICO)
            uncertain = uncertain or any(issue[1] == 'REVIEW' for issue in j['log_meta']['issues'])
            if j['reasons'] != ['SKIPPED_PICO'] or j['root_state'] != 'READABLE' or uncertain or not queue_ok:
                continue
            skipped = [url for url, state, _, _ in j['pico_rows'] if state == 'SKIPPED']
            if any(per_url_states[canonical(u)] != {'SKIPPED'} for u in skipped):
                continue
            if any(input_owners[canonical(u)] != {n} for u in skipped):
                continue
            if observed[n]['root'] != fingerprint(roots[n]) or observed[n]['log'] != fingerprint(Path(j['log_path'])):
                continue
            export_rows.extend((n, u) for u in skipped)
            export_jobs.add(n)
        # Each run has a new directory; empty export is explicitly written, never stale.
        urls = sorted(set(u for _,u in export_rows))
        target = report / ('MissingDst_' + name + '.list')
        target.write_text(''.join(u+'\n' for u in urls))
        write_tsv(report / 'export_jobs.tsv', ['jobIndex','url'], export_rows)
        say('Export: {} potvrzenych vynechanych vstupu -> {}'.format(len(urls), target))
        say('Exportovane pico-only joby se v tomto behu neposlou znovu. Seznam se sam neodesila.')
        say('Pozdeji kombinuj bud puvodni blok + doplnek, nebo novy cely blok; ne oboji.')
    submit_ids = [n for n in candidates if n not in export_jobs]
    if args.jobs and not resume:
        chosen = select_indices(args.jobs,indices)
        submit_ids = [n for n in submit_ids if n in chosen]
    if args.use_menu and not args.scan_only and submit_ids:
        labels = ['{} | {}'.format(n,','.join(jobs[n]['reasons'])) for n in submit_ids]
        marked = menu('Vyber cele joby k opakovani',labels,multiple=True,checked=range(len(labels)))
        submit_ids = [submit_ids[i] for i in marked]
    if recovery and args.use_menu and not args.scan_only and submit_ids and args.request_memory is None:
        args.request_memory = memory_menu('Pamet pro nove Condor joby')
    (report / 'resubmit_jobs.list').write_text(''.join(str(n)+'\n' for n in submit_ids))
    (report / 'scan_config.json').write_text(json.dumps({'session': str(session) if session else None, 'session_version': version,
        'project_dir':str(base), 'require_keys':args.require_key, 'manual_jobs':sorted(manual_ids),
        'source_report':str(args.resume_report) if resume else None,
        'resume_baseline_ns':args.resume_since if resume else None,
        'recover': token if recovery else None, 'request_memory':args.request_memory,
        'output_pattern': output_pattern, 'log_pattern': log_pattern, 'checks': sorted(checks),
        'rules': str(args.log_rules) if args.log_rules else None, 'schedd': args.schedd,
        'exported_jobs': sorted(export_jobs), 'structural_check_only': True}, indent=2))
    if resume and any(n in queue for n in indices):
        say('Resume: joby pritomne ve fronte se znovu neodesilaji (vcetne HELD).')
    held_failed = False if resume else manage_held(queue, queue_ok, token, base, report, jobs, args)
    if held_failed:
        say('Akce HELD nebyla dokoncena; dalsi odesilani neprovadim.')
        return 2
    if recovery:
        if not args.scan_only and not queue_ok:
            say('Resubmit zablokovan: nepodarilo se spolehlive overit frontu.')
            return 2
        return submit_condor(submit_ids,recovery,report,base,jobs,observed,roots,args)
    if args.scan_only or not submit_ids:
        say('Kontrola dokoncena; zadny novy SUMS job nebyl odeslan.')
        return 0
    if not queue_ok:
        say('Resubmit zablokovan: nepodarilo se spolehlive overit frontu.')
        return 2
    if args.request_memory:
        raise ValueError('--request-memory meni nove joby pouze v --recover rezimu. SUMS -r by pouzil pamet ulozenou v session.')
    if args.schedd:
        say('Pozor: --schedd meni jen dotaz na frontu. SUMS resubmit provadej na puvodnim submit uzlu.')
    if not args.yes and not confirm('Odeslat {} celych jobu?'.format(len(submit_ids)), args.use_menu):
        say('Nic nebylo odeslano.')
        return 0
    if fingerprint(session) != session_stamp:
        raise ValueError('Session se behem kontroly zmenila. Spust kontrolu znovu.')
    # Preserve the working workaround only when the actual versions agree.
    env = os.environ.copy()
    env.pop('no_version_check', None)
    try:
        v = run([args.submit_command, '-version'], cwd=str(base), timeout=45, env=env)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ValueError('Nelze overit verzi submitteru: ' + str(exc))
    (report / 'submit_version.log').write_text(v.stdout)
    found = set(re.findall(r'Meta Scheduler\s+(\d+\.\d+\.\d+)', v.stdout))
    if v.returncode or not version or found != {version}:
        raise ValueError('Verze submitteru a session nelze potvrdit jako shodne. Viz submit_version.log.')
    env['no_version_check'] = '1'
    shutil.copy2(str(session), str(report / session.name))
    # Prevent concurrent instances of this script from resubmitting the same session.
    # flock is released by the OS on exit/crash; the lock file itself may remain.
    import fcntl
    lock = (report_base / 'resubmit.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise ValueError('Pro tuto session jiz bezi jiny resubmit z CheckMissJobs.')
    attempts = []
    failures = False
    for n in submit_ids:
        fresh_queue, reliable, note = queue_snapshot(token, base, args.schedd)
        if not reliable:
            say('Zastavuji: nelze znovu overit frontu: ' + note)
            failures = True
            break
        if n in fresh_queue:
            say('Job {} je nyni ve fronte; vynechavam.'.format(n))
            attempts.append((n, 'SKIPPED_IN_QUEUE', '', ''))
            continue
        logpath = Path(jobs[n]['log_path']) if jobs[n]['log_path'] else None
        if observed[n]['root'] != fingerprint(roots[n]) or (logpath and observed[n]['log'] != fingerprint(logpath)):
            say('Job {}: vystup/log se od kontroly zmenil; vynechavam.'.format(n))
            attempts.append((n, 'SKIPPED_CHANGED', '', ''))
            continue
        if logpath and logpath.is_file():
            oldlogs = report / 'previous_logs'
            oldlogs.mkdir(exist_ok=True)
            shutil.copy2(str(logpath), str(oldlogs / ('{}_'.format(n) + logpath.name)))
        command = [args.submit_command, '-r', str(n), str(session)]
        say('Resubmit {}: {}'.format(n, ' '.join(shlex.quote(s) for s in command)))
        log = report / ('resubmit_{}.log'.format(n))
        # Do not kill a submitter on a short timeout: acceptance could already have occurred.
        with log.open('w') as stream:
            process = subprocess.run(command, cwd=str(base), env=env, stdout=stream, stderr=subprocess.STDOUT)
        output = log.read_text(errors='replace')
        dispatched = re.search(r'Dispatching process[^\n]*' + re.escape(token) + '_' + str(n) +
                               r'(?!\d).*?\bdone\.', output, re.S | re.I)
        success = process.returncode == 0 and bool(dispatched)
        state = 'DISPATCH_REPORTED' if success else 'SUBMIT_UNCONFIRMED'
        attempts.append((n, state, process.returncode, str(log)))
        write_tsv(report / 'resubmit_status.tsv', ['jobIndex','state','exit_code','log'], attempts)
        say('  {} (viz {})'.format(state, log.name))
        if not success:
            say('Zastavuji dalsi odesilani. Over frontu a log; nejasny pokus automaticky neopakuji.')
            failures = True
            break
    write_tsv(report / 'resubmit_status.tsv', ['jobIndex','state','exit_code','log'], attempts)
    say('Hotovo. Vystup opakovaneho jobu nahrazuje puvodni blok; nemerguj obe verze.')
    return 2 if failures else 0


try:
    sys.exit(main())
except KeyboardInterrupt:
    say('\nPreruseno. Pokud bezelo odesilani, pred opakovanim over frontu a resubmit log.')
    sys.exit(130)
except (OSError, ValueError, ET.ParseError, subprocess.TimeoutExpired) as exc:
    say('CHYBA: ' + str(exc))
    sys.exit(1)
CHECKMISS_PYTHON
