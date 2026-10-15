#!/bin/csh -f

echo "Are you sure that you want to delete all logs, reports and production root files? (y/n): "
set yn = $<
while ("$yn" != "y"&&"$yn" != "n")
echo "Hmm?"
set yn = $<
end

if ( "$yn" == "y" ) then
    set failed = 0

    find ./csh ./err ./list ./log ./out ./report ./production \
        -mindepth 1 -maxdepth 1 ! -name '.*' ! -type d \
        -exec /bin/rm -f -- {} +
    if ( $status != 0 ) set failed = 1

    find . -maxdepth 1 ! -type d \
        \( -name 'LocalLibraries.zip' \
        -o -name 'sched*.dataset' \
        -o -name 'sched*.session.xml' \
        -o -name '*dataset.tmp' \
        -o -name 'D0JetAnalysis*.dataset' \
        -o -name 'D0JetAnalysis*.session.xml' \) \
        -exec /bin/rm -f -- {} +
    if ( $status != 0 ) set failed = 1

    find . -mindepth 1 -maxdepth 1 \
        \( -name 'D0JetAnalysis*.package' \
        -o -name 'LocalLibraries.package' \) \
        -exec /bin/rm -rf -- {} +
    if ( $status != 0 ) set failed = 1

    if ( $failed == 0 ) then
        echo "It is cleaned up."
    else
        echo "Cleanup encountered errors. Check the messages above."
    endif
else
    echo "Ok, I'll leave it alone"
endif
