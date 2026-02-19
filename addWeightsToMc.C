// SlimJetsFile.C
// Použití:
//   root -l -q 'SlimJetsFile.C("Input.root","Output_slim.root")'

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TString.h"
#include "TMath.h"
#include <iostream>
#include <cstring>
#include "TH1D.h"

// seznam větví, které chceme ve stromu "jets" zachovat
static const char* neededBranches[] = {
  "centralityAlt",
  "weightCentrality",
  "mcJetPt",
  "mcJetLambda1_1",
  "mcJetLambda1_1_5",
  "mcJetLambda1_2",
  "mcJetLambda1_3",
  "mcJetLambda1_0_5",
  "mcJetMomDisp",
  "mcJetD0Z",
  "mcD0Pt",
  "mcJetEta",
  "recoJetPt",
  "recoJetPtCorr",
  "recoJetRho",
  "recoJetNConst",
  "recoJetLambda1_1",
  "recoJetLambda1_1_5",
  "recoJetLambda1_2",
  "recoJetLambda1_3",
  "recoJetLambda1_0_5",
  "recoJetMomDisp",
  "recoJetD0Z",
  "recoJetEta",
  "recoJetArea",
  "ICS_recoJetRho",
  "ICS_recoJetPt",
  "ICS_recoJetRho",
  "ICS_recoJetNConst",
  "ICS_recoJetLambda1_1",
  "ICS_recoJetLambda1_1_5",
  "ICS_recoJetLambda1_2",
  "ICS_recoJetLambda1_3",
  "ICS_recoJetLambda1_0_5",
  "ICS_recoJetMomDisp",
  "ICS_recoJetD0Z",
  "ICS_recoJetEta",
  "ICS_recoJetArea",
  "ICS_recoJetRho"
};
static const int nNeededBranches = sizeof(neededBranches)/sizeof(neededBranches[0]);

void CopyDirSlimmed(TDirectory* source, TDirectory* dest)
{
  if (!source || !dest) return;

  source->cd();
  TIter nextKey(source->GetListOfKeys());
  TKey* key = nullptr;

  while ((key = (TKey*)nextKey())) {

    // jets: jen poslední cyklus (ponechávám, jak to máš)
    if (std::strcmp(key->GetName(), "jets") == 0) {
      TKey* lastJetsKey = (TKey*)source->GetKey("jets");
      if (key != lastJetsKey) continue;
    }

    TObject* obj = key->ReadObj();
    if (!obj) continue;

    if (obj->InheritsFrom("TDirectory")) {
      TDirectory* srcSubDir = (TDirectory*)obj;
      dest->cd();
      TDirectory* dstSubDir = dest->mkdir(srcSubDir->GetName(), srcSubDir->GetTitle());
      CopyDirSlimmed(srcSubDir, dstSubDir);
    }
    else if (obj->InheritsFrom("TTree") && TString(obj->GetName()) == "jets") {
      TTree* oldTree = (TTree*)obj;

      // vypnout všechny větve, pak zapnout jen ty potřebné
      oldTree->SetBranchStatus("*", 0);
      for (int i = 0; i < nNeededBranches; ++i) {
        oldTree->SetBranchStatus(neededBranches[i], 1);
      }

      // --- čtecí proměnné pro filtr + pT ---
      Float_t mcJetEta   = 0.f;
      Float_t recoJetEta = 0.f;
      Float_t mcJetPt    = 0.f;

      if (oldTree->GetBranch("mcJetEta"))   oldTree->SetBranchAddress("mcJetEta", &mcJetEta);
      if (oldTree->GetBranch("recoJetEta")) oldTree->SetBranchAddress("recoJetEta", &recoJetEta);
      if (oldTree->GetBranch("mcJetPt"))    oldTree->SetBranchAddress("mcJetPt",  &mcJetPt);

      // ---- PASS 1: naplnit dN/dpT po filtru ----
      // Nastav si rozumný rozsah; případně uprav podle tvých dat
      const int    nPtBins = 200;
      const double ptMin   = 0.0;
      const double ptMax   = 30.0;

      TH1D hPt("hPt_tmp", "mcJetPt after eta filter;mcJetPt;counts", nPtBins, ptMin, ptMax);
      hPt.Sumw2();

      Long64_t nEntries = oldTree->GetEntries();
      for (Long64_t i = 0; i < nEntries; ++i) {
        oldTree->GetEntry(i);

        if (TMath::Abs(mcJetEta) > 0.6 || TMath::Abs(recoJetEta) > 0.6) continue;
        if (mcJetPt > 30) continue;

        // počítáme tvar mcJetPt
        hPt.Fill(mcJetPt);
      }

      // ---- z histogramu udělat váhy w(pt)=1/N(bin), normalizace na <w>=1 ----
      std::vector<double> wPt(nPtBins + 2, 1.0); // ROOT bins: 0..n+1
      double sumW = 0.0;
      int nNonZero = 0;

      for (int b = 1; b <= nPtBins; ++b) {
        const double c = hPt.GetBinContent(b);
        if (c > 0.0) {
          wPt[b] = 1.0 / c;
          sumW += wPt[b];
          nNonZero++;
        } else {
          wPt[b] = 0.0; // prázdné biny nedávat nekonečno
        }
      }

      const double meanW = (nNonZero > 0) ? (sumW / nNonZero) : 1.0;
      if (meanW > 0.0) {
        for (int b = 1; b <= nPtBins; ++b) {
          if (wPt[b] > 0.0) wPt[b] /= meanW; // <w> ~ 1
        }
      }

      // ---- PASS 2: uložit slim tree + přidat branch weightFlatPt ----
      dest->cd();
      TTree* newTree = oldTree->CloneTree(0);

      Float_t weightFlatPt = 1.f;
      TBranch* bWeightFlat = newTree->Branch("weightFlatPt", &weightFlatPt, "weightFlatPt/F");

      for (Long64_t i = 0; i < nEntries; ++i) {
        oldTree->GetEntry(i);

        if (TMath::Abs(mcJetEta) > 0.6 || TMath::Abs(recoJetEta) > 0.6) continue;

        // binning přesně jako histogram (stejné ptMin/ptMax/nPtBins)
        int bin = hPt.FindBin(mcJetPt);
        if (bin < 1) bin = 1;
        if (bin > nPtBins) bin = nPtBins;

        weightFlatPt = (Float_t)wPt[bin];
        newTree->Fill();
      }

      newTree->Write("jets", TObject::kOverwrite);

      // vrátit status větví
      oldTree->SetBranchStatus("*", 1);

    }
    else {
      dest->cd();
      obj->Write();
    }

    if (!obj->InheritsFrom("TDirectory")) delete obj;
  }
}

void addWeightsToMc(const char* inFileName  = "Output_sim_final_01022026.root",
                  const char* outFileName = "Output_sim_final_01022026_slim.root")
{
  TFile* fin = TFile::Open(inFileName, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Nemohu otevrit vstupni soubor: " << inFileName << std::endl;
    return;
  }

  TFile* fout = TFile::Open(outFileName, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Nemohu vytvorit vystupni soubor: " << outFileName << std::endl;
    fin->Close();
    return;
  }

  CopyDirSlimmed(fin, fout);

  fout->Close();
  fin->Close();

  std::cout << "Hotovo. Zmenseny soubor ulozen jako: " << outFileName << std::endl;
}

