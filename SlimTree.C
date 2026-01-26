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

// seznam větví, které chceme ve stromu "jets" zachovat
static const char* neededBranches[] = {
  "centrality",
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
  "gRefMult",
  "mcJetEta",
  "mcSmearedD0Pt",
  "mcSmearedJetPt",
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
  "recoJetRho",
  "mcSmearedD0Eta",
  "mcSmearedJetEta"
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

      // --- čtecí proměnné pro filtr ---
      Float_t mcJetEta   = 0.f;
      Float_t recoJetEta = 0.f;

      // nastav adresy jen pokud větev existuje a je aktivní
      if (oldTree->GetBranch("mcJetEta"))
        oldTree->SetBranchAddress("mcJetEta", &mcJetEta);
      if (oldTree->GetBranch("recoJetEta"))
        oldTree->SetBranchAddress("recoJetEta", &recoJetEta);

      dest->cd();
      TTree* newTree = oldTree->CloneTree(0); // klon struktury, Fill bere hodnoty z oldTree po GetEntry

      Long64_t nEntries = oldTree->GetEntries();
      for (Long64_t i = 0; i < nEntries; ++i) {
        oldTree->GetEntry(i);

        // filtr: vyhoď, pokud aspoň jeden je mimo akceptanci
        if (TMath::Abs(mcJetEta) > 0.6 || TMath::Abs(recoJetEta) > 0.6)
          continue;

        // NIC NEMĚNIT: centrality, centralityAlt, ani nic dalšího
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

void SlimTree(const char* inFileName  = "Output2HardCF_corr_23012026.root",
                  const char* outFileName = "Output_sim_2HardCF_corr_26012026.root")
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

