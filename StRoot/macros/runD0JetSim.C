//#ifndef __CINT__
#include "TSystem.h"
#include "TFile.h"
#include "TROOT.h"
#include <iostream>


// basic STAR classes
class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;
////class StRefMultCorr;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0 * TMath::Pi();

enum EJetType_t{
        kFullJet, // tracks + clusters
        kChargedJet,
        kNeutralJet
};
    
// jet algorithm
enum EJetAlgo_t    {
        kt_algorithm = 0,     // background jets
        antikt_algorithm = 1, // signal jets
        cambridge_algorithm = 2,
        genkt_algorithm = 3,
        cambridge_for_passive_algorithm = 11,
        genkt_for_passive_algorithm = 13,
        plugin_algorithm = 99,
        undefined_jet_algorithm = 999
};
// jet recombination scheme
enum ERecoScheme_t
{
        E_scheme = 0,
        pt_scheme = 1,
        pt2_scheme = 2,
        Et_scheme = 3,
        Et2_scheme = 4,
        BIpt_scheme = 5,
        BIpt2_scheme = 6,
        WTA_pt_scheme = 7,
        WTA_modp_scheme = 8,
        external_scheme = 99
};
StChain *chain;

#include <ctime> // Pro měření času


void runD0JetSim(const Char_t *inputFile = "",
                                    // char *mcfilename = "/gpfs01/star/pwg/droy1/STAR-Workspace/D0Analysis/AuAuSimFramework/Test3.list",
                                    char *mcfilename = "./mcEventLists/Run14_Sim_TestList.list",
                                    const Char_t *outputFile = "output_D0JetSimTest", Int_t nEvents = 20000, bool doTEST = kTRUE, int pYear = 2014) ///Zmenit to v xml!!!!
{
   TString out = outputFile;
    out += ".root";

    TFile *outFile = new TFile(out.Data(), "RECREATE");
    outFile->Close();

    // Load necessary libraries and macros
    LoadLibs();
    LoadMacros();
 
    // =============================================================================== //

    // input file for tests (based on Run) - updated for new Runs as needed
    if ((pYear == 2014) && doTEST)
    {
        //inputFile = "lists/testPicoDsts2014.list";
        //inputFile = "lists/test_pico_2014.list";
        inputFile = "eventLists/Run14_SL22c.list";
        
    }
    if ((pYear == 2016) && doTEST)
        inputFile = "eventLists/testPicoDsts2016.list";
        
    std::cout << "inputFileName = " << inputFile << std::endl;


    // =============================================================================== //
 
    // create chain
    StChain *chain = new StChain();

    StRefMultCorr* grefmultCorrUtil = new StRefMultCorr("grefmult", "Run14_AuAu200_VpdMB5", "P16id");

    // create the picoMaker maker
    StPicoDstMaker *picoMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile, "picoDst");

    //Load database for BEMC
    StMessMgr *msg = StMessMgr::Instance();
    msg->SwitchOff("Could not make BEMC detector");
    St_db_Maker *dbMaker = new St_db_Maker("db","MySQL:StarDb","$STAR/StarDb","StarDb");
    StEmcADCtoEMaker *adc = new StEmcADCtoEMaker();
 
    StHIOverlayAngularities *HIOverlayMaker = new StHIOverlayAngularities("StHIOverlayAngularities", picoMaker, out.Data(), mcfilename, grefmultCorrUtil);
    
    // -------------- MC events ---------------------------
    HIOverlayMaker->setAllMcSeedsToEventId(false); //True -> seed = eventId+runId
    
    // -------------- Run cuts ---------------------------
    HIOverlayMaker->setFRunBadlist(0); //0 - 2014 Hanseul's //1 - 2014 Neil's //JetInfo.h
 
    // -------------- Event cuts -------------------------
    HIOverlayMaker->SetEventZVtxRange(-6, 6); //z vertex cut
    HIOverlayMaker->SetEventRVtxMaxCut(2); // sqrt(v_x^2 + v_y^2) cut
    HIOverlayMaker->SetEventVtxVpdVzMaxCut(3); // |V_z - V_vpd| cut
    std::set<int>* triggers = new std::set<int>();
	    triggers->insert(450005);
	    triggers->insert(450015);
	    triggers->insert(450025);
	    triggers->insert(450050);
	    triggers->insert(450060);
    HIOverlayMaker->SetEventTriggers(*triggers);
    HIOverlayMaker->SetNumberOfEventsToOverLay(1);
    
    // ---------------D0 cuts-----------------------------    
    HIOverlayMaker->SetMinMcPtD0(1);
    HIOverlayMaker->setMcAbsEtaD0(false, 1); //false = not used, abs(eta) <= value
    HIOverlayMaker->SetMcD0Mass(1.86484); //D0 = 1.86484, GAMMA = 0, PION_PLUS = 0.13957
    
    //----------------Jet constituent cuts-----------------------
    //Charged tracks
    HIOverlayMaker->SetUsePrimaryTracks(kFALSE); //Primary/global tracks
    HIOverlayMaker->SetMinJetTrackPt(0.2); //Min pT
    HIOverlayMaker->SetMaxJetTrackPt(30.0);  //Max pT
    HIOverlayMaker->SetJetTrackEtaRange(-1,1); //Eta range
    HIOverlayMaker->SetChargedPart(0.13957); //D0 = 1.86484, GAMMA = 0, PION_PLUS = 0.13957
    HIOverlayMaker->SetMcChargedPart(0.13957); //D0 = 1.86484, GAMMA = 0, PION_PLUS = 0.13957
    HIOverlayMaker->SetJetTracknHitsFit(15); //=>
    HIOverlayMaker->SetJetTracknHitsRatio(0.52); //=>
    HIOverlayMaker->SetJetTrackDCAcut(3.0); // <= cm

    //Towers
    HIOverlayMaker->SetMinJetTowerET(0.2); 
    HIOverlayMaker->setTowerBadList(0); //0 - 2014 Hanseul's //1 - 2014 Neil's //JetInfo.h
    HIOverlayMaker->setTowerCalibrEnergy(true);
    HIOverlayMaker->SetNeutralPart(0); //D0 = 1.86484, GAMMA = 0, PION_PLUS = 0.13957 //TEST
    HIOverlayMaker->SetMcNeutralPart(0);
    
    //Hadronic correction
    HIOverlayMaker->SetHadronicCorrFrac(1.0); // fractional hadronic correction
    HIOverlayMaker->setMaxDcaZHadronCorr(true, 3.0); // true = dca_z, false = dca; <= [cm]
    HIOverlayMaker->SetHadronicCorrMass(0.13957); //D0 = 1.86484, GAMMA = 0, PION_PLUS = 0.13957
    //Global/Primary, Pt, eta, nHitsFit, nHitsRato are same as for charged tracks
    //---------------------------------------
    //Jet cuts
    ////HIOverlayMaker->SetJetPhiRange(0, 2.0 * pi); //not used
    Double_t jetRadius = 0.4;
    HIOverlayMaker->SetJetRad(jetRadius);        // jet radius
    HIOverlayMaker->SetJetType(kFullJet); // jetType: kFullJet, kChargedJet, kNeutralJet
    HIOverlayMaker->SetJetAlgo(antikt_algorithm); // jetAlgo
    HIOverlayMaker->SetRecombScheme(E_scheme);    // recombination scheme
    HIOverlayMaker->SetGhostArea(0.005);                                          // ghost area (Changes: 0.005)
    HIOverlayMaker->SetPrintLevel(0); 
    HIOverlayMaker->DoTrackingEfficiency(1.);
    
    //Different level cuts                       	//mc	     //mcSmeared      //reco     //value
    HIOverlayMaker->setJetFractionNeutralToJetPt(	false, 		true, 		true,	    0.95    	); //
    HIOverlayMaker->setJetMinPt(			false, 		false, 		false,	    1       	); // >= [GeV/c]
    HIOverlayMaker->setJetMinArea(			false, 		false, 		false,	    0      	); // >= value*pi*R^2
    HIOverlayMaker->setJetMinAbsEta(			false, 		false, 		false,	 1.-jetRadius   ); // abs(value) <= 

   //Strategy for jet background subtraction
   HIOverlayMaker->SetBgSubtraction(12); //1 - Area based method + jet shape method // 2 - ICS // 12 or 21 - both
   HIOverlayMaker->setJetNHardestSkipped(2, 2); // First: 0-10%; Second: 10-80% //CHANGE
   HIOverlayMaker->SetPhiBgModulation(false);
   
   //Fastjet
   HIOverlayMaker->setJetFixedSeed(false, 12345); // false = random seed, true =

   //ICS
   std::vector<Double_t> maxDeltaRs;
   maxDeltaRs.push_back(0.100);
   maxDeltaRs.push_back(0.175);

   std::vector<Double_t> alphas;
   alphas.push_back(0.0);
   alphas.push_back(0.0);
   HIOverlayMaker->setICSSubtractionParams(maxDeltaRs, alphas);

    // initialize chain
    chain->Init();
    std::cout << "chain->Init();" << std::endl;
    int total = picoMaker->chain()->GetEntries();
    std::cout << " Total entries = " << total << std::endl;

    ////nEvents = TMath::Min(total, nEvents);
    std::cout << " Actual entries = " << nEvents << std::endl;
    
    std::time_t start_time = std::time(NULL);


    for (Int_t i = 0; i < nEvents; i++)
    {
        if (doTEST)
        {
        /*
            if (i % 100 == 0)
            {
                std::cout << "****************************************** " << std::endl;
                std::cout << "Working on eventNumber " << i << std::endl;
            }*/
               if(i%200!=0) progres(i,nEvents, start_time);
        }
        else
        {
            if (i % 500 == 0)
            {
                // std::cout << "****************************************** " << std::endl;
                std::cout << "Working on eventNumber " << i << std::endl;
            }
        }

        // std::cout << "Working on eventNumber " << i << std::endl;

        chain->Clear();
        int iret = chain->Make(i);
        if (iret)
        {
            std::cout << "Bad return code!" << iret << std::endl;
            break;
        }

        total++;
    }
       if(doTEST &&i%200!=0) progres(i,nEvents, start_time);
    std::cout << std::endl;
    std::cout << "****************************************** " << std::endl;
    std::cout << "Work done... now its time to close up shop!" << std::endl;
    std::cout << "****************************************** " << std::endl;
    chain->Finish();
    std::cout << "****************************************** " << std::endl;
    std::cout << "total number of events  " << nEvents << std::endl;
    std::cout << "****************************************** " << std::endl;

    delete chain;

    // close output file if open
    if (outFile->IsOpen())
        outFile->Close();

    // StMemStat::PrintMem("load StChain");
}

void LoadLibs()
{
    // load fastjet libraries 3.x
    // gSystem->Load("libCGAL"); - not installed
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjet");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libsiscone");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libsiscone_spherical");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjetplugins");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjettools");
    gSystem->Load("/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/lib/libfastjetcontribfragile");

    // add include path to use its functionality
    //gSystem->AddIncludePath("-I/gpfs01/star/pwg/lomicond/Ondrej/Jets/Alex_install/install/fastjet-install/include");

    // load the system libraries - these were defaults
    gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
    loadSharedLibraries();

    // these are needed for new / additional classes
    gSystem->Load("libStPicoEvent");
    gSystem->Load("libStPicoDstMaker");

gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");

  gSystem->Load("St_db_Maker");
  gSystem->Load("StDaqLib");
  
    gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StTriggerUtilities");
 // gSystem->Load("StRefMultCorr"); //
//gSystem->Load("StPicoCutsBase");
  gSystem->Load("StBTofUtil");
 // gSystem->Load("StPicoD0EventMaker");
  //gSystem->Load("StPicoHFMaker");
 // gSystem->Load("StPicoCuts");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbUtilities");
  gSystem->Load("Sti");
  gSystem->Load("StiUtilities");
  gSystem->Load("StSsdDbMaker");
  gSystem->Load("StSvtDbMaker");
  gSystem->Load("StiMaker");
  gSystem->Load("StDbBroker");
  gSystem->Load("libgeometry_Tables"); //rember, order of loading makers matters
  //gSystem->Load("StPicoD0AnaMaker");

    // my libraries
    gSystem->Load("StRefMultCorr");
    gSystem->Load("StPicoD0JetAnaMaker");

    gSystem->ListLibraries();
}

void LoadMacros()
{
}

void progres(double citatel, double jmenovatel, std::time_t start_time){
      int Ndilky=50;
      int proc=floor(citatel/jmenovatel*Ndilky);
      
      
       std::time_t now = std::time(NULL);
    double elapsed_seconds = difftime(now, start_time);
    double estimated_total_time = (elapsed_seconds / citatel) * jmenovatel;
    int remaining_time = static_cast<int>(estimated_total_time - elapsed_seconds);
      
         // Převod zbývajícího času na minuty a sekundy
    int hours = remaining_time / 3600;
    int minutes = (remaining_time / 60) % 60;
    int seconds = remaining_time % 60;
      
      cout << "\r"<< flush;
      cout << "  │" << flush;
      for (int i=1; i<=proc; i++){
      cout <<"█"<<flush; 
      }
      for (int j=proc+1; j<=Ndilky;j++){
      cout <<  "░"<<flush;
      }
      cout << "│ " << flush;
    if (citatel != jmenovatel) {
       cout << Form("Completed: %.2f%% | Remaining: %02d:%02d:%02d", citatel / jmenovatel * 100.,hours, minutes, seconds) << flush;
    } else {
       cout << Form("\033[1;32m Completed: %.2f%% \033[0m", citatel / jmenovatel * 100.) << "\033[1;32m%\033[0m" << endl;
    }
    
      return;
    }
