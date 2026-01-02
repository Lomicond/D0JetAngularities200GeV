// ################################################################
// Author:  Joel Mazer for the STAR Collaboration
// Affiliation: Rutgers University
//
// ################################################################

#include "StHIOverlayAngularities.h"
#include "StMemStat.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include <THnSparse.h>
#include "TParameter.h"
#include <TProfile.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TEfficiency.h"

// STAR includes
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoMtdTrigger.h"
#include "StRoot/StPicoEvent/StPicoBEmcPidTraits.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoMtdPidTraits.h"

// jet-framework includes
#include "StFJWrapper.h"
//#include "StEmcPosition.h"
#include "StRoot/StEmcUtil/projection/StEmcPosition.h"
#include "phys_constants.h"

// centrality includes
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
//#include "StCentMaker.h"
#include "StCuts.h"

//Towers
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcTables.h"


Double_t standardPhi(const Double_t &phi)
{
  Double_t phi_standard = phi;
  if (phi_standard < 0)
    phi_standard += 2 * (TMath::Pi());
  if (phi_standard > 2 * (TMath::Pi()))
    phi_standard -= 2 * (TMath::Pi());
  return phi_standard;
}

ClassImp(StHIOverlayAngularities)

    //________________________________________________________________________
    StHIOverlayAngularities::StHIOverlayAngularities(const char *name, StPicoDstMaker *picoMaker, const char *outName, const char *filename, StRefMultCorr* grefmultCorrUtil) : /*StJetFrameworkPicoBase(name),*/ mGRefMultCorrUtil(grefmultCorrUtil) // StMaker(name),
{
  //moved here
  if (!name)return;
  
  fRunNumber = 0;
  mPicoDstMaker = 0x0;
  mPicoDst = 0x0;
  mPicoEvent = 0x0;
  //grefmultCorr = 0x0;
  mOutName = outName;
  doUsePrimTracks = kFALSE;
  fJetRad = 0.4;
  fJetAlgo = 1;      // antikt_algorithm
  fJetType = 1;      // kChargedJet
  fRecombScheme = 0; // E_scheme
  //fJetPhiMin = 0.;
 // fJetPhiMax = 2 * pi;
  //fJetEtaMin = -0.6;
  //fJetEtaMax = 0.6;
  fGhostArea = 0.01;
  fEventZVtxMinCut = -40.0;
  fEventZVtxMaxCut = 40.0;
  fEventRVtxMaxCut = -1;
  fEventVtxVpdVzMaxCut = -1;
  //fTrackPtMinCut = 0.2;
  //fTrackPtMaxCut = 30.0;
  fTracknHitsFit = 15;
  fTracknHitsRatio = 0.52;
  fCentrality = -999.;
  fCentralityAlt = -999.;
  fCentralityWeight = 0.;
  Bfield = 0.0;
  zVtx = 0.0;
  fRhoVal = 0;
  mEmcPosition = 0x0;
  fKgrefMult_uncorr = 0;
  fKrefMult_uncorr = 0;
  fMCFileListName = filename;
  fHadrCorrTrackDCAZcut = 3.0;
  fMinMcPtD0 = -1;
  fSetMcAbsEtaD0 = false;
  fMcAbsEtaD0 = 1;
  fMcD0Mass = 0;
  fMcChargedPart =0;
  fMcNeutralPart =0;
  fMassiveAll = true;
  fgAlpha1 = 0;
  fgAlpha2 = 0;
  fPhiBgModulation = false;
  fQ_1 = -999;
  fQ_2 = -999;
  fPsi_2 = -999;
  fPsi_2_shifted = -999;
  fQ_1_rec = -999;
  fQ_2_rec = -999;
  fBgSubtraction = 1;
  

  //for (Int_t i = 0; i < 8; i++) fEmcTriggerArr[i] = 0;

  for (Int_t i = 0; i < 4800; i++)
  {
    for (Int_t j = 0; j < 7; j++)
      mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;
  }

  mBemcGeom = 0x0;

  fMinJetTrackPt = 0.2;
  fMaxJetTrackPt = 30.0;
  fJetTrackEtaMin = -1.0;
  fJetTrackEtaMax = 1.0;
  fJetTrackPhiMin = 0.0;
  fJetTrackPhiMax = 2.0 * TMath::Pi();
  fJetTrackDCAcut = 3.0;
  fJetTracknHitsFit = 15;
  fJetTracknHitsRatio = 0.52;
  //fJetNeutralPtFraction = 0.95;
  fJetTowerEtaMin = -1.0;
  fJetTowerEtaMax = 1.0;
  fJetTowerPhiMin = 0.0;
  fJetTowerPhiMax = 2.0 * TMath::Pi();
  mTowerEnergyTMin = 0.2;
  mHadronicCorrFrac = 1.;

	fSetJetFracCut[0] = kTRUE;
	fSetJetFracCut[1] = kTRUE;
	fSetJetFracCut[2] = kTRUE;
	
	fSetJetMinPtCut[0] = kTRUE;
	fSetJetMinPtCut[1] = kTRUE;
	fSetJetMinPtCut[2] = kTRUE;
	
	fSetJetMinAreaCut[0] = kTRUE;
	fSetJetMinAreaCut[1] = kTRUE;
	fSetJetMinAreaCut[2] = kTRUE;
	
	fSetJetMinAbsEtaCut[0] = kTRUE;
	fSetJetMinAbsEtaCut[1] = kTRUE;
	fSetJetMinAbsEtaCut[2] = kTRUE;	
		
  fNumberOfEventsToOverLay = 1;

  for (Int_t i = 0; i < 100; i++)
  {
    fMcEventTracks[i].clear();
    fMcEventTowers[i].clear();
    fRecoMcEventTracks[i].clear();
    fRecoMcEventTowers[i].clear();
    fMcD0Information[i] = {};
    fMcRecoD0Information[i] = {};
    fOrigin[i].SetXYZ(0, 0, 0);
    fMcEventInfo[i] = {};
  }

  fPrintLevel = 0;
  fCentBin = 0;
  fTrackingEfficiency = kFALSE;
  fTrackingEfficiencyPercentage = 1.;
  

  
  SetName(name);
}

//
//________________________________________________________________________
StHIOverlayAngularities::~StHIOverlayAngularities()
{ /*  */
}

//
//________________________________________________________________________
Int_t StHIOverlayAngularities::Init()
{

  ////StJetFrameworkPicoBase::Init();
  fastjet::ClusterSequence::print_banner();
  // position object for Emc
  mBemcGeom = StEmcGeom::instance("bemc");
  mEmcPosition = new StEmcPosition();

  //Loading of BEMC tables
  //StMaker* maker = GetMaker("Eread");
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  ifstream filelistforMCEvents(fMCFileListName.Data());
/*
  // get base class pointer - this class does not inherit from base class: StJetFrameworkPicoBase, but we want to reduce redundancy
  mBaseMaker = static_cast<StJetFrameworkPicoBase *>(GetMaker("baseClassMaker"));
  if (!mBaseMaker)
  {
    LOG_WARN << " No baseMaker! Skip! " << endm;
    return kStWarn;
  }
*/
  // get bad run, dead & bad tower lists
/*
  badRuns = mBaseMaker->GetBadRuns();
  deadTowers = mBaseMaker->GetDeadTowers();
  badTowers = mBaseMaker->GetBadTowers();
*/
  if (!filelistforMCEvents.is_open())
  {
    LOG_ERROR << "No MC File List! Exiting!" << endm;
    return kStOk;
  }

  string line;

  while (getline(filelistforMCEvents, line))
  {
    TString s(line);
    filenamesforHIOverlay.push_back(s);
  }

  OutputTreeInit();

  // histograms
  //=======================================================================================================//
  	// Event histograms:
	hVtxZ = new TH1D("hVtxZ", ";PVtx.z() [cm]; Count", 100, -10, 10);
	hVtxR = new TH2D("hVtxR", ";PVtx.x() [cm]; PVtx.y() [cm]", 100, -3, 3, 100, -3, 3);
	hVzDiff = new TH1D("hVzDiff", "V_{z} - V_{z}^{VPD};V_{z} - V_{z}^{VPD} [cm];Count", 80, -4.0, 4.0);
	hCentrality = new TH1D("hCentrality", ";C_{ID}", 9, -0.5, 8.5);
	hCentralityW = new TH1D("hCentralityW", ";C_{ID}", 9, -0.5, 8.5);
	
	//Mc event
	hMcJetConstMom            = new TH2D("hMcJetConstMom", "Mc D0 jet constituent momementum, w/o D^{0};p_{T_Jet} [GeV/c];p_{T,i} [GeV/c]", 100, 0, 40, 350, 0, 35); 
	//hMcJetConstMom5GeV        = new TH2D("hMcJetConstMom5GeV", "Mc D0 jet (p_{T,Jet} > 5 GeV/c) constituent momentum, w/o D^{0};p_{z} [GeV/c];p_{T} [GeV/c]", 100, -30, 30, 350, 0, 35); 
	hMcJetConstTheta	= new TH2D("hMcJetConstTheta", "Mc D0 jet constituent theta, w/o D^{0};p_{T_Jet} [GeV/c];#theta_{i}", 100, 0, 40, 300, 0, TMath::Pi());
	hPureMcNeutralEtaPhi         = new TH2D("hPureMcNeutralEtaPhi",      "Neutral particles from MC w/o bad tower cut #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1);
	
	//Jet constituents:
	hJetTracksDedx           = new TH2D("hJetTracksDedx",           "Track dE/dx vs signed |p|; sign(q) #times |p| [GeV/c];dE/dx [keV/cm]", 200, -4, 4, 400, 0, 20);
	hJetTracksDedxAfterCuts  = new TH2D("hJetTracksDedxAfterCuts",  "Track dE/dx vs signed |p| (after cuts);sign(q) #times |p| [GeV/c];dE/dx [keV/cm]", 400, -4, 4, 200, 0, 20);
	hJetTracksPt             = new TH1D("hJetTracksPt",         "Jet charged constituent p_{T};p_{T} [GeV/c];Count", 200, 0, 40);
	hJetTracksEtaPhi         = new TH2D("hJetTracksEtaPhi",      "Jet charged constituents #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1);
	hJetTracksNHitsFit       = new TH1D("hJetTracksNHitsFit",     "Jet charged constituents nHitsFit;nHitsFit;Count", 51, -0.5, 50.5);
	hJetTracksNHitsRatio     = new TH1D("hJetTracksNHitsRatio",   "Jet charged constituents nHitsRatio;nHitsRatio;Count", 101, 0.0, 1.01);
	hJetTracksDCA            = new TH1D("hJetTracksDCA",           "Jet charged constituents DCA;DCA [cm];Count", 100, 0, 4);
	 
	hJetNeutralPt         = new TH1D("hJetNeutralPt",        "Jet neutral constituent p_{T};p_{T} [GeV/c];Count", 200, 0, 40);
	hJetNeutralEtaPhi     = new TH2D("hJetNeutralEtaPhi",    "Jet neutral constituent #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1);
	hJetNeutralEtBefAftHC = new TH2D("hJetNeutralEtBefAftHC", "Tower E_{T} before vs after HC;E_{T}^{before} [GeV];E_{T}^{after} [GeV]", 350, 0, 35, 350, 0, 35);
	hJetNeutralECalibBefAft = new TH2D("hJetNeutralECalibBefAft","Tower E before vs after calibration;E^{before calib} [GeV];E^{after calib} [GeV]", 350, 0, 35, 350, 0, 35);
	
	hJetConstCharge = new TH1D("hJetConstCharge", "Jet constituent charge;charge;Count", 3, -1.5, 1.5);
	hJetConstRapPhi = new TH2D("hJetConstRapPhi", "Jet constituent (w/o D^{0}) y vs #phi;#phi;rapidity", 120, -TMath::Pi(), TMath::Pi(), 240, -1.2, 1.2);
	hJetConstRapPhiICS = new TH2D("hJetConstRapPhiICS", "Jet constituent (w/o D^{0}) y vs #phi after ICS;#phi;rapidity", 120, -TMath::Pi(), TMath::Pi(), 240, -1.2, 1.2);
	hJetConstEtaPhi = new TH2D("hJetConstEtaPhi", "Jet constituent (w/o D^{0}) #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 200, -5, 5);
	hJetConstEtaPhiICS = new TH2D("hJetConstEtaPhiICS", "Jet constituent (w/o D^{0}) #eta vs #phi after ICS;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 200, -5, 5);
	hJetConstPt = new TH1D("hJetConstPt", "Jet constituent (w/o D^{0}); p_{T} [GeV/c];Count", 350, 0, 35);
	hJetConstPtICS = new TH1D("hJetConstPtICS", "Jet constituent (w/o D^{0}) after ICS; p_{T} [GeV/c];Count", 350, 0, 35);
	hJetMcRecoTracksPt             = new TH1D("hJetMcRecoTracksPt",         "Jet mc reco charged constituent p_{T};p_{T} [GeV/c];Count", 200, 0, 40);
	hJetMcRecoTracksEtaPhi         = new TH2D("hJetMcRecoTracksEtaPhi",      "Jet mc reco charged constituents #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1);

	hJetMcRecoTracksNHitsFit       = new TH1D("hJetMcRecoTracksNHitsFit",     "Jet mc reco charged constituents nHitsFit;nHitsFit;Count", 51, -0.5, 50.5);
	hJetMcRecoTracksNHitsRatio     = new TH1D("hJetMcRecoTracksNHitsRatio",   "Jet mc reco charged constituents nHitsRatio;nHitsRatio;Count", 101, 0.0, 1.01);
	hJetMcRecoTracksDCA            = new TH1D("hJetMcRecoTracksDCA",           "Jet mc reco charged constituents DCA;DCA [cm];Count", 100, 0, 4);
  
  	hJetMcRecoNeutralPt         = new TH1D("hJetMcRecoNeutralPt",        "Jet mc reco neutral constituent p_{T};p_{T} [GeV/c];Count", 200, 0, 40);
	hJetMcRecoNeutralEtaPhi     = new TH2D("hJetMcRecoNeutralEtaPhi",    "Jet mc reco neutral constituent #eta vs #phi;#phi;#eta", 120, -TMath::Pi(), TMath::Pi(), 100, -1.1, 1.1);
	hJetMcRecoNeutralEtBefAftHC = new TH2D("hJetMcRecoNeutralEtBefAftHC", "Mc reco tower E_{T} before vs after HC;E_{T}^{before} [GeV];E_{T}^{after} [GeV]", 350, 0, 35, 350, 0, 35);
        
        hJetMcRecoConstCharge = new TH1D("hJetMcRecoConstCharge", "Jet mc reco constituent charge;charge;Count", 3, -1.5, 1.5);
        hFractionNeutralToJetPt = new TH1D("hFractionNeutralToJetPt", "; Fraction of Neutral to Jet p_{t}; Counts", 200, 0, 1);
        
        effJetRec00_10 = new TEfficiency("effJetRec00_10","Percentage of rec. mc jets C = 00-10%; Jet p_{T} [GeV/c];Efficiency",200,0,40);
        effJetRec10_40 = new TEfficiency("effJetRec10_40","Percentage of rec. mc jets C = 10-40%; Jet p_{T} [GeV/c];Efficiency",200,0,40);
        effJetRec40_80 = new TEfficiency("effJetRec40_80","Percentage of rec. mc jets C = 40-80%; Jet p_{T} [GeV/c];Efficiency",200,0,40);


	etaRecoTrue00_10 = new TH2D("etaRecoTrue00_10","D0 jet eta reco vs true; #eta^{reco}_{Jet};#eta^{true}_{Jet}",300,-1.5,1.5,300,-1.5,1.5);
	etaRecoTrue10_40 = new TH2D("etaRecoTrue10_40","D0 jet eta reco vs true; #eta^{reco}_{Jet};#eta^{true}_{Jet}",300,-1.5,1.5,300,-1.5,1.5);
	etaRecoTrue40_80 = new TH2D("etaRecoTrue40_80","D0 jet eta reco vs true; #eta^{reco}_{Jet};#eta^{true}_{Jet}",300,-1.5,1.5,300,-1.5,1.5);

        //Neutral particles hadronic correction:
	hJetHadrCorrNHitsFit   = new TH1D("hJetHadrCorrNHitsFit",   "Hadron. corr. track nHitsFit;nHitsFit;Count", 51, -0.5, 50.5);
	hJetHadrCorrNHitsRatio = new TH1D("hJetHadrCorrNHitsRatio", "Hadron. corr. track nHitsRatio;nHitsRatio;Count", 101, 0.0, 1.01);
	hJetHadrCorrDcaZ       = new TH1D("hJetHadrCorrDcaZ",       "Hadron. corr. track DCA_{z};DCA_{z} [cm];Count", 100, 0, 4);
	hJetHadrCorrEtaVsPt    = new TH2D("hJetHadrCorrEtaVsPt",    "Hadron. corr. track #eta vs p_{T};p_{T} [GeV/c];#eta", 200, 0, 40, 150, -1.5, 1.5);
	hJetHadrCorrE          = new TH1D("hJetHadrCorrE",    "Hadron. corr. track energy;E [GeV];Count", 200, 0, 40);
	
	hResponseJetPt             = new TH2D("hResponseJetPt", "hResponseJetPt; p_{T,Jet}^{reco} [GeV/c]; p_{T,Jet}^{true} [GeV/c]", 200, 0, 50, 200, 0, 50);
	hResponseJetD0Z            = new TH2D("hResponseJetD0Z", "hResponseJetD0Z; z^{reco}; z^{true}", 200, 0, 1.01, 200, 0, 1.01);
	hResponseJetNConst         = new TH2D("hResponseJetNConst", "hResponseJetNConst; N_{const}^{reco}; N_{const}^{true}", 61, -0.5, 60.5, 61, -0.5, 60.5);
	hResponseJetLambda1_0_5    = new TH2D("hResponseJetLambda1_0_5", "hResponseJetLambda1_0_5; #lambda^{1, reco}_{0.5}; #lambda^{1, true}_{0.5}", 200, 0, 1.5, 200, 0, 1.5);
	hResponseJetLambda1_1      = new TH2D("hResponseJetLambda1_1", "hResponseJetLambda1_1; #lambda^{1, reco}_{1}; #lambda^{1, true}_{1}", 200, 0, 1.5, 200, 0, 1.5);
	hResponseJetLambda1_1_5    = new TH2D("hResponseJetLambda1_1_5", "hResponseJetLambda1_1_5; #lambda^{1, reco}_{1.5}; #lambda^{1, true}_{1.5}", 200, 0, 1.5, 200, 0, 1.5);
  	hResponseJetLambda1_2      = new TH2D("hResponseJetLambda1_2", "hResponseJetLambda1_2; #lambda^{1, reco}_{2}; #lambda^{1, true}_{2}", 200, 0, 1.5, 200, 0, 1.5);
  	hResponseJetLambda1_3      = new TH2D("hResponseJetLambda1_3", "hResponseJetLambda1_3; #lambda^{1, reco}_{3}; #lambda^{1, true}_{3}", 200, 0, 1.5, 200, 0, 1.5);
  	hResponseJetMomDisp  = new TH2D("hResponseJetMomDisp", "hResponseJetMomDisp; p^{D, reco}_{T}; p^{D, true}_{T}", 200, 0, 1.01, 200, 0, 1.01);
  	hResponseJetD0DeltaR  = new TH2D("hResponseJetD0DeltaR", "hResponseJetD0DeltaR; #Delta R_{D^{0}}^{reco}; #Delta R_{D^{0}}^{true}", 200, 0, 5, 200, 0, 5);
  	hResponseJetD0Pt     = new TH2D("hResponseJetD0Pt", "D0 response in jet;   p_{T}^{D^{0}, reco} [GeV/c]; p_{T}^{D^{0}, true} [GeV/c]", 200, 0, 12, 200, 0, 12);

  //=======================================================================================================//

//  TFile f("/star/u/droy1/Y2019/STAR/Momentum_resolution_SL16d.root");
  TFile *f;
  f = new TFile("StRoot/StPicoD0JetAnaMaker/Files/Momentum_resolution_Run14.root");

  fPionMomResolution = (TF1 *)f->Get("fPion")->Clone("fPion");
  fKaonMomResolution = (TF1 *)f->Get("fKaon")->Clone("fKaon");
  fProtonMomResolution = (TF1 *)f->Get("fProton")->Clone("fProton");

//  TFile effweight("/star/u/droy1/Y2019/STAR/EffWeightsInCentralityBins.root");
  TFile effweight("StRoot/StPicoD0JetAnaMaker/Files/EffWeightsInCentralityBins.root");
  fPionWeight[0] = (TGraph *)effweight.Get("Pion_0_10");
  fKaonWeight[0] = (TGraph *)effweight.Get("Kaon_0_10");
  fProtonWeight[0] = (TGraph *)effweight.Get("Proton_0_10");
  fAProtonWeight[0] = (TGraph *)effweight.Get("AProton_0_10");

  fPionWeight[1] = (TGraph *)effweight.Get("Pion_10_40");
  fKaonWeight[1] = (TGraph *)effweight.Get("Kaon_10_40");
  fProtonWeight[1] = (TGraph *)effweight.Get("Proton_10_40");
  fAProtonWeight[1] = (TGraph *)effweight.Get("AProton_10_40");

  fPionWeight[2] = (TGraph *)effweight.Get("Pion_40_80");
  fKaonWeight[2] = (TGraph *)effweight.Get("Kaon_40_80");
  fProtonWeight[2] = (TGraph *)effweight.Get("Proton_40_80");
  fAProtonWeight[2] = (TGraph *)effweight.Get("AProton_40_80");

  // ============================ set jet parameters for fastjet wrapper  =======================
  // recombination schemes:
  fastjet::RecombinationScheme recombScheme;
  if (fRecombScheme == 0)
    recombScheme = fastjet::E_scheme;
  if (fRecombScheme == 1)
    recombScheme = fastjet::pt_scheme;
  if (fRecombScheme == 2)
    recombScheme = fastjet::pt2_scheme;
  if (fRecombScheme == 3)
    recombScheme = fastjet::Et_scheme;
  if (fRecombScheme == 4)
    recombScheme = fastjet::Et2_scheme;
  if (fRecombScheme == 5)
    recombScheme = fastjet::BIpt_scheme;
  if (fRecombScheme == 6)
    recombScheme = fastjet::BIpt2_scheme;
  if (fRecombScheme == 7)
    recombScheme = fastjet::WTA_pt_scheme;
  if (fRecombScheme == 8)
    recombScheme = fastjet::WTA_modp_scheme;
  if (fRecombScheme == 99)
    recombScheme = fastjet::external_scheme;

  // jet algorithm
  fastjet::JetAlgorithm algorithm;
  if (fJetAlgo == 1)
    algorithm = fastjet::antikt_algorithm;
  if (fJetAlgo == 0)
    algorithm = fastjet::kt_algorithm;
  // extra algorithms
  if (fJetAlgo == 2)
    algorithm = fastjet::cambridge_algorithm;
  if (fJetAlgo == 3)
    algorithm = fastjet::genkt_algorithm;
  if (fJetAlgo == 11)
    algorithm = fastjet::cambridge_for_passive_algorithm;
  if (fJetAlgo == 13)
    algorithm = fastjet::genkt_for_passive_algorithm;
  if (fJetAlgo == 99)
    algorithm = fastjet::plugin_algorithm;
  if (fJetAlgo == 999)
    algorithm = fastjet::undefined_jet_algorithm;
  fastjet::Strategy strategy = fastjet::Best;

  // setup fj wrapper
  fjw = new StFJWrapper("HIOverlay", "HIOverlay");
  fjw->SetHJetConstRapPhiICS(hJetConstRapPhiICS);
  fjw->SetHJetConstEtaPhiICS(hJetConstEtaPhiICS);
  fjw->SetHJetConstPtICS(hJetConstPtICS);
  fjw->SetBackgroundSub(kTRUE);
  fjw->SetSubtractionMc(kTRUE);
  fjw->SetAreaType(fastjet::active_area_explicit_ghosts);
  fjw->SetStrategy(strategy);
  fjw->SetGhostArea(fGhostArea);
  fjw->SetR(fJetRad); //Checked
  fjw->SetAlgorithm(algorithm);       // fJetAlgo); //Checked
  fjw->SetRecombScheme(recombScheme); // fRecombScheme);
  fjw->SetMaxRap(10.);                // tracks


  // // ghost-area specifications
  // Double_t ghost_maxrap = 1.2;
  // fastjet::GhostedAreaSpec area_spec(ghost_maxrap);
  // fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);

////////////////////////////////////////////////////
cout << "----------------------------" << endl;
cout << (fSetMcSeed  == 0 ? "MC random seed" : fSetMcSeed  == 1 ? "MC seed = eventId + runId" : "Weird" ) << endl;
cout << "----------------------------" << endl;
cout << "Run cuts" << endl;
cout << (fRunBadlist  == 0 ? "Hanseul's bad run 2014 list" : fRunBadlist  == 1 ? "Neil's bad run 2014 list" : "None bad run listed used!!!") << endl;
cout << "----------------------------" << endl;
cout << "Event cuts" << endl;
cout << "V_z min =< " << fEventZVtxMinCut << " cm" << endl;
cout << "V_z max =< " << fEventZVtxMaxCut << " cm" << endl;
cout << "V_r =< " << fEventRVtxMaxCut << " cm" << endl;
cout << "|V_z-V_{z,VPD}| =< " << fEventVtxVpdVzMaxCut << " cm" << endl;
cout << "V_x !=0 && V_y !=0 && V_z !=0" << endl;
cout << "Triggers: {";
bool first = true;
for (const Int_t &trigger : fEventTriggers) {
    if (!first) cout << ", ";
    cout << trigger;
    first = false;
}
cout << "}" << endl;
cout << "----------------------------" << endl;
cout << "Jet constituent cuts" << endl;
cout << "********************" << endl;
cout << "  D0 meson" << endl;
cout << "MC particle mass: " << fMcD0Mass << " GeV/c^2"   << endl;
cout << "pT_min => " << fMinMcPtD0 << " GeV/c" << endl; 
TString etamaxabs;
if (fSetMcAbsEtaD0) etamaxabs=TString::Format("%.2f", fMcAbsEtaD0); else etamaxabs = "INF";
cout << "|eta| =< " << etamaxabs << endl; 
cout << "********************" << endl;
cout << "  Charged tracks" << endl;
TString GlobPrim;
if (doUsePrimTracks) GlobPrim = "Primary tracks"; else GlobPrim = "Global tracks";
cout << GlobPrim << endl;
cout << "pT_min => " << fMinJetTrackPt << " GeV/c" << endl; 
cout << "pT_max =< " << fMaxJetTrackPt << " GeV/c" << endl; 
cout << "eta_min => " << fJetTrackEtaMin << endl; 
cout << "eta_max =< " << fJetTrackEtaMax << endl; 
cout << "nHitsFit => " << fJetTracknHitsFit << endl;
cout << "nHitsRatio => " << fJetTracknHitsRatio << endl;
cout << "Particle mass: " << fChargedPart << " GeV/c^2"   << endl;
cout << "MC particle mass: " << fMcChargedPart << " GeV/c^2"   << endl;
cout << "DCA =< " << fJetTrackDCAcut << " cm" << endl;
cout << "********************" << endl;
cout << "  Towers " << endl;
cout << "ET => " << mTowerEnergyTMin << " GeV" << endl;
cout << (fSetTowerCalibrEnergy  == true ? "Hanseul's energy calibration" : fSetTowerCalibrEnergy  == false ? "Production energy calibration" : "Weird") << endl;
cout << (fTowerBadlist  == 0 ? "Hanseul's bad tower 2014 list" : fTowerBadlist  == 1 ? "Neil's bad tower 2014 list" : "None bad tower listed used!!!") << endl;
cout << "Particle mass: " << fNeutralPart << " GeV/c^2"  << endl; 
cout << "MC particle mass: " << fMcNeutralPart << " GeV/c^2"  << endl;
cout << "********************" << endl;
cout << "  Hadronic correction " << endl;
cout << "Fraction of hadr. corr. = " << mHadronicCorrFrac << endl;
cout << GlobPrim << endl;
cout << "pT_min => " << fMinJetTrackPt << " GeV/c" << endl; 
cout << "pT_max =< " << fMaxJetTrackPt << " GeV/c" << endl; 
cout << "eta_min => " << fJetTrackEtaMin << endl; 
cout << "eta_max =< " << fJetTrackEtaMax << endl; 
cout << "nHitsFit => " << fJetTracknHitsFit << endl;
cout << "nHitsRatio => " << fJetTracknHitsRatio << endl;
cout << "Particle mass: " << fHadronicCorrMass << " GeV/c^2" << endl;
//cout << "DCA_z =< " << fHadrCorrTrackDCAZcut << " cm" << endl;
cout << (fSetDcaZHadronCorr  == true ? "DCA_Z" : fSetDcaZHadronCorr  == false ? "DCA" : "Weird") << " =< " << fDcaZHadronCorr << " cm" << endl; 
cout << endl;
cout << "----------------------------" << endl;
cout << "  Jets " << endl;

cout << "**********" << endl;
cout << "Jet type: " << (fJetType == 0 ? "Full jets" : fJetType == 1 ? "Charged jets" : fJetType == 2 ? "Neutral jets" : "?")  << endl;
cout << "Jet rec. algorithm: " << (antikt_algorithm == 0 ? "kt" : antikt_algorithm == 1 ? "anti-kt" : "Other" ) << endl;
cout << "Jet R = " << fJetRad << endl;
cout << "Jet area min >= " << fJetMinAreaCut << endl;
cout << "|eta_jet| <= " << 1.-fJetRad << endl;
cout << "fJetFractionNeutralToJetPt: " <<fJetFractionNeutralToJetPt << endl;
cout << "**********" << endl;
cout << "  Backgroud estimation " << endl;
// 0 - Area based, 1 - ICS, 2 - jet shape
	// 0 - Area based, 1 - ICS, 2 - jet shape
	switch (fBgSubtraction) {
	    case 0:
		cout << "Area based method" << endl;
		break;

	    case 1:
		cout << "ICS method" << endl;
		break;

	    case 2:
		cout << "Jet shape method" << endl;
		break;

	    default:
		std::cerr << "Warning: Unknown fBgSubtraction option " << fBgSubtraction << std::endl;
		exit(1); // 1 bývá běžnější pro chybu
	}
cout << "Massive particles? " << (fMassiveAll == false ? "false" : fMassiveAll == true ? "true" : "Weird") << endl;
cout << "alpha_1: " << fgAlpha1 << " alpha_2: " << fgAlpha2 << endl;
cout << "Background phi modulation: " <<  (fPhiBgModulation == false ? "false" : fPhiBgModulation == true ? "true" : "Weird") << endl;
cout << "Number of skipped hardest jets: " << fJetNHardestSkipped_010 << " for 0-10\%" << " and " << fJetNHardestSkipped_1080 << " for 10-80\%" << endl; 
cout << (fSetJetFixedSeed  == true ? Form("Fastjet with fixed seed: %d",fJetFJSeed) : fSetJetFixedSeed  == false ? "Fastjet with random seed" : "Weird") << endl;
  return kStOK;
}
//
// Function:  write to output file and close
//________________________________________________________________________
Int_t StHIOverlayAngularities::Finish()
{
  cout << "StHIOverlayAngularities::Finish()\n";

  //  Write  to file and close it.
  if (mOutName != "")
  {
    TFile *fout = new TFile(mOutName.Data(), "UPDATE");
    fout->cd();
    
    TList *infoList = outputTree->GetUserInfo();
    ////////-------------------------------------------------------------------------------------------------------------
    TString paramStr;
    
    TString runListStr;    
    runListStr = (fRunBadlist == 0) ? "Hanseul's bad run 2014 list\n" :
	         (fRunBadlist == 1) ? "Neil's bad run 2014 list\n" :
		                     "No bad run list used!!!\n";    
            
    TString mcRandomSeed;
    mcRandomSeed = 	(fSetMcSeed  == 0) ? "MC random seed\n" : 
    			(fSetMcSeed  == 1) ? "MC seed = eventId + runId\n" :
    			"Weird";
             
                                  
paramStr.Form(
    "----------------------------\n"
    "%s"
    "----------------------------\n"
    "Run cuts\n"
    "%s"
    "----------------------------\n"
    "Event cuts\n"
    "V_z min =< %.2f cm\n"
    "V_z max =< %.2f cm\n"
    "V_r =< %.2f cm\n"
    "|V_z-V_{z,VPD}| =< %.2f cm\n"
    "V_x != 0 && V_y != 0 && V_z != 0\n"
    "Triggers: {",
    mcRandomSeed.Data(),
    runListStr.Data(),
    fEventZVtxMinCut,
    fEventZVtxMaxCut,
    fEventRVtxMaxCut,
    fEventVtxVpdVzMaxCut
);

// Přidání seznamu triggerů
bool first = true;
for (const Int_t &trigger : fEventTriggers) {
    if (!first) paramStr += ", ";
    paramStr += TString::Format("%d", trigger);
    first = false;
}
paramStr += "}\n----------------------------\n";

paramStr += TString::Format("Jet constituent cuts\n"
              "********************\n"
              "  D0 meson\n"
              "MC particle mass: %.6f GeV/c^2\n"
              "pT_min => %.2f GeV/c\n"
              "|eta| =< %s\n"
              "********************\n",
              fMcD0Mass, fMinMcPtD0, (fSetMcAbsEtaD0) ? TString::Format("%.2f", fMcAbsEtaD0).Data() : "INF");

// Charged tracks
paramStr += TString::Format("  Charged tracks\n"
              "%s\n"
              "pT_min => %.2f GeV/c\n"
              "pT_max =< %.2f GeV/c\n"
              "eta_min => %.2f\n"
              "eta_max =< %.2f\n"
              "nHitsFit => %d\n"
              "nHitsRatio => %.2f\n"
              "Particle mass: %.6f GeV/c^2\n"
              "MC particle mass: %.6f GeV/c^2\n"
              "DCA =< %.2f cm\n"
              "********************\n",
              doUsePrimTracks ? "Primary tracks" : "Global tracks",
              fMinJetTrackPt, fMaxJetTrackPt, fJetTrackEtaMin, fJetTrackEtaMax, 
              fJetTracknHitsFit, fJetTracknHitsRatio, fChargedPart, fMcChargedPart, fJetTrackDCAcut);

// Towers
TString EventListStrA;
EventListStrA = (fSetTowerCalibrEnergy  == true ? "Hanseul's energy calibration\n" : fSetTowerCalibrEnergy  == false ? "Production energy calibration\n" : "Weird\n");

TString EventListStr;
EventListStr = 	(fTowerBadlist  == 0 ? "Hanseul's bad tower 2014 list\n" : 
		fTowerBadlist  == 1 ? "Neil's bad tower 2014 list\n" : 
		"None bad tower listed used!!!\n");

paramStr += TString::Format("  Towers\n"
              "ET => %.2f GeV\n"
              "%s"
              "%s"
              "Particle mass: %.6f GeV/c^2\n"
              "MC particle mass: %.6f GeV/c^2\n"
              "********************\n",
              mTowerEnergyTMin, EventListStr.Data(), EventListStrA.Data(), fNeutralPart, fMcNeutralPart);

TString paramHadr;

paramHadr = (fSetDcaZHadronCorr  == true ? "DCA_Z\n" : fSetDcaZHadronCorr  == false ? "DCA\n" : "Weird\n");

// Hadronic correction
paramStr += TString::Format("  Hadronic correction\n"
              "Fraction of hadr. corr. = %.2f\n"
              "%s\n"
              "pT_min => %.2f GeV/c\n"
              "pT_max =< %.2f GeV/c\n"
              "eta_min => %.2f\n"
              "eta_max =< %.2f\n"
              "nHitsFit => %d\n"
              "nHitsRatio => %.2f\n"
              "Particle mass: %.6f GeV/c^2\n"
              //"DCA_z =< %.2f cm\n"
              "%s =< %.2f cm\n"
              "----------------------------\n",
              mHadronicCorrFrac,
              doUsePrimTracks ? "Primary tracks" : "Global tracks",
              fMinJetTrackPt, fMaxJetTrackPt, fJetTrackEtaMin, fJetTrackEtaMax,
              fJetTracknHitsFit, fJetTracknHitsRatio, fHadronicCorrMass, 
              paramHadr.Data(), fDcaZHadronCorr);

// Jets
paramStr += TString::Format("  Jets\n"
              "**********\n"
              "Jet type: %s\n"
              "Jet rec. algorithm: %s\n"
              "Jet R = %.2f\n"
              "Jet area min >= %.2f\n"
              "|eta_jet| <= %.2f\n"
              "fJetFractionNeutralToJetPt: %.2f\n"
              "**********\n",
              fJetType == 0 ? "Full jets" : fJetType == 1 ? "Charged jets" : fJetType == 2 ? "Neutral jets" : "?",
              antikt_algorithm == 0 ? "kt" : antikt_algorithm == 1 ? "anti-kt" : "Other",
              fJetRad, fJetMinAreaCut, 1.0 - fJetRad, fJetFractionNeutralToJetPt);

// Background estimation
//Parama
	switch (fBgSubtraction) {
	    case 0:
		paramStr +=  "Area based method\n";
		break;

	    case 1:
		paramStr +=  "ICS method\n";
		break;

	    case 2:
		paramStr +=  "Jet shape method";
		break;

	    default:
		std::cerr << "Warning: Unknown fBgSubtraction option " << fBgSubtraction << std::endl;
		exit(1);
	}
paramStr += TString::Format("  Background estimation\n"
              "Massive particles? %s\n"
              "alpha_1: %.3f alpha_2: %.3f\n",
              fMassiveAll ? "true" : "false", fgAlpha1, fgAlpha2);
 /*             
              paramStr += TString::Format("  Background estimation\n"
              "Massive particles? %s\n"
              "alpha_1: %.3f\n",
              fMassiveAll ? "true" : "false", fgAlpha1);
   */
 paramStr += "0_0: 4it, R = 0.05; 1_0: 3it, R = 0.125; 2_0: 2it, R = 0.1-0.175; 3_0: 2it, R = 0.15-0.1 \n";       
 paramStr += "Default one: 2it, 0.15, 0.2 \n";                 
 paramStr += TString::Format("Background phi modulation? %s\n",
              fPhiBgModulation ? "true" : "false");      
 paramStr += TString::Format("Number of skipped hardest jets: %d for 0-10%% and %d for 10-80%%\n", fJetNHardestSkipped_010,fJetNHardestSkipped_1080); 
 
 	switch (fSetJetFixedSeed) {
	    case 0:
		paramStr +=  "Fastjet with random seed";
		break;

	    case 1:
		paramStr +=  TString::Format("Fastjet with fixed seed: %d",fJetFJSeed);
		break;
	}            
   
              


    ////////-------------------------------------------------------------------------------------------------------------
    infoList->Add(new TObjString(paramStr));
    
    TDirectory* dirEvent = fout->mkdir("event");
    dirEvent->cd();
    hVtxZ->Write();
    hVtxR->Write();
    hVzDiff->Write();
    hCentrality->Write();
    hCentralityW->Write();
    fout->cd();

    //Mc event:
    	TDirectory* dirMcEvent = fout->mkdir("mcEvent");
    	dirMcEvent->cd();
    	hPureMcNeutralEtaPhi->Write();
    	hMcJetConstMom->Write();
    	hMcJetConstTheta->Write();
	//hMcJetConstMom5GeV->Write();
	fout->cd();

    //Jet constituents:
	TDirectory* dirJetConstituents = fout->mkdir("jetConstituents");
	dirJetConstituents->cd();
	hJetTracksDedx->Write();
	hJetTracksDedxAfterCuts->Write();
	hJetTracksPt->Write();
	hJetTracksEtaPhi->Write();
	hJetTracksNHitsFit->Write();
	hJetTracksNHitsRatio->Write();
	hJetTracksDCA->Write();
        hJetMcRecoTracksPt->Write();
	hJetMcRecoTracksEtaPhi->Write();
	hJetMcRecoTracksNHitsFit->Write();
	hJetMcRecoTracksNHitsRatio->Write();
	hJetMcRecoTracksDCA->Write();
	hJetNeutralPt->Write();
	hJetNeutralEtaPhi->Write();
	hJetNeutralEtBefAftHC->Write();
	hJetNeutralECalibBefAft->Write();
	hJetMcRecoNeutralPt->Write();
	hJetMcRecoNeutralEtaPhi->Write();
	hJetMcRecoNeutralEtBefAftHC->Write();
	hJetConstCharge->Write();
	hJetConstEtaPhi->Write();
	hJetConstEtaPhiICS->Write();  
	hJetConstRapPhi->Write();
	hJetConstRapPhiICS->Write();   
	hJetConstPt->Write();
	hJetConstPtICS->Write();  
	hJetMcRecoConstCharge->Write();
	hFractionNeutralToJetPt->Write();
	effJetRec00_10->Write();
	effJetRec10_40->Write();
	effJetRec40_80->Write();
	etaRecoTrue00_10->Write();
	etaRecoTrue10_40->Write();
	etaRecoTrue40_80->Write();
	fout->cd();

   //Neutral particles hadronic correction:
   TDirectory* dirHadronCorr = fout->mkdir("hadronCorr");
   dirHadronCorr->cd();
   hJetHadrCorrEtaVsPt->Write();
   hJetHadrCorrE->Write();
   hJetHadrCorrNHitsFit->Write();
   hJetHadrCorrNHitsRatio->Write();
   hJetHadrCorrDcaZ->Write();
   fout->cd();
   
   
   TDirectory* dirResponseMatrix = fout->mkdir("responseMatrix");
   dirResponseMatrix->cd();
   hResponseJetPt->Write();
   hResponseJetD0Z->Write();
   hResponseJetNConst->Write();
   hResponseJetLambda1_0_5->Write();
   hResponseJetLambda1_1->Write();
   hResponseJetLambda1_1_5->Write();
   hResponseJetLambda1_2->Write();
   hResponseJetLambda1_3->Write();
   hResponseJetMomDisp->Write();
   hResponseJetD0DeltaR->Write();
   hResponseJetD0Pt->Write();
   fout->cd();

    outputTree->Write();

/*
    // write all histograms
    hRecoJetPt->Write();

    hMcJetPt->Write();
    hMcRecoJetPt->Write();
    hMcD0Pt->Write();
    hMcRecoD0Pt->Write();
    hMcJetZ->Write();
    hRecoJetZ->Write();
    hMcD0PtMcRecoD0Pt->Write();
    hRecoJetPtMcJetPt->Write();
    hRecoJetZMcJetZ->Write();*/
   fout->cd();

    // fout->Write();
    fout->Close();
    
    delete infoList;
  }
  

  

  cout << "End of StHIOverlayAngularities::Finish" << endl;
  StMemStat::PrintMem("End of Finish...");

  return kStOK;
}

// OLD user code says: //  Called every event after Make().
//_____________________________________________________________________________
void StHIOverlayAngularities::Clear(Option_t *opt)
{
}
//
//  Function: This method is called every event.
//_____________________________________________________________________________

//--------Event-plane---------------
Int_t StHIOverlayAngularities::EP_IsGoodTrack(StPicoTrack *trk, TVector3 pVertex){

	enum Position {West = -1, East = 1};
	Position Side = West;
	
	TVector3 mTrkMom;
	
	//Primary tracks
   	mTrkMom = trk->pMom(); //If not exists, it equals 0, 0, 0

	//0.2 < pt < 2 GeV
	float pt = mTrkMom.Perp();
	if (pt!=pt||pt <= 0.2 || pt >= 2.0) return 0;

	//0.05 < |eta| < 1.00 (West/East)
	float eta = mTrkMom.PseudoRapidity();
	if (abs(eta) <= 0.05 || abs(eta) >= 1) return 0;
	if (eta > 0) Side = East;

	//nHitsFit > 15
	float nHitsFit = trk->nHitsFit();
	float nHitsMax = trk->nHitsMax();
	
	//nHitsFit/nHitsMax => 0.52
	float nHitsRatio = 1.0*nHitsFit/nHitsMax;
	if (nHitsFit < 15 || nHitsRatio < 0.52) return 0;
	
	
	//DCA < 1 cm
	float dca = trk->gDCA(pVertex).Mag();
	if (dca >= 1) return 0;

return Side;
}


void StHIOverlayAngularities::CalculateEventPlane(){

    //Primary vertex
    TVector3 prVertex = mPicoEvent->primaryVertex();
    //Q-vectors
    double Q_1 = 0;
    double Q_2 = 0;
    
 
    

    for (UInt_t iTrack = 0; iTrack < mPicoDst->numberOfTracks(); iTrack++){
      
        StPicoTrack *trk = static_cast<StPicoTrack *>(mPicoDst->track(iTrack));
        if (!trk) continue;
        
        double Goodtrack = EP_IsGoodTrack(trk,prVertex);
        if (!abs(Goodtrack)) continue;
        
        TVector3 trackMom = trk->pMom();
        
        double pPt = trackMom.Perp();
        double phi = trackMom.Phi();
        
        //Q-vectors calculating
        Q_1 += 1.*pPt*cos(2*phi);
        Q_2 += 1.*pPt*sin(2*phi);
    }
    
   // cout << "Q_1: " << Q_1 << " Q_2: " << Q_2 << endl;
    
    fQ_1 = Q_1;
    fQ_2 = Q_2;
    
    //Recentering
    //-------------------
    //You have to run the code for "all" events to get the mean value of Q vector
        //7.2.2025
    double Q_1rc = -0.5976;
    double Q_2rc = 1.842;
    
    double Q_1corr = Q_1 - Q_1rc;
    double Q_2corr = Q_2 - Q_2rc;
    
    fQ_1_rec = Q_1corr;
    fQ_2_rec = Q_2corr;
    //-------------------
    //7.2.2025

    double Psi_2 = 1./2 * TMath::ATan2(Q_2corr, Q_1corr);
   
    //-------------------
    //Psi shift	
    //You have to run the code for "all" events to get the mean value of Q vector (again)
    std::vector<double> A_2 = {-0.0471131, -0.0114311, -0.000157011, 0.00128572, -0.000253291, 0.000446538, 0.00223142, 0.000574021, 0.00133264, 0.00157648, 0.00169576, -0.00023281, 0.00174849, 0.00194597, 0.000343813, 0.00104523, 0.00194369, -0.000649944, 0.00073013, -0.00138147, -0.000686723};
    std::vector<double> B_2 = {0.0204715, -0.0151115, -0.00306008, -0.000966112, 0.000103265, -0.0014271, -0.00244014, 0.000801171, -0.00235772, 0.00173085, -8.08274e-05, -0.000880916, 0.000967381, 0.00182509, 0.00167076, -6.81227e-05, 8.50044e-05, -0.000525437, -0.000912995, 0.000685826, -0.000467597};

    double CorrectedPsi2 = Psi_2;
    for (int i = 1; i <= 21; i++){
         CorrectedPsi2 +=(1.0 / 2) * (2.0 / i) *(-A_2[i - 1] * cos(2 * i * Psi_2) + B_2[i - 1] * sin(2 * i * Psi_2));
    }
    
    if (Psi_2!=Psi_2 || CorrectedPsi2!=CorrectedPsi2) return;
    
    fPsi_2 = Psi_2;
    
    //Force the range (-pi/2,pi/2)
    CorrectedPsi2 = TMath::ATan2(TMath::Sin(2 * CorrectedPsi2), TMath::Cos(2 * CorrectedPsi2)) / 2.;
    Psi_2 = TMath::ATan2(TMath::Sin(2 * Psi_2), TMath::Cos(2 * Psi_2)) / 2.;

    
    //-------------------
    fPsi_2_shifted = CorrectedPsi2;
   



return;
}

Int_t StHIOverlayAngularities::Make()
{



  // ZERO these out for Double_t checking they aren't set
  for (Int_t i = 0; i < 4800; i++){
  
    for (Int_t j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    mTowerStatusArr[i] = 0;

  }

  // get PicoDstMaker
  mPicoDstMaker = static_cast<StPicoDstMaker *>(GetMaker("picoDst"));
  if (!mPicoDstMaker)
  {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }
  // construct PicoDst object from maker
  mPicoDst = static_cast<StPicoDst *>(mPicoDstMaker->picoDst());
  if (!mPicoDst)
  {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  // create pointer to PicoEvent
  mPicoEvent = static_cast<StPicoEvent *>(mPicoDst->event());
  if (!mPicoEvent)
  {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }
  
  fRunNumber = mPicoEvent->runId();
  
  fMcSeed = 0; //0 = random
  if(fSetMcSeed) fMcSeed = fRunNumber + mPicoEvent->eventId();
  
  
  // get event B (magnetic) field
  Bfield = mPicoEvent->bField();
  
  // get vertex 3-vector and z-vertex component
  mVertex = mPicoEvent->primaryVertex();
  zVtx = mVertex.z();

  Double_t zVtx_VPD = mPicoEvent->vzVpd();

  // ============================ CENTRALITY ============================== //

  //Loadig of the GRefMultCorr information
  mGRefMultCorrUtil->init(fRunNumber);
  mGRefMultCorrUtil->initEvent(mPicoDst->event()->grefMult(),zVtx,mPicoDst->event()->ZDCx());

  //Check if the event is in a bad run (due to centrality estimation)
  if (mGRefMultCorrUtil->isBadRun(fRunNumber)) return kStOK;
  
  //Loading of the centrality
  Int_t centrality  = mGRefMultCorrUtil->getCentralityBin9();
   
  //Check if the centrality is in the range  // cut on unset centrality, > 80%
  if(centrality<0) return kStOK;
  
  //Loading of the weight for the centrality correction
  fCentralityWeight = mGRefMultCorrUtil->getWeight();
  
  // centrality variables
  //Double_t refMult = mGRefMultCorrUtil->getRefMultCorr(); //krefCorr2
  fKgrefMult_uncorr = mPicoEvent->grefMult(); // grefMult
  fKrefMult_uncorr = mPicoEvent->refMult();   // refMult 
  

  if (centrality == 7 || centrality == 8) // 0-10%
    centralitybinforefficiency = 0;
  else if (centrality == 4 || centrality == 5 || centrality == 6 ) // 10-40%
    centralitybinforefficiency = 1;
  else if (centrality == 0 || centrality == 1 || centrality == 2 || centrality == 3 ) // 40-80%
    centralitybinforefficiency = 2;

  fCentrality = centrality;

  fCentralityAlt =        fCentrality==8? 2.5 //0-5%
  			: fCentrality==7? 7.5 //5-10%
  			: fCentrality==6? 15  //10-20%
  			: fCentrality==5? 25  //20-30%
  			: fCentrality==4? 35  //30-40%
  			: fCentrality==3? 45  //40-50%
  			: fCentrality==2? 55  //50-60%
  			: fCentrality==1? 65  //60-70%
  			: fCentrality==0? 75  //70-80%
  			: -999;
  			
  // ============================ end of CENTRALITY ============================== //


//CHANGE
/*
  Bool_t isBadRun = false;  
  if (AllBadRunList2014.count(fRunNumber)) {
      isBadRun = true;

  }
  */
  
/*
  Bool_t isBadRun = true;  
  if (goodRun2014.count(fRunNumber)) {
      isBadRun = false;

  }
  */
  //if (isBadRun) return kStOK;  
  
  //2014 good run list
  //if (AllBadRunList2014.count(fRunNumber)) return kStOK;
  
  Bool_t isBadRun = true;
  
  //2014 good run list 
  if (true) {
  	if (fRunBadlist == 0) //Hanseul's
        	isBadRun = mycuts::AllBadRunList2014.count(fRunNumber);
	else if (fRunBadlist == 1) //Neil's
		isBadRun = !mycuts::goodRun2014.count(fRunNumber);
  }
   
  if (true && isBadRun) return kStOK;


  //Function checks if the event is triggered by the MB trigger
  ////set<Int_t> triggers = {450005, 450015, 450025, 450050, 450060}; // Run14_AuAu200
  //// {520001, 520011, 520021, 520031, 520041, 520051, 570002, 570001}; // Run16
 
  Bool_t matchMB = kFALSE;
  for (const Int_t &trigger : fEventTriggers){
  
    if (mPicoEvent->isTrigger(trigger)) matchMB = kTRUE;
      
    if (matchMB) break;
      
  }

  if (!matchMB) return kStOk;
  
  //Check on suspicious all-0 position
  //if (abs(mVertex.x()) < 1.0e-5 || abs(mVertex.y()) < 1.0e-5 || abs(mVertex.z()) < 1.0e-5) return kStOK;
  if (not (abs(mVertex.x()) !=0 && abs(mVertex.y()) !=0 && abs(mVertex.z()) !=0 )) return kStOK;
  
  //Checking vz
  if ((zVtx < fEventZVtxMinCut) || (zVtx > fEventZVtxMaxCut)) return kStOk;
  
  
  //Checking vr = sqrt(vx^2+vy^2)
  if (TMath::Sqrt(pow(mVertex.x(), 2) + pow(mVertex.y(), 2)) > fEventRVtxMaxCut) return kStOk;
  
  //Checking vzVpdVz
  if (abs(zVtx - zVtx_VPD) > fEventVtxVpdVzMaxCut) return kStOk;
  
  //Filling events histograms
  hVtxZ->Fill(zVtx);
  hVtxR->Fill(mVertex.x(),mVertex.y());
  hVzDiff->Fill(zVtx - zVtx_VPD);
  hCentrality->Fill(centrality);
  hCentralityW->Fill(centrality,fCentralityWeight);
  
  /////////////////
  if(fPhiBgModulation) CalculateEventPlane();
  
  // ============================ end of VERTEX+TRIGGER ============================== //

  ReadTreeMc();                                // load fMCPico + tree variables
  Int_t nMcEvents = fMCPico->GetEntriesFast(); // get number of events


 
  Int_t counterEvent = 0;
  Int_t eventsconsidered = 0; // different with counterEvent due to multiple events with D0s in them
  vector<Int_t> randomEventList;

  //Clear arrays
  for (Int_t iD0Event = 0; iD0Event < kMaxNumberOfD0Events; iD0Event++){
    fRecoMcEventTracks[iD0Event].clear();
    fRecoMcEventTowers[iD0Event].clear();
    fMcEventTracks[iD0Event].clear();
    fMcEventTowers[iD0Event].clear();
    fMcD0Information[iD0Event] = {};
    fMcRecoD0Information[iD0Event] = {};
    fOrigin[iD0Event].SetXYZ(0, 0, 0);
    fMcEventInfo[iD0Event] = {};
  }
  
    TRandom3 ran(fMcSeed);
    ran.SetSeed(fMcSeed); 

  //How many times is one event used
  while (counterEvent < fNumberOfEventsToOverLay){
  
    Int_t randomMcEvent = ran.Integer(nMcEvents);
    
    if (std::find(randomEventList.begin(), randomEventList.end(), randomMcEvent) != randomEventList.end()) continue;
      
    randomEventList.push_back(randomMcEvent);
    
    if (fPrintLevel)
    {
      cout << randomEventList.back() << endl;
      cout << "Size of randomEventList = " << randomEventList.size() << endl;
    }
    
    if (randomEventList.size() > static_cast<size_t>(nMcEvents)) continue;
    

    fMCPico->GetEntry(randomMcEvent); // load this MC event into memory
    
    if (fPrintLevel)
      cout << "======================= Event Number After MC = " << Event_mEventId[0] << "\t" << counterEvent << "=======================" << endl;

    //========================================================================================//
    //================================== Find D0 ids in MC part ==============================//
    //========================================================================================//
    vertexids.clear();
    pionids.clear();
    kaonids.clear();
    
    matchedpionids.clear();
    matchedkaonids.clear();
 //cout << "McTrack_: " << McTrack_<<endl;
    //Loop over all McTracks
    for (Int_t iMcTrack = 0; iMcTrack < McTrack_; iMcTrack++) {
   
      //37 D0 //38 antiD0
      if (McTrack_mGePid[iMcTrack] != 37 && McTrack_mGePid[iMcTrack] != 38) continue;
        
      //Save the D0  
      TLorentzVector particle;
      particle.SetPxPyPzE(McTrack_mPx[iMcTrack], McTrack_mPy[iMcTrack], McTrack_mPz[iMcTrack], McTrack_mE[iMcTrack]);
      
      //Track variables
      Double_t pt = particle.Pt();
      
      // Only D0s we care about need to have pT > 1 GeV
      if (pt < fMinMcPtD0) continue; 
      
      //Phi of D0
      Double_t phi = particle.Phi();
      while (phi < 0.0)
        phi += 2.0 * pi; // force from 0-2pi
      while (phi > 2.0 * pi)
        phi -= 2.0 * pi; // force from 0-2pi
        
      //Eta of D0  
      Double_t eta = particle.PseudoRapidity();

      //Eta D0 cut //not used if fMcAbsEtaD0 == -1
      if (fSetMcAbsEtaD0 && (abs(eta) > fMcAbsEtaD0)) continue;
        
      if (fPrintLevel) cout << "MC D0 Found = " << pt << "\t" << eta << "\t" << phi << endl;
      
      Bool_t isGoodDecay = kTRUE;
      Int_t pionTrack = -99;
      Int_t kaonTrack = -99;

      for (Int_t iMcDaughterTrack = 0; iMcDaughterTrack < McTrack_; iMcDaughterTrack++){
      
	      //We only want the kaon and pion that originated from the D0 we are interested in.
	      if (McTrack_mIdVtxStart[iMcDaughterTrack] != McTrack_mIdVtxStop[iMcTrack]) continue; 

              if (McTrack_mGePid[iMcDaughterTrack] != 8 && McTrack_mGePid[iMcDaughterTrack] != 9 && McTrack_mGePid[iMcDaughterTrack] != 11 && McTrack_mGePid[iMcDaughterTrack] != 12){
	          isGoodDecay = kFALSE;
		  break; // This should never happen
	      }
		
	      // Push the pion and kaon ids into the vector
              if (McTrack_mGePid[iMcDaughterTrack] == 8 || McTrack_mGePid[iMcDaughterTrack] == 9){
	      	pionTrack = iMcDaughterTrack;	
	      } else if (McTrack_mGePid[iMcDaughterTrack] == 11 || McTrack_mGePid[iMcDaughterTrack] == 12){
		kaonTrack = iMcDaughterTrack;
	      }
      }

      if (isGoodDecay){
		if (pionTrack == -99 || kaonTrack == -99){
		  cout << "Something wrong with D0. Exiting." << endl;
		  continue;
		}

		vertexids.push_back(McTrack_mIdVtxStop[iMcTrack]);
		pionids.push_back(pionTrack);
		kaonids.push_back(kaonTrack);
		matchedpionids.push_back(GetMatchedRecoTrackFromMCTrack(pionTrack));
		matchedkaonids.push_back(GetMatchedRecoTrackFromMCTrack(kaonTrack));
      
      }
    
    } //End of loop over all McTracks
    
    assert((pionids.size() == kaonids.size()) && "Same number of kaons and pions \n");
    
    if (fPrintLevel) cout << "=============== Number of D0s =============== " << vertexids.size() << "\t" << pionids.size() << endl;

    if (vertexids.size() == 0) continue; // While loop continues.

    //========================================================================================//
    //========================================================================================//

    // fill vertices map to loop over further when using GetAllTracksFromVertex()

    //Clear vertex position
    fVertexToTracks.clear();

    for (Int_t iMcTrack = 0; iMcTrack < McTrack_; iMcTrack++){
    
      Int_t idVertexStart = McTrack_mIdVtxStart[iMcTrack] - 1;
      fVertexToTracks[idVertexStart].push_back(iMcTrack);
    
    }

    //=========================================================================================//
    //================================== Loop over D0 found in event ==========================//
    //=========================================================================================//
    for (UInt_t iD0 = 0; iD0 < vertexids.size(); iD0++){
    //cerr << iD0+1 << "/" << vertexids.size() << endl;
    
    
      if (fPrintLevel)
      {
        cout << "==================================MC part==============================" << endl;
        cout << "=======================================================================" << endl;
      }
      // I only want to use TLorentzVector to propagate the information ahead. All the processing is done within this class.
      // This means I need to find a way to make sure I can identify the D0 track when I save out the information.
      // Since only D0 needs to be identified, I am saving the mass information for the D0 as 1.865.
      // All the other tracks are saved out with pi0mass, because ultimately, jet constituents are assumed to have that mass.

      // How do I propagate charge information though?
      // Do I need it? Mostly, nope! In fact, once I separate track and tower, it should be enough.

      // This loop fills the input vector for the MC Side

      // This function prepares the input list for that event for a particular D0.
      // The event list will of course be a little different for each D0, even for the same event

      fDroppedMCTracks.clear();

      Int_t daughterVertexStop1 = McTrack_mIdVtxStop[pionids[iD0]] - 1;
      Int_t daughterVertexStop2 = McTrack_mIdVtxStop[kaonids[iD0]] - 1;

      //Get all associated tracks
      GetAllTracksFromVertex(daughterVertexStop1, fDroppedMCTracks);
      GetAllTracksFromVertex(daughterVertexStop2, fDroppedMCTracks);

      if (fPrintLevel)
      {
        cout << "Dropped Track List for D0 # = " << iD0 << endl;
        for (size_t i = 0; i < fDroppedMCTracks.size(); i++)
        {
          TVector3 particle(McTrack_mPx[fDroppedMCTracks[i]], McTrack_mPy[fDroppedMCTracks[i]], McTrack_mPz[fDroppedMCTracks[i]]);
          cout << fDroppedMCTracks[i] << "\t" << McTrack_mGePid[fDroppedMCTracks[i]] << "\t" << particle.Perp() << endl;
        }
        cout << endl;
      }

      //=========================================================================================//
      //================ Loop for preparing MC tracks to fastjet wrapper ========================//
      //=========================================================================================//

      //Loop over all mctracks
      for (Int_t iMcTrack = 0; iMcTrack < McTrack_; iMcTrack++){
      
        TLorentzVector particle;
        
        //37 D0 //38 antiD0
        if ((McTrack_mGePid[iMcTrack] == 37 || McTrack_mGePid[iMcTrack] == 38)) particle.SetXYZM(McTrack_mPx[iMcTrack], McTrack_mPy[iMcTrack], McTrack_mPz[iMcTrack], 1.86482);
        //Charged particles
        else if (McTrack_mCharge[iMcTrack] != 0) particle.SetXYZM(McTrack_mPx[iMcTrack], McTrack_mPy[iMcTrack], McTrack_mPz[iMcTrack], fMcChargedPart);
        //Neutral particles
        else particle.SetXYZM(McTrack_mPx[iMcTrack], McTrack_mPy[iMcTrack], McTrack_mPz[iMcTrack], fMcNeutralPart);

        // track variables
        Double_t pt = particle.Pt();
        Double_t phi = particle.Phi();
        
        if (phi < 0.0) phi += 2.0 * pi; // force from 0-2pi
        if (phi > 2.0 * pi)phi -= 2.0 * pi; // force from 0-2pi
        
        Double_t eta = particle.PseudoRapidity();

        if (std::find(fDroppedMCTracks.begin(), fDroppedMCTracks.end(), iMcTrack) != fDroppedMCTracks.end())
        {
          if (fPrintLevel)
            cout << "Dropped Track = " << McTrack_mGePid[iMcTrack] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
          continue; // Discard tracks which come from kaon decay after D0 decay
        }
        
	// Only Final State Particles are included in JetMaker. D0s are the only exception
        if (McTrack_mIdVtxStop[iMcTrack] != 0 && (McTrack_mGePid[iMcTrack] != 37 && McTrack_mGePid[iMcTrack] != 38)) continue; 
          

        if (McTrack_mGePid[iMcTrack] == 37 || McTrack_mGePid[iMcTrack] == 38){
        
        if (McTrack_mIdVtxStop[iMcTrack] != vertexids[iD0]) continue; // Only consider the current D0
           
        }
        
        
        // Unstable particles which shouldn't make it to the end are discarded by hand. The list provisionally includes:
        /*
          Lambda, Eta, Sigma0, Xi0, Muon, Neutrino, KS0, KL0
        */
        if (McTrack_mGePid[iMcTrack] == 4 || McTrack_mGePid[iMcTrack] == 5 ||
            McTrack_mGePid[iMcTrack] == 6 || McTrack_mGePid[iMcTrack] == 10 ||
            McTrack_mGePid[iMcTrack] == 16 || McTrack_mGePid[iMcTrack] == 17 ||
            McTrack_mGePid[iMcTrack] == 18 || McTrack_mGePid[iMcTrack] == 20 ||
            McTrack_mGePid[iMcTrack] == 22)
          continue;

        // Have to discard Kaons and Pions which come from the Current D0

        if (McTrack_mIdVtxStart[iMcTrack] == vertexids[iD0]) continue;

 //Deleted (Check later!)
        //if ((pt < fMinJetTrackPt) || (pt > fMaxJetTrackPt) ||  (phi < fJetTrackPhiMin) || (phi > fJetTrackPhiMax)) continue;
        if ((pt < fMinJetTrackPt) || (pt > fMaxJetTrackPt)) continue; //It doesn't make sense for true level?
         
        //eta cut only for non-D0 particles
        if (((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax)) && (McTrack_mGePid[iMcTrack] != 37 && McTrack_mGePid[iMcTrack] != 38)) continue; 
       
        /* 
         cout << "Type of McTrack_mCharge: " << typeid(McTrack_mCharge[iMcTrack]).name() << endl;
         cout << "McTrack_mCharge: " << McTrack_mCharge[iMcTrack] << endl;
         cout << "McTrack_mCharge: " << static_cast<int>(McTrack_mCharge[iMcTrack]) << endl;
*/
        //Charged tracks + D0
        if (McTrack_mCharge[iMcTrack] != 0 || McTrack_mGePid[iMcTrack] == 37 || McTrack_mGePid[iMcTrack] == 38){
        
          fMcEventTracks[counterEvent].push_back(particle); 
          if (fPrintLevel == 2)
            cout << "Track Constituents = " << McTrack_mGePid[iMcTrack] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
        }//Neutral particles
        else{
        
          fMcEventTowers[counterEvent].push_back(particle);
          if (fPrintLevel == 2)
            cout << "Towers = " << McTrack_mGePid[iMcTrack] << "\t" << pt << "\t" << eta << "\t" << phi << endl;
        }
      
      } //End of loop over all mctracks
      
      //=========================================================================================//
      //=========================================================================================//
      
      ////cerr << "MC" << " " << iD0 << endl;
      fjw->Clear();
      fjw->SetCentrality(fCentrality);
      fjw->SetCentralityW(fCentralityWeight);
      fjw->SetBackgroundSub(kFALSE);
      fjw->SetSubtractionMc(kTRUE);
      fjw->SetMassiveTest(fMassiveAll);
      fjw->SetAlpha1(fgAlpha1);
      fjw->SetAlpha2(fgAlpha2);
      fjw->SetEventPlane2(fPsi_2_shifted);
      fjw->SetPhiModulation(fPhiBgModulation);
      fjw->setJetNHardestSkipped(fJetNHardestSkipped_010, fJetNHardestSkipped_1080);
      fjw->setJetFixedSeed(fSetJetFixedSeed, fJetFJSeed);
      
      // Loop over all saved particles for MC in vectors and enter them into fastjet wrapper
      for (UInt_t iTrack = 0; iTrack < fMcEventTracks[counterEvent].size(); iTrack++)
      {
        TLorentzVector v = fMcEventTracks[counterEvent][iTrack];
        // track variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + fMcChargedPart * fMcChargedPart);
        
        Double_t d0Factor = 1;
        if (int(v.M())==1) {
        	energy = 1.0 * TMath::Sqrt(p * p + fMcD0Mass * fMcD0Mass);
        	d0Factor = 3;
        }
        
        fjw->AddInputVector(px, py, pz, energy, iTrack + 10000*d0Factor); // includes E
      }                                                          // track loop
	
      if ((fJetType == kFullJet) || (fJetType == kNeutralJet)){
    
    	      //Tower loop
	      for (UInt_t iTower = 0; iTower < fMcEventTowers[counterEvent].size(); iTower++){
	      
		TLorentzVector v = fMcEventTowers[counterEvent][iTower];
		// tower variables
		Double_t px = v.X();
		Double_t py = v.Y();
		Double_t pz = v.Z();
		Double_t p = v.P();
		Double_t energy = 1.0 * TMath::Sqrt(p * p + fMcNeutralPart*fMcNeutralPart); 
		fjw->AddInputVector(px, py, pz, energy, iTower + 20000); // includes E
		
	      } //End of the tower loop        
      }                          
                          
      // Jet Running
      if (fPrintLevel == 2) fjw->PrintInput();

      // This function makes the jets with the MC tracks.
      // If it doesn't find a D0 jet with pt > 5 GeV in the event,
      // we discard the whole event and go to the next case.

      // Make the jets here. If they are not what we want them to be, discard the whole event and keep going.

      StJet *jetMc = DoesItHaveAGoodD0Jet(fMcEventTracks[counterEvent], 0);
      

      if (jetMc == NULL)
      {
        fMcEventTracks[counterEvent].clear();
        fMcEventTowers[counterEvent].clear();
        eventsconsidered++;
        continue;
      }
      
      jetMc->GetMomentumOfConst(hMcJetConstMom);
      jetMc->GetThetaOfConst(hMcJetConstTheta);
      //jetMc->GetMomentumOfConst(hMcJetConstMom5GeV,5);
      fMcJet[counterEvent] = jetMc;
      

      //=========================================================================================//
      //===================== preparing Reco MC tracks to fastjet wrapper =======================//
      //=========================================================================================//
      fjw->Clear();
      fjw->SetBackgroundSub(kFALSE);
      fjw->SetSubtractionMc(kTRUE);
      fjw->SetMassiveTest(fMassiveAll);
      fjw->SetCentrality(fCentrality);
      fjw->SetCentralityW(fCentralityWeight);
      fjw->SetAlpha1(fgAlpha1);
      fjw->SetAlpha2(fgAlpha2);
      fjw->SetEventPlane2(fPsi_2_shifted);
      fjw->SetPhiModulation(fPhiBgModulation);
      fjw->setJetNHardestSkipped(fJetNHardestSkipped_010, fJetNHardestSkipped_1080);
      fjw->setJetFixedSeed(fSetJetFixedSeed, fJetFJSeed);
      ////cerr << "Reco MC" << " " << iD0 << endl;
      if (fPrintLevel)
      {
        cout << "=================================RECO part=============================" << endl;
        cout << "=======================================================================" << endl;
      }
      PrepareSetOfRecoInput(counterEvent, iD0);
      // Loop over all saved particles for MC in StHIOverlayAngularities vectors and enter them into fastjet wrapper
      for (UInt_t iTracks = 0; iTracks < fRecoMcEventTracks[counterEvent].size(); iTracks++)
      {
        TLorentzVector v = fRecoMcEventTracks[counterEvent][iTracks];
        // track variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + fMcChargedPart * fMcChargedPart);
        Double_t d0Factor = 1;
        
        if (int(v.M())==1) {
        	energy = 1.0 * TMath::Sqrt(p * p + fMcD0Mass * fMcD0Mass);
        	d0Factor = 3;
        	}
        fjw->AddInputVector(px, py, pz, energy, iTracks + 10000*d0Factor); //  10k-20k is Mc tracks range
      
        //hJetMcRecoTracksPt->Fill(v.Perp(),fCentralityWeight);
        //hJetMcRecoTracksEtaPhi->Fill(v.Phi(), v.Eta(), fCentralityWeight);
        //hJetMcTracksNHitsFit->Fill(trk->nHitsFit(), fCentralityWeight);
	//hJetMcTracksNHitsRatio->Fill(1.0*trk->nHitsFit()/trk->nHitsMax(), fCentralityWeight);
	//hJetMcTracksDCA->Fill(fabs(trk->gDCA(mVertex.x(),mVertex.y(),mVertex.z())), fCentralityWeight);
      }                                                           // track loop

     if ((fJetType == kFullJet) || (fJetType == kNeutralJet))
    {
      for (UInt_t iTowers = 0; iTowers < fRecoMcEventTowers[counterEvent].size(); iTowers++)
      {
        TLorentzVector v = fRecoMcEventTowers[counterEvent][iTowers];
        // tower variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        //Double_t energy = 1.0 * TMath::Sqrt(p * p + M_PION_PLUS * M_PION_PLUS);
	Double_t energy = 1.0 * TMath::Sqrt(p * p + fMcNeutralPart*fMcNeutralPart);
        fjw->AddInputVector(px, py, pz, energy, iTowers + 20000); //  //  20k-30k is Mc tower range
      }
      }                                                           // tower loop

	////cerr << "jetMcreco" << endl;
      StJet *jetMcReco = DoesItHaveAGoodD0Jet(fRecoMcEventTracks[counterEvent], 1);
      fMcRecoJet[counterEvent] = jetMcReco;

      if (fPrintLevel == 2)
      {
        cout << "Event Number After Reco = " << Event_mEventId[0] << "\t" << counterEvent << endl;

        for (UInt_t iRecoTrack = 0; iRecoTrack < fRecoMcEventTracks[counterEvent].size(); iRecoTrack++)
        {
          TLorentzVector v = fRecoMcEventTracks[counterEvent][iRecoTrack];
          if ((int)v.M() != 1)
            cout << "Reco Constituents = ";
          else
            cout << "Reco D0 = ";
          cout << v.Pt() << "\t" << v.PseudoRapidity() << "\t" << v.Phi() << "\t" << v.M() << endl;
        }

        for (UInt_t iRecoTower = 0; iRecoTower < fRecoMcEventTowers[counterEvent].size(); iRecoTower++)
        {
          TLorentzVector v = fRecoMcEventTowers[counterEvent][iRecoTower];
          cout << "Reco Tower = " << v.Pt() << "\t" << v.PseudoRapidity() << "\t" << v.Phi() << "\t" << v.M() << endl;
        }
      }
      counterEvent++; // Since one MC event can be invoked multiple times (due to having multiple D0s, this step is necessary.)
      ////cerr << "Event započítán" << endl;
    }
    eventsconsidered++; // loop over D0s in event
    

  }                     // loop over while
  


  fMCPico->Reset();
  f->Close();
  //========================================================================================//
  //====================================end over MC part====================================//
  //========================================================================================//
  if (fPrintLevel)
  {
    cout << "================================ Real data + embed ====================" << endl;
    cout << "=======================================================================" << endl;
  }
  //=========================================================================================//
  //===================== Loop over Real data and embed PYTHIA tracks =======================//
  //=========================================================================================//

  for (Int_t iMcD0Event = 0; iMcD0Event < counterEvent; iMcD0Event++){
  
    fFull_Event.clear();
    
      // RESET tower-track matching pro hadron corr 
    for (Int_t i = 0; i < 4800; i++){
	    mTowerStatusArr[i] = 0;
	    for (Int_t j = 0; j < 7; j++) mTowerMatchTrkIndex[i][j] = -1;
    }

    if (fPrintLevel) cout << "Event Number  = " << iMcD0Event << endl;

    fjw->Clear();
    fjw->SetCentrality(fCentrality);
    fjw->SetCentralityW(fCentralityWeight);
    fjw->SetBackgroundSub(kTRUE);
    fjw->SetSubtractionMc(kTRUE);
    fjw->SetMassiveTest(fMassiveAll);
    fjw->SetAlpha1(fgAlpha1);
    fjw->SetAlpha2(fgAlpha2);
    fjw->SetEventPlane2(fPsi_2_shifted);
    fjw->SetPhiModulation(fPhiBgModulation);
    fjw->setJetNHardestSkipped(fJetNHardestSkipped_010, fJetNHardestSkipped_1080);
    fjw->setJetFixedSeed(fSetJetFixedSeed, fJetFJSeed);

    ////cerr << "Real data + embed " << " " << iMcD0Event << endl;
   	
    // loop over ALL tracks in PicoDst and add to jet, after acceptance and quality cuts
    if ((fJetType == kFullJet) || (fJetType == kChargedJet)){
    
      for (UInt_t iTrack = 0; iTrack < mPicoDst->numberOfTracks(); iTrack++){
      
        StPicoTrack *trk = static_cast<StPicoTrack *>(mPicoDst->track(iTrack));
        
        if (!trk) continue;

        hJetTracksDedx->Fill(trk->gPtot()*trk->charge(),trk->dEdx(),fCentralityWeight);
        // acceptance and kinematic quality cuts - pt cut is also applied here currently
        if (!IsAcceptedTrackAndPt(trk, Bfield, mVertex)) continue;

        // get momentum vector of track - global or primary track
        TVector3 mTrkMom;
        //cout << "doUsePrimTracks: " << doUsePrimTracks << endl;
        if (doUsePrimTracks)
        {
          // get primary track vector
          mTrkMom = trk->pMom();
        }
        else
        {
          // get global track vector
          mTrkMom = trk->gMom(mVertex, Bfield);
        }
        // track variables
        Double_t px = mTrkMom.x();
        Double_t py = mTrkMom.y();
        Double_t pz = mTrkMom.z();
        Double_t p = mTrkMom.Mag();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + fChargedPart * fChargedPart);

	//DCA < 3
	if(fabs(trk->gDCA(mVertex.x(),mVertex.y(),mVertex.z())) <= fJetTrackDCAcut && mTrkMom.Perp() <= fMaxJetTrackPt){
        //Filling jet track histogram before after cuts
        hJetTracksDedxAfterCuts->Fill(trk->gPtot()*trk->charge(),trk->dEdx(),fCentralityWeight);
        
        //cout << "input_particles.push_back({"<<px<<","<<py<<","<<pz<<","<<energy<<"}); //charged" << endl;
        fjw->AddInputVector(px, py, pz, energy, iTrack); // includes E
 
        fastjet::PseudoJet particle_Track(px, py, pz, energy);
        
          
        particle_Track.set_user_index(iTrack);
        fFull_Event.push_back(particle_Track);
        
        //QA histograms
        hJetTracksPt->Fill(particle_Track.perp(),fCentralityWeight);
        hJetTracksEtaPhi->Fill(particle_Track.phi_std(), particle_Track.eta(), fCentralityWeight);
        hJetTracksNHitsFit->Fill(trk->nHitsFit(), fCentralityWeight);
	hJetTracksNHitsRatio->Fill(1.0*trk->nHitsFit()/trk->nHitsMax(), fCentralityWeight);
	hJetTracksDCA->Fill(fabs(trk->gDCA(mVertex.x(),mVertex.y(),mVertex.z())), fCentralityWeight);
	hJetConstCharge->Fill(trk->charge(), fCentralityWeight);
	hJetConstRapPhi->Fill(particle_Track.phi_std(),particle_Track.rap(), fCentralityWeight);
	hJetConstEtaPhi->Fill(particle_Track.phi_std(),particle_Track.eta(), fCentralityWeight);
	hJetConstPt->Fill(particle_Track.perp(),fCentralityWeight);

        }// end of DCA < 3
        
        //Old
        //Int_t matchedTowerIndex = abs(GetMatchedBtowID(trk)) - 1; // towerIndex = towerID - 1

        Int_t matchedTowerIndex = trk->bemcTowerIndex(); //returns ID

        if (matchedTowerIndex < 0) continue;

        //DCA_z < 3 cm
        Float_t dca_z = -999;
        if (fSetDcaZHadronCorr) dca_z = trk->gDCAz(mVertex.z());
        else dca_z = trk->gDCA(mVertex.x(),mVertex.y(),mVertex.z());
        
        if (fabs(dca_z) > fDcaZHadronCorr) continue;
        
        ////if(abs(trk->gDCA(mVertex).z()) > fHadrCorrTrackDCAZcut) continue;
        ////if(fabs(trk->gDCA(mVertex.x(),mVertex.y(),mVertex.z())) > fJetTrackDCAcut) continue;
        
        mTowerMatchTrkIndex[matchedTowerIndex][mTowerStatusArr[matchedTowerIndex]] = iTrack;
        mTowerStatusArr[matchedTowerIndex] = mTowerStatusArr[matchedTowerIndex] + 1; // 1+ means match, 0 for no match
      }
    }

    // full or neutral jets - get towers and apply hadronic correction
    if ((fJetType == kFullJet) || (fJetType == kNeutralJet))
    {

      fMaxTowerEtBeforeHC = -99;
      fMaxTowerEtAfterHC = -99;

      Int_t nTowers = mPicoDst->numberOfBTowHits(); // barrel tower hists, always 4800

      // loop over towers and add input vectors to fastjet
      for (Int_t itow = 0; itow < nTowers; itow++)
      {
        // get tower pointer
        StPicoBTowHit *tower = static_cast<StPicoBTowHit *>(mPicoDst->btowHit(itow));
        if (!tower)
        {
          cout << "No tower pointer... iTow = " << itow << endl;
          continue;
        }

        // tower ID - get from itow shift: which is 1 more than array index
        // itow = matchedTowerIndex
        Int_t towerID = itow + 1;
        Int_t towerIndex = towerID - 1;
        if (towerID < 0) continue; // Double_t check these aren't still in the event list

        //if (BadTowerMap[towerID-1]) continue; 
        if (fTowerBadlist == 0 && mycuts::BadTowerMap[towerID-1]) continue;
	if (fTowerBadlist == 1 && mycuts::NeilBadTowers2014.count(towerID)) continue;
 
        // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
        //TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
        StThreeVectorF vtx(mVertex.X(), mVertex.Y(), mVertex.Z());
        StThreeVectorF tmpTowerPosition = mEmcPosition->getPosFromVertex(vtx, towerID);
        TVector3 towerPosition(tmpTowerPosition.x(), tmpTowerPosition.y(), tmpTowerPosition.z());
        
        Double_t towerPhi = towerPosition.Phi();
        if (towerPhi < 0.0)
          towerPhi += 2.0 * pi; // force from 0-2pi
        if (towerPhi > 2.0 * pi)
          towerPhi -= 2.0 * pi; // force from 0-2pi
          
        Double_t towerEta = towerPosition.PseudoRapidity();

        Double_t towE = -999;
        //Calculate the tower energy
        if(fSetTowerCalibrEnergy) towE = GetTowerCalibEnergy(towerID);
        else towE = tower->energy(); 
        
        Double_t towEUncorr = tower->energy(); //energy w/o calibration
	Double_t towerEunCorr = towE; // uncorrected energy
        Double_t towerE = towE;       // corrected energy (hadronically - done below)
       // Double_t towEtunCorr = towerEunCorr / (1.0 * TMath::CosH(towerEta));

        // cut on min tower energy after filling histos
        if (towerEunCorr < mTowerEnergyTMin) continue; // if we don't have enough E to start with, why mess around

        if (towerEunCorr > fMaxTowerEtBeforeHC) fMaxTowerEtBeforeHC = towerEunCorr;

	//cout << "A: " << mTowerEnergyTMin << " " << fMaxTowerEtBeforeHC << endl;

        // =======================================================================
        // HADRONIC CORRECTION
        //Double_t maxEt = 0.;
        Double_t sumEt = 0.;

        // if tower was is matched to a track or multiple, add up the matched track energies - (mult opt.) to then subtract from the corresponding tower
        // August 15: if *have* 1+ matched trk-tow AND uncorrected energy of tower is at least your tower constituent cut, then CONTINUE
        if (mTowerStatusArr[towerIndex] > 0.5 && towerEunCorr > mTowerEnergyTMin)
        {
          Double_t maxE = 0.0;
          Double_t sumE = 0.0;
          // =======================================================================================================================
          // --- finds max E track matched to tower *AND* the sum of all matched track E and subtract from said tower
          //     USER provides readMacro.C which method to use for their analysis via SetJetHadCorrType(type);
          // loop over ALL matched tracks
          for (Int_t itrk = 0; itrk < mTowerStatusArr[towerIndex]; itrk++)
          {
            StPicoTrack *trk = static_cast<StPicoTrack *>(mPicoDst->track(mTowerMatchTrkIndex[towerIndex][itrk]));
            if (!trk)
            {
              cout << "No track pointer..." << endl;
              continue;
            }
            if (!IsAcceptedTrack(trk, Bfield, mVertex))
            {
              cout << "track matched back doesn't pass cuts" << endl;
              continue;
            }
            
            //DCA_z < 3 cm
            Float_t dca_z = -999;
            if (fSetDcaZHadronCorr) dca_z = trk->gDCAz(mVertex.z());
            else dca_z = trk->gDCA(mVertex.x(),mVertex.y(),mVertex.z());
        
            if (fabs(dca_z) > fDcaZHadronCorr)
            {
              cout << "track matched back doesn't pass cut dca/dca_z" << endl;
              continue;
            }



            // get track variables to matched tower from 3-vector
            TVector3 mTrkMom;
            if (doUsePrimTracks)
            { // get primary track vector
              mTrkMom = trk->pMom();
            }
            else
            { // get global track vector
              mTrkMom = trk->gMom(mVertex, Bfield);
            }

            double pT = mTrkMom.Perp();
            if (pT < fMinJetTrackPt) continue;


            // track variables
            Double_t p = mTrkMom.Mag();
            Double_t E = 1.0 * TMath::Sqrt(p * p + fHadronicCorrMass * fHadronicCorrMass);


	    hJetHadrCorrNHitsFit->Fill(trk->nHitsFit(), fCentralityWeight);
	    hJetHadrCorrNHitsRatio->Fill(1.*trk->nHitsFit()/trk->nHitsMax(), fCentralityWeight);
	    hJetHadrCorrDcaZ->Fill(abs(trk->gDCA(mVertex).z()), fCentralityWeight);
	    hJetHadrCorrEtaVsPt->Fill(pT, trk->gMom().Eta(), fCentralityWeight);
	    hJetHadrCorrE->Fill(E, fCentralityWeight);

            if (E > maxE) maxE = E;
            sumE = sumE + E;

       } // track loop

                // apply hadronic correction to tower
         // maxEt = (towerEunCorr - (mHadronicCorrFrac * maxE)) / (1.0 * TMath::CosH(towerEta));
          sumEt = (towerEunCorr - (mHadronicCorrFrac * sumE)) / (1.0 * TMath::CosH(towerEta));
          //=================================================================================================================
        } // have a track-tower match
        // else - no match so treat towers on their own. Must meet constituent cut

        // Et - correction comparison
        //Double_t fMaxEt = (maxEt == 0) ? towerEunCorr / (1.0 * TMath::CosH(towerEta)) : maxEt;
        Double_t fSumEt = (sumEt == 0) ? towerEunCorr / (1.0 * TMath::CosH(towerEta)) : sumEt;
        // if(mTowerStatusArr[towerIndex] > 0.5) cout<<"towerEunCorr = "<<towerEunCorr<<"  CosH: "<<1.0*TMath::CosH(towerEta)<<"   fMaxEt: "<<fMaxEt<<"   fSumEt: "<<fSumEt<<endl;

        // cut on transverse tower energy (more uniform)
        Double_t towerEt = 0.0;
        if (mTowerStatusArr[towerIndex] < 1){
           // no matches, use towers uncorrected energy
           towerEt = towerEunCorr / (1.0 * TMath::CosH(towerEta));
        }
        else{
        
	   towerEt = fSumEt;
	   towerE = fSumEt * 1.0 * TMath::CosH(towerEta);
       
        }

        if (towerEt < 0) towerEt = 0.0;

        if (towerEt < mTowerEnergyTMin) continue;


        if (towerEt > fMaxTowerEtAfterHC) fMaxTowerEtAfterHC = towerEt;


        // get components from Energy (p - momentum) - the below lines 'converts' the tower energy to momentum:
        TVector3 mom;
        GetMomentum(mom, tower, fNeutralPart, towerID, towerE); //gammas
        Double_t towerPx = mom.x();
        Double_t towerPy = mom.y();
        Double_t towerPz = mom.z();

        // add towers to fastjet - shift tower index (tracks 0+, ghosts = -1, towers < -1)
        Int_t uidTow = -(itow + 2);

        fjw->AddInputVector(towerPx, towerPy, towerPz, towerE, uidTow); // includes E
        ////cout << "input_particles.push_back({"<<towerPx<<","<<towerPy<<","<<towerPz<<","<<towerE<<"}); //neutral" << endl;
        fastjet::PseudoJet particle_Track(towerPx, towerPy, towerPz, towerE);
        particle_Track.set_user_index(uidTow);
        fFull_Event.push_back(particle_Track);
        
        hJetNeutralPt->Fill(particle_Track.perp(), fCentralityWeight);
	hJetNeutralEtaPhi->Fill(particle_Track.phi_std(), particle_Track.eta(), fCentralityWeight);
	hJetNeutralEtBefAftHC->Fill(towerEunCorr/ (1.0 * TMath::CosH(towerEta)), towerEt, fCentralityWeight);
	hJetNeutralECalibBefAft->Fill(towEUncorr, towE, fCentralityWeight);
	hJetConstCharge->Fill(0., fCentralityWeight);
	hJetConstRapPhi->Fill(particle_Track.phi_std(),particle_Track.rap(), fCentralityWeight); 
	hJetConstEtaPhi->Fill(particle_Track.phi_std(),particle_Track.eta(), fCentralityWeight); 
	hJetConstPt->Fill(particle_Track.perp(),fCentralityWeight);
      } // tower loop

    } // neutral/full jets

    for (UInt_t iMcTrack = 0; iMcTrack < fRecoMcEventTracks[iMcD0Event].size(); iMcTrack++){
    

      TLorentzVector v = fRecoMcEventTracks[iMcD0Event][iMcTrack];
      // track variables
      Double_t px = v.X();
      Double_t py = v.Y();
      Double_t pz = v.Z();
      Double_t p = v.P();
      Double_t energy = 1.0 * TMath::Sqrt(p * p + fMcChargedPart * fMcChargedPart);
      
      Double_t d0Factor = 1.;
      if ((int)v.M() == 1) {
      	energy = 1.0 * TMath::Sqrt(p * p + fMcD0Mass * fMcD0Mass);
      	d0Factor = 3;
      	} //D0 mass
      
      if ((int)v.M() == 1 && fPrintLevel) cout << "Reco D0 In Jet Loop = " << v.Pt() << "\t" << v.PseudoRapidity() << endl;



      // MC Tracks will start from 10000
      fjw->AddInputVector(px, py, pz, energy, iMcTrack + 10000*d0Factor); // includes E
      
      fastjet::PseudoJet particle_Track(px, py, pz, energy);
      particle_Track.set_user_index(iMcTrack + 10000*d0Factor);
      ////if ((int)v.M() == 1) cout << "input_particles.push_back("<<px<<","<<py<<","<<pz<<","<<energy <<"); //D0 " << iMcTrack+10000 << endl;
      ////else cout << "input_particles.push_back("<<px<<","<<py<<","<<pz<<","<<energy <<"); //Charged " << iMcTrack+10000 << endl;
      
    } 

    if (fJetType == kFullJet || (fJetType == kNeutralJet))
    {
      for (UInt_t iTowers = 0; iTowers < fRecoMcEventTowers[iMcD0Event].size(); iTowers++)
      {
        TLorentzVector v;
        v = fRecoMcEventTowers[iMcD0Event][iTowers];
        // tower variables
        Double_t px = v.X();
        Double_t py = v.Y();
        Double_t pz = v.Z();
        Double_t p = v.P();
        Double_t energy = 1.0 * TMath::Sqrt(p * p + fMcNeutralPart*fMcNeutralPart); // Towers, gammas
        // MC Towers will start from
        fjw->AddInputVector(px, py, pz, energy, -1 * iTowers - 10000); // includes E
        fastjet::PseudoJet particle_Track(px, py, pz, energy);
        particle_Track.set_user_index(-1 * iTowers - 10000);

      } // tower loop

    } // if full/charged jets

int NiterA[4] = {4,3,2,2};
double R_max1A[4] = {0.05,0.125,0.100,0.150};
double R_max2A[4] = {0.005,0.005,0.175,0.100};

int i_1 = 2;

		fjw->SetNIter(NiterA[i_1]);
		fjw->SetCentrality(fCentrality);
		fjw->SetCentralityW(fCentralityWeight);
   		fjw->SetBackgroundSub(kTRUE);
   	        fjw->SetSubtractionMc(kTRUE);
   		fjw->SetMassiveTest(fMassiveAll);
   		fjw->SetEventPlane2(fPsi_2_shifted);
   		fjw->SetPhiModulation(fPhiBgModulation);
		fjw->SetAlpha1(fgAlpha1);
		fjw->SetAlpha2(fgAlpha2);
		fjw->SetRMax1(R_max1A[i_1]);
		fjw->SetRMax2(R_max2A[i_1]);
		fjw->setJetNHardestSkipped(fJetNHardestSkipped_010, fJetNHardestSkipped_1080);	
		fjw->setJetFixedSeed(fSetJetFixedSeed, fJetFJSeed);
	        
	       /// cerr << "jetreco" << endl;
    StJet *jetReco = DoesItHaveAGoodD0Jet(fRecoMcEventTracks[iMcD0Event], 2); // run fjw and get the jets

    
    fRecoJet[iMcD0Event] = jetReco;
    

  } // loop over D0 from MC

  FillTree(counterEvent);


  return kStOK;

} // Make

//
// Function: track quality cuts
//________________________________________________________________________
Bool_t StHIOverlayAngularities::IsAcceptedTrack(StPicoTrack *trk, const Float_t &B, const TVector3 &Vert)
{
  //Double_t pi = 1.0 * TMath::Pi();
  // get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  if (doUsePrimTracks)
  {
    if (!(trk->isPrimary()))
      return kFALSE; // check if primary
    // get primary track vector
    mTrkMom = trk->pMom();
  }
  else
  {
    // get global track vector
    mTrkMom = trk->gMom(Vert, B);
  }
  // track variables
  //Double_t pt = mTrkMom.Perp();
  //Double_t phi = mTrkMom.Phi();
  Double_t eta = mTrkMom.PseudoRapidity();
  //Double_t dca = trk->gDCA(Vert).Mag(); //DCA cut done on different place
  Int_t nHitsFit = trk->nHitsFit();
  Int_t nHitsMax = trk->nHitsMax();
  Double_t nHitsRatio = 1.0 * nHitsFit / nHitsMax;
  if ((eta < fJetTrackEtaMin) || (eta > fJetTrackEtaMax)) return kFALSE;

  if (nHitsFit < fJetTracknHitsFit) return kFALSE;

  if (nHitsRatio < fJetTracknHitsRatio) return kFALSE;

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;

}
//
// Function: jet track quality cuts
// - this function should only be used for jet constituent tracks
// 	NOT when considering track-tower matches  (Aug 19')

Bool_t StHIOverlayAngularities::IsAcceptedTrackAndPt(StPicoTrack *trk, const Float_t &B, const TVector3 &Vert)
{
  if (!IsAcceptedTrack(trk, B, Vert)) return kFALSE; // first do general track QA cuts

  // get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  
  if (doUsePrimTracks) mTrkMom = trk->pMom(); // get primary track vector
  else mTrkMom = trk->gMom(Vert, B); // get global track vector

  Double_t pt = mTrkMom.Perp();
  //if (pt < fMinJetTrackPt || pt > fMaxJetTrackPt) return kFALSE;
    if (pt < fMinJetTrackPt) return kFALSE; //upper cut done later
    
  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}
//
// Tower Quality Cuts
//________________________________________________________________________
Bool_t StHIOverlayAngularities::IsAcceptedTower(StPicoBTowHit *tower, const Int_t &towerID) //not used
{ // constants:
  Double_t pi = 1.0 * TMath::Pi();
  // tower ID - passed into function: make sure some of these aren't still in event array
  if (towerID < 0)
    return kFALSE;
    
  // cluster and tower position - from vertex and ID: shouldn't need additional eta correction
  StThreeVectorF vtx(mVertex.X(), mVertex.Y(), mVertex.Z());
  StThreeVectorF tmpTowerPosition = mEmcPosition->getPosFromVertex(vtx, towerID);
  TVector3 towerPosition(tmpTowerPosition.x(), tmpTowerPosition.y(), tmpTowerPosition.z());
  
  Double_t phi = towerPosition.Phi();
  if (phi < 0.0)
    phi += 2.0 * pi; // force from 0-2pi
  if (phi > 2.0 * pi)
    phi -= 2.0 * pi; // force from 0-2pi
  Double_t eta = towerPosition.PseudoRapidity();
  // check for bad (and dead) towers
  /*
  Bool_t TowerOK = mBaseMaker->IsTowerOK(towerID);     // kTRUE means GOOD
  Bool_t TowerDead = mBaseMaker->IsTowerDead(towerID); // kTRUE means BAD
  if (!TowerOK || TowerDead)
  {
    return kFALSE;
  }
*/
  // jet track acceptance cuts njow
  if ((eta < fJetTowerEtaMin) || (eta > fJetTowerEtaMax))
    return kFALSE;
  if ((phi < fJetTowerPhiMin) || (phi > fJetTowerPhiMax))
    return kFALSE;
  // passed all above cuts - keep tower and fill input vector to fastjet
  return kTRUE;
}

void StHIOverlayAngularities::PrepareSetOfRecoInput(const Int_t &counterEvent, const Int_t &iD0){

/*
  fRecoMcEventTracks[counterEvent].clear();
  fRecoMcEventTowers[counterEvent].clear();
  */

  const Int_t numberoftowers = BTowHit_;
  Double_t towerenergy[numberoftowers];
  
  //Reset all tower energy
  for (Int_t i = 0; i < numberoftowers; i++) towerenergy[i] = 0.;

  TVector3 oVertex;
  oVertex.SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

  if (fPrintLevel == 2)
    cout << "HIOverlay Vertex " << Event_mPrimaryVertexX[0] << "\t" << Event_mPrimaryVertexY[0] << "\t" << Event_mPrimaryVertexZ[0] << endl;

  // This loop fills the input vector for the RECO side
  for (Int_t reco = 0; reco < Track_; reco++){
  

    TVector3 o;
    o.SetXYZ(Track_mOriginX[reco], Track_mOriginY[reco], Track_mOriginZ[reco]);
    TVector3 g;
    g.SetXYZ(Track_mGMomentumX[reco], Track_mGMomentumY[reco], Track_mGMomentumZ[reco]);

    // track variables
    Double_t pt = g.Perp();
    Double_t phi = g.Phi();
    if (phi < 0.0) phi += 2.0 * pi; // force from 0-2pi
    if (phi > 2.0 * pi) phi -= 2.0 * pi; // force from 0-2pi
    Double_t eta = g.PseudoRapidity();
    Double_t px = g.x();
    Double_t py = g.y();
    Double_t pz = g.z();
    Double_t p = g.Mag();
    Int_t mcid = Track_mIdTruth[reco] - 1;
//    Double_t energy = 1.0 * TMath::Sqrt(p * p + mass * mass);
    Double_t energy_hadr_corr = 1.0 * TMath::Sqrt(p * p + fHadronicCorrMass * fHadronicCorrMass);
    //Double_t energy = McTrack_mE[mcid];
    short charge = (Track_mNHitsFit[reco] > 0) ? 1 : -1;
    Double_t dca = (oVertex - o).Mag();

    Bool_t goodtrack = (dca < fJetTrackDCAcut) &&
    		       (abs(Track_mNHitsFit[reco]) >= fTracknHitsFit) && 
    		       (1.*abs(Double_t(Track_mNHitsFit[reco])) / Double_t(Track_mNHitsMax[reco]) >= fTracknHitsRatio);

	//cout << "fJetTrackDCAcut: " << fJetTrackDCAcut << endl;
	//cout << "fTracknHitsFit: " << fTracknHitsFit << endl;	
	//cout << "fTracknHitsRatio: " << fTracknHitsRatio << endl;
	//cout << "test: " << abs(Double_t(Track_mNHitsFit[reco])) / Double_t(Track_mNHitsMax[reco]) << endl;

    //// Variables For FastSim

    Double_t pt_new = pt;
    Double_t phi_new = phi;
    Double_t eta_new = eta;
    Double_t px_new = px;
    Double_t py_new = py;
    Double_t pz_new = pz;
   // Double_t p_new = p;
    //Double_t energy_new = energy;

    //short charge_new = charge;

    Bool_t mctrackavailable = kTRUE;
    //Int_t mcid = Track_mIdTruth[reco] - 1;
    if (mcid < 0) mctrackavailable = kFALSE;

    Bool_t isatrackfromD0 = kFALSE;

    // Here, we have two paths to take. If the track needs replacement, we replace it with the fastsim method that is standardised.
    // Else, the pt, eta, phi are sent as is to the final vector.

    if (mctrackavailable && goodtrack){
    
      // Kaons and Pions that come from the current D0 need to be tossed, and replaced by the fast sim version
      if (reco == matchedpionids[iD0] || reco == matchedkaonids[iD0]) isatrackfromD0 = kTRUE; 


      TVector3 mg(McTrack_mPx[mcid], McTrack_mPy[mcid], McTrack_mPz[mcid]);
      Double_t relativesmearing = TMath::Sqrt(pow(mg.Px() - px, 2) + pow(mg.Py() - py, 2)) / (mg.Pt());
      Double_t fastsimsmearing;
      Int_t pid = McTrack_mGePid[mcid];

      if (pid == 8 || pid == 9)
        fastsimsmearing = fPionMomResolution->Eval(mg.Pt()); // Pion
      else if (pid == 11 || pid == 12)
        fastsimsmearing = fKaonMomResolution->Eval(mg.Pt()); // Kaon
      else if (pid == 15 || pid == 14)
        fastsimsmearing = fProtonMomResolution->Eval(mg.Pt()); // Proton
      else
        fastsimsmearing = fPionMomResolution->Eval(mg.Pt()); // Catch all: pions

	////cout << sqrt(McTrack_mE[mcid]*McTrack_mE[mcid] - mg.Mag()*mg.Mag()) << " ate " << endl;

	
      if (relativesmearing > 3 * fastsimsmearing){
      
        TVector3 fastsimsmearedmom = FastSimMom(mg, pid);
        pt_new = fastsimsmearedmom.Perp();
        phi_new = fastsimsmearedmom.Phi();
        if (phi_new < 0.0) phi_new += 2.0 * pi; // force from 0-2pi
        if (phi_new > 2.0 * pi) phi_new -= 2.0 * pi; // force from 0-2pi

        eta_new = fastsimsmearedmom.PseudoRapidity();
        px_new = fastsimsmearedmom.x();
        py_new = fastsimsmearedmom.y();
        pz_new = fastsimsmearedmom.z();
        //p_new = fastsimsmearedmom.Mag();
      }
      
    }
    if (fPrintLevel == 2)
    {
    /*
      if (mctrackavailable)
        cout << Form("Track # %i \t %i \t %i \t %.2f \t %.2f \t %i \t %i \t %.2f \t %i \t %.2f", reco, mcid, McTrack_mGePid[mcid], pt_new, eta_new, isatrackfromD0, abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1, dca, abs(int(Track_mNHitsFit[reco])), abs(Double_t(Track_mNHitsFit[reco])) / Double_t(Track_mNHitsMax[reco])) << endl;
      else
        cout << Form("Track # %i \t %.2f \t %.2f \t %i \t %.2f \t %.2f", reco, pt_new, eta_new, abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1, dca, abs(Double_t(Track_mNHitsFit[reco])) / Double_t(Track_mNHitsMax[reco]), Track_mQATruth[reco]) << endl;
        */
    }

    if (!goodtrack) continue;
    // DCA based cuts precede everything else.

    // Track from D0 -> K Pi || D0 is in acceptance range. The KPi do not need to be in acceptance. || The KPi track is projected onto the towers.
    if (isatrackfromD0){
    
      Int_t matchedTowerIndex = abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1;
      if (matchedTowerIndex >= 0) towerenergy[matchedTowerIndex] += energy_hadr_corr;
     
    }

    Int_t particleid = -99;

    Double_t nsigpion = Track_mNSigmaPion[reco] / 1000.;
    Double_t nsigkaon = Track_mNSigmaKaon[reco] / 1000.;
    Double_t nsigproton = Track_mNSigmaProton[reco] / 1000.;

	//Co to je?
    if (abs(nsigpion) < 2 && abs(nsigkaon) > 2. && abs(nsigproton) > 2.)
      particleid = 1;
    else if (abs(nsigpion) > 2 && abs(nsigkaon) < 2. && abs(nsigproton) > 2.)
      particleid = 2;
    else if (abs(nsigpion) > 2 && abs(nsigkaon) > 2. && abs(nsigproton) < 2.)
      particleid = 3 * charge;

    Bool_t removetrack = kFALSE;
    Bool_t isD0DaugDescendant = kFALSE;
    if (!isatrackfromD0 &&
       pt > fMinJetTrackPt && pt < fMaxJetTrackPt && 
       (eta > fJetTrackEtaMin) && (eta < fJetTrackEtaMax)){
    // I am using the old pt for this because the efficiencies were derived with the old pt.
    // I am keeping things consistent with MC
    
      // KeepTrack returns False if the track is to be discarded // RemoveTrack = True if KeepTrack is False
      removetrack = !KeepTrack(particleid, centralitybinforefficiency, pt);
      // If D0 descendant, record that as well.
      isD0DaugDescendant = kFALSE;
      Int_t mctrkid = Track_mIdTruth[reco] - 1;
      if (std::find(fDroppedMCTracks.begin(), fDroppedMCTracks.end(), mctrkid) != fDroppedMCTracks.end())
        isD0DaugDescendant = kTRUE;
    }

    if (removetrack || isD0DaugDescendant)
    {
      if (fPrintLevel == 2)
      {
        if (mctrackavailable)
          cout << Form("Removed Track Input = %.2f \t %.2f \t %.2f \t %i", px_new, py_new, pz_new, McTrack_mGePid[mcid]) << endl;
        else
          cout << Form("Removed Track Input = %.2f \t %.2f \t %.2f", px_new, py_new, pz_new) << endl;
      }
    }

    // jet track acceptance cuts now
    if (pt_new < fMinJetTrackPt || pt_new > fMaxJetTrackPt) continue; // 20.0 STAR, 100.0 ALICE

    if ((eta_new < fJetTrackEtaMin) || (eta_new > fJetTrackEtaMax)) continue;

  /*  
    while (phi_new < 0.0)
      phi_new += 2.0 * pi; // force from 0-2pi
    while (phi_new > 2.0 * pi)
      phi_new -= 2.0 * pi; // force from 0-2pi
    if ((phi_new < fJetTrackPhiMin) || (phi_new > fJetTrackPhiMax))
      continue;
*/
    // additional quality cuts for tracks

    // This is questionable. Here, I don't think I should use it. But, in the QM method, I should?
    // I have now included these cuts in both. I think that's the correct thing to do.

    // This place takes care of all energy depositions due to tracks accepted for jet reco.
    // For the reco tracks that we discard, we still need to subtract their contribution from the tower
    // It also includes the kaon and pion from D0.

    Bool_t ignoretrack = removetrack || isatrackfromD0 || isD0DaugDescendant;
    if (ignoretrack) continue; // To match the efficiency, we start tossing random tracks.


    Int_t matchedTowerIndex = abs(int(Track_mBEmcMatchedTowerIndex[reco])) - 1;
    
    if (matchedTowerIndex >= 0) towerenergy[matchedTowerIndex] += energy_hadr_corr;

    TLorentzVector v;
    //v.SetXYZM(px_new, py_new, pz_new, M_PION_PLUS);
    v.SetXYZM(px_new, py_new, pz_new, fMcChargedPart); 
    fRecoMcEventTracks[counterEvent].push_back(v);
    
    
    hJetMcRecoTracksPt->Fill(v.Perp(),fCentralityWeight);
    hJetMcRecoTracksEtaPhi->Fill(v.Phi(), v.Eta(), fCentralityWeight);
    hJetMcRecoTracksNHitsFit->Fill(abs(Track_mNHitsFit[reco]), fCentralityWeight);
    hJetMcRecoTracksNHitsRatio->Fill(1.*abs(Double_t(Track_mNHitsFit[reco])) / Double_t(Track_mNHitsMax[reco]), fCentralityWeight);
    hJetMcRecoTracksDCA->Fill(dca, fCentralityWeight);
    
    hJetConstRapPhi->Fill(v.Phi(), v.Rapidity(), fCentralityWeight); 
    hJetConstEtaPhi->Fill(v.Phi(), v.Eta(), fCentralityWeight); 
    hJetConstPt->Fill(v.Perp(),fCentralityWeight);
    
    hJetMcRecoConstCharge->Fill(charge, fCentralityWeight);
    
    
    //// Fill the Tower Array with Energy Depositions here
    if (fPrintLevel == 2)
    {
      if (mctrackavailable)
        cout << Form("Reco Level Track Input = %.2f \t %.2f \t %.2f \t %i", px_new, py_new, pz_new, McTrack_mGePid[mcid]) << endl;
      else
        cout << Form("Reco Level Track Input = %.2f \t %.2f \t %.2f", px_new, py_new, pz_new) << endl;
    }
  }

  TVector3 mcKaon(McTrack_mPx[kaonids[iD0]], McTrack_mPy[kaonids[iD0]], McTrack_mPz[kaonids[iD0]]);
  TVector3 mcPion(McTrack_mPx[pionids[iD0]], McTrack_mPy[pionids[iD0]], McTrack_mPz[pionids[iD0]]);

  fMcD0Information[counterEvent] = {mcPion, mcKaon};

  TVector3 recoKaon = FastSimMom(mcKaon, 11);
  TVector3 recoPion = FastSimMom(mcPion, 8);

  fMcRecoD0Information[counterEvent] = {recoPion, recoKaon};

  fOrigin[counterEvent].SetXYZ(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);

  TVector3 recoD0;
  recoD0 = recoKaon + recoPion;

  TLorentzVector v;
  v.SetXYZM(recoD0.x(), recoD0.y(), recoD0.z(), 1.86482); // This is a D0

  if (fPrintLevel == 2)
    cout << "Reco D0 Entered Jet Loop = " << v.Pt() << "\t" << v.PseudoRapidity() << endl;

  if (fPrintLevel == 2)
    cout << Form("Reco Level Track Input D0 = %.2f \t %.2f \t %.2f", recoD0.x(), recoD0.y(), recoD0.z()) << endl;

  fRecoMcEventTracks[counterEvent].push_back(v);

  //This loop fills the input vector for the towers
  for (Int_t tower = 0; tower < numberoftowers; tower++){
  
    Int_t towerID = tower + 1;
    if (towerID < 0) continue; // Double_t check these aren't still in the event list

   // TVector3 towerPosition = mEmcPosition->getPosFromVertex(oVertex, towerID);
    StThreeVectorF vtx(oVertex.X(), oVertex.Y(), oVertex.Z());
    StThreeVectorF tmpTowerPosition = mEmcPosition->getPosFromVertex(vtx, towerID);
    TVector3 towerPosition(tmpTowerPosition.x(), tmpTowerPosition.y(), tmpTowerPosition.z());
    
    Double_t towerPhi = towerPosition.Phi();
    if (towerPhi < 0.0)
      towerPhi += 2.0 * pi; // force from 0-2pi
    if (towerPhi > 2.0 * pi)
      towerPhi -= 2.0 * pi; // force from 0-2pi
    Double_t towerEta = towerPosition.PseudoRapidity();

    if ((Double_t(BTowHit_mE[tower]) / 1000. / TMath::CosH(towerEta)) > mTowerEnergyTMin) hPureMcNeutralEtaPhi->Fill(towerPosition.Phi(),towerEta,fCentralityWeight);

    ////if (BadTowerMap[towerID-1]) continue; //Ondra
    if (fTowerBadlist == 0 && mycuts::BadTowerMap[towerID-1]) continue;
    if (fTowerBadlist == 1 && mycuts::NeilBadTowers2014.count(towerID)) continue;
  
    // jet track acceptance cuts now
    ////if ((towerEta < fJetTowerEtaMin) || (towerEta > fJetTowerEtaMax)) continue;

    ////if ((towerPhi < fJetTowerPhiMin) || (towerPhi > fJetTowerPhiMax)) continue;
      

    Double_t towerEunCorr = Double_t(BTowHit_mE[tower]) / 1000.; // uncorrected energy
    Double_t towerE = Double_t(BTowHit_mE[tower]) / 1000.;       // corrected energy (hadronically - done below)
   // Double_t towEtunCorr = towerE / (1.0 * TMath::CosH(towerEta));

    //if (fPrintLevel == 2)
    //  cout << "Tower = " << tower << "\t" << towerE << "\t" << towEtunCorr << endl;

    // cut on min tower energy after filling histos
    if (towerEunCorr < mTowerEnergyTMin) continue; // if we don't have enough E to start with, why mess around

    // =======================================================================
    // HADRONIC CORRECTION

    Double_t sumEt = (towerEunCorr - mHadronicCorrFrac*towerenergy[tower]) / (1.0 * TMath::CosH(towerEta));
    Double_t towerEt = sumEt;

    if (towerEt < mTowerEnergyTMin) continue;
    
    towerE = towerEt * 1.0 * TMath::CosH(towerEta);

    if (fPrintLevel == 2) cout << "Tower = " << tower << "\t" << towerE << "\t" << sumEt << endl;

    //Double_t p = 1.0 * TMath::Sqrt(towerE * towerE - M_PION_PLUS * M_PION_PLUS); //Gamma
    Double_t p = 1.0 * TMath::Sqrt(towerE * towerE - fMcNeutralPart*fMcNeutralPart); //Ondra
    
    Double_t posX = towerPosition.x();
    Double_t posY = towerPosition.y();
    Double_t posZ = towerPosition.z();

    Double_t r = TMath::Sqrt(posX * posX + posY * posY + posZ * posZ);

    TLorentzVector v;
    v.SetXYZM(p * posX / r, p * posY / r, p * posZ / r, fMcNeutralPart);

    fRecoMcEventTowers[counterEvent].push_back(v);
    
    hJetMcRecoNeutralPt->Fill(v.Perp(), fCentralityWeight);
    hJetMcRecoNeutralEtaPhi->Fill(v.Phi(), v.Eta(), fCentralityWeight);
    hJetMcRecoNeutralEtBefAftHC->Fill(towerEunCorr/ (1.0 * TMath::CosH(towerEta)), towerEt, fCentralityWeight);
    
    hJetConstRapPhi->Fill(v.Phi(),v.Rapidity(), fCentralityWeight); 
    hJetConstEtaPhi->Fill(v.Phi(),v.Eta(), fCentralityWeight); 
    hJetConstPt->Fill(v.Perp(),fCentralityWeight);
    hJetMcRecoConstCharge->Fill(0., fCentralityWeight);
    
    if (fPrintLevel == 2)
      cout << Form("Reco Level Tower Input = %i \t %.2f \t %.2f \t %.2f \t %.2f", towerID, v.X(), v.Y(), v.Z(), towerE) << endl;
  }
}

//
// Function: calculate momentum of a tower
//________________________________________________________________________________________________________



Bool_t StHIOverlayAngularities::GetMomentum(TVector3 &mom, const StPicoBTowHit *tower, Double_t mass, Int_t towerID, Double_t CorrectedEnergy) const
{
  // get mass, E, P, ID
  if (mass < 0) mass = 0.0;
    
  Double_t energy = CorrectedEnergy; // USE THIS, to use the corrected energy to get the momentum components
  Double_t p = 1.0 * TMath::Sqrt(energy * energy - mass * mass);
  if (p!=p) p = 0;
  // get tower position
  //TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, towerID);
  StThreeVectorF vtx(mVertex.X(), mVertex.Y(), mVertex.Z());
  StThreeVectorF tmpTowerPosition = mEmcPosition->getPosFromVertex(vtx, towerID);
  TVector3 towerPosition(tmpTowerPosition.x(), tmpTowerPosition.y(), tmpTowerPosition.z());
  Double_t posX = towerPosition.x();
  Double_t posY = towerPosition.y();
  Double_t posZ = towerPosition.z();
  // get r, set position components
  Double_t r = TMath::Sqrt(posX * posX + posY * posY + posZ * posZ);
  if (r > 1e-12)
  {
    mom.SetX(p * posX / r);
    mom.SetY(p * posY / r);
    mom.SetZ(p * posZ / r);
  }
  else
  {
    return kFALSE;
  }
  return kTRUE;
}
/*
Int_t StHIOverlayAngularities::GetMatchedBtowID(StPicoTrack *trk)
{
  Double_t bemc_radius = mBemcGeom->Radius();
  // Magnetic field in Tesla
  Double_t mBField_tesla = Bfield / 10.0; // Check this definition. Magnetic fields are minefields of error in STAR
  // Needed for projection of the track onto the barrel radius
  TVector3 bemc_pos, bemc_mom;
  // BEMC hardware indices
  Int_t h_m, h_e, h_s = 0;
  // tower index: if no tower can be matched, assign 0
  // picoTrk->setBEmcMatchedTowerIndex(0);
  Int_t tow_id = 0;
  Bool_t close_match = false;
  Int_t trkbemcid = trk->bemcTowerIndex();
  // Check if the track can be projected onto the current radius
  // if not, track can't be matched.
  // By JetCorr request the global track projection to BEMC is used.
  if (mEmcPosition->projTrack(&bemc_pos, &bemc_mom, trk, mVertex, Bfield, bemc_radius))
  {
    // First, examine track eta. If it falls in two regions:
    // 0 < |eta| < etaMin()
    // etaMax() < |eta| < 1.0
    // then shift the eta for the projection slightly into the neighboring tower
   // TVector3 towerPosition = mEmcPosition->getPosFromVertex(mVertex, trkbemcid + 1);
     StThreeVectorF vtx(mVertex.X(), mVertex.Y(), mVertex.Z());
     StThreeVectorF tmpTowerPosition = mEmcPosition->getPosFromVertex(vtx,  trkbemcid + 1);
     TVector3 towerPosition(tmpTowerPosition.x(), tmpTowerPosition.y(), tmpTowerPosition.z());
     
    if (fabs(bemc_pos.PseudoRapidity()) < mBemcGeom->EtaMin())
    {
      Double_t unsigned_eta = mBemcGeom->EtaMin() + 0.001;
      // Double_t unsigned_eta = mBemcGeom->EtaMin() + 0.000001;
      Double_t unsigned_theta = 2.0 * atan(exp(-1.0 * unsigned_eta));
      Double_t signed_theta = (bemc_pos.PseudoRapidity() >= 0 ? 1.0 : -1.0) * unsigned_theta;
      bemc_pos.SetTheta(signed_theta);
      close_match = true;
    }
    else if (fabs(bemc_pos.PseudoRapidity()) > mBemcGeom->EtaMax() &&
             fabs(bemc_pos.PseudoRapidity()) < 1.0)
    {
      Double_t unsigned_eta = mBemcGeom->EtaMax() - 0.001;
      // Double_t unsigned_eta = mBemcGeom->EtaMax() - 0.000001;
      Double_t unsigned_theta = 2.0 * atan(exp(-1.0 * unsigned_eta));
      Double_t signed_theta = (bemc_pos.PseudoRapidity() >= 0 ? 1.0 : -1.0) * unsigned_theta;
      bemc_pos.SetTheta(signed_theta);
      close_match = true;
    }

    // Get the BEMC hardware location in (m, e, s) and translate to id
    // If StEmcGeom::getBin() != 0: track was not matched to a tower.
    // Its outside of the BEMC eta range (> 1.0).

    if (mBemcGeom->getBin(bemc_pos.Phi(), bemc_pos.PseudoRapidity(), h_m, h_e, h_s) == 0)
    {
      // If StEmcGeom::getId() == 0: the track was matched successfully. Otherwise,
      // the track was not matched to a tower at this radius, the track was projected
      // into the gap between modules in phi.
      if (h_s != -1)
      {
        mBemcGeom->getId(h_m, h_e, h_s, tow_id);
        if (close_match)
        {
          return -1 * tow_id;
        }

        else
        {
          return tow_id;
        }
      }

      // Track fell in between modules in phi. We will find which module it is closer
      // to by shifting phi slightly.
      else
      {
        // Value of the "dead space" per module in phi:
        // 2*pi/60 (amount of azimuth covered per module)
        // 2*0.0495324 (active size of module)

        Double_t dphi = (TMath::Pi() / 60.0) - 0.0495324;
        // Shift the projected phi by dphi in positive and negative directions
        // if we look for the projection for both of these, only one should give
        // a tower id, and the other should still be in the inter-tower space

        TVector3 bemc_pos_shift_pos(bemc_pos);
        bemc_pos_shift_pos.SetPhi(bemc_pos_shift_pos.Phi() + dphi);
        TVector3 bemc_pos_shift_neg(bemc_pos);
        bemc_pos_shift_neg.SetPhi(bemc_pos_shift_neg.Phi() - dphi);

        if (mBemcGeom->getBin(bemc_pos_shift_pos.Phi(), bemc_pos_shift_pos.PseudoRapidity(), h_m, h_e, h_s) == 0 && h_s != -1)
        {
          mBemcGeom->getId(h_m, h_e, h_s, tow_id);
          return -1 * tow_id;
        }

        else if (mBemcGeom->getBin(bemc_pos_shift_neg.Phi(), bemc_pos_shift_neg.PseudoRapidity(), h_m, h_e, h_s) == 0 && h_s != -1)
        {
          mBemcGeom->getId(h_m, h_e, h_s, tow_id);
          return -1 * tow_id;
        }
      }
    }
  }
  return tow_id;
}
*/
Bool_t StHIOverlayAngularities::KeepTrack(const Int_t &particleid, const Int_t &centralitybin, const Double_t &pt)
{
  Bool_t keeptrack = kTRUE;
  TRandom3 *r = new TRandom3(fMcSeed);
  r->SetSeed(fMcSeed);
  Double_t rando = r->Rndm();
  if (particleid == 1 || particleid == -1)
    keeptrack = (rando > fPionWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE; // Either a pion
  else if (particleid == 2 || particleid == -2)
    keeptrack = (rando > fKaonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == 3)
    keeptrack = (rando > fProtonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  else if (particleid == -3)
    keeptrack = (rando > fAProtonWeight[centralitybin]->Eval(pt)) ? kFALSE : kTRUE;
  if (keeptrack)
  {
    if (fTrackingEfficiency)
    {
      keeptrack = (rando > fTrackingEfficiencyPercentage) ? kFALSE : kTRUE;
    }
  }
  delete r;
  return keeptrack;
}


StJet *StHIOverlayAngularities::DoesItHaveAGoodD0Jet(vector<TLorentzVector> &eventTracks, Int_t recoLevel){
//recolevel: 0 = mc; 1 = mcSmeared; 2 = reco (sim + real event) 
//Check, if the D0-jet is presented


	// 0 - Area based, 1 - ICS, 2 - jet shape
	switch (fBgSubtraction) {
	    case 0:
		fjw->Run_AreaBase();
		break;

	    case 1:
		fjw->Run();
		break;

	    case 2:
		fjw->Run_Shape();
		break;

	    default:
		std::cerr << "Warning: Unknown fBgSubtraction option " << fBgSubtraction << std::endl;
		exit(1); // 1 bývá běžnější pro chybu
	}

  //Run the jet
 
  
  //Get all jets
  std::vector<fastjet::PseudoJet> jets_incl = fjw->GetInclusiveJets();
    
  //Sort jets according to jet pt
  Int_t nJets = jets_incl.size();
  Int_t indexes[9999];
  Float_t pt[9999] = {0};
  
  //Load pT of jets
  for (Int_t i = 0; i < nJets; i++) pt[i] = jets_incl[i].perp();
  
  //Sort of jets
  TMath::Sort(nJets, pt, indexes);
  
  //Initialization
  Int_t D0JetIndex = -99;
  
  //Loop over all jets
  for (Int_t iJetUnsorted = 0; iJetUnsorted < nJets; ++iJetUnsorted){
  
    //n-th sorted jet
    Int_t iJet = indexes[iJetUnsorted];
    
    //PERFORM CUTS ON inclusive JETS before saving
    
    //Cut on min jet pt
    if (fSetJetMinPtCut[recoLevel] && (jets_incl[iJet].perp() < fJetMinPtCut)) continue;
    
    //Cut on min jet area
    if (fSetJetMinAreaCut[recoLevel]  && (fjw->GetJetArea(iJet) < fJetMinAreaCut * TMath::Pi() * fJetRad * fJetRad)) continue;

    //Cut on eta acceptance
    if (fSetJetMinAbsEtaCut[recoLevel] && (abs(jets_incl[iJet].eta()) > fJetMinAbsEtaCut)) continue; 

    //Cut on phi acceptance
    ////if ((jets_incl[iJet].phi() < fJetPhiMin) || (jets_incl[iJet].phi() > fJetPhiMax)) continue;

    //Fill jet constituents
    vector<fastjet::PseudoJet> constituents = fjw->GetJetConstituents(iJet);

    //Initialisation of D0-jet check
    Bool_t D0Jet = kFALSE;

    //Loop over all constituents
    for (UInt_t iConstituent = 0; iConstituent < constituents.size(); ++iConstituent){
    
      // get user defined index
      Int_t uid = constituents[iConstituent].user_index();

      //This is a charged MC + D0 track id range
      if ( (int(uid) >= 10000 && int(uid) < 20000) || int(uid) >= 30000){ 
      
	//Load of the particle   
	Int_t partIndex = int(uid)<30000?(uid-10000):(uid-30000);
        TLorentzVector v = eventTracks[partIndex];
        
        //Check if the particle is D0 (other particles have mass at most pi+ mass)
        if ((int)v.M() == 1){
        
          //D0-jet flag
          D0Jet = kTRUE;
          
          //D0-jet index
          D0JetIndex = iJet;
          
          //Load D0-jet parameters
          Double_t D0Pt = v.Pt();
          Double_t D0Eta = v.PseudoRapidity();
          Double_t D0Phi = v.Phi();
          
          if (fPrintLevel) cout << "D0 Found with pT eta phi = " << D0Pt << "\t" << D0Eta << "\t" << D0Phi << endl;
          
          //We do not have to check the rest of the jet constituents
          break;
          
        } //End of check if the particle is D0
        
      } //End of charged MC + D0 track condition
      
    } //End of loop over all constituents

    //If D0 is not found, continue
    if (!D0Jet) continue;
    //Load jet pt
    Double_t jet_pt = jets_incl[iJet].perp();
    
    //Initialisation of neutral_pt fraction
    Double_t neutral_pt = 0.;
    
    //Loop over all jet constituents
    for (UInt_t iConstituent = 0; iConstituent < constituents.size(); ++iConstituent){
      
	    //Get user defined index
	    Int_t uid = constituents[iConstituent].user_index();
	    
	    //This is a reco tower
	    if ((int(uid) < -1 || int(uid) >= 20000) && int(uid) < 30000) neutral_pt += constituents[iConstituent].perp();

    } //End of loop over all jet constituents

    hFractionNeutralToJetPt->Fill(1.*neutral_pt / jet_pt);

    //Neutral (D0 not cosidered) pT fraction cut
    if (fSetJetFracCut[recoLevel] && ((neutral_pt / jet_pt) > fJetFractionNeutralToJetPt) ) continue;

    //Save the D0-jet
    StJet *jet = new StJet(jets_incl[D0JetIndex].perp(), jets_incl[D0JetIndex].eta(), jets_incl[D0JetIndex].phi(), jets_incl[D0JetIndex].m());
    
    if (fPrintLevel) cout << "Jet Found with pT eta phi = " << jet->Pt() << "\t" << jet->Eta() << "\t" << jet->Phi() << endl;

    //Save D0 jet index
    jet->SetLabel(D0JetIndex);
    
    //Save in addition the jet area information
    fastjet::PseudoJet area(fjw->GetJetAreaVector(D0JetIndex));
    //jet->SetArea(fjw->GetJetArea(D0JetIndex)); //Ondra
    jet->SetArea(area.perp());
  
    jet->SetJetShapeCorrL10half(fjw->GetActualShapeL10half());
    jet->SetJetShapeCorrL11(fjw->GetActualShapeL11());
    jet->SetJetShapeCorrL11half(fjw->GetActualShapeL11half());
    jet->SetJetShapeCorrL12(fjw->GetActualShapeL12());
    jet->SetJetShapeCorrL13(fjw->GetActualShapeL13());
    jet->SetJetShapeCorrDisp(fjw->GetActualShapeDisp());
     
    jet->SetJetRho(fjw->GetJetRho());
    jet->SetJetRhoM(fjw->GetJetRhoM());
     
    jet->SetArea4Vector(fjw->GetJetAreaVector(D0JetIndex));

    //jet->SetArea(area.perp()); // same as fjw->GetJetArea(iJet)
    jet->SetAreaEta(area.eta());
    jet->SetAreaPhi(area.phi());
    jet->SetAreaE(area.E());
    
    //Save D0-jet constituents
    constituents = fjw->GetJetConstituents(D0JetIndex);
    jet->SetJetConstituents(constituents);
    
    //All information saved. Now, time to clear it.
    fjw->Clear(); 
    
    //Return the saved D0-jet
    return jet;
    
  } //End of loop over all jets
  
  
  //If no D0-jet is found
  if (fPrintLevel) cout << "No D0 Jet Found" << endl;
  //Time to clear fjw
  fjw->Clear();//POZOR
  
  //Return null pointer
  return NULL;
  
} //End of function to check, if the D0-jet is presented

void StHIOverlayAngularities::GetAllTracksFromVertex(const Int_t &vertexid, vector<Int_t> &tracksFromVertex){

  //Check if the vertex ID is valid
  if (vertexid < 0) return;
    
  if (fPrintLevel == 2){
  
    cout << "Called this function for vx " << vertexid << " with ntracks = " << fVertexToTracks[vertexid].size() << endl;
    
    for (UInt_t track = 0; track < fVertexToTracks[vertexid].size(); track++) cout << "Geant ID of tracks = " << McTrack_mGePid[fVertexToTracks[vertexid][track]] << endl;

  }
  
  // Loop over all tracks associated with this vertex
  for (UInt_t track = 0; track < fVertexToTracks[vertexid].size(); track++){
  
    //Add the current track to the output vector
    tracksFromVertex.push_back(fVertexToTracks[vertexid][track]);
    
    //Get the ID of the vertex where this track ends (stop vertex)
    Int_t idvxstop = McTrack_mIdVtxStop[fVertexToTracks[vertexid][track]];
    
    //Recursively call the function for the stop vertex, if there is one
    GetAllTracksFromVertex(idvxstop - 1, tracksFromVertex);
    
  }//End of loop over all tracks associated with this vertex
  
  return;
}

Int_t StHIOverlayAngularities::GetMatchedRecoTrackFromMCTrack(const Int_t &McTrack){
//Only for daughter pions and kaons
  Int_t recotrackmatch = -999;

  //Load vertex position
  TVector3 eventVertex(Event_mPrimaryVertexX[0], Event_mPrimaryVertexY[0], Event_mPrimaryVertexZ[0]);
  Double_t ratio = -99.;
  
  
  //Loop over all reconstructed tracks
  for (Int_t recoTrack = 0; recoTrack < Track_; recoTrack++){
  
    //Get track pointer
    TVector3 originVertex(Track_mOriginX[recoTrack], Track_mOriginY[recoTrack], Track_mOriginZ[recoTrack]);

    //Calculate parameters
    Double_t dca = (eventVertex - originVertex).Mag();
    Int_t nHitsFit = abs(Track_mNHitsFit[recoTrack]);
    Int_t nHitsMax = abs(Track_mNHitsMax[recoTrack]);
    Double_t nHitsRatio = 1.0 * nHitsFit / nHitsMax;

    // additional quality cuts for tracks
    if (dca > fJetTrackDCAcut) continue;

    if (nHitsFit < fTracknHitsFit) continue;

    if (nHitsRatio < fTracknHitsRatio) continue;

    Int_t McTrackFromReco = Track_mIdTruth[recoTrack] - 1;

    //Check, if the track corresponds to the MC particle
    if (McTrackFromReco == McTrack){
    
      // Update the best matching track if the current track has a higher hit ratio.
      if (nHitsRatio > ratio){
        recotrackmatch = recoTrack;
        ratio = nHitsRatio;
      }
      
    }
    
  } //End of loop over all reconstructed tracks
  
  //The best reconstructed track
  return recotrackmatch;
  
}

TVector3 StHIOverlayAngularities::FastSimMom(TVector3 p, Int_t pid)
{
  float pt = p.Perp();
  float pt1 = pt;
  if (pt1 > 10)
    pt1 = 10; // Used for high pt-hat bin smearing test
  float sPt = -1;

  TRandom3 *r = new TRandom3(fMcSeed);
  r->SetSeed(fMcSeed);

  if (pid == 8 || pid == 9)
    sPt = r->Gaus(pt, pt * fPionMomResolution->Eval(pt1)); // Pion
  else if (pid == 11 || pid == 12)
    sPt = r->Gaus(pt, pt * fKaonMomResolution->Eval(pt1)); // Kaon
  else if (pid == 15 || pid == 14)
    sPt = r->Gaus(pt, pt * fProtonMomResolution->Eval(pt1)); // Proton
  else
    sPt = r->Gaus(pt, pt * fPionMomResolution->Eval(pt1)); // Catch all: pions

  TVector3 smearedmom(sPt * cos(p.Phi()), sPt * sin(p.Phi()), sPt * sinh(p.PseudoRapidity()));
  delete r;
  return smearedmom;
}

void StHIOverlayAngularities::ReadTreeMc()
{

  TRandom3 *r1 = new TRandom3(fMcSeed);
  r1->SetSeed(fMcSeed);
  Int_t filenumber = r1->Integer(filenamesforHIOverlay.size());
  delete r1;
  f = new TFile(filenamesforHIOverlay.at(filenumber).Data());
  fMCPico = (TTree *)f->Get("PicoDst");
  /*if (fPrintLevel)
    cout << filenumber << "\t" << filenamesforHIOverlay.at(filenumber) << "\t" << fMCPico->GetEntriesFast() << "\t" << mCentMaker->GetWeight() << endl;*/

  fMCPico->SetMakeClass(1);
  fMCPico->SetBranchStatus("*", false);
  fMCPico->SetBranchStatus("Event", true);
  fMCPico->SetBranchStatus("Event.mRunId", true);
  fMCPico->SetBranchStatus("Event.mEventId", true);

  fMCPico->SetBranchStatus("Event.mPrimaryVertexX", true);
  fMCPico->SetBranchStatus("Event.mPrimaryVertexY", true);
  fMCPico->SetBranchStatus("Event.mPrimaryVertexZ", true);

  fMCPico->SetBranchStatus("Track", true);
  fMCPico->SetBranchStatus("Track.mGMomentumX", true);
  fMCPico->SetBranchStatus("Track.mGMomentumY", true);
  fMCPico->SetBranchStatus("Track.mGMomentumZ", true);
  fMCPico->SetBranchStatus("Track.mOriginX", true);
  fMCPico->SetBranchStatus("Track.mOriginY", true);
  fMCPico->SetBranchStatus("Track.mOriginZ", true);
  fMCPico->SetBranchStatus("Track.mNHitsFit", true);
  fMCPico->SetBranchStatus("Track.mNHitsMax", true);
  fMCPico->SetBranchStatus("Track.mNSigmaPion", true);
  fMCPico->SetBranchStatus("Track.mNSigmaKaon", true);
  fMCPico->SetBranchStatus("Track.mNSigmaProton", true);
  fMCPico->SetBranchStatus("Track.mNSigmaElectron", true);
  fMCPico->SetBranchStatus("Track.mBEmcMatchedTowerIndex", true);
  fMCPico->SetBranchStatus("Track.mIdTruth", true);
  fMCPico->SetBranchStatus("Track.mQATruth", true);

  fMCPico->SetBranchStatus("BTowHit", true);
  fMCPico->SetBranchStatus("BTowHit.mE", true);

  fMCPico->SetBranchStatus("McVertex", true);
  fMCPico->SetBranchStatus("McVertex.mId", true);
  fMCPico->SetBranchStatus("McVertex.mNoDaughters", true);
  fMCPico->SetBranchStatus("McVertex.mIdParTrk", true);
  fMCPico->SetBranchStatus("McVertex.mIsInterm", true);
  fMCPico->SetBranchStatus("McVertex.mTime", true);
  fMCPico->SetBranchStatus("McVertex.mVx", true);
  fMCPico->SetBranchStatus("McVertex.mVy", true);
  fMCPico->SetBranchStatus("McVertex.mVz", true);
  fMCPico->SetBranchStatus("McTrack", true);
  fMCPico->SetBranchStatus("McTrack.mId", true);
  fMCPico->SetBranchStatus("McTrack.mGePid", true);
  fMCPico->SetBranchStatus("McTrack.mCharge", true);
  fMCPico->SetBranchStatus("McTrack.mHits[22]", true);
  fMCPico->SetBranchStatus("McTrack.mPx", true);
  fMCPico->SetBranchStatus("McTrack.mPy", true);
  fMCPico->SetBranchStatus("McTrack.mPz", true);
  fMCPico->SetBranchStatus("McTrack.mE", true);
  fMCPico->SetBranchStatus("McTrack.mIsFromShower", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxStart", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxStop", true);
  fMCPico->SetBranchStatus("McTrack.mIdVtxItrmd", true);

  fMCPico->SetBranchAddress("Event", &Event_, &b_Event_);
  fMCPico->SetBranchAddress("Event.mRunId", Event_mRunId, &b_Event_mRunId);
  fMCPico->SetBranchAddress("Event.mEventId", Event_mEventId, &b_Event_mEventId);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexX", Event_mPrimaryVertexX, &b_Event_mPrimaryVertexX);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexY", Event_mPrimaryVertexY, &b_Event_mPrimaryVertexY);
  fMCPico->SetBranchAddress("Event.mPrimaryVertexZ", Event_mPrimaryVertexZ, &b_Event_mPrimaryVertexZ);

  fMCPico->SetBranchAddress("Track", &Track_, &b_Track_);
  fMCPico->SetBranchAddress("Track.mGMomentumX", Track_mGMomentumX, &b_Track_mGMomentumX);
  fMCPico->SetBranchAddress("Track.mGMomentumY", Track_mGMomentumY, &b_Track_mGMomentumY);
  fMCPico->SetBranchAddress("Track.mGMomentumZ", Track_mGMomentumZ, &b_Track_mGMomentumZ);
  fMCPico->SetBranchAddress("Track.mOriginX", Track_mOriginX, &b_Track_mOriginX);
  fMCPico->SetBranchAddress("Track.mOriginY", Track_mOriginY, &b_Track_mOriginY);
  fMCPico->SetBranchAddress("Track.mOriginZ", Track_mOriginZ, &b_Track_mOriginZ);
  fMCPico->SetBranchAddress("Track.mNHitsFit", Track_mNHitsFit, &b_Track_mNHitsFit);
  fMCPico->SetBranchAddress("Track.mNHitsMax", Track_mNHitsMax, &b_Track_mNHitsMax);
  fMCPico->SetBranchAddress("Track.mNSigmaPion", Track_mNSigmaPion, &b_Track_mNSigmaPion);
  fMCPico->SetBranchAddress("Track.mNSigmaKaon", Track_mNSigmaKaon, &b_Track_mNSigmaKaon);
  fMCPico->SetBranchAddress("Track.mNSigmaProton", Track_mNSigmaProton, &b_Track_mNSigmaProton);
  fMCPico->SetBranchAddress("Track.mNSigmaElectron", Track_mNSigmaElectron, &b_Track_mNSigmaElectron);
  fMCPico->SetBranchAddress("Track.mBEmcMatchedTowerIndex", Track_mBEmcMatchedTowerIndex, &b_Track_mBEmcMatchedTowerIndex);
  fMCPico->SetBranchAddress("Track.mIdTruth", Track_mIdTruth, &b_Track_mIdTruth);
  fMCPico->SetBranchAddress("Track.mQATruth", Track_mQATruth, &b_Track_mQATruth);

  fMCPico->SetBranchAddress("BTowHit", &BTowHit_, &b_BTowHit_);
  fMCPico->SetBranchAddress("BTowHit.mE", BTowHit_mE, &b_BTowHit_mE);

  fMCPico->SetBranchAddress("McVertex", &McVertex_, &b_McVertex_);
  fMCPico->SetBranchAddress("McVertex.mId", McVertex_mId, &b_McVertex_mId);
  fMCPico->SetBranchAddress("McVertex.mNoDaughters", McVertex_mNoDaughters, &b_McVertex_mNoDaughters);
  fMCPico->SetBranchAddress("McVertex.mIdParTrk", McVertex_mIdParTrk, &b_McVertex_mIdParTrk);
  fMCPico->SetBranchAddress("McVertex.mIsInterm", McVertex_mIsInterm, &b_McVertex_mIsInterm);
  fMCPico->SetBranchAddress("McVertex.mTime", McVertex_mTime, &b_McVertex_mTime);
  fMCPico->SetBranchAddress("McVertex.mVx", McVertex_mVx, &b_McVertex_mVx);
  fMCPico->SetBranchAddress("McVertex.mVy", McVertex_mVy, &b_McVertex_mVy);
  fMCPico->SetBranchAddress("McVertex.mVz", McVertex_mVz, &b_McVertex_mVz);

  fMCPico->SetBranchAddress("McTrack", &McTrack_, &b_McTrack_);
  fMCPico->SetBranchAddress("McTrack.mId", McTrack_mId, &b_McTrack_mId);
  fMCPico->SetBranchAddress("McTrack.mGePid", McTrack_mGePid, &b_McTrack_mGePid);
  fMCPico->SetBranchAddress("McTrack.mCharge", McTrack_mCharge, &b_McTrack_mCharge);
  fMCPico->SetBranchAddress("McTrack.mHits[22]", McTrack_mHits, &b_McTrack_mHits);
  fMCPico->SetBranchAddress("McTrack.mPx", McTrack_mPx, &b_McTrack_mPx);
  fMCPico->SetBranchAddress("McTrack.mPy", McTrack_mPy, &b_McTrack_mPy);
  fMCPico->SetBranchAddress("McTrack.mPz", McTrack_mPz, &b_McTrack_mPz);
  fMCPico->SetBranchAddress("McTrack.mE", McTrack_mE, &b_McTrack_mE);
  fMCPico->SetBranchAddress("McTrack.mIsFromShower", McTrack_mIsFromShower, &b_McTrack_mIsFromShower);
  fMCPico->SetBranchAddress("McTrack.mIdVtxStart", McTrack_mIdVtxStart, &b_McTrack_mIdVtxStart);
  fMCPico->SetBranchAddress("McTrack.mIdVtxStop", McTrack_mIdVtxStop, &b_McTrack_mIdVtxStop);
  fMCPico->SetBranchAddress("McTrack.mIdVtxItrmd", McTrack_mIdVtxItrmd, &b_McTrack_mIdVtxItrmd);
}

void StHIOverlayAngularities::OutputTreeInit()
{
  TString treename = "jets";
  outputTree = new TTree(treename.Data(), treename.Data());

  // Branches to save event info
  outputTree->Branch("centrality", &fRecoJetTree.centrality, "centrality/F");
  outputTree->Branch("centralityAlt", &fRecoJetTree.centralityAlt, "centralityAlt/F");
  outputTree->Branch("weightCentrality", &fRecoJetTree.weight, "weightCentrality/F");
  //outputTree->Branch("RefMult", &fRecoJetTree.refmult, "refmult/F");
  outputTree->Branch("gRefMult", &fRecoJetTree.grefmult, "grefmult/F");
  //outputTree->Branch("RefCorr2", &fRecoJetTree.refcorr2, "refcorr2/F");
  if(fPhiBgModulation){
	  outputTree->Branch("Psi2", &fRecoJetTree.Psi2, "Psi2/F");
	  outputTree->Branch("corrPsi2", &fRecoJetTree.corrPsi2, "corrPsi2/F");
	  outputTree->Branch("Q_1Vec", &fRecoJetTree.Q1_vec, "Q_1Vec/F");
	  outputTree->Branch("Q_2Vec", &fRecoJetTree.Q2_vec, "Q_2Vec/F");
	  outputTree->Branch("Q_1VecRec", &fRecoJetTree.Q1_vec_rec, "Q_1VecRec/F");
	  outputTree->Branch("Q_2VecRec", &fRecoJetTree.Q2_vec_rec, "Q_2VecRec/F");
  }
  
  
  outputTree->Branch("mcRefMult", &fMcJetTree.refmult, "mcrefmult/F");
  outputTree->Branch("mcRecoRefMult", &fMcJetTree.grefmult, "recorefmult/F");
  
  // outputTree->Branch("McPrimaryVertex", &fMcJetTree.primaryvertex);
  // outputTree->Branch("RecoPrimaryVertex", &fRecoJetTree.primaryvertex);

  outputTree->Branch("mcD0Pt", &fMcJetTree.d0pt, "mcD0Pt/F");
  outputTree->Branch("mcD0Eta", &fMcJetTree.d0eta, "mcD0Eta/F");
  outputTree->Branch("mcD0Phi", &fMcJetTree.d0phi, "mcD0Phi/F");
  outputTree->Branch("mcPionPt", &fMcJetTree.pionpt, "mcPionPt/F");
  outputTree->Branch("mcPionEta", &fMcJetTree.pioneta, "mcPionEta/F");
  outputTree->Branch("mcPionPhi", &fMcJetTree.pionphi, "mcPionPhi/F");
  outputTree->Branch("mcKaonPt", &fMcJetTree.kaonpt, "mcKaonPt/F");
  outputTree->Branch("mcKaonEta", &fMcJetTree.kaoneta, "mcKaonEta/F");
  outputTree->Branch("mcKaonPhi", &fMcJetTree.kaonphi, "mcKaonPhi/F");

  outputTree->Branch("mcSmearedD0Pt", &fRecoJetTree.d0pt, "mcSmearedD0Pt/F");
  outputTree->Branch("mcSmearedD0Eta", &fRecoJetTree.d0eta, "mcSmearedD0Eta/F");
  outputTree->Branch("mcSmearedD0Phi", &fRecoJetTree.d0phi, "mcSmearedD0Phi/F");
  outputTree->Branch("mcSmearedPionPt", &fRecoJetTree.pionpt, "mcSmearedPionPt/F");
  outputTree->Branch("mcSmearedPionEta", &fRecoJetTree.pioneta, "mcSmearedPionEta/F");
  outputTree->Branch("mcSmearedPionPhi", &fRecoJetTree.pionphi, "mcSmearedPionPhi/F");
  outputTree->Branch("mcSmearedKaonPt", &fRecoJetTree.kaonpt, "mcSmearedKaonPt/F");
  outputTree->Branch("mcSmearedKaonEta", &fRecoJetTree.kaoneta, "mcSmearedKaonEta/F");
  outputTree->Branch("mcSmearedKaonPhi", &fRecoJetTree.kaonphi, "mcSmearedKaonPhi/F");

  outputTree->Branch("mcJetPt", &fMcJetTree.jetpt, "mcJetPt/F");
  outputTree->Branch("mcJetEta", &fMcJetTree.jeteta, "mcJetEta/F");
  outputTree->Branch("mcJetPhi", &fMcJetTree.jetphi, "mcJetPhi/F");
  outputTree->Branch("mcJetArea", &fMcJetTree.jetarea, "mcJetArea/F");
  outputTree->Branch("mcJetE", &fMcJetTree.jetenergy, "mcJetE/F");
  outputTree->Branch("mcJetNConst", &fMcJetTree.numberofconstituents, "mcJetNConst/I");
  outputTree->Branch("mcJetLambda1_0_5", &fMcJetTree.lambda_1_half, "mcJetLambda1_0_5/F");
  outputTree->Branch("mcJetLambda1_1", &fMcJetTree.lambda_1_1, "mcJetLambda1_1/F");
  outputTree->Branch("mcJetLambda1_1_5", &fMcJetTree.lambda_1_1half, "mcJetLambda1_1_5/F");
  outputTree->Branch("mcJetLambda1_2", &fMcJetTree.lambda_1_2, "mcJetLambda1_2/F");
  outputTree->Branch("mcJetLambda1_3", &fMcJetTree.lambda_1_3, "mcJetLambda1_3/F");
  outputTree->Branch("mcJetMomDisp", &fMcJetTree.dispersion, "mcJetMomDisp/F");
  outputTree->Branch("mcJetD0Z", &fMcJetTree.d0z, "mcJetD0Z/F");
  outputTree->Branch("mcJetD0DeltaR", &fMcJetTree.d0DeltaR, "mcJetD0DeltaR/F");
  
  outputTree->Branch("mcSmearedJetPt", &fMcRecoJetTree.jetpt, "mcSmearedJetPt/F"); 
  outputTree->Branch("mcSmearedJetEta", &fMcRecoJetTree.jeteta, "mcSmearedJetEta/F");
  outputTree->Branch("mcSmearedJetPhi", &fMcRecoJetTree.jetphi, "mcSmearedJetPhi/F");
  outputTree->Branch("mcSmearedJetE", &fMcRecoJetTree.jetenergy, "mcSmearedJetE/F");
  outputTree->Branch("mcSmearedJetArea", &fMcRecoJetTree.jetarea, "mcSmearedJetArea/F");
  outputTree->Branch("mcSmearedJetNConst", &fMcRecoJetTree.numberofconstituents, "mcSmearedJetNConst/I");
  outputTree->Branch("mcSmearedJetLambda1_0_5", &fMcRecoJetTree.lambda_1_half, "mcSmearedJetLambda1_0_5/F");
  outputTree->Branch("mcSmearedJetLambda1_1", &fMcRecoJetTree.lambda_1_1, "mcSmearedJetLambda1_1/F");
  outputTree->Branch("mcSmearedJetLambda1_1_5", &fMcRecoJetTree.lambda_1_1half, "mcSmearedJetLambda1_1_5/F");
  outputTree->Branch("mcSmearedJetLambda1_2", &fMcRecoJetTree.lambda_1_2, "mcSmearedJetLambda1_2/F");
  outputTree->Branch("mcSmearedJetLambda1_3", &fMcRecoJetTree.lambda_1_3, "mcSmearedJetLambda1_3/F");
  outputTree->Branch("mcSmearedJetMomDisp", &fMcRecoJetTree.dispersion, "mcSmearedJetMomDisp/F");
  outputTree->Branch("mcSmearedJetD0Z", &fMcRecoJetTree.d0z, "d0z/F");
  outputTree->Branch("mcSmearedJetD0DeltaR", &fMcRecoJetTree.d0DeltaR, "mcJetSmearedD0DeltaR/F");
  outputTree->Branch("recoJetPt", &fRecoJetTree.jetpt, "jetpt/F");
  if (fBgSubtraction == 0 || fBgSubtraction == 2) outputTree->Branch("recoJetPtCorr", &fRecoJetTree.jetptcorr, "jetptcorr/F");
  outputTree->Branch("recoJetRho", &fRecoJetTree.jetrho, "jetrho/F");
  outputTree->Branch("recoJetRhoM", &fRecoJetTree.jetrhom, "jetrhom/F");
  outputTree->Branch("recoJetEta", &fRecoJetTree.jeteta, "jeteta/F");
  outputTree->Branch("recoJetPhi", &fRecoJetTree.jetphi, "jetphi/F");
  outputTree->Branch("recoJetArea", &fRecoJetTree.jetarea, "jetarea/F");
  outputTree->Branch("recoJetE", &fRecoJetTree.jetenergy, "jetenergy/F");
  //outputTree->Branch("recoJetRhoVal", &fRecoJetTree.fRhoValforjet, "fRhoValforjet/F");
  outputTree->Branch("recoJetNConst", &fRecoJetTree.numberofconstituents, "numberofconstituents/I");
  outputTree->Branch("recoJetLambda1_0_5", &fRecoJetTree.lambda_1_half, "lambda_1_half/F");
  outputTree->Branch("recoJetLambda1_1", &fRecoJetTree.lambda_1_1, "recoJetLambda1_1/F");
  outputTree->Branch("recoJetLambda1_1_5", &fRecoJetTree.lambda_1_1half, "recoJetLambda1_1_5/F");
  outputTree->Branch("recoJetLambda1_2", &fRecoJetTree.lambda_1_2, "recoJetLambda1_2/F");
  outputTree->Branch("recoJetLambda1_3", &fRecoJetTree.lambda_1_3, "recoJetLambda1_3/F");
  outputTree->Branch("recoJetMomDisp", &fRecoJetTree.dispersion, "recoJetMomDisp/F");
  outputTree->Branch("recoJetD0Z", &fRecoJetTree.d0z, "recoJetD0Z/F");
  outputTree->Branch("recoJetD0DeltaR", &fRecoJetTree.d0DeltaR, "recoJetD0DeltaR/F");
  
}

//----------------------------------
void StHIOverlayAngularities::FillTree(const Int_t &numberOfD0Events)
{
   //// cerr << "----------------------------------------------" << endl;
    ////cout << "mPicoDst->event()->eventId(): " << mPicoDst->event()->eventId() << endl;

  for (Int_t iMcD0Event = 0; iMcD0Event < numberOfD0Events; iMcD0Event++)
  {
    fMcJetTree.Clear();
    fRecoJetTree.Clear();
    fMcRecoJetTree.Clear();

    fRecoJetTree.centrality = fCentrality;
    fRecoJetTree.centralityAlt = fCentralityAlt;
    ////fRecoJetTree.weight = mCentMaker->GetWeight();
    fRecoJetTree.weight = fCentralityWeight;
    TVector3 mcPion = fMcD0Information[iMcD0Event].first;
    TVector3 mcKaon = fMcD0Information[iMcD0Event].second;

    TVector3 mcD0 = mcPion + mcKaon;

    fMcJetTree.d0pt = mcD0.Perp();
    fMcJetTree.d0eta = mcD0.PseudoRapidity();
    fMcJetTree.d0phi = standardPhi(mcD0.Phi());

    if (fPrintLevel)
      cout << "MC D0 Saved = " << fMcJetTree.d0pt << "\t" << fMcJetTree.d0eta << "\t" << fMcJetTree.d0phi << endl;

    fMcJetTree.pionpt = mcPion.Perp();
    fMcJetTree.pioneta = mcPion.PseudoRapidity();
    fMcJetTree.pionphi = standardPhi(mcPion.Phi());

    fMcJetTree.kaonpt = mcKaon.Perp();
    fMcJetTree.kaoneta = mcKaon.PseudoRapidity();
    fMcJetTree.kaonphi = standardPhi(mcKaon.Phi());

    fMcJetTree.refmult = fMcEventTracks[iMcD0Event].size();
    fMcJetTree.grefmult = fRecoMcEventTracks[iMcD0Event].size();

    TVector3 recoPion = fMcRecoD0Information[iMcD0Event].first;
    TVector3 recoKaon = fMcRecoD0Information[iMcD0Event].second;

    TVector3 recoD0 = recoPion + recoKaon;

    fRecoJetTree.d0pt = recoD0.Perp();
    fRecoJetTree.d0eta = recoD0.PseudoRapidity();
    fRecoJetTree.d0phi = standardPhi(recoD0.Phi());

    fRecoJetTree.pionpt = recoPion.Perp();
    fRecoJetTree.pioneta = recoPion.PseudoRapidity();
    fRecoJetTree.pionphi = standardPhi(recoPion.Phi());

    fRecoJetTree.kaonpt = recoKaon.Perp();
    fRecoJetTree.kaoneta = recoKaon.PseudoRapidity();
    fRecoJetTree.kaonphi = standardPhi(recoKaon.Phi());
    fRecoJetTree.grefmult = fKgrefMult_uncorr;
    
    /*
    fRecoJetTree.Q1_vec = fQ_1;
    fRecoJetTree.Q2_vec = fQ_2;
    fRecoJetTree.Q1_vec_rec = fQ_1_rec;
    fRecoJetTree.Q2_vec_rec = fQ_2_rec;
    fRecoJetTree.Psi2 = fPsi_2;
    fRecoJetTree.corrPsi2 = fPsi_2_shifted;
    */
    
    StJet *jetMc = fMcJet[iMcD0Event];
    StJet *jetMcReco = fMcRecoJet[iMcD0Event];
    StJet *jetReco = fRecoJet[iMcD0Event];


    //Pure MC
    ///cout << "STAGE I" << endl;
    FillJet(jetMc, fMcJetTree, mcD0);
    
    //Smeared MC
    ///cout << "STAGE II" << endl;
    if (fMcRecoJet[iMcD0Event] != NULL) FillJet(jetMcReco, fMcRecoJetTree, recoD0);
    
    //Embedded MC
    ///cout << "STAGE III" << endl;
    if (fRecoJet[iMcD0Event] != NULL){
    
      fRecoJetTree.fRhoValforjet = fRhoVal;
            //aditional ones (outside because I want it only for reco)
      fRecoJetTree.jetrho = jetReco->GetJetRho();
      fRecoJetTree.jetrhom = jetReco->GetJetRhoM();
      
      Double_t zD0Corr = -999;
      Double_t pt_corr = -999;
      Double_t jetPt = -999;
      Double_t lam_1_0half = -999;
      Double_t lam_1_1 = -999;
      Double_t lam_1_1half = -999;
      Double_t lam_1_2 = -999;      
      Double_t lam_1_3 = -999; 
      Double_t momDisp = -999;           
            
      if (fBgSubtraction == 0 || fBgSubtraction == 2){ 
      
		pt_corr = jetReco->Pt() - jetReco->GetJetRho() * jetReco->Area();
		///pt_corr = jetReco->Pt() - jetReco->GetJetRho() * jetReco->GetArea4Vector().pt();
	
		jetPt = jetReco->Pt();
      		Double_t px_corr = jetReco->Px() - jetReco->GetJetRho() * jetReco->GetArea4Vector().px();
            	Double_t py_corr = jetReco->Py() - jetReco->GetJetRho() * jetReco->GetArea4Vector().py();
            	zD0Corr = (px_corr*recoD0.X() + py_corr*recoD0.Y()) / (pt_corr * pt_corr);
            	
            	
            	/*
		zD0Corr = (pt_corr*recoD0.X()*TMath::Cos(jetReco->Phi())+pt_corr*recoD0.Y()*TMath::Sin(jetReco->Phi()))
		/ (pt_corr * pt_corr);*/
            	
            	fRecoJetTree.jetptcorr = pt_corr;
            	
            	
            	if(fBgSubtraction == 2){
		    	lam_1_0half = jetReco->GetActualShapeL10half();
	      		lam_1_1 = jetReco->GetActualShapeL11();
	      		lam_1_1half = jetReco->GetActualShapeL11half();
	      		lam_1_2 = jetReco->GetActualShapeL12();     
	      		lam_1_3 = jetReco->GetActualShapeL13(); 
	      		momDisp = jetReco->GetActualShapeDisp();
      		}
       
      }
      FillJet(jetReco, fRecoJetTree, recoD0);
     // cout << "Jet pT corr: " << pt_corr << endl;
      
      //In case of area based background subtraction, the observables have to be calculated differently
      if (fBgSubtraction == 0) {
      
      	fRecoJetTree.d0z = zD0Corr;
        fRecoJetTree.lambda_1_half 	*= 1.*jetPt / pt_corr;
 	fRecoJetTree.lambda_1_1 	*= 1.*jetPt / pt_corr;
  	fRecoJetTree.lambda_1_1half     *= 1.*jetPt / pt_corr;
  	fRecoJetTree.lambda_1_2         *= 1.*jetPt / pt_corr;
  	fRecoJetTree.lambda_1_3         *= 1.*jetPt / pt_corr;
  	fRecoJetTree.dispersion         *= 1.*jetPt*jetPt / pt_corr / pt_corr;
  	
      
      } else if (fBgSubtraction == 2){

        fRecoJetTree.d0z = zD0Corr;
        fRecoJetTree.lambda_1_half 	= lam_1_0half;
 	fRecoJetTree.lambda_1_1 	= lam_1_1;
  	fRecoJetTree.lambda_1_1half     = lam_1_1half;
  	fRecoJetTree.lambda_1_2         = lam_1_2;
  	fRecoJetTree.lambda_1_3         = lam_1_3;
  	fRecoJetTree.dispersion         = momDisp;
      
      }

    
    }
    

    ////exit(0);
    outputTree->Fill();

    if (fRecoJetTree.numberofconstituents > 0){

           hResponseJetPt->Fill(fRecoJetTree.jetpt, fMcJetTree.jetpt, fCentralityWeight);
	   hResponseJetD0Z->Fill(fRecoJetTree.d0z, fMcJetTree.d0z,fCentralityWeight);
	   hResponseJetNConst->Fill(fRecoJetTree.numberofconstituents, fMcJetTree.numberofconstituents,fCentralityWeight);
	   hResponseJetLambda1_0_5->Fill(fRecoJetTree.lambda_1_half, fMcJetTree.lambda_1_half,fCentralityWeight);
	   hResponseJetLambda1_1->Fill(fRecoJetTree.lambda_1_1, fMcJetTree.lambda_1_1,fCentralityWeight);
	   hResponseJetLambda1_1_5->Fill(fRecoJetTree.lambda_1_1half, fMcJetTree.lambda_1_1half,fCentralityWeight);
	   hResponseJetLambda1_2->Fill(fRecoJetTree.lambda_1_2, fMcJetTree.lambda_1_2,fCentralityWeight);
	   hResponseJetLambda1_3->Fill(fRecoJetTree.lambda_1_3, fMcJetTree.lambda_1_3,fCentralityWeight);
	   hResponseJetMomDisp->Fill(fRecoJetTree.dispersion, fMcJetTree.dispersion,fCentralityWeight);
	   hResponseJetD0DeltaR->Fill(fRecoJetTree.d0DeltaR, fMcJetTree.d0DeltaR, fCentralityWeight);
	   hResponseJetD0Pt->Fill(fRecoJetTree.d0pt, fMcJetTree.d0pt, fCentralityWeight);
	      
	   if (centralitybinforefficiency==0) etaRecoTrue00_10->Fill(fRecoJetTree.jeteta, fMcJetTree.jeteta);
	   if (centralitybinforefficiency==1) etaRecoTrue10_40->Fill(fRecoJetTree.jeteta, fMcJetTree.jeteta);
	   if (centralitybinforefficiency==2) etaRecoTrue40_80->Fill(fRecoJetTree.jeteta, fMcJetTree.jeteta);
	   
	   if (centralitybinforefficiency==0) effJetRec00_10->Fill(true, fMcJetTree.jetpt);
		   else if (centralitybinforefficiency==1) effJetRec10_40->Fill(true, fMcJetTree.jetpt);
		   else if (centralitybinforefficiency==2) effJetRec40_80->Fill(true, fMcJetTree.jetpt);
    	   } else {
    	   if (centralitybinforefficiency==0) effJetRec00_10->Fill(false, fMcJetTree.jetpt);
		   else if (centralitybinforefficiency==1) effJetRec10_40->Fill(false, fMcJetTree.jetpt);
		   else if (centralitybinforefficiency==2) effJetRec40_80->Fill(false, fMcJetTree.jetpt);
           }
    
    
  }
  if (fPrintLevel)
    cout << "Filling Tree complete" << endl;
}

Double_t StHIOverlayAngularities::GetTowerCalibEnergy(Int_t TowerId){
  //Function calculates the calibrated energy of a tower

  //Loading of the tower
  StPicoBTowHit *tower = static_cast<StPicoBTowHit*>(mPicoDst->btowHit(TowerId-1)); //ID

  //Initialization of the pedestal, rms and status
  Float_t pedestal, rms;
  Int_t status;

  //Loading of the pedestal, rms and status (it does not work, if you use root instead of root4star)
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms);
  mTables->getStatus(BTOW, TowerId, status);

  //Initialization of the tower coefficients
  const Double_t *TowerCoeff = nullptr;

  //Tower coefficients for the different runs, parameters are saved in BemcNewCalib.h
  if(fRunNumber <= 15094020) TowerCoeff = mycuts::CPre;
  else TowerCoeff = mycuts::CLowMidHigh;

  //Calculation of the calibrated energy E=C*(ADC-Pedestal)
  Double_t calibEnergy = TowerCoeff[TowerId-1]*(tower->adc() - pedestal);


  return calibEnergy;

  //Function returns the calibrated energy of the tower
}
Bool_t StHIOverlayAngularities::AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert) {
  // constants: assume neutral pion mass
  ///double pi0mass = Pico::mMass[0]; // GeV
  double pi = 1.0*TMath::Pi();

  // get momentum vector of track - global or primary track
  TVector3 mTrkMom;
  if(doUsePrimTracks) {
    if(!(trk->isPrimary())) return kFALSE; // check if primary
    // get primary track vector
    mTrkMom = trk->pMom();
  } else {
    // get global track vector
    mTrkMom = trk->gMom(Vert, B);
  }

  // track variables
  //double pt = mTrkMom.Perp();
  double phi = mTrkMom.Phi();
  double eta = mTrkMom.PseudoRapidity();
  double dca = trk->gDCA(Vert).Mag();
  int nHitsFit = trk->nHitsFit();
  int nHitsMax = trk->nHitsMax();
  double nHitsRatio = 1.0*nHitsFit/nHitsMax;

  // jet track acceptance cuts now - note difference from AcceptJetTrack()
  if((eta < -1) || (eta > 1)) return kFALSE;
  if(phi < 0.0)    phi += 2.0*pi;  // force from 0-2pi
  if(phi > 2.0*pi) phi -= 2.0*pi;  // force from 0-2pi

  // additional quality cuts for tracks
  if(dca > fJetTrackDCAcut)            return kFALSE;
  if(nHitsFit < fJetTracknHitsFit)     return kFALSE;
  if(nHitsRatio < fJetTracknHitsRatio) return kFALSE;

  // passed all above cuts - keep track and fill input vector to fastjet
  return kTRUE;
}
void StHIOverlayAngularities::FillJet(StJet *jet, StJetTreeStruct &jetTree, const TVector3 &D0)
{
  jetTree.d0z = (jet->Px() * D0.X() + jet->Py() * D0.Y()) / (jet->Pt() * jet->Pt());
  
  
  
  jetTree.d0DeltaR = sqrt(pow(D0.Eta()-jet->Eta(),2)+pow(TVector2::Phi_mpi_pi(D0.Phi() - jet->Phi()) ,2));
  jetTree.jetpt = jet->Pt();
  jetTree.jeteta = jet->Eta();
  /*
  cout << "Jet pT = " << jet->Pt() << endl;
  cout << "Jet eta = " << jet->Eta() << endl;
  cout << "D0 pT = " << D0.Perp() << endl;
  cout << "D0 Eta = " << D0.Eta() << endl;
  cout << "rho = " << fjw->GetJetRho() << endl;
  cout << "area = " << jet->Area() << endl; 
  */
  jetTree.jetphi = jet->Phi();
  jetTree.jetarea = jet->Area();
  jetTree.jetenergy = jet->E();
  jetTree.numberofconstituents = jet->GenNumberOfCons(false, false); //false = ignore ghost particles and particles with pT < 0.00001
   //   cout << "NOC = " << jet->GenNumberOfCons(false, false) << endl;
  jetTree.lambda_1_half = jet->GetAngularityLambda(1., 0.5, 0.4);
  jetTree.lambda_1_1 = jet->GetAngularityLambda(1., 1., 0.4);
  jetTree.lambda_1_1half = jet->GetAngularityLambda(1., 1.5, 0.4);
  jetTree.lambda_1_2 = jet->GetAngularityLambda(1., 2, 0.4);
  jetTree.lambda_1_3 = jet->GetAngularityLambda(1., 3, 0.4);
  jetTree.dispersion = sqrt(jet->GetAngularityLambda(2., 0., 0.4));
  
  
  

  
  
  delete jet;
}
