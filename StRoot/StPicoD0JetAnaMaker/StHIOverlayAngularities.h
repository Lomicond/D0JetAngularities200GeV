#ifndef StHIOverlayAngularities_h
#define StHIOverlayAngularities_h

#include "StEmcUtil/geometry/StEmcGeom.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TBranch.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"

#include "TGraph.h"
#include "StJetTreeStruct.h"
#include "StFJWrapper.h"
#include "FJ_includes.h"
#include "StJet.h"

#include "StMaker.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include <set>

class StFJWrapper;

// ROOT classes
class StJet;

class TString;
class TVector3;
class TLorentzVector;
class TString;
class TTree;
class TBranch;
class TGraph;

// STAR classes
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StRefMultCorr;
class StEmcPosition;
class StPicoBTowHit;

// jet-framework classes
class StCentMaker;
class StJet;

// BEMC tables
class StEmcADCtoEMaker;
class StBemcTables;

// ROOT classes
class TClonesArray;
class TF1;
class TH1;
class TH1F;
class TH2;
class TH2F;
class TH3;
class THnSparse;
class TString;


class StHIOverlayAngularities : public StMaker{

    // jet type enumerator
    enum EJetType_t {
      kFullJet,
      kChargedJet,
      kNeutralJet
    };

    // jet algorithm enumerator
    enum EJetAlgo_t {
      kt_algorithm                    = 0,
      antikt_algorithm                = 1,
      cambridge_algorithm             = 2,
      genkt_algorithm                 = 3,
      cambridge_for_passive_algorithm = 11,
      genkt_for_passive_algorithm     = 13,
      plugin_algorithm                = 99,
      undefined_jet_algorithm         = 999
    };

    // jet recombination scheme enumerator
    enum ERecoScheme_t {
      E_scheme        = 0,
      pt_scheme       = 1,
      pt2_scheme      = 2,
      Et_scheme       = 3,
      Et2_scheme      = 4,
      BIpt_scheme     = 5,
      BIpt2_scheme    = 6,
      WTA_pt_scheme   = 7,
      WTA_modp_scheme = 8,
      external_scheme = 99
    };

protected:

    //Run cuts
    Int_t                   fRunBadlist = 0;
    
    //Event cuts
    Int_t 	            fNumberOfEventsToOverLay;
    Double_t                fEventZVtxMinCut;        // min event z-vertex cut
    Double_t                fEventZVtxMaxCut;        // max event z-vertex cut
    Double_t 	            fEventRVtxMaxCut;
    Double_t		    fEventVtxVpdVzMaxCut;
    std::set<Int_t> 	    fEventTriggers; 		
    
    //MC
    Bool_t                  fSetMcSeed = false;
    Int_t                   fMcSeed = 0;
    
    //Fastjet
    Bool_t                  fSetJetFixedSeed = false;
    Int_t                   fJetFJSeed = 12345;

    //Jet Cuts
    Double_t 		    fGhostArea;  // ghost area
    Double_t                fJetRad;                 // jet radius
    Bool_t 	            fSetJetFracCut[3];
    Double_t 		    fJetFractionNeutralToJetPt = 0.95;
    Bool_t 	       	    fSetJetMinPtCut[3];
    Double_t 		    fJetMinPtCut = 1;
    Bool_t 	       	    fSetJetMinAreaCut[3];
    Double_t 		    fJetMinAreaCut = 0;    
    Bool_t 	       	    fSetJetMinAbsEtaCut[3];
    Double_t 		    fJetMinAbsEtaCut = 1;  
    
    //Jet Constituents - D0
    Double_t 		    fMinMcPtD0;
    Bool_t 	  	    fSetMcAbsEtaD0 = false;
    Double_t                fMcAbsEtaD0 = 1;
    Double_t 		    fMcD0Mass;   
    
    //Jet Constituents - Charged
    Double_t 		    fMinJetTrackPt; // min jet track transverse momentum cut
    Double_t 		    fMaxJetTrackPt; // max jet track transverse momentum cut
    Double_t                fJetTrackEtaMin;     // min jet track eta cut
    Double_t                fJetTrackEtaMax;     // max jet track eta cut
    Double_t                fJetTrackPhiMin;     // min jet track phi cut
    Double_t                fJetTrackPhiMax;     // max jet track phi cut
    Double_t                fJetTrackDCAcut;     // max jet track dca cut
    Double_t 		    fHadrCorrTrackDCAZcut;     // max jet track dca cut
    Bool_t                  doUsePrimTracks;         // primary track switch    
    Double_t 	            fChargedPart;
    Double_t 	            fMcChargedPart;  
    Int_t                   fTracknHitsFit;          // requirement for track hits  
    Double_t                fTracknHitsRatio;        // requirement for nHitsFit / nHitsMax
        
    //Jet Constituents - Towers    
    Int_t                   fTowerBadlist = 0;
    Bool_t 		    fSetTowerCalibrEnergy = false;    
    Double_t                fNeutralPart;
    Double_t                fMcNeutralPart;    
    
    //Jet Constituents - Hadronic Correction   
    Bool_t                  fSetDcaZHadronCorr = false; //false = DCA, true = DCA_z
    Double_t                fDcaZHadronCorr = 3.0;    
    Double_t                fHadronicCorrMass;    
    
    //Jet background Subtraction
    Int_t 		    fBgSubtraction;    
    Int_t                   fJetNHardestSkipped_010 = 2;
    Int_t                   fJetNHardestSkipped_1080 = 1;    
    
    //Centrality    
    Double_t		    fCentrality;
    Double_t		    fCentralityAlt;
    Double_t		    fCentralityWeight;
    Double_t		    fKgrefMult_uncorr;
    Double_t                fKrefMult_uncorr;
    
    //Event plane
    Double_t		    fPsi_2;
    Double_t 		    fQ_1;
    Double_t                fQ_2;
    Double_t                fPsi_2_shifted;
    Double_t                fQ_1_rec;
    Double_t                fQ_2_rec;

    //Event
    Double_t                Bfield;                  // event Bfield
    TVector3                mVertex;                 // event vertex 3-vector
    Double_t                zVtx;                    // z-vertex component
    
    //PicoDstMaker and PicoDst object pointer
    StPicoDstMaker          *mPicoDstMaker;
    StPicoDst               *mPicoDst;
    StPicoEvent             *mPicoEvent;

    //Position object
    StEmcPosition           *mEmcPosition;

    //Output file name string 
    TString                 mOutName;

private:
    StRefMultCorr* mGRefMultCorrUtil;
    StCentMaker           *mCentMaker;              // Centrality maker object //change

public:
  StHIOverlayAngularities(const char *name, StPicoDstMaker *picoMaker, const char *outName, const char *filename, StRefMultCorr* grefmultCorrUtil);
  virtual ~StHIOverlayAngularities();

  // class required functions
  virtual Int_t Init();
  virtual Int_t Make();
  virtual void Clear(Option_t *opt = "");
  virtual Int_t Finish();

  void SetNumberOfEventsToOverLay(Int_t a) { fNumberOfEventsToOverLay = a; }

    //Run cuts
    virtual void            SetJetRad(Double_t jrad)           { fJetRad           = jrad; } // jet radius 
    virtual void            SetUsePrimaryTracks(Bool_t P)      { doUsePrimTracks   = P; }
    virtual void            SetEventRVtxMaxCut(Double_t rVtx)  { fEventRVtxMaxCut  = rVtx; }
    virtual void            SetEventVtxVpdVzMaxCut(Double_t VpdVxVtx) {fEventVtxVpdVzMaxCut = VpdVxVtx;}
    virtual void 	    SetEventTriggers(const std::set<Int_t>& EventTriggers) { fEventTriggers = EventTriggers;}
    virtual void 	    SetMinMcPtD0(Double_t McPtD0) {fMinMcPtD0 = McPtD0;}
    
    virtual void 	    setMcAbsEtaD0(Bool_t setMcEtaD0, Double_t mcEtaD0) {
    				fSetMcAbsEtaD0 = setMcEtaD0;
    				fMcAbsEtaD0 = mcEtaD0;
    				}
    				
    virtual void            SetMcD0Mass(Double_t McD0Mass) {fMcD0Mass = McD0Mass;}
    virtual void            setAllMcSeedsToEventId(Bool_t McSeed) {fSetMcSeed = McSeed;}
    virtual void            setJetFractionNeutralToJetPt(Bool_t setJetFracCut0, Bool_t setJetFracCut1, Bool_t setJetFracCut2, Double_t jetFractionNeutralToJetPt){
    	 		    	fSetJetFracCut[0] = setJetFracCut0;
    	 		    	fSetJetFracCut[1] = setJetFracCut1;
    				fSetJetFracCut[2] = setJetFracCut2;
    				fJetFractionNeutralToJetPt = jetFractionNeutralToJetPt;
   			    } //!!!TODO
    virtual void            setJetMinPt(Bool_t setJetMinPtCut0, Bool_t setJetMinPtCut1, Bool_t setJetMinPtCut2, Double_t jetMinPtCut){
    	 		    	fSetJetMinPtCut[0] = setJetMinPtCut0;
    	 		    	fSetJetMinPtCut[1] = setJetMinPtCut1;
    				fSetJetMinPtCut[2] = setJetMinPtCut2;
    				fJetMinPtCut = jetMinPtCut;
   			    } //!!!TODO   	
   			    
   virtual void            setJetMinArea(Bool_t setJetMinAreaCut0, Bool_t setJetMinAreaCut1, Bool_t setJetMinAreaCut2, Double_t jetMinAreaCut){
    	 		    	fSetJetMinAreaCut[0] = setJetMinAreaCut0;
    	 		    	fSetJetMinAreaCut[1] = setJetMinAreaCut1;
    				fSetJetMinAreaCut[2] = setJetMinAreaCut2;
    				fJetMinAreaCut = jetMinAreaCut;
   			    } //!!!TODO      			    		    
   virtual void            setJetMinAbsEta(Bool_t setJetMinAbsEtaCut0, Bool_t setJetMinAbsEtaCut1, Bool_t setJetMinAbsEtaCut2, Double_t jetMinAbsEtaCut){
    	 		    	fSetJetMinAbsEtaCut[0] = setJetMinAbsEtaCut0;
    	 		    	fSetJetMinAbsEtaCut[1] = setJetMinAbsEtaCut1;
    				fSetJetMinAbsEtaCut[2] = setJetMinAbsEtaCut2;
    				fJetMinAbsEtaCut = jetMinAbsEtaCut;
   			    } //!!!TODO      
    
    virtual void            SetMcChargedPart(Double_t McChargedPart) {fMcChargedPart = McChargedPart;}
    virtual void            SetMcNeutralPart(Double_t McNeutralPart) {fMcNeutralPart = McNeutralPart;}
    virtual void            SetChargedPart(Double_t ChargedPart) {fChargedPart = ChargedPart;}
    virtual void            SetNeutralPart(Double_t NeutralPart) {fNeutralPart = NeutralPart;}
    virtual void            SetHadronicCorrMass(Double_t HadronicCorrMass) {fHadronicCorrMass = HadronicCorrMass;}

    
    virtual void            SetJetTracknHitsFit(Double_t JetTracknHitsFit) {fJetTracknHitsFit = JetTracknHitsFit;}
    virtual void            SetJetTracknHitsRatio(Double_t JetTracknHitsRatio) {fJetTracknHitsRatio = JetTracknHitsRatio;}
    
    virtual void            SetIBSMass(bool MassiveAll) {fMassiveAll = MassiveAll;}
    virtual void 	    SetIBSAlpha1(double ALFA_1) {fgAlpha1 = ALFA_1;}
    virtual void 	    SetIBSAlpha2(double ALFA_2) {fgAlpha2 = ALFA_2;}
    virtual void            SetPhiBgModulation(bool PhiBgModulation) {fPhiBgModulation = PhiBgModulation;} 



   virtual void SetPrintLevel(Int_t i) { fPrintLevel = i; }

  void SetJetTrackDCAcut(Double_t d) { fJetTrackDCAcut = d; }
  void SetHadrCorrTrackDCAZcut(Double_t cc) { fHadrCorrTrackDCAZcut = cc; }
  

  void DoTrackingEfficiency(const Double_t &percentage)
  {
    fTrackingEfficiency = kTRUE;
    fTrackingEfficiencyPercentage = percentage;
  } // percentage of tracks to keep
  
  void SetBgSubtraction(Int_t BgSubtraction){
  //0 - Area based, 1 - ICS, 2 - jet shape
       fBgSubtraction = BgSubtraction;
  }

  void SetJetAlgo(Int_t a) { fJetAlgo = a; }
  void SetJetType(Int_t t) { fJetType = t; }

  void SetMinJetTrackPt(Double_t min) { fMinJetTrackPt = min; }
  void SetMaxJetTrackPt(Double_t max) { fMaxJetTrackPt = max; }
  void SetJetTrackEtaRange(Double_t etaMin, Double_t etaMax)
  {
    fJetTrackEtaMin = etaMin;
    fJetTrackEtaMax = etaMax;
  }
  void SetJetTrackPhiRange(Double_t phiMin, Double_t phiMax)
  {
    fJetTrackPhiMax = phiMin;
    fJetTrackPhiMax = phiMax;
  }

  void SetRecombScheme(Int_t scheme) { fRecombScheme = scheme; }
  //void SetMinJetArea(Double_t a) { fMinJetArea = a; }
  void SetGhostArea(Double_t gharea) { fGhostArea = gharea; }
  void SetMinJetTowerET(Double_t min) { mTowerEnergyTMin = min; }
  void SetJetTowerEtaRange(Double_t etaMin, Double_t etaMax)
  {
    fJetTowerEtaMin = etaMin;
    fJetTowerEtaMax = etaMax;
  }
  void SetJetTowerPhiRange(Double_t phiMin, Double_t phiMax)
  {
    fJetTowerPhiMin = phiMin;
    fJetTowerPhiMax = phiMax;
  }

  // set hadronic correction fraction and type for matched tracks to towers
  void SetHadronicCorrFrac(float frac) { mHadronicCorrFrac = frac; }

  virtual void setFRunBadlist(Int_t tmpRunBadlist){
    fRunBadlist = tmpRunBadlist;
    //0 - Hanseul's //1 - Neil's
  }

  virtual void setTowerBadList(Int_t tmpTowerBadlist){
    fTowerBadlist = tmpTowerBadlist;
    //0 - Hanseul's //1 - Neil's
  }
  
  virtual void setTowerCalibrEnergy(Bool_t tmpSetTowerCalibrEnergy){
    fSetTowerCalibrEnergy = tmpSetTowerCalibrEnergy;
  }

  virtual void setMaxDcaZHadronCorr(Bool_t tmpSetDcaZHadronCorr, Double_t tmpDcaZHadronCorr){
    fSetDcaZHadronCorr = tmpSetDcaZHadronCorr;
    fDcaZHadronCorr = tmpDcaZHadronCorr;
  }

  virtual void setJetNHardestSkipped(Int_t tmpJetNHardestSkipped_010, Int_t tmpJetNHardestSkipped_1080) {
    fJetNHardestSkipped_010 = tmpJetNHardestSkipped_010;
    fJetNHardestSkipped_1080 = tmpJetNHardestSkipped_1080;
  }
  
  virtual void setJetFixedSeed(Bool_t tmpSetJetFixedSeed, Int_t tmpJetFJSeed){
  fSetJetFixedSeed = tmpSetJetFixedSeed;
  fJetFJSeed = tmpJetFJSeed;
  }

  virtual void SetEventZVtxRange(Double_t zmi, Double_t zma)
  {
    fEventZVtxMinCut = zmi;
    fEventZVtxMaxCut = zma;
  }

  //// Tree Variables

  // Fixed size dimensions of array or collections stored in the TTree if any.
  static const Int_t kMaxEvent = 1;
  static const Int_t kMaxTrack = 1000;
  static const Int_t kMaxEmcTrigger = 89;
  static const Int_t kMaxMtdTrigger = 1;
  static const Int_t kMaxBTowHit = 4800;

  static const Int_t kMaxMcVertex = 6000;
  static const Int_t kMaxMcTrack = 6000;


  TTree *fMCPico;
  // Declaration of leaf types
  Int_t Event_;
  Int_t Event_mRunId[kMaxEvent];   //[Event_]
  Int_t Event_mEventId[kMaxEvent]; //[Event_]
  Float_t Event_mPrimaryVertexX[kMaxEvent];
  Float_t Event_mPrimaryVertexY[kMaxEvent];
  Float_t Event_mPrimaryVertexZ[kMaxEvent];

  Int_t Track_;
  Float_t Track_mGMomentumX[kMaxTrack]; //[Track_]
  Float_t Track_mGMomentumY[kMaxTrack]; //[Track_]
  Float_t Track_mGMomentumZ[kMaxTrack]; //[Track_]
  Float_t Track_mOriginX[kMaxTrack];    //[Track_]
  Float_t Track_mOriginY[kMaxTrack];    //[Track_]
  Float_t Track_mOriginZ[kMaxTrack];    //[Track_]

  Char_t Track_mNHitsFit[kMaxTrack];  //[Track_]
  UChar_t Track_mNHitsMax[kMaxTrack]; //[Track_]

  Short_t Track_mNSigmaPion[kMaxTrack];            //[Track_]
  Short_t Track_mNSigmaKaon[kMaxTrack];            //[Track_]
  Short_t Track_mNSigmaProton[kMaxTrack];          //[Track_]
  Short_t Track_mNSigmaElectron[kMaxTrack];        //[Track_]
  Short_t Track_mBEmcMatchedTowerIndex[kMaxTrack]; //[Track_]
  UShort_t Track_mIdTruth[kMaxTrack];              //[Track_]
  UShort_t Track_mQATruth[kMaxTrack];              //[Track_]

  Int_t BTowHit_;
  Short_t BTowHit_mE[kMaxBTowHit]; //[BTowHit_]

  Int_t McVertex_;
  Int_t McVertex_mId[kMaxMcVertex];             //[McVertex_]
  UShort_t McVertex_mNoDaughters[kMaxMcVertex]; //[McVertex_]
  Int_t McVertex_mIdParTrk[kMaxMcVertex];       //[McVertex_]
  Int_t McVertex_mIsInterm[kMaxMcVertex];       //[McVertex_]
  Float_t McVertex_mTime[kMaxMcVertex];         //[McVertex_]
  Float_t McVertex_mVx[kMaxMcVertex];           //[McVertex_]
  Float_t McVertex_mVy[kMaxMcVertex];           //[McVertex_]
  Float_t McVertex_mVz[kMaxMcVertex];           //[McVertex_]
  Int_t McTrack_;
  UShort_t McTrack_mId[kMaxMcTrack];         //[McTrack_]
  Int_t McTrack_mGePid[kMaxMcTrack];         //[McTrack_]
  Char_t McTrack_mCharge[kMaxMcTrack];       //[McTrack_]
  UChar_t McTrack_mHits[kMaxMcTrack][22];    //[McTrack_]
  Float_t McTrack_mPx[kMaxMcTrack];          //[McTrack_]
  Float_t McTrack_mPy[kMaxMcTrack];          //[McTrack_]
  Float_t McTrack_mPz[kMaxMcTrack];          //[McTrack_]
  Float_t McTrack_mE[kMaxMcTrack];           //[McTrack_]
  Bool_t McTrack_mIsFromShower[kMaxMcTrack]; //[McTrack_]
  Short_t McTrack_mIdVtxStart[kMaxMcTrack];  //[McTrack_]
  Short_t McTrack_mIdVtxStop[kMaxMcTrack];   //[McTrack_]
  Short_t McTrack_mIdVtxItrmd[kMaxMcTrack];  //[McTrack_]

  // List of branches
  TBranch *b_Event_;                //!
  TBranch *b_Event_mRunId;          //!
  TBranch *b_Event_mEventId;        //!
  TBranch *b_Event_mPrimaryVertexX; //!
  TBranch *b_Event_mPrimaryVertexY; //!
  TBranch *b_Event_mPrimaryVertexZ; //!

  TBranch *b_Track_;            //!
  TBranch *b_Track_mGMomentumX; //!
  TBranch *b_Track_mGMomentumY; //!
  TBranch *b_Track_mGMomentumZ; //!
  TBranch *b_Track_mOriginX;    //!
  TBranch *b_Track_mOriginY;    //!
  TBranch *b_Track_mOriginZ;    //!

  TBranch *b_Track_mNHitsFit;              //!
  TBranch *b_Track_mNHitsMax;              //!
  TBranch *b_Track_mNSigmaPion;            //!
  TBranch *b_Track_mNSigmaKaon;            //!
  TBranch *b_Track_mNSigmaProton;          //!
  TBranch *b_Track_mNSigmaElectron;        //!
  TBranch *b_Track_mTopologyMap;           //!
  TBranch *b_Track_mBEmcMatchedTowerIndex; //!
  TBranch *b_Track_mIdTruth;               //!
  TBranch *b_Track_mQATruth;               //!

  TBranch *b_BTowHit_;   //!
  TBranch *b_BTowHit_mE; //!

  TBranch *b_McVertex_;             //!
  TBranch *b_McVertex_mId;          //!
  TBranch *b_McVertex_mNoDaughters; //!
  TBranch *b_McVertex_mIdParTrk;    //!
  TBranch *b_McVertex_mIsInterm;    //!
  TBranch *b_McVertex_mTime;        //!
  TBranch *b_McVertex_mVx;          //!
  TBranch *b_McVertex_mVy;          //!
  TBranch *b_McVertex_mVz;          //!
  TBranch *b_McTrack_;              //!
  TBranch *b_McTrack_mId;           //!
  TBranch *b_McTrack_mGePid;        //!
  TBranch *b_McTrack_mCharge;       //!
  TBranch *b_McTrack_mHits;         //!
  TBranch *b_McTrack_mPx;           //!
  TBranch *b_McTrack_mPy;           //!
  TBranch *b_McTrack_mPz;           //!
  TBranch *b_McTrack_mE;            //!
  TBranch *b_McTrack_mIsFromShower; //!
  TBranch *b_McTrack_mIdVtxStart;   //!
  TBranch *b_McTrack_mIdVtxStop;    //!
  TBranch *b_McTrack_mIdVtxItrmd;   //!

  void ReadTreeMc();
  void FillTree(const Int_t &numberOfD0Events);
  void FillJet(StJet *jet, StJetTreeStruct &jetTree, const TVector3 & D0);
  void OutputTreeInit();
  
    virtual Double_t GetTowerCalibEnergy(Int_t TowerId);
    StEmcADCtoEMaker *mADCtoEMaker;
    StBemcTables     *mTables;

  void 		GetAllTracksFromVertex(const Int_t &vertexid, vector<Int_t> &trackvec);
  StJet 	*DoesItHaveAGoodD0Jet(vector<TLorentzVector> &eventTracks, Int_t recoLevel);
  void 		PrepareSetOfRecoInput(const Int_t & counterEvent, const Int_t & iD0);
  Int_t 	GetMatchedRecoTrackFromMCTrack(const Int_t & mctrkid);
  Bool_t 	IsAcceptedTrack(StPicoTrack *trk, const Float_t &B, const TVector3 &Vert);
  Bool_t 	IsAcceptedTrackAndPt(StPicoTrack *trk, const Float_t &B, const TVector3 &Vert);
  Bool_t 	IsAcceptedTower(StPicoBTowHit *tower, const Int_t &towerID);
  Int_t 	EP_IsGoodTrack(StPicoTrack *trk, TVector3 pVertex);
  void 		CalculateEventPlane();
  Bool_t 	AcceptTrack(StPicoTrack *trk, Float_t B, TVector3 Vert);
  TVector3 	FastSimMom(TVector3 p, Int_t pid);
  Bool_t 	GetMomentum(TVector3 &mom, const StPicoBTowHit *tower, Double_t mass, Int_t towerID, Double_t CorrectedEnergy) const;
  Int_t 	GetMatchedBtowID(StPicoTrack *trk);
  Bool_t 	KeepTrack(const Int_t & particleid, const Int_t & centralitybin, const Double_t &  pt);
  Int_t 	fPrintLevel;


  // event selection types

 // Int_t fEmcTriggerArr[8];     // EMCal triggers array: used to select signal and do QA
  Int_t 			fCentBin;

  // variables
  Int_t 			fRunNumber;
  Int_t 			centralitybinforefficiency;
  Double_t 			fRhoVal;
  vector<Int_t> 		vertexids;
  vector<Int_t> 		pionids;
  vector<Int_t>		 	kaonids;
  vector<Int_t> 		matchedpionids;
  vector<Int_t> 		matchedkaonids;
  map<int, vector<Int_t>> 	fVertexToTracks;
  
  // Array of vectors which saves the MC track IDs matched to each vertex (START) // This can be private as I am not calling this function outside this class
  //  The dropped MC tracks will be different for each D0 because we only want to drop final state tracks which are from the KPi from D0 decaying.
  vector<Int_t> 		fDroppedMCTracks; // Dropped MC Tracks (This will be tracks which came from the KPi from D0 decaying. We don't want them.)


  Int_t 			fJetTracknHitsFit;      // requirement for track hits
  Double_t 			fJetTracknHitsRatio; // requirement for nHitsFit / nHitsMax
  bool				 fMassiveAll;
  Double_t 			fgAlpha1;
  Double_t 			fgAlpha2;
  Bool_t 			fPhiBgModulation;

  Double_t 			fJetTowerEtaMin; // min jet tower eta cut
  Double_t 			fJetTowerEtaMax; // max jet tower eta cut
  Double_t 			fJetTowerPhiMin; // min jet tower phi cut
  Double_t 			fJetTowerPhiMax; // max jet tower phi cut
  Double_t 			mTowerEnergyTMin; // min jet tower energy cut

  bool fTrackingEfficiency;               // track reconstruction efficiency
  Double_t fTrackingEfficiencyPercentage; // percentage of tracks to keep

  Double_t fMaxTowerEtBeforeHC;
  Double_t fMaxTowerEtAfterHC;
  Float_t mHadronicCorrFrac; // hadronic correction fraction from 0.0 to 1.0
  Double_t mTowerMatchTrkIndex[4800][7];
  Int_t mTowerStatusArr[4800];

 
  TFile *f;
  TFile *fout;

  // This is where I will store MC event by event information. All the processing will be done in this class, and once finished, we will have no access to the picodsts we got the files from.
  // All the track selection cuts for jets are done here.
  static const Int_t kMaxNumberOfD0Events = 100;
  vector<TLorentzVector> fMcEventTracks[kMaxNumberOfD0Events];         // For charged particles
  vector<TLorentzVector> fMcEventTowers[kMaxNumberOfD0Events];         // For neutral particles
  vector<TLorentzVector> fRecoMcEventTracks[kMaxNumberOfD0Events];     // For charged tracks
  vector<TLorentzVector> fRecoMcEventTowers[kMaxNumberOfD0Events];     // For neutral towers
  pair<TVector3, TVector3> fMcD0Information[kMaxNumberOfD0Events];     // Pion momenta, Kaon momenta (3 components each)
  pair<TVector3, TVector3> fMcRecoD0Information[kMaxNumberOfD0Events]; // Pion momenta, Kaon momenta (3 components each)
  TVector3 fOrigin[kMaxNumberOfD0Events];                              // MC Event Origin Information
  pair<int, int> fMcEventInfo[kMaxNumberOfD0Events];                   // RunID, EventID
  vector<TLorentzVector> fRecoTracks[kMaxNumberOfD0Events];
  vector<TLorentzVector> fRecoTowers[kMaxNumberOfD0Events];
  StJet *fMcJet[kMaxNumberOfD0Events];
  StJet *fMcRecoJet[kMaxNumberOfD0Events];
  StJet *fRecoJet[kMaxNumberOfD0Events];
  StJet *fRecoJet_test[kMaxNumberOfD0Events][5][5];
  //alpha11_t
  
  vector<fastjet::PseudoJet> All_fjinput; //! jet input vectors
  vector<fastjet::PseudoJet> fFull_Event; //! jet input vectors

  TTree *outputTree;

// Event histograms:
TH1D* hVtxZ;
TH2D* hVtxR;
TH1D* hVzDiff;
TH1D* hCentrality;
TH1D* hCentralityW;

//Mc event
TH2D* hPureMcNeutralEtaPhi;
TH2D* hMcJetConstMom;
TH2D* hMcJetConstTheta;
//TH2D* hMcJetConstMom5GeV;

//Jet constituents:
TH2D* hJetTracksDedx;
TH2D* hJetTracksDedxAfterCuts;
TH1D* hJetTracksPt;
TH2D* hJetTracksEtaPhi;
TH1D* hJetTracksNHitsFit;
TH1D* hJetTracksNHitsRatio;
TH1D* hJetTracksDCA;
TH1D* hJetNeutralPt;
TH2D* hJetNeutralEtaPhi;
TH2D* hJetNeutralEtBefAftHC;
TH2D* hJetNeutralECalibBefAft;
TH1D* hJetConstCharge;
TH2D* hJetConstRapPhi;
TH2D* hJetConstRapPhiICS;
TH2D* hJetConstEtaPhi;
TH2D* hJetConstEtaPhiICS;
TH1D* hJetConstPt;
TH1D* hJetConstPtICS;

TH1D* hJetMcRecoTracksPt;
TH2D* hJetMcRecoTracksEtaPhi;
TH1D* hJetMcRecoTracksNHitsFit;
TH1D* hJetMcRecoTracksNHitsRatio;
TH1D* hJetMcRecoTracksDCA;

TH1D* hJetMcRecoNeutralPt;
TH2D* hJetMcRecoNeutralEtaPhi;
TH2D* hJetMcRecoNeutralEtBefAftHC;

TH1D* hJetMcRecoConstCharge;
TEfficiency* effJetRec00_10;
TEfficiency* effJetRec10_40;
TEfficiency* effJetRec40_80;

TH2D* etaRecoTrue00_10;
TH2D* etaRecoTrue10_40;
TH2D* etaRecoTrue40_80;
//Neutral particles hadronic correction:
TH1D* hJetHadrCorrNHitsFit;
TH1D* hJetHadrCorrNHitsRatio;
TH1D* hJetHadrCorrDcaZ;
TH2D* hJetHadrCorrEtaVsPt;
TH1D* hJetHadrCorrE;


  StJetTreeStruct fMcJetTree;
  StJetTreeStruct fRecoJetTree;
  StJetTreeStruct fMcRecoJetTree;

TH2D *hResponseJetPt;
TH2D *hResponseJetD0Z;
TH2D *hResponseJetNConst;
TH2D *hResponseJetLambda1_0_5;
TH2D *hResponseJetLambda1_1;
TH2D *hResponseJetLambda1_1_5;
TH2D *hResponseJetLambda1_2;
TH2D *hResponseJetLambda1_3;
TH2D *hResponseJetMomDisp;
TH2D *hResponseJetD0DeltaR;
TH2D *hResponseJetD0Pt;



TH1D *hFractionNeutralToJetPt;

  // position objection
  StEmcGeom *mBemcGeom;

  TString fMCFileListName;

  vector<TString> filenamesforHIOverlay;

  TF1 *fKaonMomResolution;
  TF1 *fPionMomResolution;
  TF1 *fProtonMomResolution;

  TGraph *fPionWeight[3];
  TGraph *fKaonWeight[3];
  TGraph *fProtonWeight[3];
  TGraph *fAProtonWeight[3];

  Int_t fJetAlgo;      // jet algorithm (kt, akt, etc)
  Int_t fJetType;      // jet type (full, charged, neutral)
  Int_t fRecombScheme; // recombination scheme used by fastjet

  StFJWrapper *fjw;     //! fastjet wrapper
                        // jet attributes
 
  ClassDef(StHIOverlayAngularities, 1)
};
#endif
