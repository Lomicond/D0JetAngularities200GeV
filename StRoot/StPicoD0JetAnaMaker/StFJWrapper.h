#ifndef StFJWrapper_H
#define StFJWrapper_H

#if !defined(__CINT__)

// ROOT includes
#include <TMath.h>
#include <TList.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <vector>
#include <TString.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TF1.h"

// jet includes
#include "FJ_includes.h"
////#include "StJetShape.h"
////#include "StJetPicoDefinitions.h"
#include "StFJAngularityDefinition.h"

using namespace std;

class StFJWrapper
{
 public:
  StFJWrapper(const char *name, const char *title);
  
  virtual ~StFJWrapper();

  virtual void  AddInputVector (Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index = -99999);
  virtual void  AddInputVector (const fastjet::PseudoJet& vec,                Int_t index = -99999);
  virtual void  AddInputVectors(const std::vector<fastjet::PseudoJet>& vecs,  Int_t offsetIndex = -99999);
  virtual void  AddInputGhost  (Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index = -99999);
  virtual const char *ClassName()                            const { return "StFJWrapper";              }
  virtual void  Clear(const Option_t* /*opt*/ = "");
  virtual void  ClearMemory();
  virtual void  CopySettingsFrom (const StFJWrapper& wrapper);
  virtual void  GetMedianAndSigma(Double_t& median, Double_t& sigma, Int_t remove = 0) const;
  fastjet::ClusterSequenceArea*           GetClusterSequence() const   { return fClustSeq;                 }
  fastjet::ClusterSequence*               GetClusterSequenceSA() const { return fClustSeqSA;               }
  fastjet::ClusterSequenceActiveAreaExplicitGhosts* GetClusterSequenceGhosts() const { return fClustSeqActGhosts; }
  const std::vector<fastjet::PseudoJet>&  GetInputVectors()    const { return fInputVectors;               }
  const std::vector<fastjet::PseudoJet>&  GetInputGhosts()     const { return fInputGhosts;                }
  const std::vector<fastjet::PseudoJet>&  GetInclusiveJets()   const { return fInclusiveJets;              }
  Double_t  				  GetJetRho()   	const { return fJetRho;              }
  Double_t  				  GetJetRhoM()   	const { return fJetRhoM;              }
  const std::vector<fastjet::PseudoJet>&  GetFilteredJets()    const { return fFilteredJets;               }
  std::vector<fastjet::PseudoJet>         GetJetConstituents(UInt_t idx) const;
  std::vector<fastjet::PseudoJet>         GetFilteredJetConstituents(UInt_t idx) const;
  Double_t                                GetMedianUsedForBgSubtraction() const { return fMedUsedForBgSub; }
  const char*                             GetName()            const { return fName;                       }
  const char*                             GetTitle()           const { return fTitle;                      }
  Double_t                                GetJetArea         (UInt_t idx) const;
  fastjet::PseudoJet                      GetJetAreaVector   (UInt_t idx) const;
  Double_t                                GetFilteredJetArea (UInt_t idx) const;
  fastjet::PseudoJet                      GetFilteredJetAreaVector(UInt_t idx) const;
  Double_t                                GetJetSubtractedPt (UInt_t idx) const;
  virtual std::vector<double>             GetSubtractedJetsPts(Double_t median_pt = -1, Bool_t sorted = kFALSE);
  Bool_t                                  GetLegacyMode()            { return fLegacyMode; }
  Bool_t                                  GetDoFilterArea()          { return fDoFilterArea; }
  Double_t                                NSubjettiness(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Int_t Option=0);
  Double32_t                              NSubjettinessDerivativeSub(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Double_t JetR, fastjet::PseudoJet jet, Int_t Option=0);
#ifdef FASTJET_VERSION

  const std::vector<fastjet::PseudoJet>                      GetConstituentSubtrJets()            const { return fConstituentSubtrJets           ; }
  const std::vector<fastjet::PseudoJet>                      GetGroomedJets()                     const { return fGroomedJets                    ; }
////  Int_t CreateGenSub();          // fastjet::contrib::GenericSubtractor
  Int_t CreateConstituentSub();  // fastjet::contrib::ConstituentSubtractor
  Int_t CreateSoftDrop();
#endif
  virtual std::vector<double>                                GetGRNumerator()                     const { return fGRNumerator                    ; }
  virtual std::vector<double>                                GetGRDenominator()                   const { return fGRDenominator                  ; }
  virtual std::vector<double>                                GetGRNumeratorSub()                  const { return fGRNumeratorSub                 ; }
  virtual std::vector<double>                                GetGRDenominatorSub()                const { return fGRDenominatorSub               ; }

  virtual void RemoveLastInputVector();

  virtual Int_t Run();
  virtual Int_t Run_Shape();
  virtual Int_t Run_AreaBase();
  virtual Int_t Filter();
  virtual Int_t DoConstituentSubtraction();
  virtual Int_t DoSoftDrop();
  
  void SetName(const char* name)        { fName           = name;    }
  void SetTitle(const char* title)      { fTitle          = title;   }
  void SetStrategy(const fastjet::Strategy &strat)                 { fStrategy = strat;  }
  void SetAlgorithm(const fastjet::JetAlgorithm &algor)            { fAlgor    = algor;  }
  void SetRecombScheme(const fastjet::RecombinationScheme &scheme) { fScheme   = scheme; }
  void SetAreaType(const fastjet::AreaType &atype)                 { fAreaType = atype;  }
  void SetNRepeats(Int_t nrepeat)       { fNGhostRepeats  = nrepeat; }
  void SetGhostArea(Double_t gharea)    { fGhostArea      = gharea;  }
  void SetMaxRap(Double_t maxrap)       { fMaxRap         = maxrap;  }
  void SetR(Double_t r)                 { fR              = r;       }
  void SetGridScatter(Double_t gridSc)  { fGridScatter    = gridSc;  }
  void SetKtScatter(Double_t ktSc)      { fKtScatter      = ktSc;    }
  void SetMeanGhostKt(Double_t meankt)  { fMeanGhostKt    = meankt;  }
  void SetPluginAlgor(Int_t plugin)     { fPluginAlgor    = plugin;  }
  void SetUseArea4Vector(Bool_t useA4v) { fUseArea4Vector = useA4v;  }
  void SetCentrality(Int_t centrality)  { fCentrality = centrality; }
  void SetCentralityW(Int_t centrality)  { fCentralityWeight = centrality; }
  void SetRMax1(double rMax1) {frMax1 = rMax1;}
  void SetRMax2(double rMax2) {frMax2 = rMax2;}
  void SetNIter(double NIter) {fNIter = NIter;}
  void SetAlpha1(double alpha1) {fAlpha1 = alpha1;}
  void SetAlpha2(double alpha2) {fAlpha2 = alpha2;}
  void SetMassiveTest(bool MassiveTest) {fMassiveTest = MassiveTest;}
  void SetEventPlane2(Double_t EP_psi2) {fEP_psi2 = EP_psi2;}
  //fBackSub
  void SetBackgroundSub(Bool_t BackSub)  { fBackSub = BackSub; }
  void SetSubtractionMc(Bool_t SubtractMc)  { fSubtractMc = SubtractMc; }
  Bool_t GetBackgroundSub() const { return fBackSub; }
  Double_t GetActualShapeL10half() const { return fAngul10half; }
  Double_t GetActualShapeL11() const { return fAngul11; }
  Double_t GetActualShapeL11half() const { return fAngul11half; }
  Double_t GetActualShapeL12() const { return fAngul12; }
  Double_t GetActualShapeL13() const { return fAngul13; }
  Double_t GetActualShapeDisp() const { return fAngulDisp; }


  void SetPhiModulation(Bool_t PhiModulation) {fPhiModulation = PhiModulation;}
  void setJetNHardestSkipped(Int_t tmpJetNHardestSkipped_010, Int_t tmpJetNHardestSkipped_1080) {
  	fJetNHardestSkipped_010 = tmpJetNHardestSkipped_010;
  	fJetNHardestSkipped_1080 = tmpJetNHardestSkipped_1080;
  	}
  void setJetFixedSeed(Bool_t tmpSetJetFixedSeed, Int_t tmpJetFJSeed){
        fSetJetFixedSeed = tmpSetJetFixedSeed;
	fJetFJSeed = tmpJetFJSeed;
  }
  void SetupAlgorithmfromOpt(const char *option);
  void SetupAreaTypefromOpt(const char *option);
  void SetupSchemefromOpt(const char *option);
  void SetupStrategyfromOpt(const char *option);
  void SetLegacyMode (Bool_t mode)      { fLegacyMode ^= mode; }
  void SetLegacyFJ();
  void SetUseExternalBkg(Bool_t b, Double_t rho, Double_t rhom) { fUseExternalBkg = b; fRho = rho; fRhom = rhom;}
  void SetRMaxAndStep(Double_t rmax, Double_t dr) {fRMax = rmax; fDRStep = dr; }
  void SetRhoRhom (Double_t rho, Double_t rhom) { fUseExternalBkg = kTRUE; fRho = rho; fRhom = rhom;} // if using rho,rhom then fUseExternalBkg is true
  
  void SetHJetConstRapPhiICS(TH2D* hist) { hJetConstRapPhiICS = hist; }
  void SetHJetConstEtaPhiICS(TH2D* hist2) { hJetConstEtaPhiICS = hist2; }
  void SetHJetConstPtICS(TH1D* hist3) {hJetConstPtICS = hist3;}


  void PrintInput();
  // void PrintJetDescription();

 protected:
  Double_t fAngul10half;
  Double_t fAngul11;
  Double_t fAngul11half;
  Double_t fAngul12;
  Double_t fAngul13;
  Double_t fAngul1half;
  Double_t fAngulDisp; 
 
  Double_t fJetRho;
  Double_t fJetRhoM;
  TString                                fName;               //!
  TString                                fTitle;              //!
  std::vector<fastjet::PseudoJet>        fInputVectors;       //!
  std::vector<fastjet::PseudoJet>        fInputGhosts;        //!
  std::vector<fastjet::PseudoJet>        fInclusiveJets;      //!
  std::vector<fastjet::PseudoJet>        fFilteredJets;       //!
  std::vector<double>                    fSubtractedJetsPt;   //!
  std::vector<fastjet::PseudoJet>        fConstituentSubtrJets; //!
  std::vector<fastjet::PseudoJet>        fGroomedJets;        //!
  fastjet::AreaDefinition               *fAreaDef;            //!
  fastjet::VoronoiAreaSpec              *fVorAreaSpec;        //!
  fastjet::GhostedAreaSpec              *fGhostedAreaSpec;    //!
  fastjet::JetDefinition                *fJetDef;             //!
  fastjet::JetDefinition::Plugin        *fPlugin;             //!
#ifndef FASTJET_VERSION
  fastjet::RangeDefinition              *fRange;              //!
#else
  fastjet::Selector                     *fRange;              //!
#endif
  fastjet::ClusterSequenceArea          *fClustSeq;           //!
  fastjet::ClusterSequence              *fClustSeqSA;                //!
  fastjet::ClusterSequenceActiveAreaExplicitGhosts *fClustSeqActGhosts; //!
  fastjet::Strategy                      fStrategy;           //!
  fastjet::JetAlgorithm                  fAlgor;              //!
  fastjet::RecombinationScheme           fScheme;             //!
  fastjet::AreaType                      fAreaType;           //!
  Int_t                                  fNGhostRepeats;      //!
  Double_t                               fGhostArea;	      //!
  Double_t                               fMaxRap;	      //!
  Double_t                               fR;                  //!
  
  Double_t frMax1;
  Double_t frMax2;
  Int_t fNIter;
  Double_t fAlpha1;
  Double_t fAlpha2;
  bool fMassiveTest;
  Double_t fEP_psi2;
  
  TH2D* hJetConstRapPhiICS;
  TH2D* hJetConstEtaPhiICS;
  TH1D* hJetConstPtICS;
  
  // no setters for the moment - used default values in the constructor
  Double_t                               fGridScatter;        //!
  Double_t                               fKtScatter;	      //!
  Double_t                               fMeanGhostKt;        //!
  Int_t                                  fPluginAlgor;        //!
  Int_t                                  fCentrality;        //!
  Double_t				 fCentralityWeight;
  // extra parameters
  Double_t                               fMedUsedForBgSub;    //!
  Bool_t                                 fUseArea4Vector;     //!
  Bool_t				 fBackSub;
  Bool_t				 fSubtractMc;
  Bool_t 				 fPhiModulation;
  Int_t  		    		 fJetNHardestSkipped_010;
  Int_t    				 fJetNHardestSkipped_1080;
  Bool_t                  		 fSetJetFixedSeed;
  Int_t                                  fJetFJSeed;
  // condition to stop the grooming (rejection of soft splitting) z > fZcut theta^fBeta
  Double_t                               fZcut;               // fZcut = 0.1                
  Double_t                               fBeta;               // fBeta = 0
#ifdef FASTJET_VERSION
  fastjet::JetMedianBackgroundEstimator   *fBkrdEstimator;    //!
  //from contrib package
  fastjet::contrib::ConstituentSubtractor *fConstituentSubtractor;    //!
  fastjet::contrib::SoftDrop              *fSoftDrop;        //!

#endif
  Bool_t                                   fDoFilterArea;         //!
  Bool_t                                   fLegacyMode;           //!
  Bool_t                                   fUseExternalBkg;       //!
  Double_t                                 fRho;                  //  pT background density
  Double_t                                 fRhom;                 //  mT background density
  Double_t                                 fRMax;             //!
  Double_t                                 fDRStep;           //!
  std::vector<double>                      fGRNumerator;      //!
  std::vector<double>                      fGRDenominator;    //!
  std::vector<double>                      fGRNumeratorSub;   //!
  std::vector<double>                      fGRDenominatorSub; //!

  virtual void   SubtractBackground(const Double_t median_pt = -1);

 private:
  StFJWrapper();
  StFJWrapper(const StFJWrapper& wrapper);
  StFJWrapper& operator = (const StFJWrapper& wrapper);
};
#endif
#endif

#ifdef StFJWrapper_CXX
#undef StFJWrapper_CXX

#if defined __GNUC__
#pragma GCC system_header
#endif

namespace fj = fastjet;

//_________________________________________________________________________________________________
StFJWrapper::StFJWrapper(const char *name, const char *title)
  :
    fName              (name)
  , fTitle             (title)
  , fInputVectors      ( )
  , fAngul10half       (-999)
  , fAngul11           (-999)
  , fAngul11half       (-999)
  , fAngul12           (-999)
  , fAngul13           (-999)
  , fAngul1half        (-999)
  , fAngulDisp         (-999)
  , fInputGhosts       ( )
  , fInclusiveJets     ( )
  , fJetRho            (-999)
  , fJetRhoM           (-999)
  , fFilteredJets      ( )
  , fSubtractedJetsPt  ( )
  , fConstituentSubtrJets ( )
  , fSoftDrop          ( )
  , fAreaDef           (0)
  , fVorAreaSpec       (0)
  , fGhostedAreaSpec   (0)
  , fJetDef            (0)
  , fPlugin            (0)
  , fRange             (0)
  , fClustSeq          (0)
  , fClustSeqSA        (0)
  , fClustSeqActGhosts (0)
  , fStrategy          (fj::Best)
  , fAlgor             (fj::kt_algorithm)
  , fScheme            (fj::BIpt_scheme)
  , fAreaType          (fj::active_area)
  , fNGhostRepeats     (1)
  , fGhostArea         (0.005)
  , fMaxRap            (1.)
  , fR                 (0.4)
  , fGridScatter       (1.0)
  , fKtScatter         (0.1)
  , fMeanGhostKt       (1e-100)
  , fPluginAlgor       (0)
  , fMedUsedForBgSub   (0)
  , fUseArea4Vector    (kFALSE)
  , fZcut(0.1)
  , fBeta(0)
  , fBackSub	       (kTRUE)
  , fSubtractMc        (kTRUE)
  , fPhiModulation     (kFALSE)
  , fJetNHardestSkipped_010 (2)
  , fJetNHardestSkipped_1080 (1)
  , fSetJetFixedSeed   (kFALSE)
  , fJetFJSeed         (12345)
  , frMax1(0.15)
  , frMax2(0.20)
  , fNIter(0)
  , fAlpha1(1)
  , fAlpha2(1)
  , fMassiveTest(kTRUE)
  , fEP_psi2(-999)
#ifdef FASTJET_VERSION
  , fBkrdEstimator     (0)
  , fConstituentSubtractor (0)
#endif
  , fDoFilterArea      (false)
  , fLegacyMode        (false)
  , fUseExternalBkg    (false)
  , fRho               (0)
  , fRhom              (0)
  , fRMax(2.)
  , fDRStep(0.04)
  , fGRNumerator()
  , fGRDenominator()
  , fGRNumeratorSub()
  , fGRDenominatorSub()
  , fCentrality(-999)
  , fCentralityWeight(-999)
{
  // Constructor.
}

//_________________________________________________________________________________________________
StFJWrapper::~StFJWrapper()
{
  // Destructor.
  ClearMemory();
}

//_________________________________________________________________________________________________
void StFJWrapper::ClearMemory()
{
  // Destructor.
  if (fAreaDef)           { delete fAreaDef;           fAreaDef         = NULL; }
  if (fVorAreaSpec)       { delete fVorAreaSpec;       fVorAreaSpec     = NULL; }
  if (fGhostedAreaSpec)   { delete fGhostedAreaSpec;   fGhostedAreaSpec = NULL; }
  if (fJetDef)            { delete fJetDef;            fJetDef          = NULL; }
  if (fPlugin)            { delete fPlugin;            fPlugin          = NULL; }
  if (fRange)             { delete fRange;             fRange           = NULL; }
  if (fClustSeq)          { delete fClustSeq;          fClustSeq        = NULL; }
  if (fClustSeqSA)        { delete fClustSeqSA;        fClustSeqSA        = NULL; }
  if (fClustSeqActGhosts) { delete fClustSeqActGhosts; fClustSeqActGhosts = NULL; }
  #ifdef FASTJET_VERSION
  if (fBkrdEstimator)          { delete fBkrdEstimator; fBkrdEstimator = NULL; }
////  if (fGenSubtractor)          { delete fGenSubtractor; fGenSubtractor = NULL; }
  if (fConstituentSubtractor)  { delete fConstituentSubtractor; fConstituentSubtractor = NULL; }
  if (fSoftDrop)          { delete fSoftDrop; fSoftDrop = NULL;}
  #endif
}

//_________________________________________________________________________________________________
void StFJWrapper::CopySettingsFrom(const StFJWrapper& wrapper)
{
  // Copy some settings.
  // You very often want to keep most of the settings
  // but change only the algorithm or R - do it after call to this function

  fStrategy         = wrapper.fStrategy;
  fAlgor            = wrapper.fAlgor;
  fScheme           = wrapper.fScheme;
  fAreaType         = wrapper.fAreaType;
  fNGhostRepeats    = wrapper.fNGhostRepeats;
  fGhostArea        = wrapper.fGhostArea;
  fMaxRap           = wrapper.fMaxRap;
  fR                = wrapper.fR;
  fGridScatter      = wrapper.fGridScatter;
  fKtScatter        = wrapper.fKtScatter;
  fMeanGhostKt      = wrapper.fMeanGhostKt;
  fPluginAlgor      = wrapper.fPluginAlgor;
  fUseArea4Vector   = wrapper.fUseArea4Vector;
  fZcut             = wrapper.fZcut;
  fBeta             = wrapper.fBeta;
  fLegacyMode       = wrapper.fLegacyMode;
  fUseExternalBkg   = wrapper.fUseExternalBkg;
  fRho              = wrapper.fRho;
  fRhom             = wrapper.fRhom;
  fCentrality	    = wrapper.fCentrality;
  fCentralityWeight = wrapper.fCentralityWeight;		
}

//_________________________________________________________________________________________________
void StFJWrapper::Clear(const Option_t */*opt*/)
{
  // Simply clear the input vectors.
  // Make sure done on every event if the instance is reused
  // Reset the median to zero.

  fInputVectors.clear();
  fInputGhosts.clear();
  fMedUsedForBgSub = 0;

  // for the moment brute force delete everything
  ClearMemory();
}

//_________________________________________________________________________________________________
void StFJWrapper::RemoveLastInputVector()
{
  // Remove last input vector
  fInputVectors.pop_back();
}

// void StFJWrapper::PrintJetDescription()
// {
//   cout << "Clustering with " << fJetDef->description() << endl; 
// }

void StFJWrapper::PrintInput()
{
  cout << "        pt eta phi" << endl;
  for (int i = 0; i < fInputVectors.size(); i++){
    cout << Form("Input # %i \t %.2f \t %.2f \t %.2f \t %.2f ", i, fInputVectors[i].pt(), fInputVectors[i].eta(), fInputVectors[i].phi(), fInputVectors[i].e()) << endl;
  }
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputVector(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index)
{
  // Make the input pseudojet.
  fastjet::PseudoJet inVec(px, py, pz, E);

  // Salvatore Aiola: not sure why this was done...
  //if (index > -99999) {
  inVec.set_user_index(index);
  //} else {
  //inVec.set_user_index(fInputVectors.size());
  //}

  // add to the fj container of input vectors
  fInputVectors.push_back(inVec);
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputVector(const fj::PseudoJet& vec, Int_t index)
{
  // Add an input pseudojet.
  fj::PseudoJet inVec = vec;

  // Salvatore Aiola: not sure why this was done...
  ///if (index > -99999) {
  inVec.set_user_index(index);
  //} else {
  //inVec.set_user_index(fInputVectors.size());
  //}

  // add to the fj container of input vectors
  fInputVectors.push_back(inVec);
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputVectors(const std::vector<fj::PseudoJet>& vecs, Int_t offsetIndex)
{
  // Add the input from vector of pseudojets.
  for (UInt_t i = 0; i < vecs.size(); ++i) {
    fj::PseudoJet inVec = vecs[i];
    if (offsetIndex > -99999)
      inVec.set_user_index(fInputVectors.size() + offsetIndex);
    // add to the fj container of input vectors
    fInputVectors.push_back(inVec);
  }
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputGhost(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index)
{
  // Make the input pseudojet.
  fastjet::PseudoJet inVec(px, py, pz, E);

  if (index > -99999) {
    inVec.set_user_index(index);
  } else {
    inVec.set_user_index(fInputGhosts.size());
  }

  // add to the fj container of input vectors
  fInputGhosts.push_back(inVec);
  if (!fDoFilterArea) fDoFilterArea = kTRUE;
}

//_________________________________________________________________________________________________
Double_t StFJWrapper::GetJetArea(UInt_t idx) const
{
  // Get the jet area.
  Double_t retval = -1; // really wrong area..
  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->area(fInclusiveJets[idx]);
  } else {
    //__ERROR(Form("Wrong index: %d",idx));
    cout << (Form("Wrong index: %d",idx)) << endl;
  }
  return retval;
}

//_________________________________________________________________________________________________
Double_t StFJWrapper::GetFilteredJetArea(UInt_t idx) const
{
  // Get the filtered jet area.
  Double_t retval = -1; // really wrong area..
  if (fDoFilterArea && fClustSeqActGhosts && (idx<fFilteredJets.size())) {
    retval = fClustSeqActGhosts->area(fFilteredJets[idx]);
  } else {
    //__ERROR(Form("Wrong index: %d",idx));
    cout << (Form("Wrong index: %d",idx)) << endl;
  }
  return retval;
}

//_________________________________________________________________________________________________
fastjet::PseudoJet StFJWrapper::GetJetAreaVector(UInt_t idx) const
{
  // Get the jet area as vector.
  fastjet::PseudoJet retval;
  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->area_4vector(fInclusiveJets[idx]);
  } else {
    //__ERROR(Form("Wrong index: %d",idx));
    cout << (Form("Wrong index: %d",idx)) << endl;
  }
  return retval;
}

//_________________________________________________________________________________________________
fastjet::PseudoJet StFJWrapper::GetFilteredJetAreaVector(UInt_t idx) const
{
  // Get the jet area as vector.
  fastjet::PseudoJet retval;
  if (fDoFilterArea && fClustSeqActGhosts && (idx<fFilteredJets.size())) {
    retval = fClustSeqActGhosts->area_4vector(fFilteredJets[idx]);
  } else {
   // __ERROR(Form("Wrong index: %d",idx));
    cout << (Form("Wrong index: %d",idx)) << endl;    
  }
  return retval;
}

//_________________________________________________________________________________________________
std::vector<double> StFJWrapper::GetSubtractedJetsPts(Double_t median_pt, Bool_t sorted)
{
  // Get subtracted jets pTs, returns vector.
  SubtractBackground(median_pt);

  if (kTRUE == sorted) {
    std::sort(fSubtractedJetsPt.begin(), fSubtractedJetsPt.begin());
  }
  return fSubtractedJetsPt;
}

//_________________________________________________________________________________________________
Double_t StFJWrapper::GetJetSubtractedPt(UInt_t idx) const
{
  // Get subtracted jets pTs, returns Double_t.
  Double_t retval = -99999.; // really wrong pt..
  if ( idx < fSubtractedJetsPt.size() ) {
    retval = fSubtractedJetsPt[idx];
  }
  return retval;
}

//_________________________________________________________________________________________________
std::vector<fastjet::PseudoJet>
StFJWrapper::GetJetConstituents(UInt_t idx) const
{
  // Get jets constituents.
  std::vector<fastjet::PseudoJet> retval;

  if ( idx < fInclusiveJets.size() ) {
    retval = fClustSeq->constituents(fInclusiveJets[idx]);
  } else {
   // __ERROR(Form("Wrong index: %d",idx));
    cout << (Form("Wrong index: %d",idx)) << endl;
  }

  return retval;
}

//_________________________________________________________________________________________________
std::vector<fastjet::PseudoJet>
StFJWrapper::GetFilteredJetConstituents(UInt_t idx) const
{
  // Get jets constituents.
  std::vector<fastjet::PseudoJet> retval;

  if ( idx < fFilteredJets.size() ) {
    if (fClustSeqSA)        retval = fClustSeqSA->constituents(fFilteredJets[idx]);
    if (fClustSeqActGhosts) retval = fClustSeqActGhosts->constituents(fFilteredJets[idx]);
  } else {
   // __ERROR(Form("Wrong index: %d",idx));
   cout << (Form("Wrong index: %d",idx)) << endl;
  }

  return retval;
}
///
//_________________________________________________________________________________________________
void StFJWrapper::GetMedianAndSigma(Double_t &median, Double_t &sigma, Int_t remove) const
{
  // Get the median and sigma from fastjet.
  // User can also do it on his own because the cluster sequence is exposed (via a getter)
  if (!fClustSeq) {
   // __ERROR(Form("Run the jfinder first."));
    cout << (Form("Run the jfinder first.")) << endl;
    return;
  }

  Double_t mean_area = 0;
  try {
    if(0 == remove) {
      fClustSeq->get_median_rho_and_sigma(*fRange, fUseArea4Vector, median, sigma, mean_area);
    }  else {
      std::vector<fastjet::PseudoJet> input_jets = sorted_by_pt(fClustSeq->inclusive_jets());
      input_jets.erase(input_jets.begin(), input_jets.begin() + remove);
      fClustSeq->get_median_rho_and_sigma(input_jets, *fRange, fUseArea4Vector, median, sigma, mean_area);
      input_jets.clear();
    }
  } catch (fj::Error) {
    //__WARNING(Form("FJ Exception caught."));
    cout << Form("FJ Exception caught.") << endl;
    median = -1.;
    sigma = -1;
  }
}

//_________________________________________________________________________________________________

Int_t StFJWrapper::Run_AreaBase(){

Int_t NiterA[4] = {4,3,2,2};
Double_t R_max1A[4] = {0.05,0.125,0.100,0.150};
Double_t R_max2A[4] = {0.005,0.005,0.175,0.100};
Int_t difiter = 2;



//////________________Background_estimation_______________

    //Scaling    //v2 settings
    Double_t v2 = 0.0;
    //https://journals.aps.org/prc/pdf/10.1103/PhysRevC.77.054901
    //https://www.hepdata.net/record/ins777954
    switch (fCentrality) {

    case 0: v2 = 0.0696; break; //70-80%
    case 1: v2 = 0.0723; break; //60-70%
    case 2: v2 = 0.0744; break; //50-60%
    case 3: v2 = 0.074;  break; //40-50%
    case 4: v2 = 0.0703; break; //30-40%
    case 5: v2 = 0.0618; break; //20-30%
    case 6: v2 = 0.0476; break; //10-20%
    case 7: v2 = 0.0339; break; //5-10%
    case 8: v2 = 0.0232; break; //0-5%
 
    default: v2 = 0.0; break;
}
    
    Double_t v3 = 0;
    Double_t v4 = 0;
    Double_t psi = fEP_psi2;
    
    fj::contrib::BackgroundRescalingYPhi rescaling(v2,v3,v4,psi,0.,1.,0.,1.);
    
    //y - not used, for RHIC it is particaly flat
    rescaling.use_rap_term(false);
    //phi used - but no visible impact
    rescaling.use_phi_term(true);
    
    //--- Background Estimation ---
    // Define jet algorithm for background estimation
    fj::JetDefinition jet_def_bkgd(fj::kt_algorithm, fR, fj::E_scheme, fj::Best); 
    // Define the area for background estimation
    fj::AreaDefinition area_def_bkgd(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.005));
    
    //For tuning, fix the seed for fastjet
    if (fSetJetFixedSeed) {
	
    	Int_t seed1 = fJetFJSeed;
	Int_t seed2 = fJetFJSeed;
	std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
	
	area_def_bkgd = area_def_bkgd.with_fixed_seed(seeds);
	
    }
    
    //One hardest jet by default, two for the most central //fJetNHardestSkipped_010
    Int_t nJetsRemove = fJetNHardestSkipped_1080;
    if (fCentrality == 7 || fCentrality == 8) nJetsRemove = fJetNHardestSkipped_010;
    
    //Definition of the selector for background estimation (eta and pt cut + remove the n hardest jets)
    fj::Selector selector = (!fj::SelectorNHardest(nJetsRemove)) * fj::SelectorAbsEtaMax(0.6); // * fj::SelectorPtMin(0.01); //Comp. test
    
    // Create background estimator using the previously defined selector and jet algorithm
    fj::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    
    
	std::vector<fastjet::PseudoJet> inputRealOnly;
	inputRealOnly.reserve(fInputVectors.size());

	for (const auto &pj : fInputVectors) {
	  const int uid = pj.user_index();

	  // vyhoď MC reco tracky i MC reco towery
	  if (fSubtractMc && uid >= 10000)  continue;
	  if (fSubtractMc && uid <= -10000) continue;

	  inputRealOnly.push_back(pj);
	}

    if (fBackSub){
	    //Estimation of the background using only charged tracks
	    if (fPhiModulation) bkgd_estimator.set_rescaling_class(&rescaling);
	    
	    bkgd_estimator.set_particles(inputRealOnly); 
	    fJetRho = bkgd_estimator.rho();
	    fJetRhoM = bkgd_estimator.rho_m(); 

    }
    
//////________________Jet_estimation_______________

     //Inclusive (or track based) jet definition
     fj::JetDefinition jet_def(fj::antikt_algorithm, fR, fj::E_scheme, fj::Best);
     
     //Definition of the area for jet reconstruction
     fj::AreaDefinition area_def(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.005)); //Comp. test
     
    //For tuning, fix the seed for fastjet
    if (fSetJetFixedSeed) {
	
    	Int_t seed1 = fJetFJSeed;
	Int_t seed2 = fJetFJSeed;
	std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
	
	area_def = area_def.with_fixed_seed(seeds);
	
    }
    
     //Definition of the clustering
     //fj::ClusterSequenceArea clust_seq_hard(fInputVectors, jet_def, area_def);
     
     fClustSeq = new fastjet::ClusterSequenceArea(fInputVectors, jet_def, area_def);
	
     // Sort the jets by transverse momentum (pt)
     Double_t ptmin = 0;
     fInclusiveJets.clear();
     fInclusiveJets = sorted_by_pt(fClustSeq->inclusive_jets(ptmin));
    
     return 0;

}

Int_t StFJWrapper::Run()
{
//----------------------------------------------------------------


Int_t NiterA[4] = {4,3,2,2};
Double_t R_max1A[4] = {0.05,0.125,0.100,0.150};
Double_t R_max2A[4] = {0.005,0.005,0.175,0.100};
Int_t difiter = 2;
/*****/
   //Scaling    //v2 settings
    Double_t v2 = 0.0;
//https://journals.aps.org/prc/pdf/10.1103/PhysRevC.77.054901
//https://www.hepdata.net/record/ins777954
    switch (fCentrality) {

    case 0: v2 = 0.0696; break; //70-80%
    case 1: v2 = 0.0723; break; //60-70%
    case 2: v2 = 0.0744; break; //50-60%
    case 3: v2 = 0.074;  break; //40-50%
    case 4: v2 = 0.0703; break; //30-40%
    case 5: v2 = 0.0618; break; //20-30%
    case 6: v2 = 0.0476; break; //10-20%
    case 7: v2 = 0.0339; break; //5-10%
    case 8: v2 = 0.0232; break; //0-5%
 
    default: v2 = 0.0; break;
}
    
    Double_t v3 = 0;
    Double_t v4 = 0;
    Double_t psi = fEP_psi2;

    //BackgroundRescalingYPhiUsingVectorForY(Double_t v2, Double_t v3, Double_t v4, Double_t psi, std::vector<Double_t> values, std::vector<Double_t> rap_binning);
    //BackgroundRescalingYPhi(): _v2(0), _v3(0), _v4(0), _psi(0), _a1(1), _sigma1(1000), _a2(0), _sigma2(1000), _use_rap(false), _use_phi(false) {}
    fj::contrib::BackgroundRescalingYPhi rescaling(v2,v3,v4,psi,0.,1.,0.,1.);
    rescaling.use_rap_term(false);    // this is useful to check if the vectors with rapidity dependence have good sizes, if one wants to use also the rapidity rescaling.
    rescaling.use_phi_term(true);

	//--- Background Estimation ---
	fj::JetDefinition jet_def_bkgd(fj::kt_algorithm, fR, fj::E_scheme, fj::Best); // Define jet algorithm for background estimation
	fj::AreaDefinition area_def_bkgd(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.005)); // Define the area for background estimation
	
	//For tuning, fix the seed for fastjet
	if (fSetJetFixedSeed) {
		
	    	Int_t seed1 = fJetFJSeed;
		Int_t seed2 = fJetFJSeed;
		std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
		
		area_def_bkgd = area_def_bkgd.with_fixed_seed(seeds);
		
	}
	
	// Decide how many jets to remove based on centrality
        Int_t nJetsRemove = fJetNHardestSkipped_1080;
        if (fCentrality == 7 || fCentrality == 8) nJetsRemove = fJetNHardestSkipped_010;

	// Selector to choose the hardest jets based on centrality and eta
	fj::Selector selector = (!fj::SelectorNHardest(nJetsRemove)) * fj::SelectorAbsEtaMax(0.6);

	// Create background estimator using the previously defined selector and jet algorithm
	fj::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);

	// Set the maximum eta value for background estimation
	Double_t max_eta = 1;
        
	//--- Jet Reconstruction ---
	fj::contrib::IterativeConstituentSubtractor subtractor; // Set up the background subtraction algorithm
	subtractor.set_distance_type(fj::contrib::ConstituentSubtractor::deltaR); // Set distance type for subtraction

  	// Set parameters for the background subtraction (maximum distance and alpha values)
	vector<Double_t> max_distances;
	vector<Double_t> alphas;
	

	
	//For higher iterations always used the same R_1max and alpha
	if (difiter < 0 || difiter >= 4) {
	    std::cerr << "ERROR: difiter out of range! difiter = " << difiter << std::endl;
	    exit(1);
	}
	

	if(NiterA[difiter] > 2){
	
		for (Int_t it = 0; it < NiterA[difiter] ; it++){
			max_distances.push_back(R_max1A[difiter]);
			alphas.push_back(0);
		}
	} else //else two iterations
	{
		max_distances.push_back(R_max1A[difiter]);
		max_distances.push_back(R_max2A[difiter]);
		alphas.push_back(0);
		alphas.push_back(0);
	}
	

	
	if (fBackSub&&false){
		cout << NiterA[difiter] << "test" << endl;
		for (auto& a : alphas) {
			cout << "alpha: " << a << endl;
		}
		for (auto& rr : max_distances) {	
			cout << "max_distances: " << rr << endl;
		}
	}
	
	//Exclude MC for Background estimation
	std::vector<fastjet::PseudoJet> inputRealOnly;
	inputRealOnly.reserve(fInputVectors.size());

	for (const auto &pj : fInputVectors) {
	  const int uid = pj.user_index();

	  // exclude MC reco tracks, MC D0, and MC reco towers
	  if (fSubtractMc &&uid >= 10000)  continue;
	  if (fSubtractMc &&uid <= -10000) continue;

	  inputRealOnly.push_back(pj);
	}
	
	//exclude D0
	fj::Selector notD02 = fastjet::SelectorMassMax(1); // 1 GeV max mass
	std::vector<fastjet::PseudoJet> particles_without_D0;
	std::vector<fastjet::PseudoJet> D0_pseudojet;

	for (auto& p : fInputVectors) {
	
	    if (notD02(p) && p.user_index()<30000) {
	    	fastjet::PseudoJet temp_jet; 

	    	if (fMassiveTest) temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E());
	    	else temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()));
	    	
		temp_jet.set_user_index(p.user_index());
		
	    	particles_without_D0.push_back(temp_jet);
	    } else {
		    fastjet::PseudoJet temp_jet; 
		    
		    if (fMassiveTest) temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), p.E());
		    else temp_jet = fastjet::PseudoJet(p.px(), p.py(), p.pz(), sqrt(p.px()*p.px()+p.py()*p.py()+p.pz()*p.pz()));
		    
		    temp_jet.set_user_index(p.user_index());
		    
		    D0_pseudojet.push_back(temp_jet);
	    }
	}
	///////////////////////////////////
	
	// Apply the subtraction parameters
	subtractor.set_parameters(max_distances, alphas);
	subtractor.set_ghost_removal(true); // Enable ghost removal
	subtractor.set_ghost_area(0.005); // Set ghost area value
	subtractor.set_max_eta(max_eta); // Set maximum eta for particles
	
	
	subtractor.set_background_estimator(&bkgd_estimator); // Link the background estimator
	
	if(fMassiveTest) subtractor.set_common_bge_for_rho_and_rhom(true); // Set common background estimation for rho and rhom
	if(fMassiveTest) subtractor.set_keep_original_masses(); // Keep the original masses of particles
	subtractor.set_scale_fourmomentum();
	
        
       	// Selector for particles with mass below 1 GeV (likely excludes certain particles)
	fj::Selector notD0 = fastjet::SelectorMassMax(1); // 1 GeV max mass
	subtractor.set_particle_selector(&notD0);  

	subtractor.initialize(); // Initialize the background subtraction algorithm

	// Define jet algorithm for actual jet reconstruction (Anti-kt with radius fR)
	fj::JetDefinition jet_def(fj::antikt_algorithm, fR, fj::E_scheme, fj::Best);

	// Define area for jet reconstruction
	fj::AreaDefinition area_def_jet(fj::active_area, fj::GhostedAreaSpec(1.2, 1, 0.005));
	
	//For tuning, fix the seed for fastjet
	if (fSetJetFixedSeed) {
		
	    	Int_t seed1 = fJetFJSeed;
		Int_t seed2 = fJetFJSeed;
		std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
		
		area_def_jet = area_def_jet.with_fixed_seed(seeds);
		
	}

	// Set particles for background estimator (background estimation for the input particles)
	////bkgd_estimator.set_particles(fInputVectors); //all_vectors
	//Rescaling
	if (fPhiModulation) bkgd_estimator.set_rescaling_class(&rescaling);
	bkgd_estimator.set_particles(inputRealOnly); //MC cant be included

	
	// Minimum pt for jets (below this, jets are excluded)
	Double_t ptmin = 0;
	vector<fj::PseudoJet> corrected_event;

	// Apply background subtraction if enabled
	if (fBackSub) corrected_event = subtractor.subtract_event(particles_without_D0); 
	else corrected_event = particles_without_D0; // Otherwise, use original input vectors
	
	
	
	if (fBackSub){
	
        for (vector<fastjet::PseudoJet>::const_iterator particle = corrected_event.begin(); particle != corrected_event.end(); ++particle) {
		 hJetConstRapPhiICS->Fill(particle->phi_std(), particle->rap(),fCentralityWeight);
		 hJetConstEtaPhiICS->Fill(particle->phi_std(), particle->eta(),fCentralityWeight);
		 hJetConstPtICS->Fill(particle->perp(),fCentralityWeight);		
	}
	
	}
	
	// Return back D0
	//corrected_event.push_back(D0_pseudojet.back());
	if (!D0_pseudojet.empty()) corrected_event.push_back(D0_pseudojet.back());
	else cout << "D0 missing!!!" << endl;

	
	// Perform jet clustering with the background-subtracted (or original) particles
	if (fClustSeq) {
	    	delete fClustSeq;
	    	fClustSeq = nullptr;
	}
	fClustSeq = new fastjet::ClusterSequenceArea(corrected_event, jet_def, area_def_jet);
	
	// Sort the jets by transverse momentum (pt)
	//vector<fastjet::PseudoJet> fInclusiveJets;
	fInclusiveJets.clear();
	fInclusiveJets = sorted_by_pt(fClustSeq->inclusive_jets(ptmin));

  return 0;
}
Int_t StFJWrapper::Run_Shape()
{

//----------------------------------------------------------------
        

    //Scaling    //v2 settings
    double v2 = 0.0;
//https://journals.aps.org/prc/pdf/10.1103/PhysRevC.77.054901
//https://www.hepdata.net/record/ins777954
    switch (fCentrality) {

    case 0: v2 = 0.0696; break; //70-80%
    case 1: v2 = 0.0723; break; //60-70%
    case 2: v2 = 0.0744; break; //50-60%
    case 3: v2 = 0.074;  break; //40-50%
    case 4: v2 = 0.0703; break; //30-40%
    case 5: v2 = 0.0618; break; //20-30%
    case 6: v2 = 0.0476; break; //10-20%
    case 7: v2 = 0.0339; break; //5-10%
    case 8: v2 = 0.0232; break; //0-5%
 
    default: v2 = 0.0; break;
}
    
    double v3 = 0;
    double v4 = 0;
    double psi = fEP_psi2;

    //BackgroundRescalingYPhiUsingVectorForY(double v2, double v3, double v4, double psi, std::vector<double> values, std::vector<double> rap_binning);
    //BackgroundRescalingYPhi(): _v2(0), _v3(0), _v4(0), _psi(0), _a1(1), _sigma1(1000), _a2(0), _sigma2(1000), _use_rap(false), _use_phi(false) {}
    fj::contrib::BackgroundRescalingYPhi rescaling(v2,v3,v4,psi,0.,1.,0.,1.);
    rescaling.use_rap_term(false);    // this is useful to check if the vectors with rapidity dependence have good sizes, if one wants to use also the rapidity rescaling.
    rescaling.use_phi_term(true);
    
	

	//-----------------------------------------------------------
	
	// Define jet algorithm for actual jet reconstruction (Anti-kt with radius fR)
	fj::JetDefinition jet_def(fj::antikt_algorithm, fR, fj::E_scheme, fj::Best);

	// Define area for jet reconstruction
	fj::AreaDefinition area_def_jet(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.005));
	
	//For tuning, fix the seed for fastjet
	if (fSetJetFixedSeed) {
		
	    	Int_t seed1 = fJetFJSeed;
		Int_t seed2 = fJetFJSeed;
		std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
		
		area_def_jet = area_def_jet.with_fixed_seed(seeds);
	}
          
    	std::vector<fastjet::PseudoJet> inputRealOnly;
	inputRealOnly.reserve(fInputVectors.size());
	
	for (const auto &pj : fInputVectors) {
		  const int uid = pj.user_index();

		  if (fSubtractMc &&uid >= 10000)  continue;
		  if (fSubtractMc &&uid <= -10000) continue;

		  inputRealOnly.push_back(pj);
	}
	
	// Perform jet clustering with the background-subtracted (or original) particles
	fClustSeq = new fj::ClusterSequenceArea(fInputVectors, jet_def, area_def_jet);
	
	// Sort the jets by transverse momentum (pt)
	double ptmin = 0;
	fInclusiveJets.clear();
	fInclusiveJets = fClustSeq->inclusive_jets();
	
	//----------------------------------------------------------
	//--- Background Estimation ---
	fj::JetDefinition jet_def_bkgd(fj::kt_algorithm, fR, fj::E_scheme, fj::Best); // Define jet algorithm for background estimation
	fj::AreaDefinition area_def_bkgd(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.005)); // Define the area for background estimation
	
	//For tuning, fix the seed for fastjet
	if (fSetJetFixedSeed) {
		
	    	Int_t seed1 = fJetFJSeed;
		Int_t seed2 = fJetFJSeed;
		std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
		
		area_def_bkgd = area_def_bkgd.with_fixed_seed(seeds);
		
	}
	
	// Decide how many jets to remove based on centrality
        Int_t nJetsRemove = fJetNHardestSkipped_1080;
        if (fCentrality == 7 || fCentrality == 8) nJetsRemove = fJetNHardestSkipped_010;

	// Selector to choose the hardest jets based on centrality and eta
	fj::Selector selector = (!fj::SelectorNHardest(nJetsRemove)) * fj::SelectorAbsEtaMax(0.6);// * fj::SelectorPtMin(0.01);

	// Create background estimator using the previously defined selector and jet algorithm
	fj::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
	
	//Estimation of the background using only charged tracks
	if (fPhiModulation) bkgd_estimator.set_rescaling_class(&rescaling);
	bkgd_estimator.set_particles(inputRealOnly);
	
	if (fBackSub){
	   fJetRho = bkgd_estimator.rho();
	   fJetRhoM = bkgd_estimator.rho_m(); 
	}
	//--------------------------------------------------------------
	
	fastjet::contrib::GenericSubtractor gensub(&bkgd_estimator);
	gensub.set_common_bge_for_rho_and_rhom(true);
	fastjet::contrib::GenericSubtractorInfo info;



	double jet_R = fR;
	//kappa,alfa,jet_R
	Angularity my_angularity_10half(1,0.5,jet_R);
	Angularity my_angularity_11(1,1,jet_R);
	Angularity my_angularity_11half(1,1.5,jet_R);
	Angularity my_angularity_12(1,2,jet_R);
	Angularity my_angularity_13(1,3,jet_R);
	Angularity my_angularity_Disp(2,0,jet_R);
	
	fastjet::FunctionOfPseudoJet<double>* shape10half = &my_angularity_10half;
	fastjet::FunctionOfPseudoJet<double>* shape11 = &my_angularity_11;
	fastjet::FunctionOfPseudoJet<double>* shape11half = &my_angularity_11half;
	fastjet::FunctionOfPseudoJet<double>* shape12 = &my_angularity_12;
	fastjet::FunctionOfPseudoJet<double>* shape13 = &my_angularity_13;
	fastjet::FunctionOfPseudoJet<double>* shapeDisp = &my_angularity_Disp;
				
for (unsigned int i = 0; i < fInclusiveJets.size(); i++) {

    const fastjet::PseudoJet &jet = fInclusiveJets[i];

    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    int d0_count = 0;
    for (size_t j = 0; j < constituents.size(); j++) {
    
        // Looking for D0
        if (constituents[j].m() > 1.0 || constituents[j].user_index()>=30000) {
            d0_count++;
        }
    }

    // Pokud je v jetu právě jeden D⁰ meson, provedeme odečet angularity
    if (d0_count == 1) {
        fAngul10half = gensub(*shape10half, jet, info);
        fAngul11 = gensub(*shape11, jet, info);
        fAngul11half = gensub(*shape11half, jet, info);
        fAngul12 = gensub(*shape12, jet, info);
        fAngul13 = gensub(*shape13, jet, info);
        fAngulDisp = sqrt(gensub(*shapeDisp, jet, info));
    } 
    
}






  return 0;
}
//_________________________________________________________________________________________________
Int_t StFJWrapper::Filter()
{
//
//  StFJWrapper::Filter
//
  fJetDef = new fj::JetDefinition(fAlgor, fR, fScheme, fStrategy);

  if (fDoFilterArea) {
    if (fInputGhosts.size()>0) {
      try {
        fClustSeqActGhosts = new fj::ClusterSequenceActiveAreaExplicitGhosts(fInputVectors,
                                                                           *fJetDef,
                                                                            fInputGhosts,
                                                                            fGhostArea);
      } catch (fj::Error) {
       // __WARNING(Form("FJ Exception caught."));
            cout << Form("FJ Exception caught.") << endl;
        return -1;
      }

      fFilteredJets.clear();
      fFilteredJets =  fClustSeqActGhosts->inclusive_jets(0.0);
    } else {
      return -1;
    }
  } else {
    try {
      fClustSeqSA = new fastjet::ClusterSequence(fInputVectors, *fJetDef);
    } catch (fj::Error) {
     // __WARNING(Form("FJ Exception caught."));
          cout << Form("FJ Exception caught.") << endl;
      return -1;
    }

    fFilteredJets.clear();
    fFilteredJets = fClustSeqSA->inclusive_jets(0.0);
  }

  return 0;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetLegacyFJ()
{
  // This methods enable legacy behaviour (FastJet 2.x) when StROOT is compiled with FastJet 3.x
#ifdef FASTJET_VERSION
    std::cout << "WARNING! Setting FastJet in legacy mode" << std::endl;
    if (fGhostedAreaSpec) { fGhostedAreaSpec->set_fj2_placement(kTRUE); }
     if (fBkrdEstimator) {
      fBkrdEstimator->set_provide_fj2_sigma(kTRUE);
      fBkrdEstimator->set_use_area_4vector(kFALSE);
    }
#endif
}

//_________________________________________________________________________________________________
void StFJWrapper::SubtractBackground(Double_t median_pt)
{
  // Subtract the background (specify the value - see below the meaning).
  // Negative argument means the bg will be determined with the current algorithm
  // this is the default behavior. Zero is allowed
  // Note: user may set the switch for area4vector based subtraction.

  Double_t median    = 0;
  Double_t sigma     = 0;
  Double_t mean_area = 0;

  // clear the subtracted jet pt's vector<double>
  fSubtractedJetsPt.clear();

  // check what was specified (default is -1)
  if (median_pt < 0) {
    try {
      fClustSeq->get_median_rho_and_sigma(*fRange, fUseArea4Vector, median, sigma, mean_area);
    }

    catch (fj::Error) {
     // __WARNING(Form("FJ Exception caught."));
          cout << Form("FJ Exception caught.") << endl;
      median = -9999.;
      sigma = -1;
      fMedUsedForBgSub = median;
      return;
    }
    fMedUsedForBgSub = median;
  } else {
    // we do not know the sigma in this case
    sigma = -1;
    if (0.0 == median_pt) {
    //  __WARNING(Form("Median specified for bg subtraction is ZERO: %f \n", median_pt ));
        cout << Form("Median specified for bg subtraction is ZERO: %f \n", median_pt ) << endl;
      fMedUsedForBgSub = 0.;
    } else {
      fMedUsedForBgSub = median_pt;
    }
  }

  // subtract:
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    if ( fUseArea4Vector ) {
      // subtract the background using the area4vector
      fj::PseudoJet area4v = fClustSeq->area_4vector(fInclusiveJets[i]);
      fj::PseudoJet jet_sub = fInclusiveJets[i] - area4v * fMedUsedForBgSub;
      fSubtractedJetsPt.push_back(jet_sub.perp()); // here we put only the pt of the jet - note: this can be negative
    } else {
      // subtract the background using scalars
      // fj::PseudoJet jet_sub = fInclusiveJets[i] - area * fMedUsedForBgSub_;
      Double_t area = fClustSeq->area(fInclusiveJets[i]);
      // standard subtraction
      Double_t pt_sub = fInclusiveJets[i].perp() - fMedUsedForBgSub * area;
      fSubtractedJetsPt.push_back(pt_sub); // here we put only the pt of the jet - note: this can be negative
    }
  }
}

//_________________________________________________________________________________________________
Int_t StFJWrapper::DoConstituentSubtraction() {
  // Do constituent subtraction
#ifdef FASTJET_VERSION
  CreateConstituentSub();
  // fConstituentSubtractor->set_alpha(/* double alpha */);
  // fConstituentSubtractor->set_max_deltaR(/* double max_deltaR */);

  // clear constituent subtracted jets
  fConstituentSubtrJets.clear();
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::PseudoJet subtracted_jet(0.,0.,0.,0.);
    if(fInclusiveJets[i].perp()>0.)
      subtracted_jet = (*fConstituentSubtractor)(fInclusiveJets[i]);
    fConstituentSubtrJets.push_back(subtracted_jet);
  }
  if(fConstituentSubtractor) { delete fConstituentSubtractor; fConstituentSubtractor = NULL; }

#endif
  return 0;
}

//_________________________________________________________________________________________________
Int_t StFJWrapper::DoSoftDrop() {
  // Do grooming
#ifdef FASTJET_VERSION
  CreateSoftDrop();

  // clear groomed jets
  fGroomedJets.clear();
  //fastjet::Subtractor fjsub (fBkrdEstimator);
  //fSoftDrop->set_subtractor(&fjsub);
  //fSoftDrop->set_input_jet_is_subtracted(false); //??
  
  for (unsigned i = 0; i < fInclusiveJets.size(); i++) {
    fj::PseudoJet groomed_jet(0.,0.,0.,0.);
    if(fInclusiveJets[i].perp()>0.){
      groomed_jet = (*fSoftDrop)(fInclusiveJets[i]);
      if(groomed_jet!=0) fGroomedJets.push_back(groomed_jet);
    }
    
  }
  if(fSoftDrop) { delete fSoftDrop; fSoftDrop = NULL; }

#endif
  return 0;
}

//_________________________________________________________________________________________________
Int_t StFJWrapper::CreateSoftDrop() {
  // Do grooming
  #ifdef FASTJET_VERSION
  if (fSoftDrop) { delete fSoftDrop; } // protect against memory leaks
  
  fSoftDrop   = new fj::contrib::SoftDrop(fBeta,fZcut);
  
  #endif
  return 0;
}

/*
//_________________________________________________________________________________________________
Int_t StFJWrapper::CreateGenSub() {
  // Do generic subtraction for jet mass
  #ifdef FASTJET_VERSION
  if (fGenSubtractor) { delete fGenSubtractor; } // protect against memory leaks

  if (fUseExternalBkg)
    { fGenSubtractor   = new fj::contrib::GenericSubtractor(fRho,fRhom); }
  else
    {
    fGenSubtractor     = new fj::contrib::GenericSubtractor(fBkrdEstimator);
    #if FASTJET_VERSION_NUMBER >= 30100
    fGenSubtractor->set_common_bge_for_rho_and_rhom(); // see contrib 1.020 GenericSubtractor.hh line 62
    #endif
    }

  #endif
  return 0;
}
*/

//_________________________________________________________________________________________________
Int_t StFJWrapper::CreateConstituentSub() {
  // Do generic subtraction for jet mass
  #ifdef FASTJET_VERSION
  if (fConstituentSubtractor) { delete fConstituentSubtractor; } // protect against memory leaks

  // see ConstituentSubtractor.hh signatures
  // ConstituentSubtractor(double rho, double rhom=0, double alpha=0, double maxDeltaR=-1)
  if (fUseExternalBkg) { fConstituentSubtractor = new fj::contrib::ConstituentSubtractor(fRho,fRhom); }
  else                 { fConstituentSubtractor = new fj::contrib::ConstituentSubtractor(fBkrdEstimator); }  // FIXME Nov15, 2018 commented back in

  #endif
  return 0;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetupAlgorithmfromOpt(const char *option)
{
  // Setup algorithm from char.
  std::string opt(option);

  if (!opt.compare("kt"))                fAlgor    = fj::kt_algorithm;
  if (!opt.compare("antikt"))            fAlgor    = fj::antikt_algorithm;
  if (!opt.compare("cambridge"))         fAlgor    = fj::cambridge_algorithm;
  if (!opt.compare("genkt"))             fAlgor    = fj::genkt_algorithm;
  if (!opt.compare("cambridge_passive")) fAlgor    = fj::cambridge_for_passive_algorithm;
  if (!opt.compare("genkt_passive"))     fAlgor    = fj::genkt_for_passive_algorithm;
  if (!opt.compare("ee_kt"))             fAlgor    = fj::ee_kt_algorithm;
  if (!opt.compare("ee_genkt"))          fAlgor    = fj::ee_genkt_algorithm;
  if (!opt.compare("plugin"))            fAlgor    = fj::plugin_algorithm;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetupAreaTypefromOpt(const char *option)
{
  // Setup area type from char.
  std::string opt(option);

  if (!opt.compare("active"))                      fAreaType = fj::active_area;
  if (!opt.compare("invalid"))                     fAreaType = fj::invalid_area;
  if (!opt.compare("active_area_explicit_ghosts")) fAreaType = fj::active_area_explicit_ghosts;
  if (!opt.compare("one_ghost_passive"))           fAreaType = fj::one_ghost_passive_area;
  if (!opt.compare("passive"))                     fAreaType = fj::passive_area;
  if (!opt.compare("voronoi"))                     fAreaType = fj::voronoi_area;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetupSchemefromOpt(const char *option)
{
  //
  // setup scheme from char
  //
  std::string opt(option);

  if (!opt.compare("BIpt"))   fScheme   = fj::BIpt_scheme;
  if (!opt.compare("BIpt2"))  fScheme   = fj::BIpt2_scheme;
  if (!opt.compare("E"))      fScheme   = fj::E_scheme;
  if (!opt.compare("pt"))     fScheme   = fj::pt_scheme;
  if (!opt.compare("pt2"))    fScheme   = fj::pt2_scheme;
  if (!opt.compare("Et"))     fScheme   = fj::Et_scheme;
  if (!opt.compare("Et2"))    fScheme   = fj::Et2_scheme;
}

//_________________________________________________________________________________________________
void StFJWrapper::SetupStrategyfromOpt(const char *option)
{
  // Setup strategy from char.
  std::string opt(option);

  if (!opt.compare("Best"))            fStrategy = fj::Best;
  if (!opt.compare("N2MinHeapTiled"))  fStrategy = fj::N2MinHeapTiled;
  if (!opt.compare("N2Tiled"))         fStrategy = fj::N2Tiled;
  if (!opt.compare("N2PoorTiled"))     fStrategy = fj::N2PoorTiled;
  if (!opt.compare("N2Plain"))         fStrategy = fj::N2Plain;
  if (!opt.compare("N3Dumb"))          fStrategy = fj::N3Dumb;
  if (!opt.compare("NlnN"))            fStrategy = fj::NlnN;
  if (!opt.compare("NlnN3pi"))         fStrategy = fj::NlnN3pi;
  if (!opt.compare("NlnN4pi"))         fStrategy = fj::NlnN4pi;
  if (!opt.compare("NlnNCam4pi"))      fStrategy = fj::NlnNCam4pi;
  if (!opt.compare("NlnNCam2pi2R"))    fStrategy = fj::NlnNCam2pi2R;
  if (!opt.compare("NlnNCam"))         fStrategy = fj::NlnNCam;
  if (!opt.compare("plugin"))          fStrategy = fj::plugin_strategy;
}

//_______________________________________________________________________________________________ 
Double_t StFJWrapper::NSubjettiness(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Int_t Option){
  //Option 0=Nsubjettiness result, 1=opening angle between axes in Eta-Phi plane, 2=Distance between axes in Eta-Phi plane
  
  fJetDef = new fj::JetDefinition(fAlgor, fR*100, fScheme, fStrategy ); //the *2 is becasue of a handful of jets that end up missing a track for some reason.

  try {
    fClustSeqSA = new fastjet::ClusterSequence(fInputVectors, *fJetDef);
    // ClustSeqSA = new fastjet::ClusterSequenceArea(fInputVectors, *fJetDef, *fAreaDef);
  } catch (fj::Error) {
    //__WARNING(Form("FJ Exception caught."));
    cout << Form("FJ Exception caught.") << endl;
    return -1;
  }
  fFilteredJets.clear();
  fFilteredJets = fClustSeqSA->inclusive_jets(0); //becasue this is < not <=
  Double_t Result=-1;
  std::vector<fastjet::PseudoJet> SubJet_Axes;
  fj::PseudoJet SubJet1_Axis;
  fj::PseudoJet SubJet2_Axis;
  if (Algorithm==0){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::KT_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==1) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::CA_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==2){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::AntiKT_Axes(Radius), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==3) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::WTA_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==4) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::WTA_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==5){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==6){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==7){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_AntiKT_Axes(Radius), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==8){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_WTA_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==9){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_WTA_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==10){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::MultiPass_Axes(100), fj::contrib::NormalizedMeasure(Beta,fR));
    Result= nSub.result(fFilteredJets[0]);
    SubJet_Axes=nSub.currentAxes();
  }

  SubJet1_Axis=SubJet_Axes[0];	
  Double_t SubJet1_Eta=SubJet1_Axis.pseudorapidity();
  Double_t SubJet2_Eta;
  Double_t SubJet1_Phi=SubJet1_Axis.phi();
  if(SubJet1_Phi < -1*TMath::Pi()) SubJet1_Phi += (2*TMath::Pi());
  else if (SubJet1_Phi > TMath::Pi()) SubJet1_Phi -= (2*TMath::Pi());
  Double_t SubJet2_Phi;
  Double_t DeltaPhi=-5;
  if (SubJet_Axes.size()>1){
    SubJet2_Axis=SubJet_Axes[1];
    SubJet2_Eta=SubJet2_Axis.pseudorapidity();
    SubJet2_Phi=SubJet2_Axis.phi();
    if(SubJet2_Phi < -1*TMath::Pi()) SubJet2_Phi += (2*TMath::Pi());
    else if (SubJet2_Phi > TMath::Pi()) SubJet2_Phi -= (2*TMath::Pi());
    DeltaPhi=SubJet1_Phi-SubJet2_Phi;
    if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
    else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  }
    
  if (Option==0) return Result;
  else if (Option==1 && SubJet_Axes.size()>1 && N==2) return TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2));
  else if (Option==2 && SubJet_Axes.size()>1 && N==2) return TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2));
  else return -2;
}

//_______________________________________________________________________________________________
Double32_t StFJWrapper::NSubjettinessDerivativeSub(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Double_t JetR, fastjet::PseudoJet jet, Int_t Option){ //For derivative subtraction

  Double_t Result=-1;
  std::vector<fastjet::PseudoJet> SubJet_Axes;
  fj::PseudoJet SubJet1_Axis;
  fj::PseudoJet SubJet2_Axis;
  if (Algorithm==0){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::KT_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==1) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::CA_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==2){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::AntiKT_Axes(Radius), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==3) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::WTA_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==4) {
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::WTA_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==5){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==6){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==7){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_AntiKT_Axes(Radius), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==8){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_WTA_KT_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==9){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::OnePass_WTA_CA_Axes(), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }
  else if (Algorithm==10){
    fj::contrib::Nsubjettiness nSub(N, fj::contrib::MultiPass_Axes(100), fj::contrib::NormalizedMeasure(Beta,JetR));
    Result= nSub.result(jet);
    SubJet_Axes=nSub.currentAxes();
  }

  SubJet1_Axis=SubJet_Axes[0];	
  Double_t SubJet1_Eta=SubJet1_Axis.pseudorapidity();
  Double_t SubJet2_Eta;
  Double_t SubJet1_Phi=SubJet1_Axis.phi();
  if(SubJet1_Phi < -1*TMath::Pi()) SubJet1_Phi += (2*TMath::Pi());
  else if (SubJet1_Phi > TMath::Pi()) SubJet1_Phi -= (2*TMath::Pi());
  Double_t SubJet2_Phi;
  Double_t DeltaPhi=-5;
  if (SubJet_Axes.size()>1){
    SubJet2_Axis=SubJet_Axes[1];
    SubJet2_Eta=SubJet2_Axis.pseudorapidity();
    SubJet2_Phi=SubJet2_Axis.phi();
    if(SubJet2_Phi < -1*TMath::Pi()) SubJet2_Phi += (2*TMath::Pi());
    else if (SubJet2_Phi > TMath::Pi()) SubJet2_Phi -= (2*TMath::Pi());
    DeltaPhi=SubJet1_Phi-SubJet2_Phi;
    if(DeltaPhi < -1*TMath::Pi()) DeltaPhi += (2*TMath::Pi());
    else if (DeltaPhi > TMath::Pi()) DeltaPhi -= (2*TMath::Pi());
  }
    
  if (Option==0) return Result;
  else if (Option==1 && SubJet_Axes.size()>1 && N==2) return TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2));
  else if (Option==2 && SubJet_Axes.size()>1 && N==2) return TMath::Sqrt(TMath::Power(SubJet1_Eta-SubJet2_Eta,2)+TMath::Power(DeltaPhi,2));
  else return -2;

}
#endif
