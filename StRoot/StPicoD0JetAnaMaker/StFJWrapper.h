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
#include "StFJAngularityDefinition.h"

using namespace std;
#include "fastjet/FunctionOfPseudoJet.hh"

class ConstantRescaling : public fastjet::FunctionOfPseudoJet<double> {
public:
  explicit ConstantRescaling(double C) : fC(C) {}
  double result(const fastjet::PseudoJet&) const override { return fC; }
private:
  double fC;
};

#include "fastjet/tools/BackgroundEstimatorBase.hh"
#include <sstream>

class ScaledBGE : public fastjet::BackgroundEstimatorBase {
public:
  // wrapper, NEvlastní podkladový estimator (jen ukazatel)
  ScaledBGE(fastjet::BackgroundEstimatorBase *bge, double scale)
  : _bge(bge), _scale(scale) {}

  // povinné virtuální metody (pure virtual) --------------------------

  // FastJet může chtít kopii background estimatoru (např. interně v subtractorech)
  fastjet::BackgroundEstimatorBase * copy() const override {
    // udělej kopii podkladového estimátoru a zabal ji do nové ScaledBGE
    fastjet::BackgroundEstimatorBase *bge_copy = _bge->copy();
    return new ScaledBGE(bge_copy, _scale, /*own=*/true);
  }

  // nastavení nové události
  void set_particles(const std::vector<fastjet::PseudoJet> &particles) override {
    _bge->set_particles(particles);
  }

  void set_particles_with_seed(const std::vector<fastjet::PseudoJet> &particles,
                               const std::vector<int> &seed) override {
    _bge->set_particles_with_seed(particles, seed);
  }

  // plný odhad pozadí (globální)
  fastjet::BackgroundEstimate estimate() const override {
    fastjet::BackgroundEstimate be = _bge->estimate();
    be.apply_rescaling_factor(_scale);   // <- rho, sigma, rho_m, sigma_m
    return be;
  }

  // plný odhad pozadí v pozici referenčního jetu
  fastjet::BackgroundEstimate estimate(const fastjet::PseudoJet &jet) const override {
    fastjet::BackgroundEstimate be = _bge->estimate(jet);
    be.apply_rescaling_factor(_scale);
    return be;
  }

  // convenience API, které některé části FastJetu používají přímo
  double rho() const override { return _scale * _bge->rho(); }

  double rho(const fastjet::PseudoJet &jet) override {
    return _scale * _bge->rho(jet);
  }

  // pokud chceš (doporučuju), přepošli i sigma/rho_m, když jsou podporované
  double sigma() const override {
    return _scale * _bge->sigma();
  }

  double sigma(const fastjet::PseudoJet &jet) override {
    return _scale * _bge->sigma(jet);
  }

  bool has_sigma() const override { return _bge->has_sigma(); }

  double rho_m() const override { return _scale * _bge->rho_m(); }
  double rho_m(const fastjet::PseudoJet &jet) override { return _scale * _bge->rho_m(jet); }

  double sigma_m() const override { return _scale * _bge->sigma_m(); }
  double sigma_m(const fastjet::PseudoJet &jet) override { return _scale * _bge->sigma_m(jet); }

  bool has_rho_m() const override { return _bge->has_rho_m(); }

  std::string description() const override {
    std::ostringstream ss;
    ss << "ScaledBGE(scale=" << _scale << ") wrapping: " << _bge->description();
    return ss.str();
  }

  // -------------------------- lifecycle ----------------------------
  ~ScaledBGE() override {
    if (_own && _bge) { delete _bge; _bge = nullptr; }
  }

private:
  // konstruktor pro copy(), kde ScaledBGE vlastní kopii podkladového estimátoru
  ScaledBGE(fastjet::BackgroundEstimatorBase *bge, double scale, bool own)
  : _bge(bge), _scale(scale), _own(own) {}

  fastjet::BackgroundEstimatorBase *_bge = nullptr;
  double _scale = 1.0;
  bool _own = false;
};

class StFJWrapper
{
 public:
  StFJWrapper(const char *name, const char *title);
  
  virtual ~StFJWrapper();

  virtual void                            AddInputVector (Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index = -99999);
  virtual void                            AddInputVector (const fastjet::PseudoJet& vec,                Int_t index = -99999);
  virtual void                            AddInputVectors(const std::vector<fastjet::PseudoJet>& vecs);
  virtual const char *                    ClassName() const { return "StFJWrapper";}
  virtual void                            Clear(const Option_t* /*opt*/ = "");
  virtual void                            ClearMemory();
  fastjet::ClusterSequenceArea*           GetClusterSequence() const   { return fClustSeq;                 }
  fastjet::ClusterSequence*               GetClusterSequenceSA() const { return fClustSeqSA;               }
  fastjet::ClusterSequenceActiveAreaExplicitGhosts* GetClusterSequenceGhosts() const { return fClustSeqActGhosts; }
  const std::vector<fastjet::PseudoJet>&  GetInputVectors()    const { return fInputVectors;               }
  const std::vector<fastjet::PseudoJet>&  GetInputGhosts()     const { return fInputGhosts;                }
  const std::vector<fastjet::PseudoJet>&  GetInclusiveJets()   const { return fInclusiveJets;              }
  Double_t  				                      GetJetRho()   	const { return fJetRho;              }
  Double_t  				                      GetJetRhoM()   	const { return fJetRhoM;              }
  const std::vector<fastjet::PseudoJet>&  GetFilteredJets()    const { return fFilteredJets;               }
  std::vector<fastjet::PseudoJet>         GetJetConstituents(UInt_t idx) const;
  const char*                             GetName()            const { return fName;                       }
  const char*                             GetTitle()           const { return fTitle;                      }
  Double_t                                GetJetArea         (UInt_t idx) const;
  fastjet::PseudoJet                      GetJetAreaVector   (UInt_t idx) const;
  Bool_t                                  GetLegacyMode()            { return fLegacyMode; }
  Double_t                                NSubjettiness(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Int_t Option=0);
  Double32_t                              NSubjettinessDerivativeSub(Int_t N, Int_t Algorithm, Double_t Radius, Double_t Beta, Double_t JetR, fastjet::PseudoJet jet, Int_t Option=0);

  virtual std::vector<Double_t>                                GetGRNumerator()                     const { return fGRNumerator                    ; }
  virtual std::vector<Double_t>                                GetGRDenominator()                   const { return fGRDenominator                  ; }
  virtual std::vector<Double_t>                                GetGRNumeratorSub()                  const { return fGRNumeratorSub                 ; }
  virtual std::vector<Double_t>                                GetGRDenominatorSub()                const { return fGRDenominatorSub               ; }

  virtual void RemoveLastInputVector();

  virtual Int_t runICSMethod();
  virtual Int_t runAreaBgAndAreaShapeMethod();
  virtual Int_t Filter();
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
  //void SetMassiveTest(Bool_t MassiveTest) {fMassiveTest = MassiveTest;}
  void SetEventPlane2(Double_t EP_psi2) {fEP_psi2 = EP_psi2;}
  void SetBackgroundSub(Bool_t BackSub)  { fBackSub = BackSub; }
  void SetSubtractionMc(Bool_t SubtractMc)  { fSubtractMc = SubtractMc; }
  Bool_t GetBackgroundSub() const { return fBackSub; }
  Double_t GetActualShapeL10half() const { return fAngul10half; }
  Double_t GetActualShapeL11() const { return fAngul11; }
  Double_t GetActualShapeL11half() const { return fAngul11half; }
  Double_t GetActualShapeL12() const { return fAngul12; }
  Double_t GetActualShapeL13() const { return fAngul13; }
  Double_t GetActualShapeDisp() const { return fAngulDisp; }
  Double_t ComputeOccupancyC_fromJetsUsed(const fastjet::JetMedianBackgroundEstimator& bge);
  void SetICSSubtractionParams(const std::vector<Double_t>& maxDistances, const std::vector<Double_t>& alphas) {
    fICSMaxDistances = maxDistances;
    fICSAlphas = alphas;
  }
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

 protected:
  Double_t                               fAngul10half;
  Double_t                               fAngul11;
  Double_t                               fAngul11half;
  Double_t                               fAngul12;
  Double_t                               fAngul13;
  Double_t                               fAngul1half;
  Double_t                               fAngulDisp; 
  Double_t                               fJetRho;
  Double_t                               fJetRhoM;
  TString                                fName;               //!
  TString                                fTitle;              //!
  std::vector<fastjet::PseudoJet>        fInputVectors;       //!
  std::vector<fastjet::PseudoJet>        fInputGhosts;        //!
  std::vector<fastjet::PseudoJet>        fInclusiveJets;      //!
  std::vector<fastjet::PseudoJet>        fFilteredJets;       //!
  std::vector<Double_t>                  fSubtractedJetsPt;   //!
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
  Double_t                               fEP_psi2;
  
  TH2D*                                  hJetConstRapPhiICS;
  TH2D*                                  hJetConstEtaPhiICS;
  TH1D*                                  hJetConstPtICS;
  
  // no setters for the moment - used default values in the constructor
  Double_t                               fGridScatter;       
  Double_t                               fKtScatter;	      
  Double_t                               fMeanGhostKt;       
  Int_t                                  fPluginAlgor;        
  Int_t                                  fCentrality;       
  Double_t				                       fCentralityWeight;
  // extra parameters
  Double_t                               fMedUsedForBgSub;    
  Bool_t                                 fUseArea4Vector;    
  Bool_t				                         fBackSub;
  Bool_t				                         fSubtractMc;
  Bool_t 				                         fPhiModulation;
  Int_t  		    		                     fJetNHardestSkipped_010;
  Int_t    				                       fJetNHardestSkipped_1080;
  Bool_t                  		           fSetJetFixedSeed;
  Int_t                                  fJetFJSeed;
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
  std::vector<Double_t>                      fGRNumerator;      //!
  std::vector<Double_t>                      fGRDenominator;    //!
  std::vector<Double_t>                      fGRNumeratorSub;   //!
  std::vector<Double_t>                      fGRDenominatorSub; //!
  std::vector<Double_t>                     fICSMaxDistances;  //!
  std::vector<Double_t>                     fICSAlphas; 
 // virtual void   SubtractBackground(const Double_t median_pt = -1);

 private:
  StFJWrapper();
  StFJWrapper(const StFJWrapper& wrapper);
  StFJWrapper& operator = (const StFJWrapper& wrapper);
  fastjet::JetMedianBackgroundEstimator   CreateBackgroundEstimator();
  fastjet::contrib::BackgroundRescalingYPhi CreateBackgroundRescaling();
  std::vector<fastjet::PseudoJet>         CreateInputRealOnly() const;
  fastjet::AreaDefinition                      CreateAreaDefinition(bool explicitGhosts) const;
  void                                      SplitOutD0(const std::vector<fastjet::PseudoJet>& input,
                                                     std::vector<fastjet::PseudoJet>& out_without_D0,
                                                     std::vector<fastjet::PseudoJet>& out_D0) const;
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
fj::contrib::BackgroundRescalingYPhi StFJWrapper::CreateBackgroundRescaling()
{
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

  fj::contrib::BackgroundRescalingYPhi rescaling(v2, v3, v4, psi, 0., 1., 0., 1.);
  //y - not used, for RHIC it is particaly flat
  rescaling.use_rap_term(false);
  //phi used - but no visible impact
  rescaling.use_phi_term(true);

  return rescaling;
}

//_________________________________________________________________________________________________
fastjet::JetMedianBackgroundEstimator StFJWrapper::CreateBackgroundEstimator() {

    // Define jet algorithm and area for background estimation
    fj::JetDefinition jet_def_bkgd(fj::kt_algorithm, fR, fj::E_scheme, fj::Best);
    fj::AreaDefinition area_def_bkgd(fj::active_area_explicit_ghosts, fj::GhostedAreaSpec(1.2, 1, 0.005));

    // For tuning, fix the seed for fastjet
    if (fSetJetFixedSeed) {
        Int_t seed1 = fJetFJSeed;
        Int_t seed2 = fJetFJSeed;
        std::vector<Int_t> seeds = { static_cast<Int_t>(seed1), static_cast<Int_t>(seed2) };
        area_def_bkgd = area_def_bkgd.with_fixed_seed(seeds);
    }

    // Decide how many jets to remove based on centrality
    Int_t nJetsRemove = fJetNHardestSkipped_1080;
    if (fCentrality == 7 || fCentrality == 8) nJetsRemove = fJetNHardestSkipped_010;

    // Selector for background estimation (eta cut and remove hardest jets)
    fj::Selector selector = (!fj::SelectorNHardest(nJetsRemove)) * fj::SelectorAbsEtaMax(0.6);

    // Create and return the background estimator
    fj::JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
    return bkgd_estimator;
}
//_________________________________________________________________________________________________
Double_t StFJWrapper::ComputeOccupancyC_fromJetsUsed(const fastjet::JetMedianBackgroundEstimator& bge) {
  const Double_t etaMax  = 0.6;
  const Double_t A_range = (2.0*etaMax) * (2.0*TMath::Pi());
  if (A_range <= 0) return 1.0;

  Double_t A_occ = 0.0;
  const std::vector<fastjet::PseudoJet> ju = bge.jets_used();
  for (const auto& j : ju) {
    if (j.perp() > 0) A_occ += j.area();
  }

  Double_t C = A_occ / A_range;
  if (C < 0) C = 0;
  if (C > 1) C = 1;
  return C;
}

//_________________________________________________________________________________________________
std::vector<fastjet::PseudoJet> StFJWrapper::CreateInputRealOnly() const {

    std::vector<fastjet::PseudoJet> inputRealOnly;
    inputRealOnly.reserve(fInputVectors.size());

    for (const auto &pj : fInputVectors) {
      const int uid = pj.user_index();
      if (fSubtractMc && uid >= 10000)  continue;
      if (fSubtractMc && uid <= -10000) continue;
      inputRealOnly.push_back(pj);
    }

    return inputRealOnly;
}

//_________________________________________________________________________________________________
void StFJWrapper::SplitOutD0(const std::vector<fastjet::PseudoJet>& input,
                             std::vector<fastjet::PseudoJet>& out_without_D0,
                             std::vector<fastjet::PseudoJet>& out_D0) const
{
  out_without_D0.clear();
  out_D0.clear();

  fj::Selector notD0sel = fastjet::SelectorMassMax(1); // 1 GeV max mass

  for (size_t i = 0; i < input.size(); ++i) {
    const fastjet::PseudoJet &p = input[i];

    if (notD0sel(p) && p.user_index() < 30000) { // exclude D0 candidates (D0 mass > 1 GeV/c^2 || user_index >= 30000)
      fastjet::PseudoJet temp_jet(p.px(), p.py(), p.pz(), p.E());
      temp_jet.set_user_index(p.user_index());
      out_without_D0.push_back(temp_jet);
    } else {
      fastjet::PseudoJet temp_jet(p.px(), p.py(), p.pz(), p.E());
      temp_jet.set_user_index(p.user_index());
      out_D0.push_back(temp_jet);
    }
  }
}

//_________________________________________________________________________________________________
fj::AreaDefinition StFJWrapper::CreateAreaDefinition(bool explicitGhosts) const {

  fj::GhostedAreaSpec gas(1.2, 1, 0.005);
  fj::AreaDefinition area_def = explicitGhosts
    ? fj::AreaDefinition(fj::active_area_explicit_ghosts, gas)
    : fj::AreaDefinition(fj::active_area, gas);

  if (fSetJetFixedSeed) {
    Int_t seed1 = fJetFJSeed;
    Int_t seed2 = fJetFJSeed;
    std::vector<int> seeds = { static_cast<int>(seed1), static_cast<int>(seed2) };
    area_def = area_def.with_fixed_seed(seeds);
  }

  return area_def;
}
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
  , fBackSub	       (kTRUE)
  , fSubtractMc        (kTRUE)
  , fPhiModulation     (kFALSE)
  , fJetNHardestSkipped_010 (2)
  , fJetNHardestSkipped_1080 (1)
  , fSetJetFixedSeed   (kFALSE)
  , fJetFJSeed         (12345)
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
void StFJWrapper::Clear(const Option_t */*opt*/) {
  // Simply clear the input vectors.
  // Make sure done on every event if the instance is reused
  // Reset the median to zero.

  fInputVectors.clear();
  fInputGhosts.clear();
  fMedUsedForBgSub = 0;
  fAngul10half = -999;
  fAngul11 = -999;
  fAngul11half = -999;
  fAngul12 = -999;
  fAngul13 = -999;
  fAngulDisp = -999;

  // for the moment brute force delete everything
  ClearMemory();
}

//_________________________________________________________________________________________________
void StFJWrapper::RemoveLastInputVector() {
  // Remove last input vector
  fInputVectors.pop_back();
}
//_________________________________________________________________________________________________
void StFJWrapper::PrintInput() {
  cout << "        pt eta phi" << endl;
  for (int i = 0; i < fInputVectors.size(); i++){
    cout << Form("Input # %i \t %.2f \t %.2f \t %.2f \t %.2f ", i, fInputVectors[i].pt(), fInputVectors[i].eta(), fInputVectors[i].phi(), fInputVectors[i].e()) << endl;
  }
}

//_________________________________________________________________________________________________
void StFJWrapper::AddInputVector(Double_t px, Double_t py, Double_t pz, Double_t E, Int_t index) {
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
void StFJWrapper::AddInputVector(const fj::PseudoJet& vec, Int_t index) {
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
void StFJWrapper::AddInputVectors(const std::vector<fj::PseudoJet>& vecs) {
  // Add the input from vector of pseudojets.
  for (UInt_t i = 0; i < vecs.size(); ++i) {
    fj::PseudoJet inVec = vecs[i];
    // add to the fj container of input vectors
    fInputVectors.push_back(inVec);
  }
}

//_________________________________________________________________________________________________
Double_t StFJWrapper::GetJetArea(UInt_t idx) const {
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
fastjet::PseudoJet StFJWrapper::GetJetAreaVector(UInt_t idx) const {
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
std::vector<fastjet::PseudoJet> StFJWrapper::GetJetConstituents(UInt_t idx) const {
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

Int_t StFJWrapper::runICSMethod() {

  //--- Set up Background Estimation and Rescaling ---

  //Background scaling
  fj::contrib::BackgroundRescalingYPhi rescaling = CreateBackgroundRescaling();

  //Background estimation
  fj::JetMedianBackgroundEstimator bkgd_estimator = CreateBackgroundEstimator();

  //Exclude MC for background estimation
  std::vector<fastjet::PseudoJet> inputRealOnly = CreateInputRealOnly();
        
  //Rescaling
	if (fPhiModulation) bkgd_estimator.set_rescaling_class(&rescaling);

	// Set particles for background estimator (background estimation for the input particles)
	bkgd_estimator.set_particles(inputRealOnly); //MC cant be included

  double C = ComputeOccupancyC_fromJetsUsed(bkgd_estimator);
  C = 1.0; // Temporarily disable occupancy correction
  ScaledBGE scaled(&bkgd_estimator, C);

  //--- Event wide ICS background subtraction ---

  // Split out D0 candidates from the input vectors
  // D0 can't be affected by the background subtraction
  std::vector<fastjet::PseudoJet> inputWithoutD0;
	std::vector<fastjet::PseudoJet> onlyD0Meson;
  SplitOutD0(fInputVectors, inputWithoutD0, onlyD0Meson);

  // Set up the ICS background subtraction
	fj::contrib::IterativeConstituentSubtractor subtractor; // Set up the background subtraction algorithm

  //Check parameters
	if (fICSMaxDistances.empty() || fICSAlphas.empty() || fICSMaxDistances.size() != fICSAlphas.size()){
	    std::cerr << "ERROR: ICS subtraction parameters not set correctly!" << std::endl;
	    exit(1);
	}

 	// Apply the subtraction parameters
	subtractor.set_distance_type(fj::contrib::ConstituentSubtractor::deltaR); // Set distance type for subtraction
	subtractor.set_parameters(fICSMaxDistances, fICSAlphas); // Set ICS subtraction parameters
	subtractor.set_ghost_removal(true); // Enable ghost removal
	subtractor.set_ghost_area(0.005); // Set ghost area value
	subtractor.set_max_eta(1); // Set maximum eta for particles (all particles within |eta| < 1 are considered anyway)
	//subtractor.set_background_estimator(&bkgd_estimator); // Link the background estimator
  subtractor.set_background_estimator(&scaled);
	subtractor.set_common_bge_for_rho_and_rhom(true); // Set common background estimation for rho and rhom
	subtractor.set_keep_original_masses(); // Keep the original masses of particles
	subtractor.set_scale_fourmomentum(); // Scale four-momentum of particles after subtraction
	
  // Initialize the background subtraction algorithm
	subtractor.initialize(); 

	// Minimum pt for jets (below this, jets are excluded)
	vector<fj::PseudoJet> corrected_event;

	// Apply background subtraction if enabled
	if (fBackSub) corrected_event = subtractor.subtract_event(inputWithoutD0); 
	else corrected_event = inputWithoutD0; // Otherwise, use original input vectors

	// Fill histograms for ICS jet constituents
	if (fBackSub){

    fJetRho = scaled.rho();
    fJetRhoM = scaled.rho_m();
	  for (vector<fastjet::PseudoJet>::const_iterator particle = corrected_event.begin(); particle != corrected_event.end(); ++particle) {
      hJetConstRapPhiICS->Fill(particle->phi_std(), particle->rap(),fCentralityWeight);
      hJetConstEtaPhiICS->Fill(particle->phi_std(), particle->eta(),fCentralityWeight);
      hJetConstPtICS->Fill(particle->perp(),fCentralityWeight);		
	  }
	
	}

  // Return back D0
	if (!onlyD0Meson.empty() && onlyD0Meson.size() == 1) corrected_event.push_back(onlyD0Meson.back());
	else std::cerr << "WARNING: "<< onlyD0Meson.size() <<" D0 meson(s) found in the event!" << std::endl;

  //--- Jet reconstruction after event-wide ICS background subtraction ---

	// Define jet algorithm for actual jet reconstruction (Anti-kt with radius fR)
	fj::JetDefinition jet_def(fj::antikt_algorithm, fR, fj::E_scheme, fj::Best);

  // Define area for jet reconstruction
  fj::AreaDefinition area_def_jet = CreateAreaDefinition(false); //false - no explicit ghosts for jet clustering (due to ICS)

	// Perform jet clustering with the background-subtracted (or original) particles
  if (fClustSeq){delete fClustSeq; fClustSeq = nullptr;}
	fClustSeq = new fastjet::ClusterSequenceArea(corrected_event, jet_def, area_def_jet);
	
	// Get inclusive jets
	fInclusiveJets.clear();
	fInclusiveJets = fClustSeq->inclusive_jets();

  return 0;
}

Int_t StFJWrapper::runAreaBgAndAreaShapeMethod(){

	//--- Jet Reconstruction ---

	// Define jet algorithm for actual jet reconstruction (Anti-kt with radius fR)
	fastjet::JetDefinition jet_def(fj::antikt_algorithm, fR, fj::E_scheme, fj::Best);

  // Define area for jet reconstruction
  fastjet::AreaDefinition area_def_jet = CreateAreaDefinition(true); //true - explicit ghosts for jet clustering
	
	// Perform jet clustering with the background-subtracted (or original) particles
  if (fClustSeq){delete fClustSeq; fClustSeq = nullptr;}
	fClustSeq = new fj::ClusterSequenceArea(fInputVectors, jet_def, area_def_jet);
        
	// Get inclusive jets
	fInclusiveJets.clear();
	fInclusiveJets = fClustSeq->inclusive_jets();
	
  //--- Background Estimation ---
  
  //Background estimation
  fastjet::JetMedianBackgroundEstimator bkgd_estimator = CreateBackgroundEstimator();

  //Background scaling
  fastjet::contrib::BackgroundRescalingYPhi rescaling = CreateBackgroundRescaling();

  // Exclude MC for Background estimation
  std::vector<fastjet::PseudoJet> inputRealOnly = CreateInputRealOnly();

	//Estimation of the background using only charged tracks
	if (fPhiModulation) bkgd_estimator.set_rescaling_class(&rescaling);
  
	bkgd_estimator.set_particles(inputRealOnly);
	
  Double_t C = 0;

  fJetRho = 0;
  fJetRhoM = 0;

	if (fBackSub){
    C = ComputeOccupancyC_fromJetsUsed(bkgd_estimator);//
    C = 1.0; // Temporarily disable occupancy correction
	  fJetRho = C*bkgd_estimator.rho();
	  fJetRhoM = C*bkgd_estimator.rho_m(); 
    //cout << Form("StFJWrapper::runAreaBgAndAreaShapeMethod: Rho = %.2f, RhoM = %.2f, C = %.2f", fJetRho, fJetRhoM, C) << endl;
	}

  if (fJetRho < 0) fJetRho = 0;
  if (fJetRhoM < 0) fJetRhoM = 0;
   //--- Generic Subtractor Setup ---
	//fastjet::contrib::GenericSubtractor gensub(&bkgd_estimator);
  //gensub.set_background_density(fJetRho, fJetRhoM);//
	fastjet::contrib::GenericSubtractor gensub(fJetRho, fJetRhoM);

	/////gensub.set_common_bge_for_rho_and_rhom(true);
	fastjet::contrib::GenericSubtractorInfo info;

  //--- Jet Shapes Calculation ---

	//kappa,alfa,jet_R
	Angularity my_angularity_10half(  1,  0.5,  fR);
	Angularity my_angularity_11(      1,  1,    fR);
	Angularity my_angularity_11half(  1,  1.5,  fR);
	Angularity my_angularity_12(      1,  2,    fR);
	Angularity my_angularity_13(      1,  3,    fR);
	Angularity my_angularity_Disp(    2,  0,    fR);
	
	fastjet::FunctionOfPseudoJet<Double_t>* shape10half = &my_angularity_10half;
	fastjet::FunctionOfPseudoJet<Double_t>* shape11 = &my_angularity_11;
	fastjet::FunctionOfPseudoJet<Double_t>* shape11half = &my_angularity_11half;
	fastjet::FunctionOfPseudoJet<Double_t>* shape12 = &my_angularity_12;
	fastjet::FunctionOfPseudoJet<Double_t>* shape13 = &my_angularity_13;
	fastjet::FunctionOfPseudoJet<Double_t>* shapeDisp = &my_angularity_Disp;
				
  // Loop over jets and calculate shapes
  for (unsigned int i = 0; i < fInclusiveJets.size(); i++) {

    // Get i-th jet
    const fastjet::PseudoJet &jet = fInclusiveJets[i];

    // Get i-th jet constituents
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    
    // Count D0 mesons in the jet
    Int_t d0_count = 0;

    // Loop over constituents to identify D0 mesons
    for (size_t j = 0; j < constituents.size(); j++) {
      
      // Looking for D0 (mass should be > 1.0 GeV/c^2 (but it doesn't have to) or user_index >= 30000)
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
/*
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
}*/

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
