#ifndef CUTS_H
#define CUTS_H

/* **************************************************
 *  Cuts namespace.
 *
 *  Authors:  Xin Dong        (xdong@lbl.gov)
 *            Michael Lomnitz (mrlomnitz@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)   
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>
#include <set>
#include <vector>
#include <../../fastjet/config.h>
#include <../../fastjet/PseudoJet.hh>
#include <../../fastjet/JetDefinition.hh>
#include <../../fastjet/ClusterSequence.hh>
#include <../../fastjet/ClusterSequenceArea.hh>

#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoEvent/StPicoBTowHit.h"

#include "StEmcUtil/geometry/StEmcGeom.h"


//#include "../StPicoJetMaker/StPicoJetMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TPythia6.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TDatime.h"
#include <array>

using namespace std;
namespace mycuts
{
   extern std::string const prescalesFilesDirectoryName;
   //event
   extern UShort_t const triggerWord;
   extern float const vz;
   extern float const vzVpdVz;
    extern float const vr;

   //tracking
   extern int const nHitsFit;
   extern bool const requireHFT;
   extern float const minPt;

   //pions
   extern float const nSigmaPion;
   extern float const pTofBetaDiff;
   
   //kaons
   extern float const nSigmaKaon;
   extern float const kTofBetaDiff;

   // tree kaonPion pair cuts
   extern float const cosTheta;
   extern float const cosTheta_2014;
   extern float const cosTheta_2016;
   extern float const dcaDaughters;
   extern float const decayLength;
   extern float const minMass;
   extern float const maxMass;
   extern float const kDca;
   extern float const pDca;

   //Jet track cuts
    extern float const jetTrackPtMin;
    extern float const jetTrackPtMax;
    extern float const jetTrackEta;
    extern float const jetTracknHitsFit;
    extern float const jetTracknHitsRatio; // nHitsFit/nHitsMax >= 0.52
    extern float const jetTrackDCA;

   // histograms kaonPion pair cuts
   extern float const qaNHitsFit;
   extern float const qaNSigmaKaon;
   extern float const qaCosTheta;
   extern float const qaDcaDaughters;
   extern float const qaKDca;
   extern float const qaPDca;
   //Hadron cuts
   extern float const hadronPtMin;
   extern float const hadronPtMax;
   extern float const corDetaMin;
   extern float const corDetaMax;

   extern float const pionDCA_cut_2014[6][5];
   extern float const kaonDCA_cut_2014[6][5];
   extern float const DCA_D0_cut_2014[6][5];
   extern float const D0_decayLength_cut_2014[6][5];
   extern float const pionkaonDCA_cut_2014[6][5];

   extern float const pionDCA_cut_2016[6][5];
   extern float const kaonDCA_cut_2016[6][5];
   extern float const DCA_D0_cut_2016[6][5];
   extern float const D0_decayLength_cut_2016[6][5];
   extern float const pionkaonDCA_cut_2016[6][5];

   extern const std::set<int> mbTriggers2014;
   extern const std::set<int> mbTriggers2016;
   
   extern const std::set<int> goodRun2014;
   extern const std::set<int> AllBadRunList2014;
   
   extern const int BadTowerArrT[433];
   extern const int BadTowerArr[822];
   extern const int BadTowerMap[4800];
   extern const std::set<int> NeilBadTowers2014;
   extern const Double_t CLowMidHigh[4800];
   
   extern const Double_t CPre[4800];
   extern array<float, 4800> SumE;
}
#endif
