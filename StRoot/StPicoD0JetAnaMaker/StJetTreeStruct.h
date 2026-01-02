#ifdef __CINT__
#pragma link C++ class vector<float>+;
#pragma link C++ class vector<vector<float>>+;
#pragma link C++ class vector<int>+;
#pragma link C++ class vector<vector<int>>+;
#endif

#ifndef StJetTreeStruct_h
#define StJetTreeStruct_h

const int ConstMax = 5000;
struct StJetTreeStruct
{
    /* There's only one jet per event, due to the structure of our analysis. If an event has more than one jet of interest, that eventid will have a frequency > 1.
    Nothing else should change.
    The event level info we store are            : run id, event id, event refmult, event centrality, event triggers.
    The jet level info we store are              : jetpt, jetcorrectedpt, jeteta, jetphi, jetarea, jetradius, jetenergy, jetnef (neutral energy fraction), fRhoVal used,
                                                   highest energy track pt, number of constituents
    The constituent level info we store are      : trackid, trackpt, tracketa, trackphi, trackpx, trackpy, trackpz
    */

    // Jet Data Info
    int runid;
    int eventid;
    float refmult;
    float grefmult;
    float centrality;
    float centralityAlt;
    float refcorr2;
    float mcrefmult;
    float weight;
    vector<unsigned int> triggers;
    vector<double> primaryvertex;
    vector<double> primaryvertexerror;


    float corrPsi2;
    float Psi2;
    float Q1_vec;
    float Q2_vec;
    float Q1_vec_rec;
    float Q2_vec_rec;
    float maxtrackpt;
    float maxtowerEt;

    float jetpt;
    float jetptcorr;
    float jetrho;
    float jetrhom;
    float jetcorrectedpt;
    float jeteta;
    float jetphi;
    float jetarea;
    float jetradius;
    float jetenergy;
    float jetnef;
    float fRhoValforjet;
    float fSigValforjet;
    float jethighesttrackpt;
    int numberofconstituents;

//Parama
 float jetpt_t[4][1];
 float alpha10half_t[4][1];
 float alpha11_t[4][1];
 float alpha11half_t[4][1];
 float alpha12_t[4][1];
 float alpha13_t[4][1];
 float nOfAllConst[4][1];
 float nOfAllMcConst[4][1];
 float MomDisp[4][1];
 float D0z_t[4][1];
 

 float shapeAlpha11;
 // 0 

 
    float d0z;
    float d0mass;
    float d0pt;
    float d0phi;
    float d0eta;
    float d0DeltaR;

    float pionpt;
    float pioneta;
    float pionphi;
    float pioncharge;

    float kaonpt;
    float kaoneta;
    float kaonphi;
    float kaoncharge;

    float lambda_1_half;
    float lambda_1_1;
    float lambda_1_1half;
    float lambda_1_2;
    float lambda_1_3;
    float dispersion;
    
    float shapeLambda_1_half;
    float shapeLambda_1_1;
    float shapeLambda_1_1half;
    float shapeLambda_1_2;
    float shapeLambda_1_3;
    float shapeDispersion;

    float mTrackID[ConstMax];
    float mTrackPt[ConstMax];
    float mTrackEta[ConstMax];
    float mTrackPhi[ConstMax];
    float mTrackPx[ConstMax];
    float mTrackPy[ConstMax];
    float mTrackPz[ConstMax];
    float mTrackCharge[ConstMax];

    void Clear()
    {
        runid = -999;
        eventid = -999;
        refmult = -999;
        grefmult = -999;
        centrality = -999;
        centralityAlt = -999;
        refcorr2 = -999;
        mcrefmult = -999;
        weight = -999;
        triggers.clear();
        primaryvertex.clear();
        primaryvertexerror.clear();

        maxtrackpt = -999;
        maxtowerEt = -999;

        jetpt = -999;
        jetptcorr = -999;
        jetrho = -999;
        jetrhom = -999;
        jetcorrectedpt = -999;
        jeteta = -999;
        jetphi = -999;

        jetarea = -999;
        jetradius = -999;
        jetenergy = -999;
        jetnef = -999;
        fRhoValforjet = -999;
        fSigValforjet = -999;
        jethighesttrackpt = -999;
        shapeLambda_1_half = -999;
        shapeLambda_1_1 = -999;
        shapeLambda_1_1half = -999;
        shapeLambda_1_2 = -999;
        shapeLambda_1_3 = -999;
        shapeDispersion = -999;
        numberofconstituents = -999;

        d0z = -999;
        d0mass = -999;
        d0pt = -999;
        d0phi = -999;
        d0eta = -999;
        d0DeltaR = -999;

        pionpt = -999;
        pioneta = -999;
        pionphi = -999;
        pioncharge = -999;

        kaonpt = -999;
        kaoneta = -999;
        kaonphi = -999;
        kaoncharge = -999;

	lambda_1_half = -999;
        lambda_1_1 = -999;
        lambda_1_1half = -999;
        lambda_1_2 = -999;
        lambda_1_3 = -999;
        dispersion = -999;
        corrPsi2 = -999;
        Q1_vec = -999;
        Q2_vec = -999;
        
        shapeAlpha11 = -999;
        
        //Parama
        for (int ind1 = 0; ind1 < 4; ind1++){
                for (int ind2 = 0; ind2 < 1; ind2++){
		 	jetpt_t[ind1][ind2] = -999;
		 	D0z_t[ind1][ind2] = -999;
	 		alpha10half_t[ind1][ind2] = -999;	
	 		alpha11_t[ind1][ind2] = -999;	 		
	 		alpha11half_t[ind1][ind2] = -999;
	 		alpha12_t[ind1][ind2] = -999;	
	 		alpha13_t[ind1][ind2] = -999;		
	 		nOfAllConst[ind1][ind2] = -999;
			nOfAllMcConst[ind1][ind2] = -999;
			MomDisp[ind1][ind2] = -999;
		}
        }
        
    }
};

#endif
