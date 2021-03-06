// std includes
#include <iostream>
#include <string>
#include <vector>
#include <math.h>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "AnalysisSelector.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "RazorMultiB.hh"

RazorMultiB::RazorMultiB(TTree *tree) : Vecbos(tree) {
    _goodRunLS = false;
    _isData = false;
    _weight=1.0;
    _isSMS = false;
}

RazorMultiB::RazorMultiB(TTree *tree, string jsonFile, bool goodRunLS, bool isData) : Vecbos(tree) {
    _goodRunLS = goodRunLS;
    _isData = isData;
    _weight=1.0;
    _isSMS = false;
    
    //To read good run list!
    if (goodRunLS && isData) {
        setJsonGoodRunList(jsonFile);
        fillRunLSMap();
        InitEventFlag("/afs/cern.ch/user/w/woodson/public/WEIGHT/AllBadABCDNEWTAUID.txt");
    }
    
}


RazorMultiB::RazorMultiB(TTree *tree, string jsonFile, bool goodRunLS, bool isData, string smsName) : Vecbos(tree) {
    _goodRunLS = goodRunLS;
    _isData = isData;
    _weight=1.0;
    _isSMS = false;
    SMS_temp = smsName;
    if (smsName!="none") {
      _isSMS = true;
    }
    //To read good run list!
    if (goodRunLS && isData) {
        setJsonGoodRunList(jsonFile);
        fillRunLSMap();
        InitEventFlag("/afs/cern.ch/user/w/woodson/public/WEIGHT/AllBadABCDNEWTAUID.txt");
    }
}

RazorMultiB::~RazorMultiB() {}

void RazorMultiB::SetConditions(TTree* treeCond) {
    _treeCond = treeCond;
}

void RazorMultiB::SetWeight(double weight){
    _weight=weight;
}

void RazorMultiB::Loop(string outFileName, int start, int stop) {
    if(fChain == 0) return;
    
    // nominal triggers
    int HLT_DoubleMu;
    int HLT_DoubleEle;
    int HLT_MuEle;
    // prescaled triggers
    int HLT_Mu17;
    int HLT_Mu12;
    int HLT_Ele17;
    int HLT_Ele8;
        
    // PF block k
    double run;
    double evNum;
    double bx;
    double ls;
    double orbit;
    double PFR;
    double RSQ;
    double MR;
    double MRT;
    double pTPFHem1;
    double etaPFHem1;
    double phiPFHem1;
    double pTPFHem2;
    double etaPFHem2;
    double phiPFHem2;
    int    nBtag;
    int    nBtag_medium;
    int    nBtag_tight;
    double W=_weight;
    float mg, mchi;
    int    BOX_NUM;
    int    ss;
    int    nPV;
    
    // fast-hemispheres
    double FJR;
    double FJRSQ;
    double FJMR;
    double FJMRT;
    double pTFJHem1;
    double etaFJHem1;
    double phiFJHem1;
    double pTFJHem2;
    double etaFJHem2;
    double phiFJHem2;
    
    // gen level info
    double pT1, pT2, eta1, eta2, phi1, phi2;
    int idMcL1, idMothMcL1, idGrandMothMcL1;
    int idMcL2, idMothMcL2, idGrandMothMcL2;
    //int idMcB1, idMothMcB1, idGrandMothMcB1;
    //int idMcB2, idMothMcB2, idGrandMothMcB2;
    double genH1Mass, genH2Mass;
    double genGammaCM, genshat, genMDR, genMDR2;

    // gen-level lepton and b-quark 4-vectors
    double gpXL1, gpYL1, gpZL1, genergyL1;
    double gpXL2, gpYL2, gpZL2, genergyL2;
    double gpXB_t, gpYB_t, gpZB_t, genergyB_t;
    double gpXB_tb, gpYB_tb, gpZB_tb, genergyB_tb;

    // selected MC lepton and b-quark 4-vectors
    double pXL1, pYL1, pZL1, energyL1;
    double pXL2, pYL2, pZL2, energyL2;
    double pXB1, pYB1, pZB1, energyB1;
    double pXB2, pYB2, pZB2, energyB2;
    
    // lepton angles
    double dPhi_ll;
    // pT corrected lepton angle
    double dPhi_ll_R;
    double dPhi_ll_R_2;
    // helicity angle for bl
    double CosThetaBL1;
    double CosThetaBL2;
    // angle between bl in lab
    double dPhi_bl1;
    double dPhi_bl2;
    double dPhi_bl12;
    double dPhi_bl21;
    // pT corrected b and l angle
    double dPhi_bl1_R;
    double dPhi_bl2_R;
    double dPhi_bl12_R;
    double dPhi_bl21_R;
    //double gdPhi_ll;
    // angle between bl hemispheres in lab   
    double dPhi_bl;
    double dPhi_bl_R;
    // ttbar decay: 0 = nolep, 1 = semilep; 2 = fully lep
    int nLepTopDecay;
    int nNeutrino;
    TLorentzVector pTNeutrino;
    double pXMc;
    double pYMc;
    double pTNeutrinoMag;

    // PFElectron Block
    double pfElectron_pt;
    double pfElectron_eta;
    double pfElectron_phi;
    double pfElectron_energy;
    double pfElectron_mass = 0.511/1000;
    
    // PFMuon Block
    double pfMuon_pt;
    double pfMuon_eta;
    double pfMuon_phi;
    double pfMuon_energy;
    double pfMuon_mass = 105.7/1000;
    
    // Counters
    int nIsolatedPFJets;
    int nPFJets;
    int nBtag_lead4jets;
    int nBtag_medium_lead4jets;
    int nBtag_tight_lead4jets;
    int nBtag_TCHPT;
    int nBtag_TCHPT_lead4jets;
    double Mll;
    
    // New Razor Variables RPV 
    double MR_pTcorr;
    double gammaR;
    double shatR_bl;
    double dPhiCM;
    double EB1;
    double EB2;
    double CosThetaB1;
    double CosThetaB2;
    double TopMass1;
    double TopMass2;
    double EL1;
    double EL2;
    double CosThetaL1;
    double CosThetaL2;
    double GluinoMass1;
    double GluinoMass2;
    double VisHemMass1;
    double VisHemMass2;
    
    //Variable Arrays
    //Index 0 is Hybrid, Index 1 is pT-corrected (Rogan)
    //Index 2 is Canonical
    double TotalHemMass1List[3];
    double TotalHemMass2List[3];
    double TopHemMass1List[3];
    double TopHemMass2List[3];
    double MR1List[3];
    double MR2List[3];
    double pTNMagList[3];
    double CosThetaB1List[3];
    double CosThetaB2List[3];
    double CosThetaL1List[3];
    double CosThetaL2List[3];
    double gammaRList[3];
	
    //check btag efficiencies
    float BDiscList[2];

    //original MR using only b and l
    double MR_og;
    double pT_CM;

    //pT-corrected MTR, RSQ
    double MRT_bl;
    double RSQ_bl;

    //corrected shatR
    double rat;
    double shat_corr;
    double shat_m;
 
    //Hybrid Razor Approach
    double TotalHemMass1;
    double TotalHemMass2;
    double TopHemMass1;
    double TopHemMass2;
    double MR1;
    double MR2;
    double TotalNMag;

    //pT-corrected (Rogan) Razor Approach
    double TotalHemMass1Raz;
    double TotalHemMass2Raz;
    double TopHemMass1Raz;
    double TopHemMass2Raz;
    double MR1Raz;
    double MR2Raz;
    double TotalNMagRaz;

    //Hybrid Razor Approach + Transverse Boost
    double TotalHemMass1Trans;
    double TotalHemMass2Trans;
    double TopHemMass1Trans;
    double TopHemMass2Trans;
    double MR1Trans;
    double MR2Trans;
    double TotalNMagTrans;
    //Canonical Razor Approach
    double TotalHemMass1Star;
    double TotalHemMass2Star;
    double TopHemMass1Star;
    double TopHemMass2Star;
    double MR1Star;
    double MR2Star;
    double TotalNMagStar;
    
    
    //
    double MetMag;
    double HT;
    double MHT_x;
    double MHT_y;
    double MHT;
    
    // prepare the output tree
    TTree* outTree = new TTree("outTree", "outTree");
    outTree->Branch("run", &run, "run/D");
    outTree->Branch("evNum", &evNum, "evNum/D");
    outTree->Branch("bx", &bx, "bx/D");
    outTree->Branch("ls", &ls, "ls/D");
    outTree->Branch("orbit", &orbit, "orbit/D");
    outTree->Branch("BOX_NUM", &BOX_NUM, "BOX_NUM/I");
    outTree->Branch("ss", &ss, "ss/I");
    
    // HLT bits
    outTree->Branch("HLT_DoubleMu", &HLT_DoubleMu, "HLT_DoubleMu/I");
    outTree->Branch("HLT_DoubleEle", &HLT_DoubleEle, "HLT_DoubleEle/I");
    outTree->Branch("HLT_MuEle", &HLT_MuEle, "HLT_MuEle/I");
    outTree->Branch("HLT_Mu17", &HLT_Mu17, "HLT_Mu17/I");
    outTree->Branch("HLT_Mu12", &HLT_Mu12, "HLT_Mu12/I");
    outTree->Branch("HLT_Ele17", &HLT_Ele17, "HLT_Ele17/I");
    outTree->Branch("HLT_Ele8", &HLT_Ele8, "HLT_Ele8/I");
    
    // PF block
    outTree->Branch("PFR", &PFR, "PFR/D");
    outTree->Branch("RSQ", &RSQ, "RSQ/D");
    outTree->Branch("MR", &MR, "MR/D");
    outTree->Branch("MRT", &MRT, "MRT/D");
    outTree->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
    outTree->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
    outTree->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
    outTree->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
    outTree->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
    outTree->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
    outTree->Branch("nBtag", &nBtag, "nBtag/I");
    outTree->Branch("nBtag_medium", &nBtag_medium, "nBtag_medium/I");
    outTree->Branch("nBtag_tight", &nBtag_tight, "nBtag_tight/I");
    outTree->Branch("nBtag_TCHPT", &nBtag_TCHPT, "nBtag_TCHPT/I");
    outTree->Branch("W", &W, "W/D");
    outTree->Branch("mg", &mg, "mg/F");
    outTree->Branch("mchi", &mchi, "mchi/F");
    outTree->Branch("nPV", &nPV, "nPV/I");
    
    // original MR using only b and l
    outTree->Branch("MR_og", &MR_og, "MR_og/D");
    
    outTree->Branch("pT_CM", &pT_CM, "pT_CM/D");

    // pT-corr RSQ, MRT
    outTree->Branch("RSQ_bl", &RSQ_bl, "RSQ_bl/D");
    outTree->Branch("MRT_bl", &MRT_bl, "MRT_bl/D");

    // corrected shatR
    outTree->Branch("rat", &rat, "rat/D");
    outTree->Branch("shat_corr", &shat_corr, "shat_corr/D");
    outTree->Branch("shat_m", &shat_m, "shat_m/D");

    // fast-hemispheres
    outTree->Branch("FJR", &FJR, "FJR/D");
    outTree->Branch("FJRSQ", &FJRSQ, "FJRSQ/D");
    outTree->Branch("FJMR", &FJMR, "FJMR/D");
    outTree->Branch("FJMRT", &FJMRT, "FJMRT/D");
    outTree->Branch("pTFJHem1", &pTFJHem1, "pTFJHem1/D");
    outTree->Branch("etaFJHem1", &etaFJHem1, "etaFJHem1/D");
    outTree->Branch("phiFJHem1", &phiFJHem1, "phiPFJHem1/D");
    outTree->Branch("pTFJHem2", &pTFJHem2, "pTFJHem2/D");
    outTree->Branch("etaFJHem2", &etaFJHem2, "etaFJHem2/D");
    outTree->Branch("phiFJHem2", &phiFJHem2, "phiFJHem2/D");
    
    // Arrays of Variables
    outTree->Branch("TotalHemMass1List", TotalHemMass1List, "TotalHemMass1List[3]/D");
    outTree->Branch("TotalHemMass2List", TotalHemMass2List, "TotalHemMass2List[3]/D");
    outTree->Branch("TopHemMass1List", TopHemMass1List, "TopHemMass1List[3]/D");
    outTree->Branch("TopHemMass2List", TopHemMass2List, "TopHemMass2List[3]/D");
    outTree->Branch("MR1List", MR1List, "MR1List[3]/D");
    outTree->Branch("MR2List", MR2List, "MR2List[3]/D");
    outTree->Branch("pTNMagList", pTNMagList, "pTNMagList[3]/D");
    outTree->Branch("CosThetaB1List", CosThetaB1List, "CosThetaB1List[3]/D");
    outTree->Branch("CosThetaB2List", CosThetaB2List, "CosThetaB2List[3]/D");
    outTree->Branch("CosThetaL1List", CosThetaL1List, "CosThetaL1List[3]/D");
    outTree->Branch("CosThetaL2List", CosThetaL2List, "CosThetaL2List[3]/D");
    outTree->Branch("gammaRList", gammaRList, "gammaRList[3]/D");
    
    //b tag discriminators
    outTree->Branch("BDiscList", BDiscList, "BDiscList[2]/F");
		
    // New Razor Variables
    outTree->Branch("MR_pTcorr", &MR_pTcorr, "MR_pTcorr/D");
    outTree->Branch("gammaR", &gammaR, "gammaR/D");
    outTree->Branch("dPhiCM", &dPhiCM, "dPhiCM/D");
    outTree->Branch("shatR_bl", &shatR_bl, "shatR_bl/D");
    outTree->Branch("EB1", &EB1, "EB1/D");
    outTree->Branch("EB2", &EB2, "EB2/D");
    outTree->Branch("EL1", &EL1, "EL1/D");
    outTree->Branch("EL2", &EL2, "EL2/D");
    outTree->Branch("CosThetaB1", &CosThetaB1, "CosThetaB1/D");
    outTree->Branch("CosThetaB2", &CosThetaB2, "CosThetaB2/D");
    outTree->Branch("CosThetaL1", &CosThetaL1, "CosThetaL1/D");
    outTree->Branch("CosThetaL2", &CosThetaL2, "CosThetaL2/D");
    outTree->Branch("TopMass1", &TopMass1, "TopMass1/D");
    outTree->Branch("TopMass2", &TopMass2, "TopMass2/D");
    outTree->Branch("GluinoMass1", &GluinoMass1, "GluinoMass1/D");
    outTree->Branch("GluinoMass2", &GluinoMass2, "GluinoMass2/D");
    outTree->Branch("VisHemMass1", &VisHemMass1, "VisHemMass1/D");
    outTree->Branch("VisHemMass2", &VisHemMass2, "VisHemMass2/D");
    
    //Hybrid Razor Approach
    outTree->Branch("TotalHemMass1", &TotalHemMass1, "TotalHemMass1/D");
    outTree->Branch("TotalHemMass2", &TotalHemMass2, "TotalHemMass2/D");
    outTree->Branch("TopHemMass1", &TopHemMass1, "TopHemMass1/D");
    outTree->Branch("TopHemMass2", &TopHemMass2, "TopHemMass2/D");
    outTree->Branch("MR1", &MR1, "MR1/D");
    outTree->Branch("MR2", &MR2, "MR2/D");
    outTree->Branch("TotalNMag", &TotalNMag, "TotalNMag/D");
    //pT-corrected (Rogan) Razor Approach
    outTree->Branch("TotalHemMass1Raz", &TotalHemMass1Raz, "TotalHemMass1Raz/D");
    outTree->Branch("TotalHemMass2Raz", &TotalHemMass2Raz, "TotalHemMass2Raz/D");
    outTree->Branch("TopHemMass1Raz", &TopHemMass1Raz, "TopHemMass1Raz/D");
    outTree->Branch("TopHemMass2Raz", &TopHemMass2Raz, "TopHemMass2Raz/D");
    outTree->Branch("MR1Raz", &MR1Raz, "MR1Raz/D");
    outTree->Branch("MR2Raz", &MR2Raz, "MR2Raz/D");
    outTree->Branch("TotalNMagRaz", &TotalNMagRaz, "TotalNMagRaz/D");
    //Hybrid Razor Approach + Transverse Boost
    outTree->Branch("TotalHemMass1Trans", &TotalHemMass1Trans, "TotalHemMass1Trans/D");
    outTree->Branch("TotalHemMass2Trans", &TotalHemMass2Trans, "TotalHemMass2Trans/D");
    outTree->Branch("TopHemMass1Trans", &TopHemMass1Trans, "TopHemMass1Trans/D");
    outTree->Branch("TopHemMass2Trans", &TopHemMass2Trans, "TopHemMass2Trans/D");
    outTree->Branch("MR1Trans", &MR1Trans, "MR1Trans/D");
    outTree->Branch("MR2Trans", &MR2Trans, "MR2Trans/D");
    outTree->Branch("TotalNMagTrans", &TotalNMagTrans, "TotalNMagTrans/D");
    //Canonical Razor Approach
    outTree->Branch("TotalHemMass1Star", &TotalHemMass1Star, "TotalHemMass1Star/D");
    outTree->Branch("TotalHemMass2Star", &TotalHemMass2Star, "TotalHemMass2Star/D");
    outTree->Branch("TopHemMass1Star", &TopHemMass1Star, "TopHemMass1Star/D");
    outTree->Branch("TopHemMass2Star", &TopHemMass2Star, "TopHemMass2Star/D");
    outTree->Branch("MR1Star", &MR1Star, "MR1Star/D");
    outTree->Branch("MR2Star", &MR2Star, "MR2Star/D");
    outTree->Branch("TotalNMagStar", &TotalNMagStar, "TotalNMagStar/D");
    //
    outTree->Branch("MetMag", &MetMag, "MetMag/D");
    outTree->Branch("HT", &HT, "HT/D");
    outTree->Branch("MHT", &MHT, "MHT/D");
    outTree->Branch("MHT_x", &MHT_x, "MHT_x/D");
    outTree->Branch("MHT_y", &MHT_y, "MHT_y/D");
    
    
    // Gen-Level
    outTree->Branch("idMcL1", &idMcL1, "idMcL1/I");
    outTree->Branch("idMothMcL1", &idMothMcL1, "idMothMcL1/I");
    outTree->Branch("idGrandMothMcL1", &idGrandMothMcL1, "idGrandMothMcL1/I");
    //outTree->Branch("idMcB1", &idMcB1, "idMcB1/I");
    //outTree->Branch("idMothMcB1", &idMothMcB1, "idMothMcB1/I");
    //outTree->Branch("idGrandMothMcB1", &idGrandMothMcB1, "idGrandMothMcB1/I");    
    outTree->Branch("pT1", &pT1, "pT1/D");
    outTree->Branch("eta1", &eta1, "eta1/D");
    outTree->Branch("phi1", &phi1, "phi1/D");
    outTree->Branch("idMcL2", &idMcL2, "idMcL2/I");
    outTree->Branch("idMothMcL2", &idMothMcL2, "idMothMcL2/I");
    outTree->Branch("idGrandMothMcL2", &idGrandMothMcL2, "idGrandMothMcL2/I");
    //outTree->Branch("idMcB2", &idMcB2, "idMcB2/I");
    //outTree->Branch("idMothMcB2", &idMothMcB2, "idMothMcB2/I");
    //outTree->Branch("idGrandMothMcB2", &idGrandMothMcB2, "idGrandMothMcB2/I");
    outTree->Branch("pT2", &pT2, "pT2/D");
    outTree->Branch("eta2", &eta2, "eta2/D");
    outTree->Branch("phi2", &phi2, "phi2/D");
    outTree->Branch("nLepTopDecay",&nLepTopDecay,"nLepTopDecay/I");
    outTree->Branch("nNeutrino",&nNeutrino,"nNeutrino/I");
    outTree->Branch("pTNeutrinoMag",&pTNeutrinoMag,"pTNeutrinoMag/D");
    outTree->Branch("genH1Mass", &genH1Mass, "genH1Mass/D");
    outTree->Branch("genH2Mass", &genH2Mass, "genH2Mass/D");
    outTree->Branch("genGammaCM", &genGammaCM, "genGammaCM/D");
    outTree->Branch("genshat", &genshat, "genshat/D");
    outTree->Branch("genMDR", &genMDR, "genMDR/D");
    outTree->Branch("genMDR2", &genMDR2, "genMDR2/D");
    outTree->Branch("gpXL1", &gpXL1, "gpXL1/D");
    outTree->Branch("gpYL1", &gpYL1, "gpYL1/D");
    outTree->Branch("gpZL1", &gpZL1, "gpZL1/D");
    outTree->Branch("genergyL1", &genergyL1, "genergyL1/D");
    outTree->Branch("gpXL2", &gpXL2, "gpXL2/D");
    outTree->Branch("gpYL2", &gpYL2, "gpYL2/D");
    outTree->Branch("gpZL2", &gpZL2, "gpZL2/D");
    outTree->Branch("genergyL2", &genergyL2, "genergyL2/D");
    outTree->Branch("gpXB_t", &gpXB_t, "gpXB_t/D");
    outTree->Branch("gpYB_t", &gpYB_t, "gpYB_t/D");
    outTree->Branch("gpZB_t", &gpZB_t, "gpZB_t/D");
    outTree->Branch("genergyB_t", &genergyB_t, "genergyB_t/D");
    outTree->Branch("gpXB_tb", &gpXB_tb, "gpXB_tb/D");
    outTree->Branch("gpYB_tb", &gpYB_tb, "gpYB_tb/D");
    outTree->Branch("gpZB_tb", &gpZB_tb, "gpZB_tb/D");
    outTree->Branch("energyB_tb", &genergyB_tb, "genergyB_tb/D");
    //outTree->Branch("gdPhi_ll", &gdPhi_ll, "gdPhi_ll/D");

    //MC b's and l's
    outTree->Branch("pXL1", &pXL1, "pXL1/D");
    outTree->Branch("pYL1", &pYL1, "pYL1/D");
    outTree->Branch("pZL1", &pZL1, "pZL1/D");
    outTree->Branch("energyL1", &energyL1, "energyL1/D");
    outTree->Branch("pXL2", &pXL2, "pXL2/D");
    outTree->Branch("pYL2", &pYL2, "pYL2/D");
    outTree->Branch("pZL2", &pZL2, "pZL2/D");
    outTree->Branch("energyL2", &energyL2, "energyL2/D");
    outTree->Branch("pXB1", &pXB1, "pXB1/D");
    outTree->Branch("pYB1", &pYB1, "pYB1/D");
    outTree->Branch("pZB1", &pZB1, "pZB1/D");
    outTree->Branch("energyB1", &energyB1, "energyB1/D");
    outTree->Branch("pXB2", &pXB2, "pXB2/D");
    outTree->Branch("pYB2", &pYB2, "pYB2/D");
    outTree->Branch("pZB2", &pZB2, "pZB2/D");
    outTree->Branch("energyB2", &energyB2, "energyB2/D");
    outTree->Branch("dPhi_ll", &dPhi_ll, "dPhi_ll/D");
    outTree->Branch("dPhi_ll_R", &dPhi_ll_R, "dPhi_ll_R/D");
    outTree->Branch("dPhi_ll_R_2", &dPhi_ll_R_2, "dPhi_ll_R_2/D");
    outTree->Branch("dPhi_bl1", &dPhi_bl1, "dPhi_bl1/D");
    outTree->Branch("dPhi_bl2", &dPhi_bl2, "dPhi_bl2/D");
    outTree->Branch("dPhi_bl12", &dPhi_bl12, "dPhi_bl12/D");
    outTree->Branch("dPhi_bl21", &dPhi_bl21, "dPhi_bl21/D");
    outTree->Branch("dPhi_bl1_R", &dPhi_bl1_R, "dPhi_bl1_R/D");
    outTree->Branch("dPhi_bl2_R", &dPhi_bl2_R, "dPhi_bl2_R/D");
    outTree->Branch("dPhi_bl12_R", &dPhi_bl12_R, "dPhi_bl12_R/D");
    outTree->Branch("dPhi_bl21_R", &dPhi_bl21_R, "dPhi_bl21_R/D");
    outTree->Branch("CosThetaBL1", &CosThetaBL1, "CosThetaBL1/D");
    outTree->Branch("CosThetaBL2", &CosThetaBL2, "CosThetaBL2/D");
    outTree->Branch("dPhi_bl", &dPhi_bl, "dPhi_bl/D");
    outTree->Branch("dPhi_bl_R", &dPhi_bl_R, "dPhi_bl_R/D");

    //Selection
    outTree->Branch("nIsolatedPFJets", &nIsolatedPFJets, "nIsolatedPFJets/I");
    outTree->Branch("nPFJets", &nPFJets, "nPFJets/I");
    outTree->Branch("nBtag_lead4jets", &nBtag_lead4jets, "nBtag_lead4jets/I");
    outTree->Branch("nBtag_medium_lead4jets", &nBtag_medium_lead4jets, "nBtag_medium_lead4jets/I");
    outTree->Branch("nBtag_tight_lead4jets", &nBtag_tight_lead4jets, "nBtag_tight_lead4jets/I");
    outTree->Branch("nBtag_TCHPT_lead4jets", &nBtag_TCHPT_lead4jets, "nBtag_TCHPT_lead4jets/I");
    outTree->Branch("Mll", &Mll, "Mll/D");
    
    double Npassed_In = 0;
    double Npassed_PV = 0;
    //Jets
    double NpassedPF_4Jet = 0;
    //Leptons
    double Npassed_Lept=0;
    //B-tag
    double Npassed_1b=0;
    //Mll Cut
    double Npassed_Mll = 0;
    double Z_mass = 91.1876;
    
    double weightII = 1.;
    unsigned int lastLumi=0;
    unsigned int lastRun=0;
    
    std::vector<std::string> maskHLT_DoubleMu;
    //maskHLT_DoubleMu.push_back("HLT_Mu13_Mu8_v");
    maskHLT_DoubleMu.push_back("HLT_Mu17_Mu8_v");
    maskHLT_DoubleMu.push_back("HLT_Mu17_TkMu8_v");

    std::vector<std::string> maskHLT_DoubleEle;
    maskHLT_DoubleEle.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    
    std::vector<std::string> maskHLT_MuEle;
    maskHLT_MuEle.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
    maskHLT_MuEle.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

    std::vector<std::string> maskHLT_Mu17;
    maskHLT_Mu17.push_back("HLT_Mu17_v");

    std::vector<std::string> maskHLT_Mu12;
    maskHLT_Mu12.push_back("HLT_Mu12_v");

    std::vector<std::string> maskHLT_Ele17;
    maskHLT_Ele17.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

    std::vector<std::string> maskHLT_Ele8;
    maskHLT_Ele8.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

    Long64_t nbytes = 0;
    Long64_t nb = 0;
    cout << "Number of entries= " << stop << endl;
    for (Long64_t jentry=start;  jentry<stop;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;
        
        
        //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
        if(_isData) {
            setRequiredTriggers(maskHLT_DoubleMu); reloadTriggerMask(true); HLT_DoubleMu = hasPassedHLT();
            setRequiredTriggers(maskHLT_DoubleEle); reloadTriggerMask(true); HLT_DoubleEle = hasPassedHLT();
            setRequiredTriggers(maskHLT_MuEle); reloadTriggerMask(true); HLT_MuEle = hasPassedHLT();
            setRequiredTriggers(maskHLT_Mu17); reloadTriggerMask(true); HLT_Mu17 = hasPassedHLT();
            setRequiredTriggers(maskHLT_Mu12); reloadTriggerMask(true); HLT_Mu12 = hasPassedHLT();
            setRequiredTriggers(maskHLT_Ele17); reloadTriggerMask(true); HLT_Ele17 = hasPassedHLT();
            setRequiredTriggers(maskHLT_Ele8); reloadTriggerMask(true); HLT_Ele8 = hasPassedHLT();
        }
        
        //Good Run selection
        if (_isData && _goodRunLS && !isGoodRunLS()) {
            if ( lastRun != runNumber || lastLumi != lumiBlock) {
                lastRun = runNumber;
                lastLumi = lumiBlock;
                std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
            }
            continue;
        }
        
        if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
            lastRun = runNumber;
            lastLumi = lumiBlock;
            std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
        }
        
        // Kill bad events in Data
        if(_isData && FailFilters()) continue;
        if(_isData && isFlagged()) continue;
        
        Npassed_In += weightII;
        
        double m0=9999, m12=9999, mc=9999;
        
        
        
        if((!_isData) && _isSMS){
            std::vector<float> parameterPoint = ParseEvent();
            mg = parameterPoint[0];
            mchi = parameterPoint[1];
        }
        
        //HLT
        int passedHLT = HLT_DoubleMu + HLT_DoubleEle + HLT_MuEle + HLT_Mu17 + HLT_Mu12 + HLT_Ele17 + HLT_Ele8;
        if (_isData==true) {            
            if (passedHLT==0) continue;
        }

        
        // find highest-pT PV [replace with Sagar's code]
        int iPV = passPV();
        //if(iPV<0) continue;
        Npassed_PV += weightII;
        nPV = N_PV_EVENT;
        
        nIsolatedPFJets = 0;
        nPFJets = 0;
        vector <TLorentzVector> PFJet;
        vector <int> iPFJet;
        vector <TLorentzVector> IsolatedPFJet;
        vector <int> iIsolatedPFJet;
        bool badjet = false;
        
        
        for(int i = 0; i < nAK5PFNoPUJet; i++){
            TLorentzVector myJet;
            double px = pxAK5PFNoPUJet[i];
            double py = pyAK5PFNoPUJet[i];
            double pz = pzAK5PFNoPUJet[i];
            double E = sqrt(px*px+py*py+pz*pz);
            myJet.SetPxPyPzE(px,py,pz,E);
            
            if(myJet.Pt() > 30.0 && fabs(myJet.Eta()) < 3.0){
                if ( isLoosePFNoPUJetID(i) ) {
                    PFJet.push_back(myJet);
                    iPFJet.push_back(i);
                    nPFJets += 1;
                    // Examine if PFJet is isolated from leptons or not and Mll pass the Z mass threshold
                    
                } else {
                    PFJet.clear();
                    iPFJet.clear();
                    nPFJets = 0.;
                    break;
                }
            }
        }
        
        for (int i = 0; i < PFJet.size(); i++){
            TLorentzVector myJet = PFJet[i];
            int iJet = iPFJet[i];
            // Examine if PFJet is isolated from leptons or not and Mll pass the Z mass threshold
            bool isIsolatedJet = true;
            for(int j=0; j<nMuon; j++) {
                TLorentzVector thisMu(pxMuon[j], pyMuon[j], pzMuon[j], energyMuon[j]);
                if(isTightMuon(j) && (thisMu.Pt()>20.) && (myJet.DeltaR(thisMu)<=0.3)) isIsolatedJet = false;
            }
            for(int k=0; k<nEle; k++) {
                TLorentzVector thisEle(pxEle[k], pyEle[k], pzEle[k], energyEle[k]);
                if(isTightElectron(k) && (thisEle.Pt()>20.) && (myJet.DeltaR(thisEle)<=0.3)) isIsolatedJet = false;
            }
            if(isIsolatedJet) {
                IsolatedPFJet.push_back(myJet);
                iIsolatedPFJet.push_back(iJet);
                nIsolatedPFJets += 1;
            }
        }
        
        
        // jet ID
        if (nPFJets==0) continue;
        
        //4Jet
        /*
         if (isDeltaRIsolated) {
         //if(iIsolatedPFJet.size()<4) continue;
         NpassedPF_4Jet+=weightII;
         } else {
         //if(iPFJet.size() <4) continue;
         NpassedPF_4Jet+=weightII;
         }
         */
        
        // at least 2 jets
        if (nPFJets < 2) continue;
        
        int flag = 1;
        //bubble sort the jets by Pt()
        for (int i=0; (i < iIsolatedPFJet.size()) && flag; i++){
            TLorentzVector tempvector;
            int tempi;
            flag = 0;
            for (int j=0; j < (iIsolatedPFJet.size()-1); j++){
                if (IsolatedPFJet.at(j+1).Pt() > IsolatedPFJet.at(j).Pt()){
                    tempvector = IsolatedPFJet.at(j);
                    IsolatedPFJet.at(j) = IsolatedPFJet.at(j+1);
                    IsolatedPFJet.at(j+1) = tempvector;
                    
                    tempi = iIsolatedPFJet.at(j);
                    iIsolatedPFJet.at(j) = iIsolatedPFJet.at(j+1);
                    iIsolatedPFJet.at(j+1) = tempi;
                    flag=1;
                }
            }
        }
        
        // b from leading 4 isolated jets
        nBtag_lead4jets = 0;
        nBtag_medium_lead4jets = 0;
        nBtag_tight_lead4jets = 0;
        nBtag_TCHPT_lead4jets = 0;
        if (nIsolatedPFJets >=4) {
            for(int b=0; b< iIsolatedPFJet.size(); b++){
                int n=iIsolatedPFJet.at(b);
                if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n]) && b < 4) nBtag_lead4jets++;
                if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n]) && b < 4) nBtag_medium_lead4jets++;
                if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n]) && b < 4) nBtag_tight_lead4jets++;
                if(trackCountingHighPurBJetTagsAK5PFNoPUJet[n] > 3.41 && b < 4) nBtag_TCHPT_lead4jets++;
            }
        }
        
        //Create arrays with the b-discriminators of the jets
        nBtag = 0;
        nBtag_medium = 0;
        nBtag_tight = 0;
        nBtag_TCHPT = 0;
        for(int b=0; b< iIsolatedPFJet.size(); b++){
            int n=iIsolatedPFJet.at(b);
            if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n])) nBtag++;
            if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n])) nBtag_medium++;
            if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n])) nBtag_tight++;
            if(trackCountingHighPurBJetTagsAK5PFNoPUJet[n] > 3.41) nBtag_TCHPT++;
        }

        // 1b requirement
        //if(nBtag == 0) continue;
        Npassed_1b+=weightII;
        
        // Boxes
        vector<int> iMuTight;
        vector<TLorentzVector> MuTight;
        for(int i=0; i<nMuon; i++) {
            TLorentzVector thisMu(pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i]);
            if(isTightMuon(i) and thisMu.Pt() > 20.) {
                iMuTight.push_back(i);
                MuTight.push_back(thisMu);
            }
        }
        
        vector<int> iEleTight;
        vector<TLorentzVector> EleTight;
        for(int i=0; i<nEle; i++) {
            TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
            if(isTightElectron(i) && thisEle.Pt() > 20.) {
                iEleTight.push_back(i);
                EleTight.push_back(thisEle);
            }
        }
        
        //bubble sort the muons by Pt()
        flag = 1;
        if (iMuTight.size() > 1) {
            for (int i=0; (i < iMuTight.size()) && flag; i++){
                TLorentzVector tempvector;
                int tempi;
                flag = 0;
                for (int j=0; j < (iMuTight.size()-1); j++){
                    if (MuTight.at(j+1).Pt() > MuTight.at(j).Pt()){
                        tempvector = MuTight.at(j);
                        MuTight.at(j) = MuTight.at(j+1);
                        MuTight.at(j+1) = tempvector;
                        
                        tempi = iMuTight.at(j);
                        iMuTight.at(j) = iMuTight.at(j+1);
                        iMuTight.at(j+1) = tempi;
                        flag=1;
                    }
                }
            }
        }
        
        //bubble sort the electrons by Pt()
        flag = 1;
        if (iEleTight.size() > 1) {
            for (int i=0; (i < iEleTight.size()) && flag; i++){
                TLorentzVector tempvector;
                int tempi;
                flag = 0;
                for (int j=0; j < (iEleTight.size()-1); j++){
                    if (EleTight.at(j+1).Pt() > EleTight.at(j).Pt()){
                        tempvector = EleTight.at(j);
                        EleTight.at(j) = EleTight.at(j+1);
                        EleTight.at(j+1) = tempvector;
                        
                        tempi = iEleTight.at(j);
                        iEleTight.at(j) = iEleTight.at(j+1);
                        iEleTight.at(j+1) = tempi;
                        flag=1;
                    }
                }
            }
        }
        
        // Determine SS/OS for each dilepton events and veto on Z_mass window of width 10GeV
        int iLepton1;
        int iLepton2;
        
        vector<TLorentzVector> DiLepton;
        BOX_NUM = -99;
        if(iMuTight.size() >0 && iEleTight.size()>0){
            iLepton1 = iMuTight[0];
            iLepton2 = iEleTight[0];
            DiLepton.push_back(MuTight[0]);
            DiLepton.push_back(EleTight[0]);
            Mll = (MuTight[0]+EleTight[0]).M();
            ss = chargeMuon[iLepton1]*chargeEle[iLepton2];
            BOX_NUM = 0;
        } else if(iMuTight.size()>1) {
            iLepton1 = iMuTight[0];
            iLepton2 = iMuTight[1];
            DiLepton.push_back(MuTight[0]);
            DiLepton.push_back(MuTight[1]);
            Mll = (MuTight[0]+MuTight[1]).M();
            ss = chargeMuon[iLepton1]*chargeMuon[iLepton2];
            BOX_NUM = 1;
        } else if(iEleTight.size()>1) {
            iLepton1 = iEleTight[0];
            iLepton2 = iEleTight[1];
            DiLepton.push_back(EleTight[0]);
            DiLepton.push_back(EleTight[1]);
            Mll = (EleTight[0]+EleTight[1]).M();
            ss = chargeEle[iLepton1]*chargeEle[iLepton2];
            BOX_NUM = 2;
        }
        Npassed_Mll += weightII;
        
        // two good leptons
        if(BOX_NUM<0) continue;
        Npassed_Lept += weightII;
        
        
        // create vector of b-jets
        vector<TLorentzVector> tBJets;
        vector<float> tdiscBJets;
        vector<int> tiBJets;
        vector<TLorentzVector> otherJets;
        vector<int> iotherJets;
		
		
		
        for(int b=0; b < iIsolatedPFJet.size(); b++){
            int n = iIsolatedPFJet.at(b);
            float disc = combinedSecondaryVertexBJetTagsAK5PFNoPUJet[n];
            TLorentzVector tempvector;
            tempvector.SetPxPyPzE(pxAK5PFNoPUJet[n],pyAK5PFNoPUJet[n],pzAK5PFNoPUJet[n],energyAK5PFNoPUJet[n]);
            
            if (disc > 0.244) {
                tiBJets.push_back(n);
                tdiscBJets.push_back(disc);
                tBJets.push_back(tempvector);
            }
        }
        
        // bubble sort the b-jets by CSV b-tag discrimator values
        flag = 1;
        for(int i=0; (i < tiBJets.size()) && flag; i++){
            TLorentzVector tempvector;
            int tempi;
            float tempdisc;
            flag = 0;
            for (int j=0; j < (tiBJets.size()-1); j++){
                if (tdiscBJets.at(j+1) > tdiscBJets.at(j)){
                    tempvector  = tBJets.at(j);
                    tBJets.at(j) = tBJets.at(j+1);
                    tBJets.at(j+1) = tempvector;
                    
                    tempi = tiBJets.at(j);
                    tiBJets.at(j) = tiBJets.at(j+1);
                    tiBJets.at(j+1) = tempi;
                    
                    tempdisc = tdiscBJets.at(j);
                    tdiscBJets.at(j) = tdiscBJets.at(j+1);
                    tdiscBJets.at(j+1) = tempdisc;
                    flag = 1;
                }
            }
        }
        
		//here have a routine to just use best 2 jets if condition is set appropriately
		int UseBestBJets;
		UseBestBJets = 1;
		
		vector<TLorentzVector> BJets;
		vector<float> discBJets;
		vector<int> iBJets;
		
        if (UseBestBJets == 1){
            if (tBJets.size() > 0)BJets.push_back(tBJets.front());
            if (tBJets.size() > 1) BJets.push_back(tBJets[1]);
            if (tiBJets.size() > 0)iBJets.push_back(tiBJets.front());
            if (tiBJets.size() > 1) iBJets.push_back(tiBJets[1]);
            if (tdiscBJets.size() > 0)discBJets.push_back(tdiscBJets.front());
            if (tdiscBJets.size() > 1) discBJets.push_back(tdiscBJets[1]);         
        }
		else {
			BJets = tBJets;
			discBJets = tdiscBJets;
			iBJets = tiBJets;
		}
		
        
        // New Razor Variables
        // dummy values
        HT = -9999.;
        MHT_x = -9999.;
        MHT_y = -9999.;
        MHT = -9999.;
        MetMag = -9999.;
        EB1 = -9999.;
        EB2 = -9999.;
        EL1 = -9999.;
        EL2 = -9999.;
        dPhiCM = -9999.;
        gammaR = -9999.;
        TopMass1 = -9999.;
        TopMass2 = -9999.;
        CosThetaB1 = -9999.;
        CosThetaB2 = -9999.;
        CosThetaL1 = -9999.;
        CosThetaL2 = -9999.;
        GluinoMass1 = -9999.;
        GluinoMass2 = -9999.;
        VisHemMass1 = -9999.;
        VisHemMass2 = -9999.;
        
        //Hybrid Razor Approach
        TotalHemMass1 = -9999.;
        TotalHemMass2 = -9999.;
        TopHemMass1 = -9999.;
        TopHemMass2 = -9999.;
        MR1 = -9999.;
        MR2 = -9999.;
        TotalNMag = -9999.;
        //pT-corrected (Rogan) Razor Approach
        TotalHemMass1Raz = -9999.;
        TotalHemMass2Raz = -9999.;
        TopHemMass1Raz = -9999.;
        TopHemMass2Raz = -9999.;
        MR1Raz = -9999.;
        MR2Raz = -9999.;
        TotalNMagRaz = -9999.;
        //Canonical Razor Approach
        TotalHemMass1Star = -9999.;
        TotalHemMass2Star = -9999.;
        TopHemMass1Star = -9999.;
        TopHemMass2Star = -9999.;
        MR1Star = -9999.;
        MR2Star = -9999.;
        TotalNMagStar = -9999.;
        //Hybrid Razor Approach + Transverse Boost
        TotalNMagTrans = -9999.;
        TotalHemMass1Trans = -9999.;
        TotalHemMass2Trans = -9999.;
        TopHemMass1Trans = -9999.;
        TopHemMass2Trans = -9999.;
        MR1Trans = -9999.;
        MR2Trans = -9999.;
        TotalNMagTrans = -9999.;
		
	//BDisc initialization
        BDiscList[0] = -9999. ; 
	BDiscList[1] = -9999. ;
	
	//MR original using only b and l
	MR_og = -9999.;
	
	pT_CM = -9999.;

        //pT-corr MRT, RSQ
	MRT_bl = -9999.;
	RSQ_bl = -9999.;

	//corrected shat
	rat = -9999.;
	shat_corr = -9999.;
	shat_m = -9999.;

	// MC B and L values
	pXL1 = -9999.;
	pYL1 = -9999.;
	pZL1 = -9999.;
	energyL1 = -9999.;
	pXL2 = -9999.;
	pYL2 = -9999.;
	pZL2 = -9999.;
	energyL2 = -9999.;
	dPhi_ll = -9999.;
	dPhi_ll_R = -9999.;
	dPhi_ll_R_2 = -9999.;
	dPhi_bl1 = -9999.;
	dPhi_bl2 = -9999.;
	dPhi_bl12 = -9999.;
	dPhi_bl21 = -9999.;
	dPhi_bl1_R = -9999.;
	dPhi_bl2_R = -9999.;
	dPhi_bl12_R = -9999.;
	dPhi_bl21_R = -9999.;
	CosThetaBL1 = -9999.;
	CosThetaBL2 = -9999.;
	dPhi_bl = -9999.;
	dPhi_bl_R = -9999.;
	pXB1 = -9999.;
	pYB1 = -9999.;
	pZB1 = -9999.;
	energyB1 = -9999.;
	pXB2 = -9999.;
	pYB2 = -9999.;
	pZB2 = -9999.;
	energyB2 = -9999.;
		
	int j = 2;
	//Array var
	while (j >= 0){
	  TotalHemMass1List[j] = -9999.;
	  TotalHemMass2List[j] = -9999.;
	  TopHemMass1List[j] = -9999.;
	  TopHemMass2List[j] = -9999.;
	  MR1List[j] = -9999.;
	  MR2List[j] = -9999.;
	  pTNMagList[j] = -9999.;
	  CosThetaB1List[j] = -9999.;
	  CosThetaB2List[j] = -9999.;
	  CosThetaL1List[j] = -9999.;
	  CosThetaL2List[j] = -9999.;
	  gammaRList[j] = -9999.;
	  j = j-1;
        }
        //dummy var
        double thisjet_x;
        double thisjet_y;
        thisjet_x = 0;
        thisjet_y = 0;
        
        
        if (DiLepton.size()>=2 && BJets.size()>=2){
            //Start by using PF MET
            TVector3 MET(pxPFMet[2], pyPFMet[2], 0.);
            MetMag = MET.Mag();
            
            HT = 0.;
            MHT_x = 0.;
            MHT_y = 0.;
            
            
            TLorentzVector thisjet;
            
            
            for (int i=0; i < IsolatedPFJet.size(); i++) {
                
                thisjet = IsolatedPFJet[i];
                
                HT += sqrt(thisjet.X() * thisjet.X() + thisjet.Y() * thisjet.Y());
                MHT_x += thisjet.X();
                MHT_y += thisjet.Y();
            }
            MHT = sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
            
            //two leptons
            TLorentzVector L1 = DiLepton[0];
            TLorentzVector L2 = DiLepton[1];
            
            //and two b's
            TLorentzVector B1, B2;
            int iB1, iB2;
            //first, we choose the hemisphere pairing of the b's and leptons by minimizing invariant masses
            
            float smallestM2 = 9999999999.;
            for (int i=0; i < BJets.size(); i++) {
                TLorentzVector testB1 = BJets[i];
                for (int j=i+1; j < BJets.size(); j++) {
                    TLorentzVector testB2 = BJets[j];
                    float M2 =  std::min( (testB1+L1).M2() + (testB2+L2).M2(), (testB1+L2).M2() + (testB2+L1).M2() );
                    if (M2 <= smallestM2) {
                        smallestM2 = M2;
                        B1 = testB1;
                        iB1 = iBJets[i];
                        B2 = testB2;
                        iB2 = iBJets[j];
						BDiscList[0] = std::max( discBJets[i], discBJets[j] );
                        BDiscList[1] = std::min( discBJets[i], discBJets[j] );
                    }
                }
            }
            
            
            for (int i=0; i < IsolatedPFJet.size(); i++) {
                if (iIsolatedPFJet[i]!=iB1 && iIsolatedPFJet[i]!=iB2){
                    otherJets.push_back(IsolatedPFJet[i]);
                    iotherJets.push_back(iIsolatedPFJet[i]);
                }
            }
            
            if( (B1+L1).M2() + (B2+L2).M2() > (B1+L2).M2() + (B2+L1).M2() ){
                TLorentzVector temp = B1;
                B1 = B2;
                B2 = temp;
            }
            
	    pXL1 = L1.Px();
	    pYL1 = L1.Py();
	    pZL1 = L1.Pz();
	    energyL1 = L1.E();
	    pXL2 = L2.Px();
	    pYL2 = L2.Py();
	    pZL2 = L2.Pz();
	    energyL2 = L2.E();
	    dPhi_ll = atan2(pYL1, pXL1); 
	    dPhi_ll = dPhi_ll - atan2(pYL2, pXL2 );
	    dPhi_ll = fabs(dPhi_ll);
	    dPhi_ll = min(dPhi_ll, 2 * atan(1)*4 - dPhi_ll);
	    pXB1 = B1.Px();
	    pYB1 = B1.Py();
	    pZB1 = B1.Pz();
	    energyB1 = B1.E();
	    pXB2 = B2.Px();
	    pYB2 = B2.Py();
	    pZB2 = B2.Pz();
	    energyB2 = B2.E();
	    dPhi_bl1 = L1.Vect().DeltaPhi(B1.Vect());
	    dPhi_bl2 = L2.Vect().DeltaPhi(B2.Vect());
	    dPhi_bl12 = L1.Vect().DeltaPhi(B2.Vect());
	    dPhi_bl21 = L2.Vect().DeltaPhi(B1.Vect());
	    dPhi_bl = (L1.Vect()+B1.Vect()).DeltaPhi(L2.Vect()+B2.Vect());
            // HYBRID RAZOR APPROACH
            // Combine jets keeping the Tops (b+l) together and in separate hemispheres:
            vector<TLorentzVector> Tops;
            TLorentzVector b1, b2, lep1, lep2;
            TVector3 vBeta1, vBeta2;
            TLorentzVector Top1 = (B1+L1);
            TLorentzVector Top2 = (B2+L2);
            Tops.push_back(Top1);
            Tops.push_back(Top2);
            
            vector<TLorentzVector> H12  = CombineJetsTs(otherJets, Tops);
            TLorentzVector H1 = H12[0];
            TLorentzVector H2 = H12[1];
            
            // Now, we must perform longitudinal boost from lab frame to CMz frame
            // in order to make procedure invariant under longitundinal boosts
            TVector3 BL = H1.Vect()+H2.Vect();
            BL.SetX(0.0);
            BL.SetY(0.0);
            BL = (1./(H1.E()+H2.E()))*BL;
            
            // Boost to CMz frame
            Top1.Boost(-BL);
            Top2.Boost(-BL);
            H1.Boost(-BL);
            H2.Boost(-BL);
            
            VisHemMass1 = H1.M();
            VisHemMass2 = H2.M();
            
            // Now, we must boost to each of the respective hemisphere rest frames
            TVector3 vBETA = Boost_type1(H1,H2);
            gammaRList[0] = 1./sqrt(1.-vBETA.Mag2());
            
            H1.Boost(-vBETA);
            H2.Boost(vBETA);
            Top1.Boost(-vBETA);
            Top2.Boost(vBETA);
            
            TLorentzVector N1, N2;
            N1.SetXYZM(-H1.X(),-H1.Y(),-H1.Z(),0.0);
            N2.SetXYZM(-H2.X(),-H2.Y(),-H2.Z(),0.0);
            
            
            //TotalHemMass1 = (H1+N1).M();
            //TotalHemMass2 = (H2+N2).M();
            TotalHemMass1List[0] = (H1+N1).M();
            TotalHemMass2List[0] = (H2+N2).M();
            
            //TopHemMass1 = (Top1+N1).M();
            //TopHemMass2 = (Top2+N2).M();
            
            TopHemMass1List[0] = (Top1+N1).M();
            TopHemMass2List[0] = (Top2+N2).M();
            
            //MR1 = 2*H1.Vect().Mag();
            //MR2 = 2*H2.Vect().Mag();
            
            MR1List[0] = 2*H1.Vect().Mag();
            MR2List[0] = 2*H2.Vect().Mag();
            
            //Boost neutrinos back to CM frame
            N1.Boost(vBETA);
            N2.Boost(-vBETA);
            
            //Boost neutrinos back to lab frame
            N1.Boost(BL);
            N2.Boost(BL);
            
            //Find sum of pTs for neutrinos
            N1.SetZ(0.0);
            N2.SetZ(0.0);
            
            pTNMagList[0] = (N1+N2).Pt();
            
            b1.SetPxPyPzE(B1.X(), B1.Y(), B1.Z(), B1.E());
            b2.SetPxPyPzE(B2.X(), B2.Y(), B2.Z(), B2.E());
            b1.Boost(-BL);
            b2.Boost(-BL);
            b1.Boost(-vBETA);
            b2.Boost(vBETA);
            
            //calculate top helicity angles
            CosThetaB1List[0] = b1.Vect().Dot(vBETA)/(b1.P()*vBETA.Mag());
            CosThetaB2List[0] = b2.Vect().Dot(-vBETA)/(b2.P()*vBETA.Mag());
            
            //W helicities
            lep1.SetPxPyPzE(L1.X(), L1.Y(), L1.Z(), L1.E());
            lep2.SetPxPyPzE(L2.X(), L2.Y(), L2.Z(), L2.E());
            lep1.Boost(-BL);
            lep2.Boost(-BL);
            lep1.Boost(-vBETA);
            lep2.Boost(vBETA);
            
            vBeta1 = (-1./(lep1.E()+sqrt((lep1+b1).Vect().Mag2())))*b1.Vect();
            vBeta2 = (-1./(lep2.E()+sqrt((lep2+b2).Vect().Mag2())))*b2.Vect();
            lep1.Boost(-vBeta1);
            lep2.Boost(-vBeta2);
            
            //calculate W helicity angles
            CosThetaL1List[0] = lep1.Vect().Dot(vBeta1)/(lep1.P()*vBeta1.Mag());
            CosThetaL2List[0] = lep2.Vect().Dot(-vBeta2)/(lep2.P()*vBeta2.Mag());
            
            
            //CANONICAL RAZOR APPROACH
            
            
            // Combine jets keeping the Tops (b+l) together and in separate hemispheres:
            vector<TLorentzVector> TopsStar;
            Top1 = (B1+L1);
            Top2 = (B2+L2);
            TopsStar.push_back(Top1);
            TopsStar.push_back(Top2);
            
            H12  = CombineJetsTs(otherJets, TopsStar);
            H1 = H12[0];
            H2 = H12[1];
            
            // Now, we must perform longitudinal boost from lab frame to CMz frame
            // in order to make procedure invariant under longitundinal boosts
            BL = H1.Vect()+H2.Vect();
            BL.SetX(0.0);
            BL.SetY(0.0);
            BL = (1./(H1.E()+H2.E()))*BL;
            
            
            //Define transverse boost to Razor frame
            double BetaStarMag = CalcBetaMRStarMag(H1,H2);
            TLorentzVector H1Trans = H1;
            TLorentzVector H2Trans = H2;
            H1Trans.SetZ(0.0);
            H2Trans.SetZ(0.0);
            TLorentzVector HTransSum = H1Trans + H2Trans;
            TVector3 BetaStar = (BetaStarMag * HTransSum.Vect())*(1./HTransSum.Vect().Mag());
            gammaRList[2] = 1./sqrt(1.-BetaStar.Mag2());
            
            // Boost to CMz frame
            Top1.Boost(-BL);
            Top2.Boost(-BL);
            H1.Boost(-BL);
            H2.Boost(-BL);
            
            
            //continue as before
            //TVector3
            
            H1.Boost(-BetaStar);
            H2.Boost(BetaStar);
            Top1.Boost(-BetaStar);
            Top2.Boost(BetaStar);
            
            N1, N2;
            N1.SetXYZM(-H1.X(),-H1.Y(),-H1.Z(),0.0);
            N2.SetXYZM(-H2.X(),-H2.Y(),-H2.Z(),0.0);
            
            
            //TotalHemMass1Star = (H1+N1).M();
            //TotalHemMass2Star = (H2+N2).M();
            TotalHemMass1List[2] = (H1+N1).M();
            TotalHemMass2List[2] = (H2+N2).M();
            
            //TopHemMass1Star = (Top1+N1).M();
            //TopHemMass2Star = (Top2+N2).M();
            TopHemMass1List[2] = (Top1+N1).M();
            TopHemMass2List[2] = (Top2+N2).M();
            
            //MR1Star = 2*H1.E();
            //MR2Star = 2*H2.E();
            MR1List[2] = 2*H1.Vect().Mag();
            MR2List[2] = 2*H2.Vect().Mag();
            
            //Boost neutrinos back to CM frame
            N1.Boost(BetaStar);
            N2.Boost(-BetaStar);
            
            //Boost neutrinos back to lab frame
            N1.Boost(BL);
            N2.Boost(BL);
            
            //Find sum of pTs for neutrinos
            N1.SetZ(0.0);
            N2.SetZ(0.0);
            
            pTNMagList[2] = (N1+N2).Pt();
            
            
            b1.SetPxPyPzE(B1.X(), B1.Y(), B1.Z(), B1.E());
            b2.SetPxPyPzE(B2.X(), B2.Y(), B2.Z(), B2.E());
            b1.Boost(-BL);
            b2.Boost(-BL);
            b1.Boost(-BetaStar);
            b2.Boost(BetaStar);
            //calculate top helicity angles
            //CosThetaB1 = b1.Vect().Dot(vBETA)/(b1.P()*vBETA.Mag());
            //CosThetaB2 = b2.Vect().Dot(-vBETA)/(b2.P()*vBETA.Mag());
            CosThetaB1List[2] = b1.Vect().Dot(vBETA)/(b1.P()*vBETA.Mag());
            CosThetaB2List[2] = b2.Vect().Dot(-vBETA)/(b2.P()*vBETA.Mag());
            
            //W helicities
            lep1.SetPxPyPzE(L1.X(), L1.Y(), L1.Z(), L1.E());
            lep2.SetPxPyPzE(L2.X(), L2.Y(), L2.Z(), L2.E());
            lep1.Boost(-BL);
            lep2.Boost(-BL);
            lep1.Boost(-BetaStar);
            lep2.Boost(BetaStar);
            
            vBeta1 = (-1./(lep1.E()+sqrt((lep1+b1).Vect().Mag2())))*b1.Vect();
            vBeta2 = (-1./(lep2.E()+sqrt((lep2+b2).Vect().Mag2())))*b2.Vect();
            lep1.Boost(-vBeta1);
            lep2.Boost(-vBeta2);
            
            //calculate W helicity angles
            CosThetaL1List[2] = lep1.Vect().Dot(vBeta1)/(lep1.P()*vBeta1.Mag());
            CosThetaL2List[2] = lep2.Vect().Dot(-vBeta2)/(lep2.P()*vBeta2.Mag());
            
            
            
            // HYBRID RAZOR APPROACH + TRANSVERSE BOOST
            
            // Combine jets keeping the Tops (b+l) together and in separate hemispheres:
            /*
             vector<TLorentzVector> TopsTrans;
             Top1 = (B1+L1);
             Top2 = (B2+L2);
             
             TopsTrans.push_back(Top1);
             TopsTrans.push_back(Top2);
             
             H12  = CombineJetsTs(otherJets, TopsTrans);
             H1 = H12[0];
             H2 = H12[1];
             
             
             
             // Now, we must perform longitudinal boost from lab frame to CMz frame
             // in order to make procedure invariant under longitundinal boosts
             BL = H1.Vect()+H2.Vect();
             TVector3 BTr = H1.Vect() + H2.Vect();
             BL.SetX(0.0);
             BL.SetY(0.0);
             BL = (1./(H1.E()+H2.E()))*BL;
             BTr = (1./(H1.E()+H2.E()))*BTr;
             // Boost to CMz frame
             Top1.Boost(-BL);
             Top2.Boost(-BL);
             H1.Boost(-BL);
             H2.Boost(-BL);
             
             // Transverse Boosts go here
             
             BTr.SetZ(0.0);
             Top1.Boost(-BTr);
             Top2.Boost(-BTr);
             H1.Boost(-BTr);
             H2.Boost(-BTr);
             //
             
             
             //continue as before
             //TVector3
             vBETA = Boost_type1(H1,H2);
             
             H1.Boost(-vBETA);
             H2.Boost(vBETA);
             Top1.Boost(-vBETA);
             Top2.Boost(vBETA);
             
             N1, N2;
             N1.SetXYZM(-H1.X(),-H1.Y(),-H1.Z(),0.0);
             N2.SetXYZM(-H2.X(),-H2.Y(),-H2.Z(),0.0);
             
             
             TotalHemMass1Trans = (H1+N1).M();
             TotalHemMass2Trans = (H2+N2).M();
             
             TopHemMass1Trans = (Top1+N1).M();
             TopHemMass2Trans = (Top2+N2).M();
             
             MR1Trans = 2*H1.Vect().Mag();
             MR2Trans = 2*H2.Vect().Mag();
             
             N1.SetZ(0.0);
             N2.SetZ(0.0);
             
             TotalNMagTrans = (N1+N2).Pt();
             */
            
            
            
            
            
            
            // pT-CORRECTED RAZOR APPROACH
            // Rogan's Approach to fully leptonic ttbar
            
            // Hemispheres only b's and l's
            H1 = (B1+L1);
            H2 = (B2+L2);

	    // MR original using only b's and l's
	    MR_og = sqrt(pow(H1.E()+H2.E(),2)-pow(H1.Pz()+H2.Pz(),2));

	   
	    // temp variables used later
	    TLorentzVector H1_temp;
	    TLorentzVector H2_temp = H2;
	    H1_temp.SetPxPyPzE(H1.Px() + MET.Px(), H1.Py() + MET.Py(), H1.Pz(), H1.E());
	    double HemMass_temp = abs(H1_temp.M()) + abs(H2_temp.M());
	    rat = (HemMass_temp / (2 * 173.5));
	    H1_temp.SetE(0);
	    H1_temp.SetPz(0);
	    H2_temp.SetPz(0);
	    H2_temp.SetE(0);

	    pT_CM = abs((H1_temp+H2_temp).M());

            // Now, we must perform longitudinal boost from lab frame to CMz frame
            // in order to make procedure invariant under longitundinal boosts
            BL = H1.Vect()+H2.Vect();
            BL.SetX(0.0);
            BL.SetY(0.0);
            BL = (1./(H1.E()+H2.E()))*BL;
            
            // Boost to CMz frame
            B1.Boost(-BL);
            L1.Boost(-BL);
            B2.Boost(-BL);
            L2.Boost(-BL);
            H1.Boost(-BL);
	    H2.Boost(-BL);

            // Now, we will perform a transverse boost from the CMz frame to our
            // approximation of the CM rest frame
            
            shatR_bl = shatR(H1+H2,MET);
            TVector3 betaTR = BetaTR(H1+H2,MET);
	    MRT_bl = CalcMTR(H1,H2,MET);
	    RSQ_bl = (MRT_bl * MRT_bl) / ((shatR_bl) * (shatR_bl));

	    if (rat>1) shat_corr = shatR_bl * rat;
	    else shat_corr = shatR_bl;

            // Boost to ~CM frame
            B1.Boost(-betaTR);
            B2.Boost(-betaTR);
            L1.Boost(-betaTR);
            L2.Boost(-betaTR);
            H1.Boost(-betaTR);
            H2.Boost(-betaTR);
            
            for (int i=0; i < otherJets.size(); i++) {
                otherJets[i].Boost(-BL);
                otherJets[i].Boost(-betaTR);
            }
            
            
            //at this stage you can calculate a useful angle, the azimuthal angle between
            //the boost from the CMz to ~CM frame and the visible system of particles
            dPhiCM = (B1+B2+L1+L2).Vect().DeltaPhi(betaTR);
	    dPhi_ll_R_2 = L1.Vect().DeltaPhi(L2.Vect());
	    dPhi_bl1_R = L1.Vect().DeltaPhi(B1.Vect());
	    dPhi_bl2_R = L2.Vect().DeltaPhi(B2.Vect());
	    dPhi_bl12_R = L1.Vect().DeltaPhi(B2.Vect());
	    dPhi_bl21_R = L2.Vect().DeltaPhi(B1.Vect());
	    dPhi_bl_R = (L1.Vect()+B1.Vect()).DeltaPhi(L2.Vect()+B2.Vect());
	    // check to make sure deltaphi is right
	    dPhi_ll_R = atan2(L1.Py(), L1.Px()); 
	    dPhi_ll_R = dPhi_ll_R - atan2(L2.Py(), L2.Px() );
	    dPhi_ll_R = fabs(dPhi_ll_R);
	    dPhi_ll_R = min(dPhi_ll_R, 2 * atan(1)*4 - dPhi_ll_R);
            
            //Now, we must boost to each of the respective TR frames (or top rest frames in ttbar)
            vBETA = Boost_type1(H1,H2);
            B1.Boost(-vBETA);
            B2.Boost(vBETA);
            L1.Boost(-vBETA);
            L2.Boost(vBETA);
            H1.Boost(-vBETA);
            H2.Boost(vBETA);

            
            //Calculate neutrinos
            N1, N2;
            N1.SetXYZM(-H1.X(),-H1.Y(),-H1.Z(),0.0);
            N2.SetXYZM(-H2.X(),-H2.Y(),-H2.Z(),0.0);
            
            TotalHemMass1List[1] = (H1+N1).M();
            TotalHemMass2List[1] = (H2+N2).M();
            //TotalHemMass1Raz = (H1+N1).M();
            //TotalHemMass2Raz = (H2+N2).M();
            
            //TopHemMass1Raz = (Top1+N1).M();
            //TopHemMass2Raz = (Top2+N2).M();
            TopHemMass1List[1] = (Top1+N1).M();
            TopHemMass2List[1] = (Top2+N2).M();
            
            //MR1Raz = 2*H1.E();
            //MR2Raz = 2*H2.E();
            MR1List[1] = shatR_bl;
            MR2List[1] = shatR_bl;
            
            //Boost neutrinos back to CM frame
            N1.Boost(vBETA);
            N2.Boost(-vBETA);
            
            //Boost neutrinos back to CMz frame
            N1.Boost(betaTR);
            N2.Boost(betaTR);
            
            //Boost neutrinos back to lab frame
            N1.Boost(BL);
            N2.Boost(BL);
            
            //Find pT of neutrinos
            N1.SetZ(0.0);
            N2.SetZ(0.0);
            
            pTNMagList[1] = (N1+N2).Pt();
            
            //Useful variable is gamma associated with that boost, which corresponds to how far off threshold the ttbar event is
            //gammaR = 1./sqrt(1.-vBETA.Mag2());
            gammaRList[1] = 1./sqrt(1.-vBETA.Mag2());

	    if (MetMag/74.>1) shat_m = shatR_bl/gammaRList[1] * MetMag/74.;
	    else shat_m = shatR_bl/gammaRList[1];
            
            //In the respective top rest frames you can calculate the energy of the B's (sensitive to first mass splitting) and the helicity angle
            //energies
            EB1 = B1.E();
            EB2 = B2.E();
            
            //calculate top helicity angles
            //CosThetaB1 = B1.Vect().Dot(vBETA)/(B1.P()*vBETA.Mag());
            //CosThetaB2 = B2.Vect().Dot(-vBETA)/(B2.P()*vBETA.Mag());
	    CosThetaBL1 = (B1+L1).Vect().Dot(vBETA)/((B1+L1).P()*vBETA.Mag());
	    CosThetaBL2 = (B2+L2).Vect().Dot(vBETA)/((B2+L2).P()*vBETA.Mag());
            CosThetaB1List[1] = B1.Vect().Dot(vBETA)/(B1.P()*vBETA.Mag());
            CosThetaB2List[1] = B2.Vect().Dot(-vBETA)/(B2.P()*vBETA.Mag());
            
            //calculate top mass - energy of all objects in T-frame assuming massless neutrinos
            TopMass1 = B1.E() + L1.E() + sqrt((B1+L1).Vect().Mag2());
            TopMass2 = B2.E() + L2.E() + sqrt((B2+L2).Vect().Mag2());
            
            TLorentzVector TopJet1, TopJet2;
            TopJet1.SetXYZM(0.,0.,0.,TopMass1);
            TopJet2.SetXYZM(0.,0.,0.,TopMass2);
            
            // Boost back to CM frame
            TopJet1.Boost(vBETA);
            TopJet2.Boost(-vBETA);
            // Now COMBINE JETS on all jets in CM frame
            otherJets.push_back(TopJet1);
            otherJets.push_back(TopJet2);
            vector<TLorentzVector> CMHems  = CombineJets(otherJets);
            
            if (CMHems.size()>=2){
                TLorentzVector CMHem1 = CMHems[0];
                TLorentzVector CMHem2 = CMHems[1];
                GluinoMass1 = CMHem1.M();
                GluinoMass2 = CMHem2.M();
            }
            
            //Now, we must boost to each of the respective WR frames to evaluate the lepton energy
            vBeta1 = (-1./(L1.E()+sqrt((L1+B1).Vect().Mag2())))*B1.Vect();
            vBeta2 = (-1./(L2.E()+sqrt((L2+B2).Vect().Mag2())))*B2.Vect();
            L1.Boost(-vBeta1);
            L2.Boost(-vBeta2);
            double mygamma1 = 1./sqrt(1.-vBeta1.Mag2());
            double mygamma2 = 1./sqrt(1.-vBeta2.Mag2());
            
            //second mass splitting
            EL1 = L1.E();
            EL2 = L2.E();
            
            //calculate W helicity angles
            //CosThetaL1 = L1.Vect().Dot(vBeta1)/(L1.P()*vBeta1.Mag());
            //CosThetaL2 = L2.Vect().Dot(vBeta2)/(L2.P()*vBeta2.Mag());
            CosThetaL1List[1] = L1.Vect().Dot(vBeta1)/(L1.P()*vBeta1.Mag());
            CosThetaL2List[1] = L2.Vect().Dot(vBeta2)/(L2.P()*vBeta2.Mag());
        }
        
        // Classic Razor Variables
        // dummy values
        pTPFHem1 = -9999.;
        etaPFHem1 = -9999.;
        phiPFHem1 = -9999.;
        pTPFHem2 = -9999.;
        etaPFHem2 = -9999.;
        phiPFHem2 = -9999.;
        PFR = -99999.;
        MR = -99999.;
        MR_pTcorr = -99999.;
        
        // hemispheres
        vector<TLorentzVector> tmpJet = CombineJets(PFJet);
        
        if(tmpJet.size() >= 2) {
            TLorentzVector PFHem1 = tmpJet[0];
            TLorentzVector PFHem2 = tmpJet[1];
            // PFMET
            TVector3 MET(pxPFMet[2], pyPFMet[2], 0.);
            MRT = CalcMTR(PFHem1, PFHem2, MET);
            double variable = -999999.;
            double Rvariable = -999999.;
            variable = CalcGammaMRstar(PFHem1, PFHem2);
            if(variable >0) Rvariable = MRT/variable;
            
            // *NEW* pT-corrected version of MR
            MR_pTcorr = shatR(PFHem1+PFHem2,MET);
            
            // fill the R and hem part of the output tree
            pTPFHem1 = PFHem1.Pt();
            etaPFHem1 = PFHem1.Eta();
            phiPFHem1 = PFHem1.Phi();
            pTPFHem2 = PFHem2.Pt();
            etaPFHem2 = PFHem2.Eta();
            phiPFHem2 = PFHem2.Phi();
            PFR = Rvariable;
            RSQ = Rvariable*Rvariable;
            MR = variable;
        }
        
        
        
        // Classic Razor Variables with FastJet Hemispheres
        // int fjetsize = 9999;
        // double conesize = 0.5;
        // vector<Jet> fJet = FastJetAlgorithmForceTwo(PFJet,conesize,2.);
        // fjetsize = (int) fJet.size();
        // if (fjetsize!=2) cout << "FAILURE: fJet.size() = "  << fjetsize << endl;
        //  vector<TLorentzVector> fastJet;
        //  fastJet.push_back(fJet[0].Get4Vector());
        //  fastJet.push_back(fJet[1].Get4Vector());
        
        //  // fast-hemispheres
        //  if(fastJet.size() >= 2) {
        //    TLorentzVector FJHem1 = fastJet[0];
        //    TLorentzVector FJHem2 = fastJet[1];
        //    // PFMET
        //    TVector3 MET(pxPFMet[2], pyPFMet[2], 0.);
        //    FJMRT = CalcMTR(FJHem1, FJHem2, MET);
        //    double variable = -999999.;
        //    double Rvariable = -999999.;
        //    variable = CalcGammaMRstar(FJHem1, FJHem2);
        //    if(variable >0) Rvariable = MRT/variable;
        
        //    // fill the R and hem part of the output tree
        //    pTFJHem1 = FJHem1.Pt();
        //    etaFJHem1 = FJHem1.Eta();
        //    phiFJHem1 = FJHem1.Phi();
        //    pTFJHem2 = FJHem2.Pt();
        //    etaFJHem2 = FJHem2.Eta();
        //    phiFJHem2 = FJHem2.Phi();
        //    FJR = Rvariable;
        //    FJRSQ = Rvariable*Rvariable;
        //    FJMR = variable;
        //  }
        
        //Gen-Level
        pT1 = -9999.;
        eta1 = -9999.;
        phi1 = -9999.;
        pT2 = -9999.;
        eta2 = -9999.;
        phi2 = -9999.;
        idMcL1 = -9999.;
        idMcL2 = -9999.;

        gpXL1 = -9999.;
        gpYL1 = -9999.;
	gpZL1 = -9999.;
	genergyL1 = -9999.;
	gpXL2 = -9999.;
	gpYL2 = -9999.;
	gpZL2 = -9999.;
	genergyL2 = -9999.;
	gpXB_t = -9999.;
	gpYB_t = -9999.;
	gpZB_t = -9999.;
	genergyB_t = -9999.;
	gpXB_tb = -9999.;
	gpYB_tb = -9999.;
	gpZB_tb = -9999.;
	genergyB_tb = -9999.;
	genH1Mass = -9999.;
	genH2Mass = -9999.;
	genshat = -9999.;
	genGammaCM = -9999.;
	genMDR = -9999.;
	genMDR2 = -9999.;
        if(!_isData) {
            nLepTopDecay = 0;
            nNeutrino = 0;
            pTNeutrino.SetXYZT(0.,0.,0.,0.);
            pTNeutrinoMag = 0.;

            int iL1 = -99;
            int iL2 = -99;
            int iB_t = -99;
	    int iB_tb = -99;
	    int iLSP_t = -99;
	    int iLSP_tb = -99;
	    TLorentzVector genH1 = -99;
	    TLorentzVector genH2 = -99;
	    TLorentzVector genLSP_t = -99;
	    TLorentzVector genLSP_tb = -99;
	    TLorentzVector genB_t = -99;
	    TLorentzVector genB_tb = -99;
	    TLorentzVector genL1 = -99;
	    TLorentzVector genL2 = -99;
	    TLorentzVector genTop = -99;
	    TLorentzVector genTopBar = -99;
            
            double deltaRGen1;
            double deltaEGen1 = 999999999;
            double deltaRGen2;
            double deltaEGen2 = 999999999;


            
            for(int i=0; i<nMc; i++) {
                // Neutrinos
                if (statusMc[i] == 3 &&
                    (abs(idMc[i]) == 12 ||
                     abs(idMc[i]) == 14 ||
                     abs(idMc[i]) == 16 ||
                     abs(idMc[i]) == 18 )){
                        nNeutrino ++;
                        pXMc = pMc[i]*cos(phiMc[i])*sin(thetaMc[i]);
                        pYMc = pMc[i]*sin(phiMc[i])*sin(thetaMc[i]);
                        pTNeutrino.SetX(pXMc+pTNeutrino.X());
                        pTNeutrino.SetY(pYMc+pTNeutrino.Y());
                    }
                // TT final state
                if(abs(idMc[mothMc[i]]) == 24) {
                    if(abs(idMc[i]) >= 11 &&
                       abs(idMc[i]) <= 18) {
                        nLepTopDecay ++;
                    }
                }
       		
		// For ttbar b + top variables:
		
		// b associated with t
		if ((SMS_temp == "none")){
	 	if(idMc[mothMc[i]] == 6 &&
		   idMc[i] == 5){
		  iB_t = i;
		}
		// b associated with tbar
		if(idMc[mothMc[i]] == -6 &&
		   idMc[i] == -5){
		  iB_tb = i;
		}
		}

		// For T2bw b + stop variables:
		if ((SMS_temp == "T2bw")){
		  if(idMc[mothMc[i]] == 1000006 &&
		     idMc[i] == 5){
		    iB_t = i;
		  }
		  if(idMc[mothMc[i]] == -1000006 &&
		     idMc[i] == -5){
		    iB_tb = i;
		  }
		  if(abs(idMc[i]) == 1000022 && 
		     idMc[mothMc[mothMc[i]]] == 1000006){
		    iLSP_t = i;
		  }
		  if(abs(idMc[i]) == 1000022 && 
		     idMc[mothMc[mothMc[i]]] == -1000006){
		    iLSP_tb = i;
		  }

		}
				
				
				
		// For T2tt b + stop variables:
		if ((SMS_temp == "T2tt")){
		  if(idMc[mothMc[i]] == 6 &&
		     idMc[i] == 5){
		    iB_t = i;
		  }
		  // b associated with tbar
		  if(idMc[mothMc[i]] == -6 &&
		     idMc[i] == -5){
		    iB_tb = i;
		  }
		  if(abs(idMc[i]) == 1000022 && 
		     idMc[mothMc[i]] == 1000006){
		    iLSP_t = i;
		  }
		  if(abs(idMc[i]) == 1000022 && 
		     idMc[mothMc[i]] == -1000006){
		    iLSP_tb = i;
		  }

		}

                pTNeutrinoMag = pTNeutrino.Pt();
                // Delta R 0.2 window of selecton dileptons
                double tempDeltaEGen;
                deltaRGen1 = sqrt(pow(DiLepton[0].Eta()-etaMc[i],2) + pow(DiLepton[0].Phi()-phiMc[i],2));
                tempDeltaEGen = abs(DiLepton[0].E()-energyMc[i]);
                if(deltaRGen1<0.2 && statusMc[i]==3) {
                    if (tempDeltaEGen < deltaEGen1) {
                        iL1 = i;
                        deltaEGen1 = tempDeltaEGen;
                    }
                }
                deltaRGen2 = sqrt(pow(DiLepton[1].Eta()-etaMc[i],2) + pow(DiLepton[1].Phi()-phiMc[i],2));
                tempDeltaEGen = abs(DiLepton[1].E()-energyMc[i]);
                if(deltaRGen2<0.2 && statusMc[i]==3) {
                    if (tempDeltaEGen < deltaEGen2) {
                        iL2 = i;
                        deltaEGen2 = tempDeltaEGen;
                    }
                }
                
                
            }
            
            nLepTopDecay = nLepTopDecay/2;
            
            if(iL1>=0) {
                pT1 = pMc[iL1]*sin(thetaMc[iL1]);
                eta1 = etaMc[iL1];
                phi1 = phiMc[iL1];
                idMcL1 = idMc[iL1];
                idMothMcL1 = idMc[mothMc[iL1]];
                if (SMS_temp == "none") idGrandMothMcL1 = idMc[mothMc[mothMc[iL1]]];
				
				// For T2bw, idGrandMothMcL1 stores great grandmother , stop
				if (SMS_temp == "T2bw") idGrandMothMcL1 = idMc[mothMc[mothMc[mothMc[iL1]]]];
				
				//For T2tt, idGrandMothMcL1 stores great grandmother , stop
				if (SMS_temp == "T2tt") idGrandMothMcL1 = idMc[mothMc[mothMc[mothMc[iL1]]]];
				
				gpXL1 = pMc[iL1]*cos(phiMc[iL1])*sin(thetaMc[iL1]);
				gpYL1 = pMc[iL1]*sin(phiMc[iL1])*sin(thetaMc[iL1]);
				gpZL1 = pMc[iL1]*cos(thetaMc[iL1]);
				genergyL1 = energyMc[iL1];
				genL1.SetPxPyPzE(gpXL1,gpYL1,gpZL1,genergyL1);
            }
            if(iL2>=0) {
                pT2 = pMc[iL2]*sin(thetaMc[iL2]);
                eta2 = etaMc[iL2];
                phi2 = phiMc[iL2];
                idMcL2 = idMc[iL2];
                idMothMcL2 = idMc[mothMc[iL2]];
                if (SMS_temp == "none") idGrandMothMcL2 = idMc[mothMc[mothMc[iL2]]];
                if (SMS_temp == "T2bw") idGrandMothMcL2 = idMc[mothMc[mothMc[mothMc[iL2]]]];
				if (SMS_temp == "T2tt") idGrandMothMcL1 = idMc[mothMc[mothMc[mothMc[iL1]]]];
		gpXL2 = pMc[iL2]*cos(phiMc[iL2])*sin(thetaMc[iL2]);
		gpYL2 = pMc[iL2]*sin(phiMc[iL2])*sin(thetaMc[iL2]);
		gpZL2 = pMc[iL2]*cos(thetaMc[iL2]);
		genergyL2 = energyMc[iL2];
		genL2.SetPxPyPzE(gpXL2,gpYL2,gpZL2,genergyL2);
            }
	    if(iB_t>=0) {
	      gpXB_t = pMc[iB_t]*cos(phiMc[iB_t])*sin(thetaMc[iB_t]);
	      gpYB_t = pMc[iB_t]*sin(phiMc[iB_t])*sin(thetaMc[iB_t]);
	      gpZB_t = pMc[iB_t]*cos(thetaMc[iB_t]);
	      genergyB_t = energyMc[iB_t];
	      genB_t.SetPxPyPzE(gpXB_t,gpYB_t,gpZB_t,genergyB_t);

	    }
	    if(iB_tb>=0) {
	      gpXB_tb = pMc[iB_tb]*cos(phiMc[iB_tb])*sin(thetaMc[iB_tb]);
	      gpYB_tb = pMc[iB_tb]*sin(phiMc[iB_tb])*sin(thetaMc[iB_tb]);
	      gpZB_tb = pMc[iB_tb]*cos(thetaMc[iB_tb]);
	      genergyB_tb = energyMc[iB_tb];
	      genB_tb.SetPxPyPzE(gpXB_tb,gpYB_tb,gpZB_tb,genergyB_tb);
	    }
	    double gpXLSP_t, gpYLSP_t, gpZLSP_t, genergyLSP_t;
	    double gpXLSP_tb, gpYLSP_tb, gpZLSP_tb, genergyLSP_tb;
	 
	    if(iLSP_t>=0) {
	      gpXLSP_t = pMc[iLSP_t]*cos(phiMc[iLSP_t])*sin(thetaMc[iLSP_t]);
	      gpYLSP_t = pMc[iLSP_t]*sin(phiMc[iLSP_t])*sin(thetaMc[iLSP_t]);
	      gpZLSP_t = pMc[iLSP_t]*cos(thetaMc[iLSP_t]);
	      genergyLSP_t = energyMc[iLSP_t];

	    }
	    if(iLSP_tb>=0) {
	      gpXLSP_tb = pMc[iLSP_tb]*cos(phiMc[iLSP_tb])*sin(thetaMc[iLSP_tb]);
	      gpYLSP_tb = pMc[iLSP_tb]*sin(phiMc[iLSP_tb])*sin(thetaMc[iLSP_tb]);
	      gpZLSP_tb = pMc[iLSP_tb]*cos(thetaMc[iLSP_tb]);
	      genergyLSP_tb = energyMc[iLSP_tb];
	    }
	    // for ttbar background, calculate four-vectors of two hemispheres with B's and L's
	    // for T2bw, calculate four-vectors of two hemispheres with B's and L's, ignoring chargino 
            if ((idGrandMothMcL1 == 6 && idGrandMothMcL2 == -6)||
		(idGrandMothMcL1 == 1000006 && idGrandMothMcL2 == -1000006)) {
	      genH1.SetPxPyPzE( gpXL1 + gpXB_t,
				gpYL1 + gpYB_t,
				gpZL1 + gpZB_t,
			        genergyL1 + genergyB_t);
	      genH2.SetPxPyPzE( gpXL2 + gpXB_tb,
				gpYL2 + gpYB_tb,
				gpZL2 + gpZB_tb,
			        genergyL2 + genergyB_tb);
	    }
	    else if ((idGrandMothMcL1 == -6 && idGrandMothMcL2 == 6)||
		     (idGrandMothMcL1 == -1000006 && idGrandMothMcL2 == 1000006)) {
	      genH1.SetPxPyPzE( gpXL2 + gpXB_t,
				gpYL2 + gpYB_t,
				gpZL2 + gpZB_t,
			        genergyL2 + genergyB_t );
	      genH2.SetPxPyPzE( gpXL1 + gpXB_tb,
				gpYL1 + gpYB_tb,
				gpZL1 + gpZB_tb,
			        genergyL1 + genergyB_tb );
	    }
	    // set genH1Mass to smaller of the two masses
	    if (abs(genH2.M()) >= abs(genH1.M())){
	      genH2Mass = abs(genH2.M());
	      genH1Mass = abs(genH1.M());
	    }
	    else {
	      genH1Mass = abs(genH2.M());
	      genH2Mass = abs(genH1.M());
	    }
	    // for ttbar background or T2bw, calculate gammaCM-gen level
	    // first find four-vectors for the top and antitop
			
		if (SMS_temp == "none" || SMS_temp == "T2bw"){
			
			
			if (iB_t>0){
				
				double pX_t = pMc[mothMc[iB_t]]*cos(phiMc[mothMc[iB_t]])*sin(thetaMc[mothMc[iB_t]]);
				double pY_t = pMc[mothMc[iB_t]]*sin(phiMc[mothMc[iB_t]])*sin(thetaMc[mothMc[iB_t]]);
				double pZ_t = pMc[mothMc[iB_t]]*cos(thetaMc[mothMc[iB_t]]);
				double energy_t = energyMc[mothMc[iB_t]];
				genTop.SetPxPyPzE( pX_t, pY_t, pZ_t, energy_t);
			}
			if (iB_tb>0){
				double pX_tb = pMc[mothMc[iB_tb]]*cos(phiMc[mothMc[iB_tb]])*sin(thetaMc[mothMc[iB_tb]]);
				double pY_tb = pMc[mothMc[iB_tb]]*sin(phiMc[mothMc[iB_tb]])*sin(thetaMc[mothMc[iB_tb]]);
				double pZ_tb = pMc[mothMc[iB_tb]]*cos(thetaMc[mothMc[iB_tb]]);
				double energy_tb = energyMc[mothMc[iB_tb]];
				genTopBar.SetPxPyPzE( pX_tb, pY_tb, pZ_tb, energy_tb);
			}
		}
			
			
		//Here we do the same for T2tt, using the grandmother stop particle instead though. 	
		else if (SMS_temp == "T2tt"){
			
			if (iB_t>0){
				
				double pX_t = pMc[mothMc[mothMc[iB_t]]]*cos(phiMc[mothMc[mothMc[iB_t]]])*sin(thetaMc[mothMc[mothMc[iB_t]]]);
				double pY_t = pMc[mothMc[mothMc[iB_t]]]*sin(phiMc[mothMc[mothMc[iB_t]]])*sin(thetaMc[mothMc[mothMc[iB_t]]]);
				double pZ_t = pMc[mothMc[mothMc[iB_t]]]*cos(thetaMc[mothMc[mothMc[iB_t]]]);
				double energy_t = energyMc[mothMc[mothMc[iB_t]]];
				genTop.SetPxPyPzE( pX_t, pY_t, pZ_t, energy_t);
			}
			if (iB_tb>0){
				double pX_tb = pMc[mothMc[mothMc[iB_tb]]]*cos(phiMc[mothMc[mothMc[iB_tb]]])*sin(thetaMc[mothMc[mothMc[iB_tb]]]);
				double pY_tb = pMc[mothMc[mothMc[iB_tb]]]*sin(phiMc[mothMc[mothMc[iB_tb]]])*sin(thetaMc[mothMc[mothMc[iB_tb]]]);
				double pZ_tb = pMc[mothMc[mothMc[iB_tb]]]*cos(thetaMc[mothMc[mothMc[iB_tb]]]);
				double energy_tb = energyMc[mothMc[mothMc[iB_tb]]];
				genTopBar.SetPxPyPzE( pX_tb, pY_tb, pZ_tb, energy_tb);
			}
			
			
			
		}
			
			
	    // boost to CM frame
	    TVector3 genBoostL = genTop.Vect() + genTopBar.Vect();
	    genBoostL = (1./(genTop.E()+genTopBar.E())) * genBoostL;
	    genTop.Boost(-genBoostL);
	    genTopBar.Boost(-genBoostL);
	    
	    // boost to top rest frames
	    TVector3 genBoostCM = genTop.Vect();
	    genBoostCM = (1./genTop.E()) * genBoostCM;
	    genTop.Boost(-genBoostCM);
	    genTopBar.Boost(genBoostCM);
	    // gammaCM and shat for gen-level
	    genGammaCM = 1./sqrt(1.-genBoostCM.Mag2());
	    genshat = sqrt(4 * genGammaCM * genGammaCM * genTop.M()* genTop.M());
	    if (SMS_temp == "none"){
	      genMDR = genTop.M();
	      genMDR2 = (genTop.M()*genTop.M() + (genH1Mass+genH2Mass)*(genH1Mass+genH2Mass))/(2.*genTop.M());
	    }
	    else if ((SMS_temp == "T2bw")||(SMS_temp == "T2tt")){
	      genMDR = (genTop.M()*genTop.M() - genLSP_t.M()*genLSP_t.M())/genTop.M();
	      genMDR2 = (genTop.M()*genTop.M() - genLSP_t.M()*genLSP_t.M() + ((genL1.M()+genL2.M()+genB_t.M()+genB_tb.M())*0.5)*(0.5*(genL1.M()+genL2.M()+genB_t.M()+genB_tb.M())))/(2.*genTop.M()) ;
	      
	    }
	    			     	      
        }
        
        
        // fill output tree
        run = runNumber;
        evNum = eventNumber;
        bx = eventNumber;
        ls = lumiBlock;
        orbit = orbitNumber;
        outTree->Fill();
    }
    
    // fill efficiency tree
    TTree* effTree = new TTree("effTree", "effTree");
    effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
    effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
    effTree->Branch("NpassedPF_4Jet",  &NpassedPF_4Jet,  "NpassedPF_4Jet/D");
    effTree->Branch("Npassed_1b",      &Npassed_1b,      "Npassed_1b/D");
    effTree->Branch("Npassed_Mll",     &Npassed_Mll,     "Npassed_Mll/D");
    effTree->Branch("Npassed_Lept",    &Npassed_Lept,    "Npassed_Lept/D");
    
    effTree->Fill();
    
    TFile *file = new TFile(outFileName.c_str(),"RECREATE");
    outTree->Write();
    effTree->Write();
    //  if(isSMS)FullSMSTree->Write();
    file->Close();
}

struct RazorMultiB::JetConfig{
    fastjet::JetDefinition theJetDef;
    fastjet::ActiveAreaSpec theAreaSpec;
};

vector<Jet> RazorMultiB::FastJetAlgorithmForceTwo(vector<TLorentzVector> InputCollection, double Rparam, double thePtMin){
    string JetFinder = "antikt_algorithm";
    string Strategy = "Best";
    double theDcut = -1;
    double theNjets = -1;
    string UE_subtraction = "no";
    bool theDoSubtraction = false;
    double theGhost_EtaMax = 6.0;
    int theActive_Area_Repeats = 5;
    double theGhostArea = 0.01;
    
    theJetConfig = new JetConfig;
    theJetConfig->theAreaSpec=fastjet::ActiveAreaSpec(theGhost_EtaMax, theActive_Area_Repeats, theGhostArea);
    
    fastjet::JetFinder jet_finder = fastjet::antikt_algorithm;
    
    fastjet::Strategy strategy = fastjet::Best;
    
    int theMode = 0;
    
    if((theNjets!=-1)&&(theDcut==-1)){
        theMode = 3;
    } else if((theNjets==-1)&&(theDcut!=-1)){
        theMode = 2;
    } else if((theNjets!=-1)&&(theDcut!=-1)){
        theMode = 1;
    } else {
        theMode = 0;
    }
    //
    theJetConfig->theJetDef = fastjet::JetDefinition(jet_finder, Rparam, strategy);
    
    std::vector<fastjet::PseudoJet> input_vectors;
    int index_ = 0;
    for(int i = 0; i < InputCollection.size(); i++){
        double px = InputCollection[i].Px();
        double py = InputCollection[i].Py();
        double pz = InputCollection[i].Pz();
        double E = InputCollection[i].E();
        fastjet::PseudoJet PsJet(px,py,pz,E);
        PsJet.set_user_index(index_);
        input_vectors.push_back(PsJet);
        index_++;
    }
    
    vector<Jet> output;
    if(index_ == 0) return output;
    
    std::vector<fastjet::PseudoJet> theJets;
    
    // running without subtraction; need to add code for with subtraction
    
    fastjet::ClusterSequence clust_seq(input_vectors, theJetConfig->theJetDef);
    
    if((theNjets==-1)&&(theDcut==-1)){
        theJets=clust_seq.inclusive_jets(thePtMin);
    } else if((theNjets!=-1)&&(theDcut==-1)){
        theJets=clust_seq.exclusive_jets(theNjets);
    } else if((theNjets==-1)&&(theDcut!=-1)){
        theJets=clust_seq.exclusive_jets(theDcut);
    } else if((theNjets!=-1)&&(theDcut!=-1)){
        theJets=clust_seq.inclusive_jets(thePtMin);
    } else {
        theJets=clust_seq.inclusive_jets(thePtMin);
    }
    while (theJets.size()>2 and Rparam < 3.14*0.5) {
        theJetConfig->theJetDef = fastjet::JetDefinition(jet_finder, Rparam, strategy);
        fastjet::ClusterSequence clust_seq_temp(input_vectors, theJetConfig->theJetDef);
        theJets=clust_seq.inclusive_jets(thePtMin);
        Rparam += 0.01;
    }
    //cout << "Rparam = " << Rparam << endl;
    // setting up code shamelessly stolen from
    // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/RecoParticleFlow/PostProcessing/interface/FastJetAlgo.h?revision=1.2&view=markup
    
    std::vector<fastjet::PseudoJet> output_ = theJets;
    std::vector< std::vector< fastjet::PseudoJet > > ptrvector;
    for (int ct = 0; ct != (int) theJets.size(); ++ct ) {
        ptrvector.push_back(clust_seq.constituents(theJets[ct]));
    }
    int algorithm_ = 2;
    double cont =0;
    double distance_ = Rparam - 0.5;
    /////////// cycle to really force the production of two final jets, implemented for kt and anti-kt only
    while(output_.size() > 2 &&  cont < 100){
        // if (cont!=0) cout << "distance_ = " << distance_ << endl;
        double dimin = 1e09;
        if(algorithm_ == 2)
            dimin = 0;
        unsigned sel1 = 0;
        unsigned sel2 = 0;
        double px1 = 0.,py1 = 0.,pt1 = 0.;
        double phi1 = 0.,eta1 = 0.;
        double di = 0;
        for(int it = 0; it < output_.size(); it++){
            px1 = output_[it].px();
            py1 = output_[it].py();
            pt1 = sqrt( px1*px1 + py1*py1 );
            if(algorithm_ == 2){
                di = pow(1/pt1,2);
                if( di > dimin){
                    dimin = di;
                    sel1 = it;
                }
            }
        }
        double dijmin = 1e09;
        if(algorithm_ == 2)
            dijmin = 0;
        phi1 = output_[sel1].phi();
        eta1 = output_[sel1].eta();
        pt1 = sqrt( output_[sel1].px()*output_[sel1].px() + output_[sel1].py()*output_[sel1].py() );
        for(int it = 0; it < output_.size(); it++){
            if(it != sel1){
                double px,py,pz,pt;
                px = output_[it].px();
                py = output_[it].py();
                pz = output_[it].pz();
                pt = sqrt(px*px+py*py);
                double phi = output_[it].phi();
                double eta = output_[it].eta();
                double r = sqrt( (phi1-phi)*(phi1-phi) + (eta1 -eta)*(eta1-eta) );
                if(algorithm_ == 2 ){
                    if( 1./(pt1*pt1) < 1./(pt*pt) )
                        di = 1./(pt1*pt1);
                    else
                        di = 1./(pt*pt);
                }
                double dij = di*r*r/pow(distance_,2);
                if(algorithm_ == 2){
                    if(dijmin < dij){
                        dijmin = dij;
                        sel2 = it;
                    }
                }
            }
        }
        distance_ = sqrt(dijmin*pow(distance_,2)/dimin);
        std::vector< fastjet::PseudoJet > temp_;
        std::vector< std::vector< fastjet::PseudoJet > > tempC_;
        std::vector< fastjet::PseudoJet > tempb_;
        temp_ = output_;
        output_.clear();
        tempC_ = ptrvector;
        ptrvector.clear();
        for(int it = 0; it < temp_.size(); it++){
            if( it != sel1 && it != sel2  ){
                output_.push_back(temp_[it]);
                ptrvector.push_back(tempC_[it]);
            }
            else if(it == sel1){
                fastjet::PseudoJet psj(temp_[sel1].px() + temp_[sel2].px(),
                                       temp_[sel1].py() + temp_[sel2].py(),
                                       temp_[sel1].pz() + temp_[sel2].pz(),
                                       temp_[sel1].e() + temp_[sel2].e() );
                unsigned it1;
                for(it1 = 0; it1 < tempC_[sel1].size() ; it1++){
                    tempb_.push_back(tempC_[sel1][it1]);
                }
                for(it1 =0; it1 < tempC_[sel2].size() ; it1++){
                    tempb_.push_back(tempC_[sel2][it1]);
                }
                ptrvector.push_back(tempb_);
                output_.push_back(psj);
            }
        }
        cont++;
        if(cont == 100)
            cout << "If this is printed out, you just tried to create a final state of two jets using an inclusive algorithm a hundred times, please check the logic, something is going on!" << endl;
    }
    
    
    
    //here, for the reco jets, need to loop through constituents to get fractions
    for(std::vector<fastjet::PseudoJet>::const_iterator itJet=output_.begin();
        itJet!=output_.end();itJet++){
        
        double px = (*itJet).px();
        double py = (*itJet).py();
        double pz = (*itJet).pz();
        double E = (*itJet).E();
        TLorentzVector J(px,py,pz,E);
        output.push_back(Jet(J,0.0,0.0));
    }
    
    delete theJetConfig;
    return output;
}


TVector3 RazorMultiB::Boost_type1(TLorentzVector H1, TLorentzVector H2){
	
	TVector3 vBETA = (1./(H1.E()+H2.E()))*(H1.Vect()-H2.Vect());
	
	return vBETA;
}

TVector3 RazorMultiB::BetaTR(TLorentzVector TOT, TVector3 MET){
	double myshatR = shatR(TOT,MET);
	TVector3 BCM = TOT.Vect()+MET;
	BCM.SetZ(0.0);
	BCM = (1./(sqrt(4.*myshatR*myshatR+BCM.Dot(BCM))))*BCM;
	
	return BCM;
}

double RazorMultiB::shat3D(TLorentzVector TOT, TVector3 P){
	double E = TOT.E();
	double Pz = 0.0;
	
	float MR = sqrt(E*E-Pz*Pz);
	
	
	TVector3 vI = P;
	
	TVector3 vpt = TOT.Vect();
	
	float MR2 = 2.*(MR*MR-vpt.Dot(vI)+MR*sqrt(MR*MR+vI.Dot(vI)-2.*vI.Dot(vpt)));
	
	return sqrt(MR2);
	
}	
double RazorMultiB::CalcBetaMRStarMag(TLorentzVector ja, TLorentzVector jb){
    double A = ja.P();
    double B = jb.P();
    double az = ja.Pz();
    double bz = jb.Pz();
    TVector3 jaT, jbT;
    jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
    jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
    double ATBT = (jaT+jbT).Mag2();
    
    double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                       (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());
    
    double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));
    
    return mybeta;
}
double RazorMultiB::shatR(TLorentzVector TOT, TVector3 MET){
	double E = TOT.E();
	double Pz = TOT.Pz();
	
	float MR = sqrt(E*E-Pz*Pz);
	
	
	TVector3 vI = MET+TOT.Vect();
	vI.SetZ(0.0);
	
	TVector3 vpt = TOT.Vect();
	vpt.SetZ(0.0);
	
	float MR2 = 0.5*(MR*MR-vpt.Dot(vI)+MR*sqrt(MR*MR+vI.Dot(vI)-2.*vI.Dot(vpt)));
	
	return sqrt(MR2);
	
}



std::vector<float> RazorMultiB::ParseEvent(){
	
    std::vector<std::string>::const_iterator c_begin = commentLHE->begin();
    std::vector<std::string>::const_iterator c_end = commentLHE->end();
    
    float mg, mchi;
    for( std::vector<std::string>::const_iterator cit=c_begin; cit!=c_end; ++cit) {
        size_t found = (*cit).find("model");
        if( found != std::string::npos)   {    
            size_t foundLength = (*cit).size();
            found = (*cit).find(" ");
            std::string smaller = (*cit).substr(found+1,foundLength);
            found = smaller.find("_");
            smaller = smaller.substr(found+1,smaller.size());
            //cout << smaller;
            
            std::istringstream iss(smaller);
            iss >> mg;
            iss.clear();
            //cout << mg << endl;
            
            found = smaller.find("_");
            smaller = smaller.substr(found+1,smaller.size());
            
            iss.str(smaller);
            iss >> mchi;
            iss.clear();
            //cout << mchi << endl;
        }
    }
    
    std::vector<float> parameterPoint;
    parameterPoint.push_back(mg);
    parameterPoint.push_back(mchi);
    return parameterPoint;
}

bool RazorMultiB::FailFilters(){
    bool FAIL = false;
    /*
     This is how METFlags is defined in Vecbos:
     http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsToWW2e/src/CmsMetFiller.cc?revision=1.13&view=markup
     
     *(privateData_->filterBits) = 
     (eeBadScFilterFlag << 8) | (hcalLaserEventFilterFlag << 7) | (HBHENoiseFilterResultFlag << 6) | 
     (isNotDeadEcalCluster << 5) | (trackerFailureFilterFlag << 4) | (CSCHaloFilterFlag << 3) | 
     ( drDead << 2 ) | ( drBoundary << 1 ) | ECALTPFilterFlag;
     */
    //cout << METFlags << endl;
    
    if((METFlags >> 0)%2 == 0) FAIL = true; //ecal dead cell tp
    //if((METFlags >> 1)%2 == 0) FAIL = true; // dr boundary
    //if((METFlags >> 2)%2 == 0) FAIL = true; // dr dead
    if((METFlags >> 3)%2 == 0) FAIL = true; // csc hal
    if((METFlags >> 4)%2 == 0) FAIL = true; // tracker failure
    //if((METFlags >> 5)%2 == 0) FAIL = true; // ecal dead cluster
    if((METFlags >> 6)%2 == 0) FAIL = true; // hbhe noise
    if((METFlags >> 7)%2 == 0) FAIL = true; // hcal laser
    if((METFlags >> 8)%2 == 0) FAIL = true; // bad ee sc
    //if((METFlags >> 9)%2 == 0) FAIL = true; //ecal laser
    
    return FAIL;
}


struct RazorMultiB::EventIndex {
	int RunNumber;
	Long64_t EventNumber;
	
	EventIndex() : RunNumber(0), EventNumber(0) {}
	
	bool operator <(const EventIndex &other) const
	{
		if(RunNumber < other.RunNumber)
			return true;
		if(RunNumber > other.RunNumber)
			return false;
		
		if(EventNumber < other.EventNumber)
			return true;
		if(EventNumber > other.EventNumber)
			return false;
		
		return false;
	}
};


void RazorMultiB::InitEventFlag(char *s_Event){
    ifstream inputfile(s_Event);
    
    int RUN_NUMBER;
    int LS_NUMBER;
    Long64_t EVENT_NUMBER;
    
    EventIndex index;
    
    cout << "Reading bad event list" << endl;
	
    while(!inputfile.eof()){
        inputfile >> RUN_NUMBER >> LS_NUMBER >> EVENT_NUMBER;
        
        index.RunNumber = RUN_NUMBER;
        index.EventNumber = EVENT_NUMBER;
        
        if(index.RunNumber < 0 || index.EventNumber < 0)
            continue;
        
        //Is this event/run-number combo alreading in the map?
        if(EventCounts.find(index) == EventCounts.end()){ //no
            EventCounts.insert(pair<EventIndex, int>(index, 1));
        } else { //yes
            EventCounts[index] = EventCounts[index] + 1;
        }
        index.RunNumber = -1;
        index.EventNumber = -1;
    }
}


bool RazorMultiB::isFlagged(){
    EventIndex index;
    index.EventNumber = eventNumber;
    index.RunNumber = runNumber;
    
    if(EventCounts.find(index) == EventCounts.end()){ //yes
        //cout << "Event not in list" << endl;
        return false;
    } else {
        //cout << "Event IS in list" << endl;
        return true;
        if(EventCounts[index] == 1){
            //cout << "Found the event - all is good in the world" << endl;
            return true;
        }
    } 
}



vector<TLorentzVector> RazorMultiB::CombineJetsTs(vector<TLorentzVector> myjets, vector<TLorentzVector> Ts){
    
    vector<TLorentzVector> mynewjets;
    TLorentzVector j1, j2;
    bool foundGood = false;
    
    int N_comb = 1;
    for(int i = 0; i < myjets.size(); i++){
        N_comb *= 2;
    }
    
    double M_min = 9999999999.0;
    int j_count;
    for(int i = 1; i < N_comb-1; i++){
        TLorentzVector j_temp1, j_temp2;
        int itemp = i;
        j_count = N_comb/2;
        int count = 0;
        while(j_count > 0){
            if(itemp/j_count == 1){
                j_temp1 += myjets[count];
            } else {
                j_temp2 += myjets[count];
            }
            itemp -= j_count*(itemp/j_count);
            j_count /= 2;
            count++;
        }    
        
        
        
        double M_temp = std::min((j_temp1 + Ts[0]).M2()+(j_temp2 + Ts[1]).M2() , (j_temp1 + Ts[1]).M2()+(j_temp2 + Ts[0]).M2() ) ;
        // smallest mass
        if (M_temp < M_min){
            M_min = M_temp;
            if ((j_temp1 + Ts[0]).M2()+(j_temp2 + Ts[1]).M2() == M_temp ){
                j1 = j_temp1 + Ts[0];
                j2 = j_temp2 + Ts[1];  
            }
            else{ 
                j1 = j_temp1 + Ts[1];
                j2 = j_temp2 + Ts[0];
            }
        }
    }
    
    // set masses to 0
    // j1.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),0.0);
    // j2.SetPtEtaPhiM(j2.Pt(),j2.Eta(),j2.Phi(),0.0);
    
    if(j2.Pt() > j1.Pt()){
        TLorentzVector temp = j1;
        j1 = j2;
        j2 = temp;
    }
    
    mynewjets.push_back(j1);
    mynewjets.push_back(j2);
    return mynewjets;  
}



