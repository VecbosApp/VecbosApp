// std includes
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream>
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
#include "RazorDMAnalysis.hh"

//Float_t Jet_Min_Pt = 70.0;
Float_t Jet_Min_Pt = 80.0;//at least 2 jest PT>80 GeV

RazorDMAnalysis::RazorDMAnalysis(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight=1.0;  
  _isSMS = false;
}

RazorDMAnalysis::RazorDMAnalysis(TTree *tree, string jsonFile,bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(jsonFile);
    fillRunLSMap();
  }
  
}

RazorDMAnalysis::~RazorDMAnalysis() {}

void RazorDMAnalysis::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorDMAnalysis::SetWeight(double weight){
  _weight=weight;
}

void RazorDMAnalysis::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;
  int pdgID = -99;
  if(outFileName.find("DYJetsHT200To400") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/DYJetsHT200To400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("QCD") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/DYJetsHT200To400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("DYJetsHT400") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/DYJetsHT400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("TTJetsFullyLeptMGDecays") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/TTj_Lep.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 6;
  }else if(outFileName.find("TTJetsSemiLeptMGDecays") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/TTj_Semilep.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 6;
  }else if(outFileName.find("TTJetsHadMGDecays") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/TTj_Had.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 6;
  }else if(outFileName.find("WJetsToLNu_150_HT_200") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/Wpj_150_HT_200.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("WJetsToLNu_200_HT_250") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/Wpj_200_HT_250.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("WJetsToLNu_250_HT_300") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/Wpj_250_HT_300.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("WJetsToLNu_300_HT_400") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/Wpj_300_HT_400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("WJetsToLNu_400_HT_Inf") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/Wpj_400_HT_inf.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("ZJetsToNuNu_50_HT_100") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/ZJetsToNuNu_50_HT_100.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("ZJetsToNuNu_100_HT_200") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/ZJetsToNuNu_100_HT_200.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("ZJetsToNuNu_200_HT_400") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/ZJetsToNuNu_200_HT_400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("ZJetsToNuNu_400_HT_inf") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/Files/ZJetsToNuNu_400_HT_inf.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("mDm") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/monoBt.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 18;
  }else if(outFileName.find("DMm") != string::npos){
    pu_f = new TFile("/afs/cern.ch/work/c/cpena/scratch_btagEff/CMSSW_5_2_3/src/PileUpCorrection/DMmTotal.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 18;
  }else{
    std::cout << "-----------INVALID MC NAME-------------------" << std::endl;
  }

  mu_corr_f = new TFile("/afs/cern.ch/user/w/woodson/public/WEIGHT/MuScaleFactorMap_MC53X_2012HZZ4L.root");
  mu_corr_h = (TH2F*)mu_corr_f->Get("TH2D_ALL_2012");

  //JEC Uncertainty
  /*
  JetCorrectionUncertainty *jec_un = 
    new JetCorrectionUncertainty(*(new JetCorrectorParameters("JEC_Uncertainty/Summer13_V5_MC_Uncertainty_AK5PF.txt", "Total")));
  */
  JetCorrectionUncertainty *jec_un =
    new JetCorrectionUncertainty(*(new JetCorrectorParameters("JEC_Uncertainty/Summer13_V5_MC_Uncertainty_AK5PF.txt")));
    
  TFile* file = new TFile(outFileName.c_str(),"RECREATE");

  double w_isr = 0.0;
  double w_isr_up = 0.0;
  double w_isr_down = 0.0;

  double w_pdf_CTEQ66[100];
  double w_pdf_MRST2006NNLO[100];
  double w_pdf_NNPDF10100[150];

  double w_pdf_CTEQ66_isr[100];
  double w_pdf_MRST2006NNLO_isr[100];
  double w_pdf_NNPDF10100_isr[150];

  double w_pdf_CTEQ66_isr_up = 0.0;
  double w_pdf_CTEQ66_isr_down = 0.0;
  double w_pdf_MRST2006NNLO_isr_up = 0.0;
  double w_pdf_MRST2006NNLO_isr_down = 0.0;
  double w_pdf_NNPDF10100_isr_up = 0.0;
  double w_pdf_NNPDF10100_isr_down = 0.0;

  double pu_w;
  double jes_up_w;
  double jes_down_w;
  double mu_w;
  double mu_w_up;
  double mu_w_down;
  double ISR;
  double ISR_up;
  double ISR_down;

  int HLT_Razor;
  int HLT_Razor_prescaled;
  int passedHLT;
  
  bool ECALTPFilterFlag;
  bool ecalLaserFilter;
  bool hcalLaserFilter;
  bool drBoundary;
  bool drDead;
  bool CSCHaloFilterFlag;
  bool trackerFailureFilterFlag;
  bool BEECALFlag; 
  bool eeBadScFilterFlag;//New filter, Cristian.(number 8 )
  bool HBHENoiseFilterResultFlag;//New filter, Cristian.(number 6 )

  // PF block
  double run;
  double evNum;
  double bx;
  double ls;
  double orbit;
  double R[4];
  double R_up[4];//JEC Uncertaninty
  double R_down[4];
  double RSQ[4];
  double RSQ_up[4];//JEC Uncertaninty
  double RSQ_down[4];
  double MR[4];
  double MR_up[4];// JEC Uncertainty
  double MR_down[4];
  double MRT[4];
  double MRT_up[4];//JEC Uncertainty
  double MRT_down[4];
  double pTHem1;
  double pTHem1_up;//JEC Uncertainty
  double pTHem1_down;
  double etaHem1;
  double etaHem1_up;//JEC Uncertainty
  double etaHem1_down;
  double phiHem1;
  double phiHem1_up;//JEC Uncertainty
  double phiHem1_down;
  double pTHem2;
  double pTHem2_up;//JEC Uncertainty
  double pTHem2_down;
  double etaHem2;
  double etaHem2_up;//JEC Uncertainty
  double etaHem2_down;
  double phiHem2;
  double phiHem2_up;//JEC Uncertainty
  double phiHem2_down;
  double Jet_PT[20];
  double Jet_PT_up[20];//JEC Uncertainty
  double Jet_PT_down[20];
  double Jet_Eta[20];
  double Jet_Eta_up[20];//JEC Uncertainty
  double Jet_Eta_down[20];
  double Jet_Phi[20];
  double Jet_Phi_up[20];//JEC Uncertainty
  double Jet_Phi_down[20];
  double CSV[20];
  int    nBtag;
  int    nBtagMed; 
  int    nBtagTight;
  double W = _weight;
  double mst, mchi;
  int    BOX_NUM;
  int BOX[3];
  Double_t Mu_Px_[2], Mu_Py_[2], Mu_Pz_[2], Mu_E_[2];
  
  double metX[4], metY[4], metCorrX[4], metCorrY[4], JetMetX, JetMetY, JetMetX_up, JetMetX_down, JetMetY_up, JetMetY_down, ht;
  double metX_up[4], metX_down[4], metY_up[4], metY_down[4], metCorrX_up[4], metCorrX_down[4], metCorrY_up[4], metCorrY_down[4];
  double mht[3];//xyz
  
  //Cristian MC information
  
  //Int_t           nMC;
  //Float_t         pMC[2001];   //[nMC]
  //Float_t         thetaMC[2001];   //[nMC]
  //Float_t         etaMC[2001];   //[nMC]
  //Float_t         phiMC[2001];   //[nMC]
  //Float_t         energyMC[2001];   //[nMC]
  //Int_t           idMC[2001];   //[nMC]
  //Int_t           mothMC[2001];   //[nMC]
  //Int_t           statusMC[2001];   //[nMC]
  
  float mssm[3];
  //int    ss;
  int    nPV;
  int    nPu;
  int Jet_Multiplicity;
  int N_Jets;
  // gen level info
  double pT1, pT2, eta1, eta2, phi1, phi2;
  int i1, i2;
  // ttbar decay: 0 = nolep, 1 = semilep; 2 = fully lep
  int iTopDecay;
  
  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("run", &run, "run/D");
  outTree->Branch("evNum", &evNum, "evNum/D");
  outTree->Branch("bx", &bx, "bx/D");
  outTree->Branch("ls", &ls, "ls/D");
  outTree->Branch("orbit", &orbit, "orbit/D");
  outTree->Branch("nPU", &nPu, "nPU/I");
  outTree->Branch("pu_w", &pu_w, "pu_w/D");
  outTree->Branch("mu_w", &mu_w, "mu_w/D");
  outTree->Branch("mu_w_up", &mu_w_up, "mu_w_up/D");
  outTree->Branch("mu_w_down", &mu_w_down, "mu_w_down/D");
  outTree->Branch("ISR", &ISR, "ISR/D");
  outTree->Branch("ISR_up", &ISR_up, "ISR_up/D");
  outTree->Branch("ISR_down", &ISR_down, "ISR_down/D");
  outTree->Branch("BOX_NUM", &BOX_NUM, "BOX_NUM/I");
  outTree->Branch("BOX", BOX, "BOX[3]/I");
  //outTree->Branch("ss", &ss, "ss/I");
  // HLT bits
  outTree->Branch("HLT_Razor", &HLT_Razor, "HLT_Razor/I"); 
  outTree->Branch("HLT_Razor_prescaled", &HLT_Razor_prescaled, "HLT_Razor_prescaled/I"); 
  outTree->Branch("passedHLT", &passedHLT, "passedHLT/I");
  //  block
  outTree->Branch("R", R, "R[4]/D");
  outTree->Branch("R_up", R_up, "R_up[4]/D");
  outTree->Branch("R_down", R_down, "R_down[4]/D");
  outTree->Branch("RSQ", RSQ, "RSQ[4]/D");
  outTree->Branch("RSQ_up", RSQ_up, "RSQ_up[4]/D");
  outTree->Branch("RSQ_down", RSQ_down, "RSQ_down[4]/D");
  outTree->Branch("MR", MR, "MR[4]/D");
  outTree->Branch("MR_up", MR_up, "MR_up[4]/D");
  outTree->Branch("MR_down", MR_down, "MR_down[4]/D");
  outTree->Branch("MRT", MRT, "MRT[4]/D");
  outTree->Branch("MRT_up", MRT_up, "MRT_up[4]/D");
  outTree->Branch("MRT_down", MRT_down, "MRT_down[4]/D");
  
  outTree->Branch("pTHem1", &pTHem1, "pTHem1/D");
  outTree->Branch("pTHem1_up", &pTHem1_up, "pTHem1_up/D");
  outTree->Branch("pTHem1_down", &pTHem1_down, "pTHem1_down/D");
  outTree->Branch("etaHem1", &etaHem1, "etaHem1/D");
  outTree->Branch("etaHem1_up", &etaHem1_up, "etaHem1_up/D");
  outTree->Branch("etaHem1_down", &etaHem1_down, "etaHem1_down/D");
  outTree->Branch("phiHem1", &phiHem1, "phiHem1/D");
  outTree->Branch("phiHem1_up", &phiHem1_up, "phiHem1_up/D");
  outTree->Branch("phiHem1_down", &phiHem1_down, "phiHem1_down/D");
  outTree->Branch("pTHem2", &pTHem2, "pTHem2/D");
  outTree->Branch("pTHem2_up", &pTHem2_up, "pTHem2_up/D");
  outTree->Branch("pTHem2_down", &pTHem2_down, "pTHem2_down/D");
  outTree->Branch("etaHem2", &etaHem2, "etaHem2/D");
  outTree->Branch("etaHem2_up", &etaHem2_up, "etaHem2_up/D");
  outTree->Branch("etaHem2_down", &etaHem2_down, "etaHem2_down/D");
  outTree->Branch("phiHem2", &phiHem2, "phiHem2/D");
  outTree->Branch("phiHem2_up", &phiHem2_up, "phiHem2_up/D");
  outTree->Branch("phiHem2_down", &phiHem2_down, "phiHem2_down/D");
  outTree->Branch("N_Jets", &N_Jets, "N_Jets/I");
  outTree->Branch("Jet_PT", Jet_PT, "Jet_PT[N_Jets]/D");
  outTree->Branch("Jet_PT_up", Jet_PT_up, "Jet_PT_up[N_Jets]/D");
  outTree->Branch("Jet_PT_down", Jet_PT_down, "Jet_PT_down[N_Jets]/D");
  outTree->Branch("Jet_Eta", Jet_Eta, "Jet_Eta[N_Jets]/D");
  outTree->Branch("Jet_Eta_up", Jet_Eta_up, "Jet_Eta_up[N_Jets]/D");
  outTree->Branch("Jet_Eta_down", Jet_Eta_down, "Jet_Eta_down[N_Jets]/D");
  outTree->Branch("Jet_Phi", Jet_Phi, "Jet_Phi[N_Jets]/D");
  outTree->Branch("Jet_Phi_up", Jet_Phi_up, "Jet_Phi_up[N_Jets]/D");
  outTree->Branch("Jet_Phi_down", Jet_Phi_down, "Jet_Phi_down[N_Jets]/D");
  outTree->Branch("CSV", CSV, "CSV[N_Jets]/D");
  outTree->Branch("nBtag", &nBtag, "nBtag/I");
  outTree->Branch("nBtagMed", &nBtagMed, "nBtagMed/I");
  outTree->Branch("nBtagTight", &nBtagTight, "nBtagTight/I");
  outTree->Branch("W", &W, "W/D");
  outTree->Branch("mst", &mst, "mst/D");
  outTree->Branch("mchi", &mchi, "mchi/D");
  outTree->Branch("nPV", &nPV, "nPV/I");
  outTree->Branch("i1", &i1, "i1/I");
  outTree->Branch("pT1", &pT1, "pT1/D");
  outTree->Branch("eta1", &eta1, "eta1/D");
  outTree->Branch("phi1", &phi1, "phi1/D");
  outTree->Branch("i2", &i2, "i2/I");
  outTree->Branch("pT2", &pT2, "pT2/D");
  outTree->Branch("eta2", &eta2, "eta2/D");
  outTree->Branch("phi2", &phi2, "phi2/D");
  //outTree->Branch("N_Jets", &N_Jets, "N_Jets/I");
  outTree->Branch("Mu_Px", Mu_Px_,"Mu_Px_[2]/D");
  outTree->Branch("Mu_Py", Mu_Py_,"Mu_Py_[2]/D");
  outTree->Branch("Mu_Pz", Mu_Pz_,"Mu_Pz_[2]/D");
  outTree->Branch("Mu_E", Mu_E_,"Mu_E_[2]/D");
  outTree->Branch("iTopDecay", &iTopDecay, "iTopDecay/I");
  outTree->Branch("mssm", mssm, "mssm[3]/F");
  //MET Info
  outTree->Branch("metX", metX, "metX[4]/D");
  outTree->Branch("metY", metY, "metY[4]/D");
  outTree->Branch("metCorrX", metCorrX, "metCorrX[4]/D");
  outTree->Branch("metCorrY", metCorrY, "metCorrY[4]/D");
  outTree->Branch("metX_up", metX_up, "metX_up[4]/D");
  outTree->Branch("metX_down", metX_down, "metX_down[4]/D");
  outTree->Branch("metY_up", metY_up, "metY_up[4]/D");
  outTree->Branch("metY_down", metY_down, "metY_down[4]/D");
  outTree->Branch("metCorrX_up", metCorrX_up, "metCorrX_up[4]/D");
  outTree->Branch("metCorrX_down", metCorrX_down, "metCorrX_down[4]/D");
  outTree->Branch("metCorrY_up", metCorrY_up, "metCorrY_up[4]/D");
  outTree->Branch("metCorrY_down", metCorrY_down, "metCorrY_down[4]/D");
  outTree->Branch("JetMetX", &JetMetX, "JetMetX/D");
  outTree->Branch("JetMetY", &JetMetY, "JetMetY/D");
  outTree->Branch("ht", &ht, "ht/D");
  outTree->Branch("mht", mht, "mht[3]/D");//xyz->[0,1,2]
  //MC GEN LEVEL INFO
  if( _isData == 0 ){
    outTree->Branch("nMC", &nMc, "nMC/I");
    outTree->Branch("pMC", pMc, "pMC[nMC]/F");
    outTree->Branch("thetaMC", thetaMc, "thetaMC[nMC]/F");
    outTree->Branch("etaMC", etaMc, "etaMC[nMC]/F");
    outTree->Branch("phiMC", phiMc, "phiMC[nMC]/F");
    outTree->Branch("energyMC", energyMc, "energyMC[nMC]/F");
    outTree->Branch("vxMC", vxMc, "vxMC[nMC]/F");
    outTree->Branch("vyMC", vyMc, "vyMC[nMC]/F");
    outTree->Branch("vzMC", vzMc, "vzMC[nMC]/F");
    outTree->Branch("idMC", idMc, "idMC[nMC]/I");
    outTree->Branch("mothMC", mothMc, "mothMC[nMC]/I");
    outTree->Branch("statusMC", statusMc, "statusMC[nMC]/I");
    if(outFileName.find("DMm") != string::npos || outFileName.find("mDm") != string::npos){
      outTree->Branch("nCTEQ66", &nCTEQ66, "nCTEQ66/I");
      outTree->Branch("wCTEQ66", wCTEQ66, "wCTEQ66[nCTEQ66]/D");
      outTree->Branch("nMRST2006NNLO", &nMRST2006NNLO, "nMRST2006NNLO/I");
      outTree->Branch("wMRST2006NNLO", wMRST2006NNLO, "wMRST2006NNLO[nMRST2006NNLO]/D");
      outTree->Branch("nNNPDF10100", &nNNPDF10100, "nNNPDF10100/I");
      outTree->Branch("wNNPDF10100", wNNPDF10100, "wNNPDF10100[nNNPDF10100]/D");
    }
  }
  
  double Npassed_In = 0;
  double Npassed_PV = 0;
  //Jets
  double Npassed_2Jet = 0;
  //Leptons
  double Npassed_LepVeto=0;
  //B-tag
  double Npassed_0btag=0;

  double weightII = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  
  std::vector<std::string> maskHLT_Razor; 
  maskHLT_Razor.push_back("HLT_RsqMR55_Rsq0p09_MR150");
  maskHLT_Razor.push_back("HLT_RsqMR60_Rsq0p09_MR150");
  maskHLT_Razor.push_back("HLT_RsqMR65_Rsq0p09_MR150");
  
  std::vector<std::string> maskHLT_Razor_prescaled; 
  maskHLT_Razor_prescaled.push_back("HLT_RsqMR40_Rsq0p04");
  //maskHLT_Razor_prescaled.push_back("HLT_RsqMR45_Rsq0p09");
    
  /*
  std::vector<std::string> maskHLT_Razor;                                         
  maskHLT_Razor.push_back("HLT_Mu17_Mu8");                                         
   
  std::vector<std::string> maskHLT_Razor_prescaled;                                 
  maskHLT_Razor_prescaled.push_back("HLT_Mu17_TkMu8");
  */

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Number of entries = " << stop << endl;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;

    //Cristian ISR correction
    TVector3 SYS(0.0, 0.0, 0.0);
    if(!_isData){
      TVector3 aux;
      for(int i=0; i<nMc; i++) {
        if(abs(idMc[i]) == pdgID && statusMc[i] == 3){
	  //std::cout << "PDG: " << idMc[i] << "  Pt: " << pMc[i]/cosh(etaMc[i]) << std::endl;
          aux.SetPtEtaPhi(pMc[i]/cosh(etaMc[i]), etaMc[i], phiMc[i]);
          SYS += aux;
        }
      }
    }
    //std::cout << "PT SYS: " << SYS.Pt() << std::endl;
    ISR = GetISR(SYS.Pt(), "central");
    ISR_up = GetISR(SYS.Pt(), "up");
    ISR_down = GetISR(SYS.Pt(), "down");
    
    w_isr += ISR;
    w_isr_up += ISR_up;
    w_isr_down += ISR_down;
    
    if(outFileName.find("DMm") != string::npos || outFileName.find("mDm") != string::npos){
      w_pdf_CTEQ66_isr_up += ISR_up*wCTEQ66[0];
      w_pdf_CTEQ66_isr_down += ISR_down*wCTEQ66[0];

      w_pdf_MRST2006NNLO_isr_up += ISR_up*wMRST2006NNLO[0];
      w_pdf_MRST2006NNLO_isr_down += ISR_down*wMRST2006NNLO[0];

      w_pdf_NNPDF10100_isr_up += ISR_up*wNNPDF10100[0];
      w_pdf_NNPDF10100_isr_down += ISR_down*wNNPDF10100[0];
      
      for(int l = 0; l < nCTEQ66; l++){
        w_pdf_CTEQ66[l] += wCTEQ66[l];
        w_pdf_CTEQ66_isr[l] += ISR*wCTEQ66[l];
      }
      for(int l = 0; l < nMRST2006NNLO; l++){
        w_pdf_MRST2006NNLO[l] += wMRST2006NNLO[l];
        w_pdf_MRST2006NNLO_isr[l] += ISR*wMRST2006NNLO[l];
      }
      for(int l = 0; l < nNNPDF10100; l ++){
        w_pdf_NNPDF10100[l] += wNNPDF10100[l];
        w_pdf_NNPDF10100_isr[l] += ISR*wNNPDF10100[l];
      }
    }

    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
      setRequiredTriggers(maskHLT_Razor); reloadTriggerMask(true); HLT_Razor = hasPassedHLT();
      setRequiredTriggers(maskHLT_Razor_prescaled); reloadTriggerMask(true); HLT_Razor_prescaled = hasPassedHLT();
      
      ECALTPFilterFlag = (METFlags >> 0)%2;
      drBoundary = (METFlags >> 1)%2;
      drDead = (METFlags >> 2)%2;
      CSCHaloFilterFlag = (METFlags >> 3)%2;
      trackerFailureFilterFlag = (METFlags >> 4)%2;
      BEECALFlag = (METFlags >> 5)%2; 
      HBHENoiseFilterResultFlag =  (METFlags >> 6)%2;
      ecalLaserFilter = (METFlags >> 7)%2;
      eeBadScFilterFlag = (METFlags >> 8)%2;
      hcalLaserFilter = (METFlags >> 9)%2;
    }
    
    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	lastRun = runNumber;
	lastLumi = lumiBlock;
	//std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    Npassed_In += weightII;
    
    double m0=9999, m12=9999, mc=9999;
    _isSMS = false;//Chage to true when running on SMS sample
    if( _isData == 0  && _isSMS ){
      int ctr_s = 0;
      //Read and store MSSM paramenters
      mssm[0] = mssm[1] = mssm[2] = 0;
      for (std::vector<string>::iterator it = commentLHE->begin() ; it != commentLHE->end(); ++it){
	istringstream iss (*it,istringstream::in);
	string val;
	for (int n = 0; n < 5; n++){
	  iss >> val;
	  if( n == 1){
	    if( val.compare( "model" ) )break;
	  }else if( n == 2 ){
	    cout << val.substr( 0,val.find("_") ) << endl;
	    string aa1 = val.substr( val.find("_")+ 1 );
	    string sm0 = aa1.substr( 0, val.find("_")-1 );
	    string sm1 = aa1.substr( val.find("_") );
	    std::cout << "sm0: " << sm0 << " sm1:  " << sm1 << std::endl;
	    mssm[0] = atof(sm0.c_str());
	    mssm[1] = atof(sm1.c_str());
	  }else if( n == 3 )mssm[2] = atof( val.c_str() );
	}
	ctr_s++;
      }
      
    }
    
    mst=m0;
    mchi=mc;
    
    //HLT and Data Filter
    //passedHLT = HLT_Razor + HLT_Razor_prescaled;
    passedHLT = HLT_Razor;

    if ( _isData == true ) {
      if ( passedHLT == 0 ) continue;//Comment out for getting trigger turn-ons
      if ((ECALTPFilterFlag==0) || (drBoundary==0) || (drDead==0) || (CSCHaloFilterFlag==0) || (trackerFailureFilterFlag==0) || (BEECALFlag==0) || ( HBHENoiseFilterResultFlag ==0 ) || (ecalLaserFilter == 0) || (eeBadScFilterFlag == 0) || (hcalLaserFilter == 0)) continue;
    }
    
    // find highest-pT PV [replace with Sagar's code]
    int iPV = passPV();
    if( iPV < 0 ) continue;//If negative no PV found to pass the cuts
    Npassed_PV += weightII;
    nPV = N_PV_EVENT;
    
    //////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////
    ////////////////////////// Calo Jets + JetID /////////////////////// 
    ////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////
    
    /*
    vector<TLorentzVector> Jet;
    vector <int> iJet;
    bool badjet = false;
    std::cout << "calo jets number " << nAK5Jet << std::endl; 
    for( int i = 0; i < nAK5Jet; i++ ) {
      TLorentzVector jet;
      double px = pxAK5Jet[i];
      double py = pyAK5Jet[i];
      double pz = pzAK5Jet[i];
      double E = sqrt( px*px + py*py + pz*pz );//massless Jet
      jet.SetPxPyPzE( px, py, pz, E );
      double pt = sqrt( px*px + py*py );
      if( fabs( jet.Eta() ) >= 3.0 ) continue;
      if( nHit90AK5Jet[i] <= 1 || fHPDAK5Jet[i] >= 0.98 ){//What are this variables for???
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) < 2.55 && emFracAK5Jet[i] <= 0.01 ){
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) >= 2.55 && jet.Pt() > 80. && emFracAK5Jet[i] >= 1. ){
	badjet = true;
	break;
      }
      if( fabs( jet.Eta() ) >= 2.55 && emFracAK5Jet[i] <= -0.9 ){
	badjet = true;
	break;
      }
      
      if ( jet.Pt() > 40. && fabs( jet.Eta() ) < 3.0 ) {
	Jet.push_back(jet);
	iJet.push_back(i);
      }
    }
    
    // Number of Jets                                                                                             
    Jet_Multiplicity[0] = Jet.size();
    
    // jet ID                                                                                                 
    if (badjet == true) continue;// If any Jet is bad (see loop before) event is rejected                    
    
    // >= 2Jets with pT> 70 GeV                                                       
    if( int( Jet.size() ) < 2 ) continue;//At least 2 Jets                                                     
    int iFirstJet = HighestPt(Jet, -99);
    int iSecondJet = HighestPt(Jet, iFirstJet);
    if( Jet[iSecondJet].Pt() < 70. ) continue;//First and second most energetic Jets Must have Pt > 70 GeV        
    Npassed_2Jet+=weightII;
    
    //count btagged jets                                                                                          
    nBtag[0] = 0;
    nBtagTight[0] = 0;
    
    for( int b = 0; b < iJet.size(); b++ ){
      if( pfJetPassCSVL( combinedSecondaryVertexBJetTagsAK5Jet[iJet[b]] ) ) nBtag[0]++; //Loose              
      if( pfJetPassCSVT( combinedSecondaryVertexBJetTagsAK5Jet[ iJet[b] ] ) ) nBtagTight[0]++; //Tight        
    }
    //if(nBtag>0) continue;                                                                                       
    Npassed_0btag+=weightII;//What is the meaning of this variable  
    
    /////////////////////////////                                                                                 
    ////////////HT///////////////                                                                                 
    /////////////////////////////                                                                            

    ht =  mht[0] = mht[1] = mht[2] = 0;//initialize HT and MHT                                                  
    for(  std::vector<TLorentzVector>::iterator it = Jet.begin(); it != Jet.end(); ++it){
      ht += (*it).Pt();//Compute ht, scalar sum of Pt
      mht[0] += (*it).Px();//Vectorial Sum of vec(P) of the jets                                                  
      mht[1] += (*it).Py();//Vectorial Sum of vec(P) of the jets                                                  
      mht[2] += (*it).Pz();//Vectorial Sum of vec(P) of the jets                                                  
    }

    */

    /////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //////////////////////// PF JETS + JetID ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    
    vector<TLorentzVector> pfJets;
    vector<int> i_pfJets;
    vector<double> pfJets_f_photon, pfJets_f_electron, pfJets_f_muon, pfJets_f_neutralhad, pfJets_f_chargedhad,	\
      pfJets_f_hfhad, pfJets_f_hfem;
    vector<double> pfJets_mult_photon, pfJets_mult_electron, pfJets_mult_muon, pfJets_mult_neutralhad, \
      pfJets_mult_chargedhad, pfJets_mult_hfhad, pfJets_mult_hfem;
    bool bad_pfJet = false;
    int N_pfJets = 0, pfBtags = 0;
    double pfHT, pfMHTx, pfMHTy;
    
    bool good_pfjet = false;
    for(int i = 0; i < nAK5PFNoPUJet; i++){
      TLorentzVector jet;
      double px = pxAK5PFNoPUJet[i];
      double py = pyAK5PFNoPUJet[i];
      double pz = pzAK5PFNoPUJet[i];
      double E = sqrt(px*px+py*py+pz*pz);
      double scale = 1.;
      jet.SetPxPyPzE(scale*px,scale*py,scale*pz,scale*E);
      
      good_pfjet = false;
      double EU = uncorrEnergyAK5PFNoPUJet[i];
      
      if(jet.Pt() > 40.0 && fabs(jet.Eta()) < 3.0){
	double fHAD = (neutralHadronEnergyAK5PFNoPUJet[i]+chargedHadronEnergyAK5PFNoPUJet[i])/EU;
	if(fHAD > 0.99){
	  N_pfJets = 0;
	  break;// clean NOISY event
	}
	
	int nConstituents = chargedHadronMultiplicityAK5PFNoPUJet[i]+neutralHadronMultiplicityAK5PFNoPUJet[i]+photonMultiplicityAK5PFNoPUJet[i]+electronMultiplicityAK5PFNoPUJet[i]+muonMultiplicityAK5PFNoPUJet[i]+HFHadronMultiplicityAK5PFNoPUJet[i]+HFEMMultiplicityAK5PFNoPUJet[i];
	int chargedMult = chargedHadronMultiplicityAK5PFNoPUJet[i]+electronMultiplicityAK5PFNoPUJet[i]+muonMultiplicityAK5PFNoPUJet[i];
	
	float photonFrac = photonEnergyAK5PFNoPUJet[i]/EU;
	float electronFrac = electronEnergyAK5PFNoPUJet[i]/EU;
	float muonFrac = muonEnergyAK5PFNoPUJet[i]/EU;
	float neutralHadFrac = neutralHadronEnergyAK5PFNoPUJet[i]/EU;
	float chargedHadFrac = chargedHadronEnergyAK5PFNoPUJet[i]/EU;
	float HFHadFrac = HFHadronEnergyAK5PFNoPUJet[i]/EU;
	float HFEMFrac = HFEMEnergyAK5PFNoPUJet[i]/EU;
	int photonMult = photonMultiplicityAK5PFNoPUJet[i];
	int electronMult = electronMultiplicityAK5PFNoPUJet[i];
	int muonMult = muonMultiplicityAK5PFNoPUJet[i];
	int neutralHadMult = neutralHadronMultiplicityAK5PFNoPUJet[i];
	int chargedHadMult = chargedHadronMultiplicityAK5PFNoPUJet[i];
	int HFHadMult = HFHadronMultiplicityAK5PFNoPUJet[i];
	int HFEMMult = HFEMMultiplicityAK5PFNoPUJet[i];
	
	if((neutralHadFrac < 0.99) && (photonFrac < 0.99) && (nConstituents > 1)) {
	  //outside of tracker acceptance, these are the only requirements
	  if (fabs(jet.Eta())>=2.4) good_pfjet = true;
	  //inside of the tracker acceptance, there are extra requirements     
	  else {
	    if ((chargedHadFrac > 0.0) && (chargedMult > 0) && (electronFrac < 0.99)) good_pfjet = true;
	  }
	}
	
	if(good_pfjet){
	  N_pfJets++;
	  pfJets.push_back(jet);
	  i_pfJets.push_back(i);
	  pfHT += jet.Pt();
	  pfMHTx -= jet.Px();
	  pfMHTy -= jet.Py();
	  //pfJets_btag.push_back(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[i]);
	  if( pfJetPassCSVL( combinedSecondaryVertexBJetTagsAK5PFNoPUJet[i] ) )pfBtags++;
	  pfJets_f_photon.push_back(photonFrac);
	  pfJets_f_electron.push_back(electronFrac);
	  pfJets_f_muon.push_back(muonFrac);
	  pfJets_f_neutralhad.push_back(neutralHadFrac);
	  pfJets_f_chargedhad.push_back(chargedHadFrac);
	  pfJets_f_hfhad.push_back(HFHadFrac);
	  pfJets_f_hfem.push_back(HFEMFrac);
	  pfJets_mult_photon.push_back(photonMult);
	  pfJets_mult_electron.push_back(electronMult);
	  pfJets_mult_muon.push_back(muonMult);
	  pfJets_mult_neutralhad.push_back(neutralHadMult);
	  pfJets_mult_chargedhad.push_back(chargedHadMult);
	  pfJets_mult_hfhad.push_back(HFHadMult);
	  pfJets_mult_hfem.push_back(HFEMMult);
	}else {
	  cout << "clean NOISY event JET ID" << endl;
	  N_pfJets = 0;
	  break;//Only takes out the pfJets loop! But good_jet = false
	}
      }
    }
    // jet ID                                                                 
    if (N_pfJets <= 0 )  continue;// If any Jet is bad (see loop before) event is rejected
    //std::cout << "JET REQUIREMENT" << std::endl;
    //////////////////////////////////////////////////////////////
    /////////////////////Create Muon Collection///////////////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    vector<int> iMuLoose;
    vector<TLorentzVector> MuLoose;
    vector<int> iMuLoose_StdIso;
    vector<TLorentzVector> MuLoose_StdIso;
    vector<int> iMuTight;
    vector<TLorentzVector> MuTight;

    for( int i = 0; i < nMuon; i++ ) {
      TLorentzVector thisMu( pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i] );
      /*if( ( isTightMuon(i, true) ) && ( thisMu.Pt() > 15. ) ) {
        iMuTight.push_back(i);
        MuTight.push_back(thisMu);
	iMuLoose.push_back(i);
        MuLoose.push_back(thisMu);
      }else if( ( isLooseMuon(i, true) ) && ( thisMu.Pt() > 15. ) ) {
        iMuLoose.push_back(i);
        MuLoose.push_back(thisMu);
	}*/

      if( ( isLooseMuon(i, true) ) && ( thisMu.Pt() > 15. ) ) {                                            
	iMuLoose.push_back(i);
	MuLoose.push_back(thisMu);
      }
      //Standard Isolation
      if(( isLooseMuon(i, false) ) && ( thisMu.Pt() > 15. )){
	iMuLoose_StdIso.push_back(i);
        MuLoose_StdIso.push_back(thisMu);
      }
    }
    
    mu_w = 1.0;
    double mu_err_sqr = 0.0;
    for(int j = 0; j < MuLoose.size(); j++){
      if(MuLoose[j].Pt() <= 100 && fabs(MuLoose[j].PseudoRapidity()) <= 2.4){
	mu_w *= mu_corr_h->GetBinContent( mu_corr_h->FindBin( MuLoose[j].Pt(), MuLoose[j].PseudoRapidity() ) );
	mu_err_sqr += pow(mu_corr_h->GetBinError( mu_corr_h->FindBin( MuLoose[j].Pt(), MuLoose[j].PseudoRapidity() ) ), 2);
      }else if( MuLoose[j].Pt() > 100 && fabs(MuLoose[j].PseudoRapidity()) <= 2.4 ){
	mu_w *= mu_corr_h->GetBinContent( mu_corr_h->FindBin( 99.9, MuLoose[j].PseudoRapidity() ) );
        mu_err_sqr += pow(mu_corr_h->GetBinError( mu_corr_h->FindBin( 99.9, MuLoose[j].PseudoRapidity() ) ), 2);
      }else if( MuLoose[j].Pt() <= 100 && fabs(MuLoose[j].PseudoRapidity()) > 2.4 ){
	int mu_eta;
	if(MuLoose[j].PseudoRapidity()>0){
	  mu_eta = 2.39;
	}else{
	  mu_eta = -2.39;
	}
        mu_w *= mu_corr_h->GetBinContent( mu_corr_h->FindBin( MuLoose[j].Pt(), mu_eta ) );
        mu_err_sqr += pow(mu_corr_h->GetBinError( mu_corr_h->FindBin( MuLoose[j].Pt(), mu_eta ) ), 2);
      }
    }
    
    mu_w_up = mu_w + sqrt(mu_err_sqr);
    mu_w_down = mu_w - sqrt(mu_err_sqr);
    
    ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    ///////////////// Create Collection  /////////////////    
    ////////////////  pfJets muon subtracted ////////////
    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    
    vector<TLorentzVector> pfJets_noMu;
    int ctr = 0;
    vector<int> i_pfJets_noMu;
    for(std::vector<TLorentzVector>::iterator pfJet_it = pfJets.begin(); pfJet_it != pfJets.end(); ++pfJet_it){
      bool isMuon = false;
      for (std::vector<TLorentzVector>::iterator Mu_it = MuLoose.begin(); Mu_it != MuLoose.end(); ++Mu_it){
	if ( (*Mu_it).Pt() < 15. ) continue;
	if( (*pfJet_it).DeltaR( *Mu_it ) <= 0.3 ){
	  isMuon = true;
	  break;
	}
      }
      if (!isMuon){
	pfJets_noMu.push_back( *pfJet_it );//If no muon found push it back into the collection  
	i_pfJets_noMu.push_back(i_pfJets[ctr]);
      }
      ctr++;
    }
    
    //JEC Up and Down
    vector<TLorentzVector> pfJets_noMu_Up;
    vector<TLorentzVector> pfJets_noMu_Down;
    for(int j = 0; j < pfJets_noMu.size(); j++){
      //JEC Up
      jec_un->setJetEta(pfJets_noMu[j].Eta());
      jec_un->setJetPt(pfJets_noMu[j].Pt());
      //std::cout << "Getting JEC"<< std::endl;
      double deltaPt = fabs(jec_un->getUncertainty(true));//Fractional Error
      //std::cout << "------deltaPt: " << pfJets_noMu[j].Pt()*deltaPt << std::endl;
      //std::cout << "Got JEC"<< std::endl;
      
      //JEC Up Now
      //double scl = (pfJets_noMu[j].Pt() + deltaPt)/(pfJets_noMu[j].Pt());
      double scl = 1.0 + deltaPt;
      TLorentzVector aux_J;
      aux_J.SetPxPyPzE(scl*pfJets_noMu[j].Px(), scl*pfJets_noMu[j].Py(), scl*pfJets_noMu[j].Pz(), scl*pfJets_noMu[j].E());
      pfJets_noMu_Up.push_back(aux_J);
      //JEC Down
      //scl = (pfJets_noMu[j].Pt() - deltaPt)/(pfJets_noMu[j].Pt());
      scl = 1.0 - deltaPt;
      aux_J.SetPxPyPzE(scl*pfJets_noMu[j].Px(), scl*pfJets_noMu[j].Py(), scl*pfJets_noMu[j].Pz(), scl*pfJets_noMu[j].E());
      pfJets_noMu_Down.push_back(aux_J);
    }
    //std::cout << "===>OUT OF JEC" << std::endl;
    // Number of Jets                                                                                             
    Jet_Multiplicity = pfJets_noMu.size();
    N_Jets = pfJets_noMu.size();

    //Create Maps with Pt as the Key to order the Jet from Highest to Lowest PT
    std::map<double, double> JetCVS_Map;
    
    std::map<double, double> JetEtaMap;
    std::map<double, double> JetEtaMap_up;
    std::map<double, double> JetEtaMap_down;
    std::map<double, double> JetPhiMap;
    std::map<double, double> JetPhiMap_up;
    std::map<double, double> JetPhiMap_down;
    
    std::vector<double> JetPT;
    std::map<double, double> JetPtMap_up;
    std::map<double, double> JetPtMap_down;
    
    //JetMet Init
    JetMetX = 0;
    JetMetX_up = 0;
    JetMetX_down = 0;
    JetMetY = 0;
    JetMetY_up = 0;
    JetMetY_down = 0;
    for(int j = 0; j < pfJets_noMu.size(); j++){
      JetMetX += pfJets_noMu[j].Px();
      JetMetX_up += pfJets_noMu_Up[j].Px();
      JetMetX_down += pfJets_noMu_Down[j].Px();
      JetMetY += pfJets_noMu[j].Py();
      JetMetY_up += pfJets_noMu_Up[j].Py();
      JetMetY_down += pfJets_noMu_Down[j].Py();
      
      
      JetCVS_Map[pfJets_noMu[j].Pt()] = combinedSecondaryVertexBJetTagsAK5PFNoPUJet[  i_pfJets_noMu[j]  ];
      
      JetEtaMap[pfJets_noMu[j].Pt()] = pfJets_noMu[j].Eta();
      JetEtaMap_up[pfJets_noMu[j].Pt()] = pfJets_noMu_Up[j].Eta();
      JetEtaMap_down[pfJets_noMu[j].Pt()] = pfJets_noMu_Down[j].Eta();
      
      JetPhiMap[pfJets_noMu[j].Pt()] = pfJets_noMu[j].Phi();
      JetPhiMap_up[pfJets_noMu[j].Pt()] = pfJets_noMu_Up[j].Phi();
      JetPhiMap_down[pfJets_noMu[j].Pt()] = pfJets_noMu_Down[j].Phi();
      
      JetPT.push_back( pfJets_noMu[j].Pt() );
      JetPtMap_up[pfJets_noMu[j].Pt()] = pfJets_noMu_Up[j].Pt();
      JetPtMap_down[pfJets_noMu[j].Pt()] = pfJets_noMu_Down[j].Pt();
    }
    //std::cout << "===>OUT OF MAPS" << std::endl;
    std::sort(JetPT.begin(), JetPT.end());
    std::reverse(JetPT.begin(), JetPT.end());
    for(int j = 0; j < pfJets_noMu.size(); j++){
      Jet_PT[j] = JetPT[j];
      Jet_PT_up[j] = JetPtMap_up[JetPT[j]];
      Jet_PT_down[j] = JetPtMap_down[JetPT[j]];
      
      Jet_Eta[j] = JetEtaMap[JetPT[j]];
      Jet_Eta_up[j] = JetEtaMap_up[JetPT[j]];
      Jet_Eta_down[j] = JetEtaMap_down[JetPT[j]];
      
      Jet_Phi[j] = JetPhiMap[JetPT[j]];
      Jet_Phi_up[j] = JetPhiMap_up[JetPT[j]];
      Jet_Phi_down[j] = JetPhiMap_down[JetPT[j]];
      
      CSV[j] = JetCVS_Map[JetPT[j]];
    }
    
    // >= 2Jets with pT > 80 GeV                                                            
    if( int( pfJets_noMu.size() ) < 2 ) continue;//At least 2 Jets                                         
    int iFirst_pfJet = HighestPt(pfJets_noMu, -99);
    int iSecond_pfJet = HighestPt(pfJets_noMu, iFirst_pfJet);
    if( pfJets_noMu[iSecond_pfJet].Pt() < Jet_Min_Pt ) continue;//First and second most energetic Jets
    //Must have Pt > 80 GeV  
    Npassed_2Jet+=weightII;
    //count btagged jets                                                                                          
    nBtag = 0;
    nBtagMed = 0;
    nBtagTight = 0;
        
    for( int b = 0; b < i_pfJets_noMu.size(); b++ ){
      //if( pfJetPassCSVL( combinedSecondaryVertexBJetTagsAK5Jet[ i_pfJets_noMu[b] ] ) ) nBtag++;//Loose
      //if( pfJetPassCSVM( combinedSecondaryVertexBJetTagsAK5Jet[ i_pfJets_noMu[b] ] ) ) nBtagMed++;//Med
      //if( pfJetPassCSVT( combinedSecondaryVertexBJetTagsAK5Jet[ i_pfJets_noMu[b] ] ) ) nBtagTight++;//Tight  
      
      if( pfJetPassCSVL( combinedSecondaryVertexBJetTagsAK5PFNoPUJet[ i_pfJets_noMu[b] ] ) ) nBtag++;//Loose
      if( pfJetPassCSVM( combinedSecondaryVertexBJetTagsAK5PFNoPUJet[ i_pfJets_noMu[b] ] ) ) nBtagMed++;//Med
      if( pfJetPassCSVT( combinedSecondaryVertexBJetTagsAK5PFNoPUJet[ i_pfJets_noMu[b] ] ) ) nBtagTight++;//Tight
    }
        
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    ///// Muon MET Correct pfJet_noMuon ///////
    ///////////////////////////////////////////
    //////////////////////////////////////////
    //Di-muon Transverse momentum
    TVector3 Muon_MET_Correction(0,0,0);
    // Boxes
    Mu_Px_[0] = Mu_Py_[0] = Mu_Pz_[0] = Mu_E_[0] = Mu_Px_[1] = Mu_Py_[1] = Mu_Pz_[1] = Mu_E_[1] =  -9999;

    // Correct MET for Presence of Loose Muons with pT>15 GeV
    int iFirstMuon = HighestPt(MuLoose, -99);
    if(iFirstMuon >= 0) {
      Muon_MET_Correction.SetX(Muon_MET_Correction.X() + MuLoose[iFirstMuon].Px());
      Muon_MET_Correction.SetY(Muon_MET_Correction.Y() + MuLoose[iFirstMuon].Py());

      Mu_Px_[0] =  MuLoose[iFirstMuon].Px();
      Mu_Py_[0] =  MuLoose[iFirstMuon].Py();
      Mu_Pz_[0] =  MuLoose[iFirstMuon].Pz();
      Mu_E_[0] =  MuLoose[iFirstMuon].E();
    }
    int iSecondMuon = HighestPt(MuLoose, iFirstMuon);
    if(iSecondMuon >= 0) {
      Muon_MET_Correction.SetX(Muon_MET_Correction.X() +  MuLoose[iSecondMuon].Px()); 
      Muon_MET_Correction.SetY(Muon_MET_Correction.Y() +  MuLoose[iSecondMuon].Py());
      
      Mu_Px_[1] = MuLoose[iSecondMuon].Px();
      Mu_Py_[1] = MuLoose[iSecondMuon].Py();
      Mu_Pz_[1] = MuLoose[iSecondMuon].Pz();
      Mu_E_[1] = MuLoose[iSecondMuon].E();
    }

    // Tight Electron ID
    vector<int> iEleTight;
    for(int i=0; i < nEle; i++) {
      TLorentzVector thisEle(pxEle[i], pyEle[i], pzEle[i], energyEle[i]);
      if(isTrigElectron(i) && thisEle.Pt() > 15.) {
	iEleTight.push_back(i);//Fill if the electron is tight and Pt > 15 GeV
      }
    }
    
    // Tight Tau ID
    
    vector<int> iTauTight;
    for ( int i = 0; i < nPFTau; i++) {
      TLorentzVector thisTau( pxPFTau[i], pyPFTau[i], pzPFTau[i], energyPFTau[i] );
      if ( isTightTau(i) && thisTau.Pt() > 20. ) {
	iTauTight.push_back(i);//Fill if the tau is tight and Pt > 20. GeV
      }
    }
    
    
    // Electron VETO
    if(iEleTight.size()>0) continue;//REMOVE ONLY FOR Trigger TURN_ON    
    //std::cout << "Pass ELe Veto" << std::endl;

    // TAU VETO
    //Using it for QCD sample
    if (iTauTight.size()>0) continue;//removed after a bug was found in the tau ID

    //////////////////////////////////////////////                                                                      
    //////////Indirect Lepton Veto (taus)/////////                                                                      
    //////////////////////////////////////////////  
    bool IsoPF = false;
    for(int kk = 0; kk < nPFCand; kk++){
      TLorentzVector pfCand4V(pxPFCand[kk], pyPFCand[kk], pzPFCand[kk], energyPFCand[kk]);
      bool isMuon = false;
      for(std::vector<TLorentzVector>::iterator Mu_it = MuLoose.begin(); Mu_it != MuLoose.end(); ++Mu_it){
	if( (*Mu_it).Pt() < 15. ) continue;
	if( (*Mu_it).DeltaR( pfCand4V ) <= 0.3){
	  isMuon = true;
          break;//out of the muon iterator loop
        }
      }
            
      if( ILV( kk ) < 0.15 && isMuon == false){
        IsoPF = true;
        break;//out of the PFCand loop. Ths means we found an isolated tau or electron 
      }
      
    }
    
    //Remove only for QCD sample
    //if(IsoPF)continue;//Applying ILV
    
    Npassed_LepVeto += weightII;//Record how many of th
    //std::cout << "===>OUT OF ILV" << std::endl;
        
    // BOX NUMBER
    BOX_NUM = -99;
    if(iMuLoose.size() == 0){
      BOX_NUM = 0;
    } else if(iMuLoose.size() == 1) {
      BOX_NUM = 1;
    } else if(iMuLoose.size() >= 2) {
      BOX_NUM = 2;
    }
    
    BOX[0] = BOX[1] = BOX[2] = -99;
    BOX[0] = iMuLoose.size();
    BOX[1] = iMuLoose_StdIso.size();
    BOX[2] = iMuTight.size();

    for(int l = 0; l < 4; l++){
      R[l] = -99999.;
      RSQ[l] = -99999.;
      MR[l] = -99999.;
      MRT[l] = -99999.;
    }
    
    pTHem1 = -9999.;                                                                                           
    etaHem1 = -9999.;                                                                                          
    phiHem1 = -9999.;                                                                                          
    pTHem2 = -9999.;                                                                                           
    etaHem2 = -9999.;                                                                                          
    phiHem2 = -9999.;
    
    // hemispheres PfJets Muons subtracted                                                                       
    //std::cout << pxPFMet[0] << " " << pxPFMet[1] << " " << pxPFMet[2] << " " << pxPFMet[3] << std::endl;       
    vector<TLorentzVector> hem_pfJet = CombineJets(pfJets_noMu);
    vector<TLorentzVector> hem_pfJet_Up = CombineJets(pfJets_noMu_Up);
    vector<TLorentzVector> hem_pfJet_Down = CombineJets(pfJets_noMu_Down);
    if( hem_pfJet.size() >= 2 ) {
      TLorentzVector Hem1 = hem_pfJet[0];
      TLorentzVector Hem2 = hem_pfJet[1];
      TLorentzVector Hem1_Up = hem_pfJet_Up[0];
      TLorentzVector Hem2_Up = hem_pfJet_Up[1];
      TLorentzVector Hem1_Down = hem_pfJet_Down[0];
      TLorentzVector Hem2_Down = hem_pfJet_Down[1];
      
      // PFMET + Correction                                                                                     
      TVector3 MET[4];
      TVector3 MET_up[4];
      TVector3 MET_down[4];
      for(int k = 0; k < 4; k++ ){
        MET[k] = TVector3(pxPFMet[k], pyPFMet[k], 0.);
	MET_up[k] = TVector3(pxPFMet[k]+JetMetX-JetMetX_up, pyPFMet[k]+JetMetY-JetMetY_up, 0.);
	MET_down[k] = TVector3(pxPFMet[k]+JetMetX-JetMetX_down, pyPFMet[k]+JetMetY-JetMetY_down, 0.);
	
        metX[k] = pxPFMet[k];
	metX_up[k] = pxPFMet[k]+JetMetX-JetMetX_up;
	metX_down[k] = pxPFMet[k]+JetMetX-JetMetX_down;
        metY[k] = pyPFMet[k];
	metY_up[k] = pyPFMet[k]+JetMetY-JetMetY_up;
	metY_down[k] = pyPFMet[k]+JetMetY-JetMetY_down;
	
	//MET with muon as neutrinos METCorr = -sum(Pt_pfCand)+Sum(Pt_pfCand==Muon)
	metCorrX[k] = pxPFMet[k] + Muon_MET_Correction.Px();
	metCorrX_up[k] = pxPFMet[k] + Muon_MET_Correction.Px() +JetMetX-JetMetX_up;
	metCorrX_down[k] = pxPFMet[k] + Muon_MET_Correction.Px() +JetMetX-JetMetX_down;
        metCorrY[k] = pyPFMet[k] + Muon_MET_Correction.Py();
	metCorrY_up[k] = pyPFMet[k] + Muon_MET_Correction.Px() +JetMetX-JetMetX_up;
	metCorrY_down[k] = pyPFMet[k] + Muon_MET_Correction.Py() +JetMetY-JetMetY_down;
	
        MRT[k] = CalcMTR(Hem1, Hem2, MET[k] + Muon_MET_Correction);
	MRT_up[k] = CalcMTR(Hem1_Up, Hem2_Up, MET_up[k] + Muon_MET_Correction);
	MRT_down[k] = CalcMTR(Hem1_Down, Hem2_Down, MET_down[k] + Muon_MET_Correction);
        //MRT = CalcMTR(Hem1, Hem2, MET);                                                                        
        double gammaMRstar = -999999.;
        double  r = -999999.;
        gammaMRstar = CalcGammaMRstar(Hem1, Hem2);
        if(gammaMRstar >0) r = MRT[k]/gammaMRstar;
	// fill the R and hem part of the output tree                                                           
	R[k] = r;
        RSQ[k] = r*r;
        MR[k] = gammaMRstar;
	//JEC Up
	gammaMRstar = CalcGammaMRstar(Hem1_Up, Hem2_Up);
        if(gammaMRstar >0) r = MRT_up[k]/gammaMRstar;
	R_up[k] = r;
        RSQ_up[k] = r*r;
        MR_up[k] = gammaMRstar;
	//JEC Up
	gammaMRstar = CalcGammaMRstar(Hem1_Down, Hem2_Down);
	if(gammaMRstar >0) r = MRT_down[k]/gammaMRstar;
        R_down[k] = r;
        RSQ_down[k] = r*r;
        MR_down[k] = gammaMRstar;
      }
      
      pTHem1 = Hem1.Pt();
      pTHem1_up = Hem1_Up.Pt();
      pTHem1_down = Hem1_Down.Pt();
      
      etaHem1 = Hem1.Eta();
      etaHem1_up = Hem1_Up.Eta();
      etaHem1_down = Hem1_Down.Eta();
      
      phiHem1 = Hem1.Phi();
      phiHem1_up = Hem1_Up.Phi();
      phiHem1_down = Hem1_Down.Phi();
      
      pTHem2 = Hem2.Pt();
      pTHem2_up = Hem2_Up.Pt();
      pTHem2_down = Hem2_Down.Pt();
      
      etaHem2 = Hem2.Eta();
      etaHem2_up = Hem2_Up.Eta();
      etaHem2_down = Hem2_Down.Eta();
      
      phiHem2 = Hem2.Phi();
      phiHem2_up = Hem2_Up.Phi();
      phiHem2_down = Hem2_Down.Phi();
    }
    
    //std::cout << "===>OUT OF MR and R2" << std::endl;
    
    //gen-level info
    pT1 = -999;
    eta1 = -999;
    phi1 = -999;     
    pT2 = -999;
    eta2 = -999;
    phi2 = -999;
    i1 = -99;
    i2 = -99;
    if(!_isData) {
      iTopDecay = 0;
      int iL1 = -99;
      int iL2 = -99;
      for(int i=0; i<nMc; i++) {
	// TT final state
	if(abs(idMc[mothMc[i]]) == 24) {
	  if(idMc[i] >= 11 &&
	     idMc[i] <= 18) {
	    iTopDecay ++;
	  }
	}
	// Z daughters
	if(idMc[mothMc[i]] == 23) {
	  if(iL1 <0) iL1 = i;
	  else if(iL2<0) iL2 = i;
	}
      }
      iTopDecay = iTopDecay/2;
      
      if(iL1>=0) {
	pT1 = pMc[iL1]*sin(thetaMc[iL1]);
	eta1 = etaMc[iL1];
	phi1 = phiMc[iL1];
	i1 = idMc[iL1];
      } 
      if(iL2>=0) {
	pT2 = pMc[iL2]*sin(thetaMc[iL2]);
	eta2 = etaMc[iL2];
	phi2 = phiMc[iL2];
	i2 = idMc[iL2];
      } 
    }

    // fill output tree
    run = runNumber;
    evNum = eventNumber;
    bx = eventNumber;
    ls = lumiBlock;
    orbit = orbitNumber;
    nPu = nPU[0];
    if(!_isData){
      pu_w = pu_h->GetBinContent(pu_h->FindBin(nPU[0]));
    }else{
      pu_w = 1.0;
      mu_w = 1.0;
    }
    outTree->Fill();
  }
  
  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
  effTree->Branch("Npassed_In", &Npassed_In, "Npassed_In/D");
  effTree->Branch("Npassed_ISR", &w_isr, "Npassed_ISR/D");
  effTree->Branch("Npassed_ISR_up", &w_isr_up, "Npassed_ISR_up/D");
  effTree->Branch("Npassed_ISR_down", &w_isr_down, "Npassed_ISR_down/D");

  effTree->Branch("nCTEQ66", &nCTEQ66, "nCTEQ66/I");
  effTree->Branch("N_pdf_CTEQ66", w_pdf_CTEQ66, "N_pdf_CTEQ66[nCTEQ66]/D");
  effTree->Branch("N_pdf_CTEQ66_isr", w_pdf_CTEQ66_isr, "N_pdf_CTEQ66_isr[nCTEQ66]/D");
  effTree->Branch("N_pdf_CTEQ66_isr_up", &w_pdf_CTEQ66_isr_up, "N_pdf_CTEQ66_isr_up/D");
  effTree->Branch("N_pdf_CTEQ66_isr_down", &w_pdf_CTEQ66_isr_down, "N_pdf_CTEQ66_isr_down/D");

  effTree->Branch("nMRST2006NNLO", &nMRST2006NNLO, "nMRST2006NNLO/I");
  effTree->Branch("N_pdf_MRST2006NNLO", w_pdf_MRST2006NNLO, "N_pdf_MRST2006NNLO[nMRST2006NNLO]/D");
  effTree->Branch("N_pdf_MRST2006NNLO_isr", w_pdf_MRST2006NNLO_isr, "N_pdf_MRST2006NNLO_isr[nMRST2006NNLO]/D");
  effTree->Branch("N_pdf_MRST2006NNLO_isr_up", &w_pdf_MRST2006NNLO_isr_up, "N_pdf_MRST2006NNLO_isr_up/D");
  effTree->Branch("N_pdf_MRST2006NNLO_isr_down", &w_pdf_MRST2006NNLO_isr_down, "N_pdf_MRST2006NNLO_isr_down/D");

  effTree->Branch("nNNPDF10100", &nNNPDF10100, "nNNPDF10100/I");
  effTree->Branch("N_pdf_NNPDF10100", w_pdf_NNPDF10100, "N_pdf_NNPDF10100[nNNPDF10100]/D");
  effTree->Branch("N_pdf_NNPDF10100_isr", w_pdf_NNPDF10100_isr, "N_pdf_NNPDF10100_isr[nNNPDF10100]/D");
  effTree->Branch("N_pdf_NNPDF10100_isr_up", &w_pdf_NNPDF10100_isr_up, "N_pdf_NNPDF10100_isr_up/D");
  effTree->Branch("N_pdf_NNPDF10100_isr_down", &w_pdf_NNPDF10100_isr_down, "N_pdf_NNPDF10100_isr_down/D");
  
  effTree->Branch("Npassed_PV", &Npassed_PV, "Npassed_PV/D");
  effTree->Branch("Npassed_2Jet", &Npassed_2Jet, "Npassed_2Jet/D");
  effTree->Branch("Npassed_0btag", &Npassed_0btag, "Npassed_0btag/D");
  effTree->Branch("Npassed_LepVeto", &Npassed_LepVeto, "Npassed_LepVeto/D");

  effTree->Fill();
  
  //TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  outTree->Write();
  effTree->Write();
  file->Close();
}

int RazorDMAnalysis::HighestPt(vector<TLorentzVector> p, int iHIGHEST) {
  
  int iH = -99;
  double highestPT = 0.;
  for(int i=0; i<p.size();i++) {
    if((p[i].Pt()>= highestPT) && (i != iHIGHEST)) {
      iH = i;
      highestPT = p[i].Pt();
    }
  }
  return iH;
}
