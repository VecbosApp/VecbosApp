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
#include "RazorRunTwo.hh"
#include "ControlSampleEvents.hh"

Float_t Jet_Min_Pt = 80.0;//at least 2 jest PT>80 GeV

const double ele_mass = 0.000511;
const int ele_pdgID = 11;
const double muon_mass = 0.1057;
const int muon_pdgID = 13;
const double tau_mass = 1.777;
const int etau_pdgID = 15;

//Default cone matching
const double genLepton_DR = 0.1;

RazorRunTwo::RazorRunTwo(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight=1.0;  
  _isSMS = false;
}

RazorRunTwo::RazorRunTwo(TTree *tree, string jsonFile,bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight=1.0;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(jsonFile);
    fillRunLSMap();
  }
  
}

RazorRunTwo::~RazorRunTwo() {}

void RazorRunTwo::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorRunTwo::SetWeight(double weight){
  _weight=weight;
}

void RazorRunTwo::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;
  int pdgID = -99;
  if(outFileName.find("DYJetsHT200To400") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/DYJetsHT200To400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("QCD") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/DYJetsHT200To400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("DYJetsHT400") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/DYJetsHT400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("TTJetsFullyLeptMGDecays") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/TTj_Lep.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 6;
  }else if(outFileName.find("TTJetsSemiLeptMGDecays") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/TTj_Semilep.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 6;
  }else if(outFileName.find("TTJetsHadMGDecays") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/TTj_Had.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 6;
  }else if(outFileName.find("WJetsToLNu_150_HT_200") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/Wpj_150_HT_200.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("WJetsToLNu_200_HT_250") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/Wpj_200_HT_250.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("WJetsToLNu_250_HT_300") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/Wpj_250_HT_300.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("WJetsToLNu_300_HT_400") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/Wpj_300_HT_400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("WJetsToLNu_400_HT_Inf") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/Wpj_400_HT_inf.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 24;
  }else if(outFileName.find("ZJetsToNuNu_50_HT_100") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/ZJetsToNuNu_50_HT_100.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("ZJetsToNuNu_100_HT_200") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/ZJetsToNuNu_100_HT_200.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("ZJetsToNuNu_200_HT_400") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/ZJetsToNuNu_200_HT_400.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("ZJetsToNuNu_400_HT_inf") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/ZJetsToNuNu_400_HT_inf.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 23;
  }else if(outFileName.find("mDm") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/monoBt.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 18;
  }else if(outFileName.find("DMm") != string::npos){
    pu_f = new TFile("/afs/cern.ch/user/c/cpena/public/PU_Files/DMmTotal.root");
    pu_h = (TH1D*)pu_f->Get("pileup");
    pdgID = 18;
  }else{
    std::cout << "-----------INVALID MC NAME-------------------" << std::endl;
  }

  mu_corr_f = new TFile("/afs/cern.ch/user/w/woodson/public/WEIGHT/MuScaleFactorMap_MC53X_2012HZZ4L.root");
  mu_corr_h = (TH2F*)mu_corr_f->Get("TH2D_ALL_2012");

  JetCorrectionUncertainty *jec_un =
    new JetCorrectionUncertainty(*(new JetCorrectorParameters("JEC_Uncertainty/Summer13_V5_MC_Uncertainty_AK5PF.txt")));
  
  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  
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
  
  events = new ControlSampleEvents;
  events->CreateTree();
  events->tree_->SetAutoFlush(0);

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
    
    //Filling Normalization Histogram
    NEvents->Fill(1.0);
    //Filling Event Info
    events->weight = 1.0;
    events->run = runNumber;
    events->lumi = lumiBlock;
    events->event = eventNumber;
    events->processID = 1;
    events->NPU_0 = nPU[0];
    
    //Init GenLepton Variables
    InitGenLeptonVariables();
    //Init Lepton Variables
    InitLeptonVariables();
    
    //Set gen level lepton indix
    ResetGenLeptonIndex();//Resets Lepton indixes
    SetGenElectronIndex();//Set electron indixes
    SetGenMuonIndex();//Set muon indixes
    SetGenTauIndex();//Set tau indexes
    SetGenLeptonVector();//Set TLorentz vector for the two leading leptons
    
    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if( _isData ) 
      {
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
    if ( _isData && _goodRunLS && !isGoodRunLS() ) 
      {
	if ( lastRun != runNumber || lastLumi != lumiBlock) 
	  {
	    lastRun = runNumber;
	    lastLumi = lumiBlock;
	  }
	continue;
      }
    
    if ( _isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) 
      {
	lastRun = runNumber;
	lastLumi = lumiBlock;
      }
    
    Npassed_In += weightII;
    
    //HLT and Data Filter
    //passedHLT = HLT_Razor + HLT_Razor_prescaled;
    passedHLT = HLT_Razor;
    
    if ( _isData == true ) 
      {
	if ( passedHLT == 0 ) continue;//Comment out for getting trigger turn-ons
	if ( (ECALTPFilterFlag==0) || (drBoundary==0) || (drDead==0) || (CSCHaloFilterFlag==0) 
	     || (trackerFailureFilterFlag==0) || (BEECALFlag==0) || ( HBHENoiseFilterResultFlag ==0 )
	     || (ecalLaserFilter == 0) || (eeBadScFilterFlag == 0) || (hcalLaserFilter == 0)) continue;
      }
    
    // find highest-pT PV [replace with Sagar's code]
    int iPV = passPV();
    if( iPV < 0 ) continue;//If negative no PV found to pass the cuts
    Npassed_PV += weightII;
    nPV = N_PV_EVENT;
    
    /////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    //////////////////////// PF JETS + JetID ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    
    vector<TLorentzVector> pfJets;
    vector<int> i_pfJets;
    int N_pfJets = DoPfSelection(pfJets, i_pfJets);
    // jet ID                                                                 
    if (N_pfJets <= 0 )  continue;// If any Jet is bad (see loop before) event is rejected
    
    //////////////////////////////////////////////////////////////
    /////////////////////Selecting Muons//////////////////////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    std::vector< VecbosLepton > LooseLepton;
    std::vector< VecbosLepton > TightLepton;
    for( int i = 0; i < nMuon; i++ ) {
      VecbosLepton tmp;
      TLorentzVector thisMu( pxMuon[i], pyMuon[i], pzMuon[i], energyMuon[i] );
      if ( thisMu.Pt() < 15.0 || fabs( thisMu.Eta() ) > 2.4 ) continue;
      tmp.index = i;
      tmp.lepton = thisMu;
      tmp.charge = chargeMuon[i];
      tmp.mass = muon_mass;
      tmp.pdgID = muon_pdgID;
      tmp._isLoose = isLooseMuon(i, true);
      tmp._isTight = isTightMuon(i, true);
      if ( tmp._isLoose ) 
	{                                            
	  LooseLepton.push_back(tmp);
	}
      
      if ( tmp._isTight )
	{
	  TightLepton.push_back(tmp);
	}
    }//end muon loop

    //////////////////////////////////////////////////////////////
    /////////////////////Selecting Electrons///////////////////////
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    for( int i = 0; i < nEle; i++ ) {
      VecbosLepton tmp;
      TLorentzVector thisEle( pxEle[i], pyEle[i], pzEle[i], energyEle[i] );
      if ( thisEle.Pt() < 15.0 || fabs( thisEle.Eta() ) > 2.5 ) continue;
      //Look for matching muon already store in the collection                                                                      
      bool matchMuon = false;
      for ( auto lep : LooseLepton ) {
	if ( lep.lepton.DeltaR( thisEle ) < 0.1 )
	  {
	    matchMuon = true;
	    break;
	  }
      }//end loose lepton loop
      if( matchMuon ) break;//Discard electron if matches a muon
      
      tmp.index = i;
      tmp.lepton = thisEle;
      tmp.charge = chargeEle[i];
      tmp.mass = ele_mass;
      tmp.pdgID= ele_pdgID;
      tmp._isLoose = isLooseElectron(i);
      tmp._isTight = isTightElectron(i);
      if ( tmp._isLoose )
        {
	  LooseLepton.push_back( tmp );
	}
      
      if ( tmp._isTight )
	{
          TightLepton.push_back( tmp );
        }
    }//end electron loop
    
    //Veto Event if there are no muons or electrons
    if ( LooseLepton.size() == 0 ) continue;
    //Filling Letptons. PT order taken into account automatically
    FillLeptons( LooseLepton );
    SortByPt( LooseLepton );
    
    /*
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
    */
    
    ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    ///////////////// Create Collection  /////////////////    
    ////////////////  pfJets muon subtracted ////////////
    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    /*
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
    */
    /*
    //JEC Up and Down Systematic
    vector<TLorentzVector> pfJets_noMu_Up;
    vector<TLorentzVector> pfJets_noMu_Down;
    for(int j = 0; j < pfJets_noMu.size(); j++){
      //JEC Up
      jec_un->setJetEta(pfJets_noMu[j].Eta());
      jec_un->setJetPt(pfJets_noMu[j].Pt());
      //Getting JEC;
      double deltaPt = fabs(jec_un->getUncertainty(true));//Fractional Error
      //JEC Up Now
      double scl = 1.0 + deltaPt;
      TLorentzVector aux_J;
      aux_J.SetPxPyPzE(scl*pfJets_noMu[j].Px(), scl*pfJets_noMu[j].Py(), scl*pfJets_noMu[j].Pz(), scl*pfJets_noMu[j].E());
      pfJets_noMu_Up.push_back(aux_J);
      //JEC Down
      scl = 1.0 - deltaPt;
      aux_J.SetPxPyPzE(scl*pfJets_noMu[j].Px(), scl*pfJets_noMu[j].Py(), scl*pfJets_noMu[j].Pz(), scl*pfJets_noMu[j].E());
      pfJets_noMu_Down.push_back(aux_J);
    }
    */

    /*
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
    //if (iTauTight.size()>0) continue;//removed after a bug was found in the tau ID

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
    if(IsoPF)continue;//Applying ILV
    
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
    */
    
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
    //outTree->Fill();
    events->tree_->Fill();
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
  
  events->tree_->Write();
  effTree->Write();
  file->Close();
}

int RazorRunTwo::HighestPt(vector<TLorentzVector> p, int iHIGHEST) {
  
  int iH = -99;
  double highestPT = 0.;
  for(int i=0; i<p.size();i++)
    {
      if((p[i].Pt()>= highestPT) && (i != iHIGHEST))
	{
	  iH = i;
	  highestPT = p[i].Pt();
	}
    }
  return iH;
};

bool RazorRunTwo::SetGenElectronIndex()
{
  bool _findLepton = false;
  for(int i = 0; i < nMc; i++)
    {
      float pt = pMc[i]/cosh(etaMc[i]);
      if( abs(idMc[i]) == 11 && pt > 5.0 && fabs(etaMc[i]) < 2.5 && statusMc[i] == 1 )
	{
	  if( ( abs(idMc[mothMc[i]]) == 23 || abs(idMc[mothMc[i]]) == 24 ) || 
	      ( ( abs(idMc[mothMc[i]]) == 11 || abs(idMc[mothMc[i]]) == 13 || abs(idMc[mothMc[i]]) == 15 ) 
		&& statusMc[mothMc[i]] == 3 && mothMc[mothMc[i]] >= 0 && 
		( abs(idMc[mothMc[mothMc[i]]]) == 23 || abs(idMc[mothMc[mothMc[i]]]) == 24 ) 
		) )
	    {
	      genLeptonIndex.push_back(i);
	      _findLepton =  true;
	    }
	}//big if ends
  }//loop over gen particles
  
  return _findLepton;
  
};

bool RazorRunTwo::SetGenMuonIndex()
{
  bool _findLepton = false;
  for(int i = 0; i < nMc; i++)
    {
      float pt = pMc[i]/cosh(etaMc[i]);
      if( abs(idMc[i]) == 13 && pt > 5.0 && fabs(etaMc[i]) < 2.4 && statusMc[i] == 1 )
	{
	  if( ( abs(idMc[mothMc[i]]) == 23 || abs(idMc[mothMc[i]]) == 24 ) || 
	      ( ( abs(idMc[mothMc[i]]) == 11 || abs(idMc[mothMc[i]]) == 13 || abs(idMc[mothMc[i]]) == 15 ) 
		&& statusMc[mothMc[i]] == 3 && mothMc[mothMc[i]] >= 0 && 
		( abs(idMc[mothMc[mothMc[i]]]) == 23 || abs(idMc[mothMc[mothMc[i]]]) == 24 ) 
		) )
	    {
	      genLeptonIndex.push_back(i);
	      _findLepton = true;
	    }
	}//big if ends
    }//loop over gen particles
  
  return _findLepton;
  
};

bool RazorRunTwo::SetGenTauIndex()
{
  bool _findLepton = false;
  for(int i = 0; i < nMc; i++)
    {
      float pt = pMc[i]/cosh(etaMc[i]);
      if( abs(idMc[i]) == 15 && pt > 2.0 && fabs(etaMc[i]) < 2.4 && statusMc[i] == 2 )
	{
	  if( ( abs(idMc[mothMc[i]]) == 23 || abs(idMc[mothMc[i]]) == 24 ) || 
	      ( ( abs(idMc[mothMc[i]]) == 11 || abs(idMc[mothMc[i]]) == 13 || abs(idMc[mothMc[i]]) == 15 ) 
		&& statusMc[mothMc[i]] == 3 && mothMc[mothMc[i]] >= 0 && 
		( abs(idMc[mothMc[mothMc[i]]]) == 23 || abs(idMc[mothMc[mothMc[i]]]) == 24 ) 
		) )
	    {
	      genLeptonIndex.push_back(i);
	      _findLepton = true;
	    } 
	}//big if end
    }//loop over gen particles
  
  return _findLepton;
  
};

void RazorRunTwo::SetGenLeptonVector()
{
  int i_ledLep  = -99;//Indices for leading lepton
  int i_subLedLep = -99;//Indices for subleading lepton 
  float max_pt = -99.0;//store max pt in each iteration when appropiate
  //finds leading PT lepton
  for( int i : genLeptonIndex )
    {
      float pt = pMc[i]/cosh(etaMc[i]);
      if( pt > max_pt )
	{
	  max_pt = pt;
	  i_ledLep = i;
	}
    }
  //find subleading PT lepton 
  max_pt = -99.0;//Resetting max_pt
  for( int i : genLeptonIndex )
    {
      float pt = pMc[i]/cosh(etaMc[i]);
      if( pt > max_pt && i != i_ledLep )
        {
          max_pt = pt;
          i_subLedLep = i;
        }
    }
  
  //Fill leading lepton information
  if( i_ledLep  >= 0 )
    {
      float pt = pMc[i_ledLep]/cosh(etaMc[i_ledLep]);
      double mass = 0.000511;
      if( abs(idMc[i_ledLep]) == 13 )mass = 0.1057;
      if( abs(idMc[i_ledLep]) == 15 )mass = 1.777;
      events->genlep1.SetPtEtaPhiM( pt, etaMc[i_ledLep], phiMc[i_ledLep], mass );
      events->genlep1Type = idMc[i_ledLep];
    }
  //Fill subleading lepton information
  if( i_subLedLep  >= 0 )
    {
      float pt = pMc[i_subLedLep]/cosh(etaMc[i_subLedLep]);
      double mass = 0.000511;
      if( abs(idMc[i_subLedLep]) == 13 )mass = 0.1057;
      if( abs(idMc[i_subLedLep]) == 15 )mass = 1.777;
      events->genlep2.SetPtEtaPhiM( pt, etaMc[i_subLedLep], phiMc[i_subLedLep], mass );
      events->genlep2Type = idMc[i_subLedLep];
    }
  
};

void RazorRunTwo::ResetGenLeptonIndex()
{
  genLeptonIndex.clear();
};

int RazorRunTwo::DoPfSelection(std::vector<TLorentzVector>& pfJets, std::vector<int>& i_pfJets)
{
  vector<double> pfJets_f_photon, pfJets_f_electron, pfJets_f_muon, pfJets_f_neutralhad, pfJets_f_chargedhad, \
    pfJets_f_hfhad, pfJets_f_hfem;
  vector<double> pfJets_mult_photon, pfJets_mult_electron, pfJets_mult_muon, pfJets_mult_neutralhad, \
    pfJets_mult_chargedhad, pfJets_mult_hfhad, pfJets_mult_hfem;
  bool bad_pfJet = false;
  int N_pfJets = 0, pfBtags = 0;
  double pfHT, pfMHTx, pfMHTy;
  
  struct pfJetStruct{
    TLorentzVector Jet;
    int index;
  };
  
  std::map< double, pfJetStruct > map_pt;
  std::vector< double > PtVec;
  
  bool good_pfjet = false;
  for(int i = 0; i < nAK5PFNoPUJet; i++){                                                                                     
    TLorentzVector jet;
    pfJetStruct aux_pfJetStruct;
    
    double px = pxAK5PFNoPUJet[i];                                                                                            
    double py = pyAK5PFNoPUJet[i];                                                                                            
    double pz = pzAK5PFNoPUJet[i];                                                                                            
    double E = sqrt(px*px+py*py+pz*pz);                                                                                       
    double scale = 1.;                                                                                                        
    jet.SetPxPyPzE(scale*px,scale*py,scale*pz,scale*E);
    good_pfjet = false;                                                                                                       
    double EU = uncorrEnergyAK5PFNoPUJet[i];
    if( jet.Pt() > 40.0 && fabs(jet.Eta()) < 3.0 )
      {
	double fHAD = (neutralHadronEnergyAK5PFNoPUJet[i]+chargedHadronEnergyAK5PFNoPUJet[i])/EU;                              
	if(fHAD > 0.99)
	  {                                                                                                       
	    N_pfJets = 0;                                                                                                      
	    break;// clean NOISY event
	  }
	
	int nConstituents = chargedHadronMultiplicityAK5PFNoPUJet[i]+neutralHadronMultiplicityAK5PFNoPUJet[i]+
	  photonMultiplicityAK5PFNoPUJet[i]+electronMultiplicityAK5PFNoPUJet[i]+muonMultiplicityAK5PFNoPUJet[i]+
	  HFHadronMultiplicityAK5PFNoPUJet[i]+HFEMMultiplicityAK5PFNoPUJet[i];
	int chargedMult = chargedHadronMultiplicityAK5PFNoPUJet[i]+
	  electronMultiplicityAK5PFNoPUJet[i]+
	  muonMultiplicityAK5PFNoPUJet[i];
    
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
    
	if( (neutralHadFrac < 0.99) && (photonFrac < 0.99) && (nConstituents > 1) ) 
	  {
	    //outside of tracker acceptance, these are the only requirements
	    if (fabs(jet.Eta())>=2.4) good_pfjet = true;
	    //inside of the tracker acceptance, there are extra requirements              
	    else {
	      if ((chargedHadFrac > 0.0) && (chargedMult > 0) && (electronFrac < 0.99)) good_pfjet = true;
	    }
	  }
	
	if(good_pfjet)
	  {
	    N_pfJets++;
	    //pfJets.push_back(jet);
	    //i_pfJets.push_back(i);
	    aux_pfJetStruct.Jet = jet;
	    aux_pfJetStruct.index = i;
	    if( map_pt.find(jet.Pt()) == map_pt.end() )
	      {
		map_pt[jet.Pt()] = aux_pfJetStruct;
		PtVec.push_back(jet.Pt());
	      }
	    else{
	      std::cerr << "Identical Jet PT!" << std::endl;
	    }
	    pfHT += jet.Pt();
	    pfMHTx -= jet.Px();
	    pfMHTy -= jet.Py();
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
	  }
	else 
	  {
	    N_pfJets = 0;
	    break;//Only takes out the pfJets loop! But good_jet = false
	  }
      }
  }
  
  std::sort( PtVec.begin(), PtVec.end() );
  std::reverse( PtVec.begin(), PtVec.end() );
  for( double tmp : PtVec){
    pfJets.push_back( map_pt[tmp].Jet );
    i_pfJets.push_back( map_pt[tmp].index );
  }
  return N_pfJets;
  
};

void RazorRunTwo::FillJetInfo(vector<TLorentzVector> GoodJets, vector<int> GoodJetIndices, vector<TLorentzVector> GoodLeptons){
    //NOTE: GoodJets should be sorted by jet pT!

    //reset event variables
    events->HT = 0;
    events->NJets40 = 0;
    events->NBJetsLoose = 0;
    events->NBJetsMedium = 0;
    events->NBJetsTight = 0;
    events->jet1 = TLorentzVector();
    events->jet2 = TLorentzVector();
    events->bjet1 = TLorentzVector();
    events->bjet2 = TLorentzVector();
    events->jet1PassCSVLoose = false;
    events->jet1PassCSVMedium = false;
    events->jet1PassCSVTight = false;
    events->jet2PassCSVLoose = false;
    events->jet2PassCSVMedium = false;
    events->jet2PassCSVTight = false;
    events->bjet1PassLoose = false;
    events->bjet1PassMedium = false;
    events->bjet1PassTight = false;
    events->bjet2PassLoose = false;
    events->bjet2PassMedium = false;
    events->bjet2PassTight = false;
    events->MR = -1; events->MR_NoDilepton = -1; 
    events->Rsq = -1; events->Rsq_NoDilepton = -1;
    events->MET = -1; events->MET_NoDilepton = -1;
    events->minDPhi = -1; events->minDPhiN = -1;

    if(GoodJets.size() == 0) return;

    //get the PF MET as a TLorentzVector
    TLorentzVector PFMET(pxPFMet[2], pyPFMet[2], 0, sqrt(pxPFMet[2]*pxPFMet[2] + pyPFMet[2]*pyPFMet[2]));
    TVector3 PFMET3(pxPFMet[2], pyPFMet[2], 0);
    events->MET = PFMET3.Pt();
    TLorentzVector PFMETWithLeptons = PFMET; //will add leptons to this

    bool gotLeadJet = false; bool gotSubLeadJet = false;
    bool gotLeadBJet = false; bool gotSubLeadBJet = false;
    vector<TLorentzVector> GoodJetsWithoutLeptons;
    for(int j = 0; j < GoodJets.size(); j++){
        //check if this jet matches a good lepton
        //(if it does, skip it)
        double dR = -1;
        for(auto& lep : GoodLeptons){
            double thisDR = GoodJets[j].DeltaR(lep);
            if(dR < 0 || thisDR < dR) dR = thisDR;
        }
        if(dR > 0 && dR < 0.5){ //a selected lepton is inside this jet -- add to PFMETWithLeptons
            PFMETWithLeptons = PFMETWithLeptons + GoodJets[j];
            continue; 
        }

        GoodJetsWithoutLeptons.push_back(GoodJets[j]);

        if(GoodJets[j].Pt() > 40) events->NJets40++;
        events->HT += GoodJets[j].Pt();
        if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->NBJetsLoose++;
        if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->NBJetsMedium++;
        if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->NBJetsTight++;
        //minDPhiN: min phi angle between jet and MET
        double thisDPhi = GoodJets[j].DeltaPhi(PFMET);
        if(events->minDPhi < 0 || thisDPhi < events->minDPhi) events->minDPhi = thisDPhi;

        //fill info on first two jets
        if(!gotLeadJet){
            events->jet1 = GoodJets[j];
            if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->jet1PassCSVLoose = true;
            else events->jet1PassCSVLoose = false;
            if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->jet1PassCSVMedium = true;
            else events->jet1PassCSVMedium = false;
            if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->jet1PassCSVTight = true;
            else events->jet1PassCSVTight = false;
            gotLeadJet = true;
        }
        else if(!gotSubLeadJet){
            events->jet2 = GoodJets[j];
            if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->jet2PassCSVLoose = true;
            else events->jet2PassCSVLoose = false;
            if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->jet2PassCSVMedium = true;
            else events->jet2PassCSVMedium = false;
            if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->jet2PassCSVTight = true;
            else events->jet2PassCSVTight = false;
            gotSubLeadJet = true;
        }

        //fill info on first two b-jets (CSVM)
        if(!gotLeadBJet && pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])){
            events->bjet1 = GoodJets[j];
            if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->bjet1PassLoose = true;
            else events->bjet1PassLoose = false;
            if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->bjet1PassMedium = true;
            else events->bjet1PassMedium = false;
            if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->bjet1PassTight = true;
            else events->bjet1PassTight = false;
            gotLeadBJet = true;
        }
        else if(!gotSubLeadBJet && pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])){
            events->bjet2 = GoodJets[j];
            if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->bjet2PassLoose = true;
            else events->bjet2PassLoose = false;
            if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->bjet2PassMedium = true;
            else events->bjet2PassMedium = false;
            if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->bjet2PassTight = true;
            else events->bjet2PassTight = false;
            gotSubLeadBJet = true;
        }
    }

    TVector3 PFMET3WithLeptons(PFMETWithLeptons.Px(), PFMETWithLeptons.Py(), 0);
    events->MET_NoDilepton = PFMET3WithLeptons.Pt();

    //compute MR and Rsq
    vector<TLorentzVector> hemispheres = CombineJets(GoodJets);
    events->MR = CalcGammaMRstar(hemispheres[0], hemispheres[1]);
    double MTR = CalcMTR(hemispheres[0], hemispheres[1], PFMET3);
    double R = -999;
    if(events->MR > 0) R = MTR/events->MR;
    events->Rsq = R*R;

    //MR and Rsq with leptons added to MET
    vector<TLorentzVector> hemispheresNoLeps = CombineJets(GoodJetsWithoutLeptons);
    events->MR_NoDilepton = CalcGammaMRstar(hemispheresNoLeps[0], hemispheresNoLeps[1]);
    double MTRNoLeps = CalcMTR(hemispheresNoLeps[0], hemispheresNoLeps[1], PFMET3WithLeptons);
    double RNoLeps = -999;
    if(events->MR_NoDilepton > 0) RNoLeps = MTRNoLeps/events->MR_NoDilepton;
    events->Rsq_NoDilepton = RNoLeps*RNoLeps;

    //compute minDeltaPhiN (see Appendix D in CMS note AN-2011-409)
    for(int i = 0; i < GoodJetsWithoutLeptons.size(); i++){
        //compute deltaPhi between jet and MET
        double thisDPhi = GoodJetsWithoutLeptons[i].DeltaPhi(PFMET);
        //compute estimated 'perpendicular missing energy' change due to possible jet mismeasurements
        //first compute sum((px_i*py_j - py_i*px_j)^2), summing over all other jets
        double sum = 0;
        for(int j = 0; j < GoodJetsWithoutLeptons.size(); j++){
            if(i == j) continue;
            sum = sum + pow(GoodJetsWithoutLeptons[i].Px()*GoodJetsWithoutLeptons[j].Py() - GoodJetsWithoutLeptons[i].Py()*GoodJetsWithoutLeptons[j].Px(), 2);
        }
        double thisDeltaT = 0.1*sqrt(sum)/GoodJetsWithoutLeptons[i].Pt();
        //compute the normalized delta phi
        double thisDeltaPhiN = thisDPhi/atan(thisDeltaT/PFMET.Pt());
        if(events->minDPhiN < 0 || thisDeltaPhiN < events->minDPhiN) events->minDPhiN = thisDeltaPhiN;
    }
}

void RazorRunTwo::SortByPt(std::vector<VecbosLepton>& lepton)
{
  std::vector< VecbosLepton > aux;
  //Pt key map for lepton
  std::map< double, VecbosLepton > lepton_map;
  std::vector< double > lepton_pt;
  for ( auto tmp : lepton ) {
    if ( lepton_map.find( tmp.lepton.Pt() ) == lepton_map.end() )
      {
	lepton_map[ tmp.lepton.Pt() ] = tmp;
	lepton_pt.push_back( tmp.lepton.Pt() );
      }
    else
      {
	std::cerr << "[ERROR]: Identical lepton PT" << std::endl;
      }
  }//end lepton loop
  
  //sorting by pt 
  std::sort( lepton_pt.begin(), lepton_pt.end() );
  std::reverse( lepton_pt.begin(), lepton_pt.end() );
  
  //Cleaning lepton vector
  lepton.clear();
  //Filling PT ordered lepton vector
  for ( double tmp : lepton_pt ){
    if ( lepton_map.find( tmp ) != lepton_map.end() )
      {
	lepton.push_back( lepton_map[tmp] );
      }
    else
      {
	std::cerr << "[ERROR]: lepton PT not found" << std::endl;
      }
  }
  
};

void RazorRunTwo::FillLeptons(std::vector<VecbosLepton> lepton)
{
  SortByPt( lepton );//Sorting Leptons by PT
  int n_lepton = 0;
  //Fill Two Leading PT leptons
  for ( auto tmp : lepton ){
    if( n_lepton >= 2 ) break;
    if( n_lepton == 0 )
      {
	events->lep1.SetPtEtaPhiM( tmp.lepton.Pt(), tmp.lepton.Eta(), tmp.lepton.Phi(), tmp.mass );
	events->lep1PassVeto = tmp._isLoose;
	events->lep1PassLoose = tmp._isLoose;
	events->lep1PassTight = tmp._isTight;
	events->lep1Type = -1.0*tmp.pdgID*tmp.charge;
	events->lep1MatchedGenLepIndex = MatchLeptonGenLevel( tmp.lepton );
      }
    else if ( n_lepton == 1 ) 
      {
	events->lep2.SetPtEtaPhiM( tmp.lepton.Pt(), tmp.lepton.Eta(), tmp.lepton.Phi(), tmp.mass );
	events->lep2PassVeto = tmp._isLoose;
        events->lep2PassLoose = tmp._isLoose;
        events->lep2PassTight = tmp._isTight;
        events->lep2Type = -1.0*tmp.pdgID*tmp.charge;
	events->lep2MatchedGenLepIndex = MatchLeptonGenLevel( tmp.lepton );
      }
    n_lepton++;//increase lepton counter
  }//end lepton loop
};

int RazorRunTwo::MatchLeptonGenLevel(TLorentzVector lepton)
{
  double deltaR_gen1 = 9999.0;
  double deltaR_gen2 = 9999.0;
  //first gen lepton
  if ( events->genlep1.M() != 0.0 )
    {
      if ( lepton.DeltaR( events->genlep1 ) < genLepton_DR )
	{
	  deltaR_gen1 = lepton.DeltaR( events->genlep1 );
	}
    }
  //second gen lepton
  if ( events->genlep2.M() != 0.0 )
    {
      if ( lepton.DeltaR( events->genlep2 ) < genLepton_DR )
        {
          deltaR_gen2 = lepton.DeltaR( events->genlep2 );
        }
    }
  
  //return index to closes gen lepton
  
  //no match
  if( deltaR_gen1 == 9999.0 && deltaR_gen2 == 9999.0 ) return -1;
  
  if( deltaR_gen1 < deltaR_gen2 )
    {
      return 1;
    }
  else if ( deltaR_gen1 > deltaR_gen2 )
    {
      return 2;
    }
  else if ( deltaR_gen1 == deltaR_gen2 && deltaR_gen1 != 9999.0 )
    {
      return 3;
    } 

  //if everything fails return no match
  return -1;  
};

void RazorRunTwo::InitGenLeptonVariables()
{
  events->genlep1.SetPtEtaPhiM(0,0,0,0);
  events->genlep2.SetPtEtaPhiM(0,0,0,0);;
  events->genlep1Type = 0;
  events->genlep2Type = 0;
};
void RazorRunTwo::InitLeptonVariables()
{
  events->lep1.SetPtEtaPhiM(0,0,0,0);
  events->lep2.SetPtEtaPhiM(0,0,0,0);
  events->lep1Type = 0;
  events->lep2Type = 0;
  events->lep1MatchedGenLepIndex = -1;
  events->lep2MatchedGenLepIndex = -1;
  events->lep1PassVeto = false;
  events->lep1PassLoose = false;
  events->lep1PassTight = false;
  events->lep2PassVeto = false;
  events->lep2PassLoose = false;
  events->lep2PassTight = false;
};

float RazorRunTwo::GetMTLep()
{
  return sqrt(events->lep1.M2() + 2*metPt*events->lep1.Pt()*(1 - cos(deltaPhi(metPhi,events->lep1.Phi()))));
};
