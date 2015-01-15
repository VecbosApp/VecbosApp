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
  double Npassed_LepVeto = 0;
  //B-tag
  double Npassed_0btag = 0;

  double weightII = 1.;
  unsigned int lastLumi = 0;
  unsigned int lastRun = 0;
  
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
    FillJetInfo(pfJets, i_pfJets, LooseLepton);
    
    //Filling MTlepton
    FillMTLep();
    events->tree_->Fill();
  }
  

  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
  effTree->Branch("Npassed_In", &Npassed_In, "Npassed_In/D");
  effTree->Branch("Npassed_ISR", &w_isr, "Npassed_ISR/D");
  effTree->Branch("Npassed_ISR_up", &w_isr_up, "Npassed_ISR_up/D");
  effTree->Branch("Npassed_ISR_down", &w_isr_down, "Npassed_ISR_down/D");

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

void RazorRunTwo::FillJetInfo(vector<TLorentzVector> GoodJets, vector<int> GoodJetIndices, vector<VecbosLepton> GoodLeptons){
    //NOTE: GoodJets should be sorted by jet pT!

    //reset event variables
    events->HT = 0;
    events->NJets40 = 0;
    events->NJets80 = 0;
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
    events->MR = -1; events->MR_LeptonsAsMET = -1; events->MR_LeadLeptonsAsMET = -1;
    events->Rsq = -1; events->Rsq_LeptonsAsMET = -1; events->Rsq_LeadLeptonsAsMET = -1;
    events->MET = -1; events->MET_LeptonsAsMET = -1; events->MET_LeadLeptonsAsMET = -1;
    events->minDPhi = -1; events->minDPhiN = -999; 
    events->minDPhiFolded = -1; events->minDPhiNFolded = -999;
    events->dPhiHemHem = -99;
    events->dPhiHemHem_LeptonsAsMET = -99;
    events->dPhiHemHem_LeadLeptonsAsMET = -99;

    if(GoodJets.size() == 0) return;

    //get the PF MET as a TVector3
    TVector3 PFMET(pxPFMet[2], pyPFMet[2], 0);
    events->MET = PFMET.Pt();
    TVector3 PFMETWithLeptons = PFMET; //add leptons to this
    TVector3 PFMETWithLeadingLeptons = PFMET; //add the two highest pt leptons to this
    int nLepsAddedToMET = 0;
    for(auto& lep : GoodLeptons){
        //add to PFMETWithLeptons
        PFMETWithLeptons = PFMETWithLeptons + lep.lepton.Vect();
        PFMETWithLeptons.SetZ(0.0);
        //add first two leptons to PFMETWithLeadingLeptons
        if(nLepsAddedToMET < 2){
            PFMETWithLeadingLeptons = PFMETWithLeadingLeptons + lep.lepton.Vect();
            PFMETWithLeadingLeptons.SetZ(0.0);
            nLepsAddedToMET++;
        }
    }
    events->MET_LeptonsAsMET = PFMETWithLeptons.Pt();
    events->MET_LeadLeptonsAsMET = PFMETWithLeadingLeptons.Pt();

    bool gotLeadJet = false; bool gotSubLeadJet = false;
    bool gotLeadBJet = false; bool gotSubLeadBJet = false;
    vector<TLorentzVector> GoodJetsWithoutLeptons;
    vector<TLorentzVector> GoodJetsWithoutLeadingLeptons;
    for(int j = 0; j < GoodJets.size(); j++){
        //check if this jet matches a good lepton
        //(if it does, skip it)
        double dR = -1;
        int nLepsFoundInJetCollection = 0;
        for(auto& lep : GoodLeptons){
            double thisDR = GoodJets[j].DeltaR(lep.lepton);
            if(dR < 0 || thisDR < dR) dR = thisDR;
        }
        if(dR > 0 && dR < 0.5){ //a selected lepton is inside this jet
            //if we already have identified two leptons in the jet collection, add the jet to GoodJetsWithoutLeadingLeptons
            if(nLepsFoundInJetCollection >= 2) GoodJetsWithoutLeadingLeptons.push_back(GoodJets[j]);
            nLepsFoundInJetCollection++;
            continue; 
        }

        GoodJetsWithoutLeptons.push_back(GoodJets[j]);
        GoodJetsWithoutLeadingLeptons.push_back(GoodJets[j]);

        if(GoodJets[j].Pt() > 40) events->NJets40++;
        if(GoodJets[j].Pt() > 80) events->NJets80++;
        events->HT += GoodJets[j].Pt();
        if(pfJetPassCSVL(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->NBJetsLoose++;
        if(pfJetPassCSVM(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->NBJetsMedium++;
        if(pfJetPassCSVT(combinedSecondaryVertexBJetTagsAK5PFNoPUJet[GoodJetIndices[j]])) events->NBJetsTight++;
        //minDPhiN: min phi angle between jet and MET
        double thisDPhi = fabs(GoodJets[j].Vect().DeltaPhi(PFMET));
        double thisFoldedDPhi = min(thisDPhi, fabs(3.14159 - thisDPhi)); // we want the jet closest to being aligned OR anti-aligned with the MET
        if(events->minDPhiFolded < 0 || thisFoldedDPhi < events->minDPhiFolded) events->minDPhiFolded = thisFoldedDPhi;
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

    //compute MR and Rsq
    vector<TLorentzVector> hemispheres = CombineJets(GoodJets);
    events->MR = CalcGammaMRstar(hemispheres[0], hemispheres[1]);
    double MTR = CalcMTR(hemispheres[0], hemispheres[1], PFMET);
    double R = -999;
    if(events->MR > 0) R = MTR/events->MR;
    events->Rsq = R*R;
    //compute the transverse angle between the two hemispheres
    events->dPhiHemHem = hemispheres[0].DeltaPhi(hemispheres[1]);

    //MR and Rsq with leptons added to MET
    vector<TLorentzVector> hemispheresNoLeps = CombineJets(GoodJetsWithoutLeptons);
    events->MR_LeptonsAsMET = CalcGammaMRstar(hemispheresNoLeps[0], hemispheresNoLeps[1]);
    double MTRNoLeps = CalcMTR(hemispheresNoLeps[0], hemispheresNoLeps[1], PFMETWithLeptons);
    double RNoLeps = -999;
    if(events->MR_LeptonsAsMET > 0) RNoLeps = MTRNoLeps/events->MR_LeptonsAsMET;
    events->Rsq_LeptonsAsMET = RNoLeps*RNoLeps;
    //compute the transverse angle between the two hemispheres
    events->dPhiHemHem_LeptonsAsMET = hemispheresNoLeps[0].DeltaPhi(hemispheresNoLeps[1]);

    //MR and Rsq with leading leptons added to MET
    vector<TLorentzVector> hemispheresNoLeadLeps = CombineJets(GoodJetsWithoutLeadingLeptons);
    events->MR_LeadLeptonsAsMET = CalcGammaMRstar(hemispheresNoLeadLeps[0], hemispheresNoLeadLeps[1]);
    double MTRNoLeadLeps = CalcMTR(hemispheresNoLeadLeps[0], hemispheresNoLeadLeps[1], PFMETWithLeadingLeptons);
    double RNoLeadLeps = -999;
    if(events->MR_LeadLeptonsAsMET > 0) RNoLeadLeps = MTRNoLeadLeps/events->MR_LeadLeptonsAsMET;
    events->Rsq_LeadLeptonsAsMET = RNoLeadLeps*RNoLeadLeps;
    //compute the transverse angle between the two hemispheres
    events->dPhiHemHem_LeadLeptonsAsMET = hemispheresNoLeadLeps[0].DeltaPhi(hemispheresNoLeadLeps[1]);

    //compute minDeltaPhiN (see Appendix D in CMS note AN-2011-409)
    for(int i = 0; i < GoodJetsWithoutLeptons.size(); i++){
        //compute deltaPhi between jet and MET
        double thisDPhi = GoodJetsWithoutLeptons[i].Vect().DeltaPhi(PFMET);
        double thisFoldedDPhi = min(thisDPhi, fabs(3.14159 - thisDPhi));
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
        double thisDeltaPhiNFolded = thisFoldedDPhi/atan(thisDeltaT/PFMET.Pt());
        if(events->minDPhiNFolded < 0 || thisDeltaPhiNFolded < events->minDPhiNFolded) events->minDPhiNFolded = thisDeltaPhiNFolded;
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
  TVector3 lepton( events->lep1.Px(), events->lep1.Py(), events->lep1.Pz() );
  TVector3 met( pxPFMet[2], pyPFMet[2], 0.0 );
  double deltaPhi = lepton.DeltaPhi( met );
  double LeptonM2 = events->lep1.M2();
  return sqrt( LeptonM2 + 2*met.Pt()*lepton.Pt()*( 1.0 - cos( deltaPhi ) ) );
};

void RazorRunTwo::FillMTLep()
{
  events->lep1MT = GetMTLep();
};
