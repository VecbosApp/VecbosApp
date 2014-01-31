// std includes
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

// ROOT includes
#include <TTree.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>

// local includes
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CoolTools.hh"
#include "CaloTower.hh"
#include "AnalysisSelector.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/ElectronEffectiveArea.h"
#include "CommonTools/include/MuonEffectiveArea.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "Timing.hh"

Timing::Timing(TTree *tree) : Vecbos(tree) {  
  
}

Timing::~Timing(){
}

void Timing::Loop(string outFileName, int ISDATA) {
  if(fChain == 0) return;
  
  s_output = outFileName;
  is_DATA = ISDATA;

  //Initiate event dump tree
  InitEventTree();

  int NTOT = 0; //total number of events written to output
  
  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "here" << endl;
  Long64_t nentries = fChain->GetEntries();
  cout << "and here" << endl;
  cout << "Number of entries = " << nentries << endl;
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    fail2011Filter = 0;
    failHPDHits = 0;
    failHPDNoOtherHits = 0;
    failMaxZeros = 0;

    nb = fChain->GetEntry(jentry);   nbytes += nb;
      
    if (jentry%10000 == 0)
      cout << ">>> Processing event # " << jentry << endl;
    
    if(is_DATA) {
      //Check LS/RUN 
      if(isGoodRunLS() == false){
	continue;
      }
      if((fail2011Filter == 1) || (failHPDHits == 1) || 
	 (failHPDNoOtherHits == 1) || (failMaxZeros == 1)){
	continue;
      }
    }
      
    if(PassLeptonSelection() == false) continue;
    if(PassPVSelection() == false) continue;
    if(PassSCSelection() == false) continue;

    InitTracks();
    if(PassJetSelection() == false) continue;

    AddEvent();
    NTOT++;
  } 
  
  if(NTOT > 0){
    cout << "Events selected - Writing tree to " << s_output.c_str() << endl;
      f_out = new TFile(s_output.c_str(),"RECREATE");
      WriteEventTree();
      
      f_out->Close();
  } else {
    cout << "No events selected - not writing tree to " << s_output.c_str() << endl;
    f_out = new TFile(s_output.c_str(),"RECREATE");
    WriteEventTree();
    
    f_out->Close();
  }
}

void Timing::InitEventTree(){

  myEventT = (TTree*) new TTree("Event","Event");

  //first common event stuff
  myEventT->Branch("MET_FLAGS", &MET_FLAGS, "MET_FLAGS/I");
  myEventT->Branch("TRACKER_FLAGS", &TRACKER_FLAGS, "TRACKER_FLAGS/I");
  myEventT->Branch("RUN_NUM", &RUN_NUM, "RUN_NUM/I");
  myEventT->Branch("LS_NUM", &LS_NUM, "LS_NUM/I");
  myEventT->Branch("EVENT_NUM", &EVENT_NUM, "EVENT_NUM/L");

  myEventT->Branch("rhoFJ", &rhoFJ, "rhoFJ/F");
  myEventT->Branch("rhoJetsFJ", &rhoJetsFJ, "rhoJetsFJ/F");
  myEventT->Branch("rhoJetsCentralFJ", &rhoJetsCentralFJ, "rhoJetsCentralFJ/F");
  myEventT->Branch("rhoJetsNoPUFJ", &rhoJetsNoPUFJ, "rhoJetsNoPUFJ/F");
  
  myEventT->Branch("N_PU", N_PU, "N_PU[3]/I");

  //Lepton info
  myEventT->Branch("MU_pt", MU_pt, "MU_pt[2]/F");
  myEventT->Branch("MU_eta", MU_eta, "MU_eta[2]/F");
  myEventT->Branch("MU_phi", MU_phi, "MU_phi[2]/F");
  myEventT->Branch("MU_vx", MU_pt, "MU_vx[2]/F");
  myEventT->Branch("MU_vy", MU_pt, "MU_vy[2]/F");
  myEventT->Branch("MU_vz", MU_pt, "MU_vz[2]/F");

  //true Z vertex
  myEventT->Branch("vZx", &vZx, "vZx/F");
  myEventT->Branch("vZy", &vZy, "vZy/F");
  myEventT->Branch("vZz", &vZz, "vZz/F");

  //reco vertex information
  myEventT->Branch("N_PV", &N_PV, "N_PV/I");
  myEventT->Branch("PV_x", PV_x, "Pv_x[5]/F");
  myEventT->Branch("PV_y", PV_y, "Pv_y[5]/F");
  myEventT->Branch("PV_z", PV_z, "Pv_z[5]/F");
  myEventT->Branch("PV_xerr", PV_xerr, "Pv_xerr[5]/F");
  myEventT->Branch("PV_yerr", PV_yerr, "Pv_yerr[5]/F");
  myEventT->Branch("PV_zerr", PV_zerr, "Pv_zerr[5]/F");
  myEventT->Branch("PV_SumPt", PV_SumPt, "Pv_SumPt[5]/F");
  myEventT->Branch("PV_ndof", PV_ndof, "Pv_ndof[5]/F");
  myEventT->Branch("PV_chi2", PV_chi2, "Pv_chi2[5]/F");

  //supercluster info
  myEventT->Branch("N_SC", &N_SC, "N_SC/I");
  myEventT->Branch("SC_nBC", SC_nBC, "SC_nBC[40]/I");
  myEventT->Branch("SC_nCrystal", SC_nCrystal, "SC_nCrystal[40]/I");
  myEventT->Branch("SC_rawEnergy", SC_rawEnergy, "SC_rawEnergy[40]/F");
  myEventT->Branch("SC_energy", SC_energy, "SC_energy[40]/F");
  myEventT->Branch("SC_eta", SC_eta, "SC_eta[40]/F");
  myEventT->Branch("SC_phi", SC_phi, "SC_phi[40]/F");
  myEventT->Branch("SC_phiWidth", SC_phiWidth, "SC_phiWidth[40]/F");
  myEventT->Branch("SC_etaWidth", SC_etaWidth, "SC_etaWidth[40]/F");
  myEventT->Branch("SC_time", SC_time, "SC_time[40]/F");
  myEventT->Branch("SC_chi2", SC_chi2, "SC_chi2[40]/F");
  myEventT->Branch("SC_x", SC_x, "SC_x[40]/F");
  myEventT->Branch("SC_y", SC_y, "SC_y[40]/F");
  myEventT->Branch("SC_z", SC_z, "SC_z[40]/F");

  //jet info

  //CaloJet
  myEventT->Branch("N_CaloJet", &N_CaloJet, "N_CaloJet/I");
  myEventT->Branch("CaloJet_pt", CaloJet_pt, "CaloJet_pt[20]/F");
  myEventT->Branch("CaloJet_eta", CaloJet_eta, "CaloJet_eta[20]/F");
  myEventT->Branch("CaloJet_phi", CaloJet_phi, "CaloJet_phi[20]/F");
  myEventT->Branch("CaloJet_energy", CaloJet_energy, "CaloJet_energy[20]/F");
  myEventT->Branch("CaloJet_uncorrEnergy", CaloJet_uncorrEnergy, "CaloJet_uncorrEnergy[20]/F");
  myEventT->Branch("CaloJet_vx", CaloJet_vx, "CaloJet_vx[20]/F");
  myEventT->Branch("CaloJet_vy", CaloJet_vy, "CaloJet_vy[20]/F");
  myEventT->Branch("CaloJet_vz", CaloJet_vz, "CaloJet_vz[20]/F");
  myEventT->Branch("CaloJet_area", CaloJet_area, "CaloJet_area[20]/F");
  myEventT->Branch("CaloJet_EMFrac", CaloJet_EMFrac, "CaloJet_EMFrac[20]/F");
  myEventT->Branch("CaloJet_Ntrack_match", CaloJet_Ntrack_match, "CaloJet_Ntrack_match[20]/I");
  myEventT->Branch("CaloJet_Ntrack_nomatch", CaloJet_Ntrack_nomatch, "CaloJet_Ntrack_nomatch[20]/I");
  myEventT->Branch("CaloJet_pt_match", CaloJet_pt_match, "CaloJet_pt_match[20]/F");
  myEventT->Branch("CaloJet_pt_nomatch", CaloJet_pt_nomatch, "CaloJet_pt_nomatch[20]/F");

  myEventT->Branch("CaloJet_Id", CaloJet_Id, "CaloJet_Id[20]/I");
  myEventT->Branch("CaloJet_covEtaEta", CaloJet_covEtaEta, "CaloJet_covEtaEta[20]/F");
  myEventT->Branch("CaloJet_covPhiPhi", CaloJet_covPhiPhi, "CaloJet_covPhiPhi[20]/F");

  //PFNoPU
  myEventT->Branch("N_PFNoPU", &N_PFNoPU, "N_PFNoPU/I");
  myEventT->Branch("PFNoPU_pt", PFNoPU_pt, "PFNoPU_pt[20]/F");
  myEventT->Branch("PFNoPU_eta", PFNoPU_eta, "PFNoPU_eta[20]/F");
  myEventT->Branch("PFNoPU_phi", PFNoPU_phi, "PFNoPU_phi[20]/F");
  myEventT->Branch("PFNoPU_energy", PFNoPU_energy, "PFNoPU_energy[20]/F");
  myEventT->Branch("PFNoPU_uncorrEnergy", PFNoPU_uncorrEnergy, "PFNoPU_uncorrEnergy[20]/F");
  myEventT->Branch("PFNoPU_vx", PFNoPU_vx, "PFNoPU_vx[20]/F");
  myEventT->Branch("PFNoPU_vy", PFNoPU_vy, "PFNoPU_vy[20]/F");
  myEventT->Branch("PFNoPU_vz", PFNoPU_vz, "PFNoPU_vz[20]/F");
  myEventT->Branch("PFNoPU_area", PFNoPU_area, "PFNoPU_area[20]/F");
  myEventT->Branch("PFNoPU_EMFrac", PFNoPU_EMFrac, "PFNoPU_EMFrac[20]/F");
  myEventT->Branch("PFNoPU_Ntrack_match", PFNoPU_Ntrack_match, "PFNoPU_Ntrack_match[20]/I");
  myEventT->Branch("PFNoPU_Ntrack_nomatch", PFNoPU_Ntrack_nomatch, "PFNoPU_Ntrack_nomatch[20]/I");
  myEventT->Branch("PFNoPU_pt_match", PFNoPU_pt_match, "PFNoPU_pt_match[20]/F");
  myEventT->Branch("PFNoPU_pt_nomatch", PFNoPU_pt_nomatch, "PFNoPU_pt_nomatch[20]/F");
  myEventT->Branch("PFNoPU_MVAID", PFNoPU_MVAID, "PFNoPU_MVAID[20]/F");
  myEventT->Branch("PFNoPU_LoosePFJetID", PFNoPU_LoosePFJetID, "PFNoPU_LoosePFJetID[20]/I");
  myEventT->Branch("PFNoPU_MediumPFJetID", PFNoPU_MediumPFJetID, "PFNoPU_MediumPFJetID[20]/I");
  myEventT->Branch("PFNoPU_TightPFJetID", PFNoPU_TightPFJetID, "PFNoPU_TightPFJetID[20]/I");

  //PFPUcorr
  myEventT->Branch("N_PFPUcorr", &N_PFPUcorr, "N_PFPUcorr/I");
  myEventT->Branch("PFPUcorr_pt", PFPUcorr_pt, "PFPUcorr_pt[20]/F");
  myEventT->Branch("PFPUcorr_eta", PFPUcorr_eta, "PFPUcorr_eta[20]/F");
  myEventT->Branch("PFPUcorr_phi", PFPUcorr_phi, "PFPUcorr_phi[20]/F");
  myEventT->Branch("PFPUcorr_energy", PFPUcorr_energy, "PFPUcorr_energy[20]/F");
  myEventT->Branch("PFPUcorr_uncorrEnergy", PFPUcorr_uncorrEnergy, "PFPUcorr_uncorrEnergy[20]/F");
  myEventT->Branch("PFPUcorr_vx", PFPUcorr_vx, "PFPUcorr_vx[20]/F");
  myEventT->Branch("PFPUcorr_vy", PFPUcorr_vy, "PFPUcorr_vy[20]/F");
  myEventT->Branch("PFPUcorr_vz", PFPUcorr_vz, "PFPUcorr_vz[20]/F");
  myEventT->Branch("PFPUcorr_area", PFPUcorr_area, "PFPUcorr_area[20]/F");
  myEventT->Branch("PFPUcorr_EMFrac", PFPUcorr_EMFrac, "PFPUcorr_EMFrac[20]/F");
  myEventT->Branch("PFPUcorr_Ntrack_match", PFPUcorr_Ntrack_match, "PFPUcorr_Ntrack_match[20]/I");
  myEventT->Branch("PFPUcorr_Ntrack_nomatch", PFPUcorr_Ntrack_nomatch, "PFPUcorr_Ntrack_nomatch[20]/I");
  myEventT->Branch("PFPUcorr_pt_match", PFPUcorr_pt_match, "PFPUcorr_pt_match[20]/F");
  myEventT->Branch("PFPUcorr_pt_nomatch", PFPUcorr_pt_nomatch, "PFPUcorr_pt_nomatch[20]/F");
  myEventT->Branch("PFPUcorr_MVAID", PFPUcorr_MVAID, "PFPUcorr_MVAID[20]/F");
  myEventT->Branch("PFPUcorr_LoosePFJetID", PFPUcorr_LoosePFJetID, "PFPUcorr_LoosePFJetID[20]/I");
  myEventT->Branch("PFPUcorr_MediumPFJetID", PFPUcorr_MediumPFJetID, "PFPUcorr_MediumPFJetID[20]/I");
  myEventT->Branch("PFPUcorr_TightPFJetID", PFPUcorr_TightPFJetID, "PFPUcorr_TightPFJetID[20]/I");

  //GEN
  myEventT->Branch("N_GEN", &N_GEN, "N_GEN/I");
  myEventT->Branch("GEN_pt", GEN_pt, "GEN_pt[20]/F");
  myEventT->Branch("GEN_eta", GEN_eta, "GEN_eta[20]/F");
  myEventT->Branch("GEN_phi", GEN_phi, "GEN_phi[20]/F");
  myEventT->Branch("GEN_energy", GEN_energy, "GEN_energy[20]/F");
  myEventT->Branch("GEN_Ntrack_match", GEN_Ntrack_match, "GEN_Ntrack_match[20]/I");
  myEventT->Branch("GEN_Ntrack_nomatch", GEN_Ntrack_nomatch, "GEN_Ntrack_nomatch[20]/I");
  myEventT->Branch("GEN_pt_match", GEN_pt_match, "GEN_pt_match[20]/F");
  myEventT->Branch("GEN_pt_nomatch", GEN_pt_nomatch, "GEN_pt_nomatch[20]/F");
}

void Timing::AddEvent(){
  //set info that is not already set by other functions
  RUN_NUM = runNumber;
  LS_NUM = lumiBlock;
  EVENT_NUM = eventNumber;
 
  MET_FLAGS = METFlags;
  TRACKER_FLAGS = tooManyTrackerFailures;

  rhoFJ = rhoFastjet;
  rhoJetsFJ = rhoJetsFastJet;
  rhoJetsCentralFJ = rhoJetsCentralNeutralFastJet;
  rhoJetsNoPUFJ = rhoJetsFastJet_nopu;

  for(int i = 0; i < 3; i++){
    N_PU[i] = nPU[i];
  }
 
  myEventT->Fill();
}

void Timing::WriteEventTree(){
  f_out->cd();
  myEventT->Write();
  myEventT->Delete();
}

bool Timing::PassLeptonSelection(){

  vMU.clear();

  vector<int> index_MU;

  if(nMuon < 2) return false;

  for(int i = 0; i < nMuon; i++){
    double pt = sqrt(pxMuon[i]*pxMuon[i]+pyMuon[i]*pyMuon[i]);
    if(pt > 20. && fabs(etaMuon[i]) < 2.1 && isTightMuon(i,false)){
      TLorentzVector mu;
      mu.SetPxPyPzE(pxMuon[i],pyMuon[i],pzMuon[i],energyMuon[i]);
      vMU.push_back(mu);
      index_MU.push_back(i);
      continue;
    } 
    if(pt > 15. && fabs(etaMuon[i]) < 2.1 && isLooseMuon(i,false)){
      TLorentzVector mu;
      mu.SetPxPyPzE(pxMuon[i],pyMuon[i],pzMuon[i],energyMuon[i]);
      vMU.push_back(mu);
      index_MU.push_back(i);
      continue;
    }
  }
  
  //exactly two muons
  if(vMU.size() != 2) return false;

  //opposite signs
  if(chargeMuon[index_MU[0]]+chargeMuon[index_MU[1]] != 0) return false;

  //Z window
  if( fabs((vMU[0]+vMU[1]).M()-91.) > 15.) return false;
  
  //fill lepton info
  for(int i = 0; i < 2; i++){
    double pt = sqrt(pxMuon[index_MU[i]]*pxMuon[index_MU[i]]+
		     pyMuon[index_MU[i]]*pyMuon[index_MU[i]]);
    MU_pt[i] = pt;
    MU_eta[i] = etaMuon[index_MU[i]];
    MU_phi[i] = phiMuon[index_MU[i]];
    MU_vx[i] = vertexXMuon[index_MU[i]];
    MU_vy[i] = vertexYMuon[index_MU[i]];
    MU_vz[i] = vertexZMuon[index_MU[i]];
  }

  bool goodZ = false;

  //Get Z gen-level info
  int iBos = 23;
  for(int i = 0; i < nMc; i++){
    if(abs(statusMc[i]) != 3) continue;
    if(mothMc[i] >= 0 && mothMc[i] < nMc){
      if(abs(idMc[mothMc[i]]) == iBos){

	goodZ = true;
	vZx = vxMc[mothMc[i]];
	vZy = vyMc[mothMc[i]];
	vZz = vzMc[mothMc[i]];

	break;
      }
    }
  }

  if(goodZ == false) return false;
  
  return true;
}

bool Timing::PassPVSelection(){
  if(nPV < 1) return false;

  //get Z-vertex
  TVector3 vZ;
  vZ.SetXYZ(vZx,vZy,vZz);

  N_PV = 0;

  //Now we find vertex with lowest/highest z, the one closest to the real Z and the ones 
  //below/above in z this one

  float min_z = 1000000000.;
  int i_min_z = -1;
  float max_z = -1000000000.;
  int i_max_z = -1;

  float min_dist = 1000000000000.;
  int i_min_dist = -1;

  Vertices.clear();

  for(int i = 0; i < nPV; i++){
    float rhoVTX = sqrt(PVxPV[i]*PVxPV[i] + PVyPV[i]*PVyPV[i]);
    if(fabs(PVzPV[i]) <= 24. && ndofPV[i] > 4 && isFakePV[i] == 0 && rhoVTX <= 2.){
      N_PV++;

      TVector3 vtemp;
      vtemp.SetXYZ(PVxPV[i],PVyPV[i],PVzPV[i]);

      if( (vtemp-vZ).Mag() < min_dist){
	min_dist = (vtemp-vZ).Mag();
	i_min_dist = i;
	myV = vtemp;
      }
      if(PVzPV[i] < min_z){
	min_z = PVzPV[i];
	i_min_z = i;
      }
      if(PVzPV[i] > max_z){
	max_z = PVzPV[i];
	i_max_z = i;
      }
    }
  }

  if(N_PV < 1) return false;

  float low_z = 1000000000.;
  int i_low_z = -1;
  float high_z = 1000000000.;
  int i_high_z = -1;
  for(int i = 0; i < nPV; i++){
    float rhoVTX = sqrt(PVxPV[i]*PVxPV[i] + PVyPV[i]*PVyPV[i]);
    if(fabs(PVzPV[i]) <= 24. && ndofPV[i] > 4 && isFakePV[i] == 0 && rhoVTX <= 2.){

      if(i == i_min_dist) continue;

      TVector3 vtemp;
      vtemp.SetXYZ(PVxPV[i],PVyPV[i],PVzPV[i]);
      Vertices.push_back(vtemp);

      if(PVzPV[i] < PVzPV[i_min_dist] && fabs(PVzPV[i]-PVzPV[i_min_dist]) < low_z){
	low_z = fabs(PVzPV[i]-PVzPV[i_min_dist]);
	i_low_z = i;
      }
      if(PVzPV[i] > PVzPV[i_min_dist] && fabs(PVzPV[i]-PVzPV[i_min_dist]) < high_z){
	high_z = fabs(PVzPV[i]-PVzPV[i_min_dist]);
	i_high_z = i;
      }
    }
  }

  //Now we fill the info
  vector<int> Vindex;
  if(i_min_z > -1){
    Vindex.push_back(i_min_z);
  } else {
    Vindex.push_back(i_min_dist);
  }
  if(i_low_z > -1){
    Vindex.push_back(i_low_z);
  } else {
    Vindex.push_back(i_min_dist);
  }
  Vindex.push_back(i_min_dist);
  if(i_high_z > -1){
    Vindex.push_back(i_high_z);
  } else {
    Vindex.push_back(i_min_dist);
  }
  if(i_max_z > -1){
    Vindex.push_back(i_max_z);
  } else {
    Vindex.push_back(i_min_dist);
  }

  for(int i = 0; i < 5; i++){
    PV_x[i] = PVxPV[Vindex[i]];
    PV_y[i] = PVyPV[Vindex[i]];
    PV_z[i] = PVzPV[Vindex[i]];
    PV_xerr[i] = PVErrxPV[Vindex[i]];
    PV_yerr[i] = PVErryPV[Vindex[i]];
    PV_zerr[i] = PVErrzPV[Vindex[i]];
    PV_SumPt[i] = SumPtPV[Vindex[i]];
    PV_ndof[i] = ndofPV[Vindex[i]];
    PV_chi2[i] = chi2PV[Vindex[i]];
  }

  return true;

}

void Timing::InitTracks(){

  matched_tracks.clear();
  unmatched_tracks.clear();

  //get the list of tracks matched/unmatched to the PV
 
  for(int i = 0; i < nTrack; i++){
    
    TVector3 vT(trackVxTrack[i],trackVyTrack[i],trackVzTrack[i]);
    TVector3 p;
    p.SetXYZ(pxTrack[i],pyTrack[i],pzTrack[i]);

    if(p.Pt() < 0.5) continue;
    if(p.Pt() > 500.) continue;
    if(trackNormalizedChi2Track[i] > 20.0) continue;
    if(fabs(transvImpactParTrack[i]/transvImpactParErrorTrack[i]) > 3.5) continue;
    if(fabs(p.Eta()) > 2.4) continue;
    if(trackValidHitsTrack[i] < 5) continue;

    bool unmatched = false;

    double D = (vT-myV).Mag();

    for(int iv = 0; iv < Vertices.size(); iv++){
      if( (vT-Vertices[iv]).Mag() < D){
	unmatched = true;
	break;
      }
    }

    if(fabs((vT-myV).Mag()) > 0.1) unmatched = true;
    
    if(unmatched){
      unmatched_tracks.push_back(p);
    } else {
      matched_tracks.push_back(p);
    }
  }
}

bool Timing::PassSCSelection(){

  if(nSC < 1) return false;
  

  N_SC = 0;

  for(int i = 0; i < nSC && i < 40; i++){
    SC_nBC[N_SC] = nBCSC[i];
    SC_nCrystal[N_SC] = nCrystalsSC[i];
    SC_rawEnergy[N_SC] = rawEnergySC[i];
    SC_energy[N_SC] = energySC[i];
    SC_eta[N_SC] = etaSC[i];
    SC_phi[N_SC] = phiSC[i];
    SC_phiWidth[N_SC] = phiWidthSC[i];
    SC_etaWidth[N_SC] = etaWidthSC[i];
    SC_time[N_SC] = timeSC[i];
    SC_chi2[N_SC] = chi2SC[i];
    SC_x[N_SC] = xPosSC[i];
    SC_y[N_SC] = yPosSC[i];
    SC_z[N_SC] = zPosSC[i];
    N_SC++;
  }

  return true;
}



bool Timing::PassJetSelection(){

  double min_pt = 20.;


  /////////////////////////////////
  // fill the different jet collections
  /////////////////////////////////

  /////////////////////////////////////
  // calo jets w/ PU+L2L3 correction
  /////////////////////////////////////
  N_CaloJet = 0;
  for(int i = 0; i < nAK5Jet && N_CaloJet < 20; i++){
    TLorentzVector jet;
    
    double px = pxAK5Jet[i];
    double py = pyAK5Jet[i];
    double pz = pzAK5Jet[i];
    double E = sqrt(px*px+py*py+pz*pz);
    double scale = 1.;
    jet.SetPxPyPzE(scale*px,scale*py,scale*pz,scale*E);

    if(jet.Pt() < min_pt) continue;

   
    if(nHit90AK5Jet[i] <= 1 || fHPDAK5Jet[i] >= 0.98 && fabs(jet.Eta()) < 3.0){
      //Noise jet - reject event
      return false;
    }
    if(fabs(jet.Eta()) < 2.55 && emFracAK5Jet[i] <= 0.01){
      //Noise jet - reject event
      return false;
    }
    if(fabs(jet.Eta()) >= 2.55 && jet.Pt() > 80. && emFracAK5Jet[i] >= 1. ){
      //Noise jet - reject event
      return false;
    }
    if(fabs(jet.Eta()) >= 2.55 && emFracAK5Jet[i] <= -0.9){
      //Noise jet - reject event
      return false;
    }
  

    //check muons' deltaR
    for(int imu = 0; imu < 2; imu++){
      if(jet.DeltaR(vMU[imu]) < 0.5) continue;
    }

    //Add jet info
    CaloJet_pt[N_CaloJet] = jet.Pt();
    CaloJet_eta[N_CaloJet] = jet.Eta();
    CaloJet_phi[N_CaloJet] = jet.Phi();
    CaloJet_energy[N_CaloJet] = energyAK5Jet[i];
    CaloJet_uncorrEnergy[N_CaloJet] = uncorrEnergyAK5Jet[i];
    CaloJet_vx[N_CaloJet] = vertexXAK5Jet[i];
    CaloJet_vy[N_CaloJet] = vertexYAK5Jet[i];
    CaloJet_vz[N_CaloJet] = vertexZAK5Jet[i];
    CaloJet_area[N_CaloJet] = areaAK5Jet[i];
    CaloJet_EMFrac[N_CaloJet] = emFracAK5Jet[i];
    CaloJet_Id[N_CaloJet] = IdAK5Jet[i];
    CaloJet_covEtaEta[N_CaloJet] = covEtaEtaAK5Jet[i];
    CaloJet_covPhiPhi[N_CaloJet] = covPhiPhiAK5Jet[i];
    
    //Now we do track/jet matching
    int N_match = 0;
    int N_unmatch = 0;
    float pt_match = 0.0;
    float pt_unmatch = 0.0;
    for(int it = 0; it < matched_tracks.size(); it++){
      if(matched_tracks[it].DeltaR(jet.Vect()) < 0.5){
	N_match++;
	pt_match += matched_tracks[it].Pt();
      }
    }
    for(int it = 0; it < unmatched_tracks.size(); it++){
      if(unmatched_tracks[it].DeltaR(jet.Vect()) < 0.5){
	N_unmatch++;
	pt_unmatch += unmatched_tracks[it].Pt();
      }
    }
    CaloJet_Ntrack_match[N_CaloJet] = N_match;
    CaloJet_Ntrack_nomatch[N_CaloJet] = N_unmatch;
    CaloJet_pt_match[N_CaloJet] = pt_match;
    CaloJet_pt_nomatch[N_CaloJet] = pt_unmatch;

    N_CaloJet++;
  }
  
  ///////////////////////////////////
  //Now, we do PFNoPU jets with PU corrections 
  ///////////////////////////////////
  N_PFNoPU = 0;
  for(int i = 0; i < nAK5PFNoPUJet && N_PFNoPU < 20; i++){
    TLorentzVector jet;
    double px = pxAK5PFNoPUJet[i];
    double py = pyAK5PFNoPUJet[i];
    double pz = pzAK5PFNoPUJet[i];
    double E = sqrt(px*px+py*py+pz*pz);
    double scale = 1.;
    jet.SetPxPyPzE(scale*px,scale*py,scale*pz,scale*E);
       
    if(jet.Pt() < min_pt) continue;

    //check muons' deltaR
    for(int imu = 0; imu < 2; imu++){
      if(jet.DeltaR(vMU[imu]) < 0.5) continue;
    }

    PFNoPU_pt[N_PFNoPU] = jet.Pt();
    PFNoPU_eta[N_PFNoPU] = jet.Eta();
    PFNoPU_phi[N_PFNoPU] = jet.Phi();
    PFNoPU_energy[N_PFNoPU] = energyAK5PFNoPUJet[i];
    PFNoPU_uncorrEnergy[N_PFNoPU] = uncorrEnergyAK5PFNoPUJet[i];
    PFNoPU_vx[N_PFNoPU] = vertexXAK5PFNoPUJet[i];
    PFNoPU_vy[N_PFNoPU] = vertexYAK5PFNoPUJet[i];
    PFNoPU_vz[N_PFNoPU] = vertexZAK5PFNoPUJet[i];
    PFNoPU_area[N_PFNoPU] = areaAK5PFNoPUJet[i];
    PFNoPU_EMFrac[N_PFNoPU] = (HFEMEnergyAK5PFNoPUJet[i]+
			       chargedEmEnergyAK5PFNoPUJet[i]+
			       neutralEmEnergyAK5PFNoPUJet[i])/uncorrEnergyAK5PFNoPUJet[i];
    PFNoPU_MVAID[N_PFNoPU] = jetIdMvaPhilV1AK5PFNoPUJet[i];
    PFNoPU_LoosePFJetID[N_PFNoPU] = isLoosePFNoPUJetID(i);
    PFNoPU_MediumPFJetID[N_PFNoPU] = isMediumPFNoPUJetID(i);
    PFNoPU_TightPFJetID[N_PFNoPU] = isTightPFNoPUJetID(i);
    
    //Now we do track/jet matching
    int N_match = 0;
    int N_unmatch = 0;
    float pt_match = 0.0;
    float pt_unmatch = 0.0;
    for(int it = 0; it < matched_tracks.size(); it++){
      if(matched_tracks[it].DeltaR(jet.Vect()) < 0.5){
	N_match++;
	pt_match += matched_tracks[it].Pt();
      }
    }
    for(int it = 0; it < unmatched_tracks.size(); it++){
      if(unmatched_tracks[it].DeltaR(jet.Vect()) < 0.5){
	N_unmatch++;
	pt_unmatch += unmatched_tracks[it].Pt();
      }
    }
    PFNoPU_Ntrack_match[N_PFNoPU] = N_match;
    PFNoPU_Ntrack_nomatch[N_PFNoPU] = N_unmatch;
    PFNoPU_pt_match[N_PFNoPU] = pt_match;
    PFNoPU_pt_nomatch[N_PFNoPU] = pt_unmatch;

    N_PFNoPU++;
  }

  ///////////////////////////////////
  //Now, we do PFPUcorr jets with PU corrections 
  ///////////////////////////////////
  N_PFPUcorr = 0;
  for(int i = 0; i < nAK5PFPUcorrJet && N_PFPUcorr < 20; i++){
    TLorentzVector jet;
    double px = pxAK5PFPUcorrJet[i];
    double py = pyAK5PFPUcorrJet[i];
    double pz = pzAK5PFPUcorrJet[i];
    double E = sqrt(px*px+py*py+pz*pz);
    double scale = 1.;
    jet.SetPxPyPzE(scale*px,scale*py,scale*pz,scale*E);
       
    if(jet.Pt() < min_pt) continue;

    //check muons' deltaR
    for(int imu = 0; imu < 2; imu++){
      if(jet.DeltaR(vMU[imu]) < 0.5) continue;
    }

    PFPUcorr_pt[N_PFPUcorr] = jet.Pt();
    PFPUcorr_eta[N_PFPUcorr] = jet.Eta();
    PFPUcorr_phi[N_PFPUcorr] = jet.Phi();
    PFPUcorr_energy[N_PFPUcorr] = energyAK5PFPUcorrJet[i];
    PFPUcorr_uncorrEnergy[N_PFPUcorr] = uncorrEnergyAK5PFPUcorrJet[i];
    PFPUcorr_vx[N_PFPUcorr] = vertexXAK5PFPUcorrJet[i];
    PFPUcorr_vy[N_PFPUcorr] = vertexYAK5PFPUcorrJet[i];
    PFPUcorr_vz[N_PFPUcorr] = vertexZAK5PFPUcorrJet[i];
    PFPUcorr_area[N_PFPUcorr] = areaAK5PFPUcorrJet[i];
    PFPUcorr_EMFrac[N_PFPUcorr] = (HFEMEnergyAK5PFPUcorrJet[i]+
			       chargedEmEnergyAK5PFPUcorrJet[i]+
			       neutralEmEnergyAK5PFPUcorrJet[i])/uncorrEnergyAK5PFPUcorrJet[i];
    PFPUcorr_MVAID[N_PFPUcorr] = jetIdMvaPhilV1AK5PFPUcorrJet[i];
    PFPUcorr_LoosePFJetID[N_PFPUcorr] = isLoosePFPUcorrJetID(i);
    PFPUcorr_MediumPFJetID[N_PFPUcorr] = isMediumPFPUcorrJetID(i);
    PFPUcorr_TightPFJetID[N_PFPUcorr] = isTightPFPUcorrJetID(i);
    
    //Now we do track/jet matching
    int N_match = 0;
    int N_unmatch = 0;
    float pt_match = 0.0;
    float pt_unmatch = 0.0;
    for(int it = 0; it < matched_tracks.size(); it++){
      if(matched_tracks[it].DeltaR(jet.Vect()) < 0.5){
	N_match++;
	pt_match += matched_tracks[it].Pt();
      }
    }
    for(int it = 0; it < unmatched_tracks.size(); it++){
      if(unmatched_tracks[it].DeltaR(jet.Vect()) < 0.5){
	N_unmatch++;
	pt_unmatch += unmatched_tracks[it].Pt();
      }
    }
    PFPUcorr_Ntrack_match[N_PFPUcorr] = N_match;
    PFPUcorr_Ntrack_nomatch[N_PFPUcorr] = N_unmatch;
    PFPUcorr_pt_match[N_PFPUcorr] = pt_match;
    PFPUcorr_pt_nomatch[N_PFPUcorr] = pt_unmatch;

    N_PFPUcorr++;
  }


  ///////////////////////////////////
  //Now, we do GEN jets:
  ///////////////////////////////////
     
  //GEN 
  N_GEN = 0;
  for(int i = 0; i < nAK5GenJet; i++){
    TLorentzVector jet;
    double px = pxAK5GenJet[i];
    double py = pyAK5GenJet[i];
    double pz = pzAK5GenJet[i];
    double E = sqrt(px*px+py*py+pz*pz);
    jet.SetPxPyPzE(px,py,pz,E);
    
    if(jet.Pt() < min_pt/2.) continue;

    GEN_pt[N_GEN] = jet.Pt();
    GEN_eta[N_GEN] = jet.Eta();
    GEN_phi[N_GEN] = jet.Phi();
    GEN_energy[N_GEN] = energyAK5GenJet[i];
    
    //Now we do track/jet matching
    int N_match = 0;
    int N_unmatch = 0;
    float pt_match = 0.0;
    float pt_unmatch = 0.0;
    for(int it = 0; it < matched_tracks.size(); it++){
      if(matched_tracks[it].DeltaR(jet.Vect()) < 0.5){
	N_match++;
	pt_match += matched_tracks[it].Pt();
      }
    }
    for(int it = 0; it < unmatched_tracks.size(); it++){
      if(unmatched_tracks[it].DeltaR(jet.Vect()) < 0.5){
	N_unmatch++;
	pt_unmatch += unmatched_tracks[it].Pt();
      }
    }
    GEN_Ntrack_match[N_GEN] = N_match;
    GEN_Ntrack_nomatch[N_GEN] = N_unmatch;
    GEN_pt_match[N_GEN] = pt_match;
    GEN_pt_nomatch[N_GEN] = pt_unmatch;

    N_GEN++;
  }
	
  if(N_CaloJet+N_PFNoPU+N_PFPUcorr < 1) return false;
 
  return true;
}
