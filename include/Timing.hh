//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

/// The Timing class can be used to perform fast check
/// on input ntuples (in a format compatible to VecbosBase)

#ifndef Timing_h
#define Timing_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"

#include "CommonTools/include/PUWeight.h"
#include "VecbosBase.hh"
#include "Vecbos.hh"

using namespace std;

class Timing : public Vecbos{
public:

  Timing(TTree *tree=0); /// Class Constructor
  virtual ~Timing();     /// Class Destructor
  /// The function to run on each events
  void Loop(string outFileName, int ISDATA);
  
private:

  //
  bool PassLeptonSelection();
  bool PassPVSelection();
  bool PassSCSelection();
  bool PassJetSelection();
  
  void InitTracks();

  void InitEventTree();
  void AddEvent();
  void WriteEventTree();
 
  //

  //Event tree for events of interest
  TTree *myEventT;

  string s_output;
 
  bool is_DATA;
 
  TFile *f_out;

  ////////////////////////////
  // Event stuff to be cleared
  ////////////////////////////
  //the two muons;
  vector<TLorentzVector> vMU;
  //reco vertex closest to Z vertex
  TVector3 myV;
  vector<TVector3> Vertices;
  vector<TVector3> matched_tracks;
  vector<TVector3> unmatched_tracks;

  ////////////////////////
  //Event tree content
  ////////////////////////
  //Run info
  int RUN_NUM;
  int LS_NUM;
  Long64_t EVENT_NUM;
  int MET_FLAGS;
  int TRACKER_FLAGS;
 
  //Lepton info
  float MU_pt[2];
  float MU_eta[2];
  float MU_phi[2];
  float MU_vx[2];
  float MU_vy[2];
  float MU_vz[2];

  //true Z vertex positon
  float vZx, vZy, vZz;
  
  //reco vertex information
  int N_PV;
  float PV_x[5];
  float PV_y[5];
  float PV_z[5];
  float PV_xerr[5];
  float PV_yerr[5];
  float PV_zerr[5];
  float PV_SumPt[5];
  float PV_ndof[5];
  float PV_chi2[5];

  //supercluster info
  int N_SC;
  int SC_nBC[40];
  int SC_nCrystal[40];
  float SC_rawEnergy[40];
  float SC_energy[40];
  float SC_eta[40];
  float SC_phi[40];
  float SC_phiWidth[40];
  float SC_etaWidth[40];
  float SC_time[40];
  float SC_chi2[40];
  float SC_x[40];
  float SC_y[40];
  float SC_z[40];
  
  //jet info
  //CaloJet
  int N_CaloJet;
  float CaloJet_pt[20];
  float CaloJet_eta[20];
  float CaloJet_phi[20];
  float CaloJet_energy[20];
  float CaloJet_uncorrEnergy[20];
  float CaloJet_vx[20];
  float CaloJet_vy[20];
  float CaloJet_vz[20];
  float CaloJet_area[20];
  float CaloJet_EMFrac[20];
  
  int CaloJet_Ntrack_match[20];
  int CaloJet_Ntrack_nomatch[20];
  float CaloJet_pt_match[20];
  float CaloJet_pt_nomatch[20];

  int CaloJet_Id[20];
  float CaloJet_covEtaEta[20];
  float CaloJet_covPhiPhi[20];

  //PFNoPU
  int N_PFNoPU;
  float PFNoPU_pt[20];
  float PFNoPU_eta[20];
  float PFNoPU_phi[20];
  float PFNoPU_energy[20];
  float PFNoPU_uncorrEnergy[20];
  float PFNoPU_vx[20];
  float PFNoPU_vy[20];
  float PFNoPU_vz[20];
  float PFNoPU_area[20];
  float PFNoPU_EMFrac[20];
  int PFNoPU_LoosePFJetID[20];
  int PFNoPU_MediumPFJetID[20];
  int PFNoPU_TightPFJetID[20];

  int PFNoPU_Ntrack_match[20];
  int PFNoPU_Ntrack_nomatch[20];
  float PFNoPU_pt_match[20];
  float PFNoPU_pt_nomatch[20];

  float PFNoPU_MVAID[20];

  //PFPUcorr
  int N_PFPUcorr;
  float PFPUcorr_pt[20];
  float PFPUcorr_eta[20];
  float PFPUcorr_phi[20];
  float PFPUcorr_energy[20];
  float PFPUcorr_uncorrEnergy[20];
  float PFPUcorr_vx[20];
  float PFPUcorr_vy[20];
  float PFPUcorr_vz[20];
  float PFPUcorr_area[20];
  float PFPUcorr_EMFrac[20];

  int PFPUcorr_Ntrack_match[20];
  int PFPUcorr_Ntrack_nomatch[20];
  float PFPUcorr_pt_match[20];
  float PFPUcorr_pt_nomatch[20];
  
  float PFPUcorr_MVAID[20];
  int PFPUcorr_LoosePFJetID[20];
  int PFPUcorr_MediumPFJetID[20];
  int PFPUcorr_TightPFJetID[20];

  //other info for the event
  int N_PU[3];
  float rhoFJ;
  float rhoJetsFJ;
  float rhoJetsCentralFJ;
  float rhoJetsNoPUFJ;

  //GEN
  int N_GEN;
  float GEN_pt[20];
  float GEN_eta[20];
  float GEN_phi[20];
  float GEN_energy[20];

  int GEN_Ntrack_match[20];
  int GEN_Ntrack_nomatch[20];
  float GEN_pt_match[20];
  float GEN_pt_nomatch[20];
};
#endif
