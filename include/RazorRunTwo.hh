//-------------------------------------------------------
// Description:
// Razor analysis for 2Jets+MET final state.
// Keep information for btagging and leptons.
// Aimed to create the input of the data-driven
// prediction for RunII including the possibility of
// defining control and signal regions
// Authors:Cristian Pen~a
//
//-------------------------------------------------------

#ifndef RazorRunTwo_h
#define RazorRunTwo_h

//LOCAL INCLUDES
#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "JetCorrectionUncertainty.h"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "JetCorrectorParameters.h"
#include "ControlSampleEvents.hh"
#include "DefineStruct.hh"

using namespace std;

class RazorRunTwo : public Vecbos{
public:

  RazorRunTwo(TTree *tree = 0); /// Class Constructor
  RazorRunTwo(TTree *tree = 0, string jsonFile = string("none"), bool goodRunLS = false,
	      bool isData = false); /// Class Constructor
  RazorRunTwo(TTree *tree = 0, string jsonFile = string("none"), bool goodRunLS = false,
	      bool isData = false, bool keepNfiles = false);
  virtual ~RazorRunTwo();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void SetWeight(double);
  bool SetGenElectronIndex();
  bool SetGenMuonIndex();
  bool SetGenTauIndex();
  void SetGenLeptonVector();
  void SetGenPhotonVector();
  void ResetGenLeptonIndex();
  bool  DoPfSelection(std::vector<TLorentzVector>& pfJets, std::vector<int>& i_pfJets, std::vector< VecbosLepton > LooseLepton);
  void FillJetInfo(std::vector<TLorentzVector> GoodJets, std::vector<int> GoodJetIndices, std::vector<VecbosLepton> GoodLeptons);
  int FillPhotonInfo(int iPV);
  void SortByPt(std::vector<VecbosLepton>& lepton);
  void FillLeptons(std::vector<VecbosLepton> lepton);
  int  MatchLeptonGenLevel(TLorentzVector lepton);
  int  MatchPhotonGenLevel(TLorentzVector photon);
  void InitGenLeptonVariables();
  void InitLeptonVariables();
  float GetMTLep();
  void FillMTLep();
  double _weight;
  ControlSampleEvents* events;

private:
  int HighestPt(vector<TLorentzVector> p, int iHIGHEST);
  bool _isSMS;
  bool _isData;
  bool _goodRunLS;
  bool _keepNfiles;
  TTree* _treeCond;
  TFile* mu_corr_f;
  TH2F* mu_corr_h;
  TFile* pu_f;
  TH1D* pu_h;
  std::vector<int> genLeptonIndex;
};
#endif
