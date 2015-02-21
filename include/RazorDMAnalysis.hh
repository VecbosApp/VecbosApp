//-------------------------------------------------------
// Description:
//    Razor analysis for 2Jets+MET final state
// Authors:Cristian Pen~a
//
//-------------------------------------------------------

#ifndef RazorDMAnalysis_h
#define RazorDMAnalysis_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "JetCorrectionUncertainty.h"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "JetCorrectorParameters.h"

using namespace std;

class RazorDMAnalysis : public Vecbos{
public:

  RazorDMAnalysis(TTree *tree=0); /// Class Constructor
  RazorDMAnalysis(TTree *tree=0, string jsonFile=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~RazorDMAnalysis();     /// Class Destructor
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);
  void SetWeight(double);
  void SetLambdaEntries();
  void CreateHistoEntries();
  float GetQtr();
  double _weight;

private:
  int HighestPt(vector<TLorentzVector> p, int iHIGHEST);
  bool _isSMS;
  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  TFile* mu_corr_f;
  TH2F* mu_corr_h;
  TH1F* h_entries;//Lambda cut on Nentries
  TFile* pu_f;
  TH1D* pu_h;
};
#endif
