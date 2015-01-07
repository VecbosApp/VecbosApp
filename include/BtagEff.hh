//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

/// The TestAnalysis class can be used to perform fast check
/// on input ntuples (in a format compatible to VecbosBase)

#ifndef BtagEff_h
#define BtagEff_h

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"

#include "VecbosBase.hh"
#include "Vecbos.hh"

#define JetFlavor_BJet 0
#define JetFlavor_CJet 1
#define JetFlavor_LJet 2
#define PI 3.14159265358979323846264338327950288479716939937510

using namespace std;

class BtagEff : public Vecbos{

public:
  BtagEff(TTree *tree = 0); /// Class Constructor
  virtual ~BtagEff();      /// Class Destructor
  /// The function to run on each events
  void Loop(string outFileName);
  
private:
  // B-tagging stuff
  
  int DetermineJetFlavor(double JetEta, double JetPhi);
  double GetDR(double Eta1, double Phi1, double Eta2, double Phi2);
  bool SwitchTag(bool CurrentTag, double ScaleFactor, double Efficiency);
  bool SwitchJetTag(bool CurrentTag, int JetFlavor, double JetEta, double JetPT,
		    bool DoFastSim, double Shift, int Model);
  int FindBin(double Value, const vector<double> &Bins);
  
  map<string, int> ModelIndex;
  
  vector<double> B_SF_PTBins, B_SF_EtaBins;
  vector<double> C_SF_PTBins, C_SF_EtaBins;
  vector<double> L_SF_PTBins, L_SF_EtaBins;
  double B_SF[500][5], C_SF[500][5], L_SF[500][5];   // Lazy coder here!!
  double B_SF_Error[500][5], C_SF_Error[500][5], L_SF_Error[500][5];
  
  vector<double> B_Eff_PTBins[1], B_Eff_EtaBins[1];
  vector<double> C_Eff_PTBins[1], C_Eff_EtaBins[1];
  vector<double> L_Eff_PTBins[1], L_Eff_EtaBins[1];
  double B_Eff[1][15000][5], C_Eff[1][15000][5], L_Eff[1][15000][5];   // Lazy coder here!!
  
  vector<double> B_SFFast_PTBins, B_SFFast_EtaBins;
  vector<double> C_SFFast_PTBins, C_SFFast_EtaBins;
  vector<double> L_SFFast_PTBins, L_SFFast_EtaBins;
  double B_SFFast[500][5], C_SFFast[500][5], L_SFFast[500][5];   // Lazy...
  double B_SFFast_Error[500][5], C_SFFast_Error[500][5], L_SFFast_Error[500][5];
  
  vector<double> B_EffFast_PTBins[1], B_EffFast_EtaBins[1];
  vector<double> C_EffFast_PTBins[1], C_EffFast_EtaBins[1];
  vector<double> L_EffFast_PTBins[1], L_EffFast_EtaBins[1];
  double B_EffFast[1][15000][5], C_EffFast[1][15000][5], L_EffFast[1][15000][5];   // Lazy...
  
  void ReadFile(vector<double> &PT, vector<double> &Eta,
		double Array[15000][5], double Error[15000][5], string FileName,
		bool DoubleError = false);
  void InitializeAllArrays();
  
};

#endif

