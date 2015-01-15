//-------------------------------------------------------
// Description:
// Define Structs used by RazorRunTwo App 
// Authors:Cristian Pen~a
//
//-------------------------------------------------------

#ifndef DEFINESTRUCT_HH
#define DEFINESTRUCT_HH

//ROOT INCLUDES
#include "TLorentzVector.h"

struct VecbosMuon{
  int index;
  TLorentzVector muon;
  bool _isLoose;
  bool _isTight;
};

struct VecbosEle{
  int index;
  TLorentzVector ele;
  bool _isLoose;
  bool _isTight;
};

struct VecbosLepton{
  int index;
  TLorentzVector lepton;
  int charge;
  double mass;
  int pdgID;//absolute value
  bool _isLoose;
  bool _isTight;
};

#endif
