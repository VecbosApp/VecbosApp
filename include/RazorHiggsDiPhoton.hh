//-------------------------------------------------------
// Description:
//    Class for Test search analyses
// Authors:
//
//-------------------------------------------------------

#ifndef RazorHiggsDiPhoton_h
#define RazorHiggsDiPhoton_h

#include "VecbosBase.hh"
#include "Vecbos.hh"
#include "Jet.hh"
#include "CaloTower.hh"
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "CMSHemisphere.hh"

using namespace std;

class RazorHiggsDiPhoton : public Vecbos{
public:

  RazorHiggsDiPhoton(TTree *tree=0); /// Class Constructor
  RazorHiggsDiPhoton(TTree *tree=0, string json=string("none"), bool goodRunLS=false, bool isData=false); /// Class Constructor
  virtual ~RazorHiggsDiPhoton();     /// Class Destructor
  void SetWeight(double weight);
  void Loop(string outFileName, int start, int stop);
  void SetConditions(TTree* treeCond);

private:
  bool PhotonIdBarrel(int iPh, TString Selector);
  bool PhotonIdBarrelisEM(int i);
  bool PhotonIdBarrelisLoose(int i);
  bool PhotonIdBarrelisTight(int i);
  int  AntiPhotonIdBarrel(int iPh);
  bool IsPhotonBarrel(int iPh);

  bool _isData;
  bool _goodRunLS;
  TTree* _treeCond;
  double _weight;

};
#endif
