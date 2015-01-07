//-------------------------------------------------------
// Description:  CMS Hemisphere calculator
// Authors:      Maurizio Pierini, CERN
// Authors:      Christopher Rogan, Caltech
//-------------------------------------------------------

///  The CMSHemisphere class implements the 
///  algorithm to combine jets in hemispheres
///  as defined in CMS SUSY  searches

#ifndef CMSHemisphere_h
#define CMSHemisphere_h

#include <vector>
#include <iostream>
#include <TLorentzVector.h>

using namespace std;

class CMSHemisphere {

public:
  /// Constructor
  CMSHemisphere(vector<TLorentzVector> jets);

  /// Destructor
  ~CMSHemisphere();

  /// Return the hemisphere TLorentzVector;
  vector<TLorentzVector> GetHemispheres();

  vector<int> GetHem1Constituents();
  vector<int> GetHem2Constituents();

  /// Combining the jets in two hemispheres minimizing the 
  /// sum of the invariant masses of the two hemispheres
  void CombineMinMass();

  /// Combining the jets in two hemispheres minimizing the 
  /// difference of HT for the two hemispheres
  void CombineMinHT();

  /// Combining the jets in two hemispheres by minimizing m1^2/E1 + m2^2/E2
  void CombineMinEnergyMass();

  /// Combining the jets in two hemispheres by maximizing (E1-Beta*m1^2/E1 + E2-Beta*m1^2/E2)
  void CombineGeorgi();

private:

  /// Combine the jets in all the possible pairs of hemispheres
  void Combine();

  void CombineSaveConstituents();

  vector<TLorentzVector> jIN;
  vector<TLorentzVector> jOUT;
  vector<TLorentzVector> j1;
  vector<TLorentzVector> j2;
  int hem[20000][40];
  int hem2[20000][40];
  int chosen_perm;
  bool no_switch;
};

#endif
