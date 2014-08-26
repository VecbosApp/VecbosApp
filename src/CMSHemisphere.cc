#include <vector>
#include <math.h>
#include <TLorentzVector.h>
#include "CMSHemisphere.hh"

using namespace std;

CMSHemisphere::CMSHemisphere(vector<TLorentzVector> jets){ 
  if(jets.size() < 2) cout << "Error in CMSHemisphere: you should provide at least two jets to form Henispheres" << endl;
  jIN = jets;
  CombineSaveConstituents();
  //Combine();
}

CMSHemisphere::~CMSHemisphere() {
}

vector<TLorentzVector> CMSHemisphere::GetHemispheres() {
  return jOUT;
}

vector<int> CMSHemisphere::GetHem1Constituents() {
  vector<int> hem_temp;
  if(no_switch) {
    for (int i; i < 40; i++){
      hem_temp.push_back(hem[chosen_perm][i]);
    }
    return hem_temp;
    delete hem;
  }
  else {
    for (int i; i < 40; i++){
      hem_temp.push_back(hem2[chosen_perm][i]);
    }
    return hem_temp;
    delete hem2;
  }
}
vector<int> CMSHemisphere::GetHem2Constituents() {
  vector<int> hem_temp;
  if(no_switch) {
    for (int i; i < 40; i++){
      hem_temp.push_back(hem2[chosen_perm][i]);
    }
    return hem_temp;
    delete hem2;
  }
  else {
    for(int i;i < 40;i++){
      hem_temp.push_back(hem[chosen_perm][i]);
    }   
    return hem_temp;
    delete hem;
  }
}
void CMSHemisphere::Combine() {
  int N_JETS = jIN.size();

  int N_comb = 1;
  for(int i = 0; i < N_JETS; i++){
    N_comb *= 2;
  }
    
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
	j_temp1 += jIN[count];
      } else {
	j_temp2 += jIN[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }
    j1.push_back(j_temp1);
    j2.push_back(j_temp2);
  }
}

void CMSHemisphere::CombineMinMass() {
  double M_min = -1;
  // default value (in case none is found)
  TLorentzVector myJ1 = j1[0];
  TLorentzVector myJ2 = j2[0];
  for(int i=0; i< j1.size(); i++) {
    double M_temp = j1[i].M2()+j2[i].M2();
    if(M_min < 0 || M_temp < M_min){
      M_min = M_temp;
      myJ1 = j1[i];
      myJ2 = j2[i];
      chosen_perm = i;
    }
  }
  //  myJ1.SetPtEtaPhiM(myJ1.Pt(),myJ1.Eta(),myJ1.Phi(),0.0);
  //  myJ2.SetPtEtaPhiM(myJ2.Pt(),myJ2.Eta(),myJ2.Phi(),0.0);

  jOUT.clear();
  if(myJ1.Pt() > myJ2.Pt()){
    no_switch = true;
    jOUT.push_back(myJ1);
    jOUT.push_back(myJ2);
  } else {
    no_switch = false;
    jOUT.push_back(myJ2);
    jOUT.push_back(myJ1);
  }
}

void CMSHemisphere::CombineMinEnergyMass() {
  double M_min = -1;
  // default value (in case none is found)
  TLorentzVector myJ1 = j1[0];
  TLorentzVector myJ2 = j2[0];
  for(int i=0; i< j1.size(); i++) {
    double M_temp = j1[i].M2()/j1[i].E()+j2[i].M2()/j2[i].E();
    if(M_min < 0 || M_temp < M_min){
      M_min = M_temp;
      myJ1 = j1[i];
      myJ2 = j2[i];
      chosen_perm = i;
    }
  }
  
  //  myJ1.SetPtEtaPhiM(myJ1.Pt(),myJ1.Eta(),myJ1.Phi(),0.0);
  //  myJ2.SetPtEtaPhiM(myJ2.Pt(),myJ2.Eta(),myJ2.Phi(),0.0);

  jOUT.clear();
  if(myJ1.Pt() > myJ2.Pt()){
    no_switch = true;
    jOUT.push_back(myJ1);
    jOUT.push_back(myJ2);
  } else {
    no_switch = false;
    jOUT.push_back(myJ2);
    jOUT.push_back(myJ1);
  }
}

void CMSHemisphere::CombineGeorgi(){
  double M_max = -10000;
  // default value (in case none is found)
  TLorentzVector myJ1 = j1[0];
  TLorentzVector myJ2 = j2[0];
  for(int i=0; i< j1.size(); i++) {
    int myBeta = 2;
    double M_temp = (j1[i].E()-myBeta*j1[i].M2()/j1[i].E())+(j2[i].E()-myBeta*j2[i].M2()/j2[i].E());
    if(M_max < -9999 || M_temp > M_max){
      M_max = M_temp;
      myJ1 = j1[i];
      myJ2 = j2[i];
      chosen_perm = i;
    }
  }
  
  //  myJ1.SetPtEtaPhiM(myJ1.Pt(),myJ1.Eta(),myJ1.Phi(),0.0);
  //  myJ2.SetPtEtaPhiM(myJ2.Pt(),myJ2.Eta(),myJ2.Phi(),0.0);

  jOUT.clear();
  if(myJ1.Pt() > myJ2.Pt()){
    no_switch = true;
    jOUT.push_back(myJ1);
    jOUT.push_back(myJ2);
  } else {
    no_switch = false;
    jOUT.push_back(myJ2);
    jOUT.push_back(myJ1);
  }
}

void CMSHemisphere::CombineSaveConstituents() {
  int N_JETS = jIN.size();
  // jets NEED to be ordered by pT
  // saves 2D array, first index = # of combination of jets
  // saved number = jet # that goes into hemisphere
  int N_comb = 1;
  for(int i = 0; i < N_JETS; i++){
    N_comb *= 2;
  }
  for (int i = 0; i < 20000; i++){
    for (int j = 0; j < 40; j++){
      hem[i][j] = -1;
      hem2[i][j] = -1;
    }
  }
  int j_count;
  int array_count=0;
  for(int i = 1; i < N_comb-1; i++){
    TLorentzVector j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    int hem_count = 0; //loops through # of jets in first hem
    int hem_count2 = 0; // loops through # of jets in second hem
    while(j_count > 0){
      if(itemp/j_count == 1){
        j_temp1 += jIN[count];
	hem[array_count][hem_count] = count;
	hem_count++;
      } else {
        j_temp2 += jIN[count];
	hem2[array_count][hem_count2]=count;
	hem_count2++; 
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }
    j1.push_back(j_temp1);
    j2.push_back(j_temp2);
    array_count++;
  }
}

void CMSHemisphere::CombineMinHT() {
  double dHT_min = 999999999999999.0;
  // default value (in case none is found)
  TLorentzVector myJ1 = j1[0];
  TLorentzVector myJ2 = j2[0];
  for(int i=0; i< j1.size(); i++) {
    double dHT_temp = fabs(j1[i].E()-j2[i].E());
    if(dHT_temp < dHT_min){  
      dHT_min = dHT_temp;
      myJ1 = j1[i];
      myJ2 = j2[i];
    }
  }
  
  jOUT.clear();
  if(myJ1.Pt() > myJ2.Pt()){
    jOUT.push_back(myJ1);
    jOUT.push_back(myJ2);
  } else {
    jOUT.push_back(myJ2);
    jOUT.push_back(myJ1);
  }
}
