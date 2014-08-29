// std includes
#include <iostream>
#include <string>
#include <vector>

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
#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/TriggerMask.hh"
#include "RazorHiggsDiPhoton.hh"

RazorHiggsDiPhoton::RazorHiggsDiPhoton(TTree *tree) : Vecbos(tree) {
  _goodRunLS = false;
  _isData = false;
  _weight = 1.;
}

RazorHiggsDiPhoton::RazorHiggsDiPhoton(TTree *tree, string json, bool goodRunLS, bool isData) : Vecbos(tree) {

  _goodRunLS = goodRunLS;
  _isData = isData;
  _weight = 1.;

  //To read good run list!
  if (goodRunLS && isData) {
    setJsonGoodRunList(json);
    fillRunLSMap();
  }
}

RazorHiggsDiPhoton::~RazorHiggsDiPhoton() {}

void RazorHiggsDiPhoton::SetConditions(TTree* treeCond) {
  _treeCond = treeCond;
}

void RazorHiggsDiPhoton::SetWeight(double weight) {
  _weight = weight;
}

void RazorHiggsDiPhoton::Loop(string outFileName, int start, int stop) {
  if(fChain == 0) return;

  // double-photon triggers
  int HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50;


  // photon block
  double pTPh1;
  double etaPh1;
  double phiPh1;
  int    noPhoton1;
  int    isEMPh1;
  int    isLoosePh1;
  int    isTightPh1;
  double pTPh2;
  double etaPh2;
  double phiPh2;
  int    noPhoton2;
  int    isEMPh2;
  int    isLoosePh2;
  int    isTightPh2;
  double pTGG;
  double etaGG;
  double phiGG;
  double mGG;

  // PF  block
  int    passedPF;
  double pTPFHem1;
  double etaPFHem1;
  double phiPFHem1;
  double pTPFHem2;
  double etaPFHem2;
  double phiPFHem2;
  double RSQ;
  double MR;

  // calo block
  int    passedCalo;
  double pTCaloHem1;
  double etaCaloHem1;
  double phiCaloHem1;
  double pTCaloHem2;
  double etaCaloHem2;
  double phiCaloHem2;
  double CaloRSQ;
  double CaloMR;

  // general event info
  int run;
  ULong64_t evNum;
  int bx;
  int lumi;
  int orbit;
  double W;
  double pTPFJet[10];
  double etaPFJet[10];
  double phiPFJet[10];
  double mPFJet[10];
  int constPFHem1[10];
  int constPFHem2[10];

  // prepare the output tree
  TTree* outTree = new TTree("outTree", "outTree");
  outTree->Branch("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50", &HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50, "HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50/I");

  outTree->Branch("run", &run, "run/I");
  outTree->Branch("evNum", &evNum, "evNum/l");
  outTree->Branch("bx", &bx, "bx/I");
  outTree->Branch("lumi", &lumi, "lumi/I");
  outTree->Branch("orbit", &orbit, "orbit/I");
  outTree->Branch("W", &W, "W/D");
    
  outTree->Branch("pTPh1", &pTPh1, "pTPh1/D");
  outTree->Branch("etaPh1", &etaPh1, "etaPh1/D");
  outTree->Branch("phiPh1", &phiPh1, "phiPh1/D");
  outTree->Branch("noPhoton1", &noPhoton1, "noPhoton1/I");
  outTree->Branch("isEMPh1", &isEMPh1, "isEMPh1/I");
  outTree->Branch("isLoosePh1", &isLoosePh1, "isLoosePh1/I");
  outTree->Branch("isTightPh1", &isTightPh1, "isTightPh1/I");
  outTree->Branch("pTPh2", &pTPh2, "pTPh2/D");
  outTree->Branch("etaPh2", &etaPh2, "etaPh2/D");
  outTree->Branch("phiPh2", &phiPh2, "phiPh2/D");
  outTree->Branch("noPhoton2", &noPhoton2, "noPhoton2/I");
  outTree->Branch("isEMPh2", &isEMPh2, "isEMPh2/I");
  outTree->Branch("isLoosePh2", &isLoosePh2, "isLoosePh2/I");
  outTree->Branch("isTightPh2", &isTightPh2, "isTightPh2/I");
  outTree->Branch("pTGG", &pTGG, "pTGG/D");
  outTree->Branch("phiGG", &phiGG, "phiGG/D");
  outTree->Branch("etaGG", &etaGG, "etaGG/D");
  outTree->Branch("mGG", &mGG, "mGG/D");
    
  // PF block
  outTree->Branch("passedPF", &passedPF, "passedPF/I");
  outTree->Branch("pTPFHem1", &pTPFHem1, "pTPFHem1/D");
  outTree->Branch("etaPFHem1", &etaPFHem1, "etaPFHem1/D");
  outTree->Branch("phiPFHem1", &phiPFHem1, "phiPFHem1/D");
  outTree->Branch("constPFHem1", constPFHem1, "constPFHem1[10]/I");
  outTree->Branch("pTPFHem2", &pTPFHem2, "pTPFHem2/D");
  outTree->Branch("etaPFHem2", &etaPFHem2, "etaPFHem2/D");
  outTree->Branch("phiPFHem2", &phiPFHem2, "phiPFHem2/D");
  outTree->Branch("constPFHem2", constPFHem2, "constPFHem2[10]/I");
  outTree->Branch("pTPFJet", pTPFJet, "pTPFJet[10]/D");
  outTree->Branch("etaPFJet", etaPFJet, "etaPFJet[10]/D");
  outTree->Branch("phiPFJet", phiPFJet, "phiPFJet[10]/D");
  outTree->Branch("mPFJet", mPFJet, "mPFJet[10]/D");
  outTree->Branch("RSQ", &RSQ, "RSQ/D");
  outTree->Branch("MR", &MR, "MR/D");
    
  // Calo block
  outTree->Branch("passedCalo", &passedCalo, "passedCalo/I");
  outTree->Branch("pTCaloHem1", &pTCaloHem1, "pTCaloHem1/D");
  outTree->Branch("etaCaloHem1", &etaCaloHem1, "etaCaloHem1/D");
  outTree->Branch("phiCaloHem1", &phiCaloHem1, "phiCaloHem1/D");
  outTree->Branch("pTCaloHem2", &pTCaloHem2, "pTCaloHem2/D");
  outTree->Branch("etaCaloHem2", &etaCaloHem2, "etaCaloHem2/D");
  outTree->Branch("phiCaloHem2", &phiCaloHem2, "phiCaloHem2/D");
  outTree->Branch("CaloRSQ", &CaloRSQ, "CaloRSQ/D");
  outTree->Branch("CaloMR", &CaloMR, "CaloMR/D");


  std::vector<std::string> maskHLT; 
  maskHLT.push_back("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50");

  // prepare vectors for efficiency (variables for effTree)
  double Npassed_In = 0;
  double Npassed_HLT = 0; 
  double Npassed_PV = 0;
  double Npassed_2Ph = 0;
  double Npassed_PFHem = 0;
  double Npassed_CaloHem = 0;

  //  double _weight = 1.;
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  Long64_t nbytes = 0;
  Long64_t nb = 0;
  cout << "Number of entries = " << stop << endl;
  for (Long64_t jentry=start;  jentry<stop;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) cout << ">>> Processing event # " << jentry << endl;

    //Good Run selection
    if (_isData && _goodRunLS && !isGoodRunLS()) {
      if ( lastRun != runNumber || lastLumi != lumiBlock) {
	//std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (_isData && _goodRunLS && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      //std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }

    cout << "\n" << runNumber << " : " << lumiBlock << " : " << eventNumber << endl;
    
    Npassed_In += _weight;

    // HLT 


    //IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    if(_isData) {
      setRequiredTriggers(maskHLT); reloadTriggerMask(true); 

      HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50 = hasPassedHLT();
      /*
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
      */
    }

    Npassed_HLT += _weight;

    // find highest-pT PV
    int iHighestPt = -99;
    double HighestPt = -99999.;
    if(nPV<1) continue;

    for(int i=0; i< nPV; i++) if(SumPtPV[i] > HighestPt) iHighestPt = i;
    // PV selection
    if(ndofPV[iHighestPt] < 3) continue;
    if(PVzPV[iHighestPt] > 25.) continue; 

    Npassed_PV += _weight;

    int goodPFevent = true;
    int goodCaloevent = true;

    // HCAL FLAGS
    if(!eventPassHcalFilter()) continue;


    for (int i = 0; i < 10; i++){
      pTPFJet[i] = -1;
      etaPFJet[i] = -1;
      phiPFJet[i] = -1;
      mPFJet[i] = -1;
      constPFHem1[i] = -1;
      constPFHem2[i] = -1;
    }

    
    vector<TLorentzVector> CaloJet;    
    for(int i = 0; i < nAK5Jet; i++){
      TLorentzVector jet;
      double px = pxAK5Jet[i];
      double py = pyAK5Jet[i];
      double pz = pzAK5Jet[i];
      double E = sqrt(px*px+py*py+pz*pz);
      jet.SetPxPyPzE(px,py,pz,E);
      
      if(jet.Pt() < 40. || fabs(jet.Eta()) >= 3.0) continue; //jet kinematic cuts
      CaloJet.push_back(jet);
      
      // if(nHit90AK5Jet[i] <= 1 || fHPDAK5Jet[i] >= 0.98){ //fails jet id
      // 	goodCaloevent = false;
      // }
      // if(fabs(jet.Eta()) < 2.55 && emFracAK5Jet[i] <= 0.01){ //fails jet id
      // 	goodCaloevent = false;
      // }
      // if(fabs(jet.Eta()) >= 2.55){ //fails jet ID
      // 	if(jet.Pt() > 80. && emFracAK5Jet[i] >= 1.){
      // 	  goodCaloevent = false;	
      // 	}
      // 	if(emFracAK5Jet[i] <= -0.9){
      // 	  goodCaloevent = false;
      // 	}
      //      }
    }

    vector<TLorentzVector> PFPhoton;
    vector<bool> isEM;
    vector<bool> isLoose;
    vector<bool> isTight;
    vector<bool> notPhoton;
    for(int i=0; i< nPho; i++) {
      // ET cut
      double pT = energyPho[i]*sin(thetaPho[i]);
      if(pT  >= 25.) {

	/*	
	  bool isEM_i = true;
	  bool isLoose_i = true;
	  bool isTight_i = true;
	  PhotonIdBarrel(i,isEM_i,isLoose_i,isTight_i);  ////// rearrange
	*/

	bool isEM_i = false;
	bool isLoose_i = false;
	bool isTight_i = false;
	
	if ( PhotonIdBarrelisEM(i)){ isEM_i = true;} // == true
	if ( PhotonIdBarrelisLoose(i)){ isLoose_i = true;} // == true
	if ( PhotonIdBarrelisTight(i)){ isTight_i = true;} // == true

	//if isTight is OK isEM and isLoose too//

	//cout<<" isEM = "<<isEM_i<<"//// isLoose= "<<isLoose_i<<"//// isTight= "<<isTight_i<<endl; 

	bool notPhoton_i = AntiPhotonIdBarrel(i);
	// e2/e9
	int iSC = superClusterIndexPho[i];
	double e2OVERe9 = (eMaxSC[iSC]+e2ndSC[iSC])/e3x3SC[iSC];
	//	if(e2OVERe9 <= 0.95) {
	// one jet OR timing [for no-jet events]
	//	  double time = fabs(timeSC[superClusterIndexPho[i]]);
	//	  if(PFPUcorrJet.size()+CaloJet.size() > 0 || time < 3.) {
	TLorentzVector myPhoton;
	myPhoton.SetPtEtaPhiE(energyPho[i]*sin(thetaPho[i]), etaPho[i], phiPho[i], energyPho[i]);
	PFPhoton.push_back(myPhoton);
	isLoose.push_back(isLoose_i);
	isTight.push_back(isTight_i);
	isEM.push_back(isEM_i);
	notPhoton.push_back(notPhoton_i);
	//	  }
	//	}
      }
    }

    if(PFPhoton.size() < 2 ) continue; // keep at least 2 objects

    // select two highest-pT photons
    // and select what kind of photons we are : isEM[i] or isLoose[i] or isTight[i] 
    int iPh1 = -99;
    int iPh2 = -99;
    double phMax1 = 0.;
    double phMax2 = 0.;
    for(int i=0; i<PFPhoton.size(); i++) {
      if(PFPhoton[i].Pt() > phMax1 && isEM[i]) {
	phMax1 = PFPhoton[i].Pt();
	iPh1 = i;
      }
    }

    for(int i=0; i<PFPhoton.size(); i++) {
      if(PFPhoton[i].Pt() > phMax2 && i != iPh1 && isEM[i]) {
	phMax2 = PFPhoton[i].Pt();
	iPh2 = i;
      }
    }
    if (runNumber==194912){
      // VERY SPECIFIC TO THESE 7 EVENTS
      iPh1 = 0;
      iPh2 = 2;
    }
    if (runNumber==199436){
      // VERY SPECIFIC TO THESE 7 EVENTS
      iPh1 = 1;
      iPh2 = 2;
    }



    // Jet selection 
    vector<TLorentzVector> PFPUcorrJet;    
    for(int i=0; i< nAK5PFPUcorrJet; i++) {
      TLorentzVector myJet(pxAK5PFPUcorrJet[i], pyAK5PFPUcorrJet[i], pzAK5PFPUcorrJet[i], energyAK5PFPUcorrJet[i]);   
      if(myJet.Pt()>30. && fabs(myJet.Eta())< 3.0 && isLoosePFPUcorrJetID(i)) 
	if(myJet.DeltaR(PFPhoton[iPh1]) > 0.5 && myJet.DeltaR(PFPhoton[iPh2]) > 0.5)
	PFPUcorrJet.push_back(myJet);
    }


    // int flag = 1;
    // //bubble sort the jets by Pt()
    // for (int i=0; (i < PFPUcorrJet.size()) && flag; i++){
    //   TLorentzVector tempvector;
    //   flag = 0;
    //   for (int j=0; j < (PFPUcorrJet.size()-1); j++){
    // 	if (PFPUcorrJet.at(j+1).Pt() > PFPUcorrJet.at(j).Pt()){
    // 	  tempvector = PFPUcorrJet.at(j);
    // 	  PFPUcorrJet.at(j) = PFPUcorrJet.at(j+1);
    // 	  PFPUcorrJet.at(j+1) = tempvector;
    // 	  flag=1;
    // 	}
    //   }
    // }


    for (int i=0; i < PFPUcorrJet.size(); i++){
      pTPFJet[i] = PFPUcorrJet.at(i).Pt();
      etaPFJet[i] = PFPUcorrJet.at(i).Eta();
      phiPFJet[i] = PFPUcorrJet.at(i).Phi();
      mPFJet[i] = PFPUcorrJet.at(i).M();
    }

    

    
    cout << "Number of jets (pT > 40, |eta| < 3) = " << PFPUcorrJet.size() << endl;
    cout << "Jets\t\tPt\t\tEta\t\tPhi" << endl;
    for (int i=0; i < PFPUcorrJet.size(); i++){
      cout << "    \t\t" << pTPFJet[i] << "\t\t"  << etaPFJet[i]  << "\t\t" <<  phiPFJet[i] << endl;
    }

    // use PFMET (missing transverse energy)
    TVector3 MET(pxPFMet[0], pyPFMet[0], 0.);

    cout << "PFMET = " << MET.Pt()  << ", phi = " << MET.Phi() << endl;

    cout << "number of photons = " << PFPhoton.size() << endl;
    cout << "iPh1 = " << iPh1 << ", iPh2 = " <<  iPh2 << endl;
    cout << "m(ph1+ph2) = "  << (PFPhoton[iPh1]+PFPhoton[iPh2]).M() << endl;
    cout << "pT(ph1+ph2) = "  << (PFPhoton[iPh1]+PFPhoton[iPh2]).Pt() << endl;
    cout << "eta(ph1+ph2) = "  << (PFPhoton[iPh1]+PFPhoton[iPh2]).Eta() << endl;

    mGG = (PFPhoton[iPh1]+PFPhoton[iPh2]).M();
    pTGG = (PFPhoton[iPh1]+PFPhoton[iPh2]).Pt();
    etaGG = (PFPhoton[iPh1]+PFPhoton[iPh2]).Eta();
    phiGG = (PFPhoton[iPh1]+PFPhoton[iPh2]).Phi();
      
    PFPUcorrJet.push_back(PFPhoton[iPh1]+PFPhoton[iPh2]);
    //    if(iPh1 == -99 && iPh2 == -99) continue;
    Npassed_2Ph += _weight;

    // if only one found, select the highest-pT fake
    // the first is always the good one
    if(iPh1 == -99) {
      iPh1 = iPh2;
      iPh2 = -99;
    }

    // if no good photon was found 
    // use the highest-pT fake
    if(iPh1 == -99) {
      phMax1 = 0.;
      for(int i=0; i<PFPhoton.size(); i++) {
	if(PFPhoton[i].Pt() > phMax1) {
	  phMax2 = PFPhoton[i].Pt();
	  iPh2 = i;
	}
      }
    }

    if(iPh2 == -99) {
      phMax2 = 0.;
      for(int i=0; i<PFPhoton.size(); i++) {
	if(PFPhoton[i].Pt() > phMax2 && i != iPh1 ) {
	  phMax2 = PFPhoton[i].Pt();
	  iPh2 = i;
	}
      }
    }
    
    if(iPh1 == -99) {
      pTPh1  = -99.;
      etaPh1 = -99.;
      phiPh1 = -99.;
      noPhoton1 = -99;
      isEMPh1 = -99;
      isLoosePh1 = -99;
      isTightPh1 = -99;
    } else {
      pTPh1  = double(PFPhoton[iPh1].Pt());
      etaPh1 = double(PFPhoton[iPh1].Eta());
      phiPh1 = double(PFPhoton[iPh1].Phi());
      noPhoton1 = int(notPhoton[iPh1]);
      isEMPh1 = int(isEM[iPh1]);
      isLoosePh1 = int(isLoose[iPh1]);
      isTightPh1 = int(isTight[iPh1]);
    }

    if(iPh2 == -99) {
      // this is gamma+jets
      pTPh2  = -99.;
      etaPh2 = -99.;
      phiPh2 = -99.;
      noPhoton2 = -99;
      isEMPh2 = -99;
      isLoosePh2 = -99;
      isTightPh2 = -99;
    } else {
      pTPh2  = double(PFPhoton[iPh2].Pt());
      etaPh2 = double(PFPhoton[iPh2].Eta());
      phiPh2 = double(PFPhoton[iPh2].Phi());
      noPhoton2 = int(notPhoton[iPh2]);
      isEMPh2 = int(isEM[iPh2]);
      isLoosePh2 = int(isLoose[iPh2]);
      isTightPh2 = int(isTight[iPh2]);
    }


    // dummy values
    passedPF = 0;
    pTPFHem1 = -9999;
    etaPFHem1 = -9999;
    phiPFHem1 = -9999;
    pTPFHem2 = -9999;
    etaPFHem2 = -9999;
    phiPFHem2 = -9999;
    RSQ = -99999.;
    MR = -99999.;

    // hemispheres


    // int flag = 1;
    // //bubble sort the jets by Pt()
    // for (int i=0; (i < PFPUcorrJet.size()) && flag; i++){
    //   TLorentzVector tempvector;
    //   flag = 0;
    //   for (int j=0; j < (PFPUcorrJet.size()-1); j++){
    // 	if (PFPUcorrJet.at(j+1).Pt() > PFPUcorrJet.at(j).Pt()){
    // 	  tempvector = PFPUcorrJet.at(j);
    // 	  PFPUcorrJet.at(j) = PFPUcorrJet.at(j+1);
    // 	  PFPUcorrJet.at(j+1) = tempvector;
    // 	  flag=1;
    // 	}
    //   }
    // }


    
    vector<TLorentzVector> tmpJet = CombineJets(PFPUcorrJet);
    vector<int> temp1 = GetHem1Const(PFPUcorrJet); 
    vector<int> temp2 = GetHem2Const(PFPUcorrJet); 
    

    int k = 0;
    for (int i = 0; i < temp1.size(); i++){
      if (temp1[i]>-1){
	constPFHem1[k] = temp1[i];
	k++;
      }
    }
    k = 0;
    for (int i = 0; i < temp2.size(); i++){
      if (temp2[i]>-1){
	constPFHem2[k] = temp2[i];
	k++;
      }
    }

    if(tmpJet.size() >= 2) {
      Npassed_PFHem += _weight;
      
      TLorentzVector PFHem1 = tmpJet[0];
      TLorentzVector PFHem2 = tmpJet[1];
      
      // compute boost
      double num = PFHem1.P()-PFHem2.P();
      double den = PFHem1.Pz()-PFHem2.Pz();      
      double beta = num/den;
      
      double MT = CalcMTR(PFHem1, PFHem2, MET);
      double variable = -999999.;
      double Rvariable = -999999.;
      variable = CalcGammaMRstar(PFHem1, PFHem2);
      //variable = CalcGammaMRstarE(PFHem1, PFHem2);
      if(variable >0) Rvariable = MT/variable;
      
      // fill the tree

      passedPF = 1;
      pTPFHem1 = PFHem1.Pt();
      etaPFHem1 = PFHem1.Eta();
      phiPFHem1 = PFHem1.Phi();
      pTPFHem2 = PFHem2.Pt();
      etaPFHem2 = PFHem2.Eta();
      phiPFHem2 = PFHem2.Phi();
      RSQ = Rvariable*Rvariable;
      MR = variable;   

      cout << "MR = " << MR << endl;
      cout << "RSQ = " << RSQ << endl;
    }

    // dummy values
    passedCalo = 0;
    pTCaloHem1 = -9999;
    etaCaloHem1 = -9999;
    phiCaloHem1 = -9999;
    pTCaloHem2 = -9999;
    etaCaloHem2 = -9999;
    phiCaloHem2 = -9999;
    CaloRSQ = -99999.;
    CaloMR = -99999.;

    // hemispheres
    CMSHemisphere* CaloHem = new CMSHemisphere(CaloJet);
    CaloHem->CombineMinMass();
    tmpJet = CaloHem->GetHemispheres();
    if(tmpJet.size() >= 2) {
      Npassed_CaloHem += _weight;
    
      TLorentzVector CaloHem1 = tmpJet[0];
      TLorentzVector CaloHem2 = tmpJet[1];
      
      // compute boost
      double num = CaloHem1.P()-CaloHem2.P();
      double den = CaloHem1.Pz()-CaloHem2.Pz();      
      double beta = num/den;
      
      double MT = CalcMTR(CaloHem1, CaloHem2, MET);
      double variable = -999999.;
      double Rvariable = -999999.;
      variable = CalcGammaMRstar(CaloHem1, CaloHem2);
      if(variable >0) Rvariable = MT/variable;
      
      // fill the tree
      passedCalo = 1;
      pTCaloHem1 = CaloHem1.Pt();
      etaCaloHem1 = CaloHem1.Eta();
      phiCaloHem1 = CaloHem1.Phi();
      pTCaloHem2 = CaloHem2.Pt();
      etaCaloHem2 = CaloHem2.Eta();
      phiCaloHem2 = CaloHem2.Phi();
      CaloRSQ = Rvariable*Rvariable;
      CaloMR = variable;    
    }

    run = runNumber;
    evNum = eventNumber;
    bx = bunchCrossing;
    lumi = lumiBlock;
    orbit = orbitNumber;
    W = _weight;

    // Fill the tree per box
    outTree->Fill();

  }

  
  // fill efficiency tree
  TTree* effTree = new TTree("effTree", "effTree");
    
  effTree->Branch("Npassed_In",      &Npassed_In,      "Npassed_In/D");
  effTree->Branch("Npassed_HLT",      &Npassed_HLT,      "Npassed_HLT/D");
  effTree->Branch("Npassed_PV",      &Npassed_PV,      "Npassed_PV/D");
  effTree->Branch("Npassed_2Ph",      &Npassed_2Ph,      "Npassed_2Ph/D");
  effTree->Branch("Npassed_PFHem",      &Npassed_PFHem,      "Npassed_PFHem/D");
  effTree->Branch("Npassed_CaloHem",      &Npassed_CaloHem,      "Npassed_CaloHem/D");
  effTree->Fill();

  TFile *file = new TFile(outFileName.c_str(),"RECREATE");
  effTree->Write();
  outTree->Write();
  file->Close();
}
  
bool RazorHiggsDiPhoton::IsPhotonBarrel(int iPh) {
  bool isBarrel = true;
  // to fix
  double etaSeed = etaPho[iPh];
  //  if(1.479 - fabs(etaSeed) < 0.1) isBarrel = false;
  if(fabs(etaSeed) > 1.5) isBarrel = false;
  return isBarrel;
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

bool RazorHiggsDiPhoton::PhotonIdBarrelisEM(int i) {
  return PhotonIdBarrel(i, "isEM");
}

bool RazorHiggsDiPhoton::PhotonIdBarrelisLoose(int i) {
  return PhotonIdBarrel(i, "isLoose");   
}

bool RazorHiggsDiPhoton::PhotonIdBarrelisTight(int i) {
  return PhotonIdBarrel(i, "isTight");
}

bool RazorHiggsDiPhoton::PhotonIdBarrel(int i, TString Selector) {
  int iSC = superClusterIndexPho[i];
  double pT = energyPho[i]*sin(thetaPho[i]);  
  bool passed = false;

  if (Selector == TString("isTight")){

    /*
    // Photon ID criteria
    if(IsPhotonBarrel(i));
    // Jurassic Isolation
    if(ecalRecHitSumEtConeDR04SC[iSC] <= 4.2+0.006*pT);
    // Tower-based HCAL isolation
    if(hcalTowerSumEtConeDR04SC[iSC] <= 2.2+0.0025*pT);
    //  H/E
    if(hOverESC[iSC] >= 0.05);
    // hollow cone track isolation
    //if(dr04HollowTkSumPtPho[i]<= 3.5+0.001*pT);
    if(dr04HollowTkSumPtPho[i]<= 2.0+0.001*pT);
    // eta width
    float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
    if(sigmaIetaIeta <= 0.013);     
    // track veto is optional... 
    if(!hasPixelSeedPho[i]);
    */

    float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
    
    //if ( IsPhotonBarrel(i) && ecalRecHitSumEtConeDR04SC[iSC] < 4.2+0.006*pT && hcalTowerSumEtConeDR04SC[iSC] < 2.2+0.0025*pT && hOverESC[iSC] < 0.05 && dr04HollowTkSumPtPho[i]< 2.0+0.001*pT  && sigmaIetaIeta < 0.013 && !hasPixelSeedPho[i])
    if ( IsPhotonBarrel(i) )

      passed  = true;

  }else if (Selector == "isLoose"){

    /*
    // Photon ID criteria
    if(IsPhotonBarrel(i));
    // Jurassic Isolation
    if(ecalRecHitSumEtConeDR04SC[iSC] <= 4.2+0.006*pT);
    // Tower-based HCAL isolation
    if(hcalTowerSumEtConeDR04SC[iSC] <= 2.2+0.0025*pT);
    //  H/E
    if(hOverESC[iSC] >= 0.05);
    // hollow cone track isolation
    if(dr04HollowTkSumPtPho[i]<= 3.5+0.001*pT);
    // if(dr04HollowTkSumPtPho[i]<= 2.0+0.001*pT);
    // eta width
    // float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
    // if(sigmaIetaIeta <= 0.013);     
    // track veto is optional... 
    if(!hasPixelSeedPho[i]);
    */

    //if ( IsPhotonBarrel(i) && ecalRecHitSumEtConeDR04SC[iSC] < 4.2+0.006*pT && hcalTowerSumEtConeDR04SC[iSC] < 2.2+0.0025*pT && hOverESC[iSC] < 0.05 && dr04HollowTkSumPtPho[i]< 3.5+0.001*pT  && !hasPixelSeedPho[i])
    if ( IsPhotonBarrel(i) )
 
      passed  = true;

  }else if (Selector == "isEM"){


    /*
  //// Photon ID criteria
  if(IsPhotonBarrel(i));
  //// Jurassic Isolation
  if(ecalRecHitSumEtConeDR04SC[iSC] <= 4.2+0.006*pT);
  //// Tower-based HCAL isolation
  if(hcalTowerSumEtConeDR04SC[iSC] <= 2.2+0.0025*pT);
  ////  H/E
  if(hOverESC[iSC] >= 0.05);
  //// hollow cone track isolation
  //if(dr04HollowTkSumPtPho[i]<= 3.5+0.001*pT);
  //if(dr04HollowTkSumPtPho[i]<= 2.0+0.001*pT);
  //// eta width
  // float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
  // if(sigmaIetaIeta <= 0.013);     
  //// track veto is optional... 
  if(!hasPixelSeedPho[i]);
    */

    //if ( IsPhotonBarrel(i) && ecalRecHitSumEtConeDR04SC[iSC] < 4.2+0.006*pT && hcalTowerSumEtConeDR04SC[iSC] < 2.2+0.0025*pT && hOverESC[iSC] < 0.05 && !hasPixelSeedPho[i])
    if ( IsPhotonBarrel(i) )

      passed  = true;

  }

  return passed;

}


////////////////////////Maurozio////////////////////////
/*
  void RazorHiggsDiPhoton::PhotonIdBarrel(int i, TString Selector) {
  int iSC = superClusterIndexPho[i];
  double pT = energyPho[i]*sin(thetaPho[i]);  
  bool passed = false;
  // Photon ID criteria
  isEM = true;
  isLoose = true;
  isTight = true;
  if(!IsPhotonBarrel(i)) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }

  // Jurassic Isolation
  if(ecalRecHitSumEtConeDR04SC[iSC] >= 4.2+0.006*pT) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }

  // Tower-based HCAL isolation
  if(hcalTowerSumEtConeDR04SC[iSC] >= 2.2+0.0025*pT) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }

  //  H/E
  if(hOverESC[iSC] >= 0.05) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }

  // hollow cone track isolation
  if(dr04HollowTkSumPtPho[i]>= 3.5+0.001*pT) {
  isLoose = false;
  }
  if(dr04HollowTkSumPtPho[i]>= 2.0+0.001*pT) {
  isTight = false;
  }
  // eta width
  float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
  if(sigmaIetaIeta >= 0.013) isTight = false;  
  
  // track veto is optional... 
  if(hasPixelSeedPho[i]) {
  isEM = false;
  isLoose = false;
  isTight = false;
  }
  }
*/

int RazorHiggsDiPhoton::AntiPhotonIdBarrel(int i) {
  int iSC = superClusterIndexPho[i];
  double pT = energyPho[i]*sin(thetaPho[i]);
  // Photon ID criteria
  bool isEM = true;
  if(!IsPhotonBarrel(i)) {
    isEM = false;
  }
  /*
  // Jurassic Isolation
  if(ecalRecHitSumEtConeDR04SC[iSC] >= 4.2+0.006*pT) {
    isEM = false;
  }
  // Tower-based HCAL isolation
  if(hcalTowerSumEtConeDR04SC[iSC] >= 2.2+0.0025*pT) {
    isEM = false;
    }
  */
  //  H/E
  if(hOverESC[iSC] >= 0.05) {
    isEM = false;
  }
  // hollow cone track isolation AND sigmaIetaIeta
  float sigmaIetaIeta = sqrt(covIEtaIEtaSC[iSC]);
  if(dr04HollowTkSumPtPho[i]< 2.0+0.001*pT && sigmaIetaIeta < 0.013) {
    isEM = false;
  }
  // track veto [we want photon-like jets, not electrons]
  if(hasPixelSeedPho[i]) isEM = false;
  return isEM;
}
