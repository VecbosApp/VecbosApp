#ifndef ControlSampleEvents_HH
#define ControlSampleEvents_HH

//C++ INCLUDES
#include <cmath>
#include "assert.h"
#include <Rtypes.h>
//ROOT INCLUDES
#include "TFile.h"
#include "TTree.h"
#include "TError.h"
#include "TLorentzVector.h"

  class ControlSampleEvents {

    public:

      /// bit map
      /// DON'T CHANGE ORDER

      //*******************************************
      //=== Process IDs  ====
      //*******************************************
      enum BkgProcessId { kData = 0,
			  kQCD = 1,
			  kWJets = 2,
			  kZJets = 3,
			  kTTJets = 4,
			  kSingleT = 5,
			  kVV = 6,
			  kHiggs = 7,
			  kSUSY = 99,
			  kUnknown = 999
      };

      /// variables
      Float_t                 weight;
      UInt_t                  run;
      UInt_t                  lumi;
      UInt_t                  event;
      UInt_t                  processID;
      UInt_t                  NPU_0;
      UInt_t                  NPU_Minus1;
      UInt_t                  NPU_Plus1;
      Bool_t                  Trigger;
      Bool_t                  NoiseFilter;
      Float_t                 lep1Iso;
      Float_t                 lep1Iso_2;
      //float                 lep1KinkMuon;
      Float_t                 lep2Iso;
      Float_t                 lep2Iso_2;
      //float                 lep2KinkMuon;
      TLorentzVector          genlep1;
      TLorentzVector          genlep2;
      Int_t                   genlep1Type;
      Int_t                   genlep2Type;
      TLorentzVector          genPhoton1;
      TLorentzVector          genPhoton2;
      Bool_t                  foundGenPhoton1;
      Bool_t                  foundGenPhoton2;
      TLorentzVector          lep1;
      TLorentzVector          lep2;
      Int_t                   lep1Type;
      Int_t                   lep2Type;
      Int_t                   lep1MatchedGenLepIndex;
      Int_t                   lep2MatchedGenLepIndex;
      Bool_t                  lep1PassVeto;
      Bool_t                  lep1PassLoose;
      Bool_t                  lep1PassTight;
      Bool_t                  lep2PassVeto;
      Bool_t                  lep2PassLoose;
      Bool_t                  lep2PassTight;
      TLorentzVector          bjet1;
      TLorentzVector          bjet2;
      Bool_t                  bjet1PassLoose;
      Bool_t                  bjet1PassMedium;
      Bool_t                  bjet1PassTight;
      Bool_t                  bjet2PassLoose;
      Bool_t                  bjet2PassMedium;
      Bool_t                  bjet2PassTight;
      TLorentzVector          jet1;
      TLorentzVector          jet2;      
      Bool_t                  bad_jet;
      Bool_t                  jet1PassCSVLoose;
      Bool_t                  jet1PassCSVMedium;
      Bool_t                  jet1PassCSVTight;
      Bool_t                  jet2PassCSVLoose;
      Bool_t                  jet2PassCSVMedium;
      Bool_t                  jet2PassCSVTight;
      Float_t                 MR;
      Float_t                 Rsq;
      Float_t                 MR_LeadLeptonsAsMET;
      Float_t                 Rsq_LeadLeptonsAsMET;
      Float_t                 MR_LeptonsAsMET;
      Float_t                 Rsq_LeptonsAsMET;
      Float_t                 MR_lep;
      Float_t                 MRT_lep;
      Float_t                 Rsq_lep;
      Float_t                 MR_res;
      Float_t                 MRT_res;
      Float_t                 Rsq_res;
      Float_t                 MET;
      Float_t                 MET_LeadLeptonsAsMET;
      Float_t                 MET_LeptonsAsMET;
      Float_t                 minDPhi;
      Float_t                 minDPhiN;
      Float_t                 minDPhiFolded;
      Float_t                 minDPhiNFolded;
      Float_t                 dPhiHemHem;
      Float_t                 dPhiHemHem_LeptonsAsMET;
      Float_t                 dPhiHemHem_LeadLeptonsAsMET;
      UInt_t                  NJets30;
      UInt_t                  NJets40;
      UInt_t                  NJets80;
      UInt_t                  NBJetsLoose;
      UInt_t                  NBJetsMedium;
      UInt_t                  NBJetsTight;
      Float_t                 HT;
      Float_t                 lep1MT;
      UInt_t                  nPhotonsAbove20;
      TLorentzVector          photon1;
      Float_t                 photon1SCEta;
      Int_t                   photon1MatchedGenIndex;
      TLorentzVector          photon2;
      Float_t                 photon2SCEta;
      Int_t                   photon2MatchedGenIndex;

    public:
      /// this is the main element
      TTree *tree_;
      TFile *f_;
  
      /// hold the names of variables to facilitate things (filled during Init)
      std::vector<std::string> variables_;

      /// default constructor  
      ControlSampleEvents()  {
	genlep1Ptr  = &genlep1;
	genlep2Ptr  = &genlep2;
        genPhoton1Ptr = &genPhoton1;
        genPhoton2Ptr = &genPhoton2;
	lep1Ptr     = &lep1;
	lep2Ptr     = &lep2;
        photon1Ptr  = &photon1;
        photon2Ptr  = &photon2;
	bjet1Ptr    = &bjet1;
	bjet2Ptr    = &bjet2;       
	jet1Ptr     = &jet1;
	jet2Ptr     = &jet2;       
      };

      /// default destructor
      ~ControlSampleEvents(){ 
        if (f_) f_->Close();  
      };
    
      /// initialize varibles and fill list of available variables
      void InitVariables() {
	weight               = 0.0;
	run                  = 0.0;
	lumi                 = 0.0;
	event                = 0.0;
	processID            = ControlSampleEvents::kUnknown;
	NPU_0                = 0.0;
	NPU_Minus1           = 0.0;
	NPU_Plus1            = 0.0;
        Trigger              = false;
        NoiseFilter          = false;
	lep1Iso              = -999.0;
	lep1Iso_2            = -999.0;
	//lep1KinkMuon         = -999.0;
	lep2Iso              = -999.0;
        lep2Iso_2            = -999.0;
        //lep2KinkMuon         = -999.0;
	genlep1              = TLorentzVector();
	genlep2              = TLorentzVector();
	genlep1Type          = 0.0;
	genlep2Type          = 0.0;
        genPhoton1           = TLorentzVector();
        genPhoton2           = TLorentzVector();
        foundGenPhoton1      = false;
        foundGenPhoton2      = false;
	lep1                 = TLorentzVector();
	lep2                 = TLorentzVector();
	lep1Type             = 0.0;
	lep2Type             = 0.0;
	lep1MatchedGenLepIndex = -1;
	lep2MatchedGenLepIndex = -1;
	lep1PassVeto         = 0.0;
	lep1PassLoose        = 0.0;
	lep1PassTight        = 0.0;
	lep2PassVeto         = 0.0;
	lep2PassLoose        = 0.0;
	lep2PassTight        = 0.0;
	bjet1                = TLorentzVector();
	bjet2                = TLorentzVector();
	bjet1PassLoose       = 0.0;
	bjet1PassMedium      = 0.0;
	bjet1PassTight       = 0.0;
	bjet2PassLoose       = 0.0;
	bjet2PassMedium      = 0.0;
	bjet2PassTight       = 0.0;
	jet1                 = TLorentzVector();
	jet2                 = TLorentzVector();
	bad_jet              = false;
	jet1PassCSVLoose     = 0.0;
	jet1PassCSVMedium    = 0.0;
	jet1PassCSVTight     = 0.0;
	jet2PassCSVLoose     = 0.0;
	jet2PassCSVMedium    = 0.0;
	jet2PassCSVTight     = 0.0;
	MR                   = 0.0;
	Rsq                  = 0.0;
	MR_LeadLeptonsAsMET        = 0.0;
	Rsq_LeadLeptonsAsMET       = 0.0;
	MR_LeptonsAsMET         = 0.0;
	Rsq_LeptonsAsMET        = 0.0;
	MR_lep                  = 0.0;
	MRT_lep                 = 0.0;
	Rsq_lep                 = 0.0;
	MR_res                  = 0.0;
	MRT_res                 = 0.0;
	Rsq_res                 = 0.0;
	MET                  = 0.0;
	MET_LeadLeptonsAsMET       = 0.0;
	MET_LeptonsAsMET        = 0.0;
	minDPhi              = 0.0;
	minDPhiN             = 0.0;
        dPhiHemHem           = -1.0;
        dPhiHemHem_LeptonsAsMET          = -1.0;
        dPhiHemHem_LeadLeptonsAsMET           = -1.0;
	NJets30              = 0.0;
	NJets40              = 0.0;
	NJets80              = 0.0;
	NBJetsLoose          = 0.0;
	NBJetsMedium         = 0.0;
	NBJetsTight          = 0.0;
	HT                   = 0.0;      
	lep1MT               = 0.0;      
        nPhotonsAbove20      = 0;
        photon1              = TLorentzVector();
        photon1SCEta         = -999;
        photon1MatchedGenIndex = -1;
        photon2              = TLorentzVector();
        photon2SCEta         = -999;
        photon2MatchedGenIndex = -1;
      }
    
      /// load a ControlSampleEvents
      void LoadTree(const char* file){
        f_ = TFile::Open(file);
        assert(f_);
        tree_ = dynamic_cast<TTree*>(f_->Get("ControlSampleEvent"));
	InitTree();
        assert(tree_);
      }
    
      /// create a ControlSampleEvents
      void CreateTree(){
        tree_ = new TTree("ControlSampleEvent","ControlSampleEvent");
        f_ = 0;

        //book the branches
	tree_->Branch("weight",&weight,"weight/F");
	tree_->Branch("run",&run,"run/i");
	tree_->Branch("lumi",&lumi,"lumi/i");
	tree_->Branch("event",&event,"event/i");
	tree_->Branch("processID",&processID,"processID/i");
	tree_->Branch("NPU_0",&NPU_0,"NPU_0/i");
	tree_->Branch("NPU_Minus1",&NPU_Minus1,"NPU_Minus1/i");
	tree_->Branch("NPU_Plus1",&NPU_Plus1,"NPU_Plus1/i");
	tree_->Branch("Trigger", &Trigger, "Trigger/O");
	tree_->Branch("NoiseFilter", &NoiseFilter, "NoiseFilter/O");
	
	//tree_->Branch("lep1KinkMuon", &lep1KinkMuon, "lep1KinkMuon/F");
	tree_->Branch("lep1Iso", &lep1Iso, "lep1Iso/F");
	tree_->Branch("lep1Iso_2", &lep1Iso_2, "lep1Iso_2/F");
	//tree_->Branch("lep2KinkMuon", &lep2KinkMuon, "lep2KinkMuon/F");
        tree_->Branch("lep2Iso", &lep2Iso, "lep2Iso/F");
        tree_->Branch("lep2Iso_2", &lep2Iso_2, "lep2Iso_2/F");
	
	tree_->Branch("genlep1Type",&genlep1Type,"genlep1Type/I");
	tree_->Branch("genlep2Type",&genlep2Type,"genlep2Type/I");
        tree_->Branch("foundGenPhoton1", &foundGenPhoton1, "foundGenPhoton1/O");
        tree_->Branch("foundGenPhoton2", &foundGenPhoton2, "foundGenPhoton2/O");
	tree_->Branch("lep1Type",&lep1Type,"lep1Type/I");
	tree_->Branch("lep2Type",&lep2Type,"lep2Type/I");
	tree_->Branch("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex,"lep1MatchedGenLepIndex/I");
	tree_->Branch("lep2MatchedGenLepIndex",&lep2MatchedGenLepIndex,"lep2MatchedGenLepIndex/I");
	tree_->Branch("lep1PassVeto",&lep1PassVeto,"lep1PassVeto/O");
	tree_->Branch("lep1PassLoose",&lep1PassLoose,"lep1PassLoose/O");
	tree_->Branch("lep1PassTight",&lep1PassTight,"lep1PassTight/O");
	tree_->Branch("lep2PassVeto",&lep2PassVeto,"lep2PassVeto/O");
	tree_->Branch("lep2PassLoose",&lep2PassLoose,"lep2PassLoose/O");
	tree_->Branch("lep2PassTight",&lep2PassTight,"lep2PassTight/O");
        tree_->Branch("photon1SCEta", &photon1SCEta, "photon1SCEta/F");
        tree_->Branch("photon2SCEta", &photon2SCEta, "photon2SCEta/F");
        tree_->Branch("photon1MatchedGenIndex", &photon1MatchedGenIndex, "photon1MatchedGenIndex/I");
        tree_->Branch("photon2MatchedGenIndex", &photon2MatchedGenIndex, "photon2MatchedGenIndex/I");
	tree_->Branch("bjet1PassLoose",&bjet1PassLoose,"bjet1PassLoose/O");
	tree_->Branch("bjet1PassMedium",&bjet1PassMedium,"bjet1PassMedium/O");
	tree_->Branch("bjet1PassTight",&bjet1PassTight,"bjet1PassTight/O");
	tree_->Branch("bjet2PassLoose",&bjet2PassLoose,"bjet2PassLoose/O");
	tree_->Branch("bjet2PassMedium",&bjet2PassMedium,"bjet2PassMedium/O");
	tree_->Branch("bjet2PassTight",&bjet2PassTight,"bjet2PassTight/O");
	tree_->Branch("bad_jet",&bad_jet,"bad_jet/O");
	tree_->Branch("jet1PassCSVLoose",&jet1PassCSVLoose,"jet1PassCSVLoose/O");
	tree_->Branch("jet1PassCSVMedium",&jet1PassCSVMedium,"jet1PassCSVMedium/O");
	tree_->Branch("jet1PassCSVTight",&jet1PassCSVTight,"jet1PassCSVTight/O");
	tree_->Branch("jet2PassCSVLoose",&jet2PassCSVLoose,"jet2PassCSVLoose/O");
	tree_->Branch("jet2PassCSVMedium",&jet2PassCSVMedium,"jet2PassCSVMedium/O");
	tree_->Branch("jet2PassCSVTight",&jet2PassCSVTight,"jet2PassCSVTight/O");
	tree_->Branch("MR",&MR,"MR/F");
	tree_->Branch("Rsq",&Rsq,"Rsq/F");
	tree_->Branch("MR_LeadLeptonsAsMET",&MR_LeadLeptonsAsMET,"MR_LeadLeptonsAsMET/F");
	tree_->Branch("Rsq_LeadLeptonsAsMET",&Rsq_LeadLeptonsAsMET,"Rsq_LeadLeptonsAsMET/F");
	tree_->Branch("MR_LeptonsAsMET",&MR_LeptonsAsMET,"MR_LeptonsAsMET/F");
	tree_->Branch("Rsq_LeptonsAsMET",&Rsq_LeptonsAsMET,"Rsq_LeptonsAsMET/F");
	tree_->Branch("MR_lep",&MR_lep,"MR_lep/F");
	tree_->Branch("MRT_lep",&MRT_lep,"MRT_lep/F");
        tree_->Branch("Rsq_lep",&Rsq_lep,"Rsq_lep/F");
	tree_->Branch("MR_res",&MR_res,"MR_res/F");
	tree_->Branch("MRT_res",&MRT_res,"MRT_res/F");
	tree_->Branch("Rsq_res",&Rsq_res,"Rsq_res/F");
	tree_->Branch("MET",&MET,"MET/F");
	tree_->Branch("MET_LeadLeptonsAsMET",&MET_LeadLeptonsAsMET,"MET_LeadLeptonsAsMET/F");
	tree_->Branch("MET_LeptonsAsMET",&MET_LeptonsAsMET,"MET_LeptonsAsMET/F");
	tree_->Branch("minDPhi",&minDPhi,"minDPhi/F");
	tree_->Branch("minDPhiN",&minDPhiN,"minDPhiN/F");
	tree_->Branch("dPhiHemHem",&dPhiHemHem,"dPhiHemHem/F");
	tree_->Branch("dPhiHemHem_LeptonsAsMET",&dPhiHemHem_LeptonsAsMET,"dPhiHemHem_LeptonsAsMET/F");
	tree_->Branch("dPhiHemHem_LeadLeptonsAsMET",&dPhiHemHem_LeadLeptonsAsMET,"dPhiHemHem_LeadLeptonsAsMET/F");
	tree_->Branch("NJets30",&NJets30,"NJets30/i");
	tree_->Branch("NJets40",&NJets40,"NJets40/i");
	tree_->Branch("NJets80",&NJets80,"NJets80/i");
	tree_->Branch("NBJetsLoose",&NBJetsLoose,"NBJetsLoose/i");
	tree_->Branch("NBJetsMedium",&NBJetsMedium,"NBJetsMedium/i");
	tree_->Branch("NBJetsTight",&NBJetsTight,"NBJetsTight/i");
	tree_->Branch("HT",&HT,"HT/F");
	tree_->Branch("lep1MT",&lep1MT,"lep1MT/F");
	tree_->Branch("genlep1", "TLorentzVector", &genlep1Ptr);
	tree_->Branch("genlep2", "TLorentzVector", &genlep2Ptr);
        tree_->Branch("genPhoton1", "TLorentzVector", &genPhoton1Ptr);
        tree_->Branch("genPhoton2", "TLorentzVector", &genPhoton2Ptr);
	tree_->Branch("lep1",    "TLorentzVector", &lep1Ptr);
	tree_->Branch("lep2",    "TLorentzVector", &lep2Ptr);
        tree_->Branch("photon1", "TLorentzVector", &photon1Ptr);
        tree_->Branch("photon2", "TLorentzVector", &photon2Ptr);
	tree_->Branch("bjet1",   "TLorentzVector", &bjet1Ptr);
	tree_->Branch("bjet2",   "TLorentzVector", &bjet2Ptr);
	tree_->Branch("jet1",    "TLorentzVector", &jet1Ptr);
	tree_->Branch("jet2",    "TLorentzVector", &jet2Ptr);
      } 

      // initialze a ControlSampleEvents
      void InitTree(){
        assert(tree_);
        // don't forget to set pointers to zero before you set address
        // or you will fully appreciate that "ROOT sucks" :)
        InitVariables();
        //Set branch address
        Int_t currentState = gErrorIgnoreLevel;

	tree_->SetBranchAddress("weight",&weight);
	tree_->SetBranchAddress("run",&run);
	tree_->SetBranchAddress("lumi",&lumi);
	tree_->SetBranchAddress("event",&event);
	tree_->SetBranchAddress("processID",&processID);
	tree_->SetBranchAddress("NPU_0",&NPU_0);
	tree_->SetBranchAddress("NPU_Minus1",&NPU_Minus1);
	tree_->SetBranchAddress("NPU_Plus1",&NPU_Plus1);
	tree_->SetBranchAddress("Trigger", &Trigger);
	tree_->SetBranchAddress("NoiseFilter", &NoiseFilter);
	
	//tree_->SetBranchAddress("lep1KinkMuon", &lep1KinkMuon);
	tree_->SetBranchAddress("lep1Iso", &lep1Iso);
	tree_->SetBranchAddress("lep1Iso_2", &lep1Iso_2);
	
	//tree_->SetBranchAddress("lep2KinkMuon", &lep2KinkMuon);
        tree_->SetBranchAddress("lep2Iso", &lep2Iso);
        tree_->SetBranchAddress("lep2Iso_2", &lep2Iso_2);

	tree_->SetBranchAddress("genlep1Type",&genlep1Type);
	tree_->SetBranchAddress("genlep2Type",&genlep2Type);
        tree_->SetBranchAddress("foundGenPhoton1", &foundGenPhoton1);
        tree_->SetBranchAddress("foundGenPhoton2", &foundGenPhoton2);
	tree_->SetBranchAddress("lep1Type",&lep1Type);
	tree_->SetBranchAddress("lep2Type",&lep2Type);
	tree_->SetBranchAddress("lep1MatchedGenLepIndex",&lep1MatchedGenLepIndex);
	tree_->SetBranchAddress("lep2MatchedGenLepIndex",&lep2MatchedGenLepIndex);
	tree_->SetBranchAddress("lep1PassVeto",&lep1PassVeto);
	tree_->SetBranchAddress("lep1PassLoose",&lep1PassLoose);
	tree_->SetBranchAddress("lep1PassTight",&lep1PassTight);
	tree_->SetBranchAddress("lep2PassVeto",&lep2PassVeto);
	tree_->SetBranchAddress("lep2PassLoose",&lep2PassLoose);
	tree_->SetBranchAddress("lep2PassTight",&lep2PassTight);
        tree_->SetBranchAddress("photon1SCEta", &photon1SCEta);
        tree_->SetBranchAddress("photon2SCEta", &photon2SCEta);
        tree_->SetBranchAddress("photon1MatchedGenIndex", &photon1MatchedGenIndex);
        tree_->SetBranchAddress("photon2MatchedGenIndex", &photon2MatchedGenIndex);
	tree_->SetBranchAddress("bjet1PassLoose",&bjet1PassLoose);
	tree_->SetBranchAddress("bjet1PassMedium",&bjet1PassMedium);
	tree_->SetBranchAddress("bjet1PassTight",&bjet1PassTight);
	tree_->SetBranchAddress("bjet2PassLoose",&bjet2PassLoose);
	tree_->SetBranchAddress("bjet2PassMedium",&bjet2PassMedium);
	tree_->SetBranchAddress("bjet2PassTight",&bjet2PassTight);
	tree_->SetBranchAddress("bad_jet", &bad_jet);
	tree_->SetBranchAddress("jet1PassCSVLoose",&jet1PassCSVLoose);
	tree_->SetBranchAddress("jet1PassCSVMedium",&jet1PassCSVMedium);
	tree_->SetBranchAddress("jet1PassCSVTight",&jet1PassCSVTight);
	tree_->SetBranchAddress("jet2PassCSVLoose",&jet2PassCSVLoose);
	tree_->SetBranchAddress("jet2PassCSVMedium",&jet2PassCSVMedium);
	tree_->SetBranchAddress("jet2PassCSVTight",&jet2PassCSVTight);
	tree_->SetBranchAddress("MR",&MR);
	tree_->SetBranchAddress("Rsq",&Rsq);
	tree_->SetBranchAddress("MR_LeadLeptonsAsMET",&MR_LeadLeptonsAsMET);
	tree_->SetBranchAddress("Rsq_LeadLeptonsAsMET",&Rsq_LeadLeptonsAsMET);
	tree_->SetBranchAddress("MR_LeptonsAsMET",&MR_LeptonsAsMET);
	tree_->SetBranchAddress("Rsq_LeptonsAsMET",&Rsq_LeptonsAsMET);
	tree_->SetBranchAddress("MR_lep",&MR_lep);
	tree_->SetBranchAddress("MRT_lep",&MRT_lep);
        tree_->SetBranchAddress("Rsq_lep",&Rsq_lep);
	tree_->SetBranchAddress("MR_res",&MR_res);
	tree_->SetBranchAddress("MRT_res",&MRT_res);
        tree_->SetBranchAddress("Rsq_res",&Rsq_res);
	tree_->SetBranchAddress("MET",&MET);
	tree_->SetBranchAddress("MET_LeadLeptonsAsMET",&MET_LeadLeptonsAsMET);
	tree_->SetBranchAddress("MET_LeptonsAsMET",&MET_LeptonsAsMET);
	tree_->SetBranchAddress("minDPhi",&minDPhi);
	tree_->SetBranchAddress("minDPhiN",&minDPhiN);
	tree_->SetBranchAddress("dPhiHemHem",&dPhiHemHem);
	tree_->SetBranchAddress("dPhiHemHem_LeptonsAsMET",&dPhiHemHem_LeptonsAsMET);
	tree_->SetBranchAddress("dPhiHemHem_LeadLeptonsAsMET",&dPhiHemHem_LeadLeptonsAsMET);
	tree_->SetBranchAddress("NJets30",&NJets30);
	tree_->SetBranchAddress("NJets40",&NJets40);
	tree_->SetBranchAddress("NJets80",&NJets80);
	tree_->SetBranchAddress("NBJetsLoose",&NBJetsLoose);
	tree_->SetBranchAddress("NBJetsMedium",&NBJetsMedium);
	tree_->SetBranchAddress("NBJetsTight",&NBJetsTight);
	tree_->SetBranchAddress("HT",&HT);
	tree_->SetBranchAddress("lep1MT",&lep1MT);	
	tree_->SetBranchAddress("genlep1",&genlep1Ptr);
	tree_->SetBranchAddress("genlep2",&genlep2Ptr);
        tree_->SetBranchAddress("genPhoton1", &genPhoton1Ptr);
        tree_->SetBranchAddress("genPhoton2", &genPhoton2Ptr);
	tree_->SetBranchAddress("lep1",&lep1Ptr);
	tree_->SetBranchAddress("lep2",&lep2Ptr);
        tree_->SetBranchAddress("photon1", &photon1Ptr);
        tree_->SetBranchAddress("photon2", &photon2Ptr);
	tree_->SetBranchAddress("bjet1",&bjet1Ptr);
	tree_->SetBranchAddress("bjet2",&bjet2Ptr);
	tree_->SetBranchAddress("jet1",&jet1Ptr);
	tree_->SetBranchAddress("jet2",&jet2Ptr);
	
        gErrorIgnoreLevel = currentState;
      }

    private:
      TLorentzVector* genlep1Ptr;
      TLorentzVector* genlep2Ptr;
      TLorentzVector* genPhoton1Ptr;
      TLorentzVector* genPhoton2Ptr;
      TLorentzVector* lep1Ptr;
      TLorentzVector* lep2Ptr;
      TLorentzVector* photon1Ptr;
      TLorentzVector* photon2Ptr;
      TLorentzVector* bjet1Ptr;
      TLorentzVector* bjet2Ptr;
      TLorentzVector* jet1Ptr;
      TLorentzVector* jet2Ptr;
      
  }; 


#endif

