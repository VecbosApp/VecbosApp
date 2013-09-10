// std includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
using namespace std;

// ROOT includes
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"

// local includes
#include "Vecbos.hh"
#include "BtagEff.hh"

BtagEff::BtagEff(TTree *tree) : Vecbos(tree) {}

BtagEff::~BtagEff(){}

//void BtagEff::Loop(string OutputFileName){}


void BtagEff::Loop(string OutputFileName)
{
  
   if(fChain == NULL)
      return;

   int Low = 0;
   int High = 3000;

   InitializeAllArrays();

   Long64_t nbytes = 0;
   Long64_t nb = 0;
   Long64_t nentries = fChain->GetEntries();
   cout << "Number of entries = " << nentries << endl;

   int BJetCount = 0;
   int CJetCount = 0;
   int LJetCount = 0;

   int RawBTagCount = 0, Pass1BTagCount = 0, Pass2BTagCount = 0;
   int RawCTagCount = 0, Pass1CTagCount = 0, Pass2CTagCount = 0;
   int RawLTagCount = 0, Pass1LTagCount = 0, Pass2LTagCount = 0;

   TFile OutputFile(OutputFileName.c_str(), "RECREATE");
   TTree OutputTree("OutputTree", "Output tree");

   char ModelName[1024] = "";
   int JetCount = 0;
   double JetPT[50] = {0};
   double JetEta[50] = {0};
   int JetFlavor[50] = {0};
   bool JetTagged[50] = {0};
   bool JetTaggedPass1[50] = {0};
   bool JetTaggedPass2[50] = {0};
   bool JetTaggedPass2Up[50] = {0};
   bool JetTaggedPass2Down[50] = {0};
   int RunNumber = 0;
   long long EventNumber = 0;

   OutputTree.Branch("ModelName", &ModelName, "ModelName/C");
   OutputTree.Branch("JetCount", &JetCount, "JetCount/I");
   OutputTree.Branch("JetPT", JetPT, "JetPT[50]/D");
   OutputTree.Branch("JetEta", JetEta, "JetEta[50]/D");
   OutputTree.Branch("JetFlavor", JetFlavor, "JetFlavor[50]/I");
   OutputTree.Branch("JetTagged", JetTagged, "JetTagged[50]/O");
   OutputTree.Branch("JetTaggedPass1", JetTaggedPass1, "JetTaggedPass1[50]/O");
   OutputTree.Branch("JetTaggedPass2", JetTaggedPass2, "JetTaggedPass2[50]/O");
   OutputTree.Branch("JetTaggedPass2Up", JetTaggedPass2Up, "JetTaggedPass2Up[50]/O");
   OutputTree.Branch("JetTaggedPass2Down", JetTaggedPass2Down, "JetTaggedPass2Down[50]/O");
   OutputTree.Branch("RunNumber", &RunNumber, "RunNumber/I");
   OutputTree.Branch("EventNumber", &EventNumber, "EventNumber/LL");

   for(long iEntry = 0; iEntry < nentries; iEntry++)
   {
      long ientry = LoadTree(iEntry);
      if(ientry < 0)
         break;

      nb = fChain->GetEntry(iEntry);
      nbytes += nb;

      if(iEntry % 500 == 0)
         cout << ">>> Processing event # " << iEntry << endl;

      // PYTHIA STRING
      //string Model = (*commentLHE)[0].substr(8);
      /*
      string Model = (*commentLHE)[1].substr(8);
      cout << "model is " << Model << endl;
      for(int i = 0; i < (int)Model.size(); i++)
      {
         if(Model[i] == '\n')
         {
            Model.erase(Model.begin() + i);
            i = i - 1;
         }
      }
      for(int i = 0; i < (int)Model.size(); i++)
      {
         if(Model[i] == ' ')
         {
	   Model.erase(Model.begin() + i, Model.end());
	   break;
         }
      }

      cout << "model is " << Model << endl;
      strcpy(ModelName, Model.c_str());
      */
      string Model = "TTFullLep";
      strcpy(ModelName, Model.c_str());
      RunNumber = runNumber;
      EventNumber = eventNumber;

      if(ModelIndex.find(Model) == ModelIndex.end())
      {
         if(ModelIndex.size() == 50)   // Max!  Clear everything and restart
            ModelIndex.clear();

         int Current = ModelIndex.size();
         ModelIndex.insert(pair<string, int>(Model, Current));

         double DummyError[500][5] = {{0}};

         ReadFile(B_Eff_PTBins[Current], B_Eff_EtaBins[Current], B_Eff[Current], DummyError,
            "BTagging/8TeVNew/" + Model + "_BEff_Pass1.txt");
         ReadFile(C_Eff_PTBins[Current], C_Eff_EtaBins[Current], C_Eff[Current], DummyError,
            "BTagging/8TeVNew/" + Model + "_CEff_Pass1.txt");
         ReadFile(L_Eff_PTBins[Current], L_Eff_EtaBins[Current], L_Eff[Current], DummyError,
            "BTagging/8TeVNew/" + Model + "_LEff_Pass1.txt");
         
         ReadFile(B_EffFast_PTBins[Current], B_EffFast_EtaBins[Current], B_EffFast[Current], DummyError,
            "BTagging/8TeVNew/" + Model + "_BEff_Raw.txt");
         ReadFile(C_EffFast_PTBins[Current], C_EffFast_EtaBins[Current], C_EffFast[Current], DummyError,
            "BTagging/8TeVNew/" + Model + "_CEff_Raw.txt");
         ReadFile(L_EffFast_PTBins[Current], L_EffFast_EtaBins[Current], L_EffFast[Current], DummyError,
            "BTagging/8TeVNew/" + Model + "_LEff_Raw.txt");
      }
      int CurrentModel = ModelIndex[Model];
      std::cout << "debug 0" << std::endl;
      JetCount = nAK5PFNoPUJet;
      for(int i = 0; i < nAK5PFNoPUJet; i++)
      {
         JetPT[i] = 0;
         JetEta[i] = 0;
         JetFlavor[i] = -1;
         JetTagged[i] = false;
      }

      bool OriginalTag[500] = {false};
      bool Pass1Tag[500] = {false};
      bool Pass2TagCentral[500] = {false};
      bool Pass2TagUp[500] = {false};
      bool Pass2TagDown[500] = {false};
      std::cout << "debug 1" <<std::endl;
      
      for(int i = 0; i < nAK5PFNoPUJet; i++)
      {
         double eta = etaAK5PFNoPUJet[i];
         double phi = phiAK5PFNoPUJet[i];
         double csv = combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[i];
         double px = pxAK5PFNoPUJet[i];
         double py = pyAK5PFNoPUJet[i];
         double pt = sqrt(px * px + py * py);

         JetPT[i] = pt;
         JetEta[i] = eta;
	 std::cout << "debug 2" <<std::endl;
         JetFlavor[i] = DetermineJetFlavor(eta, phi);
	 std::cout << "debug 2.1" <<std::endl;
         if(eta < 2.4 && eta > -2.4 && pt > 40)
         {
            OriginalTag[i] = (csv > 0.679);
	    std::cout << "debug 2.2" <<std::endl;
            Pass1Tag[i] = SwitchJetTag(OriginalTag[i], JetFlavor[i], eta, pt, true, 0, CurrentModel);
	    std::cout << "debug 2.3" <<std::endl;
            Pass2TagCentral[i] = SwitchJetTag(Pass1Tag[i], JetFlavor[i], eta, pt, false, 0, CurrentModel);
	    std::cout << "debug 2.4" <<std::endl;
            Pass2TagUp[i] = SwitchJetTag(Pass1Tag[i], JetFlavor[i], eta, pt, false, 1, CurrentModel);
	    std::cout << "debug 2.5" <<std::endl;
            Pass2TagDown[i] = SwitchJetTag(Pass1Tag[i], JetFlavor[i], eta, pt, false, -1, CurrentModel);
	    std::cout << "debug 2.6" <<std::endl;
            JetTagged[i] = OriginalTag[i];
	    std::cout << "debug 2.7" <<std::endl;
            JetTaggedPass1[i] = Pass1Tag[i];
	    std::cout << "debug 2.8" <<std::endl;
            JetTaggedPass2[i] = Pass2TagCentral[i];
	    std::cout << "debug 2.9" <<std::endl;
            JetTaggedPass2Up[i] = Pass2TagUp[i];
	    std::cout << "debug 2.10" <<std::endl;
            JetTaggedPass2Down[i] = Pass2TagDown[i];
	    std::cout << "debug 2.11" <<std::endl;
         }
	 
	 std::cout << "debug 3" <<std::endl;

         if(JetFlavor[i] == JetFlavor_BJet)   BJetCount = BJetCount + 1;
         if(JetFlavor[i] == JetFlavor_CJet)   CJetCount = CJetCount + 1;
         if(JetFlavor[i] == JetFlavor_LJet)   LJetCount = LJetCount + 1;

         if(OriginalTag[i] == true)
         {
            if(JetFlavor[i] == JetFlavor_BJet)   RawBTagCount = RawBTagCount + 1;
            if(JetFlavor[i] == JetFlavor_CJet)   RawCTagCount = RawCTagCount + 1;
            if(JetFlavor[i] == JetFlavor_LJet)   RawLTagCount = RawLTagCount + 1;
         }

	 std::cout << "debug 4" <<std::endl;
         if(Pass1Tag[i] == true)
         {
            if(JetFlavor[i] == JetFlavor_BJet)   Pass1BTagCount = Pass1BTagCount + 1;
            if(JetFlavor[i] == JetFlavor_CJet)   Pass1CTagCount = Pass1CTagCount + 1;
            if(JetFlavor[i] == JetFlavor_LJet)   Pass1LTagCount = Pass1LTagCount + 1;
         }
	 std::cout << "debug 5" <<std::endl;
         if(Pass2TagCentral[i] == true)
         {
            if(JetFlavor[i] == JetFlavor_BJet)   Pass2BTagCount = Pass2BTagCount + 1;
            if(JetFlavor[i] == JetFlavor_CJet)   Pass2CTagCount = Pass2CTagCount + 1;
            if(JetFlavor[i] == JetFlavor_LJet)   Pass2LTagCount = Pass2LTagCount + 1;
         }
      }

      std::cout << "debug 7" <<std::endl;
      OutputTree.Fill();
      std::cout << "debug 8" <<std::endl;
   }

   OutputTree.Write();
   OutputFile.Close();

   cout << "RAW:" << endl;
   cout << "   B = " << RawBTagCount << "/" << BJetCount << endl;
   cout << "   C = " << RawCTagCount << "/" << CJetCount << endl;
   cout << "   L = " << RawLTagCount << "/" << LJetCount << endl;
   cout << endl;
   
   cout << "Pass1:" << endl;
   cout << "   B = " << Pass1BTagCount << "/" << BJetCount << endl;
   cout << "   C = " << Pass1CTagCount << "/" << CJetCount << endl;
   cout << "   L = " << Pass1LTagCount << "/" << LJetCount << endl;
   cout << endl;
   
   cout << "Pass2 (central):" << endl;
   cout << "   B = " << Pass2BTagCount << "/" << BJetCount << endl;
   cout << "   C = " << Pass2CTagCount << "/" << CJetCount << endl;
   cout << "   L = " << Pass2LTagCount << "/" << LJetCount << endl;
   cout << endl;
  
}

int BtagEff::DetermineJetFlavor(double JetEta, double JetPhi)
{
   // Not sure if this is correct.  We should probably match b/c-hadrons,
   //    not bare quarks.  But this is what's in the MultiTop analysis.

  
   double dR = 0.5;

   bool HasB = false;
   bool HasC = false;

   for(int i = 0; i < nMc; i++)
   {
      if(GetDR(JetEta, JetPhi, etaMc[i], phiMc[i]) < dR)
      {
         if(idMc[i] == 5 || idMc[i] == -5)
            HasB = true;
         if(idMc[i] == 4 || idMc[i] == -4)
            HasC = true;
      }
   }

   if(HasB == true)
      return JetFlavor_BJet;
   if(HasC == true)
      return JetFlavor_CJet;
   return JetFlavor_LJet;

  
}

double BtagEff::GetDR(double Eta1, double Phi1, double Eta2, double Phi2)
{
  
   double DEta = Eta1 - Eta2;
   double DPhi = Phi1 - Phi2;

   while(DPhi < 0)
      DPhi = DPhi + 2 * PI;
   while(DPhi >= 2 * PI)
      DPhi = DPhi - 2 * PI;

   DPhi = min(DPhi, 2 * PI - DPhi);

   return sqrt(DEta * DEta + DPhi * DPhi);
  
}

bool BtagEff::SwitchTag(bool CurrentTag, double ScaleFactor, double Efficiency)
{
  
   static TRandom3 Random(time(NULL));   // lazy coder here~

   if(ScaleFactor == 1)
      return CurrentTag;
   if(Efficiency == 0)
      return false;

   if(ScaleFactor > 1)
   {
      if(CurrentTag == false)   // Chance of upgrade!
      {
         float MistagPercentage = (ScaleFactor - 1) * Efficiency / (1 - Efficiency);
         // float MistagPercentage = (ScaleFactor - 1) * Efficiency;

         if(Random.Uniform(1) < MistagPercentage)
            return true;
      }
   }
   else
   {
      if(CurrentTag == true)   // Maybe downgrade!
      {
         if(Random.Uniform(1) > ScaleFactor)
            return false;
         // if(Random.Uniform(1) > (1 - ScaleFactor) * Efficiency)
         //    return false;
      }
   }

   return CurrentTag;

  
}

bool BtagEff::SwitchJetTag(bool CurrentTag, int JetFlavor,
   double JetEta, double JetPT, bool DoFastSim, double Shift, int Model)
{
  
  std::cout << "SwitchJetTag 0" << std::endl; 
  if(JetEta < 0)
    JetEta = -JetEta;
  if(JetEta > 2.4)   // anything greater than 2.4 is not tagged
    return false;
  
  double ScaleFactor = 0, Efficiency = 0;
  std::cout << "SwitchJetTag 1"<< std::endl;
  if(DoFastSim == false){
    int PTBin = -1, EtaBin = -1;
    
    if(JetFlavor == JetFlavor_BJet)   PTBin = FindBin(JetPT, B_SF_PTBins);
    if(JetFlavor == JetFlavor_CJet)   PTBin = FindBin(JetPT, C_SF_PTBins);
    if(JetFlavor == JetFlavor_LJet)   PTBin = FindBin(JetPT, L_SF_PTBins);
    
    if(JetFlavor == JetFlavor_BJet)   EtaBin = FindBin(JetEta, B_SF_EtaBins);
    if(JetFlavor == JetFlavor_CJet)   EtaBin = FindBin(JetEta, C_SF_EtaBins);
    if(JetFlavor == JetFlavor_LJet)   EtaBin = FindBin(JetEta, L_SF_EtaBins);
    std::cout << "SwitchJetTag 2"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)
      ScaleFactor = B_SF[PTBin][EtaBin] + Shift * B_SF_Error[PTBin][EtaBin];
    if(JetFlavor == JetFlavor_CJet)
      ScaleFactor = C_SF[PTBin][EtaBin] + Shift * C_SF_Error[PTBin][EtaBin];
    if(JetFlavor == JetFlavor_LJet)
      ScaleFactor = L_SF[PTBin][EtaBin] + Shift * L_SF_Error[PTBin][EtaBin];
    
    if(JetFlavor == JetFlavor_BJet)   PTBin = FindBin(JetPT, B_Eff_PTBins[Model]);
    if(JetFlavor == JetFlavor_CJet)   PTBin = FindBin(JetPT, C_Eff_PTBins[Model]);
    if(JetFlavor == JetFlavor_LJet)   PTBin = FindBin(JetPT, L_Eff_PTBins[Model]);
    std::cout << "SwitchJetTag 3"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)   EtaBin = FindBin(JetEta, B_Eff_EtaBins[Model]);
    if(JetFlavor == JetFlavor_CJet)   EtaBin = FindBin(JetEta, C_Eff_EtaBins[Model]);
    if(JetFlavor == JetFlavor_LJet)   EtaBin = FindBin(JetEta, L_Eff_EtaBins[Model]);
    std::cout << "SwitchJetTag 4"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)   Efficiency = B_Eff[Model][PTBin][EtaBin];
    if(JetFlavor == JetFlavor_CJet)   Efficiency = C_Eff[Model][PTBin][EtaBin];
    if(JetFlavor == JetFlavor_LJet)   Efficiency = L_Eff[Model][PTBin][EtaBin];
  }else{//Do fastsim
    int PTBin = -1, EtaBin = -1;
    std::cout << "SwitchJetTag 5"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)   PTBin = FindBin(JetPT, B_SFFast_PTBins);
    if(JetFlavor == JetFlavor_CJet)   PTBin = FindBin(JetPT, C_SFFast_PTBins);
    if(JetFlavor == JetFlavor_LJet)   PTBin = FindBin(JetPT, L_SFFast_PTBins);
    std::cout << "SwitchJetTag 6"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)   EtaBin = FindBin(JetEta, B_SFFast_EtaBins);
    if(JetFlavor == JetFlavor_CJet)   EtaBin = FindBin(JetEta, C_SFFast_EtaBins);
    if(JetFlavor == JetFlavor_LJet)   EtaBin = FindBin(JetEta, L_SFFast_EtaBins);
    std::cout << "SwitchJetTag 7"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)
      ScaleFactor = B_SFFast[PTBin][EtaBin] + Shift * B_SFFast_Error[PTBin][EtaBin];
    if(JetFlavor == JetFlavor_CJet)
      ScaleFactor = C_SFFast[PTBin][EtaBin] + Shift * C_SFFast_Error[PTBin][EtaBin];
    if(JetFlavor == JetFlavor_LJet)
      ScaleFactor = L_SFFast[PTBin][EtaBin] + Shift * L_SFFast_Error[PTBin][EtaBin];
    std::cout << "SwitchJetTag 8"<< std::endl;
    
    // if(JetFlavor == JetFlavor_BJet)
    //{
    //   PTBin = FindBin(JetPT, B_SF_PTBins);
    //   cout << JetPT << " " << PTBin << " " << B_SFFast_PTBins[PTBin] << endl;
    //   cout << JetEta << " " << EtaBin << " " << B_SFFast_EtaBins[EtaBin] << endl;
    //  cout << B_SFFast[PTBin][EtaBin] << endl;
    //}
    std::cout << "Model: " << Model << std::endl; 
    std::cout << "SwitchJetTag 9"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)   PTBin = FindBin(JetPT, B_EffFast_PTBins[Model]);
    if(JetFlavor == JetFlavor_CJet)   PTBin = FindBin(JetPT, C_EffFast_PTBins[Model]);
    if(JetFlavor == JetFlavor_LJet)   PTBin = FindBin(JetPT, L_EffFast_PTBins[Model]);
    std::cout << "SwitchJetTag 10"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)   EtaBin = FindBin(JetEta, B_EffFast_EtaBins[Model]);
    if(JetFlavor == JetFlavor_CJet)   EtaBin = FindBin(JetEta, C_EffFast_EtaBins[Model]);
    if(JetFlavor == JetFlavor_LJet)   EtaBin = FindBin(JetEta, L_EffFast_EtaBins[Model]);
    std::cout << "SwitchJetTag 11"<< std::endl;
    if(JetFlavor == JetFlavor_BJet)   Efficiency = B_EffFast[Model][PTBin][EtaBin];
    if(JetFlavor == JetFlavor_CJet)   Efficiency = C_EffFast[Model][PTBin][EtaBin];
    if(JetFlavor == JetFlavor_LJet)   Efficiency = L_EffFast[Model][PTBin][EtaBin];
  }
  std::cout << "SwitchJetTag 12"<< std::endl;
  return SwitchTag(CurrentTag, ScaleFactor, Efficiency);
  
  
}

int BtagEff::FindBin(double Value, const vector<double> &Bins)
{
  
  std::cout << "find bin 0" << std::endl;
  if(Value < Bins[0])
    return 0;
  std::cout << "find bin 1" << std::endl;
  for(int i = 1; i < (int)Bins.size(); i++)
    {
      if(Value < Bins[i])
	return i - 1;
    }
  std::cout << "find bin 2" << std::endl;
  return Bins.size() - 1;
  
}

void BtagEff::ReadFile(vector<double> &PT, vector<double> &Eta,
   double Array[500][5], double Error[500][5], string FileName, bool DoubleError)
{
  
   cout << "Start reading in file \"" << FileName << "\"" << endl;

   for(int i = 0; i < 500; i++)
   {
      for(int j = 0; j < 5; j++)
      {
         Array[i][j] = 0;
         Error[i][j] = 0;
      }
   }

   map<double, string> PTStrings;
   map<double, string> EtaStrings;
   map<pair<string, string>, double> Values;
   map<pair<string, string>, double> Errors;

   ifstream in(FileName.c_str());

   while(in)
   {
      char line[10240] = "";
      in.getline(line, 10239, '\n');
      if(line[0] == '\0')
         continue;

      stringstream str(line);

      string pt = "", eta = "";
      double value = -1, error = 0;

      str >> pt >> eta >> value >> error;

      if(value < 0 || pt == "" || eta == "")
         continue;

      PTStrings.insert(pair<double, string>(atof(pt.c_str()), pt));
      EtaStrings.insert(pair<double, string>(atof(eta.c_str()), eta));

      pair<string, string> index(pt, eta);

      Values.insert(pair<pair<string, string>, double>(index, value));
      Errors.insert(pair<pair<string, string>, double>(index, error));
   }

   PT.clear();
   for(map<double, string>::iterator iter = PTStrings.begin();
      iter != PTStrings.end(); iter++)
      PT.push_back(iter->first);
   
   Eta.clear();
   for(map<double, string>::iterator iter = EtaStrings.begin();
      iter != EtaStrings.end(); iter++)
      Eta.push_back(iter->first);
   
   int indexpt = 0, indexeta = 0;
   for(map<double, string>::iterator iter1 = PTStrings.begin();
      iter1 != PTStrings.end(); iter1++)
   {
      indexeta = 0;

      for(map<double, string>::iterator iter2 = EtaStrings.begin();
         iter2 != EtaStrings.end(); iter2++)
      {
         pair<string, string> Index(iter1->second, iter2->second);

         Array[indexpt][indexeta] = Values[Index];

         if(DoubleError == true)
            Error[indexpt][indexeta] = Errors[Index] * 2;
         else
            Error[indexpt][indexeta] = Errors[Index];

         indexeta = indexeta + 1;
      }

      indexpt = indexpt + 1;
   }
   
   in.close();

   if(PT.size() == 0 || Eta.size() == 0)
      cerr << "Warning!  Nothing is read from file " << FileName << endl;

  
}

void BtagEff::InitializeAllArrays()
{
  
   double DummyError[500][5];
   
   ReadFile(B_SF_PTBins, B_SF_EtaBins, B_SF, B_SF_Error,
      "BTagging/8TeVNew/BEff_SF_CSVM.txt", false);
   ReadFile(C_SF_PTBins, C_SF_EtaBins, C_SF, C_SF_Error,
      "BTagging/8TeVNew/BEff_SF_CSVM.txt", true);
   ReadFile(L_SF_PTBins, L_SF_EtaBins, L_SF, L_SF_Error,
      "BTagging/8TeVNew/Mistag_SF_CSVM.txt", false);

   
   //ReadFile(B_Eff_PTBins, B_Eff_EtaBins, B_Eff, DummyError,
   //   "BTagging/8TeVNew/T1bbbb_BEff_Pass1.txt");
   //ReadFile(C_Eff_PTBins, C_Eff_EtaBins, C_Eff, DummyError,
   //   "BTagging/8TeVNew/T1bbbb_CEff_Pass1.txt");
   //ReadFile(L_Eff_PTBins, L_Eff_EtaBins, L_Eff, DummyError,
   //   "BTagging/8TeVNew/T1bbbb_LEff_Pass1.txt");
   
   
   ReadFile(B_SFFast_PTBins, B_SFFast_EtaBins, B_SFFast, B_SFFast_Error,
      "BTagging/8TeVNew/BEff_SF_FastSim.txt", false);
   ReadFile(C_SFFast_PTBins, C_SFFast_EtaBins, C_SFFast, C_SFFast_Error,
      "BTagging/8TeVNew/CEff_SF_FastSim.txt", false);
   ReadFile(L_SFFast_PTBins, L_SFFast_EtaBins, L_SFFast, L_SFFast_Error,
      "BTagging/8TeVNew/Mistag_SF_FastSim.txt", false);
   
   
   //ReadFile(B_EffFast_PTBins, B_EffFast_EtaBins, B_EffFast, DummyError,
   //   "BTagging/8TeVNew/T1bbbb_BEff_Raw.txt");
   //ReadFile(C_EffFast_PTBins, C_EffFast_EtaBins, C_EffFast, DummyError,
   //   "BTagging/8TeVNew/T1bbbb_CEff_Raw.txt");
   //ReadFile(L_EffFast_PTBins, L_EffFast_EtaBins, L_EffFast, DummyError,
   //   "BTagging/8TeVNew/T1bbbb_LEff_Raw.txt");
     
  
}

