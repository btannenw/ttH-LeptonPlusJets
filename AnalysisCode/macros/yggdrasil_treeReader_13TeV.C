// Description: Quickly move from Yggdrasil Trees or Express Trees to BDTs



#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TMath.h"
#include "TRandom3.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "TGraphAsymmErrors.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "Math/Interpolator.h"

#include <time.h>


#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"


#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "ttH-LeptonPlusJets/AnalysisCode/interface/YggdrasilEventVars.h"

#endif




typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;


void fillCSVhistos(TFile *fileHF, TFile *fileLF);
double get_csv_wgt( vvdouble jets, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );

// CSV reweighting
TH1D* h_csv_wgt_hf[9][6];
TH1D* c_csv_wgt_hf[9][6];
TH1D* h_csv_wgt_lf[9][4][3];





void yggdrasil_treeReader_13TeV(  int SampleType=2,int SplitType =2, int maxNentries=-1, int Njobs=1, int jobN=1 ) {

int insample=1;


clock_t t;
t = clock();

  TFile* f_CSVwgt_HF = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/csv_rwt_hf_IT_FlatSF_2015_07_27.root").c_str());
  TFile* f_CSVwgt_LF = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/csv_rwt_lf_IT_FlatSF_2015_07_27.root").c_str());

  fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);
  
  std::cout << "   ===> load the root files! " << std::endl;

std::string sampleName = "TTJets";
if(SampleType == 2){
	if(SplitType == 1) sampleName = "ttb";
	if(SplitType == 2) sampleName = "ttbb";
	if(SplitType == 3) sampleName = "tt2b";
	if(SplitType == 4) sampleName = "ttcc";
	if(SplitType == 5) sampleName = "ttlf";
 }
if(SampleType == 1){
	if(SplitType == 9) sampleName = "ttH";
	if(SplitType == 1) sampleName = "ttH_hbb_err"; //Sometimes hbb is not seperated from ttH samples. 
	if(SplitType == 2) sampleName = "ttH_hww";
	if(SplitType == 3) sampleName = "ttH_htt";
	if(SplitType == 4) sampleName = "ttH_hgluglu";
	if(SplitType == 5) sampleName = "ttH_hcc";
	if(SplitType == 6) sampleName = "ttH_hzz";
	if(SplitType == 7) sampleName = "ttH_hzg";
	if(SplitType == 8) sampleName = "ttH_hgg";
	if(SplitType == 10) sampleName = "ttH_hss";
	if(SplitType == 11) sampleName = "ttH_hmumu";
} 
if(SampleType == 3) sampleName ="ttH_hbb";
if(SampleType == 4) sampleName = "ttZ";
if(SampleType == 5) sampleName = "ttW";

 
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();

  std::string treefilename ="/eos/uscms/store/user/sflowers/treeMaker/Spring15_Sep3/TT_TuneCUETP8M1_13TeV-powheg-pythia8/Spring15DR-PU20bx25_MCRUN2_74_V9-v1_yggdrasilTree_v1/150911_144050/0000/*.root";


  //  std::string __my_eos__("/eos/uscms/store/user/satoshi/");
  std::string __my_eos__("/eos/uscms/store/user/sflowers/treeMaker/Spring15_Sep3/");


  treefilename.assign( __my_eos__ + "TT_TuneCUETP8M1_13TeV-powheg-pythia8/Validation_v1/150915_132026/0000" + "/*.root" );
  if(SampleType == 1)treefilename =  __my_eos__ + "ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8/Validation_v1/150911_174759/0000" + "/*.root";
  if(SampleType == 3)treefilename =  __my_eos__ + "" + "/*.root";
  if(SampleType == 4)treefilename =  __my_eos__ + "" + "/*.root";

  std::string s_end = "_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = ".root";
  
  cout<<sampleName<<" *** "<<endl;

  std::string histofilename = "/uscms/home/satoshi/work/tth/osu/data/" + sampleName + s_end; ### Change This For your environment.

  std::cout << "  treefilename  = " << treefilename.c_str() << std::endl;
  std::cout << "  histofilename = " << histofilename.c_str() << std::endl;

  TChain *chain = new TChain("ttHTreeMaker/worldTree");
  chain->Add(treefilename.c_str());
  
  ////////////////////////////////////
  
  
  std::vector<std::string> cat_labels;
  cat_labels.push_back("4j2t"); //0
  cat_labels.push_back("5j2t"); //1
  cat_labels.push_back("6j2t"); //2
  cat_labels.push_back("4j3t"); //3
  cat_labels.push_back("5j3t"); //4
  cat_labels.push_back("6j3t"); //5
  cat_labels.push_back("4j4t"); //6
  cat_labels.push_back("5j4t"); //7
  cat_labels.push_back("6j4t"); //8
  
  int NumCat = int(cat_labels.size());

  TMVA::Reader *reader_Spring15_KITreoptimized_BDT[NumCat];
  TMVA::Reader *reader_Splitting_1_1[NumCat];
  TMVA::Reader *reader_Splitting_1_2[NumCat];
  
  for( int c=0; c<NumCat; c++ ){
	  reader_Spring15_KITreoptimized_BDT[c] = new TMVA::Reader("!Color:Silent");
	  reader_Splitting_1_1[c] = new TMVA::Reader("!Color:Silent");
	  reader_Splitting_1_2[c] = new TMVA::Reader("!Color:Silent");
  }
	  
  ///BDT Variables
 
  Float_t numJets_float; // due to a bug the reader only takes floats
  Float_t numTags_float;

  Float_t tight_lepton_pt;
  Float_t tight_lepton_phi;
  Float_t first_jet_pt;
  Float_t min_dr_tagged_jets;
  Float_t avg_dr_tagged_jets;
  Float_t aplanarity;
  Float_t sphericity;
  Float_t avg_btag_disc_non_btags;
  Float_t MET;
  Float_t second_jet_pt;
  Float_t dr_between_lep_and_closest_jet;
  Float_t h0;
  Float_t avg_btag_disc_btags;
  Float_t dev_from_avg_disc_btags;
  Float_t third_jet_pt;
  Float_t fourth_jet_pt;
  Float_t avg_tagged_dijet_mass;
  Float_t h1;
  Float_t h2;
  Float_t lowest_btag;
  Float_t all_sum_pt_with_met;
  Float_t all_sum_pt_incl_met;
  Float_t all_sum_pt;
  Float_t avg_untagged_dijet_mass;
  Float_t closest_tagged_dijet_mass;
  Float_t h3;
  Float_t h4;
  Float_t first_highest_btag;
  Float_t second_highest_btag;
  Float_t third_highest_btag;
  Float_t fourth_highest_btag;
  Float_t fifth_highest_CSV;
  Float_t sixth_highest_CSV;
  
  Float_t CSV_Average;
  Float_t Deta_JetsAverage;
  Float_t Deta_JetsAverageIC;
  
  //Float_t dijet_mass_of_everything;

  Float_t MHT;
  Float_t MHT_phi;

  Float_t invariant_mass_of_everything;

  Float_t MT_mht;

  Float_t HT;

  Float_t M3;
  Float_t M3_1tag;
  Float_t Mlb;

  Float_t best_higgs_mass;
  Float_t minChi2;
  Float_t dRbb;

  Float_t tagged_dijet_mass_closest_to_125;

  
  Float_t bhmv2 = 0.;
  Float_t maxeta_jet_jet = 0.;			 
  Float_t maxeta_jet_tag = 0.;			 
  Float_t maxeta_tag_tag = 0.;			 
  Float_t abs_dEta_leptop_bb = 0.;		 
  Float_t abs_dEta_hadtop_bb = 0.;		 
  Float_t dEta_fn = 0.;					 
  Float_t abs_bb_eta = 0.;				 
  Float_t angle_tops_bb = 0.;			 
  Float_t median_bb_mass = 0.;			 
  Float_t pt_all_jets_over_E_all_jets = 0.;

  
  //Weights Locations
  
  TString BDTWgtDIR_Spring15      = string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/V3_KIT_weights/";
  TString BDTWgtDIR_Splitting_1_1 = string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/Train_1_1/";
  TString BDTWgtDIR_Splitting_1_2 = string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/Train_1_2/";
  
  
  /////////////

//KITBDTS


  //4j3t
  
  //reader_Spring15_KITreoptimized_BDT[3]->AddVariable("CSV_Average", &CSV_Average);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("BDTOhio_v2_input_h1", &h1);
  //reader_Spring15_KITreoptimized_BDT[3]->AddVariable("HT", &HT);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("BDTOhio_v2_input_avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("BDTOhio_v2_input_sphericity", &sphericity);
  //reader_Spring15_KITreoptimized_BDT[3]->AddVariable("avg_btag_disc_btags", &avg_btag_disc_btags);
  //reader_Spring15_KITreoptimized_BDT[3]->AddVariable("dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("BDTOhio_v2_input_third_highest_btag", &third_highest_btag);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("BDTOhio_v2_input_HT",&HT);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("BDTOhio_v2_input_dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("BDTOhio_v2_input_M3", &M3);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("BDTOhio_v2_input_min_dr_tagged_jets", &min_dr_tagged_jets);
  reader_Spring15_KITreoptimized_BDT[3]->AddVariable("Evt_CSV_Average", &CSV_Average);
  reader_Spring15_KITreoptimized_BDT[3]->BookMVA("BDT method", BDTWgtDIR_Spring15 + "weights_Final_43_KITV3.xml");
  
    
  //4j4t
  
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_h1", &h1);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_sphericity", &sphericity);
  //reader_Spring15_KITreoptimized_BDT[6]->AddVariable("dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_avg_btag_disc_btags", &avg_btag_disc_btags);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_h2", &h2);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_closest_tagged_dijet_mass", &closest_tagged_dijet_mass);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_second_highest_btag", &second_highest_btag);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_fourth_highest_btag", &fourth_highest_btag);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_maxeta_jet_jet", &maxeta_jet_jet);
  reader_Spring15_KITreoptimized_BDT[6]->AddVariable("BDTOhio_v2_input_pt_all_jets_over_E_all_jets", &pt_all_jets_over_E_all_jets);
  reader_Spring15_KITreoptimized_BDT[6]->BookMVA("BDT method", BDTWgtDIR_Spring15 + "weights_Final_44_KITV3.xml");
  
  
  //5j3t
  
  reader_Spring15_KITreoptimized_BDT[4]->AddVariable("BDTOhio_v2_input_h1", &h1);
  //reader_Spring15_KITreoptimized_BDT[4]->AddVariable("HT", &HT);
  reader_Spring15_KITreoptimized_BDT[4]->AddVariable("BDTOhio_v2_input_avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Spring15_KITreoptimized_BDT[4]->AddVariable("BDTOhio_v2_input_sphericity", &sphericity);
  //reader_Spring15_KITreoptimized_BDT[4]->AddVariable("dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  //reader_Spring15_KITreoptimized_BDT[4]->AddVariable("h3", &h3);
  reader_Spring15_KITreoptimized_BDT[4]->AddVariable("BDTOhio_v2_input_third_highest_btag", &third_highest_btag);
  reader_Spring15_KITreoptimized_BDT[4]->AddVariable("BDTOhio_v2_input_h3", &h3);
  reader_Spring15_KITreoptimized_BDT[4]->AddVariable("BDTOhio_v2_input_HT", &HT);
  reader_Spring15_KITreoptimized_BDT[4]->AddVariable("BDTOhio_v2_input_dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  reader_Spring15_KITreoptimized_BDT[4]->AddVariable("BDTOhio_v2_input_fourth_highest_btag", &fourth_highest_btag);
  reader_Spring15_KITreoptimized_BDT[4]->BookMVA("BDT method", BDTWgtDIR_Spring15 + "weights_Final_53_KITV3.xml");
  
  
  //5j4t
    reader_Spring15_KITreoptimized_BDT[7]->AddVariable("BDTOhio_v2_input_avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Spring15_KITreoptimized_BDT[7]->AddVariable("BDTOhio_v2_input_HT", &HT);
 reader_Spring15_KITreoptimized_BDT[7]->AddVariable("BDTOhio_v2_input_M3", &M3);
   reader_Spring15_KITreoptimized_BDT[7]->AddVariable("BDTOhio_v2_input_fifth_highest_CSV", &fifth_highest_CSV);
//  reader_Spring15_KITreoptimized_BDT[7]->AddVariable("avg_btag_disc_btags", &avg_btag_disc_btags);
  reader_Spring15_KITreoptimized_BDT[7]->AddVariable("BDTOhio_v2_input_fourth_highest_btag", &fourth_highest_btag);

  reader_Spring15_KITreoptimized_BDT[7]->AddVariable("Evt_Deta_JetsAverage", &Deta_JetsAverage);
 reader_Spring15_KITreoptimized_BDT[7]->AddVariable("BDTOhio_v2_input_avg_btag_disc_btags", &avg_btag_disc_btags);
  reader_Spring15_KITreoptimized_BDT[7]->BookMVA("BDT method", BDTWgtDIR_Spring15 + "weights_Final_54_KITV3.xml");
  

  
  //6j2t
  
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_h1", &h1);
  //reader_Spring15_KITreoptimized_BDT[2]->AddVariable("HT", &HT);
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_sphericity", &sphericity);
  //reader_Spring15_KITreoptimized_BDT[2]->AddVariable("h3", &h3);
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_third_highest_btag", &third_highest_btag);
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_h3", &h3);
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_HT", &HT);
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_Mlb", &Mlb);
  //reader_Spring15_KITreoptimized_BDT[2]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_fifth_highest_CSV", &fifth_highest_CSV);
  reader_Spring15_KITreoptimized_BDT[2]->AddVariable("BDTOhio_v2_input_fourth_highest_btag", &fourth_highest_btag);
  //reader_Spring15_KITreoptimized_BDT[2]->AddVariable("Mlb", &Mlb);
  reader_Spring15_KITreoptimized_BDT[2]->BookMVA("BDT method", BDTWgtDIR_Spring15 + "weights_Final_62_KITV3.xml");
  
  
  //6j3t
  
  reader_Spring15_KITreoptimized_BDT[5]->AddVariable("BDTOhio_v2_input_h1", &h1);
  //reader_Spring15_KITreoptimized_BDT[5]->AddVariable("HT", &HT);
  reader_Spring15_KITreoptimized_BDT[5]->AddVariable("BDTOhio_v2_input_avg_dr_tagged_jets", &avg_dr_tagged_jets);
  
  //reader_Spring15_KITreoptimized_BDT[5]->AddVariable("h3", &h3);
  reader_Spring15_KITreoptimized_BDT[5]->AddVariable("BDTOhio_v2_input_third_highest_btag", &third_highest_btag);
  reader_Spring15_KITreoptimized_BDT[5]->AddVariable("BDTOhio_v2_input_HT", &HT);
  reader_Spring15_KITreoptimized_BDT[5]->AddVariable("BDTOhio_v2_input_pt_all_jets_over_E_all_jets", &pt_all_jets_over_E_all_jets);
  reader_Spring15_KITreoptimized_BDT[5]->AddVariable("BDTOhio_v2_input_fifth_highest_CSV", &fifth_highest_CSV);
  
  reader_Spring15_KITreoptimized_BDT[5]->AddVariable("BDTOhio_v2_input_fourth_highest_btag", &fourth_highest_btag);
  reader_Spring15_KITreoptimized_BDT[5]->AddVariable("BDTOhio_v2_input_min_dr_tagged_jets", &min_dr_tagged_jets);
  //reader_Spring15_KITreoptimized_BDT[5]->AddVariable("pt_all_jets_over_E_all_jets", &pt_all_jets_over_E_all_jets);
  reader_Spring15_KITreoptimized_BDT[5]->BookMVA("BDT method", BDTWgtDIR_Spring15 + "weights_Final_63_KITV3.xml");
  
  
  //6j4t
  
 // reader_Spring15_KITreoptimized_BDT[8]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Spring15_KITreoptimized_BDT[8]->AddVariable("BDTOhio_v2_input_h2", &h2);
  reader_Spring15_KITreoptimized_BDT[8]->AddVariable("BDTOhio_v2_input_avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Spring15_KITreoptimized_BDT[8]->AddVariable("BDTOhio_v2_input_third_highest_btag", &third_highest_btag);
  reader_Spring15_KITreoptimized_BDT[8]->AddVariable("BDTOhio_v2_input_M3", &M3);
   reader_Spring15_KITreoptimized_BDT[8]->AddVariable("BDTOhio_v2_input_fifth_highest_CSV", &fifth_highest_CSV);
  reader_Spring15_KITreoptimized_BDT[8]->AddVariable("BDTOhio_v2_input_fourth_highest_btag", &fourth_highest_btag);
 
  //reader_Spring15_KITreoptimized_BDT[8]->AddVariable("M3", &M3);
  reader_Spring15_KITreoptimized_BDT[8]->AddVariable("BDTOhio_v2_input_best_higgs_mass", &best_higgs_mass);
  reader_Spring15_KITreoptimized_BDT[8]->AddVariable("BDTOhio_v2_input_dEta_fn", &dEta_fn);
  reader_Spring15_KITreoptimized_BDT[8]->BookMVA("BDT method", BDTWgtDIR_Spring15 + "weights_Final_64_KITV3.xml");
  


  //BDTs for checking even/odd splitting. 
  

  //4j3t
  reader_Splitting_1_1[3]->AddVariable("CSV_Average", &CSV_Average);
  reader_Splitting_1_1[3]->AddVariable("h1", &h1);
  reader_Splitting_1_1[3]->AddVariable("HT",&HT);
  reader_Splitting_1_1[3]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_1[3]->AddVariable("sphericity", &sphericity);
  reader_Splitting_1_1[3]->AddVariable("dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  reader_Splitting_1_1[3]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_1[3]->AddVariable("M3", &M3);
  reader_Splitting_1_1[3]->AddVariable("min_dr_tagged_jets", &min_dr_tagged_jets);
  reader_Splitting_1_1[3]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_1 + "433/TMVAClassification_BDT.weights.xml");
  //4j4t
  reader_Splitting_1_1[6]->AddVariable("h1", &h1);
  reader_Splitting_1_1[6]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_1[6]->AddVariable("sphericity", &sphericity);
  reader_Splitting_1_1[6]->AddVariable("avg_btag_disc_btags", &avg_btag_disc_btags);
  reader_Splitting_1_1[6]->AddVariable("h2", &h2);
  reader_Splitting_1_1[6]->AddVariable("closest_tagged_dijet_mass", &closest_tagged_dijet_mass);
  reader_Splitting_1_1[6]->AddVariable("second_highest_btag", &second_highest_btag);
  reader_Splitting_1_1[6]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_1[6]->AddVariable("maxeta_jet_jet", &maxeta_jet_jet);
  reader_Splitting_1_1[6]->AddVariable("pt_all_jets_over_E_all_jets", &pt_all_jets_over_E_all_jets);
  reader_Splitting_1_1[6]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_1 + "443/TMVAClassification_BDT.weights.xml");
  //5j3t
  reader_Splitting_1_1[4]->AddVariable("h1", &h1);
  reader_Splitting_1_1[4]->AddVariable("HT", &HT);
  reader_Splitting_1_1[4]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_1[4]->AddVariable("sphericity", &sphericity);
  reader_Splitting_1_1[4]->AddVariable("dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  reader_Splitting_1_1[4]->AddVariable("h3", &h3);
  reader_Splitting_1_1[4]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_1[4]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_1[4]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_1 + "533/TMVAClassification_BDT.weights.xml");
  //5j4t
  reader_Splitting_1_1[7]->AddVariable("HT", &HT);
  reader_Splitting_1_1[7]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_1[7]->AddVariable("avg_btag_disc_btags", &avg_btag_disc_btags);
  reader_Splitting_1_1[7]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_1[7]->AddVariable("fifth_highest_CSV", &fifth_highest_CSV);
  reader_Splitting_1_1[7]->AddVariable("Deta_JetsAverage", &Deta_JetsAverage);
  reader_Splitting_1_1[7]->AddVariable("M3", &M3);
  reader_Splitting_1_1[7]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_1 + "543/TMVAClassification_BDT.weights.xml");
  //6j2t
  reader_Splitting_1_1[2]->AddVariable("h1", &h1);
  reader_Splitting_1_1[2]->AddVariable("HT", &HT);
  reader_Splitting_1_1[2]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_1[2]->AddVariable("sphericity", &sphericity);
  reader_Splitting_1_1[2]->AddVariable("h3", &h3);
  reader_Splitting_1_1[2]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_1[2]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_1[2]->AddVariable("fifth_highest_CSV", &fifth_highest_CSV);
  reader_Splitting_1_1[2]->AddVariable("Mlb", &Mlb);
  reader_Splitting_1_1[2]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_1 + "623/TMVAClassification_BDT.weights.xml");
  //6j3t
  reader_Splitting_1_1[5]->AddVariable("h1", &h1);
  reader_Splitting_1_1[5]->AddVariable("HT", &HT);
  reader_Splitting_1_1[5]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_1[5]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_1[5]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_1[5]->AddVariable("fifth_highest_CSV", &fifth_highest_CSV);
  reader_Splitting_1_1[5]->AddVariable("min_dr_tagged_jets", &min_dr_tagged_jets);
  reader_Splitting_1_1[5]->AddVariable("pt_all_jets_over_E_all_jets", &pt_all_jets_over_E_all_jets);
  reader_Splitting_1_1[5]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_1 + "633/TMVAClassification_BDT.weights.xml");
  //6j4t
  reader_Splitting_1_1[8]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_1[8]->AddVariable("h2", &h2);
  reader_Splitting_1_1[8]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_1[8]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_1[8]->AddVariable("fifth_highest_CSV", &fifth_highest_CSV);
  reader_Splitting_1_1[8]->AddVariable("M3", &M3);
  reader_Splitting_1_1[8]->AddVariable("best_higgs_mass", &best_higgs_mass);
  reader_Splitting_1_1[8]->AddVariable("dEta_fn", &dEta_fn);
  reader_Splitting_1_1[8]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_1 + "643/TMVAClassification_BDT.weights.xml");
  
  //4j3t
  reader_Splitting_1_2[3]->AddVariable("CSV_Average", &CSV_Average);
  reader_Splitting_1_2[3]->AddVariable("h1", &h1);
  reader_Splitting_1_2[3]->AddVariable("HT",&HT);
  reader_Splitting_1_2[3]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_2[3]->AddVariable("sphericity", &sphericity);
  reader_Splitting_1_2[3]->AddVariable("dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  reader_Splitting_1_2[3]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_2[3]->AddVariable("M3", &M3);
  reader_Splitting_1_2[3]->AddVariable("min_dr_tagged_jets", &min_dr_tagged_jets);
  reader_Splitting_1_2[3]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_2 + "433/TMVAClassification_BDT.weights.xml");
  //4j4t
  reader_Splitting_1_2[6]->AddVariable("h1", &h1);
  reader_Splitting_1_2[6]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_2[6]->AddVariable("sphericity", &sphericity);
  reader_Splitting_1_2[6]->AddVariable("avg_btag_disc_btags", &avg_btag_disc_btags);
  reader_Splitting_1_2[6]->AddVariable("h2", &h2);
  reader_Splitting_1_2[6]->AddVariable("closest_tagged_dijet_mass", &closest_tagged_dijet_mass);
  reader_Splitting_1_2[6]->AddVariable("second_highest_btag", &second_highest_btag);
  reader_Splitting_1_2[6]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_2[6]->AddVariable("maxeta_jet_jet", &maxeta_jet_jet);
  reader_Splitting_1_2[6]->AddVariable("pt_all_jets_over_E_all_jets", &pt_all_jets_over_E_all_jets);
  reader_Splitting_1_2[6]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_2 + "443/TMVAClassification_BDT.weights.xml");
  //5j3t
  reader_Splitting_1_2[4]->AddVariable("h1", &h1);
  reader_Splitting_1_2[4]->AddVariable("HT", &HT);
  reader_Splitting_1_2[4]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_2[4]->AddVariable("sphericity", &sphericity);
  reader_Splitting_1_2[4]->AddVariable("dev_from_avg_disc_btags", &dev_from_avg_disc_btags);
  reader_Splitting_1_2[4]->AddVariable("h3", &h3);
  reader_Splitting_1_2[4]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_2[4]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_2[4]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_2 + "533/TMVAClassification_BDT.weights.xml");
  //5j4t
  reader_Splitting_1_2[7]->AddVariable("HT", &HT);
  reader_Splitting_1_2[7]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_2[7]->AddVariable("avg_btag_disc_btags", &avg_btag_disc_btags);
  reader_Splitting_1_2[7]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_2[7]->AddVariable("fifth_highest_CSV", &fifth_highest_CSV);
  reader_Splitting_1_2[7]->AddVariable("Deta_JetsAverage", &Deta_JetsAverage);
  reader_Splitting_1_2[7]->AddVariable("M3", &M3);
  reader_Splitting_1_2[7]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_2 + "543/TMVAClassification_BDT.weights.xml");
  //6j2t
  reader_Splitting_1_2[2]->AddVariable("h1", &h1);
  reader_Splitting_1_2[2]->AddVariable("HT", &HT);
  reader_Splitting_1_2[2]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_2[2]->AddVariable("sphericity", &sphericity);
  reader_Splitting_1_2[2]->AddVariable("h3", &h3);
  reader_Splitting_1_2[2]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_2[2]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_2[2]->AddVariable("fifth_highest_CSV", &fifth_highest_CSV);
  reader_Splitting_1_2[2]->AddVariable("Mlb", &Mlb);
  reader_Splitting_1_2[2]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_2 + "623/TMVAClassification_BDT.weights.xml");
  //6j3t
  reader_Splitting_1_2[5]->AddVariable("h1", &h1);
  reader_Splitting_1_2[5]->AddVariable("HT", &HT);
  reader_Splitting_1_2[5]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_2[5]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_2[5]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_2[5]->AddVariable("fifth_highest_CSV", &fifth_highest_CSV);
  reader_Splitting_1_2[5]->AddVariable("min_dr_tagged_jets", &min_dr_tagged_jets);
  reader_Splitting_1_2[5]->AddVariable("pt_all_jets_over_E_all_jets", &pt_all_jets_over_E_all_jets);
  reader_Splitting_1_2[5]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_2 + "633/TMVAClassification_BDT.weights.xml");
  //6j4t
  reader_Splitting_1_2[8]->AddVariable("avg_dr_tagged_jets", &avg_dr_tagged_jets);
  reader_Splitting_1_2[8]->AddVariable("h2", &h2);
  reader_Splitting_1_2[8]->AddVariable("third_highest_btag", &third_highest_btag);
  reader_Splitting_1_2[8]->AddVariable("fourth_highest_btag", &fourth_highest_btag);
  reader_Splitting_1_2[8]->AddVariable("fifth_highest_CSV", &fifth_highest_CSV);
  reader_Splitting_1_2[8]->AddVariable("M3", &M3);
  reader_Splitting_1_2[8]->AddVariable("best_higgs_mass", &best_higgs_mass);
  reader_Splitting_1_2[8]->AddVariable("dEta_fn", &dEta_fn);
  reader_Splitting_1_2[8]->BookMVA("BDT method", BDTWgtDIR_Splitting_1_2 + "643/TMVAClassification_BDT.weights.xml");
  



 
  
 
  /////////////
 
  
	yggdrasilEventVars *eve=0;
    chain->SetBranchAddress("eve.", &eve );
	
    TFile histofile(histofilename.c_str(),"recreate");
    histofile.cd();
	
	///////Systematics/////
	
    std::vector<std::string> sys_cat_labels;
    sys_cat_labels.push_back("");               //0
    sys_cat_labels.push_back("_lepIdSFUp");     //1
    sys_cat_labels.push_back("_lepIdSFDown");   //2
    sys_cat_labels.push_back("_PUUp");          //3
    sys_cat_labels.push_back("_PUDown");        //4
    sys_cat_labels.push_back("_JERUp");         //5
    sys_cat_labels.push_back("_JERDown");       //6
    sys_cat_labels.push_back("_JESUp");         //7
    sys_cat_labels.push_back("_JESDown");       //8
    sys_cat_labels.push_back("_CSVLFUp");       //9
    sys_cat_labels.push_back("_CSVLFDown");     //10
    sys_cat_labels.push_back("_CSVHFUp");       //11
    sys_cat_labels.push_back("_CSVHFDown");     //12
    sys_cat_labels.push_back("_CSVHFStats1Up");     //13
    sys_cat_labels.push_back("_CSVHFStats1Down");   //14
    sys_cat_labels.push_back("_CSVHFStats2Up");     //15
    sys_cat_labels.push_back("_CSVHFStats2Down");   //16
    sys_cat_labels.push_back("_CSVLFStats1Up");     //17
    sys_cat_labels.push_back("_CSVLFStats1Down");   //18
    sys_cat_labels.push_back("_CSVLFStats2Up");     //19
    sys_cat_labels.push_back("_CSVLFStats2Down");   //20
 
    sys_cat_labels.push_back("_CSVCErr1Up");     //21
    sys_cat_labels.push_back("_CSVCErr1Down");   //22
    sys_cat_labels.push_back("_CSVCErr2Up");     //23
    sys_cat_labels.push_back("_CSVCErr2Down");   //24
     sys_cat_labels.push_back("_Q2scaleUp");     //25
    sys_cat_labels.push_back("_Q2scaleDown");   //26
    sys_cat_labels.push_back("_topPtcorrUp");   //27
    sys_cat_labels.push_back("_topPtcorrDown"); //28
	
	 int NumSysCat = int(sys_cat_labels.size());
	 
	 //Histos//
	 
	 TH1::SetDefaultSumw2();

	 TH1D* h_BDT_Output_IncludeWgt[NumCat][NumSysCat];
	 TH1D* h_BDT_Output_IncludeWgt_even[NumCat][NumSysCat];
	 TH1D* h_BDT_Output_IncludeWgt_odd[NumCat][NumSysCat];
	 
	 TH1D* h_LeadJetPt_IncludeWgt[NumCat][NumSysCat];
	 TH1D* h_LeadJetCSV_IncludeWgt[NumCat][NumSysCat];
	
	 TH1D* h_BDT_Output_NoWgt[NumCat][NumSysCat];
	 TH1D* h_LeadJetPt_NoWgt[NumCat][NumSysCat];
	 TH1D* h_LeadJetCSV_NoWgt[NumCat][NumSysCat];
	 
	 
     TH1D* h_SplittingOutput_1_1[NumCat][NumSysCat];
     TH1D* h_SplittingOutput_1_2[NumCat][NumSysCat];
     TH1D* h_SplittingOutput_1_1_Train[NumCat][NumSysCat];
     TH1D* h_SplittingOutput_1_2_Train[NumCat][NumSysCat];
	 
	 
	 

	 //Variable Histos//

     TH1D* h_first_jet_pt[NumCat][NumSysCat]; 
     TH1D* h_min_dr_tagged_jets[NumCat][NumSysCat];
     TH1D* h_avg_dr_tagged_jets[NumCat][NumSysCat];
     TH1D* h_aplanarity[NumCat][NumSysCat];
     TH1D* h_sphericity[NumCat][NumSysCat];
     //TH1D* h_avg_btag_disc_non_btags[NumCat][NumSysCat];
     TH1D* h_second_jet_pt[NumCat][NumSysCat];
    // TH1D* h_dr_between_lep_and_closest_jet[NumCat][NumSysCat];
     TH1D* h_h0[NumCat][NumSysCat];
     TH1D* h_avg_btag_disc_btags[NumCat][NumSysCat];
     TH1D* h_dev_from_avg_disc_btags[NumCat][NumSysCat];
     TH1D* h_third_jet_pt[NumCat][NumSysCat];
     TH1D* h_fourth_jet_pt[NumCat][NumSysCat];
     //TH1D* h_avg_tagged_dijet_mass[NumCat][NumSysCat];
     TH1D* h_h1[NumCat][NumSysCat];
     TH1D* h_h2[NumCat][NumSysCat];
     TH1D* h_all_sum_pt_with_met[NumCat][NumSysCat];
     TH1D* h_closest_tagged_dijet_mass[NumCat][NumSysCat];
     TH1D* h_h3[NumCat][NumSysCat];
    // TH1D* h_h4[NumCat][NumSysCat];
    // TH1D* h_first_highest_btag[NumCat][NumSysCat];
     TH1D* h_second_highest_btag[NumCat][NumSysCat];
     TH1D* h_third_highest_btag[NumCat][NumSysCat];
     TH1D* h_fourth_highest_btag[NumCat][NumSysCat];
     TH1D* h_fifth_highest_CSV[NumCat][NumSysCat];
     //TH1D* h_sixth_highest_CSV[NumCat][NumSysCat];
     TH1D* h_invariant_mass_of_everything[NumCat][NumSysCat];
     TH1D* h_Mlb[NumCat][NumSysCat];
     TH1D* h_tagged_dijet_mass_closest_to_125[NumCat][NumSysCat];
     TH1D* h_maxeta_jet_jet[NumCat][NumSysCat];		 
     TH1D* h_maxeta_jet_tag[NumCat][NumSysCat];		 
     TH1D* h_maxeta_tag_tag[NumCat][NumSysCat];			 
     TH1D* h_dEta_fn[NumCat][NumSysCat];			 
     TH1D* h_pt_all_jets_over_E_all_jets[NumCat][NumSysCat];  
     TH1D* h_HT[NumCat][NumSysCat];
	 
	 
  /////Special histos for KIT binning
	
 	TH1F *h1_6j2t_KIT[NumSysCat];
 	TH1F *avg_dr_tagged_jets_6j2t_KIT[NumSysCat];
 	TH1F *sphericity_6j2t_KIT[NumSysCat];
 	TH1F *third_highest_btag_6j2t_KIT[NumSysCat];
 	TH1F *h3_6j2t_KIT[NumSysCat];
 	TH1F *HT_6j2t_KIT[NumSysCat];
 	TH1F *Mlb_6j2t_KIT[NumSysCat];
 	TH1F *fifth_highest_CSV_6j2t_KIT[NumSysCat];
 	TH1F *fourth_highest_btag_6j2t_KIT[NumSysCat];
 	TH1F *h1_4j3t_KIT[NumSysCat];
 	TH1F *avg_dr_tagged_jets_4j3t_KIT[NumSysCat];
 	TH1F *sphericity_4j3t_KIT[NumSysCat];
 	TH1F *third_highest_btag_4j3t_KIT[NumSysCat];
 	TH1F *HT_4j3t_KIT[NumSysCat];
 	TH1F *dev_from_avg_disc_btags_4j3t_KIT[NumSysCat];
 	TH1F *M3_4j3t_KIT[NumSysCat];
 	TH1F *min_dr_tagged_jets_4j3t_KIT[NumSysCat];
 	TH1F *Evt_CSV_Average_4j3t_KIT[NumSysCat];
 	TH1F *h1_5j3t_KIT[NumSysCat];
 	TH1F *avg_dr_tagged_jets_5j3t_KIT[NumSysCat];
 	TH1F *sphericity_5j3t_KIT[NumSysCat];
 	TH1F *third_highest_btag_5j3t_KIT[NumSysCat];
 	TH1F *h3_5j3t_KIT[NumSysCat];
 	TH1F *HT_5j3t_KIT[NumSysCat];
 	TH1F *dev_from_avg_disc_btags_5j3t_KIT[NumSysCat];
 	TH1F *fourth_highest_btag_5j3t_KIT[NumSysCat];
 	TH1F *h1_6j3t_KIT[NumSysCat];
 	TH1F *avg_dr_tagged_jets_6j3t_KIT[NumSysCat];
 	TH1F *third_highest_btag_6j3t_KIT[NumSysCat];
 	TH1F *HT_6j3t_KIT[NumSysCat];
 	TH1F *pt_all_jets_over_E_all_jets_6j3t_KIT[NumSysCat];
 	TH1F *fifth_highest_CSV_6j3t_KIT[NumSysCat];
 	TH1F *fourth_highest_btag_6j3t_KIT[NumSysCat];
 	TH1F *min_dr_tagged_jets_6j3t_KIT[NumSysCat];
 	TH1F *h1_4j4t_KIT[NumSysCat];
 	TH1F *avg_dr_tagged_jets_4j4t_KIT[NumSysCat];
 	TH1F *sphericity_4j4t_KIT[NumSysCat];
 	TH1F *avg_btag_disc_btags_4j4t_KIT[NumSysCat];
 	TH1F *h2_4j4t_KIT[NumSysCat];
 	TH1F *closest_tagged_dijet_mass_4j4t_KIT[NumSysCat];
 	TH1F *second_highest_btag_4j4t_KIT[NumSysCat];
 	TH1F *fourth_highest_btag_4j4t_KIT[NumSysCat];
 	TH1F *maxeta_jet_jet_4j4t_KIT[NumSysCat];
 	TH1F *pt_all_jets_over_E_all_jets_4j4t_KIT[NumSysCat];
 	TH1F *avg_dr_tagged_jets_5j4t_KIT[NumSysCat];
 	TH1F *HT_5j4t_KIT[NumSysCat];
 	TH1F *M3_5j4t_KIT[NumSysCat];
 	TH1F *fifth_highest_CSV_5j4t_KIT[NumSysCat];
 	TH1F *fourth_highest_btag_5j4t_KIT[NumSysCat];
 	TH1F *Evt_Deta_JetsAverage_5j4t_KIT[NumSysCat];
 	TH1F *avg_btag_disc_btags_5j4t_KIT[NumSysCat];
 	TH1F *h2_6j4t_KIT[NumSysCat];
 	TH1F *avg_dr_tagged_jets_6j4t_KIT[NumSysCat];
 	TH1F *third_highest_btag_6j4t_KIT[NumSysCat];
 	TH1F *M3_6j4t_KIT[NumSysCat];
 	TH1F *fifth_highest_CSV_6j4t_KIT[NumSysCat];
 	TH1F *fourth_highest_btag_6j4t_KIT[NumSysCat];
 	TH1F *best_higgs_mass_6j4t_KIT[NumSysCat];
 	TH1F *dEta_fn_6j4t_KIT[NumSysCat];
     
	////// 
	
	TH1D *h_NumCount;
	TH1D *h_NumCountRaw;
	
	h_NumCount = new TH1D("h_NumCount","h_NumCount",10,0,10);
	h_NumCountRaw = new TH1D("h_NumCountRaw","h_NumCountRaw",10,0,10);
	 

	 ///
	 
	 double jetptmax = 500.;
	 int NjetptBins = int ( jetptmax/10. + 0.0001);
	 double massmax = 500. ;
	 int NmassBins = int( massmax/10. + 0.0001);
	 int NcsvBins = 44;
	 double htmax = 2000;
	 int NhtBins = int ( htmax/50. + 0.0001 );
	 
	 
     for( int b=0; b<NumSysCat; b++ ){
		 for( int c=0; c<NumCat; c++ ){ 
         std::string suffix = "_" + cat_labels[c];
         suffix += sys_cat_labels[b];
  	
	h_BDT_Output_IncludeWgt[c][b] = new TH1D((std::string("h_BDT_Output_IncludeWgt" + suffix)).c_str(),(std::string("BDT_Output" + suffix)).c_str(),20,-1,1);
	h_BDT_Output_IncludeWgt_even[c][b] = new TH1D((std::string("h_BDT_Output_IncludeWgt_even" + suffix)).c_str(),(std::string("BDT_Output_even" + suffix)).c_str(),20,-1,1);
	h_BDT_Output_IncludeWgt_odd[c][b] = new TH1D((std::string("h_BDT_Output_IncludeWgt_odd" + suffix)).c_str(),(std::string("BDT_Output_odd" + suffix)).c_str(),20,-1,1);
  	
    h_SplittingOutput_1_1[c][b] = new TH1D((std::string("h_SplittingOutput_1_1" + suffix)).c_str(),"h_SplittingOutput_1_1",20,-1,1);
    h_SplittingOutput_1_2[c][b] = new TH1D((std::string("h_SplittingOutput_1_2" + suffix)).c_str(),"h_SplittingOutput_1_2",20,-1,1);
    
	h_SplittingOutput_1_1_Train[c][b] = new TH1D((std::string("h_SplittingOutput_1_1_Train" + suffix)).c_str(),"h_SplittingOutput_1_1_Train",20,-1,1);
    h_SplittingOutput_1_2_Train[c][b] = new TH1D((std::string("h_SplittingOutput_1_2_Train" + suffix)).c_str(),"h_SplittingOutput_1_2_Train",20,-1,1);
	
	
	h_LeadJetPt_IncludeWgt[c][b] = new TH1D((std::string("h_LeadJetPt_IncludeWgt" + suffix)).c_str(),(std::string("LeadJetPt" + suffix)).c_str(),20,0.0,500.0);
  	h_LeadJetCSV_IncludeWgt[c][b] = new TH1D((std::string("h_LeadJetCSV_IncludeWgt" + suffix)).c_str(),(std::string("LeadjetCSV" + suffix)).c_str(),20,0.0,1.0);
	
	h_BDT_Output_NoWgt[c][b] = new TH1D((std::string("h_BDT_Output_NoWgt" + suffix)).c_str(),(std::string("BDT_Output" + suffix)).c_str(),20,-1,1);
  	h_LeadJetPt_NoWgt[c][b] = new TH1D((std::string("h_LeadJetPt_NoWgt" + suffix)).c_str(),(std::string("LeadJetPt" + suffix)).c_str(),20,0.0,500.0);
  	h_LeadJetCSV_NoWgt[c][b] = new TH1D((std::string("h_LeadJetCSV_NoWgt" + suffix)).c_str(),(std::string("LeadjetCSV" + suffix)).c_str(),20,0.0,1.0);
	
	h_first_jet_pt[c][b] = new TH1D((std::string("h_first_jet_pt" + suffix)).c_str(),";first jet pt", NjetptBins, 0, jetptmax );
      	h_second_jet_pt[c][b] = new TH1D((std::string("h_second_jet_pt" + suffix)).c_str(),";second jet pt", NjetptBins, 0, jetptmax );
        h_third_jet_pt[c][b] = new TH1D((std::string("h_third_jet_pt" + suffix)).c_str(),";third jet pt", NjetptBins, 0, jetptmax );
        h_fourth_jet_pt[c][b] = new TH1D((std::string("h_fourth_jet_pt" + suffix)).c_str(),";fourth jet pt", NjetptBins, 0, jetptmax );

      h_aplanarity[c][b] = new TH1D((std::string("h_aplanarity" + suffix)).c_str(),";aplanarity", 25, 0, 0.5 );
      h_sphericity[c][b] = new TH1D((std::string("h_sphericity" + suffix)).c_str(),";sphericity", 25, 0, 1.0 );
     // h_avg_btag_disc_non_btags[c][b] = new TH1D((std::string("h_avg_btag_disc_non_btags" + suffix)).c_str(),";ave btag disc non btags", 230, -10.5, 1.0 );
      h_avg_btag_disc_btags[c][b] = new TH1D((std::string("h_avg_btag_disc_btags" + suffix)).c_str(),";ave btag disc btags", 100, -0.01, 1.01 );
      h_dev_from_avg_disc_btags[c][b] = new TH1D((std::string("h_dev_from_avg_disc_btags" + suffix)).c_str(),";dev from avg disc btags", 25, -0.0001, 0.04 );
      h_closest_tagged_dijet_mass[c][b] = new TH1D((std::string("h_closest_tagged_dijet_mass" + suffix)).c_str(),";mass of closest tagged jets", NmassBins, 0, massmax );

      h_h0[c][b] = new TH1D((std::string("h_h0" + suffix)).c_str(),";h0", 25, 0, 0.5 );
      h_h1[c][b] = new TH1D((std::string("h_h1" + suffix)).c_str(),";h1", 25, -0.3, 0.5 );
      h_h2[c][b] = new TH1D((std::string("h_h2" + suffix)).c_str(),";h2", 25, -0.2, 0.5 );
      h_h3[c][b] = new TH1D((std::string("h_h3" + suffix)).c_str(),";h3", 25, -0.2, 1.1 );
   //   h_h4[c][b] = new TH1D((std::string("h_h4" + suffix)).c_str(),";h4", 25, -0.2, 0.3 );

      h_maxeta_jet_jet[c][b] = new TH1D((std::string("h_maxeta_jet_jet" + suffix)).c_str(),";max #Delta#eta(jet,ave jet #eta)", 50, 0., 6. );
      h_maxeta_jet_tag[c][b] = new TH1D((std::string("h_maxeta_jet_tag" + suffix)).c_str(),";max #Delta#eta(tag,ave jet #eta)", 50, 0., 6. );
      h_maxeta_tag_tag[c][b] = new TH1D((std::string("h_maxeta_tag_tag" + suffix)).c_str(),";max #Delta#eta(tag,ave tag #eta)", 50, 0., 6. );
      h_dEta_fn[c][b] = new TH1D((std::string("h_dEta_fn" + suffix)).c_str(),";product(#Delta#eta(leptonic top, bb),#Delta#eta(hadronic top, bb))", 20, 0., 6.);
      h_pt_all_jets_over_E_all_jets[c][b] = new TH1D((std::string("h_pt_all_jets_over_E_all_jets" + suffix)).c_str(),";(#Sigma jet p_{T})/(#Sigma jet E)",50,0.,1.);

    //  h_jet_csv_1[c][b] = new TH1D((std::string("h_jet_csv_1" + suffix)).c_str(),";jet 1 CSV", NcsvBins, -1.1, 1.1 );
      h_second_highest_btag[c][b] = new TH1D((std::string("h_second_highest_btag" + suffix)).c_str(),";jet 2 CSV", NcsvBins, -1.1, 1.1 );
      h_third_highest_btag[c][b] = new TH1D((std::string("h_third_highest_btag" + suffix)).c_str(),";jet 3 CSV", NcsvBins, -1.1, 1.1 );
      h_fourth_highest_btag[c][b] = new TH1D((std::string("h_fourth_highest_btag" + suffix)).c_str(),";jet 4 CSV", NcsvBins, -1.1, 1.1 );
      h_fifth_highest_CSV[c][b] = new TH1D((std::string("h_fifth_highest_CSV" + suffix)).c_str(),";jet 5 CSV", NcsvBins, -1.1, 1.1 );
    
      h_min_dr_tagged_jets[c][b] = new TH1D((std::string("h_min_dr_tagged_jets" + suffix)).c_str(),";min #DeltaR(tag,tag)", 30, 0, 6. );
      h_avg_dr_tagged_jets[c][b] = new TH1D((std::string("h_avg_dr_tagged_jets" + suffix)).c_str(),";ave #DeltaR(tag,tag)", 30, 0, 6. );
	
      h_HT[c][b] = new TH1D((std::string("h_HT" + suffix)).c_str(),";H_{T} (jets)", NhtBins, 0, htmax );
      h_all_sum_pt_with_met[c][b] = new TH1D((std::string("h_all_sum_pt_with_met" + suffix)).c_str(),";sum p_{T} (jets,lepton,MET)", NhtBins, 0, htmax );
      h_invariant_mass_of_everything[c][b] = new TH1D((std::string("h_invariant_mass_of_everything" + suffix)).c_str(),";M_{inv} (jets,lepton,MET)", NhtBins, 0, htmax );
     
      h_Mlb[c][b] = new TH1D((std::string("h_Mlb" + suffix)).c_str(),";M(lep,tag)", NmassBins, 0, massmax );
      h_tagged_dijet_mass_closest_to_125[c][b] = new TH1D((std::string("h_tagged_dijet_mass_closest_to_125" + suffix)).c_str(),";tagged dijet mass closest to 125", NmassBins, 0, massmax );

	 }
	 
	  std::string short_suffix = sys_cat_labels[b];
	  
	      //SpecialKIT binning 
    
	  h1_6j2t_KIT[b] = new TH1F((std::string("h1_6j2t_KIT" + short_suffix)).c_str(),"h1 (6j2t)",20,-0.175,0.325);    
	  avg_dr_tagged_jets_6j2t_KIT[b] = new TH1F((std::string("avg_dr_tagged_jets_6j2t_KIT" + short_suffix)).c_str(),"avg_dr_tagged_jets (6j2t)",19,0.0,4.75);
	  sphericity_6j2t_KIT[b] = new TH1F((std::string("sphericity_6j2t_KIT" + short_suffix)).c_str(),"sphericity (6j2t)",18,0.0,0.9);
	  third_highest_btag_6j2t_KIT[b] = new TH1F((std::string("third_highest_btag_6j2t_KIT" + short_suffix)).c_str(),"third_highest_btag (6j2t)",18,0.0,0.9);
	  h3_6j2t_KIT[b] = new TH1F((std::string("h3_6j2t_KIT" + short_suffix)).c_str(),"h3 (6j2t)",19,0.0,0.95);
	  HT_6j2t_KIT[b] = new TH1F((std::string("HT_6j2t_KIT" + short_suffix)).c_str(),"HT (6j2t)",17,0.0,1700.0);
	  Mlb_6j2t_KIT[b] = new TH1F((std::string("Mlb_6j2t_KIT" + short_suffix)).c_str(),"Mlb (6j2t)",20,0.0,400.0);
	  fifth_highest_CSV_6j2t_KIT[b] = new TH1F((std::string("fifth_highest_CSV_6j2t_KIT" + short_suffix)).c_str(),"fifth_highest_CSV (6j2t)",13,0.0,0.65);
	  fourth_highest_btag_6j2t_KIT[b] = new TH1F((std::string("fourth_highest_btag_6j2t_KIT" + short_suffix)).c_str(),"fourth_highest_btag (6j2t)",17,0.0,0.85);
	  h1_4j3t_KIT[b] = new TH1F((std::string("h1_4j3t_KIT" + short_suffix)).c_str(),"h1 (4j3t)",11,-0.2,0.35);
	  avg_dr_tagged_jets_4j3t_KIT[b] = new TH1F((std::string("avg_dr_tagged_jets_4j3t_KIT" + short_suffix)).c_str(),"avg_dr_tagged_jets (4j3t)",19,0.0,3.8);
	  sphericity_4j3t_KIT[b] = new TH1F((std::string("sphericity_4j3t_KIT" + short_suffix)).c_str(),"sphericity (4j3t)",18,0.0,0.9);
	  third_highest_btag_4j3t_KIT[b] = new TH1F((std::string("third_highest_btag_4j3t_KIT" + short_suffix)).c_str(),"third_highest_btag (4j3t)",11,0.88,0.99);
	  HT_4j3t_KIT[b] = new TH1F((std::string("HT_4j3t_KIT" + short_suffix)).c_str(),"HT (4j3t)",19,0.0,950.0);
	  dev_from_avg_disc_btags_4j3t_KIT[b] = new TH1F((std::string("dev_from_avg_disc_btags_4j3t_KIT" + short_suffix)).c_str(),"dev_from_avg_disc_btags (4j3t)",12,0.0,0.0024);
	  M3_4j3t_KIT[b] = new TH1F((std::string("M3_4j3t_KIT" + short_suffix)).c_str(),"M3 (4j3t)",17,0.0,850.0);
	  min_dr_tagged_jets_4j3t_KIT[b] = new TH1F((std::string("min_dr_tagged_jets_4j3t_KIT" + short_suffix)).c_str(),"min_dr_tagged_jets (4j3t)",16,0.0,3.2);
	  Evt_CSV_Average_4j3t_KIT[b] = new TH1F((std::string("Evt_CSV_Average_4j3t_KIT" + short_suffix)).c_str(),"Evt_CSV_Average (4j3t)",13,0.68,0.94);
	  h1_5j3t_KIT[b] = new TH1F((std::string("h1_5j3t_KIT" + short_suffix)).c_str(),"h1 (5j3t)",20,-0.2,0.3);
	  avg_dr_tagged_jets_5j3t_KIT[b] = new TH1F((std::string("avg_dr_tagged_jets_5j3t_KIT" + short_suffix)).c_str(),"avg_dr_tagged_jets (5j3t)",18,0.0,3.6);
	  sphericity_5j3t_KIT[b] = new TH1F((std::string("sphericity_5j3t_KIT" + short_suffix)).c_str(),"sphericity (5j3t)",18,0.0,0.9);
	  third_highest_btag_5j3t_KIT[b] = new TH1F((std::string("third_highest_btag_5j3t_KIT" + short_suffix)).c_str(),"third_highest_btag (5j3t)",12,0.88,1.0);
	  h3_5j3t_KIT[b] = new TH1F((std::string("h3_5j3t_KIT" + short_suffix)).c_str(),"h3 (5j3t)",20,-0.05,0.95);
	  HT_5j3t_KIT[b] = new TH1F((std::string("HT_5j3t_KIT" + short_suffix)).c_str(),"HT (5j3t)",12,0.0,1200.0);
	  dev_from_avg_disc_btags_5j3t_KIT[b] = new TH1F((std::string("dev_from_avg_disc_btags_5j3t_KIT" + short_suffix)).c_str(),"dev_from_avg_disc_btags (5j3t)",12,0.0,0.0024);
	  fourth_highest_btag_5j3t_KIT[b] = new TH1F((std::string("fourth_highest_btag_5j3t_KIT" + short_suffix)).c_str(),"fourth_highest_btag (5j3t)",17,0.0,0.85);
	  h1_6j3t_KIT[b] = new TH1F((std::string("h1_6j3t_KIT" + short_suffix)).c_str(),"h1 (6j3t)",20,-0.175,0.325);
	  avg_dr_tagged_jets_6j3t_KIT[b] = new TH1F((std::string("avg_dr_tagged_jets_6j3t_KIT" + short_suffix)).c_str(),"avg_dr_tagged_jets (6j3t)",19,0.0,3.8);
	  third_highest_btag_6j3t_KIT[b] = new TH1F((std::string("third_highest_btag_6j3t_KIT" + short_suffix)).c_str(),"third_highest_btag (6j3t)",12,0.88,1.0);
	  HT_6j3t_KIT[b] = new TH1F((std::string("HT_6j3t_KIT" + short_suffix)).c_str(),"HT (6j3t)",19,0.0,1900.0);
	  pt_all_jets_over_E_all_jets_6j3t_KIT[b] = new TH1F((std::string("pt_all_jets_over_E_all_jets_6j3t_KIT" + short_suffix)).c_str(),"pt_all_jets_over_E_all_jets (6j3t)",14,0.2,0.9);
	  fifth_highest_CSV_6j3t_KIT[b] = new TH1F((std::string("fifth_highest_CSV_6j3t_KIT" + short_suffix)).c_str(),"fifth_highest_CSV (6j3t)",13,0.0,0.65);
	  fourth_highest_btag_6j3t_KIT[b] = new TH1F((std::string("fourth_highest_btag_6j3t_KIT" + short_suffix)).c_str(),"fourth_highest_btag (6j3t)",18,0.0,0.9);
	  min_dr_tagged_jets_6j3t_KIT[b] = new TH1F((std::string("min_dr_tagged_jets_6j3t_KIT" + short_suffix)).c_str(),"min_dr_tagged_jets (6j3t)",17,0.0,3.4);
	  h1_4j4t_KIT[b] = new TH1F((std::string("h1_4j4t_KIT" + short_suffix)).c_str(),"h1 (4j4t)",19,-0.2,0.275);
	  avg_dr_tagged_jets_4j4t_KIT[b] = new TH1F((std::string("avg_dr_tagged_jets_4j4t_KIT" + short_suffix)).c_str(),"avg_dr_tagged_jets (4j4t)",15,1.4,2.9);
	  sphericity_4j4t_KIT[b] = new TH1F((std::string("sphericity_4j4t_KIT" + short_suffix)).c_str(),"sphericity (4j4t)",17,0.0,0.85);
	  avg_btag_disc_btags_4j4t_KIT[b] = new TH1F((std::string("avg_btag_disc_btags_4j4t_KIT" + short_suffix)).c_str(),"avg_btag_disc_btags (4j4t)",16,0.915,0.995);
	  h2_4j4t_KIT[b] = new TH1F((std::string("h2_4j4t_KIT" + short_suffix)).c_str(),"h2 (4j4t)",20,-0.16,0.24);
	  closest_tagged_dijet_mass_4j4t_KIT[b] = new TH1F((std::string("closest_tagged_dijet_mass_4j4t_KIT" + short_suffix)).c_str(),"closest_tagged_dijet_mass (4j4t)",11,0.0,220.0);
	  second_highest_btag_4j4t_KIT[b] = new TH1F((std::string("second_highest_btag_4j4t_KIT" + short_suffix)).c_str(),"second_highest_btag (4j4t)",16,0.92,1.0);
	  fourth_highest_btag_4j4t_KIT[b] = new TH1F((std::string("fourth_highest_btag_4j4t_KIT" + short_suffix)).c_str(),"fourth_highest_btag (4j4t)",11,0.88,0.99);
	  maxeta_jet_jet_4j4t_KIT[b] = new TH1F((std::string("maxeta_jet_jet_4j4t_KIT" + short_suffix)).c_str(),"maxeta_jet_jet (4j4t)",15,0.0,1.5);
	  pt_all_jets_over_E_all_jets_4j4t_KIT[b] = new TH1F((std::string("pt_all_jets_over_E_all_jets_4j4t_KIT" + short_suffix)).c_str(),"pt_all_jets_over_E_all_jets (4j4t)",14,0.2,0.9);
	  avg_dr_tagged_jets_5j4t_KIT[b] = new TH1F((std::string("avg_dr_tagged_jets_5j4t_KIT" + short_suffix)).c_str(),"avg_dr_tagged_jets (5j4t)",18,1.3,3.1);
	  HT_5j4t_KIT[b] = new TH1F((std::string("HT_5j4t_KIT" + short_suffix)).c_str(),"HT (5j4t)",11,0.0,1100.0);
	  M3_5j4t_KIT[b] = new TH1F((std::string("M3_5j4t_KIT" + short_suffix)).c_str(),"M3 (5j4t)",17,0.0,850.0);
	  fifth_highest_CSV_5j4t_KIT[b] = new TH1F((std::string("fifth_highest_CSV_5j4t_KIT" + short_suffix)).c_str(),"fifth_highest_CSV (5j4t)",19,-0.05,0.9);
	  fourth_highest_btag_5j4t_KIT[b] = new TH1F((std::string("fourth_highest_btag_5j4t_KIT" + short_suffix)).c_str(),"fourth_highest_btag (5j4t)",11,0.88,0.99);
	  Evt_Deta_JetsAverage_5j4t_KIT[b] = new TH1F((std::string("Evt_Deta_JetsAverage_5j4t_KIT" + short_suffix)).c_str(),"Evt_Deta_JetsAverage (5j4t)",12,0.0,2.4);
	  avg_btag_disc_btags_5j4t_KIT[b] = new TH1F((std::string("avg_btag_disc_btags_5j4t_KIT" + short_suffix)).c_str(),"avg_btag_disc_btags (5j4t)",18,0.9,0.99);
	  h2_6j4t_KIT[b] = new TH1F((std::string("h2_6j4t_KIT" + short_suffix)).c_str(),"h2 (6j4t)",20,-0.12,0.28);
	  avg_dr_tagged_jets_6j4t_KIT[b] = new TH1F((std::string("avg_dr_tagged_jets_6j4t_KIT" + short_suffix)).c_str(),"avg_dr_tagged_jets (6j4t)",16,0.0,3.2);
	  third_highest_btag_6j4t_KIT[b] = new TH1F((std::string("third_highest_btag_6j4t_KIT" + short_suffix)).c_str(),"third_highest_btag (6j4t)",11,0.89,1.0);
	  M3_6j4t_KIT[b] = new TH1F((std::string("M3_6j4t_KIT" + short_suffix)).c_str(),"M3 (6j4t)",13,0.0,1300.0);
	  fifth_highest_CSV_6j4t_KIT[b] = new TH1F((std::string("fifth_highest_CSV_6j4t_KIT" + short_suffix)).c_str(),"fifth_highest_CSV (6j4t)",19,0.0,0.95);
	  fourth_highest_btag_6j4t_KIT[b] = new TH1F((std::string("fourth_highest_btag_6j4t_KIT" + short_suffix)).c_str(),"fourth_highest_btag (6j4t)",11,0.88,0.99);
	  best_higgs_mass_6j4t_KIT[b] = new TH1F((std::string("best_higgs_mass_6j4t_KIT" + short_suffix)).c_str(),"best_higgs_mass (6j4t)",12,0.0,1200.0);
	  dEta_fn_6j4t_KIT[b] = new TH1F((std::string("dEta_fn_6j4t_KIT" + short_suffix)).c_str(),"dEta_fn (6j4t)",18,0.0,4.5);

	  ///
	
	
	
	}

    //////////////////////////////////////////////////////////////////////////
    ///  Weights
    //////////////////////////////////////////////////////////////////////////
 

    int nentries = chain->GetEntries();
    std::cout << "\n\t Number of entries = " << nentries << std::endl;
    std::cout << "\t Max number of entries = " << maxNentries << std::endl;
    std::cout << "\n" << std::endl;

    int use_nentries = std::max( maxNentries, nentries);

    int NeventsPerJob = int( double(use_nentries)/double(Njobs) + 0.000001 ) + 1;

    int firstEvent = (jobN-1)*NeventsPerJob + 1;
    int lastEvent  = firstEvent + NeventsPerJob;
    if( jobN==Njobs ) lastEvent = -1;
    if( jobN==1 ) firstEvent = 0;
  
  
    double MyWgt = 1;
    double ttJet_xSec = 832;
    //double ttH_xSec = 0.5085;
	double ttHbb_xSec = 0.2934;
	double ttHnobb_xSec = 0.2151;
    double ttZ_xSec = 2.232;
    double ttW_xSec = 1.152;
    double nGen = nentries;
    double nGenNeg=1580920; //for amc@nlo sample

	double MyttHwgtNeg = (ttHnobb_xSec/nGenNeg)*10000;
	double MyttHwgt = ttHbb_xSec/nGen*10000; 
	double MyttJetwgt = (ttJet_xSec/nGen)*10000;
	double MyttZwgt = (ttZ_xSec/nGen)*10000;
	double MyttWwgt = (ttW_xSec/nGen)*10000;




  cout<<" BEGIN wgt ="<<MyttHwgtNeg<<endl;


  if(SampleType ==1)MyWgt=MyttHwgtNeg;
  if(SampleType ==2)MyWgt=MyttJetwgt;
  if(SampleType ==3)MyWgt=MyttHwgt;
  if(SampleType ==4)MyWgt=MyttZwgt;
  if(SampleType ==5)MyWgt=MyttWwgt;
  
  int cut[7];
	
    std::cout << "========  Starting Event Loop  ========" << std::endl;
    for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    //Long64_t

	    if( ievt<firstEvent ) continue;
	    if( ievt==lastEvent ) break;

      	if( ievt==0 )        std::cout << "     Event " << ievt+1 << std::endl;
      	if( ievt%10000==0 && ievt!=1) std::cout << "           " << ievt << "\t" 
  					     << int(double(ievt-1)/double(nentries)*100) << "% done" << std::endl;
	  
	 	if( ievt==(maxNentries+1) && ievt!=0 ) break;
	  
		cut[0]=0;
		cut[1]=0;
		cut[2]=0;
		cut[3]=0;
		cut[4]=0;
		cut[5]=0;
		cut[6]=0;
		
		cut[0]++;
		chain->GetEntry(ievt);  
		long evt = eve->evt_;
	  	
		
		
		bool GoodFirstPV = eve->GoodFirstPV_;
		if(GoodFirstPV)cut[1]++;
		
		// LJ Event selection
		
		int passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1 = eve->passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1_;
		int passHLT_IsoMu24_eta2p1_v = eve->passHLT_IsoMu24_eta2p1_v_;
		if((passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1==1 || passHLT_IsoMu24_eta2p1_v==1) && cut[1]==1)cut[2]++;
		
		int numTightElectrons = eve->numTightElectrons_;
		int numLooseElectrons = eve->numLooseElectrons_;
		int numTightMuons = eve->numTightMuons_;
		int numLooseMuons = eve->numLooseMuons_;
		int numLooselep=-1;
		int numTightlep=-1;

		if((numLooseMuons + numLooseElectrons)==1)numLooselep=1;
		if(((numLooseMuons==1 && passHLT_IsoMu24_eta2p1_v==1) || (numLooseElectrons==1 && passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1==1)) && cut[2]==1)cut[3]++;
		if(((numTightElectrons==1 && passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1==1 && numLooseElectrons==1) || (numTightMuons==1 && passHLT_IsoMu24_eta2p1_v==1 && numLooseMuons==1)) && cut[3]==1)cut[4]++;
		if(((numTightElectrons==1 && passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1==1) || (numTightMuons==1 && passHLT_IsoMu24_eta2p1_v==1)) && numLooselep==1 && cut[4]==1)cut[5]++;
		
		
		// Cat Event selection
		
		int numJets = int( eve->numJets_float_[0] + 0.001 );
		int numTags = int( eve->numTags_float_[0] + 0.001 );
		
		if(((numJets>=4 && numTags>=3) || (numJets>=6 && numTags>=2)) && cut[5]==1)cut[6]++;
		
		for(int i=0;i<7;i++){
		  if(cut[i]==1) h_NumCountRaw->Fill(i+1);
	    }
		
		if(cut[5]==0)continue; // Lepton + Jet selection.
		
		
	    // Filter by decay

		if(SampleType == 2){
		int additionalJetEventId			 = eve->additionalJetEventId_;
		if(SplitType ==1 && !(additionalJetEventId == 51) )continue; //ttb
		if(SplitType ==3 && !(additionalJetEventId == 52) )continue; //ttB
		if(SplitType ==2 && !(additionalJetEventId == 53 || additionalJetEventId == 54 || additionalJetEventId == 55))continue; //ttbb
		if(SplitType ==4 && !(additionalJetEventId == 41 || additionalJetEventId == 42 || additionalJetEventId == 43 || additionalJetEventId == 44 || additionalJetEventId == 45))continue; //ttc
		//if(SplitType ==5 && !(additionalJetEventId == 43 || additionalJetEventId == 44 || additionalJetEventId == 45))continue; //ttcc
		if(SplitType ==5 && !(additionalJetEventId == 0))continue; //ttlf
	    }
	
		if(SampleType == 1){
			int Hdecay = eve->higgsDecayType_;
			if(SplitType==1 && !(Hdecay == 1))continue; //bb
			if(SplitType==2 && !(Hdecay == 2))continue; //WW
			if(SplitType==3 && !(Hdecay == 3))continue; //TauTau
			if(SplitType==4 && !(Hdecay == 4))continue; //gluglu
			if(SplitType==5 && !(Hdecay == 5))continue; //cc
			if(SplitType==6 && !(Hdecay == 6))continue; //ZZ
			if(SplitType==7 && !(Hdecay == 7))continue; //Zgamma
			if(SplitType==8 && !(Hdecay == 8))continue; //gammagamma
			if(SplitType==10 && !(Hdecay == 10))continue; //ss
			if(SplitType==11 && !(Hdecay == 11))continue; //mumu
		}
	
	
		
		
		
	    int iSys=-1;
	    for(int treeSys=0; treeSys<rNumSys; treeSys++){
	    HT=0;
	      iSys++;
		

  		// Cat Event selection
		
  		numJets = int( eve->numJets_float_[treeSys] + 0.001 );
  		numTags = int( eve->numTags_float_[treeSys] + 0.001 );
		
		
    




		
	
		// Grab variables from YggdrasilTrees
		
	    int this_category = -1;
	    if( numJets==4 && numTags==2) this_category=0;
	    if( numJets==5 && numTags==2) this_category=1;
	    if( numJets>=6 && numTags==2) this_category=2;	
	    if( numJets==4 && numTags==3) this_category=3;
	    if( numJets==5 && numTags==3) this_category=4;
	    if( numJets>=6 && numTags==3) this_category=5;
	    if( numJets==4 && numTags>=4) this_category=6;
	    if( numJets==5 && numTags>=4) this_category=7;
	    if( numJets>=6 && numTags>=4) this_category=8;
	    
	   
  
	
		if(this_category<2)continue;
		
		int run = eve->run_;
		int lumi = eve->lumi_;
	
		//Associate BDT variables

	    aplanarity             	 	 = eve->aplanarity_[treeSys];
	    sphericity              	 	 = eve->sphericity_[treeSys];
	    second_highest_btag   			 = eve->second_highest_btag_[treeSys];
	    third_highest_btag     		 = eve->third_highest_btag_[treeSys];
	    fourth_highest_btag     		 = eve->fourth_highest_btag_[treeSys]; 
	    Mlb   							 = eve->Mlb_[treeSys];
		first_jet_pt          	  		 = eve->first_jet_pt_[treeSys];
		second_jet_pt          		 = eve->second_jet_pt_[treeSys];
	    third_jet_pt           	  	 = eve->third_jet_pt_[treeSys];
	    fourth_jet_pt           		 = eve->fourth_jet_pt_[treeSys];
		avg_btag_disc_btags    		 = eve->avg_btag_disc_btags_[treeSys]; 
		h0								 = eve->h0_[treeSys];
		h1								 = eve->h1_[treeSys];
		h2								 = eve->h2_[treeSys];
		h3								 = eve->h3_[treeSys];
		avg_dr_tagged_jets      		 = eve->avg_dr_tagged_jets_[treeSys];
		dev_from_avg_disc_btags  		 = eve->dev_from_avg_disc_btags_[treeSys];
		closest_tagged_dijet_mass		 = eve->closest_tagged_dijet_mass_[treeSys];			
		invariant_mass_of_everything 	 = eve->invariant_mass_of_everything_[treeSys];
		min_dr_tagged_jets      		 = eve->min_dr_tagged_jets_[treeSys];
		tagged_dijet_mass_closest_to_125= eve->tagged_dijet_mass_closest_to_125_[treeSys];
		maxeta_jet_jet					 = eve->maxeta_jet_jet_[treeSys];  
		maxeta_jet_tag					 = eve->maxeta_jet_tag_[treeSys];  
		maxeta_tag_tag					 = eve->maxeta_tag_tag_[treeSys]; 
		dEta_fn							 = eve->dEta_f_[treeSys];
		pt_all_jets_over_E_all_jets	 = eve->pt_all_jets_over_E_all_jets_[treeSys];
		fifth_highest_CSV				 = eve->fifth_highest_CSV_[treeSys];
		best_higgs_mass			= eve->best_higgs_mass_[treeSys];
		M3				= eve->M3_[treeSys];
		//sixth_highest_CSV_f                = eve->sixth_highest_CSV_[treeSys];
		//HT_f 								=eve->HT_[treeSys];
		//all_sum_pt_with_met_f			= eve->all_sum_pt_with_met_[treeSys];
		
	
		//If running over yggdrasil_trees uncomment this.

	
		//need to get proper HT and KIT variables
		HT=0;
		vvdouble jet_vect_TLV = eve->jet_vect_TLV_[treeSys];
   		vdouble jet_CSV = eve->jet_CSV_[treeSys];
   		vint    jet_flavour = eve->jet_flavour_[treeSys];
		
		double maxJetPt=-1;
		double maxJetPtCSV=-1;
		
		double sumJetEta = 0;
		double sumTagEta = 0;
		double cntJetEta = 0;
		double cntTagEta = 0;
		double sumCSV =0;
		double cntCSV =0;
		double sumDeta =0;
		double cntDeta =0;
		
		for( int iJet=0; iJet<int(jet_vect_TLV.size()); iJet++ ){
		  TLorentzVector myJet;
		  myJet.SetPxPyPzE( jet_vect_TLV[iJet][0], jet_vect_TLV[iJet][1], jet_vect_TLV[iJet][2], jet_vect_TLV[iJet][3] );
		  double myJetPT = myJet.Pt();
		  double myCSV = jet_CSV[iJet];
		  sumCSV = sumCSV + myCSV;
		  cntCSV +=1;
		  HT += myJetPT;
		  if(myJetPT>maxJetPt){
		  	maxJetPt=myJetPT;
		  	maxJetPtCSV=myCSV;
		  }
 
          for( int jJet =iJet;jJet<int(jet_vect_TLV.size()); jJet++){
                   if(iJet==jJet)continue;
                   TLorentzVector myJetJ;
                   myJetJ.SetPxPyPzE(jet_vect_TLV[jJet][0],jet_vect_TLV[jJet][1],jet_vect_TLV[jJet][2],jet_vect_TLV[jJet][3]);
                  // sumDabseta = abs(myJet.Eta())-abs(myJetJ.Eta());
                   sumDeta = sumDeta + abs(myJet.Eta()-myJetJ.Eta());
                   cntDeta +=1;
			   }
	 	 }
	  
	   double MET                       = eve->MET_[treeSys];
	   double tight_lepton_pt           = eve->tight_lepton_pt_;
	   
	  all_sum_pt_with_met =HT + tight_lepton_pt + MET;
	  CSV_Average = sumCSV/cntCSV;
	  Deta_JetsAverage = sumDeta/cntDeta;
	  
		
	
	
		
		double csvWgtHF, csvWgtLF, csvWgtCF;
		double newCSVwgt = ( insample<0 ) ? 1 : get_csv_wgt(jet_vect_TLV,jet_CSV,jet_flavour,treeSys, csvWgtHF, csvWgtLF, csvWgtCF);
		double wgt=MyWgt;
		
		//Dealing with negWgts
		double GenWgt = eve->wgt_generator_;
		if(SampleType ==1)wgt=MyWgt*GenWgt;
		//
		
		wgt=wgt*newCSVwgt;
		
	
		
		
		double wgthalf=wgt*2;  //halving the amount of events

	
		double Reader_Output = reader_Spring15_KITreoptimized_BDT[this_category]->EvaluateMVA( "BDT method" );
		double SplittingOutput_1_1 = reader_Splitting_1_1[this_category]->EvaluateMVA("BDT method");
		double SplittingOutput_1_2 = reader_Splitting_1_2[this_category]->EvaluateMVA("BDT method");
		
		h_BDT_Output_IncludeWgt[this_category][iSys]->Fill(Reader_Output,wgt);
	  	h_LeadJetPt_IncludeWgt[this_category][iSys]->Fill(maxJetPt,wgt);
	  	h_LeadJetCSV_IncludeWgt[this_category][iSys]->Fill(maxJetPtCSV,wgt);
		
		if(int(evt) % 2 !=0){
			h_BDT_Output_IncludeWgt_odd[this_category][iSys]->Fill(Reader_Output,wgthalf);
			h_SplittingOutput_1_1[this_category][iSys]->Fill(SplittingOutput_1_1,wgthalf);
			h_SplittingOutput_1_2_Train[this_category][iSys]->Fill(SplittingOutput_1_2,wgthalf);
		}
		if(int(evt) % 2 ==0){
			h_BDT_Output_IncludeWgt_even[this_category][iSys]->Fill(Reader_Output,wgthalf);
			h_SplittingOutput_1_2[this_category][iSys]->Fill(SplittingOutput_1_2,wgthalf);
			h_SplittingOutput_1_1_Train[this_category][iSys]->Fill(SplittingOutput_1_1,wgthalf);
		}
		
		h_BDT_Output_NoWgt[this_category][iSys]->Fill(Reader_Output,1);
	  	h_LeadJetPt_NoWgt[this_category][iSys]->Fill(maxJetPt,1);
	  	h_LeadJetCSV_NoWgt[this_category][iSys]->Fill(maxJetPtCSV,1);
		
		  h_first_jet_pt[this_category][iSys]->Fill(first_jet_pt,wgt); 
      h_min_dr_tagged_jets[this_category][iSys]->Fill(min_dr_tagged_jets,wgt);
      h_avg_dr_tagged_jets[this_category][iSys]->Fill(avg_dr_tagged_jets,wgt);
      h_aplanarity[this_category][iSys]->Fill(aplanarity,wgt);
      h_sphericity[this_category][iSys]->Fill(sphericity,wgt);
     // h_avg_btag_disc_non_btags[this_category][iSys]->Fill(avg_btag_disc_non_btags,wgt);
      h_second_jet_pt[this_category][iSys]->Fill(second_jet_pt,wgt);
    //  h_dr_between_lep_and_closest_jet[this_category][iSys]->Fill(dr_between_lep_and_closest_jet,wgt);
      h_h0[this_category][iSys]->Fill(h0,wgt);
      h_avg_btag_disc_btags[this_category][iSys]->Fill(avg_btag_disc_btags,wgt);
      h_dev_from_avg_disc_btags[this_category][iSys]->Fill(dev_from_avg_disc_btags,wgt);
      h_third_jet_pt[this_category][iSys]->Fill(third_jet_pt,wgt);
      h_fourth_jet_pt[this_category][iSys]->Fill(fourth_jet_pt,wgt);
     // h_avg_tagged_dijet_mass[this_category][iSys]->Fill(avg_tagged_dijet_mass,wgt);
      h_h1[this_category][iSys]->Fill(h1,wgt);
      h_h2[this_category][iSys]->Fill(h2,wgt);
      h_all_sum_pt_with_met[this_category][iSys]->Fill(all_sum_pt_with_met,wgt);
      h_closest_tagged_dijet_mass[this_category][iSys]->Fill(closest_tagged_dijet_mass,wgt);
      h_h3[this_category][iSys]->Fill(h3,wgt);
    //  h_h4[this_category][iSys]->Fill(h4,wgt);
    //  h_first_highest_btag[this_category][iSys]->Fill(first_highest_btag,wgt);
      h_second_highest_btag[this_category][iSys]->Fill(second_highest_btag,wgt);
      h_third_highest_btag[this_category][iSys]->Fill(third_highest_btag,wgt);
      h_fourth_highest_btag[this_category][iSys]->Fill(fourth_highest_btag,wgt);
      h_fifth_highest_CSV[this_category][iSys]->Fill(fifth_highest_CSV,wgt);
     // h_sixth_highest_CSV[this_category][iSys]->Fill(sixth_highest_CSV,wgt);
      h_invariant_mass_of_everything[this_category][iSys]->Fill(invariant_mass_of_everything,wgt);
      h_Mlb[this_category][iSys]->Fill(Mlb,wgt);
      h_tagged_dijet_mass_closest_to_125[this_category][iSys]->Fill(tagged_dijet_mass_closest_to_125,wgt);
      h_maxeta_jet_jet[this_category][iSys]->Fill(maxeta_jet_jet,wgt);		 
      h_maxeta_jet_tag[this_category][iSys]->Fill(maxeta_jet_tag,wgt);		 
      h_maxeta_tag_tag[this_category][iSys]->Fill(maxeta_tag_tag,wgt);			 
      h_dEta_fn[this_category][iSys]->Fill(dEta_fn,wgt);			 
      h_pt_all_jets_over_E_all_jets[this_category][iSys]->Fill(pt_all_jets_over_E_all_jets,wgt);  
      h_HT[this_category][iSys]->Fill(HT,wgt);
     
	//if(iSys==0)		cout<<run<<" "<<lumi<<" "<<evt<<" "<<aplanarity_f<<" "<<sphericity_f<<" "<<tight_lepton_pt<<" "<<MET<<endl;
	
	
	  //Fill KIT histos
	  //KIT uses odd events on training and even events on testing and evaluating limits
	  if(int(evt)%2==0){
	  if(this_category==2){
	    	h1_6j2t_KIT[iSys]->Fill(h1,wgthalf);
			avg_dr_tagged_jets_6j2t_KIT[iSys]->Fill(avg_dr_tagged_jets,wgthalf);
			sphericity_6j2t_KIT[iSys]->Fill(sphericity,wgthalf);
			third_highest_btag_6j2t_KIT[iSys]->Fill(third_highest_btag,wgthalf);
			h3_6j2t_KIT[iSys]->Fill(h3,wgthalf);
			HT_6j2t_KIT[iSys]->Fill(HT,wgthalf);
			Mlb_6j2t_KIT[iSys]->Fill(Mlb,wgthalf);
			fifth_highest_CSV_6j2t_KIT[iSys]->Fill(fifth_highest_CSV,wgthalf);
			fourth_highest_btag_6j2t_KIT[iSys]->Fill(fourth_highest_btag,wgthalf);
		}
		if(this_category==3){
		    h1_4j3t_KIT[iSys]->Fill(h1,wgthalf);
			avg_dr_tagged_jets_4j3t_KIT[iSys]->Fill(avg_dr_tagged_jets,wgthalf);
			sphericity_4j3t_KIT[iSys]->Fill(sphericity,wgthalf);
			third_highest_btag_4j3t_KIT[iSys]->Fill(third_highest_btag,wgthalf);
			HT_4j3t_KIT[iSys]->Fill(HT,wgthalf);
			dev_from_avg_disc_btags_4j3t_KIT[iSys]->Fill(dev_from_avg_disc_btags,wgthalf);
			M3_4j3t_KIT[iSys]->Fill(M3,wgthalf);
			min_dr_tagged_jets_4j3t_KIT[iSys]->Fill(min_dr_tagged_jets,wgthalf);
			Evt_CSV_Average_4j3t_KIT[iSys]->Fill(CSV_Average,wgthalf);
		}
		
		if(this_category==4){
		
			h1_5j3t_KIT[iSys]->Fill(h1,wgthalf);
			avg_dr_tagged_jets_5j3t_KIT[iSys]->Fill(avg_dr_tagged_jets,wgthalf);
			sphericity_5j3t_KIT[iSys]->Fill(sphericity,wgthalf);
			third_highest_btag_5j3t_KIT[iSys]->Fill(third_highest_btag,wgthalf);
			h3_5j3t_KIT[iSys]->Fill(h3,wgthalf);
			HT_5j3t_KIT[iSys]->Fill(HT,wgthalf);
			dev_from_avg_disc_btags_5j3t_KIT[iSys]->Fill(dev_from_avg_disc_btags,wgthalf);
			fourth_highest_btag_5j3t_KIT[iSys]->Fill(fourth_highest_btag,wgthalf);
		}
			if(this_category==5){
			h1_6j3t_KIT[iSys]->Fill(h1,wgthalf);
			avg_dr_tagged_jets_6j3t_KIT[iSys]->Fill(avg_dr_tagged_jets,wgthalf);
			third_highest_btag_6j3t_KIT[iSys]->Fill(third_highest_btag,wgthalf);
			HT_6j3t_KIT[iSys]->Fill(HT,wgthalf);
			pt_all_jets_over_E_all_jets_6j3t_KIT[iSys]->Fill(pt_all_jets_over_E_all_jets,wgthalf);
			fifth_highest_CSV_6j3t_KIT[iSys]->Fill(fifth_highest_CSV,wgthalf);
			fourth_highest_btag_6j3t_KIT[iSys]->Fill(fourth_highest_btag,wgthalf);
			min_dr_tagged_jets_6j3t_KIT[iSys]->Fill(min_dr_tagged_jets,wgthalf);
		}
		
		if(this_category==6){	
			h1_4j4t_KIT[iSys]->Fill(h1,wgthalf);
			avg_dr_tagged_jets_4j4t_KIT[iSys]->Fill(avg_dr_tagged_jets,wgthalf);
			sphericity_4j4t_KIT[iSys]->Fill(sphericity,wgthalf);
			avg_btag_disc_btags_4j4t_KIT[iSys]->Fill(avg_btag_disc_btags,wgthalf);
			h2_4j4t_KIT[iSys]->Fill(h2,wgthalf);
			closest_tagged_dijet_mass_4j4t_KIT[iSys]->Fill(closest_tagged_dijet_mass,wgthalf);
			second_highest_btag_4j4t_KIT[iSys]->Fill(second_highest_btag,wgthalf);
			fourth_highest_btag_4j4t_KIT[iSys]->Fill(fourth_highest_btag,wgthalf);
			maxeta_jet_jet_4j4t_KIT[iSys]->Fill(maxeta_jet_jet,wgthalf);
			pt_all_jets_over_E_all_jets_4j4t_KIT[iSys]->Fill(pt_all_jets_over_E_all_jets,wgthalf);
		}
		
		if(this_category==7){
			avg_dr_tagged_jets_5j4t_KIT[iSys]->Fill(avg_dr_tagged_jets,wgthalf);
			HT_5j4t_KIT[iSys]->Fill(HT,wgthalf);
			M3_5j4t_KIT[iSys]->Fill(M3,wgthalf);
			fifth_highest_CSV_5j4t_KIT[iSys]->Fill(fifth_highest_CSV,wgthalf);
			fourth_highest_btag_5j4t_KIT[iSys]->Fill(fourth_highest_btag,wgthalf);
			Evt_Deta_JetsAverage_5j4t_KIT[iSys]->Fill(Deta_JetsAverage,wgthalf);
			avg_btag_disc_btags_5j4t_KIT[iSys]->Fill(avg_btag_disc_btags,wgthalf);
		}
			if(this_category==8){
				h2_6j4t_KIT[iSys]->Fill(h2,wgthalf);
				avg_dr_tagged_jets_6j4t_KIT[iSys]->Fill(avg_dr_tagged_jets,wgthalf);
				third_highest_btag_6j4t_KIT[iSys]->Fill(third_highest_btag,wgthalf);
				M3_6j4t_KIT[iSys]->Fill(M3,wgthalf);
				fifth_highest_CSV_6j4t_KIT[iSys]->Fill(fifth_highest_CSV,wgthalf);
				fourth_highest_btag_6j4t_KIT[iSys]->Fill(fourth_highest_btag,wgthalf);
				best_higgs_mass_6j4t_KIT[iSys]->Fill(best_higgs_mass,wgthalf);
				dEta_fn_6j4t_KIT[iSys]->Fill(dEta_fn,wgthalf);
		    }
			
		}//end KIT evenevts
		}//end loop over systematics
	
		
		
	
  	} // end loop over evts
	
	
	
	
	
	double RawCutEvts[7];
	double CutEvts[7];
	double CatEvts[9];
	double RawCatEvts[9];
	for(int i=0;i<9;i++){
	if(i<7)RawCutEvts[i] = h_NumCountRaw->GetBinContent(i+2);
	//if(i<7)CutEvts[i] = h_NumCount->GetBinContent(i+1);
	if(i>1)CatEvts[i] = h_BDT_Output_IncludeWgt[i][0]->Integral();
	if(i>1)RawCatEvts[i] = h_BDT_Output_NoWgt[i][0]->Integral();
    }
    
	std::string t_end = "_" + str_jobN + "_Info.txt";
    if( Njobs==1 ) t_end = "_Info.txt";
	std::string textfilename = "/eos/uscms/store/user/sflowers/TreeReader/Spring15_Sep15Trees/Test/" + sampleName + t_end;
	
	ofstream textfile;
	  textfile.open (textfilename);
	  
	  //textfile << "Cut Evts   Raw : Weighted"<<endl;
	  textfile<<"Cut Evts Raw"<<endl;
	  for(int i=0;i<7;i++){
		  //textfile<<"Cut #"<<i<<" "<<RawCutEvts[i]<<"   :   "<<CutEvts<<endl;
		  textfile<<"Cut #"<<i<<" "<<RawCutEvts[i]<<endl;
	  }
	  textfile<<endl;
	  textfile<< "Category Yields    Raw : Weighted"<<endl;
	  for(int i=2;i<9;i++){
		  textfile<<cat_labels[i]<<" "<<RawCatEvts[i]<<"     :     "<<CatEvts[i]<<endl;
	  }
	 
	 
	 
	  textfile.close();
		
   	 t=clock()-t;
	
	
		std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;
	
		
    
	

    histofile.Write();
    histofile.Close();

	

	 

}






void fillCSVhistos(TFile* fileHF, TFile* fileLF){

  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c = "final";
    TString syst_csv_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<6; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<6; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

  return;
}

double get_csv_wgt( vvdouble jets, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){

  int iSysHF = 0;
  switch(iSys){
  case 7: iSysHF=1; break;
  case 8: iSysHF=2; break;
  case 9: iSysHF=3; break;
  case 10: iSysHF=4; break;
  case 13: iSysHF=5; break;
  case 14: iSysHF=6; break;
  case 15: iSysHF=7; break;
  case 16: iSysHF=8; break;
  default : iSysHF = 0; break;
  }

  int iSysC = 0;
  switch(iSys){
  case 21: iSysC=1; break;
  case 22: iSysC=2; break;
  case 23: iSysC=3; break;
  case 24: iSysC=4; break;
  default : iSysC = 0; break;
  }

  int iSysLF = 0;
  switch(iSys){
  case 7: iSysLF=1; break;
  case 8: iSysLF=2; break;
  case 11: iSysLF=3; break;
  case 12: iSysLF=4; break;
  case 17: iSysLF=5; break;
  case 18: iSysLF=6; break;
  case 19: iSysLF=7; break;
  case 20: iSysLF=8; break;
  default : iSysLF = 0; break;
  }

  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;


  for( int iJet=0; iJet<int(jets.size()); iJet++ ){
    TLorentzVector myJet;
    myJet.SetPxPyPzE( jets[iJet][0], jets[iJet][1], jets[iJet][2], jets[iJet][3] );
	  
    double csv = jetCSV[iJet];
    double jetPt = myJet.Pt();
    double jetAbsEta = fabs( myJet.Eta() );
    int flavor = jetFlavor[iJet];

    int iPt = -1; int iEta = -1;
    if (jetPt >=19.99 && jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100 && jetPt<160) iPt = 4;
    else if (jetPt >=160 && jetPt<10000) iPt = 5;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

    // tentative solution...
    if (iPt < 0 || iEta < 0){
      if( iPt < 0 ){ iPt = 0  ;} 
      if( iEta < 0 ){ iEta = 0;} 
    }

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;

      // if( iSysHF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			     flavor, jetPt, iJet->eta, csv, iCSVWgtHF );
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
      // if( iSysC==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      //      flavor, jetPt, iJet->eta, csv, iCSVWgtC );
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;

      // if( iSysLF==0 ) printf(" iJet,\t flavor=%d,\t pt=%.1f,\t eta=%.2f,\t csv=%.3f,\t wgt=%.2f \n",
      // 			     flavor, jetPt, iJet->eta, csv, iCSVWgtLF );
    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  return csvWgtTotal;
}


