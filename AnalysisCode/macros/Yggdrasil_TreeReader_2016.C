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
#include "ttH-LeptonPlusJets/AnalysisCode/interface/YggHelper.h"
#include "TTH/CommonClassifier/interface/BDTClassifier.h"
//#include "TTH/CommonClassifier/interface/MEMClassifier.h"
# include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilEventSelection.h"



#endif




typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;




void Yggdrasil_TreeReader_2016(  int SampleType=2,int SplitType =2, int usedSlimmedTrees=0, int maxNentries=-1, int Njobs=1, int jobN=1) {

clock_t t;
t = clock();

//Setup Input and Output  

int insample = 0;
  
  std::cout << "===> Loading root files" << std::endl;

  std::string sampleName = "TTJets";
  if(SampleType == 2){
	if(SplitType == 1) sampleName = "ttb";
	if(SplitType == 2) sampleName = "ttbb";
	if(SplitType == 3) sampleName = "tt2b";
	if(SplitType == 4) sampleName = "ttcc";
	if(SplitType == 5) sampleName = "ttlf";
  }
  if(SampleType == 1){
	if(SplitType == 0) sampleName = "ttH_Nonbb";
	if(SplitType == 1) sampleName = "ttH_hbb";
	if(SplitType == 2) sampleName = "ttH_hww";
	if(SplitType == 3) sampleName = "ttH_htt";
	if(SplitType == 4) sampleName = "ttH_hgluglu";
	if(SplitType == 5) sampleName = "ttH_hcc";
	if(SplitType == 6) sampleName = "ttH_hzz";
	if(SplitType == 7) sampleName = "ttH_hzg";
	if(SplitType == 8) sampleName = "ttH_hgg";
	if(SplitType == 9) sampleName = "ttH_hss";
	if(SplitType == 10) sampleName = "ttH_hmumu";
  } 
  if(SampleType == 3) sampleName = "ttHTobb";
  if(SampleType == 6) sampleName = "SingleElectronB";
  if(SampleType == 7) sampleName = "SingleElectronC";
  if(SampleType == 8) sampleName = "SingleElectronD";
  if(SampleType == 9) sampleName = "SingleMuonB";
  if(SampleType == 10) sampleName = "SingleMuonC";
  if(SampleType == 11) sampleName = "SingleMuonD";

 
  std::string str_jobN;
  std::stringstream stream;
  stream << jobN;
  str_jobN = stream.str();
 
 
  std::string treefilename = "";
  std::string treefilename2 = "";
  std::string treefilename3 = "";
  
  if(SampleType == 2 && usedSlimmedTrees==0){
  	 treefilename = "/eos/uscms/store/user/lpctthrun2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/ttbar_Satoshi_Yggdra_Aug01/160801_165115/0000/*.root";
	 treefilename2 = "/eos/uscms/store/user/lpctthrun2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/ttbar_Satoshi_Yggdra_Aug01/160801_165115/0001/*.root";
	 treefilename3 = "/eos/uscms/store/user/lpctthrun2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/ttbar_Satoshi_Yggdra_Aug01/160801_165115/0002/*.root";
  }
  
  if(usedSlimmedTrees == 1){
  	if(SampleType == 1)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_ttHToNonbb_*.root";
  	if(SampleType == 2)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_ttJets_*.root";
	if(SampleType == 3)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_ttHTobb_*.root";
	if(SampleType == 6)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_SingleElectronB_*.root";
	if(SampleType == 7)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_SingleElectronC_*.root";
	if(SampleType == 8)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_SingleElectronD_*.root";
	if(SampleType == 9)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_SingleMuonB_*.root";
	if(SampleType == 10)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_SingleMuonC_*.root";
	if(SampleType == 11)treefilename = "/eos/uscms/store/user/sflowers/SlimTrees/SlimTrees_Aug5th_SingleMuonD_*.root";
  }
  
  std::string OutputDir ="/eos/uscms/store/user/sflowers/TreeReader/YT2016_Histos_Aug9/";
 
  std::string s_end = "_" + str_jobN + ".root";
  if( Njobs==1 ) s_end = ".root";
  std::string histofilename = OutputDir + "YggdrasilHistoFile_" + sampleName +s_end;
  
  std::cout << "  treefilename  = " << treefilename.c_str() << std::endl;
  if(SampleType ==2 && usedSlimmedTrees==0)std::cout<< "  +  "<<treefilename2.c_str()<<std::endl;
  if(SampleType ==2 && usedSlimmedTrees==0)std::cout<< "  +  "<<treefilename3.c_str()<<std::endl;
  
  std::cout << "  histofilename = " << histofilename.c_str() << std::endl;

  TChain *chain = new TChain("ttHTreeMaker/worldTree");
  chain->Add(treefilename.c_str());
  if(SampleType ==2 && usedSlimmedTrees ==0)chain->Add(treefilename3.c_str());   // for longer samples
  if(SampleType ==2 && usedSlimmedTrees ==0)chain->Add(treefilename2.c_str());
 
  yggdrasilEventVars *eve=0;
  chain->SetBranchAddress("eve.", &eve );
  
  

  
  ////////////////////////////////////

	 
	 
	cout<<"===> Loading classifier information"<<endl;
	 

	BDTClassifier kitBDTclass((string(getenv("CMSSW_BASE"))+"/src/TTH/CommonClassifier/data/bdtweights_v5/"));
	YggHelper yGGhelper;
	
	
	 
  	TFile* f_CSVwgt_HF = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/csv_rwt_fit_hf_76x_2016_02_08.root").c_str());
 	TFile* f_CSVwgt_LF = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/csv_rwt_fit_lf_76x_2016_02_08.root").c_str());
 	
	TFile* f_LEPwgt_Ele_Trig = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/eleTrig_SF.root").c_str());
 	TFile* f_LEPwgt_Ele_Iso  = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/eleIsolation_SF.root").c_str());
 	TFile* f_LEPwgt_Ele_ID   = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/ScaleFactor_GsfElectronToRECO_passingTrigWP80.txt.egamma_SF2D.root").c_str());

	TFile* f_LEPwgt_Mu_Trig = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root").c_str());
    TFile* f_LEPwgt_Mu_Iso  = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/MuonID_Z_RunCD_Reco76X_Feb15.root").c_str());
    TFile* f_LEPwgt_Mu_ID   = new TFile ((string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/AnalysisCode/data/MuonIso_Z_RunCD_Reco76X_Feb15.root").c_str());


 	yGGhelper.fillCSVhistos(f_CSVwgt_HF, f_CSVwgt_LF);	 
	yGGhelper.FillLeptonWgtHistos(f_LEPwgt_Ele_Trig, f_LEPwgt_Ele_ID, f_LEPwgt_Ele_Iso, f_LEPwgt_Mu_Trig, f_LEPwgt_Mu_ID,f_LEPwgt_Mu_Iso);
	 

    //////////////////////////////////////////////////////////////////////////
    ///  Histograms
    //////////////////////////////////////////////////////////////////////////
  
  
     TFile histofile(histofilename.c_str(),"recreate");
     histofile.cd();
	 
 
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
	 /*
	 sys_cat_labels.push_back("_CSVCErr1Up");     //21
     sys_cat_labels.push_back("_CSVCErr1Down");   //22
     sys_cat_labels.push_back("_CSVCErr2Up");     //23
     sys_cat_labels.push_back("_CSVCErr2Down");   //24
     sys_cat_labels.push_back("_Q2scaleUp");     //25
     sys_cat_labels.push_back("_Q2scaleDown");   //26
     sys_cat_labels.push_back("_topPtcorrUp");   //27
     sys_cat_labels.push_back("_topPtcorrDown"); //28
	*/
   	 int NumSysCat = int(sys_cat_labels.size());

	 
	 TH1::SetDefaultSumw2();
	 
	 

	 //Event Information//
	 TH1D* h_additionalJetEventId;
	 
	 TH1D* h_numJets[NumSysCat];
	 TH1D* h_numTags[NumSysCat];
	 TH1D* h_category_yield[NumSysCat];
	 TH1D* h_category_yield_NoWgt[NumSysCat];
	 TH1D* h_CutFlow[NumSysCat];
	 TH1D* h_CutFlow_NoWgt[NumSysCat];
	 TH1D* h_CutFlow_NoCSV[NumSysCat];
	 
	 TH1D* h_CutFlow_NoWgt_KIT[NumSysCat];
	 
	 //ttH, Higgs specific//
	 
	 TH1D* h_higgsDecayType[NumCat][NumSysCat];
	 TH2D* h_category_yield_higgsDecayType[NumSysCat];
	 
	 //Lepton Plots//
	 
	 TH1D* h_tight_lepton_pt[NumCat];
	 TH1D* h_tight_lepton_phi[NumCat];
	 TH1D* h_tight_lepton_eta[NumCat];
	 
	 //Jet Plots//
	 
	 TH1D* h_jet_pt[NumCat][NumSysCat];
	 TH1D* h_jet_pt_NoWgt[NumCat][NumSysCat];
	 TH1D* h_jet_eta[NumCat][NumSysCat];
	 TH1D* h_jet_eta_NoWgt[NumCat][NumSysCat];
	 TH1D* h_jet_csv[NumCat][NumSysCat];
	 TH1D* h_jet_csv_NoWgt[NumCat][NumSysCat];
	 
	 TH1D* h_jet_tag_pt[NumCat][NumSysCat];
	 TH1D* h_jet_tag_pt_NoWgt[NumCat][NumSysCat];
	 TH1D* h_jet_tag_eta[NumCat][NumSysCat];
	 TH1D* h_jet_tag_eta_NoWgt[NumCat][NumSysCat];
	 TH1D* h_jet_tag_csv[NumCat][NumSysCat];
	 TH1D* h_jet_tag_csv_NoWgt[NumCat][NumSysCat];
	 
	 //Individual Jet Plots//
	 
	 TH1D* h_first_jet_pt[NumCat][NumSysCat];
	 TH1D* h_second_jet_pt[NumCat][NumSysCat];
	 TH1D* h_third_jet_pt[NumCat][NumSysCat];
	 TH1D* h_fourth_jet_pt[NumCat][NumSysCat];
	 
	 TH1D* h_first_jet_eta[NumCat][NumSysCat];
	 TH1D* h_second_jet_eta[NumCat][NumSysCat];
	 TH1D* h_third_jet_eta[NumCat][NumSysCat];
	 TH1D* h_fourth_jet_eta[NumCat][NumSysCat];
	 
	 TH1D* h_first_jet_csv[NumCat][NumSysCat];
	 TH1D* h_second_jet_csv[NumCat][NumSysCat];
	 TH1D* h_third_jet_csv[NumCat][NumSysCat];
	 TH1D* h_fourth_jet_csv[NumCat][NumSysCat];
	 
	 TH1D* h_first_jet_pt_NoWgt[NumCat][NumSysCat];
	 TH1D* h_first_jet_csv_NoWgt[NumCat][NumSysCat];
	 
	 TH1D* h_first_highest_btag[NumCat][NumSysCat];
	 TH1D* h_second_highest_btag[NumCat][NumSysCat];
	 TH1D* h_third_highest_btag[NumCat][NumSysCat];
	 TH1D* h_fourth_highest_btag[NumCat][NumSysCat];
	 TH1D* h_fifth_highest_CSV[NumCat][NumSysCat];
	 TH1D* h_sixth_highest_CSV[NumCat][NumSysCat];
	 
	 TH1D* h_lowest_btag[NumCat][NumSysCat];
	 
	 //Multi Object Plots//
	 
	 TH1D* h_CSV_Average[NumCat][NumSysCat];
	 TH1D* h_Deta_JetsAverage[NumCat][NumSysCat];
	 TH1D* h_min_dr_tagged_jets[NumCat][NumSysCat];
	 TH1D* h_avg_dr_tagged_jets[NumCat][NumSysCat];
	 TH1D* h_dr_between_lep_and_closest_jet[NumCat][NumSysCat];
	 TH1D* h_avg_btag_disc_btags[NumCat][NumSysCat];
	 TH1D* h_dev_from_avg_disc_btags[NumCat][NumSysCat];
	 
	 TH1D* h_HT[NumCat][NumSysCat];
	 TH1D* h_all_sum_pt_with_met[NumCat][NumSysCat];
	 TH1D* h_MET_pt[NumCat][NumSysCat];
	 TH1D* h_MHT[NumCat][NumSysCat];
	 TH1D* h_invariant_mass_of_everything[NumCat][NumSysCat];
	 TH1D* h_avg_untagged_dijet_mass[NumCat][NumSysCat];;
	 TH1D* h_closest_tagged_dijet_mass[NumCat][NumSysCat];;
	 TH1D* h_M3[NumCat][NumSysCat];
	 TH1D* h_Mlb[NumCat][NumSysCat];
	 TH1D* h_tagged_dijet_mass_closest_to_125[NumCat][NumSysCat];
	 
	 //BDTvars variables//
	 
	 TH1D* h_aplanarity[NumCat][NumSysCat];
	 TH1D* h_sphericity[NumCat][NumSysCat];
	
	 TH1D* h_h0[NumCat][NumSysCat];
	 TH1D* h_h1[NumCat][NumSysCat];
	 TH1D* h_h2[NumCat][NumSysCat];
	 TH1D* h_h3[NumCat][NumSysCat];
	 TH1D* h_h4[NumCat][NumSysCat];
	 
	 TH1D* h_maxeta_jet_jet[NumCat][NumSysCat];
	 TH1D* h_maxeta_jet_tag[NumCat][NumSysCat];
	 TH1D* h_maxeta_tag_tag[NumCat][NumSysCat];
	 
	 TH1D* h_best_higgs_mass[NumCat][NumSysCat];
	 TH1D* h_minChi2[NumCat][NumSysCat];
	 TH1D* h_dRbb[NumCat][NumSysCat];
	 TH1D* h_abs_dEta_leptop_bb[NumCat][NumSysCat];
	 TH1D* h_abs_dEta_hadtop_bb[NumCat][NumSysCat];
	 TH1D* h_dEta_fn[NumCat][NumSysCat];
	 TH1D* h_median_bb_mass[NumCat][NumSysCat];
	 TH1D* h_pt_all_jets_over_E_all_jets[NumCat][NumSysCat];
	 
	 //BDT Histos//
	 
	 TH1D* h_bdtOutput[NumCat][NumSysCat];
	
	 //Histogram Parameters

     double NpvMin = 0-0.5;
     double NpvMax = 50+0.5;
     int NpvBins = int( NpvMax - NpvMin + 0.0001 );

     int NjetMin_full = 4;
     int NjetMax_full = 10;
     int NjetBins_full = NjetMax_full - NjetMin_full + 1;

     int NjetMin = 4;
     int NjetMax = 6;
     int NjetBins = NjetMax - NjetMin + 1;
     int NtagMin = 2;
     int NtagMax = 4;
     int NtagBins = NtagMax - NtagMin + 1;

     int NtagMin_full = 2;
     int NtagMax_full = 7;
     int NtagBins_full = NtagMax_full - NtagMin_full + 1;


     int NjetMin_loose = 0;
     int NjetMax_loose = 5;
     int NjetBins_loose = NjetMax_loose - NjetMin_loose + 1;

     int NtagMin_loose = 0;
     int NtagMax_loose = 2;
     int NtagBins_loose = NtagMax_loose - NtagMin_loose + 1;

     double metmax   = 500.;
     double lepPtMax = 300.;
     double jetptmax = 500.;
     double htmax    = 2000.;
     double massmax  = 500.;
     double m3max    = 1000.;
     double massmax_sumTop = 2000.;

     int NmetBins   = int( metmax/20. + 0.0001 );
     int NlepPtBins = int( lepPtMax/10. + 0.0001 );
     int NjetptBins = int( jetptmax/10. + 0.0001 );
     int NhtBins    = int( htmax/50. + 0.0001 );
     int NmassBins  = int( massmax/10. + 0.0001 );
     int Nm3Bins    = int( m3max/10. + 0.0001 );
     int NmassBins_sumTop  = int( massmax/10. + 0.0001 );

     int NcsvBins = 44;
	 
	 h_additionalJetEventId = new TH1D("h_additionalJetEventId","h_additionalJetEventId",100,0,100);
	 
	 for( int b=0; b<NumSysCat; b++ ){
		 std::string short_suffix = sys_cat_labels[b];
		 
	     h_numJets[b] = new TH1D((std::string("h_numJets" + short_suffix)).c_str(), "Number of Jets", NjetBins_full, NjetMin_full-0.5, NjetMax_full+0.5 );
	     h_numTags[b] = new TH1D((std::string("h_numTags" + short_suffix)).c_str(), "Number of Tags", NtagBins_full, NtagMin_full-0.5, NtagMax_full+0.5 );
	     h_category_yield[b] = new TH1D((std::string("h_category_yield" + short_suffix)).c_str(), "Category Yield", NumCat, 0, NumCat );
	     h_category_yield_NoWgt[b] = new TH1D((std::string("h_category_yield_NoWgt" + short_suffix)).c_str(), "Category Yield (unweighted)", NumCat, 0, NumCat );
	     h_category_yield_higgsDecayType[b] = new TH2D((std::string("h_category_yield_higgsDecayType" + short_suffix)).c_str(), "Higgs decay", NumCat, 0, NumCat, 11, 0, 11 );
		 
   
         h_CutFlow[b] = new TH1D((std::string("h_CutFlow" + short_suffix)).c_str(),"CutFlow",10,0,10);
		 h_CutFlow_NoCSV[b] = new TH1D((std::string("h_CutFlow_NoCSV" + short_suffix)).c_str(),"CutFlow (no csv weight)",10,0,10);
		 h_CutFlow_NoWgt[b] = new TH1D((std::string("h_CutFlow_NoWgt" + short_suffix)).c_str(),"CutFlow (unweighted)",10,0,10);
		 
		 h_CutFlow_NoWgt_KIT[b] = new TH1D((std::string("h_CutFlow_NoWgt_KIT" + short_suffix)).c_str(),"CutFlow (unweighted)",10,0,10);
		 
		
         h_tight_lepton_pt[b] = new TH1D((std::string("h_tight_lepton_pt" + short_suffix)).c_str(),"Tight Lepton p_{T}", NlepPtBins, 0, lepPtMax );
         h_tight_lepton_phi[b] = new TH1D((std::string("h_tight_lepton_phi" + short_suffix)).c_str(),"Tight Lepton #phi", 34, -3.4, 3.4 );
         h_tight_lepton_eta[b] = new TH1D((std::string("h_tight_lepton_eta" + short_suffix)).c_str(),"Tight Lepton #eta", 25, -2.5, 2.5 );
		 
		 
		 for( int c=0; c<NumCat; c++ ){ 
        	std::string suffix = "_" + cat_labels[c];
         	suffix += sys_cat_labels[b];
			
   		    h_higgsDecayType[c][b] = new TH1D((std::string("h_higgsDecayType" + suffix)).c_str(), ";Higgs decay", 11, 0, 11 );
	     
			
            if(b==0)h_tight_lepton_pt[c] = new TH1D((std::string("h_tight_lepton_pt" + suffix)).c_str(),"Tight Lepton p_{T}", NlepPtBins, 0, lepPtMax );
            if(b==0)h_tight_lepton_phi[c] = new TH1D((std::string("h_tight_lepton_phi" + suffix)).c_str(),"Tight Lepton #phi", 34, -3.4, 3.4 );
            if(b==0)h_tight_lepton_eta[c] = new TH1D((std::string("h_tight_lepton_eta" + suffix)).c_str(),"Tight Lepton #eta", 25, -2.5, 2.5 );
  		 
		 	h_jet_pt[c][b] = new TH1D((std::string("h_jet_pt" + suffix)).c_str(),"Jet p_{T}", NjetptBins, 0, jetptmax );
		 	h_jet_pt_NoWgt[c][b] = new TH1D((std::string("h_jet_pt_NoWgt" + suffix)).c_str(),"Jet p_{T} (unweighted)", NjetptBins, 0, jetptmax );
			h_jet_eta[c][b] = new TH1D((std::string("h_jet_eta" + suffix)).c_str(),"Jet #eta", 25, -2.5, 2.5 );
			h_jet_eta_NoWgt[c][b] = new TH1D((std::string("h_jet_eta_NoWgt" + suffix)).c_str(),"Jet #eta (unweighted)", 25, -2.5, 2.5 );
			h_jet_csv[c][b] = new TH1D((std::string("h_jet_csv" + suffix)).c_str(),"Jet CSV", NcsvBins, -1.1, 1.1 );
			h_jet_csv_NoWgt[c][b] = new TH1D((std::string("h_jet_csv_NoWgt" + suffix)).c_str(),"Jet CSV (unweighted)", NcsvBins, -1.1, 1.1 );
		 
		 	h_jet_tag_pt[c][b] = new TH1D((std::string("h_jet_tag_pt" + suffix)).c_str(),"Tagged Jet p_{T}", NjetptBins, 0, jetptmax );
		 	h_jet_tag_pt_NoWgt[c][b] = new TH1D((std::string("h_jet_tag_pt_NoWgt" + suffix)).c_str(),"Tagged Jet p_{T} (unweighted)", NjetptBins, 0, jetptmax );
		 	h_jet_tag_eta[c][b] = new TH1D((std::string("h_jet_tag_eta" + suffix)).c_str(),"Tagged Jet #eta", 25, -2.5, 2.5 );
		 	h_jet_tag_eta_NoWgt[c][b] = new TH1D((std::string("h_jet_tag_eta_NoWgt" + suffix)).c_str(),"Tagged Jet #eta (unweighted)", 25, -2.5, 2.5 );
		 	h_jet_tag_csv[c][b] = new TH1D((std::string("h_jet_tag_csv" + suffix)).c_str(),"Tagged Jet CSV", NcsvBins, -1.1, 1.1 );
		 	h_jet_tag_csv_NoWgt[c][b] = new TH1D((std::string("h_jet_tag_csv_NoWgt" + suffix)).c_str(),"Tagged Jet CSV (unweighted)", NcsvBins, -1.1, 1.1 );
		 
		 	h_first_jet_pt[c][b] = new TH1D((std::string("h_first_jet_pt" + suffix)).c_str(),"First Jet p_{T}", NjetptBins, 0, jetptmax );
		 	h_second_jet_pt[c][b] = new TH1D((std::string("h_second_jet_pt" + suffix)).c_str(),"Second Jet p_{T}", NjetptBins, 0, jetptmax );
		 	h_third_jet_pt[c][b] = new TH1D((std::string("h_third_jet_pt" + suffix)).c_str(),"Third Jet p_{T}", NjetptBins, 0, jetptmax );
		 	h_fourth_jet_pt[c][b] = new TH1D((std::string("h_fourth_jet_pt" + suffix)).c_str(),"Fourth Jet p_{T}", NjetptBins, 0, jetptmax );
		 
		 	h_first_jet_eta[c][b] = new TH1D((std::string("h_first_jet_eta" + suffix)).c_str(),"First Jet #eta", 25, -2.5, 2.5 );
		 	h_second_jet_eta[c][b] = new TH1D((std::string("h_second_jet_eta" + suffix)).c_str(),"Second Jet #eta", 25, -2.5, 2.5 );
		 	h_third_jet_eta[c][b] = new TH1D((std::string("h_third_jet_eta" + suffix)).c_str(),"Third Jet #eta", 25, -2.5, 2.5 );
		 	h_fourth_jet_eta[c][b] = new TH1D((std::string("h_fourth_jet_eta" + suffix)).c_str(),"Fourth Jet #eta", 25, -2.5, 2.5 );
		
		 	h_first_jet_csv[c][b] = new TH1D((std::string("h_first_jet_csv" + suffix)).c_str(),"First Jet CSV", NcsvBins, -1.1, 1.1 );
		 	h_second_jet_csv[c][b] = new TH1D((std::string("h_second_jet_csv" + suffix)).c_str(),"Second Jet CSV", NcsvBins, -1.1, 1.1 );
		 	h_third_jet_csv[c][b] = new TH1D((std::string("h_third_jet_csv" + suffix)).c_str(),"Third Jet CSV", NcsvBins, -1.1, 1.1 );
		 	h_fourth_jet_csv[c][b] = new TH1D((std::string("h_fourth_jet_csv" + suffix)).c_str(),"Fourth Jet CSV", NcsvBins, -1.1, 1.1 ); 
		 
		 	h_first_jet_pt_NoWgt[c][b] = new TH1D((std::string("h_first_jet_pt_NoWgt" + suffix)).c_str(),"First Jet p_{T} (unweighted)", NjetptBins, 0, jetptmax );
		 	h_first_jet_csv_NoWgt[c][b] = new TH1D((std::string("h_first_jet_csv_NoWgt" + suffix)).c_str(),"First Jet CSV (unweighted)", NcsvBins, -1.1, 1.1 );
		
         	h_first_highest_btag[c][b] = new TH1D((std::string("h_first_highest_btag" + suffix)).c_str(),"First Highest CSV", NcsvBins, -1.1, 1.1 );
         	h_second_highest_btag[c][b] = new TH1D((std::string("h_second_highest_btag" + suffix)).c_str(),"Second Highest CSV", NcsvBins, -1.1, 1.1 );
         	h_third_highest_btag[c][b] = new TH1D((std::string("h_third_highest_btag" + suffix)).c_str(),"Third Highest CSV", NcsvBins, -1.1, 1.1 );
         	h_fourth_highest_btag[c][b] = new TH1D((std::string("h_fourth_highest_btag" + suffix)).c_str(),"Fourth Highest CSV", NcsvBins, -1.1, 1.1 );
         	h_fifth_highest_CSV[c][b] = new TH1D((std::string("h_fifth_highest_CSV" + suffix)).c_str(),"Fifth Highest CSV", NcsvBins, -1.1, 1.1 );
         	h_sixth_highest_CSV[c][b] = new TH1D((std::string("h_sixth_highest_CSV" + suffix)).c_str(),"Sixth Highest CSV", NcsvBins, -1.1, 1.1 );
		 
         	h_lowest_btag[c][b] = new TH1D((std::string("h_lowest_btag" + suffix)).c_str(),"Lowest Tagged CSV", NcsvBins, -1.1, 1.1 );
		 
        	h_CSV_Average[c][b] = new TH1D((std::string("h_CSV_Average" + suffix)).c_str(),"CSV Average",13,0.68,0.94);
		 	h_Deta_JetsAverage[c][b] = new TH1D((std::string("h_Deta_JetsAverage" + suffix)).c_str(),"Deta Jets Average",40,-2,2);
        	h_min_dr_tagged_jets[c][b] = new TH1D((std::string("h_min_dr_tagged_jets" + suffix)).c_str(),"Min #DeltaR(tag,tag)", 30, 0, 6. );
         	h_avg_dr_tagged_jets[c][b] = new TH1D((std::string("h_ave_dr_tagged_jets" + suffix)).c_str(),"Avg #DeltaR(tag,tag)", 30, 0, 6. );
		 	h_dr_between_lep_and_closest_jet[c][b] = new TH1D((std::string("h_dr_between_lep_and_closest_jet" + suffix)).c_str(),"Min #DeltaR(lep,jet)", 30, 0, 6. );
         	h_avg_btag_disc_btags[c][b] = new TH1D((std::string("h_avg_btag_disc_btags" + suffix)).c_str(),"Avg btag disc btags", 100, -0.01, 1.01 );
         	h_dev_from_avg_disc_btags[c][b] = new TH1D((std::string("h_dev_from_avg_disc_btags" + suffix)).c_str(),"Dev from avg disc btags", 25, -0.0001, 0.04 );
         
		 	h_HT[c][b] = new TH1D((std::string("h_HT" + suffix)).c_str(),"H_{T} (jets)", NhtBins, 0, htmax );
         	h_all_sum_pt_with_met[c][b] = new TH1D((std::string("h_all_sum_pt_with_met" + suffix)).c_str(),"Sum p_{T} (jets,lepton,MET)", NhtBins, 0, htmax );
         	h_MET_pt[c][b]  = new TH1D((std::string("h_MET_pt" + suffix)).c_str(),"MET p_{T}", NmetBins, 0, metmax );
   	     	h_MHT[c][b]  = new TH1D((std::string("h_MHT" + suffix)).c_str(),"MHT", NmetBins, 0, metmax );
         	h_invariant_mass_of_everything[c][b] = new TH1D((std::string("h_invariant_mass_of_everything" + suffix)).c_str(),"M_{inv} (jets,lepton,MET)", NhtBins, 0, htmax );
         	h_avg_untagged_dijet_mass[c][b] = new TH1D((std::string("h_avg_untagged_dijet_mass" + suffix)).c_str(),"Avg mass(untag,untag)", NmassBins, 0, massmax );
         	h_closest_tagged_dijet_mass[c][b] = new TH1D((std::string("h_closest_tagged_dijet_mass" + suffix)).c_str(),"Mass of closest tagged jets", NmassBins, 0, massmax );
         	h_M3[c][b] = new TH1D((std::string("h_M3" + suffix)).c_str(),"M3", Nm3Bins, 0, m3max );
         	h_Mlb[c][b] = new TH1D((std::string("h_Mlb" + suffix)).c_str(),"M(lep,tag)", NmassBins, 0, massmax );
         	h_tagged_dijet_mass_closest_to_125[c][b] = new TH1D((std::string("h_tagged_dijet_mass_closest_to_125" + suffix)).c_str(),"Tagged dijet mass closest to 125", NmassBins, 0, massmax );
		 
		 	h_aplanarity[c][b] = new TH1D((std::string("h_aplanarity" + suffix)).c_str(),"Aplanarity", 25, 0, 0.5 );
         	h_sphericity[c][b] = new TH1D((std::string("h_sphericity" + suffix)).c_str(),"Sphericity", 25, 0, 1.0 );
		 
         	h_h0[c][b] = new TH1D((std::string("h_h0" + suffix)).c_str(),"h0", 25, 0, 0.5 );
         	h_h1[c][b] = new TH1D((std::string("h_h1" + suffix)).c_str(),"h1", 20, -0.3, 0.5 );
         	h_h2[c][b] = new TH1D((std::string("h_h2" + suffix)).c_str(),"h2", 25, -0.2, 0.5 );
         	h_h3[c][b] = new TH1D((std::string("h_h3" + suffix)).c_str(),"h3", 25, -0.2, 1.1 );
         	h_h4[c][b] = new TH1D((std::string("h_h4" + suffix)).c_str(),"h4", 25, -0.2, 0.3 );
		 
         	h_maxeta_jet_jet[c][b] = new TH1D((std::string("h_maxeta_jet_jet" + suffix)).c_str(),"max #Delta#eta(jet,ave jet #eta)", 50, 0., 6. );
         	h_maxeta_jet_tag[c][b] = new TH1D((std::string("h_maxeta_jet_tag" + suffix)).c_str(),"max #Delta#eta(tag,ave jet #eta)", 50, 0., 6. );
         	h_maxeta_tag_tag[c][b] = new TH1D((std::string("h_maxeta_tag_tag" + suffix)).c_str(),"max #Delta#eta(tag,ave tag #eta)", 50, 0., 6. );
      
         	h_best_higgs_mass[c][b] = new TH1D((std::string("h_best_higgs_mass" + suffix)).c_str(),"Best Higgs Mass", NmassBins, 0, massmax );
         	h_minChi2[c][b] = new TH1D((std::string("h_minChi2" + suffix)).c_str(),"minChi2", 500, 0, 500 );
         	h_dRbb[c][b] = new TH1D((std::string("h_dRbb" + suffix)).c_str(),"Best #DeltaR(b,b)", 30, 0, 6. );
		 	h_abs_dEta_leptop_bb[c][b] = new TH1D((std::string("h_abs_dEta_leptop_bb" + suffix)).c_str(),"abs(#Delta#eta(leptonic top, bb))", 20, 0., 6. );
         	h_abs_dEta_hadtop_bb[c][b] = new TH1D((std::string("h_abs_dEta_hadtop_bb" + suffix)).c_str(),"abs(#Delta#eta(hadronic top, bb))", 20, 0., 6. );
         	h_dEta_fn[c][b] = new TH1D((std::string("h_dEta_fn" + suffix)).c_str(),"product(#Delta#eta(leptonic top, bb),#Delta#eta(hadronic top, bb))", 20, 0., 6.);
         	h_median_bb_mass[c][b] = new TH1D((std::string("h_median_bb_mass" + suffix)).c_str(),"median inv. mass, all btag pairs",20, 0, massmax );
         	h_pt_all_jets_over_E_all_jets[c][b] = new TH1D((std::string("h_pt_all_jets_over_E_all_jets" + suffix)).c_str(),"(#Sigma jet p_{T})/(#Sigma jet E)",50,0.,1.);
      
		 	h_bdtOutput[c][b] = new TH1D((std::string("h_bdtOutput" + suffix)).c_str(),"h_bdtOutput",20,-1,1);
		
	
		}
	}
	
	for( int b=0; b<NumSysCat; b++ ){
	    for( int c=0; c<NumCat; c++ ){
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(1+1,"bb");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(2+1,"WW");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(3+1,"#tau#tau");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(4+1,"gg");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(5+1,"cc");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(6+1,"ZZ");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(7+1,"Z#gamma");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(8+1,"#gamma#gamma");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(9+1,"ss");
	      h_higgsDecayType[c][b]->GetXaxis()->SetBinLabel(10+1,"#mu#mu");
		  
	  	  h_category_yield[b]->GetXaxis()->SetBinLabel(1,cat_labels[c].c_str());
	  	  h_category_yield_NoWgt[b]->GetXaxis()->SetBinLabel(1,cat_labels[c].c_str());
	   }
          h_CutFlow[b]->GetXaxis()->SetBinLabel(1,"All Events ");
		  h_CutFlow[b]->GetXaxis()->SetBinLabel(2,"Good 1st PV");
		  h_CutFlow[b]->GetXaxis()->SetBinLabel(3,"HLT Trigger");
		  h_CutFlow[b]->GetXaxis()->SetBinLabel(4,"1 Loose lep");
		  h_CutFlow[b]->GetXaxis()->SetBinLabel(5,"1 Tight lep");
		  h_CutFlow[b]->GetXaxis()->SetBinLabel(6,"No add. lep");
		  h_CutFlow[b]->GetXaxis()->SetBinLabel(7,"4 Jet 2 Tag");
          h_CutFlow_NoWgt[b]->GetXaxis()->SetBinLabel(1,"All");
		  h_CutFlow_NoWgt[b]->GetXaxis()->SetBinLabel(2,"Good 1st PV");
		  h_CutFlow_NoWgt[b]->GetXaxis()->SetBinLabel(3,"HLT Trigger");
		  h_CutFlow_NoWgt[b]->GetXaxis()->SetBinLabel(4,"1 Loose lep");
		  h_CutFlow_NoWgt[b]->GetXaxis()->SetBinLabel(5,"1 Tight lep");
		  h_CutFlow_NoWgt[b]->GetXaxis()->SetBinLabel(6,"No add. lep");
		  h_CutFlow_NoWgt[b]->GetXaxis()->SetBinLabel(7,"4 Jet 2 Tag");
          h_CutFlow_NoCSV[b]->GetXaxis()->SetBinLabel(1,"All");
		  h_CutFlow_NoCSV[b]->GetXaxis()->SetBinLabel(2,"Good 1st PV");
		  h_CutFlow_NoCSV[b]->GetXaxis()->SetBinLabel(3,"HLT Trigger");
		  h_CutFlow_NoCSV[b]->GetXaxis()->SetBinLabel(4,"1 Loose lep");
		  h_CutFlow_NoCSV[b]->GetXaxis()->SetBinLabel(5,"1 Tight lep");
		  h_CutFlow_NoCSV[b]->GetXaxis()->SetBinLabel(6,"No add. lep");
		  h_CutFlow_NoCSV[b]->GetXaxis()->SetBinLabel(7,"4 Jet 2 Tag");
		  
	  }
	  
	
	  

    //////////////////////////////////////////////////////////////////////////
    ///  Event Info
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
  
	//Event Weighting
  
    double MyWgt = 1;
    double ttJet_xSec = 832;
	double ttHnobb_xSec = 0.5071-0.2934;
	double ttHbb_xSec = 0.2934; //Hbb
    double ttZ_xSec = 2.232;
    double ttW_xSec = 1.152;
    double nGen = nentries;
    double nGenNeg=89671.4;

   // double MyttHwgtNeg = (ttH_xSec/nGenNeg)*10000;
	if(SampleType==1)MyWgt = (ttHnobb_xSec/nGen)*12890;
	if(SampleType==2)MyWgt = (ttJet_xSec/nGen)*12890;
    if(SampleType==3)MyWgt = (ttHbb_xSec/nGen)*12890;
   

	cout<<" BEGIN wgt ="<<MyWgt<<endl;

  
    //////////////////////////////////////////////////////////////////////////
    ///  Begin Looping over Events
    //////////////////////////////////////////////////////////////////////////
	
    
    std::cout << "========  Starting Event Loop  ========" << std::endl;
    int ievt=0;
  
    for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    //Long64_t

	    if( ievt<firstEvent ) continue;
	    if( ievt==lastEvent ) break;

      	if( ievt==0 )        std::cout << "     Event " << ievt+1 << std::endl;
      	if( ievt%10000==0 && ievt!=1) std::cout << "           " << ievt << "\t" 
  					     << int(double(ievt-1)/double(nentries)*100) << "% done" << std::endl;
	  
	 	if( ievt==(maxNentries+1) && ievt!=0 ) break;
	  
	  
	  	chain->GetEntry(ievt);
		
		long evt = eve->evt_;
		
	
		if(SampleType == 2){
		int additionalJetEventId			 = eve->additionalJetEventId_;
		if(SplitType ==1 && !(additionalJetEventId == 51) )continue; //ttb
		if(SplitType ==3 && !(additionalJetEventId == 52) )continue; //ttB
		if(SplitType ==2 && !(additionalJetEventId == 53 || additionalJetEventId == 54 || additionalJetEventId == 55))continue; //ttbb
		if(SplitType ==4 && !(additionalJetEventId == 41 || additionalJetEventId == 42 || additionalJetEventId == 43 || additionalJetEventId == 44 || additionalJetEventId == 45))continue; //ttc
		//if(SplitType ==5 && !(additionalJetEventId == 43 || additionalJetEventId == 44 || additionalJetEventId == 45))continue; //ttcc
		if(SplitType ==5 && !(additionalJetEventId == 0))continue; //ttlf
		h_additionalJetEventId->Fill(additionalJetEventId);
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
		
/////Event Selection
		
		//Set Thresholds;
		
			ThresholdList Yggdrasil_Threshold_Cuts;
			yGGhelper.SetCuts(Yggdrasil_Threshold_Cuts,25.0,15.0,2.1,2.4,0.15,0.25,30,15,2.1,2.4,0.15,0.15,30,2.4,0.8);
		
		
		
		// LEPTON SELECTION
		
			LeptonCollection uncleanedleptons;
			LeptonCollection selectedTightleptons;
			LeptonCollection selectedLooseleptons;
		
			uncleanedleptons.SetLeptonCollection(eve->lepton_pt_,eve->lepton_eta_,eve->lepton_phi_,eve->lepton_e_,eve->lepton_charge_,eve->lepton_isMuon_,eve->lepton_relIso_,eve->lepton_isTight_,eve->lepton_isLoose_);
	       	if(!usedSlimmedTrees)yGGhelper.SelectLeptons(uncleanedleptons,selectedTightleptons,selectedLooseleptons, Yggdrasil_Threshold_Cuts);
		    if(usedSlimmedTrees)selectedTightleptons=uncleanedleptons;
			
			std::vector<TLorentzVector> vTLlepton;
			vint lepnums;
		
			yGGhelper.makeVTLVleptons(vTLlepton, lepnums, selectedTightleptons, selectedLooseleptons);
			
			if(usedSlimmedTrees)lepnums=eve->lepnums_;
		
		//Begin looping on systematics
		
		
	    	int iSys=-1;
	    	for(int treeSys=0; treeSys<NumSysCat; treeSys++){
				iSys++;
				//if(treeSys!=0)continue;
		
				bool passEvent = false;	
		 
		  
		// JET SELECTION
		 
		 	   	JetCollection uncleanedJets;
		  	 	JetCollection selectedJets;

		  	  	uncleanedJets.SetJetCollection(eve->jet_pt_[treeSys],eve->jet_eta_[treeSys],eve->jet_phi_[treeSys],eve->jet_m_[treeSys],eve->jet_combinedInclusiveSecondaryVertexV2BJetTags_[treeSys]);
		  	  	if(!usedSlimmedTrees)yGGhelper.SelectJets(uncleanedJets,selectedJets,selectedLooseleptons, Yggdrasil_Threshold_Cuts);
		  		if(usedSlimmedTrees)selectedJets=uncleanedJets;
		  
		  	  	int numJets, numTags;
		 	   	std::vector<TLorentzVector> vTLjet;
		  	  	vdouble jet_CSV;
		  
		  	  	yGGhelper.makeVTLVjets(vTLjet,selectedJets,numJets,numTags, Yggdrasil_Threshold_Cuts);
		  	  	jet_CSV = selectedJets.jet_bDiscriminant;
		 
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
		
		  // Get CSV Weights		
			
				double csvWgtHF, csvWgtLF, csvWgtCF;
				double newCSVwgt = ( SampleType>5 ) ? 1 : yGGhelper.get_csv_wgt(vTLjet,jet_CSV,eve->jet_flavour_[treeSys],iSys, csvWgtHF, csvWgtLF, csvWgtCF);
				double wgt=MyWgt;
				double GenWgt = eve->wgt_generator_;
				//if(SampleType ==1)wgt=MyWgt*GenWgt; for negative weight samples
				double nocsvwgt=wgt;
				wgt=wgt*newCSVwgt;
		
				bool GoodFirstPV = eve->GoodFirstPV_;
				bool passSingleElectronTrigger = eve->passHLT_Ele27_eta2p1_WPTight_Gsf_v_;
				bool passSingleMuonTrigger = eve->passHLT_IsoMu22_v_;
				bool passSingleMuonTriggerTk = eve->passHLT_IsoTkMu22_v_;
				bool passSingleMuTriggers = (passSingleMuonTrigger || passSingleMuonTriggerTk);
				
				
				
		// Perform Event Selection
		
				int cutflowvar = yGGhelper.passCuts(GoodFirstPV,passSingleElectronTrigger,passSingleMuTriggers,lepnums,numJets,numTags);
				
				for(int i=0;i<cutflowvar;i++){
					h_CutFlow[iSys]->Fill(i,wgt);
					h_CutFlow_NoWgt[iSys]->Fill(i);
					h_CutFlow_NoCSV[iSys]->Fill(i,nocsvwgt);
				}
		
				if(iSys==0 && cutflowvar==7)passEvent=true;
				if(!passEvent)continue;
		
		        
		
	
	  		
			
			double newLepWgt = (SampleType>5) ? 1 : yGGhelper.getLeptonWgt(lepnums, vTLlepton);
			wgt=wgt*newLepWgt;
			
	  		TLorentzVector TLMET;
	  		TLMET.SetPtEtaPhiE(eve->MET_[treeSys],0,eve->MET_phi_[treeSys],eve->MET_[treeSys]);
			
			h_category_yield[iSys]->Fill(this_category,wgt);
			h_category_yield_NoWgt[iSys]->Fill(this_category);
			if(this_category<2)continue; //I want to see weighted category yields for the 4j2t and 5j2t category, but thats it.
	
		//Execute Classifier
			
			double BDToutput = kitBDTclass.GetBDTOutput(vTLlepton,vTLjet,jet_CSV,vTLjet,jet_CSV,TLMET);
			std::map<std::string,float> varsList = kitBDTclass.GetVariablesOfLastEvaluation();
			
			
	    //Fill Histograms
	       //Most Variables are calculated only within the BDTclass and can only be accessed using the variableMap
			
			h_numJets[iSys]->Fill(numJets,wgt);
			h_numTags[iSys]->Fill(numTags,wgt);
		
			

			//if(SampleType==1 || SampleType==3)h_higgsDecayType[this_category][iSys]->Fill(eve->higgsDecayType_,wgt);
			//if(SampleType==1 || SampleType==3)h_category_yield_higgsDecayType[iSys]->Fill(this_category,eve->higgsDecayType_,wgt);
		
			if(iSys==0)h_tight_lepton_pt[this_category]->Fill(selectedTightleptons.lep_pt.at(0),wgt);
			if(iSys==0)h_tight_lepton_phi[this_category]->Fill(selectedTightleptons.lep_phi.at(0),wgt);
			if(iSys==0)h_tight_lepton_eta[this_category]->Fill(selectedTightleptons.lep_eta.at(0),wgt);
			
		
		// Jet Histograms
		
  			for(int i=0;i<numJets;i++){
  				
				h_jet_pt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i),wgt);
				h_jet_pt_NoWgt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i));
				h_jet_eta[this_category][iSys]->Fill(selectedJets.jet_eta.at(i),wgt);
				h_jet_eta_NoWgt[this_category][iSys]->Fill(selectedJets.jet_eta.at(i));
				h_jet_csv[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i),wgt);
				h_jet_csv_NoWgt[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i));
			
				if(jet_CSV[i]>0.80){
				
					h_jet_tag_pt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i),wgt);
					h_jet_tag_pt_NoWgt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i));
					h_jet_tag_eta[this_category][iSys]->Fill(selectedJets.jet_eta.at(i),wgt);
					h_jet_tag_eta_NoWgt[this_category][iSys]->Fill(selectedJets.jet_eta.at(i));
					h_jet_tag_csv[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i),wgt);
					h_jet_tag_csv_NoWgt[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i));
				
				}
			
				if(i==0){
					h_first_jet_pt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i),wgt);
					h_first_jet_eta[this_category][iSys]->Fill(selectedJets.jet_eta.at(i),wgt);
					h_first_jet_csv[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i),wgt);
				
					h_first_jet_pt_NoWgt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i));
					h_first_jet_csv_NoWgt[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i));
				}
				if(i==1){
					h_second_jet_pt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i),wgt);
					h_second_jet_eta[this_category][iSys]->Fill(selectedJets.jet_eta.at(i),wgt);
					h_second_jet_csv[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i),wgt);
				}
				if(i==2){
					h_third_jet_pt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i),wgt);
					h_third_jet_eta[this_category][iSys]->Fill(selectedJets.jet_eta.at(i),wgt);
					h_third_jet_csv[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i),wgt);
				}
				if(i==3){
					h_fourth_jet_pt[this_category][iSys]->Fill(selectedJets.jet_pt.at(i),wgt);
					h_fourth_jet_eta[this_category][iSys]->Fill(selectedJets.jet_eta.at(i),wgt);
					h_fourth_jet_csv[this_category][iSys]->Fill(selectedJets.jet_bDiscriminant.at(i),wgt);
				}
  		
			}
		
		
			//h_first_highest_btag[this_category][iSys]->Fill(bdt_osu_30_2p4.GetBDTvariable(9),wgt);
			h_second_highest_btag[this_category][iSys]->Fill(varsList["second_highest_btag"],wgt);
			h_third_highest_btag[this_category][iSys]->Fill(varsList["third_highest_btag"],wgt);
			h_fourth_highest_btag[this_category][iSys]->Fill(varsList["fourth_highest_btag"],wgt);
			h_fifth_highest_CSV[this_category][iSys]->Fill(varsList["fifth_highest_CSV"],wgt);
			h_sixth_highest_CSV[this_category][iSys]->Fill(varsList["sixth_highest_CSV"],wgt);
		
			h_lowest_btag[this_category][iSys]->Fill(varsList["lowest_btag"],wgt);
		
		//Multi Object Plots//
		
			h_CSV_Average[this_category][iSys]->Fill(varsList["CSV_Average"],wgt);
			h_Deta_JetsAverage[this_category][iSys]->Fill(varsList["Deta_JetsAverage"],wgt);
			h_min_dr_tagged_jets[this_category][iSys]->Fill(varsList["min_dr_tagged_jets"],wgt);
			h_avg_dr_tagged_jets[this_category][iSys]->Fill(varsList["avg_dr_tagged_jets"],wgt);
			h_dr_between_lep_and_closest_jet[this_category][iSys]->Fill(varsList["dr_between_lep_and_closest_jet"],wgt);
			h_avg_btag_disc_btags[this_category][iSys]->Fill(varsList["avg_btag_disc_btags"],wgt);
			h_dev_from_avg_disc_btags[this_category][iSys]->Fill(varsList["dev_from_avg_disc_btags"],wgt);
		
			h_HT[this_category][iSys]->Fill(varsList["HT"],wgt);
			h_all_sum_pt_with_met[this_category][iSys]->Fill(varsList["all_sum_pt_with_met"],wgt);
			h_MET_pt[this_category][iSys]->Fill(varsList["MET_pt"],wgt);
			h_MHT[this_category][iSys]->Fill(varsList["MHT"],wgt);
			h_invariant_mass_of_everything[this_category][iSys]->Fill(varsList["invariant_mass_of_everything"],wgt);
			h_avg_untagged_dijet_mass[this_category][iSys]->Fill(varsList["avg_untagged_dijet_mass"],wgt);
			h_closest_tagged_dijet_mass[this_category][iSys]->Fill(varsList["closest_tagged_dijet_mass"],wgt);
			h_M3[this_category][iSys]->Fill(varsList["M3"],wgt);
			h_Mlb[this_category][iSys]->Fill(varsList["Mlb"],wgt);
			h_tagged_dijet_mass_closest_to_125[this_category][iSys]->Fill(varsList["tagged_dijet_mass_closest_to_125"],wgt);
			
		//BDTvars variables//
	
			h_aplanarity[this_category][iSys]->Fill(varsList["aplanarity"],wgt);
			h_sphericity[this_category][iSys]->Fill(varsList["sphericity"],wgt);
		
			h_h0[this_category][iSys]->Fill(varsList["h0"],wgt);
			h_h1[this_category][iSys]->Fill(varsList["h1"],wgt);
			h_h2[this_category][iSys]->Fill(varsList["h2"],wgt);
			h_h3[this_category][iSys]->Fill(varsList["h3"],wgt);
			h_h4[this_category][iSys]->Fill(varsList["h4"],wgt);
		
			h_maxeta_jet_jet[this_category][iSys]->Fill(varsList["maxeta_jet_jet"],wgt);
			h_maxeta_jet_tag[this_category][iSys]->Fill(varsList["maxeta_jet_tag"],wgt);
			h_maxeta_tag_tag[this_category][iSys]->Fill(varsList["maxeta_tag_tag"],wgt);
			
			h_best_higgs_mass[this_category][iSys]->Fill(varsList["best_higgs_mass"],wgt);
			h_minChi2[this_category][iSys]->Fill(varsList["minChi2"],wgt);
			h_dRbb[this_category][iSys]->Fill(varsList["dRbb"],wgt);
			h_abs_dEta_leptop_bb[this_category][iSys]->Fill(varsList["abs_dEta_leptop_bb"],wgt);
			h_abs_dEta_hadtop_bb[this_category][iSys]->Fill(varsList["abs_dEta_hadtop_bb"],wgt);
			h_dEta_fn[this_category][iSys]->Fill(varsList["dEta_fn"],wgt);
			h_median_bb_mass[this_category][iSys]->Fill(varsList["median_bb_mass"],wgt);
			h_pt_all_jets_over_E_all_jets[this_category][iSys]->Fill(varsList["pt_all_jets_over_E_all_jets"],wgt);
		
		//BDT Histograms
		
			h_bdtOutput[this_category][iSys]->Fill(BDToutput);
		
  		} // end loop over systematics
	}// end event loop

//Print Concluding Information
	
	//Cut Flow
	
	int cutfflow[7][NumSysCat];
	string cutflowlabel[7];
	double percentchange=-99;
   
	for(int i=0;i<7;i++){
		cutflowlabel[i]=h_CutFlow[0]->GetXaxis()->GetBinLabel(i+1);
		cutfflow[i][0]=h_CutFlow_NoWgt[0]->GetBinContent(i+1);
		cout<<cutflowlabel[i]<<" - "<<cutfflow[i][0]<<endl;
	}
	cout<<endl;
	 
	//Weighted Event yields by category 
	 
	float yyields[7][NumSysCat];
   
	for(int i=0;i<9;i++){
		yyields[i][0]=h_category_yield[0]->GetBinContent(i+1);
	    cout<<cat_labels[i]<<" - "<<yyields[i][0]<<endl;
	}
	cout<<endl;
 
    //Timer
	t=clock()-t;
	std::cout << " Done! " <<((float)t)/CLOCKS_PER_SEC<< std::endl;
	

 histofile.Write();
 histofile.Close();

}









