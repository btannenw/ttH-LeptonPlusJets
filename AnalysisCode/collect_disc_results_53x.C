#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TAxis.h"
#include "TList.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TKey.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>


void collect_disc_results_53x( int lepSelection=-1, int ANNtype_=0 ) {


  ///////////////////////////////////////////////////////

  int Nsample = 6;//38;
  int Nsignal = 1; //9+72;
  int Nsample_all = Nsample + Nsignal;
  TFile *file[7];

  TString lepType = "";
  if( lepSelection==0 )      lepType = "mu_sel";
  else if( lepSelection==1 ) lepType = "ele_sel";
  else                       lepType = "lep_sel";

  TString baseDir = "/eos/uscms/store/user/satoshi/testout/" ;
  //  TString baseDir = "/afs/cern.ch/work/l/lwming/";
  TString basePrefix = "yggdrasil_treeReader";
  //  TString baseSuffix = "_" + lepType + "_histo.root";
  TString baseSuffix =  ".root";

  std::vector<TString> datasets(Nsample_all);
  datasets[0]  = "SingleMu_2012ABCD_BEAN_53xOn53x";

  //datasets[1] = "ttZ";
  //datasets[2] = "ttW";

  datasets[1] = "ttlf";
  datasets[2] = "ttcc";
  datasets[3] = "ttbb";
  datasets[4] = "ttb";

  datasets[5] = "ttH";

  datasets[6] = "tt2b";



  if( lepSelection==0 )      datasets[0]  = "SingleMu_2012ABCD_BEAN_53xOn53x";
  else if( lepSelection==1 ) datasets[0]  = "SingleElectron_2012ABCD_BEAN_53xOn53x";
  else                       datasets[0]  = "data";

  // TString beanVer = "_beanVer2";

  for( int iSample=0; iSample<Nsample_all; iSample++ ){
    // TString typePrefix = ( datasets[iSample].Contains("Single") ) ? "_data_" : "_mc_";
    
    // baseDir = "/home/puigh/work/physics/8TeV/CMSSW_5_2_6/src/Output/HistoFiles/treeReader/v25_newBEANs_newCharmCSVSFUnc/";

    file[iSample] = new TFile(baseDir + datasets[iSample] + baseSuffix);
  }



  std::vector<std::string> sample;
  sample.push_back("data_obs");
  // sample.push_back("diboson");
  // sample.push_back("zjets");
  // sample.push_back("wjets");
  // sample.push_back("singlet");
  sample.push_back("ttbar");
  sample.push_back("ttbarPlusCCbar");
  sample.push_back("ttbarPlusBBbar");
  sample.push_back("ttbarPlusB");
  //sample.push_back("ttbarW");
  //sample.push_back("ttbarZ");

  sample.push_back("ttH125");
  sample.push_back("ttbarPlus2B");




  TString file_suffix = "_final10v1_8TeV_CFMlpANN.root";
  if( ANNtype_==10 )      file_suffix = "_final10v10_8TeV_BDT.root";
  else if( ANNtype_==15 ) file_suffix = "_final10v15_8TeV_BDT.root";
  else if( ANNtype_==16 ) file_suffix = "_final10v16_8TeV_BDT.root";
  else if( ANNtype_==17 ) file_suffix = "_final10v17_8TeV_BDT.root";
  else if( ANNtype_==18 ) file_suffix = "_final10v18_8TeV_BDT.root";
  else if( ANNtype_==19 ) file_suffix = "_final10v19_8TeV_BDT.root";

  file_suffix = "LJ_13TeV_BDT.root";
  TString histofilename =  file_suffix;


  TFile histofile(histofilename,"recreate");

  histofile.cd();



  std::vector<std::string> histoname1;
  if( ANNtype_==10)      histoname1.push_back("h_disc_final10v10_8TeV_CFMlpANN_BDT");
  else if( ANNtype_==15) histoname1.push_back("h_disc_final10v15_8TeV_CFMlpANN_BDT");
  else if( ANNtype_==16) histoname1.push_back("h_disc_final10v16_8TeV_CFMlpANN_BDT");
  else if( ANNtype_==17) histoname1.push_back("h_disc_final10v17_8TeV_CFMlpANN_BDT");
  else if( ANNtype_==18) histoname1.push_back("h_disc_final10v16_8TeV_CFMlpANN_BDT");
  else if( ANNtype_==19) histoname1.push_back("h_disc_final10v16_8TeV_CFMlpANN_BDT");
  else histoname1.push_back("h_BDT_Output_IncludeWgt_even_");
//else histoname1.push_back("h_SplittingOutput_1_1_");

  std::vector<std::string> old_catNames;
  // old_catNames.push_back("4j2t");
  old_catNames.push_back("4j3t");
  old_catNames.push_back("4j4t");
  // old_catNames.push_back("5j2t");
  old_catNames.push_back("5j3t");
  old_catNames.push_back("5j4t");
  old_catNames.push_back("6j2t");
  old_catNames.push_back("6j3t");
  old_catNames.push_back("6j4t");

  std::vector<std::string> old_sysLabels;
  old_sysLabels.push_back("");
  // old_sysLabels.push_back("_PUUp");
  // old_sysLabels.push_back("_PUDown");
  // old_sysLabels.push_back("_hfSFUp");
  // old_sysLabels.push_back("_hfSFDown");
  // old_sysLabels.push_back("_lfSFUp");
  // old_sysLabels.push_back("_lfSFDown");
   old_sysLabels.push_back("_JERUp");
   old_sysLabels.push_back("_JERDown");
  old_sysLabels.push_back("_JESUp");
  old_sysLabels.push_back("_JESDown");
  // old_sysLabels.push_back("_topPtcorrUp");
  // old_sysLabels.push_back("_topPtcorrDown");
  // old_sysLabels.push_back("_Q2scaleUp");
  // old_sysLabels.push_back("_Q2scaleDown");
  //
  old_sysLabels.push_back("_CSVLFUp");     //17
  old_sysLabels.push_back("_CSVLFDown");   //18
  old_sysLabels.push_back("_CSVHFUp");     //19
  old_sysLabels.push_back("_CSVHFDown");   //20
  old_sysLabels.push_back("_CSVHFStats1Up");     //21
  old_sysLabels.push_back("_CSVHFStats1Down");   //22
  old_sysLabels.push_back("_CSVLFStats1Up");     //23
  old_sysLabels.push_back("_CSVLFStats1Down");   //24
  old_sysLabels.push_back("_CSVHFStats2Up");     //25
  old_sysLabels.push_back("_CSVHFStats2Down");   //26
  old_sysLabels.push_back("_CSVLFStats2Up");     //27
  old_sysLabels.push_back("_CSVLFStats2Down");   //28
   old_sysLabels.push_back("_CSVCErr1Up");     //29
   old_sysLabels.push_back("_CSVCErr1Down");   //30
   old_sysLabels.push_back("_CSVCErr2Up");     //31
   old_sysLabels.push_back("_CSVCErr2Down");   //32

  std::vector<std::string> new_sysLabels;
  new_sysLabels.push_back("");
  // new_sysLabels.push_back("_CMS_ttH_puUp");
  // new_sysLabels.push_back("_CMS_ttH_puDown");
  // new_sysLabels.push_back("_CMS_eff_bUp");
  // new_sysLabels.push_back("_CMS_eff_bDown");
  // new_sysLabels.push_back("_CMS_fake_bUp");
  // new_sysLabels.push_back("_CMS_fake_bDown");
   new_sysLabels.push_back("_CMS_res_jUp");
   new_sysLabels.push_back("_CMS_res_jDown");
  new_sysLabels.push_back("_CMS_scale_jUp");
  new_sysLabels.push_back("_CMS_scale_jDown");
  // new_sysLabels.push_back("_CMS_ttH_topPtcorrUp");
  // new_sysLabels.push_back("_CMS_ttH_topPtcorrDown");
  // new_sysLabels.push_back("_Q2scale_ttH_ttbarUp");
  // new_sysLabels.push_back("_Q2scale_ttH_ttbarDown");
  //
  new_sysLabels.push_back("_CMS_ttH_CSVLFUp");     //17
  new_sysLabels.push_back("_CMS_ttH_CSVLFDown");   //18
  new_sysLabels.push_back("_CMS_ttH_CSVHFUp");     //19
  new_sysLabels.push_back("_CMS_ttH_CSVHFDown");   //20
  new_sysLabels.push_back("_CMS_ttH_CSVHFStats1Up");     //21
  new_sysLabels.push_back("_CMS_ttH_CSVHFStats1Down");   //22
  new_sysLabels.push_back("_CMS_ttH_CSVLFStats1Up");     //23
  new_sysLabels.push_back("_CMS_ttH_CSVLFStats1Down");   //24
  new_sysLabels.push_back("_CMS_ttH_CSVHFStats2Up");     //25
  new_sysLabels.push_back("_CMS_ttH_CSVHFStats2Down");   //26
  new_sysLabels.push_back("_CMS_ttH_CSVLFStats2Up");     //27
  new_sysLabels.push_back("_CMS_ttH_CSVLFStats2Down");   //28
   new_sysLabels.push_back("_CMS_ttH_CSVCErr1Up");     //29
   new_sysLabels.push_back("_CMS_ttH_CSVCErr1Down");   //30
   new_sysLabels.push_back("_CMS_ttH_CSVCErr2Up");     //31
   new_sysLabels.push_back("_CMS_ttH_CSVCErr2Down");   //32



  std::vector<std::string> new_catNames;
  // new_catNames.push_back("j4_t2");
  new_catNames.push_back("j4_t3");
  new_catNames.push_back("j4_t4");
  // new_catNames.push_back("j5_t2");
  new_catNames.push_back("j5_t3");
  new_catNames.push_back("j5_tge4");
  new_catNames.push_back("jge6_t2");
  new_catNames.push_back("jge6_t3");
  new_catNames.push_back("jge6_tge4");




  // zjets
  // int bin_zjets = 1;

  // wjets
  // int bin_wjets = 7;

  // diboson

  // single top

  // ttV
//  int bin_ttZ     = 1;
//  int bin_ttW     = 2;

  // tt + lf
  int bin_ttbar_lf    = 1;

  // tt + cc
  int bin_ttbar_cc    = 2;

  // tt + bb
  int bin_ttbar_bb    = 3;

  // tt + b
  int bin_ttbar_b = 4;

  // ttH
  int bin_ttH125 = 5;

  int bin_ttbar_B = 6;


  std::vector<int> use_data;
  use_data.push_back(0);

  std::vector<int> bkg_ttbar_b;
  bkg_ttbar_b.push_back(bin_ttbar_b);

  std::vector<int> bkg_ttbar_B;
  bkg_ttbar_B.push_back(bin_ttbar_B);

  std::vector<int> bkg_ttbar_bbbar;
  bkg_ttbar_bbbar.push_back(bin_ttbar_bb);

  std::vector<int> bkg_ttbar_ccbar;
  bkg_ttbar_ccbar.push_back(bin_ttbar_cc);

  std::vector<int> bkg_ttbar_other;
  bkg_ttbar_other.push_back(bin_ttbar_lf);


  //std::vector<int> bkg_ttZ;
  //bkg_ttZ.push_back(bin_ttZ);
  //std::vector<int> bkg_ttW;
  //bkg_ttW.push_back(bin_ttW);


  std::vector<int> sig_ttH125;
  sig_ttH125.push_back(bin_ttH125);




  std::vector<std::vector<int> > samples_to_use;
  samples_to_use.push_back(use_data);
  // samples_to_use.push_back(bkg_diboson);
  // samples_to_use.push_back(bkg_zjets);
  // samples_to_use.push_back(bkg_wjets);
  // samples_to_use.push_back(bkg_singlet);
  samples_to_use.push_back(bkg_ttbar_other);
  samples_to_use.push_back(bkg_ttbar_ccbar);
  samples_to_use.push_back(bkg_ttbar_bbbar);
  samples_to_use.push_back(bkg_ttbar_b);
  //samples_to_use.push_back(bkg_ttW);
  //samples_to_use.push_back(bkg_ttZ);

  samples_to_use.push_back(sig_ttH125);
  samples_to_use.push_back(bkg_ttbar_B);

  //////////

  for( int iNSamp=0; iNSamp<int(sample.size()); iNSamp++ ){

    for( int i=0; i<int(histoname1.size()); i++ ){
      std::string hname_tmp = histoname1[i];

      for( int j=0; j<int(old_catNames.size()); j++ ){

	// if( (old_catNames[j].find("6j3t")!=std::string::npos) ){
	//   if( ANNtype_==18 )      hname_tmp = "h_disc_final10v15_8TeV_CFMlpANN_BDT";
	// }
	// if( (old_catNames[j].find("5j4t")!=std::string::npos) ){
	//   if( ANNtype_==19 )      hname_tmp = "h_disc_final10v15_8TeV_CFMlpANN_BDT";
	// }

	// h_Reader_Output_4j3t_final
	for( int k=0; k<int(old_sysLabels.size()); k++ ){
	  if( (sample[iNSamp].find("data")!=std::string::npos) && k!=0 ) continue;

	  //std::string hname = hname_tmp + "_" + old_catNames[j] + "_finall_May10_optimized" + old_sysLabels[k];
	  std::string hname = hname_tmp + old_catNames[j] + old_sysLabels[k];
	  std::string hname_use = hname_tmp + "_use_" + old_catNames[j] + old_sysLabels[k];

	  bool firstSample = true;
	  TH1D* hist;
	  for( int iSample=0; iSample<Nsample_all; iSample++ ){
	    bool useSample = false;
	    for( int jSample=0; jSample<int(samples_to_use[iNSamp].size()); jSample++ ){
	      if( iSample==samples_to_use[iNSamp][jSample] ){
		useSample = true;
		break;
	      }
	    }
	    if( !useSample ) continue;

	    TH1D* hist_temp = (TH1D*)file[iSample]->Get(hname.c_str());

	    if( firstSample ){
	      hist = (TH1D*)hist_temp->Clone();
	      firstSample = false;
	    }
	    else hist->Add(hist_temp);
	  }

	  std::string temp = hname;
	  int rebin = 1;

	  //TH1D* hist_use;
	  //	  if( (hname.find("disc")!=std::string::npos) && (hname.find("CFMlpANN")!=std::string::npos) ){
	    if( (hname.find("4t")!=std::string::npos) )      rebin = 2;
	    if( (hname.find("TTHTTBB")!=std::string::npos) ) rebin = 2;

	  TH1D*  hist_use = (TH1D*)hist->Clone(hname_use.c_str());
	    hist_use->Rebin(rebin);
	  // }
	  // else if( (hname.find("avg_btag_disc_btags")!=std::string::npos)){
	  //   if( (hname.find("4j1t")!=std::string::npos) )     rebin = 1;
	  //   else if( (hname.find("j2t")!=std::string::npos) ) rebin = 1;
	  //   else if( (hname.find("j3t")!=std::string::npos) ) rebin = 2;
	  //   else if( (hname.find("j4t")!=std::string::npos) ) rebin = 2;

	  //   hist_use = (TH1D*)hist->Clone(hname_use.c_str());
	  //   hist_use->Rebin(rebin);
	  // }


	  int nbins = hist_use->GetNbinsX();
	  for( int bin=0; bin<nbins; bin++ ){
	    if( !(hist_use->GetBinContent(bin+1)>0.) ) hist_use->SetBinContent(bin+1,0.00001);
	  }

	  // std::string new_prefix = "CFMlpANN";
	  std::string new_prefix = "MVA";

	  std::string new_syst_label = new_sysLabels[k];
	  if( (new_sysLabels[k].find("Q2scale_ttH_ttbarUp")!=std::string::npos)  ){
	    if( (sample[iNSamp].find("ttbarPlusBB")!=std::string::npos)  )          new_syst_label = "_Q2scale_ttH_ttbar_bbUp";
	    else if( (sample[iNSamp].find("ttbarPlusB")!=std::string::npos)  )    new_syst_label = "_Q2scale_ttH_ttbar_bUp";
	    else if( (sample[iNSamp].find("ttbarPlusCCbar")!=std::string::npos)  ) new_syst_label = "_Q2scale_ttH_ttbar_ccUp";
	  }
	  else if( (new_sysLabels[k].find("Q2scale_ttH_ttbarDown")!=std::string::npos)  ){
	    if( (sample[iNSamp].find("ttbarPlusBB")!=std::string::npos)  )          new_syst_label = "_Q2scale_ttH_ttbar_bbDown";
	    else if( (sample[iNSamp].find("ttbarPlusB")!=std::string::npos)  )    new_syst_label = "_Q2scale_ttH_ttbar_bDown";
	    else if( (sample[iNSamp].find("ttbarPlusCCbar")!=std::string::npos)  ) new_syst_label = "_Q2scale_ttH_ttbar_ccDown";
	  }

	  std::string new_hname = sample[iNSamp] + "_" + new_prefix + "_ljets_" + new_catNames[j] + new_syst_label;

	  hist_use->Write(new_hname.c_str());
	}
      }
    }
  }


  histofile.Write();
  histofile.Close();

  // close all files
  for( int j=0; j<Nsample_all; j++ ) file[j]->Close();

}


/*

TH1D* h_data = (TH1D*)data_obs_CFMlpANN_ljets_jge6_tge4->Clone();

TH1D* h_bkg = (TH1D*)diboson_CFMlpANN_ljets_jge6_tge4->Clone();
h_bkg->Add(zjets_CFMlpANN_ljets_jge6_tge4);
h_bkg->Add(wjets_CFMlpANN_ljets_jge6_tge4);
h_bkg->Add(singlet_CFMlpANN_ljets_jge6_tge4);
h_bkg->Add(ttbar_CFMlpANN_ljets_jge6_tge4);
h_bkg->Add(ttbarPlusCCbar_CFMlpANN_ljets_jge6_tge4);
h_bkg->Add(ttbarPlusBBbar_CFMlpANN_ljets_jge6_tge4);
h_bkg->Add(ttbarW_CFMlpANN_ljets_jge6_tge4);
h_bkg->Add(ttbarZ_CFMlpANN_ljets_jge6_tge4);

TH1D* h_bkg_CMS_ttH_CSVLFUp = (TH1D*)diboson_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp->Clone();
h_bkg_CMS_ttH_CSVLFUp->Add(zjets_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp);
h_bkg_CMS_ttH_CSVLFUp->Add(wjets_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp);
h_bkg_CMS_ttH_CSVLFUp->Add(singlet_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp);
h_bkg_CMS_ttH_CSVLFUp->Add(ttbar_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp);
h_bkg_CMS_ttH_CSVLFUp->Add(ttbarPlusCCbar_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp);
h_bkg_CMS_ttH_CSVLFUp->Add(ttbarPlusBBbar_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp);
h_bkg_CMS_ttH_CSVLFUp->Add(ttbarW_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp);
h_bkg_CMS_ttH_CSVLFUp->Add(ttbarZ_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFUp);

TH1D* h_bkg_CMS_ttH_CSVLFDown = (TH1D*)diboson_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown->Clone();
h_bkg_CMS_ttH_CSVLFDown->Add(zjets_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown);
h_bkg_CMS_ttH_CSVLFDown->Add(wjets_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown);
h_bkg_CMS_ttH_CSVLFDown->Add(singlet_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown);
h_bkg_CMS_ttH_CSVLFDown->Add(ttbar_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown);
h_bkg_CMS_ttH_CSVLFDown->Add(ttbarPlusCCbar_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown);
h_bkg_CMS_ttH_CSVLFDown->Add(ttbarPlusBBbar_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown);
h_bkg_CMS_ttH_CSVLFDown->Add(ttbarW_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown);
h_bkg_CMS_ttH_CSVLFDown->Add(ttbarZ_CFMlpANN_ljets_jge6_tge4_CMS_ttH_CSVLFDown);

h_data->SetMarkerStyle(20);
h_bkg->SetLineWidth(2);
h_bkg_CMS_ttH_CSVLFUp->SetLineWidth(2);
h_bkg_CMS_ttH_CSVLFDown->SetLineWidth(2);

h_bkg->SetLineColor(kBlack);
h_bkg_CMS_ttH_CSVLFUp->SetLineColor(kRed);
h_bkg_CMS_ttH_CSVLFDown->SetLineColor(kBlue);


h_data->SetStats(0);

h_data->Draw("pe1");
h_bkg->Draw("histsame");
h_bkg_CMS_ttH_CSVLFUp->Draw("histsame");
h_bkg_CMS_ttH_CSVLFDown->Draw("histsame");


TFile* file_v16 = new TFile("HistoFiles/collect_disc_results_53x_histo_3rd40_2013_06_24_v1_LJ_lep_sel_final10v16_8TeV_BDT.root");
TFile* file_v19 = new TFile("HistoFiles/collect_disc_results_53x_histo_3rd40_2013_06_24_v1_LJ_lep_sel_final10v19_8TeV_BDT.root");

TH1D* ttH_v16 = (TH1D*)file_v16->Get("ttH125_MVA_ljets_jge6_t3");
TH1D* ttH_v19 = (TH1D*)file_v19->Get("ttH125_MVA_ljets_jge6_t3");

ttH_v16->Draw();
ttH_v19->Draw("histsame");






   1 & avg\_dr\_tagged\_jets               & 1.068e-01 \\
   2 & maxeta\_tag\_tag                   & 8.248e-02 \\
   3 & pt\_all\_jets\_over\_E\_all\_jets      & 7.428e-02 \\
   4 & h1                               & 6.709e-02 \\
   5 & tagged\_dijet\_mass\_closest\_to\_125 & 6.667e-02 \\
   6 & all\_sum\_pt\_with\_met              & 6.415e-02 \\
   7 & fourth\_highest\_btag              & 6.177e-02 \\
   8 & aplanarity                       & 6.070e-02 \\
   9 & h3                               & 5.450e-02 \\
  10 & MET                              & 3.341e-02 \\


   1 & h3                               & 2.854e-02 \\
   2 & maxeta\_tag\_tag                   & 2.824e-02 \\
   3 & dEta\_fn                          & 2.784e-02 \\
   4 & pt\_all\_jets\_over\_E\_all\_jets      & 2.664e-02 \\
   5 & abs\_dEta\_hadtop\_bb               & 2.374e-02 \\
   6 & tagged\_dijet\_mass\_closest\_to\_125 & 2.349e-02 \\
   7 & maxeta\_jet\_jet                   & 2.259e-02 \\
   8 & M3                               & 2.216e-02 \\
   9 & maxeta\_jet\_tag                   & 2.161e-02 \\
  10 & abs\_dEta\_leptop\_bb               & 2.098e-02 \\
  11 & h1                               & 1.991e-02 \\
  12 & sphericity                       & 1.923e-02 \\
  13 & aplanarity                       & 1.845e-02 \\
  14 & min\_dr\_tagged\_jets               & 1.579e-02 \\
  15 & third\_jet\_pt                     & 1.127e-02 \\



 */
