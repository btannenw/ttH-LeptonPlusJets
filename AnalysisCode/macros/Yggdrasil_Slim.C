#include "TFile.h"
#include "TChain.h"

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


#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif


#if !defined(__CINT__) && !defined(__MAKECINT__)

#include "ttH-LeptonPlusJets/AnalysisCode/interface/YggHelper.h"
#include "ttH-LeptonPlusJets/AnalysisCode/interface/YggdrasilEventVars.h"

#endif




typedef std::vector< TLorentzVector >          vecTLorentzVector;
typedef std::vector<int>                       vint;
typedef std::vector<double>                    vdouble;
typedef std::vector<std::vector<double> >      vvdouble;

void make_the_tree(std::string newfilename = "newfilename.root", std::string treefilename = "readfile.root", int maxNentries=-1, int Njobs=1, int jobN=1);

void Yggdrasil_Slim(int SampleType =2, int maxNentries=-1, int Njobs=1, int jobN=1 )
{
    std::string str_jobN;
    std::stringstream stream;
    stream << jobN;
    str_jobN = stream.str();
    
    std::vector<string> treefilename;
    
    std::string SampleName;
    
    
    if(SampleType == 1){
    	SampleName = "ttHToNonbb";
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/ttHToNonbb_M125_13TeV_powheg_pythia8/tthnonbb_Satoshi_Yggdra_Aug01/160801_165219/0000/*.root");
    }
    
    
    if(SampleType == 2){
    	 SampleName = "ttJets";
  	 	treefilename.push_back("/eos/uscms/store/user/lpctthrun2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/ttbar_Satoshi_Yggdra_Aug01/160801_165115/0000/*.root");
	 	treefilename.push_back("/eos/uscms/store/user/lpctthrun2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/ttbar_Satoshi_Yggdra_Aug01/160801_165115/0001/*.root");
	 	treefilename.push_back("/eos/uscms/store/user/lpctthrun2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/ttbar_Satoshi_Yggdra_Aug01/160801_165115/0002/*.root");
    }
    
    
    if(SampleType == 3){
    	SampleName = "ttHTobb";
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/ttHTobb_M125_13TeV_powheg_pythia8/tthbb_Satoshi_Yggdra_Aug01/160801_165143/0000/*.root");
    }
    
    if(SampleType == 6){
        SampleName = "SingleElectronB";
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleElectron/DataMuB_Satoshi_Yggdra_Aug01/160801_164829/0000/*.root");
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleElectron/DataMuB_Satoshi_Yggdra_Aug01/160801_164829/0001/*.root");
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleElectron/DataMuB_Satoshi_Yggdra_Aug01/160801_164829/0002/*.root");
    }
    
    
    if(SampleType == 7){
    	SampleName = "SingleElectronC";
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleElectron/DataMuC_Satoshi_Yggdra_Aug01/160801_164858/0000/*.root");
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleElectron/DataMuC_Satoshi_Yggdra_Aug01/160801_164858/0001/*.root");
    }
    
    if(SampleType == 8){
        SampleName = "SingleElectronD";
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleElectron/DataMuD_Satoshi_Yggdra_Aug01/160801_164939/0000/*.root");
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleElectron/DataMuD_Satoshi_Yggdra_Aug01/160801_164939/0001/*.root");
    }
    
    if(SampleType == 9){
    	SampleName = "SingleMuonB";
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleMuon/DataElB_Satoshi_Yggdra_Aug01/160801_164709/0000/*.root");
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleMuon/DataElB_Satoshi_Yggdra_Aug01/160801_164709/0001/*.root");
    }
    
    if(SampleType == 10){
        SampleName = "SingleMuonC";
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleMuon/DataElC_Satoshi_Yggdra_Aug01/160801_164732/0000/*.root");
    }
    
    if(SampleType == 11){
    	SampleName = "SingleMuonD";
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleMuon/DataElD_Satoshi_Yggdra_Aug01/160801_164758/0000/*.root");
		treefilename.push_back("/eos/uscms/store/user/lpctthrun2/SingleMuon/DataElD_Satoshi_Yggdra_Aug01/160801_164758/0001/*.root");
    }
    
  
 
    std::string OutputDir = "";
    //std::string OutputDir ="/eos/uscms/store/user/sflowers/SlimTrees/";
    for(int iSample=0; iSample<int(treefilename.size());iSample++){
    	make_the_tree(OutputDir + "condortestSlimTrees_Aug5th_" + SampleName + "_" + std::to_string(iSample) + "_" + str_jobN + ".root",treefilename[iSample],maxNentries,Njobs,jobN); cout<<SampleName<<" done "<<iSample+1<<"/"<<int(treefilename.size())<<endl;
    }
}


void make_the_tree(std::string newfilename, std::string treefilename, int maxNentries, int Njobs, int jobN)
{
    std::cout << "   ===> load the root file! " << std::endl;
	std::cout << "  treefilename  = " << treefilename.c_str() << std::endl;
	std::cout << "  newfilename   = " << newfilename.c_str() << std::endl;
	TChain *chain = new TChain("ttHTreeMaker/worldTree");
	
        chain->Add(treefilename.c_str());
	
	yggdrasilEventVars *eve=0;
	yggdrasilEventVars *neve=0;
       chain->SetBranchAddress("eve.", &eve );
	
	////// Create new shaved down tree (Shrub) and add variables
	
	TFile *newtreefile = new TFile(newfilename.c_str(),"RECREATE","testfile");
	TDirectory *ttHTreeMaker = newtreefile->mkdir("ttHTreeMaker");
	ttHTreeMaker->cd();
	TTree *YggdrasilSlim = new TTree("worldTree", "YggdrasilSlim");
	YggdrasilSlim->Branch("eve.", "yggdrasilEventVars" , &neve);
	
	YggHelper yGGhelper;
	
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
    

	
    std::cout << "========  Starting Event Loop  ========" << std::endl;
    for (Long64_t ievt=0; ievt<chain->GetEntries();ievt++) {    //Long64_t

	    if( ievt<firstEvent ) continue;
	    if( ievt==lastEvent ) break;

      	if( ievt==0 )        std::cout << "     Event " << ievt+1 << std::endl;
      	if( ievt%10000==0 && ievt!=1) std::cout << "           " << ievt << "\t" 
  					     << int(double(ievt-1)/double(nentries)*100) << "% done" << std::endl;
	  
	 	if( ievt==(maxNentries+1) && ievt!=0 ) break;
	  
	  
	  	chain->GetEntry(ievt);
		
		
		neve->evt_ = eve->evt_;
		neve->additionalJetEventId_ = eve->additionalJetEventId_;
		neve->higgsDecayType_ = eve->higgsDecayType_;
		
		ThresholdList Yggdrasil_Threshold_Cuts;
		yGGhelper.SetCuts(Yggdrasil_Threshold_Cuts,25.0,15.0,2.1,2.4,0.15,0.25,30,15,2.1,2.4,0.15,0.15,30,2.4,0.8);
		
		//Slim Leptons
		
		LeptonCollection uncleanedleptons;
		LeptonCollection selectedTightleptons;
		LeptonCollection selectedLooseleptons;
		
		uncleanedleptons.SetLeptonCollection(eve->lepton_pt_,eve->lepton_eta_,eve->lepton_phi_,eve->lepton_e_,eve->lepton_charge_,eve->lepton_isMuon_,eve->lepton_relIso_,eve->lepton_isTight_,eve->lepton_isLoose_);
	        yGGhelper.SelectLeptons(uncleanedleptons,selectedTightleptons,selectedLooseleptons, Yggdrasil_Threshold_Cuts);
		
		std::vector<TLorentzVector> vTLlepton;
		vint lepnums;
		
		yGGhelper.makeVTLVleptons(vTLlepton, lepnums, selectedTightleptons, selectedLooseleptons);
		
		
		
		neve->lepton_pt_ = selectedTightleptons.lep_pt;
		neve->lepton_eta_ = selectedTightleptons.lep_eta;
		neve->lepton_phi_ = selectedTightleptons.lep_phi;
		neve->lepton_e_ = selectedTightleptons.lep_e;
		neve->lepton_charge_ = selectedTightleptons.lep_charge;
		neve->lepton_isMuon_ = selectedTightleptons.lep_isMuon;
		neve->lepton_relIso_ = selectedTightleptons.lep_relIso;
		neve->lepton_isTight_ = selectedTightleptons.lep_POGTight;
		neve->lepton_isLoose_ = selectedTightleptons.lep_POGLoose;
		neve->lepnums_ = lepnums;
		
		//Loop of Systematics
		
		int numSys = 0;
		int treeSys=0;
		
		//Slim Jets
		
		JetCollection uncleanedJets;
		JetCollection selectedJets;
		uncleanedJets.SetJetCollection(eve->jet_pt_[treeSys],eve->jet_eta_[treeSys],eve->jet_phi_[treeSys],eve->jet_m_[treeSys],eve->jet_combinedInclusiveSecondaryVertexV2BJetTags_[treeSys]);
		yGGhelper.SelectJets(uncleanedJets,selectedJets,selectedLooseleptons, Yggdrasil_Threshold_Cuts);
		  
		
		neve->jet_pt_[numSys] = selectedJets.jet_pt;
		neve->jet_eta_[numSys] = selectedJets.jet_eta;
		neve->jet_phi_[numSys] = selectedJets.jet_phi;
		neve->jet_m_[numSys] = selectedJets.jet_m;
		neve->jet_combinedInclusiveSecondaryVertexV2BJetTags_[numSys] = selectedJets.jet_bDiscriminant;
		neve->jet_flavour_[numSys] = eve->jet_flavour_[numSys];
		
		neve->MET_[numSys] = eve->MET_[numSys];
		neve->MET_phi_[numSys] = eve->MET_phi_[numSys];
		
		
		//End Systematic Loop;
		
		neve->wgt_generator_ = eve->wgt_generator_;
		neve->GoodFirstPV_ = eve->GoodFirstPV_;
		
		neve->passHLT_Ele27_eta2p1_WPTight_Gsf_v_ = eve->passHLT_Ele27_eta2p1_WPTight_Gsf_v_;
		neve->passHLT_IsoMu22_v_ = eve->passHLT_IsoMu22_v_;
		neve->passHLT_IsoTkMu22_v_ = eve->passHLT_IsoTkMu22_v_;
		neve->passHLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_ = eve->passHLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_;
		neve->passHLT_IsoMu24_eta2p1_v_ = eve->passHLT_IsoMu24_eta2p1_v_;

	  
	  YggdrasilSlim->Fill();
  	}

   	newtreefile->Write();
   	newtreefile->Close();

	

	 

}

