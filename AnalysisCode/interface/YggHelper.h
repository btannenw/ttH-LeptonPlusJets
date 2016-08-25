#ifndef _YggHelper_h
#define _YggHelper_h

#include <iostream>
#include <vector>
#include <map>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <algorithm>
#include "TVector.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2D.h"


#ifdef __MAKECINT__
#pragma link C++ class std::vector< TLorentzVector >+; 
#endif

#if !defined(__CINT__) && !defined(__MAKECINT__)


#endif

typedef std::map<std::string, std::string> mparams;
typedef std::vector< TLorentzVector > vecTLorentzVector;
typedef std::vector<std::vector<double> > vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<std::vector<int> > vvint;
typedef std::vector<std::string> vstring;
typedef std::vector<double> vdouble;
typedef std::vector<int> vint;
typedef std::vector<float> vfloat;

struct ThresholdList {
	double Thre_TightMu_PT;
	double Thre_LooseMu_PT;
	double Thre_TightMu_Eta;
	double Thre_LooseMu_Eta;
	double Thre_TightMu_Iso;
	double Thre_LooseMu_Iso;
	
	double Thre_TightEl_PT;
	double Thre_LooseEl_PT;
	double Thre_TightEl_Eta;
	double Thre_LooseEl_Eta;
	double Thre_TightEl_Iso;
	double Thre_LooseEl_Iso;
	
	double Thre_Jet_PT;
	double Thre_Jet_Eta;
	double Thre_Jet_Btag;
};

struct JetCollection{

	std::vector<double> jet_pt;
	std::vector<double> jet_eta;
	std::vector<double> jet_phi;
	std::vector<double> jet_m;
	std::vector<double> jet_bDiscriminant;
	
	void SetJetCollection(vdouble j_pt, vdouble j_eta, vdouble j_phi, vdouble j_m, vdouble j_csv);
	void AddJet(JetCollection JetCollectionB, int iJet);
};

struct LeptonCollection {
			std::vector<double> lep_pt;
			std::vector<double> lep_eta;
			std::vector<double> lep_phi;
			std::vector<double> lep_e;
			std::vector<int> lep_charge;
			std::vector<int> lep_isMuon;
			std::vector<double> lep_relIso;
			std::vector<int> lep_POGLoose;
			std::vector<int> lep_POGTight;
			
			void SetLeptonCollection(vdouble l_pt, vdouble l_eta, vdouble l_phi, vdouble l_e, vint l_charge, vint l_isMuon, vdouble l_relIso, vint l_isTight, vint l_isLoose);
			void AddLep(LeptonCollection LeptonCollectionB,int iLep);
			void ClearLep();
			void SortLeptons();
};
		





using namespace std;

class YggHelper{

	// === Functions === //
	public: 
		// Constructor(s) and destructor
		YggHelper();
		virtual ~YggHelper();
		
		// Set up BDTvars
   
		void Test();
	

	//Algorithms 
		
		
		
		void remakeJetV(std::vector<TLorentzVector> JetV, vdouble csvV, double JetPTcut, double JetEtacut, std::vector<TLorentzVector> &NewJetV, vdouble &NewJetcsvV, int &NewnumJets, int &NewtagJets );
		void remakeLeptonV(std::vector<TLorentzVector> LeptonV, double eleCut, double muonCut, std::vector<TLorentzVector> &NewLeptonV);
		
		int passCuts(bool goodfirstPV, int ElectronTrigger, int MuonTrigger, vint lepnums, int numJet, int numTag);
		void fillCSVhistos(TFile *fileHF, TFile *fileLF);
		double get_csv_wgt( std::vector<TLorentzVector>, vdouble jetCSV, vint jetFlavor, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );
		void FillLeptonWgtHistos(TFile *fileEleTrig, TFile *fileEleID, TFile *fileEleIso, TFile *fileMuonTrig, TFile *fileMuonID, TFile *fileMuonIso);
		double getLeptonWgt(vint lepnums, std::vector<TLorentzVector> vlepton);
		double reweightPU( int nPU, int iSys );
		
		void SetCuts(ThresholdList &cuts, double tightmupt, double loosemupt, double tightmueta, double loosemueta, double tightmuiso, double loosemuiso, double tightelpt, double looseelpt, double tighteleta, double looseeleta, double tighteliso, double looseeliso, double jetpt, double jeteta, double jetbtag);
		void SelectLeptons(LeptonCollection inputLeptons, LeptonCollection &outputTightLeptons, LeptonCollection &outputLooseLeptons, ThresholdList cuts);
		void ElectronSelection(LeptonCollection inputLeptons, LeptonCollection &outputTightLeptons, LeptonCollection &outputLooseLeptons, ThresholdList cuts);
		void MuonSelection(LeptonCollection inputLeptons, LeptonCollection &outputTightLeptons, LeptonCollection &outputLooseLeptons, ThresholdList cuts);
		void SelectJets(JetCollection inputJets, JetCollection &outputJets, LeptonCollection inputLooseLeptons, ThresholdList cuts);
		bool OverlapWithLooseLeptons(double jeta, double jphi, LeptonCollection inputLoose);
		void makeVTLVjets(std::vector<TLorentzVector> &NewJetV, JetCollection inputJets, int &NewNumJets, int &NewNumTags, ThresholdList cuts);		
		void makeVTLVleptons(std::vector<TLorentzVector> &NewLeptonV, vint &lepnums, LeptonCollection inputTightLeptons, LeptonCollection inputLooseLeptons);
		
	
		
	
	

	// === Variables === //
	public:

	protected:
	    
		TH1D* h_csv_wgt_hf[9][6];
		TH1D* c_csv_wgt_hf[9][6];
		TH1D* h_csv_wgt_lf[9][4][3];
		
		TH2D* h_lep_wgt_ele_trig;
		TH2D* h_lep_wgt_ele_iso;
		TH2D* h_lep_wgt_ele_id;
		TH2D* h_lep_wgt_mu_trig1;
		TH2D* h_lep_wgt_mu_trig2;
		TH2D* h_lep_wgt_mu_iso;
		TH2D* h_lep_wgt_mu_id;
		
		double csvwp;
		

	private:

}; // End of class prototype






#endif // _YggHelper_h







