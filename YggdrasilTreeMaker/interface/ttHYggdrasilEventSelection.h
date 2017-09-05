
#ifndef TTHYGGDRASILEVENTSELECTION
#define TTHYGGDRASILEVENTSELECTION

#include <vector>
#include <TLorentzVector.h>

class ttHYggdrasilEventSelection{

 public :

  ttHYggdrasilEventSelection();
  ~ttHYggdrasilEventSelection();

  void SetElTrigger( const int * trigFlag );
  void SetMuTrigger( const int * trigFlag );

  void SetGoodVtx( const bool * _goodvtx );

  void SetElElTrigger( const int * trigFlag );
  void SetMuMuTrigger( const int * trigFlag );
  void SetElMuTrigger( const int * trigFlag );

  void SetEl_ORTrigger  ( const int n, ... );
  void SetMu_ORTrigger  ( const int n, ... );
  void SetElEl_ORTrigger( const int n, ... );
  void SetMuMu_ORTrigger( const int n, ... );
  void SetElMu_ORTrigger( const int n, ... );

  void SetLeptons( const std::vector<double> * pt, 
		   const std::vector<double> * eta, 
		   const std::vector<double> * scEta,
		   const std::vector<double> * phi,
		   const std::vector<double> * e,
		   const std::vector<int>    * charge,
		   const std::vector<int>    * isMuon, 
		   const std::vector<double> * relIso,
		   const std::vector<int> * POGLoose,
		   const std::vector<int> * POGTight );

  void SetLepTrigDR( const std::vector<double> * trig_dr );

  void SetJets( const std::vector<double> * pt, 
		const std::vector<double> * eta, 
		const std::vector<double> * phi, 
		const std::vector<double> * m,
		const std::vector<double> * bDiscriminant ,
		const std::vector<int> * flav );


  void SetJetsPTBeforRecorrection( const std::vector<double> * pt ,  const std::vector<double> * phi );

  void SetMet( const float * _met_pt, const float * _met_phi );

  void EnableInfoDumpForDebug();

  void doEventSelection();

  // ttH Tight lepton for SingleLepton
  std::vector<const TLorentzVector*> leptons();
  std::vector<double>                leptonsRelIso();
  std::vector<double>                leptonsSCEta();
  std::vector<int>                   leptonsIsMuon();
  std::vector<int>                   leptonsCharge();

  // ttH Loose leoton
  std::vector<const TLorentzVector*> looseLeptons();
  std::vector<double>                looseLeptonsRelIso();
  std::vector<double>                looseLeptonsScEta();
  std::vector<int>                   looseLeptonsIsMuon();
  std::vector<int>                   looseLeptonsCharge();

  // ttH Jet 
  std::vector<const TLorentzVector*> jets();
  std::vector<double> jetsBdiscriminant();
  std::vector<int> jetsFlav();

  // passed bTagging
  std::vector<const TLorentzVector*> bjets();
  std::vector<double> bjetsBdiscriminant();
  std::vector<int> bjetsFlav();

  bool PassSingleMuCh();
  bool PassSingleElCh();
  bool PassElEl();
  bool PassMuMu();
  bool PassElMu();

  float metAbs();
  float metPhi();

  // ** For DiLepton channel study **
  // ttH Tight lepton for DL
  std::vector<const TLorentzVector*> DLTightLeptons();
  std::vector<double>                DLTightLeptonsRelIso();
  std::vector<int>                   DLTightLeptonsIsMuon();
  std::vector<int>                   DLTightLeptonsCharge();

  // ** For DiLepton channel study **
  // ttH Jets LoosePTCut
  std::vector<const TLorentzVector*> DLSofterjets();
  std::vector<double>                DLSofterjetsBdiscriminant();

  // ** For DiLepton channel study **
  // ttH Jets LoosePTCut
  // passed bTagging
  std::vector<const TLorentzVector*> DLSofterbjets();
  std::vector<double> DLSofterbjetsBdiscriminant();

  bool passElTrig() ; 	
  bool passMuTrig() ; 	
  bool passElElTrig() ; 
  bool passMuMuTrig() ; 
  bool passElMuTrig() ; 

  double getLeptonDR( unsigned int idx );

  int getnNonIsoEl();
  int getnNonIsoMu() ;

  void SetFakeEstimationModeON();

  void SetTightMuonRelIsoCut( double cut );
  void SetLooseMuonRelIsoCut( double cut );

  double MinDRLepJet();


 private :

  bool b_InfoDumpForDebug;

  bool b_FakeEstimation;

  void _InitInternalVariables();

  void _ElectronSelection();
  void _MuonSelection();
  void _JetSelection();
  void _SortChargedLepton();
  void _SortChargedLepton( std::vector<const TLorentzVector*> * v_TLV ,
			   std::vector<double>                * v_iso ,
			   std::vector<double>                * v_scEta ,
			   std::vector<int>                   * v_isMuon ,
			   std::vector<int>                   * v_charge );
  bool _OverlapWithLooseLeptons( double eta, double phi);
  double _calcDR2( double eta1, double eta2, double phi1, double phi2 );

  void _MetCorrection();

  std::vector< const int * > ElTrig ; 
  std::vector< const int * > MuTrig ; 
  std::vector< const int * > ElElTrig ; 
  std::vector< const int * > MuMuTrig ; 
  std::vector< const int * > ElMuTrig ; 

  bool _passTrig( const std::vector< const int * > & triggerSets );

  double Thre_TightMu_PT ;
  double Thre_TightMu_Eta ;
  double Thre_TightMu_Iso ;

  double Thre_LooseMu_PT ;
  double Thre_LooseMu_Eta ;
  double Thre_LooseMu_Iso ;

  double Thre_TightEl_PT ;
  double Thre_TightEl_Eta ;
  double Thre_TightEl_Iso[2] ;

  double Thre_LooseEl_PT ;
  double Thre_LooseEl_Eta ;
  double Thre_LooseEl_Iso[2] ;

  double Thre_Jet_PT ;
  double Thre_Jet_Eta ;
  double Thre_Jet_Btag ;


  const std::vector<double> * lep_pt;
  const std::vector<double> * lep_eta; 
  const std::vector<double> * lep_phi;
  const std::vector<double> * lep_e;
  const std::vector<double> * lep_scEta;
  const std::vector<int>    * lep_charge;
  const std::vector<int>    * lep_isMuon; 
  const std::vector<double> * lep_relIso;
  const std::vector<int>    * lep_POGLoose;
  const std::vector<int>    * lep_POGTight;

  const std::vector<double> * lep_trigDR ; // extra info.

  const std::vector<double> * jet_pt; 
  const std::vector<double> * jet_eta; 
  const std::vector<double> * jet_phi; 
  const std::vector<double> * jet_m;
  const std::vector<double> * jet_bDiscriminant ;
  const std::vector<int>    * jet_flav; 

  const std::vector<double> * jet_pt_beforeRecorrection; 
  const std::vector<double> * jet_phi_beforeRecorrection; 

  const float * met_preRecorrection_pt , *met_preRecorrection_phi ;
  const bool * goodvtx;

  std::vector<const TLorentzVector*> selected_tightLeptons;
  std::vector<double>                selected_tightLeptonsRelIso;
  std::vector<int>                   selected_tightLeptonsIsMuon;
  std::vector<int>                   selected_tightLeptonsCharge;
  std::vector<double>                selected_tightLeptonsScEta;

  std::vector<const TLorentzVector*> selected_looseLeptons;
  std::vector<double>                selected_looseLeptonsRelIso;
  std::vector<double>                selected_looseLeptonsScEta;
  std::vector<int>                   selected_looseLeptonsIsMuon;
  std::vector<int>                   selected_looseLeptonsCharge;

  std::vector<const TLorentzVector*> selected_jets;
  std::vector<double>                selected_jetsBdiscriminant;
  std::vector<int>                selected_jetsFlav;
  std::vector<double>                selected_jetPt_preRecorrection;
  std::vector<double>                selected_jetPhi_preRecorrection;

  std::vector<const TLorentzVector*> selected_bjets;
  std::vector<double>                selected_bjetsBdiscriminant;
  std::vector<int>                selected_bjetsFlav;

  float met_pt , met_phi ;

  // ** for DiLepton channel Study **
  std::vector<const TLorentzVector*> DLselected_tightLeptons;
  std::vector<double>                DLselected_tightLeptonsRelIso;
  std::vector<double>                DLselected_tightLeptonsScEta;
  std::vector<int>                   DLselected_tightLeptonsIsMuon;
  std::vector<int>                   DLselected_tightLeptonsCharge;

  // ** for DiLepton channel Study **
  std::vector<const TLorentzVector*> DLsofterselected_jets;
  std::vector<double>                DLsofterselected_jetsBdiscriminant;
  std::vector<int>                DLsofterselected_jetsFlav;

  // ** for DiLepton channel Study **
  std::vector<const TLorentzVector*> DLsofterselected_bjets;
  std::vector<double>                DLsofterselected_bjetsBdiscriminant;
  std::vector<int>                DLsofterselected_bjetsFlav;


  std::vector<long> nEvent_passSingleMuCh;
  std::vector<long> nEvent_passSingleElCh;
  std::vector<long> nEvent_passElElCh;
  std::vector<long> nEvent_passMuMuCh;
  std::vector<long> nEvent_passElMuCh;

  int    _GetOriginalIdxOfTightLepton( unsigned int idx );

  int nNonIsolatedElectron ;
  int nNonIsolatedMuon ;

};



#endif
