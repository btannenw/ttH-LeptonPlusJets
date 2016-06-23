
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

  void SetLeptons( const std::vector<double> * pt, 
		   const std::vector<double> * eta, 
		   const std::vector<double> * phi,
		   const std::vector<double> * e,
		   const std::vector<int>    * charge,
		   const std::vector<int>    * isMuon, 
		   const std::vector<double> * relIso,
		   const std::vector<int> * POGLoose,
		   const std::vector<int> * POGTight );

  void SetJets( const std::vector<double> * pt, 
		const std::vector<double> * eta, 
		const std::vector<double> * phi, 
		const std::vector<double> * m,
		const std::vector<double> * bDiscriminant );

  void SetMet( const float * _met_pt, const float * _met_phi );

  void doEventSelection();

  // ttH Tight lepton for SingleLepton
  std::vector<const TLorentzVector*> leptons();
  std::vector<double>                leptonsRelIso();
  std::vector<int>                   leptonsIsMuon();
  std::vector<int>                   leptonsCharge();

  // ttH Loose leoton
  std::vector<const TLorentzVector*> looseLeptons();
  std::vector<double>                looseLeptonsRelIso();
  std::vector<int>                   looseLeptonsIsMuon();
  std::vector<int>                   looseLeptonsCharge();

  // ttH Jet 
  std::vector<const TLorentzVector*> jets();
  std::vector<double> jetsBdiscriminant();

  // passed bTagging
  std::vector<const TLorentzVector*> bjets();
  std::vector<double> bjetsBdiscriminant();

  bool PassSingleMuCh();
  bool PassSingleElCh();
  bool PassElEl();
  bool PassMuMu();
  bool PassElMu();


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


 private :
   
  void _InitInternalVariables();

  void _ElectronSelection();
  void _MuonSelection();
  void _JetSelection();
  void _SortChargedLepton();
  void _SortChargedLepton( std::vector<const TLorentzVector*> * v_TLV ,
			   std::vector<double>                * v_iso ,
			   std::vector<int>                   * v_isMuon ,
			   std::vector<int>                   * v_charge );
  bool _OverlapWithLooseLeptons( double eta, double phi);
  double _calcDR2( double eta1, double eta2, double phi1, double phi2 );


  const int * ElTrig ; 
  const int * MuTrig ; 
  const int * ElElTrig ; 
  const int * MuMuTrig ; 
  const int * ElMuTrig ; 

  double Thre_TightMu_PT ;
  double Thre_TightMu_Eta ;
  double Thre_TightMu_Iso ;

  double Thre_LooseMu_PT ;
  double Thre_LooseMu_Eta ;
  double Thre_LooseMu_Iso ;

  double Thre_TightEl_PT ;
  double Thre_TightEl_Eta ;
  double Thre_TightEl_Iso ;

  double Thre_LooseEl_PT ;
  double Thre_LooseEl_Eta ;
  double Thre_LooseEl_Iso ;

  double Thre_Jet_PT ;
  double Thre_Jet_Eta ;
  double Thre_Jet_Btag ;


  const std::vector<double> * lep_pt;
  const std::vector<double> * lep_eta; 
  const std::vector<double> * lep_phi;
  const std::vector<double> * lep_e;
  const std::vector<int>    * lep_charge;
  const std::vector<int>    * lep_isMuon; 
  const std::vector<double> * lep_relIso;
  const std::vector<int>    * lep_POGLoose;
  const std::vector<int>    * lep_POGTight;

  const std::vector<double> * jet_pt; 
  const std::vector<double> * jet_eta; 
  const std::vector<double> * jet_phi; 
  const std::vector<double> * jet_m;
  const std::vector<double> * jet_bDiscriminant ;

  const float * met_pt , *met_phi ;
  const bool * goodvtx;

  std::vector<const TLorentzVector*> selected_tightLeptons;
  std::vector<double>                selected_tightLeptonsRelIso;
  std::vector<int>                   selected_tightLeptonsIsMuon;
  std::vector<int>                   selected_tightLeptonsCharge;

  std::vector<const TLorentzVector*> selected_looseLeptons;
  std::vector<double>                selected_looseLeptonsRelIso;
  std::vector<int>                   selected_looseLeptonsIsMuon;
  std::vector<int>                   selected_looseLeptonsCharge;

  std::vector<const TLorentzVector*> selected_jets;
  std::vector<double>                selected_jetsBdiscriminant;

  std::vector<const TLorentzVector*> selected_bjets;
  std::vector<double>                selected_bjetsBdiscriminant;


  // ** for DiLepton channel Study **
  std::vector<const TLorentzVector*> DLselected_tightLeptons;
  std::vector<double>                DLselected_tightLeptonsRelIso;
  std::vector<int>                   DLselected_tightLeptonsIsMuon;
  std::vector<int>                   DLselected_tightLeptonsCharge;

  // ** for DiLepton channel Study **
  std::vector<const TLorentzVector*> DLsofterselected_jets;
  std::vector<double>                DLsofterselected_jetsBdiscriminant;

  // ** for DiLepton channel Study **
  std::vector<const TLorentzVector*> DLsofterselected_bjets;
  std::vector<double>                DLsofterselected_bjetsBdiscriminant;


};



#endif