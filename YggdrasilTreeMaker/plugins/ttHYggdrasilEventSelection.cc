
#ifdef STANDALONECOMPILE
#include "./ttHYggdrasilEventSelection.h"
#else
#include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilEventSelection.h"
#endif

#include <iostream>
#include <stdarg.h>

ttHYggdrasilEventSelection::ttHYggdrasilEventSelection(){

 lep_trigDR = 0 ;

 Thre_TightMu_PT  = 26.0 ;
 Thre_LooseMu_PT  = 15.0;

 Thre_TightMu_Eta = 2.1 ;
 Thre_LooseMu_Eta = 2.4 ;

 Thre_TightMu_Iso = 0.15 ;
 Thre_LooseMu_Iso = 0.25 ;

 Thre_TightEl_PT  = 30;
 Thre_LooseEl_PT  = 15;

 Thre_TightEl_Eta = 2.1;
 Thre_LooseEl_Eta = 2.4;

 Thre_TightEl_Iso = 0.15;
 Thre_LooseEl_Iso = 0.15;

 Thre_Jet_PT      = 30.0;
 Thre_Jet_Eta     =  2.4;

 Thre_Jet_Btag    =  0.8484;

 b_InfoDumpForDebug = false ;

 b_FakeEstimation = false ;

 for( int i = 0 ; i < 20 ; i ++ ){
   nEvent_passSingleMuCh .push_back(0);
   nEvent_passSingleElCh .push_back(0);
   nEvent_passElElCh .push_back(0);
   nEvent_passMuMuCh .push_back(0);
   nEvent_passElMuCh .push_back(0);
 }

}


ttHYggdrasilEventSelection::~ttHYggdrasilEventSelection(){


  if( b_InfoDumpForDebug ){

    std::cout <<" * * * Cut flow : Mu ch * * *"  << std::endl ; 
    for( int i = 0 ; i < 20 && nEvent_passSingleMuCh[i] !=0 ; i ++ ){
      std::cout << "Cut_" << i << " " << nEvent_passSingleMuCh [i]<< std::endl;
    }
    std::cout <<" * * * Cut flow : El ch * * *"  << std::endl ; 
    for( int i = 0 ; i < 20 && nEvent_passSingleElCh[i]!=0 ; i ++ ){
      std::cout << "Cut_" << i << " " << nEvent_passSingleElCh[i]<< std::endl ;
    }
    std::cout <<" * * * Cut flow : ElEl ch * * *"  << std::endl ; 
    for( int i = 0 ; i < 20 && nEvent_passElElCh[i]!=0 ; i ++ ){
      std::cout << "Cut_" << i << " " << nEvent_passElElCh[i] << std::endl ;
    }
    std::cout <<" * * * Cut flow : MuMu ch * * *"  << std::endl ; 
    for( int i = 0 ; i < 20 && nEvent_passMuMuCh[i]!=0 ; i ++ ){
      std::cout << "Cut_" << i << " " << nEvent_passMuMuCh[i] << std::endl ;
    }
    std::cout <<" * * * Cut flow : ElMu ch * * *"  << std::endl ; 
    for( int i = 0 ; i < 20 && nEvent_passElMuCh[i]!=0 ; i ++ ){
      std::cout << "Cut_" << i << " " << nEvent_passElMuCh[i]  << std::endl ;
    }
    
  }

  _InitInternalVariables();

}

void ttHYggdrasilEventSelection::SetElTrigger( const int* trigFlag){

  ElTrig .clear();
  ElTrig .push_back( trigFlag ) ;

}

void ttHYggdrasilEventSelection::SetMuTrigger( const int * trigFlag){

  MuTrig . clear();
  MuTrig . push_back( trigFlag );

}

void ttHYggdrasilEventSelection::SetElElTrigger( const int * trigFlag){

  ElElTrig .clear();
  ElElTrig .push_back( trigFlag );

}
void ttHYggdrasilEventSelection::SetMuMuTrigger( const int * trigFlag){

  MuMuTrig . clear() ;
  MuMuTrig .push_back( trigFlag );

}
void ttHYggdrasilEventSelection::SetElMuTrigger( const int * trigFlag){

  ElMuTrig . clear();
  ElMuTrig .push_back( trigFlag );

}


void ttHYggdrasilEventSelection::SetEl_ORTrigger( int n , ... ){

  ElTrig .clear();
  va_list args;
  va_start(args , n );
  for(int i = 0 ; i < n ; i ++){
    ElTrig .push_back(  va_arg(args , const int *) ) ;
  }
  va_end(args);
  return ; 
}

void ttHYggdrasilEventSelection::SetMu_ORTrigger( int n , ... ){


  MuTrig .clear();
  va_list args;
  va_start(args , n);
  for(int i = 0 ; i < n ; i ++){
    MuTrig .push_back(  va_arg(args , const int *) ) ;
  }
  va_end(args);
  return ; 

}

void ttHYggdrasilEventSelection::SetElEl_ORTrigger( int n , ... ){


  ElElTrig .clear();
  va_list args;
  va_start(args , n);
  for(int i = 0 ; i < n ; i ++){
    ElElTrig .push_back(  va_arg(args , const int *) ) ;
  }
  va_end(args);
  return ; 


}
void ttHYggdrasilEventSelection::SetMuMu_ORTrigger( int n , ... ){

  
  MuMuTrig .clear();
  va_list args;
  va_start(args , n);
  for(int i = 0 ; i < n ; i ++){
    MuMuTrig .push_back(  va_arg(args , const int *) ) ;
  }
  va_end(args);
  return ; 


}
void ttHYggdrasilEventSelection::SetElMu_ORTrigger( int n , ... ){

  ElMuTrig .clear();
  va_list args;
  va_start(args , n);
  for(int i = 0 ; i < n ; i ++){
    ElMuTrig .push_back(  va_arg(args , const int *) ) ;
  }
  va_end(args);
  return ; 

}


void ttHYggdrasilEventSelection::SetGoodVtx( const bool * _goodvtx ){
  goodvtx = _goodvtx ;
}

void ttHYggdrasilEventSelection::SetLeptons( const std::vector<double> * pt, 
					     const std::vector<double> * eta, 
					     const std::vector<double> * scEta,
					     const std::vector<double> * phi,
					     const std::vector<double> * e,
					     const std::vector<int>    * charge, 
					     const std::vector<int>    * isMuon, 
					     const std::vector<double> * relIso,
					     const std::vector<int> * POGLoose,
					     const std::vector<int> * POGTight ){
  
  lep_pt  = pt ;  
  lep_eta = eta ;  
  lep_scEta = scEta ;  
  lep_phi = phi ;
  lep_e   = e ;
  lep_charge = charge ;
  lep_isMuon = isMuon ; 
  lep_relIso = relIso ; 
  lep_POGLoose = POGLoose ; 
  lep_POGTight = POGTight ;

}

void ttHYggdrasilEventSelection::SetLepTrigDR( const std::vector<double> * trig_dr ){
  lep_trigDR = trig_dr ;
} 


void ttHYggdrasilEventSelection::SetJets( const std::vector<double> * pt, 
					  const std::vector<double> * eta, 
					  const std::vector<double> * phi, 
					  const std::vector<double> * m,
					  const std::vector<double> * bDiscriminant ,
					  const std::vector<int> * flav ){

  jet_pt = pt ; 
  jet_eta = eta ; 
  jet_phi = phi ;  
  jet_m = m ;
  jet_bDiscriminant = bDiscriminant ;
  jet_flav = flav;
}

void ttHYggdrasilEventSelection::SetMet( const float * _met_pt, const float * _met_phi ){

  met_pt  = _met_pt ;
  met_phi = _met_phi ;

}

void ttHYggdrasilEventSelection::doEventSelection(){

  _InitInternalVariables();

  // Object Definition 
  _ElectronSelection();
  _MuonSelection();
  _SortChargedLepton();
  _JetSelection();

  
}

void  ttHYggdrasilEventSelection::_SortChargedLepton(){

  _SortChargedLepton( & selected_tightLeptons,
		      & selected_tightLeptonsRelIso,
		      & selected_tightLeptonsScEta,
		      & selected_tightLeptonsIsMuon,
		      & selected_tightLeptonsCharge );

  _SortChargedLepton( & selected_looseLeptons,
		      & selected_looseLeptonsRelIso,
		      & selected_looseLeptonsScEta,
		      & selected_looseLeptonsIsMuon,
		      & selected_looseLeptonsCharge );

  _SortChargedLepton( & DLselected_tightLeptons,
		      & DLselected_tightLeptonsRelIso,
		      & DLselected_tightLeptonsScEta,
		      & DLselected_tightLeptonsIsMuon,
		      & DLselected_tightLeptonsCharge
		      );

}


void  ttHYggdrasilEventSelection::_SortChargedLepton( std::vector<const TLorentzVector*> * v_TLV ,
						      std::vector<double>                * v_iso ,
						      std::vector<double>                * v_scEta ,
						      std::vector<int>                   * v_isMuon ,
						      std::vector<int>                   * v_charge ){


  for( unsigned int i = 0 ; i < v_TLV -> size() ; i++  ){
    for( unsigned int j = i+1 ; j < v_TLV -> size() ; j++  ){
      
      if( v_TLV->at(i)->Pt() > v_TLV->at(j)->Pt() ){ continue ; }

      // Perform swap
      const TLorentzVector* swap_TLV    = v_TLV    -> at(i);
      double                swap_iso    = v_iso    -> at(i);
      double                swap_scEta  = v_scEta  -> at(i);
      int                   swap_isMuon = v_isMuon -> at(i);
      int                   swap_charge = v_charge -> at(i);

       v_TLV    -> at(i)  = v_TLV    -> at(j);
       v_iso    -> at(i)  = v_iso    -> at(j);
       v_scEta  -> at(i)  = v_scEta  -> at(j);
       v_isMuon -> at(i)  = v_isMuon -> at(j);
       v_charge -> at(i)  = v_charge -> at(j);

       v_TLV    -> at(j)  = swap_TLV   ;
       v_iso    -> at(j)  = swap_iso   ;
       v_scEta  -> at(j)  = swap_scEta ;
       v_isMuon -> at(j)  = swap_isMuon;
       v_charge -> at(j)  = swap_charge;

    }
  }

}

void ttHYggdrasilEventSelection::_InitInternalVariables(){

  // Clear pointers
  // ( Remember TLorentzVector* of TightLeptons are included in the Loose,
  //   and does not need "delete" for them.
  //   Similarly, bjet are in jet.)

  for( std::vector<const TLorentzVector*>::iterator vec = selected_looseLeptons . begin();
       vec != selected_looseLeptons.end();
       vec ++
       ){
    if( (*vec) == 0 ) continue ;
    delete (*vec);
    (*vec) = 0 ;
  }

  for( std::vector<const TLorentzVector*>::iterator vec = selected_jets . begin();
       vec != selected_jets.end();
       vec ++
       ){
    if( (*vec) == 0 ) continue ;
    delete (*vec);
    (*vec) = 0 ;
  }


  // ** For DiLepton channel study **
  for( std::vector<const TLorentzVector*>::iterator vec = DLselected_tightLeptons . begin();
       vec != DLselected_tightLeptons.end();
       vec ++
       ){
    if( (*vec) == 0 ) continue ;
    delete (*vec);
    (*vec) = 0 ;
  }

  // ** For DiLepton channel study **
  for( std::vector<const TLorentzVector*>::iterator vec = DLsofterselected_jets . begin();
       vec != DLsofterselected_jets.end();
       vec ++
       ){
    if( (*vec) == 0 ) continue ;
    delete (*vec);
    (*vec) = 0 ;
  }

  selected_tightLeptons.clear();
  selected_tightLeptonsRelIso.clear();
  selected_tightLeptonsScEta.clear();
  selected_tightLeptonsIsMuon.clear();
  selected_tightLeptonsCharge.clear();
  
  selected_looseLeptons.clear();
  selected_looseLeptonsRelIso.clear();
  selected_looseLeptonsScEta.clear();
  selected_looseLeptonsIsMuon.clear();
  selected_looseLeptonsCharge.clear();
  
  selected_jets.clear();
  selected_jetsBdiscriminant.clear();
  selected_jetsFlav.clear();
  
  selected_bjets.clear();
  selected_bjetsBdiscriminant.clear();
  selected_bjetsFlav.clear();

  // ** For DiLepton channel study **
  DLselected_tightLeptons.clear();
  DLselected_tightLeptonsRelIso.clear();
  DLselected_tightLeptonsScEta.clear();
  DLselected_tightLeptonsIsMuon.clear();
  DLselected_tightLeptonsCharge.clear();

  // ** For DiLepton channel study **  
  DLsofterselected_jets.clear();
  DLsofterselected_jetsBdiscriminant.clear();
  DLsofterselected_jetsFlav.clear();

  // ** For DiLepton channel study **  
  DLsofterselected_bjets.clear();
  DLsofterselected_bjetsBdiscriminant.clear();
  DLsofterselected_bjetsFlav.clear();

}

void ttHYggdrasilEventSelection::_ElectronSelection(){

  for( unsigned int idx = 0 ; idx < lep_pt->size() ; idx++ ){

    if( lep_isMuon -> at(idx) != 0 ) continue ;
    
    if(     lep_POGTight -> at(idx)  != 1  ) continue ;
    if(        lep_pt    -> at(idx)   < Thre_LooseEl_PT  ) continue ;
    if( fabs( lep_eta    -> at(idx) ) > Thre_LooseEl_Eta ) continue ;
    if(       lep_relIso -> at(idx)   > Thre_LooseEl_Iso ) continue ;

    TLorentzVector * vec = new TLorentzVector;
    vec->SetPtEtaPhiE( lep_pt  -> at(idx),
		       lep_eta -> at(idx),
		       lep_phi -> at(idx),
		       lep_e   -> at(idx) );

    selected_looseLeptons.push_back( vec );
    selected_looseLeptonsRelIso.push_back( lep_relIso -> at(idx) );
    selected_looseLeptonsScEta.push_back( lep_scEta -> at(idx) );
    selected_looseLeptonsIsMuon.push_back( 0 );
    selected_looseLeptonsCharge.push_back( lep_charge -> at(idx) );


  }

  nNonIsolatedElectron = 0 ; 

  for( unsigned int idx = 0 ; idx < lep_pt->size() ; idx++ ){
    
    if( lep_isMuon -> at(idx) != 0 ) continue ;
    if(        lep_pt    -> at(idx)   < Thre_TightEl_PT  ) continue ;
    if( fabs( lep_eta    -> at(idx) ) > Thre_TightEl_Eta ) continue ;

    if(     lep_POGLoose -> at(idx)  != 1  ) continue ;
    nNonIsolatedElectron ++ ;

    if( ! b_FakeEstimation  &&  ( lep_POGTight -> at(idx)  != 1  ||  fabs( lep_relIso -> at(idx) ) > Thre_TightEl_Iso ) ) continue ;
    if(   b_FakeEstimation  &&    lep_POGTight -> at(idx)  == 1  ) continue ;
    // Note : electron ID requires relIso in ID, and no need to check iso after requiring tight!=1.
	  
    TLorentzVector * vec = new TLorentzVector;
    vec->SetPtEtaPhiE( lep_pt  -> at(idx),
		       lep_eta -> at(idx),
		       lep_phi -> at(idx),
		       lep_e   -> at(idx) );

    selected_tightLeptons.push_back( vec );
    selected_tightLeptonsRelIso.push_back( lep_relIso -> at(idx) );
    selected_tightLeptonsScEta.push_back( lep_scEta -> at(idx) );
    selected_tightLeptonsIsMuon.push_back( 0 );
    selected_tightLeptonsCharge.push_back( lep_charge -> at(idx) );

  }


  // Independent loop for DiLepton Tight Lepton.
  for( unsigned int idx = 0 ; idx < lep_pt->size() ; idx++ ){

    if( lep_isMuon -> at(idx) != 0 ) continue ;
    
    if(     lep_POGTight -> at(idx)  != 1     ) continue ;
    if(        lep_pt    -> at(idx)   < 25.0  ) continue ;
    if( fabs( lep_eta    -> at(idx) ) >  2.4  ) continue ;
    if(       lep_relIso -> at(idx)   >  0.15 ) continue ;

    TLorentzVector * vec = new TLorentzVector;
    vec->SetPtEtaPhiE( lep_pt  -> at(idx),
		       lep_eta -> at(idx),
		       lep_phi -> at(idx),
		       lep_e   -> at(idx) );

    DLselected_tightLeptons.push_back( vec );
    DLselected_tightLeptonsRelIso.push_back( lep_relIso -> at(idx) );
    DLselected_tightLeptonsScEta.push_back( lep_scEta -> at(idx) );
    DLselected_tightLeptonsIsMuon.push_back( 0 );
    DLselected_tightLeptonsCharge.push_back( lep_charge -> at(idx) );

  }


}

void ttHYggdrasilEventSelection::_MuonSelection(){

  for( unsigned int idx = 0 ; idx < lep_pt->size() ; idx++ ){

    if( lep_isMuon -> at(idx) != 1 ) continue ;
    
    if(     lep_POGTight -> at(idx)  != 1  ) continue ;
    if(        lep_pt    -> at(idx)   < Thre_LooseMu_PT  ) continue ;
    if( fabs( lep_eta    -> at(idx) ) > Thre_LooseMu_Eta ) continue ;
    if(       lep_relIso -> at(idx)   > Thre_LooseMu_Iso ) continue ;

    TLorentzVector * vec = new TLorentzVector;
    vec->SetPtEtaPhiE( lep_pt  -> at(idx),
		       lep_eta -> at(idx),
		       lep_phi -> at(idx),
		       lep_e   -> at(idx) );

    selected_looseLeptons.push_back( vec );
    selected_looseLeptonsRelIso.push_back( lep_relIso -> at(idx) );
    selected_looseLeptonsScEta.push_back( lep_scEta -> at(idx) );
    selected_looseLeptonsIsMuon.push_back( 1 );
    selected_looseLeptonsCharge.push_back( lep_charge -> at(idx) );
  }

  nNonIsolatedMuon = 0 ; 

  for( unsigned int idx = 0 ; idx < lep_pt->size() ; idx++ ){

    if( lep_isMuon -> at(idx) != 1 ) continue ;

    if(     lep_POGTight -> at(idx)  != 1  ) continue ;
    if(        lep_pt    -> at(idx)   < Thre_TightMu_PT  ) continue ;
    if( fabs( lep_eta    -> at(idx) ) > Thre_TightMu_Eta ) continue ;

    nNonIsolatedMuon ++ ; 
  
    if( ! b_FakeEstimation && fabs( lep_relIso -> at(idx) ) > Thre_TightMu_Iso ) continue ;
    if(   b_FakeEstimation && fabs( lep_relIso -> at(idx) ) < Thre_TightMu_Iso ) continue ;


    TLorentzVector * vec = new TLorentzVector;
    vec->SetPtEtaPhiE( lep_pt  -> at(idx),
		       lep_eta -> at(idx),
		       lep_phi -> at(idx),
		       lep_e   -> at(idx) );

    selected_tightLeptons.push_back( vec );
    selected_tightLeptonsRelIso.push_back( lep_relIso -> at(idx) );
    selected_tightLeptonsScEta.push_back( lep_scEta -> at(idx) );
    selected_tightLeptonsIsMuon.push_back( 1 );
    selected_tightLeptonsCharge.push_back( lep_charge -> at(idx) );

  }

  // Independent loop for DiLepton Tight Lepton.
  for( unsigned int idx = 0 ; idx < lep_pt->size() ; idx++ ){

    if( lep_isMuon -> at(idx) != 1 ) continue ;
    
    if(     lep_POGTight -> at(idx)  != 1     ) continue ;
    if(        lep_pt    -> at(idx)   < 25.0  ) continue ;
    if( fabs( lep_eta    -> at(idx) ) >  2.4  ) continue ;
    if(       lep_relIso -> at(idx)   >  0.25 ) continue ;

    TLorentzVector * vec = new TLorentzVector;
    vec->SetPtEtaPhiE( lep_pt  -> at(idx),
		       lep_eta -> at(idx),
		       lep_phi -> at(idx),
		       lep_e   -> at(idx) );

    DLselected_tightLeptons.push_back( vec );
    DLselected_tightLeptonsRelIso.push_back( lep_relIso -> at(idx) );
    DLselected_tightLeptonsScEta.push_back( lep_scEta -> at(idx) );
    DLselected_tightLeptonsIsMuon.push_back( 1 );
    DLselected_tightLeptonsCharge.push_back( lep_charge -> at(idx) );

  }

}

void ttHYggdrasilEventSelection::_JetSelection(){


  for( unsigned int idx = 0 ; idx < jet_pt->size() ; idx++ ){

    if(        jet_pt    -> at(idx)   < Thre_Jet_PT  ) continue ;
    if( fabs( jet_eta    -> at(idx) ) > Thre_Jet_Eta ) continue ;

    if( _OverlapWithLooseLeptons( jet_eta -> at(idx), jet_phi -> at(idx) ) ) continue ; 

    TLorentzVector * vec = new TLorentzVector;
    vec->SetPtEtaPhiM( jet_pt  -> at(idx),
		       jet_eta -> at(idx),
		       jet_phi -> at(idx),
		       jet_m   -> at(idx) );

    selected_jets.push_back( vec );
    selected_jetsBdiscriminant.push_back( jet_bDiscriminant ->at(idx) );
    selected_jetsFlav         .push_back( jet_flav ->at(idx) );

    if( jet_bDiscriminant ->at(idx) < Thre_Jet_Btag ) continue ;

    selected_bjets.push_back( vec );
    selected_bjetsBdiscriminant.push_back( jet_bDiscriminant ->at(idx) );
    selected_bjetsFlav.push_back( jet_flav ->at(idx) );
  }




  // ** independent loop for Di-Lepton channel **
  for( unsigned int idx = 0 ; idx < jet_pt->size() ; idx++ ){

    if(        jet_pt    -> at(idx)   < 20  ) continue ;
    if( fabs( jet_eta    -> at(idx) ) > 2.4 ) continue ;

    if( _OverlapWithLooseLeptons( jet_eta -> at(idx), jet_phi -> at(idx) ) ) continue ; 

    TLorentzVector * vec = new TLorentzVector;
    vec->SetPtEtaPhiM( jet_pt  -> at(idx),
		       jet_eta -> at(idx),
		       jet_phi -> at(idx),
		       jet_m   -> at(idx) );

    DLsofterselected_jets.push_back( vec );
    DLsofterselected_jetsBdiscriminant.push_back( jet_bDiscriminant ->at(idx) );
    DLsofterselected_jetsFlav.push_back( jet_flav ->at(idx) );

    if( jet_bDiscriminant ->at(idx) < Thre_Jet_Btag ) continue ;

    DLsofterselected_bjets.push_back( vec );
    DLsofterselected_bjetsBdiscriminant.push_back( jet_bDiscriminant ->at(idx) );
    DLsofterselected_bjetsFlav.push_back( jet_flav ->at(idx) );

  }

}

bool ttHYggdrasilEventSelection::_OverlapWithLooseLeptons( double eta, double phi){

  for( std::vector<const TLorentzVector*>::iterator lep = selected_looseLeptons.begin(); 
       lep != selected_looseLeptons.end();
       lep++){
    if( _calcDR2( eta, (*lep)->Eta(), phi, (*lep)->Phi() ) < 0.4 * 0.4 ){ return true ; }
  }
  return false ;
}


bool ttHYggdrasilEventSelection::PassSingleMuCh(){

  int i_step = 0 ;
  nEvent_passSingleMuCh[ i_step++ ]++ ; 

  // Trig Requirement.
  if( ! passMuTrig() ){ return false ;}
  nEvent_passSingleMuCh[ i_step++ ]++ ; 

  // (no vertex requirement inside the tool.)
  if( ! (*goodvtx) ){return false;}
  nEvent_passSingleMuCh[ i_step++ ]++ ; 

  // One tight muon.
  if( selected_tightLeptons.size() != 1 || selected_tightLeptonsIsMuon[0] != 1 ){ return false ; }
  nEvent_passSingleMuCh[ i_step++ ]++ ; 

  // No additional loose leptons.
  if( selected_looseLeptons.size() != 1 ){ return false ; } ;
  nEvent_passSingleMuCh[ i_step++ ]++ ; 

  // jet
  if( selected_jets . size() < 4 ){ return false ;}
  nEvent_passSingleMuCh[ i_step++ ]++ ; 

  // btagjet
  // jet
  if( selected_bjets . size() < 2 ){ return false ;}
  nEvent_passSingleMuCh[ i_step++ ]++ ; 

  return true; 
}

bool ttHYggdrasilEventSelection::PassSingleElCh(){

  int i_step = 0 ;
  nEvent_passSingleElCh[ i_step++ ]++ ; 

  // Trig Requirement.
  if( ! passElTrig() ){ return false ;}
  nEvent_passSingleElCh[ i_step++ ]++ ; 

  // (no vertex requirement inside the tool.)
  if( ! (*goodvtx) ){ return false;}
  nEvent_passSingleElCh[ i_step++ ]++ ; 

  // One tight muon.
  if( selected_tightLeptons.size() != 1 || selected_tightLeptonsIsMuon[0] == 1 ){ return false ; }
  nEvent_passSingleElCh[ i_step++ ]++ ; 

  // No additional loose leptons.
  if( selected_looseLeptons.size() != 1 ){ return false ; } ;
  nEvent_passSingleElCh[ i_step++ ]++ ; 

  // jet
  if( selected_jets . size() < 4 ){ return false ;}
  nEvent_passSingleElCh[ i_step++ ]++ ; 

  // btagjet
  // jet
  if( selected_bjets . size() < 2 ){ return false ;}
  nEvent_passSingleElCh[ i_step++ ]++ ; 

  return true; 

}

bool ttHYggdrasilEventSelection::PassElEl(){

  int i_step = 0 ;
  nEvent_passElElCh[ i_step++ ]++ ; 

  // Trigger 
  if( ! passElElTrig()  ){ return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 

  // GoodVtx
  if( ! (*goodvtx) ) { return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 

  // 2 OS-leptons (in DL definition) && FlavMatching with triggers.
  if( DLselected_tightLeptons.size() < 1 ){ return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 
  if(   selected_looseLeptons.size() < 2 ){ return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 
  if(   selected_looseLeptonsIsMuon.at(0) != 0 ){ return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 
  if(   selected_looseLeptonsIsMuon.at(1) != 0 ){ return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 

  // Veto additional leptons
  if(   selected_looseLeptons.size() != 2 ){ return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 

  const double invMass = ( (*selected_looseLeptons.at(0)) + (*selected_looseLeptons.at(1)) ).M() ;

  // mll>20 GeV.
  if( invMass < 20 ){ return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 

  // Exclude Z window (mll >76 GeV and <106 GeV) for ee and mumu.
  if( 76 < invMass &&  invMass < 106 ){ return false ;}
  nEvent_passElElCh[ i_step++ ]++ ; 

  // >=2 jets with pT>30 GeV, other jets with pT>20 GeV.
  if(         selected_jets . size() < 2 ){ return false ;} 
  nEvent_passElElCh[ i_step++ ]++ ; 
  if( DLsofterselected_jets . size() < 3 ){ return false ;} 
  nEvent_passElElCh[ i_step++ ]++ ; 

  // MET>=40 GeV for ee and mumu (study removal or adjustment of cut, perhaps also 2D: MET vs. mll).
  if( (*met_pt) < 40 ) { return false ;} 
  nEvent_passElElCh[ i_step++ ]++ ; 

  // >=1 b-tags.
  if( DLsofterselected_bjets . size() < 1 ){ return false ;} 
  nEvent_passElElCh[ i_step++ ]++ ; 

  return true;

}

bool ttHYggdrasilEventSelection::PassMuMu(){

  int i_step = 0 ;
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  // Trigger 
  if( ! passMuMuTrig() ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  // GoodVtx
  if( ! (*goodvtx) ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  // 2 OS-leptons (in DL definition) && FlavMatching with triggers.
  if( DLselected_tightLeptons.size() < 1 ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 
  if(   selected_looseLeptons.size() < 2 ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 
  if(   selected_looseLeptonsIsMuon.at(0) != 1 ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 
  if(   selected_looseLeptonsIsMuon.at(1) != 1 ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  // Veto additional leptons
  if(   selected_looseLeptons.size() != 2 ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  const double invMass = ( (*selected_looseLeptons.at(0)) + (*selected_looseLeptons.at(1)) ).M() ;

  // mll>20 GeV.
  if( invMass < 20 ) { return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  // Exclude Z window (mll >76 GeV and <106 GeV) for ee and mumu.
  if( 76 < invMass &&  invMass < 106 ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  // >=2 jets with pT>30 GeV, other jets with pT>20 GeV.
  if(         selected_jets . size() < 2 ){ return false ;} 
  nEvent_passMuMuCh[ i_step++ ]++ ; 
  if( DLsofterselected_jets . size() < 3 ){ return false ;} 
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  // MET>=40 GeV for ee and mumu (study removal or adjustment of cut, perhaps also 2D: MET vs. mll).
  if( (*met_pt) < 40 ){ return false ;}
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  // >=1 b-tags.
  if( DLsofterselected_bjets . size() < 1 ){ return false ;} 
  nEvent_passMuMuCh[ i_step++ ]++ ; 

  return true;
}

bool ttHYggdrasilEventSelection::PassElMu(){

  int i_step = 0 ;
  nEvent_passElMuCh[ i_step++ ]++ ; 

  // Trigger 
  if( ! passElMuTrig() ){ return false ;}
  nEvent_passElMuCh[ i_step++ ]++ ; 

  // GoodVtx
  if( ! (*goodvtx) ){ return false ;}
  nEvent_passElMuCh[ i_step++ ]++ ; 


  // 2 OS-leptons (in DL definition) && FlavMatching with triggers.
  if( DLselected_tightLeptons.size() < 1 ){ return false ;}
  nEvent_passElMuCh[ i_step++ ]++ ; 

  if(   selected_looseLeptons.size() < 2 ){ return false ;}
  nEvent_passElMuCh[ i_step++ ]++ ; 

  if(  !
       ( ( selected_looseLeptonsIsMuon.at(0) == 1 && selected_looseLeptonsIsMuon.at(1) == 0 )
	 ||
	 ( selected_looseLeptonsIsMuon.at(0) == 0 && selected_looseLeptonsIsMuon.at(1) == 1 )
	 )
       ){ return false ;}
  nEvent_passElMuCh[ i_step++ ]++ ; 

  // Veto additional leptons
  if(   selected_looseLeptons.size() != 2 ){ return false ;}
  nEvent_passElMuCh[ i_step++ ]++ ; 

  const double invMass = ( (*selected_looseLeptons.at(0)) + (*selected_looseLeptons.at(1)) ).M() ;

  // mll>20 GeV.
  if( invMass < 20 ){ return false ;}
  nEvent_passElMuCh[ i_step++ ]++ ; 

  // Exclude Z window (mll >76 GeV and <106 GeV) for ee and mumu.
  // --> Not required in El_Mu ch

  // >=2 jets with pT>30 GeV, other jets with pT>20 GeV.
  if(         selected_jets . size() < 2 ){ return false ;} 
  nEvent_passElMuCh[ i_step++ ]++ ; 
  if( DLsofterselected_jets . size() < 3 ){ return false ;} 
  nEvent_passElMuCh[ i_step++ ]++ ; 

  // MET>=40 GeV for ee and mumu (study removal or adjustment of cut, perhaps also 2D: MET vs. mll).
  // --> Not required in El_Mu ch

  // >=1 b-tags.
  if( DLsofterselected_bjets . size() < 1 ){ return false ;} 
  nEvent_passElMuCh[ i_step++ ]++ ; 

  return true;
}


double ttHYggdrasilEventSelection::_calcDR2( double eta1, double eta2, double phi1, double phi2 ){
  
  double d_eta = eta1 - eta2 ;
  double d_phi = fabs( phi1 - phi2 ) ; 
  d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ; 
  
  return  d_eta*d_eta + d_phi*d_phi ;
  
}

std::vector<const TLorentzVector*> ttHYggdrasilEventSelection::leptons()      {return selected_tightLeptons;}
std::vector<double>                ttHYggdrasilEventSelection::leptonsRelIso(){return selected_tightLeptonsRelIso; }
std::vector<int>                   ttHYggdrasilEventSelection::leptonsIsMuon(){return selected_tightLeptonsIsMuon; }
std::vector<int>                   ttHYggdrasilEventSelection::leptonsCharge(){return selected_tightLeptonsCharge; }
std::vector<double>                ttHYggdrasilEventSelection::leptonsSCEta (){return selected_tightLeptonsScEta; }

std::vector<const TLorentzVector*> ttHYggdrasilEventSelection::looseLeptons()      {return selected_looseLeptons;}
std::vector<double>                ttHYggdrasilEventSelection::looseLeptonsRelIso(){return selected_looseLeptonsRelIso;}
std::vector<int>                   ttHYggdrasilEventSelection::looseLeptonsIsMuon(){return selected_looseLeptonsIsMuon;}
std::vector<int>                   ttHYggdrasilEventSelection::looseLeptonsCharge(){return selected_looseLeptonsCharge;}
std::vector<double>                ttHYggdrasilEventSelection::looseLeptonsScEta (){return selected_looseLeptonsScEta;}

std::vector<const TLorentzVector*> ttHYggdrasilEventSelection::jets()             { return selected_jets;}
std::vector<double> ttHYggdrasilEventSelection::jetsBdiscriminant(){ return selected_jetsBdiscriminant;}
std::vector<int> ttHYggdrasilEventSelection::jetsFlav(){ return selected_jetsFlav;}

std::vector<const TLorentzVector*> ttHYggdrasilEventSelection::bjets()             { return selected_bjets;}
std::vector<double> ttHYggdrasilEventSelection::bjetsBdiscriminant(){ return selected_bjetsBdiscriminant;}
std::vector<int> ttHYggdrasilEventSelection::bjetsFlav(){ return selected_bjetsFlav;}


// ** for DL channel study ***
std::vector<const TLorentzVector*> ttHYggdrasilEventSelection::DLTightLeptons()      {return DLselected_tightLeptons;}
std::vector<double>                ttHYggdrasilEventSelection::DLTightLeptonsRelIso(){return DLselected_tightLeptonsRelIso; }
std::vector<int>                   ttHYggdrasilEventSelection::DLTightLeptonsIsMuon(){return DLselected_tightLeptonsIsMuon; }
std::vector<int>                   ttHYggdrasilEventSelection::DLTightLeptonsCharge(){return DLselected_tightLeptonsCharge; }
std::vector<const TLorentzVector*> ttHYggdrasilEventSelection::DLSofterjets() { return DLsofterselected_jets;}
std::vector<double> ttHYggdrasilEventSelection::DLSofterjetsBdiscriminant()   { return DLsofterselected_jetsBdiscriminant;}
std::vector<const TLorentzVector*> ttHYggdrasilEventSelection::DLSofterbjets(){ return DLsofterselected_bjets;}
std::vector<double> ttHYggdrasilEventSelection::DLSofterbjetsBdiscriminant()  { return DLsofterselected_bjetsBdiscriminant;}


void ttHYggdrasilEventSelection::EnableInfoDumpForDebug(){ b_InfoDumpForDebug = true;}


bool ttHYggdrasilEventSelection::_passTrig( const std::vector< const int * > & triggerSets ){

  for( unsigned int i = 0 ; i < triggerSets . size() ; i++) {
    if( *(triggerSets[i]) == 1 ) return true;
  }
  return false ; 
}

bool ttHYggdrasilEventSelection::passElTrig(){
  return _passTrig( ElTrig );
}

bool ttHYggdrasilEventSelection::passMuTrig(){
  return _passTrig( MuTrig );
}

bool ttHYggdrasilEventSelection::passElElTrig(){
  return _passTrig( ElElTrig );
}

bool ttHYggdrasilEventSelection::passMuMuTrig(){
  return _passTrig( MuMuTrig );
}

bool ttHYggdrasilEventSelection::passElMuTrig(){
  return _passTrig( ElMuTrig );
}

float ttHYggdrasilEventSelection::metAbs(){
  return * met_pt ;
}
float ttHYggdrasilEventSelection::metPhi(){
  return * met_phi ;
}

double ttHYggdrasilEventSelection::getLeptonDR( unsigned int idx ){

  if( lep_trigDR -> size() == 0 ){ return -1; } // trigDR flag has not passed to this tool by user.

  if( idx >= selected_tightLeptons.size() ){ return -1 ; }

  const int OriginalIdx = _GetOriginalIdxOfTightLepton( idx );
  if( OriginalIdx < 0 ) return 100;// failed to find  original idx (should not).

  if( (unsigned int)OriginalIdx >= lep_trigDR -> size() ){
    return -1 ;
  }
  return  lep_trigDR->at( OriginalIdx );
};

int ttHYggdrasilEventSelection::_GetOriginalIdxOfTightLepton( unsigned int idx ){

  for( unsigned int iLep = 0 ; iLep < lep_pt->size() ; iLep ++ ){
    if( fabs(  selected_tightLeptons[idx]->Pt() - lep_pt->at(iLep)) > 0.1 ) continue ;

    if( _calcDR2( selected_tightLeptons[idx]->Eta(), lep_eta ->at(iLep),
		  selected_tightLeptons[idx]->Phi(), lep_phi ->at(iLep) ) > 0.0001 ) continue ;
    return iLep ;
  }
  std::cout <<"[ERROR] _GetOriginalIdxOfTightLepton can not find original lepton" << std::endl ;
  return -1;
}

int ttHYggdrasilEventSelection::getnNonIsoEl(){ return  nNonIsolatedElectron  ;}
int ttHYggdrasilEventSelection::getnNonIsoMu(){ return  nNonIsolatedMuon  ; }


void ttHYggdrasilEventSelection::SetFakeEstimationModeON(  ){ 
  b_FakeEstimation = true ; 
  std::cout <<"ttHYggdrasilEventSelection [debug] : b_FakeEstimation is set to true." << std::endl ; 
}


void ttHYggdrasilEventSelection::SetTightMuonRelIsoCut( double cut ){
  Thre_TightMu_Iso = cut ; 
  std::cout <<"[ttHYggdrasilEventSelection.cc] : Tight muon isolation is set to " << Thre_TightMu_Iso << std::endl ; 
}
void ttHYggdrasilEventSelection::SetLooseMuonRelIsoCut( double cut ){
  Thre_LooseMu_Iso = cut ;
  std::cout <<"[ttHYggdrasilEventSelection.cc] : Loose muon isolation is set to " << Thre_LooseMu_Iso << std::endl ; 
}
