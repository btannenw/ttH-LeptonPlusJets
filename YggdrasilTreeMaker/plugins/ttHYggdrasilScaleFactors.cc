
#ifdef STANDALONECOMPILE
#include "./ttHYggdrasilScaleFactors.h"
#else
#include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilScaleFactors.h"
#endif

#include <iostream>
#include <stdlib.h>
#include <cassert>

#include <TFile.h>


ttHYggdrasilScaleFactors::ttHYggdrasilScaleFactors() :
  initialized( false )
  , period ( period_all )
{

#ifdef STANDALONECOMPILE
  SFfileDir =
    std::string("ttHYggdrasilScaleFactors_data/" );
#else
  SFfileDir =
    (std::string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/" );
#endif

  PileupHistogram . assign( "PileupHistogram_2017data.root" );

}

ttHYggdrasilScaleFactors::ttHYggdrasilScaleFactors( char * sf_file_directory )
  : ttHYggdrasilScaleFactors()
{

  SFfileDir .assign( sf_file_directory );

}


void ttHYggdrasilScaleFactors::init_all(){

  initialized = true ; 
  
  MC_PU_DISTRIBUTION_CHANNEL = MC_PU_DEFAULT ;

  init_btagSF();
  init_Pileup();
  init_ElectronSF();
  init_MuonSF();
  init_TrigMuSF();
  init_TrigElSF();
}

void ttHYggdrasilScaleFactors::init_TrigElSF(){

  {
    std::string input = SFfileDir +"/" + "trig/TriggerSF_Run2016All_v1.root";
    h_EleSF_Trig_SF    = (TH2D*) getTH2HistogramFromFile( input , std::string ("Ele27_WPTight_Gsf") );
  }

}
void ttHYggdrasilScaleFactors::init_TrigMuSF(){

  std::string input = SFfileDir +"/" + "trig/EfficienciesAndSF_RunBtoF_Nov17Nov2017.root";
  TH2D  h_ (* (TH2D*) getTH2HistogramFromFile( input , std::string ("IsoMu27_PtEtaBins/pt_abseta_ratio") ) ) ;
  h_MuSF_Trig_SF  = new TH2D( h_ );
  
}

void ttHYggdrasilScaleFactors::init_ElectronSF(){
  
  {
    std::string input ;

    if(      period == period_B  ){  input  = SFfileDir +"/el/RecoSF2017/egammaEffi.txt_EGM2D_runB_passingRECO.root" ;}
    else if( period == period_C  ){  input  = SFfileDir +"/el/RecoSF2017/egammaEffi.txt_EGM2D_runC_passingRECO.root" ;} 
    else if( period == period_D  ){  input  = SFfileDir +"/el/RecoSF2017/egammaEffi.txt_EGM2D_runD_passingRECO.root" ;} 
    else if( period == period_E  ){  input  = SFfileDir +"/el/RecoSF2017/egammaEffi.txt_EGM2D_runE_passingRECO.root" ;} 
    else if( period == period_F  ){  input  = SFfileDir +"/el/RecoSF2017/egammaEffi.txt_EGM2D_runF_passingRECO.root" ;} 
    else{                           input  = SFfileDir +"/el/RecoSF2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root" ; }    
    
    std::cout << __FILE__ << " Electron tracking SF file : " << input << std::endl ; 
    h_EleSF_ID = (TH2F*) getTH2HistogramFromFile( input , std::string ("EGamma_SF2D") );
    // x : Super cluster eta : -2.5 to 2.5
    // y : ET.
  }
  { 

    std::string input ;
    if(      period == period_B  ){  input  = SFfileDir +"/el/2017CutBaseTight_SF/egammaEffi.txt_EGM2D_runB_passingTight94X.root" ;}
    else if( period == period_C  ){  input  = SFfileDir +"/el/2017CutBaseTight_SF/egammaEffi.txt_EGM2D_runC_passingTight94X.root" ;} 
    else if( period == period_D  ){  input  = SFfileDir +"/el/2017CutBaseTight_SF/egammaEffi.txt_EGM2D_runD_passingTight94X.root" ;} 
    else if( period == period_E  ){  input  = SFfileDir +"/el/2017CutBaseTight_SF/egammaEffi.txt_EGM2D_runE_passingTight94X.root" ;} 
    else if( period == period_F  ){  input  = SFfileDir +"/el/2017CutBaseTight_SF/egammaEffi.txt_EGM2D_runF_passingTight94X.root" ;} 
    else{                            input  = SFfileDir +"/el/2017CutBaseTight_SF/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root" ; }    

    std::cout << __FILE__ << " Electron Reco SF file : " << input << std::endl ; 
    h_EleSF_Reco = (TH2F*) getTH2HistogramFromFile( input , std::string ("EGamma_SF2D") );
    // x : super cluster -2.5 to 2.5
    // y : PT
  }

}

void ttHYggdrasilScaleFactors::init_MuonSF(){


  // Muon ID : setup histogram for 2 files (Moriond17)
  h_MuSF_ID      .clear();
  h_MuSF_ID_Lumi .clear();
  {
    std::string input = SFfileDir +"/" + "muon/ID/EfficienciesAndSF_BCDEF.root"; 
    h_MuSF_ID . push_back( (TH2D*) getTH2HistogramFromFile( input , std::string ("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")) );
    h_MuSF_ID_Lumi . push_back( 19255482132.199 ); // amount of data in the period
  }
  {
    std::string input = SFfileDir +"/" + "muon/ID/EfficienciesAndSF_GH.root";
    h_MuSF_ID . push_back( (TH2D*) getTH2HistogramFromFile( input , std::string ("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")) );
    h_MuSF_ID_Lumi . push_back( 16290016931.807 ); // amount of data in the period
  }
  h_MuSF_ID_LumiTotal = 0 ;
  for( unsigned int i = 0 ; i < h_MuSF_ID_Lumi.size() ; i++ ){
    h_MuSF_ID_LumiTotal += h_MuSF_ID_Lumi[i];
  }

  // Muon track 
  {
    // abs eta 
    std::string input = SFfileDir +"/" + "muon/track/bf/fits.root";
    TGraphAsymmErrors * e_bf = getTGraphFromFile( input , std::string( "ratio_eff_aeta_dr030e030_corr" ) );
    h_muTrack . push_back( ConvertIlldefinedTGraphToTH2D( e_bf ) ) ;
    h_muTrack_down . push_back( ConvertIlldefinedTGraphToTH2D( e_bf , -1 ) ) ;
  }
  {
    std::string input_gh = SFfileDir +"/" + "muon/track/gh/fits.root";
    TGraphAsymmErrors * e_gh = getTGraphFromFile( input_gh , std::string( "ratio_eff_aeta_dr030e030_corr" ) );
    h_muTrack . push_back( ConvertIlldefinedTGraphToTH2D( e_gh ) ) ;
    h_muTrack_down . push_back( ConvertIlldefinedTGraphToTH2D( e_gh , -1 ) ) ;
  }

  
  // Muon Iso : setup histogram with 2 giles (Moridon17)
  { 
    std::string input = SFfileDir +"/" + "muon/ISO/EfficienciesAndSF_BCDEF.root";
    h_MuSF_Iso .push_back(  (TH2D*) getTH2HistogramFromFile( input , std::string ("TightISO_TightID_pt_eta/abseta_pt_ratio") ) );
    h_MuSF_Iso_Lumi . push_back( 19255482132.199 );
  }
  {
    std::string input = SFfileDir +"/" + "muon/ISO/EfficienciesAndSF_GH.root";
    h_MuSF_Iso .push_back( (TH2D*) getTH2HistogramFromFile( input , std::string ("TightISO_TightID_pt_eta/abseta_pt_ratio") ) ) ;
    h_MuSF_Iso_Lumi . push_back( 16290016931.807 );
  }
  h_MuSF_Iso_LumiTotal = 0 ;
  for( unsigned int i = 0 ; i < h_MuSF_Iso_Lumi.size() ; i++ ){
    h_MuSF_Iso_LumiTotal += h_MuSF_Iso_Lumi[i];
  }

}



TH2D * ttHYggdrasilScaleFactors::ConvertIlldefinedTGraphToTH2D( TGraphAsymmErrors * g , int syst ){

  const long nBins = g->GetN();

  double x[ nBins + 1  ];


  x[ 0 ] = g->GetX()[0] - ( g -> GetErrorXlow(0) ) ;
  for( int i = 1 ; i <= nBins ;i++ ){
    x[ i ] = g->GetX()[i-1] + ( g -> GetErrorXhigh(i-1) ) ;
  }
  
  double y[2]={0,100000000};
  TH2D * h = new TH2D( "","",nBins , x , 1, y );
  
  for(  int i = 1 ; i <= nBins ;i++ ){
    h->SetBinContent( i , 1 ,   g->GetY()[i-1] );

    if( syst < 0 ){
      h->SetBinError( i , 1 ,   g->GetEYlow()[i-1] );
    }else{
      h->SetBinError( i , 1 ,   g->GetEYhigh()[i-1] );
    }

  }

  return h ; 

}


TGraphAsymmErrors* ttHYggdrasilScaleFactors::getTGraphFromFile( std::string input , std::string histoname ){

  TDirectory * main_directory = gDirectory;
  TFile * f = TFile::Open( input.c_str() );
  main_directory->cd();

  TGraphAsymmErrors * h_2d_tmp = 0 ;
  f-> GetObject ( histoname.c_str(), h_2d_tmp );
  if( h_2d_tmp != 0 ){
    TGraphAsymmErrors * h_2d = new TGraphAsymmErrors(*h_2d_tmp) ;
    f->Close();
    return h_2d ; 
  }
  std::cout <<"Failed to obtain histogarm named " << histoname<< " from file " << input << std::endl ; 
  assert( false );
  return 0 ; 

}


TH2* ttHYggdrasilScaleFactors::getTH2HistogramFromFile( std::string input , std::string histoname ){

  TDirectory * main_directory = gDirectory;
  TFile * f = TFile::Open( input.c_str() );
  main_directory->cd();
  
  TH2D * h_2d_tmp = 0 ;
  f-> GetObject ( histoname.c_str(), h_2d_tmp );
  if( h_2d_tmp != 0 ){
    TH2D * h_2d = new TH2D(*h_2d_tmp) ;
    f->Close();
    return h_2d ; 
  }

  TH2F * h_2f_tmp = 0 ;
  f-> GetObject ( histoname.c_str(), h_2f_tmp );
  if( h_2f_tmp != 0 ){
    TH2F * h_2f = new TH2F(*h_2f_tmp) ;
    f->Close();
    return h_2f ; 
  }
  std::cout <<"Failed to obtain histogarm named " << histoname<< " from file " << input << std::endl ; 
  assert( false );
  return 0 ; 
}


double ttHYggdrasilScaleFactors::GetBinValueFromXYValues( TH2 * h , double xVal , double yVal , int syst , bool useOveflowBinForX , bool useOveflowBinForY ){

  int bin_x = h->GetXaxis()->FindBin( xVal );
  if( ! useOveflowBinForX && bin_x < 0 ){ bin_x = 1 ;}
  if( ! useOveflowBinForX && bin_x > h->GetXaxis()->GetNbins() ){ bin_x = h->GetXaxis()->GetNbins() ;}

  int bin_y = h->GetYaxis()->FindBin( yVal );
  if(! useOveflowBinForY && bin_y < 0 ){ bin_x = 1 ;}
  if(! useOveflowBinForY && bin_y > h->GetYaxis()->GetNbins() ){ bin_y = h->GetYaxis()->GetNbins() ;}

  const double center_val = h->GetBinContent( bin_x , bin_y );
  if( syst == 0 ){
    return center_val ;
  }

  const double systematic = h->GetBinError( bin_x , bin_y );
  if( syst > 0 ){
    return  center_val +  systematic ;
  }
  return  center_val -  systematic ;
}

double ttHYggdrasilScaleFactors::getTightMuonSF( ttHYggdrasilEventSelection * event ){

  assert( initialized );
  return getTightMuon_IDSF(event ) * getTightMuon_IsoSF(event );

}


double ttHYggdrasilScaleFactors::getTightMuon_IDSF( ttHYggdrasilEventSelection * event , int syst ){

  assert( initialized );
  
  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 

    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
    const double pt      =            event->leptons().at( i )->Pt()  ; 

    double wgt_fot_this_mu = 0 ;
    double wgt_for_this_mu_tracking = 0 ; 
    for( unsigned int iSF = 0 ; iSF < h_MuSF_ID.size() ; iSF ++ ){
      wgt_fot_this_mu +=
	GetBinValueFromXYValues( h_MuSF_ID[iSF] , abs_eta , pt , syst )
	*
	( h_MuSF_ID_Lumi[iSF] / h_MuSF_ID_LumiTotal ) ; //<- Weight based on the int_lumi in the period.

      
      wgt_for_this_mu_tracking +=
	GetBinValueFromXYValues( ( syst >= 0 ? h_muTrack[iSF] : h_muTrack_down[iSF] ) 
				 , abs_eta , 10.0 , syst )// the second value is a dummpy.
	*
	( h_MuSF_ID_Lumi[iSF] / h_MuSF_ID_LumiTotal ) ; //<- Weight based on the int_lumi in the period.
    }

    weight *= wgt_fot_this_mu * wgt_for_this_mu_tracking ;
    
  }
  return weight ;

}

double ttHYggdrasilScaleFactors::getTightMuon_IsoSF( ttHYggdrasilEventSelection * event , int syst ){

  assert( initialized );
  
  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 

    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
    const double pt      =            event->leptons().at( i )->Pt()  ; 

    double wgt_fot_this_mu = 0 ;
    for( unsigned int iSF = 0 ; iSF < h_MuSF_Iso.size() ; iSF ++ ){
      wgt_fot_this_mu +=
	GetBinValueFromXYValues( h_MuSF_Iso[iSF] , abs_eta , pt , syst )
	*
	( h_MuSF_Iso_Lumi[iSF] / h_MuSF_Iso_LumiTotal ) ; //<- Weight based on the int_lumi in the period.
    }
    weight *= wgt_fot_this_mu ;
    
  }
  return weight ;

}


double ttHYggdrasilScaleFactors::getTightElectron_IDSF( ttHYggdrasilEventSelection * event , int syst ){

  assert( initialized );
  
  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 0 ) continue ; 
    
    const double sc_eta =  event->leptonsSCEta().at(i); 
    const double pt     =  event->leptons().at( i )->Pt() ; 
    
    weight *= GetBinValueFromXYValues( h_EleSF_ID , sc_eta , pt , syst );
  }
  return weight ;

}

double ttHYggdrasilScaleFactors::getTightElectron_RecoSF( ttHYggdrasilEventSelection * event , int syst ){

  assert( initialized );

  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 0 ) continue ; 
    
    const double sc_eta =  event->leptonsSCEta().at(i); 
    const double pt     =  event->leptons().at( i )->Pt() ; 
    
    weight *= GetBinValueFromXYValues( h_EleSF_Reco  , sc_eta , pt , syst );
    
  }
  return weight ;
}


double ttHYggdrasilScaleFactors::getTightElectronSF( ttHYggdrasilEventSelection * event ){

  assert( initialized );
    
  return getTightElectron_IDSF( event ) * getTightElectron_RecoSF( event );

}


void ttHYggdrasilScaleFactors::init_btagSF(){


  std::string inputFileHF = SFfileDir +"/" + "btag/csv_rwt_fit_hf_v0_final_2018_2_13.root";
  std::string inputFileLF = SFfileDir +"/" + "btag/csv_rwt_fit_lf_v0_final_2018_2_13.root";

  TFile* fileHF = new TFile ( inputFileHF .c_str());
  TFile* fileLF = new TFile ( inputFileLF .c_str());

  for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting /// only care about the nominal ones
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

    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

  return;


}


ttHYggdrasilScaleFactors::~ttHYggdrasilScaleFactors(){


}


double ttHYggdrasilScaleFactors::get_csv_wgt( std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs, std::vector<int> jetFlavors, 
					      int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){

  assert( initialized );
  
  int iSysHF = 0;
  switch(iSys){
  case 7:  iSysHF=1; break; //JESUp
  case 8:  iSysHF=2; break; //JESDown
  case 9:  iSysHF=3; break; //LFUp
  case 10: iSysHF=4; break; //LFDown
  case 13: iSysHF=5; break; //Stats1Up
  case 14: iSysHF=6; break; //Stats1Down
  case 15: iSysHF=7; break; //Stats2Up
  case 16: iSysHF=8; break; //Stats2Down
  default : iSysHF = 0; break; //NoSys
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
  case 7:  iSysLF=1; break; //JESUp
  case 8:  iSysLF=2; break; //JESDown
  case 11: iSysLF=3; break; //HFUp
  case 12: iSysLF=4; break; //HFDown
  case 17: iSysLF=5; break; //Stats1Up
  case 18: iSysLF=6; break; //Stats1Down
  case 19: iSysLF=7; break; //Stats2Up
  case 20: iSysLF=8; break; //Stats2Down
  default : iSysLF = 0; break; //NoSys
  }

  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;

  for( int iJet=0; iJet<int(jetPts.size()); iJet++ ){

    double csv = jetCSVs[iJet];
    double jetPt = jetPts[iJet];
    double jetAbsEta = fabs(jetEtas[iJet]);
    int flavor = jetFlavors[iJet];

    int iPt = -1; int iEta = -1;
    if (jetPt >=19.99 && jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100) iPt = 4;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;
    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  return csvWgtTotal;


}



double ttHYggdrasilScaleFactors::get_csv_wgt( ttHYggdrasilEventSelection * event,
					      int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF ){

  assert( initialized );
  
  std::vector<double> pt, eta  ; 
  for( unsigned int i = 0 ; i < event->jets().size() ; i++ ){
    pt  .push_back( event -> jets().at(i) -> Pt()  );
    eta .push_back( event -> jets().at(i) -> Eta() );
  }

  return get_csv_wgt(
		     pt , 
		     eta , 
		     event -> jetsBdiscriminant() , 
		     event -> jetsFlav() , 
		     iSys,
		     csvWgtHF,
		     csvWgtLF,
		     csvWgtCF 
		     );

}


void ttHYggdrasilScaleFactors::init_Pileup(){

  // Setting numbers here is just a temporal workaround.
  double PU_MC[ NBINS_PU_REWEIGHTING ] = {
    // RunIIFall17MiniAOD-94X :
    // taken from https://raw.githubusercontent.com/cms-sw/cmssw/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py
    // From 0 to 98 
    3.39597497605e-05,
    6.63688402133e-06,
    1.39533611284e-05,
    3.64963078209e-05,
    6.00872171664e-05,
    9.33932578027e-05,
    0.000120591524486,
    0.000128694546198,
    0.000361697233219,
    0.000361796847553,
    0.000702474896113,
    0.00133766053707,
    0.00237817050805,
    0.00389825605651,
    0.00594546732588,
    0.00856825906255,
    0.0116627396044,
    0.0148793350787,
    0.0179897368379,
    0.0208723871946,
    0.0232564170641,
    0.0249826433945,
    0.0262245860346,
    0.0272704617569,
    0.0283301107549,
    0.0294006137386,
    0.0303026836965,
    0.0309692426278,
    0.0308818046328,
    0.0310566806228,
    0.0309692426278,
    0.0310566806228,
    0.0310566806228,
    0.0310566806228,
    0.0307696426944,
    0.0300103336052,
    0.0288355370103,
    0.0273233309106,
    0.0264343533951,
    0.0255453758796,
    0.0235877272306,
    0.0215627588047,
    0.0195825559393,
    0.0177296309658,
    0.0160560731931,
    0.0146022004183,
    0.0134080690078,
    0.0129586991411,
    0.0125093292745,
    0.0124360740539,
    0.0123547104433,
    0.0123953922486,
    0.0124360740539,
    0.0124360740539,
    0.0123547104433,
    0.0124360740539,
    0.0123387597772,
    0.0122414455005,
    0.011705203844,
    0.0108187105305,
    0.00963985508986,
    0.00827210065136,
    0.00683770076341,
    0.00545237697118,
    0.00420456901556,
    0.00367513566191,
    0.00314570230825,
    0.0022917978982,
    0.00163221454973,
    0.00114065309494,
    0.000784838366118,
    0.000533204105387,
    0.000358474034915,
    0.000238881117601,
    0.0001984254989,
    0.000157969880198,
    0.00010375646169,
    6.77366175538e-05,
    4.39850477645e-05,
    2.84298066026e-05,
    1.83041729561e-05,
    1.17473542058e-05,
    7.51982735129e-06,
    6.16160108867e-06,
    4.80337482605e-06,
    3.06235473369e-06,
    1.94863396999e-06,
    1.23726800704e-06,
    7.83538083774e-07,
    4.94602064224e-07,
    3.10989480331e-07,
    1.94628487765e-07,
    1.57888581037e-07,
    1.2114867431e-07,
    7.49518929908e-08,
    4.6060444984e-08,
    2.81008884326e-08,
    1.70121486128e-08,
    1.02159894812e-08 } ; 


  // Overwrite the PU distribution by the information form the actual PU distribution
  if( MC_PU_DISTRIBUTION_CHANNEL != MC_PU_DEFAULT ){


    std::string path = SFfileDir + "/2017_MC_PU/" + get_MCPUDistributionFileName( MC_PU_DISTRIBUTION_CHANNEL ) +".root" ;
    std::cout << "MC PU distribution overwriting takes the histogram from file : " << path << std::endl ;
    
    TFile * f = TFile::Open( path . c_str() ) ; 
			    
    TH1D * h ; 
    f->GetObject("h_pileup_noWeight" , h ); 

    
    for( int i = 0 ; i < NBINS_PU_REWEIGHTING ; i ++ ){
      PU_MC[i] = h->GetBinContent( i + 3 );
    }
    std::cout <<"MC PU distribution overwriting : done." << std::endl; 

    f->Close();
    
  }
  
  { 
      double total = 0 ;
      for( int i = 0 ; i < NBINS_PU_REWEIGHTING ; i ++ ){
	total += PU_MC[ i ];
      }
      for( int i = 0 ; i < NBINS_PU_REWEIGHTING ; i ++ ){
	PU_MC[ i ] /= total;
      }
    }

    double PU_DATA[ NBINS_PU_REWEIGHTING ] ;
    {
      TH1D * h;
      TFile * f = TFile::Open( (SFfileDir + "/" + PileupHistogram ).c_str() );
      std::cout << "DEBUG : DATA_PU file path = "  << (SFfileDir + "/" + PileupHistogram) << std::endl ; 
      // --> The file has been produced with following command : 
      //  pileupCalc.py -i ./../data/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt  --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 71300 --maxPileupBin 1000 --numPileupBins 1000    MyDataPileupHistogram.root
      
      f -> GetObject( "pileup" , h ) ;
      for( int i = 0 ; i < NBINS_PU_REWEIGHTING ; i++ ){
	PU_DATA[i] = h->GetBinContent( i + 1 );
      }
      f -> Close();
    }



    { 
      double total = 0 ;
      for( int i = 0 ; i < NBINS_PU_REWEIGHTING ; i ++ ){
	total += PU_DATA[ i ];
      }
      // normalization.
      for( int i = 0 ; i < NBINS_PU_REWEIGHTING ; i ++ ){
	PU_DATA[ i ] /= total;
      }
    }

    
    for(int i = 0 ; i < NBINS_PU_REWEIGHTING ; i ++){
      if( PU_MC[i] == 0 ){
	PU_weight[ i ] = 0 ; 
      }else{
	PU_weight[ i ] = PU_DATA[i] / PU_MC[i];
      }
    }
    
}


double ttHYggdrasilScaleFactors::get_pu_wgt( int mc_pu ){

  assert( initialized );

  if( mc_pu >= NBINS_PU_REWEIGHTING ){ mc_pu = NBINS_PU_REWEIGHTING -1 ; }
  if( mc_pu <   0 ){
    mc_pu =   0 ;
  }
  return PU_weight[ mc_pu ];
  
}


double ttHYggdrasilScaleFactors::get_TrigMuEfficiency( ttHYggdrasilEventSelection * event ){

//  double totalInefficiency = 1 ; 
//
//  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
//    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 
//    
//    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
//    const double pt      =  event->leptons().at( i )->Pt() ; 
//
//    double in_efficiency  = 1.0 - GetBinValueFromXYValues( h_MUEff_SingleMuonTrig , abs_eta , pt );
//
//    totalInefficiency *= in_efficiency ;
//  }
//
//  return 1.0 - totalInefficiency ;

  return 1 ; 

}


double ttHYggdrasilScaleFactors::get_TrigElEfficiency( ttHYggdrasilEventSelection * event ){

  assert( initialized );
  return 1 ;
}



double ttHYggdrasilScaleFactors::get_TrigMuSF( ttHYggdrasilEventSelection * event , int syst ){

  assert( initialized );
  
 double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 
    
    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
    const double pt      =  event->leptons().at( i )->Pt() ; 
    
    //    const double trigdr = event->getLeptonDR( i );
    //    const bool isTriggered = trigdr < 0.1 ;

    const double sf = GetBinValueFromXYValues( h_MuSF_Trig_SF  , pt,  abs_eta , syst );

    weight *= sf; // This is not right calculation for the envets with 2 or more electrons. Just for single electron analysis.

  }
  return weight ;
}


double ttHYggdrasilScaleFactors::get_TrigElSF( ttHYggdrasilEventSelection * event , int syst ){

  assert( initialized );
 double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 0 ) continue ; 
    
    const double sc_eta =  event->leptonsSCEta().at(i); 
    const double pt     =  event->leptons().at( i )->Pt() ; 
    
    //    const double trigdr = event->getLeptonDR( i );
    //    const bool isTriggered = trigdr < 0.1 ;

    const bool UseOverflowBinForHighPT = true ; 

    const double sf = GetBinValueFromXYValues( h_EleSF_Trig_SF  , pt , sc_eta , syst, UseOverflowBinForHighPT );

    weight *= sf; // This is not right calculation for the envets with 2 or more electrons. Just for single electron analysis.

  }
  return weight ;

  
}

void ttHYggdrasilScaleFactors::SetupDataPileupFile( std::string filename ){
  PileupHistogram = filename ; 
}



std::string ttHYggdrasilScaleFactors::get_MCPUDistributionFileName( _MC_PU_DIST_CH_NAME ch ){

  
  if( ch == MC_PU_TT_2L   ){ return std::string( "ttto2l2nu"  ); } 
  if( ch == MC_PU_TT_1L   ){ return std::string( "tttosemilep"); } 
  if( ch == MC_PU_Z       ){ return std::string( "zjetsincl"  ); } 
  if( ch == MC_PU_ZLOWMASS){ return std::string( "ZjetLowMass"); } 
  if( ch == MC_PU_W       ){ return std::string( "wjetsincl"  ); } 
  if( ch == MC_PU_WW      ){ return std::string( "ww"         ); } 
  if( ch == MC_PU_WZ      ){ return std::string( "wz"         ); } 
  if( ch == MC_PU_ZZ      ){ return std::string( "zz"         ); } 


  std::cout <<"[FATAL] ttHYggdrasilScaleFactors::get_MCPUDistributionFileName : MC name for the given channel name is not defined." << std::endl ;
  assert ( false );

  return std::string("CAN_NOT_REACH_HERE");

}

void ttHYggdrasilScaleFactors::SetMCPileupChannel( std::string name ){

  if(name == "ttto2l2nu"   ){ MC_PU_DISTRIBUTION_CHANNEL =   MC_PU_TT_2L   ; }
  if(name == "tttosemilep" ){ MC_PU_DISTRIBUTION_CHANNEL =   MC_PU_TT_1L   ; }
  if(name == "zjetsincl"   ){ MC_PU_DISTRIBUTION_CHANNEL =   MC_PU_Z       ; }
  if(name == "ZjetLowMass" ){ MC_PU_DISTRIBUTION_CHANNEL =   MC_PU_ZLOWMASS; }
  if(name == "wjetsincl"   ){ MC_PU_DISTRIBUTION_CHANNEL =   MC_PU_W       ; }
  if(name == "ww"          ){ MC_PU_DISTRIBUTION_CHANNEL =   MC_PU_WW      ; }
  if(name == "wz"          ){ MC_PU_DISTRIBUTION_CHANNEL =   MC_PU_WZ      ; }
  if(name == "zz"          ){ MC_PU_DISTRIBUTION_CHANNEL =   MC_PU_ZZ      ; }

}


void ttHYggdrasilScaleFactors::setAnalysisPeriod( analysis_period p ){

  assert( ! initialized );

  period = p ; 


}


double ttHYggdrasilScaleFactors::get_ElectronZVertexSF(){

  assert( initialized ) ;


  // Values https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations . r12
  if( period == period_B ){ return 0.934 ; }
  // error    0.005

  if( period == period_C ){return 0.992 ; }
  // error   0.001

  if( period == period_D
      ||
      period == period_E
      ||
      period == period_F ) { return 1 ; } // no error is assigned

  // for whole period, 
  return 0.991; 
  // error : 0.001

}
