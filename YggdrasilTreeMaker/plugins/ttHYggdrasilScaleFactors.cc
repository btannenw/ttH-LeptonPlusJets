
#ifdef STANDALONECOMPILE
#include "./ttHYggdrasilScaleFactors.h"
#else
#include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilScaleFactors.h"
#endif

#include <iostream>
#include <stdlib.h>
#include <cassert>

#include <TFile.h>

ttHYggdrasilScaleFactors::ttHYggdrasilScaleFactors(){

#ifdef STANDALONECOMPILE
  SFfileDir =
    std::string("ttHYggdrasilScaleFactors_data/" );
#else
  SFfileDir =
    (std::string(getenv("CMSSW_BASE")) + "/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/" );
#endif
  PileupHistogram . assign( "PileupHistogram_Moriond17.root" );

  init_all();

}

ttHYggdrasilScaleFactors::ttHYggdrasilScaleFactors( char * sf_file_directory ){

  SFfileDir .assign( sf_file_directory );
  PileupHistogram . assign( "PileupHistogram_Moriond17.root" );

  init_all();

}


void ttHYggdrasilScaleFactors::init_all(){

  init_btagSF();
  init_Pileup();
  init_ElectronSF();
  init_MuonSF();
  init_TrigMuSF();
  init_TrigElSF();
}

void ttHYggdrasilScaleFactors::init_TrigElSF(){

  {
    std::string input = SFfileDir +"/" + "trig/ElTriggerPerformance_Jan30.root";
    h_EleSF_Trig_SF    = (TH2F*) getTH2HistogramFromFile( input , std::string ("electrontrig_sf_eta_pt") );
    h_EleSF_TrigEff_MC = (TH2F*) getTH2HistogramFromFile( input , std::string ("electrontrig_eff_mc_eta_pt") );
  }

}
void ttHYggdrasilScaleFactors::init_TrigMuSF(){

  {
    std::string input = SFfileDir +"/" + "trig/MuonTriggerPerformance_Jan30.root";
    h_MuSF_Trig_SF    = (TH2D*) getTH2HistogramFromFile( input , std::string ("muontrig_sf_abseta_pt") );
    h_MuSF_TrigEff_MC = (TH2D*) getTH2HistogramFromFile( input , std::string ("muontrig_eff_mc_abseta_pt") );
  }

  // Root file taken from https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults?rev=15
  {
    
    std::string input = SFfileDir +"/" + "SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root";
    h_MUEff_SingleMuonTrig = (TH2D*) getTH2HistogramFromFile( input , std::string ("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097/efficienciesDATA/abseta_pt_DATA") );
  }

}

void ttHYggdrasilScaleFactors::init_ElectronSF(){
  
  {
    std::string input = SFfileDir +"/" + "el/gsf_trackingEff/egammaEffi.txt_SF2D.root" ;
    h_EleSF_ID = (TH2F*) getTH2HistogramFromFile( input , std::string ("EGamma_SF2D") );
  }
  { 
    std::string input = SFfileDir +"/" + "el/CutBase76xMedium__SFfor80x/egammaEffi.txt_SF2D.root";
    h_EleSF_Reco = (TH2F*) getTH2HistogramFromFile( input , std::string ("EGamma_SF2D") );
  }

}

void ttHYggdrasilScaleFactors::init_MuonSF(){


  // Muon ID : setup histogram for 2 files (Moriond17)
  h_MuSF_ID      .clear();
  h_MuSF_ID_Lumi .clear();
  {
    std::string input = SFfileDir +"/" + "muon/ID/EfficienciesAndSF_BCDEF.root"; 
    h_MuSF_ID . push_back( (TH2D*) getTH2HistogramFromFile( input , std::string ("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")) );
    h_MuSF_ID_Lumi . push_back( 20236179680.471 ); // amount of data in the period
  }
  {
    std::string input = SFfileDir +"/" + "muon/ID/EfficienciesAndSF_GH.root";
    h_MuSF_ID . push_back( (TH2D*) getTH2HistogramFromFile( input , std::string ("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio")) );
    h_MuSF_ID_Lumi . push_back( 16578401309.707 ); // amount of data in the period
  }
  h_MuSF_ID_LumiTotal = 0 ;
  for( unsigned int i = 0 ; i < h_MuSF_ID_Lumi.size() ; i++ ){
    h_MuSF_ID_LumiTotal += h_MuSF_ID_Lumi[i];
  }

  
  // Muon Iso : setup histogram with 2 giles (Moridon17)
  { 
    std::string input = SFfileDir +"/" + "muon/ISO/EfficienciesAndSF_BCDEF.root";
    h_MuSF_Iso .push_back(  (TH2D*) getTH2HistogramFromFile( input , std::string ("TightISO_TightID_pt_eta/abseta_pt_ratio") ) );
    h_MuSF_Iso_Lumi . push_back( 20236179680.471 );
  }
  {
    std::string input = SFfileDir +"/" + "muon/ISO/EfficienciesAndSF_GH.root";
    h_MuSF_Iso .push_back( (TH2D*) getTH2HistogramFromFile( input , std::string ("TightISO_TightID_pt_eta/abseta_pt_ratio") ) ) ;
    h_MuSF_Iso_Lumi . push_back( 16578401309.707 );
  }
  h_MuSF_Iso_LumiTotal = 0 ;
  for( unsigned int i = 0 ; i < h_MuSF_Iso_Lumi.size() ; i++ ){
    h_MuSF_Iso_LumiTotal += h_MuSF_Iso_Lumi[i];
  }


}


TH2* ttHYggdrasilScaleFactors::getTH2HistogramFromFile( std::string input , std::string histoname ){

  TFile * f = TFile::Open( input.c_str() );
  TH2F * h_2f = 0 ;
  TH2D * h_2d = 0 ;
  f-> GetObject ( histoname.c_str(), h_2d );
  f-> GetObject ( histoname.c_str(), h_2f );
  if( h_2d != 0 ){
    return h_2d ; 
  }
  if( h_2f != 0 ){
    return h_2f ; 
  }
  std::cout <<"Failed to obtain histogarm named " << histoname<< " from file " << input << std::endl ; 
  assert( false );
  return 0 ; 
}


double ttHYggdrasilScaleFactors::GetBinValueFromXYValues( TH2 * h , double xVal , double yVal ){

  int bin_x = h->GetXaxis()->FindBin( xVal );
  if( bin_x < 0 ){ bin_x = 1 ;}
  if( bin_x > h->GetXaxis()->GetNbins() ){ bin_x = h->GetXaxis()->GetNbins() ;}

  int bin_y = h->GetYaxis()->FindBin( yVal );
  if( bin_y < 0 ){ bin_x = 1 ;}
  if( bin_y > h->GetYaxis()->GetNbins() ){ bin_y = h->GetYaxis()->GetNbins() ;}

  return h->GetBinContent( bin_x , bin_y );

}

double ttHYggdrasilScaleFactors::getTightMuonSF( ttHYggdrasilEventSelection * event ){

  return getTightMuon_IDSF(event ) * getTightMuon_IsoSF(event );

}


double ttHYggdrasilScaleFactors::getTightMuon_IDSF( ttHYggdrasilEventSelection * event ){

  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 

    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
    const double pt      =            event->leptons().at( i )->Pt()  ; 

    double wgt_fot_this_mu = 0 ;
    for( unsigned int iSF = 0 ; iSF < h_MuSF_ID.size() ; iSF ++ ){
      wgt_fot_this_mu +=
	GetBinValueFromXYValues( h_MuSF_ID[iSF] , abs_eta , pt )
	*
	( h_MuSF_ID_Lumi[iSF] / h_MuSF_ID_LumiTotal ) ; //<- Weight based on the int_lumi in the period.
    }
    weight *= wgt_fot_this_mu ;
    
  }
  return weight ;

}

double ttHYggdrasilScaleFactors::getTightMuon_IsoSF( ttHYggdrasilEventSelection * event ){

  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 

    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
    const double pt      =            event->leptons().at( i )->Pt()  ; 

    double wgt_fot_this_mu = 0 ;
    for( unsigned int iSF = 0 ; iSF < h_MuSF_Iso.size() ; iSF ++ ){
      wgt_fot_this_mu +=
	GetBinValueFromXYValues( h_MuSF_Iso[iSF] , abs_eta , pt )
	*
	( h_MuSF_Iso_Lumi[iSF] / h_MuSF_Iso_LumiTotal ) ; //<- Weight based on the int_lumi in the period.
    }
    weight *= wgt_fot_this_mu ;
    
  }
  return weight ;

}


double ttHYggdrasilScaleFactors::getTightElectron_IDSF( ttHYggdrasilEventSelection * event ){

  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 0 ) continue ; 
    
    const double sc_eta =  event->leptonsSCEta().at(i); 
    const double pt     =  event->leptons().at( i )->Pt() ; 
    
    weight *= GetBinValueFromXYValues( h_EleSF_ID , sc_eta , pt );
  }
  return weight ;

}

double ttHYggdrasilScaleFactors::getTightElectron_RecoSF( ttHYggdrasilEventSelection * event ){

  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 0 ) continue ; 
    
    const double sc_eta =  event->leptonsSCEta().at(i); 
    const double pt     =  event->leptons().at( i )->Pt() ; 
    
    weight *= GetBinValueFromXYValues( h_EleSF_Reco  , sc_eta , pt );
    
  }
  return weight ;
}


double ttHYggdrasilScaleFactors::getTightElectronSF( ttHYggdrasilEventSelection * event ){

  return getTightElectron_IDSF( event ) * getTightElectron_RecoSF( event );

}


void ttHYggdrasilScaleFactors::init_btagSF(){


  std::string inputFileHF = SFfileDir +"/" + "csv_rwt_fit_hf_v2_final_2017_1_10test.root";
  std::string inputFileLF = SFfileDir +"/" + "csv_rwt_fit_lf_v2_final_2017_1_10test.root";

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
  double PU_MC[75] = {
    // values from https://github.com/cms-sw/cmssw/blob/7c04e20a023bb47d86f417f6e13bc41dee6890ee/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py#L25
    1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05
  };


    { 
      double total = 0 ;
      for( int i = 0 ; i < 75 ; i ++ ){
	total += PU_MC[ i ];
      }
      for( int i = 0 ; i < 75 ; i ++ ){
	PU_MC[ i ] /= total;
      }
    }

    double PU_DATA[75] ;
    {
      TH1D * h;
      TFile * f = TFile::Open( (SFfileDir + "/" + PileupHistogram ).c_str() );
      std::cout << "DEBUG : DATA_PU file path = "  << (SFfileDir + "/" + PileupHistogram) << std::endl ; 
      // --> The file has been produced with following command : 
      //  pileupCalc.py -i ./../data/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt  --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 71300 --maxPileupBin 1000 --numPileupBins 1000    MyDataPileupHistogram.root
      
      f -> GetObject( "pileup" , h ) ;
      for( int i = 0 ; i < 75 ; i++ ){
	PU_DATA[i] = h->GetBinContent( i + 1 );
      } 
    }



    { 
      double total = 0 ;
      for( int i = 0 ; i < 75 ; i ++ ){
	total += PU_DATA[ i ];
      }
      // normalization.
      for( int i = 0 ; i < 75 ; i ++ ){
	PU_DATA[ i ] /= total;
      }
    }

    
    for(int i = 0 ; i < 75 ; i ++){
      if( PU_MC[i] == 0 ){
	PU_weight[ i ] = 0 ; 
      }else{
	PU_weight[ i ] = PU_DATA[i] / PU_MC[i];
      }
    }
    
}


double ttHYggdrasilScaleFactors::get_pu_wgt( int mc_pu ){

  return PU_weight[ mc_pu ];
  
}


double ttHYggdrasilScaleFactors::get_TrigMuEfficiency( ttHYggdrasilEventSelection * event ){

  double totalInefficiency = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 
    
    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
    const double pt      =  event->leptons().at( i )->Pt() ; 

    double in_efficiency  = 1.0 - GetBinValueFromXYValues( h_MUEff_SingleMuonTrig , abs_eta , pt );

    totalInefficiency *= in_efficiency ;
  }

  return 1.0 - totalInefficiency ;

}


double ttHYggdrasilScaleFactors::get_TrigElEfficiency( ttHYggdrasilEventSelection * event ){

  return 1 ;
}



double ttHYggdrasilScaleFactors::get_TrigMuSF( ttHYggdrasilEventSelection * event ){

 double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 
    
    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
    const double pt      =  event->leptons().at( i )->Pt() ; 
    
    const double trigdr = event->getLeptonDR( i );
    const bool isTriggered = trigdr < 0.1 ;

    const double sf = GetBinValueFromXYValues( h_MuSF_Trig_SF  , abs_eta , pt );
    if( isTriggered ){
      weight *= sf;
    }else{
      const double mc_eff = GetBinValueFromXYValues( h_MuSF_TrigEff_MC , abs_eta , pt );
      weight *= ( 1.0 - sf * mc_eff )/( 1.0 - mc_eff );
    }
  }
  return weight ;
}


double ttHYggdrasilScaleFactors::get_TrigElSF( ttHYggdrasilEventSelection * event ){

 double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 0 ) continue ; 
    
    const double sc_eta =  event->leptonsSCEta().at(i); 
    const double pt     =  event->leptons().at( i )->Pt() ; 
    
    const double trigdr = event->getLeptonDR( i );
    const bool isTriggered = trigdr < 0.1 ;

    const double sf = GetBinValueFromXYValues( h_EleSF_Trig_SF  , sc_eta , pt );

    if( isTriggered ){
      weight *= sf;
    }else{
      const double mc_eff = GetBinValueFromXYValues( h_EleSF_TrigEff_MC , sc_eta , pt );
      weight *= ( 1.0 - sf * mc_eff )/( 1.0 - mc_eff );
    }
    
  }
  return weight ;

  
}

void ttHYggdrasilScaleFactors::SetupDataPileupFile( std::string filename ){
  PileupHistogram = filename ; 
}
