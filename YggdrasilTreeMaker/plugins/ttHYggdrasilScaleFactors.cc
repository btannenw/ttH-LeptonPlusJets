
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
  PileupHistogram . assign( "Data_E_and_F.root" );

  init_all();

}

ttHYggdrasilScaleFactors::ttHYggdrasilScaleFactors( char * sf_file_directory ){

  SFfileDir .assign( sf_file_directory );
  PileupHistogram . assign( "Data_E_and_F.root" );

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
    std::string input = SFfileDir +"/" + "eleTrig_SF.root";
    h_EleSF_Trig = (TH2F*) getTH2HistogramFromFile( input , std::string ("h_eleTrig_SF") );
  }

}
void ttHYggdrasilScaleFactors::init_TrigMuSF(){

  {
    std::string input = SFfileDir +"/" + "SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root";
    h_MuSF_Trig_HLTv4p2 = (TH2D*) getTH2HistogramFromFile( input , std::string ("runD_IsoMu20_OR_IsoTkMu20_HLTv4p2_PtEtaBins/abseta_pt_ratio") );
    h_MuSF_Trig_HLTv4p3 = (TH2D*) getTH2HistogramFromFile( input , std::string ("runD_IsoMu20_OR_IsoTkMu20_HLTv4p3_PtEtaBins/abseta_pt_ratio") );
  }

  // Root file taken from https://twiki.cern.ch/twiki/bin/view/CMS/MuonWorkInProgressAndPagResults?rev=15
  {
    
    std::string input = SFfileDir +"/" + "SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root";
    h_MUEff_SingleMuonTrig = (TH2D*) getTH2HistogramFromFile( input , std::string ("IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097/efficienciesDATA/abseta_pt_DATA") );
  }

}

void ttHYggdrasilScaleFactors::init_ElectronSF(){
  
  {
    std::string input = SFfileDir +"/" + "eleRECO.txt.egamma_SF2D.root";
    h_EleSF_ID = (TH2F*) getTH2HistogramFromFile( input , std::string ("EGamma_SF2D") );
  }
  { 
    std::string input = SFfileDir +"/" + "ScaleFactor_GsfElectronToRECO_passingTrigWP80.txt.egamma_SF2D.root";
    h_EleSF_Reco = (TH2F*) getTH2HistogramFromFile( input , std::string ("EGamma_SF2D") );
  }

}

void ttHYggdrasilScaleFactors::init_MuonSF(){

  {
    std::string input = SFfileDir +"/" + "MuonID_Z_2016runB_2p6fb.root";
    h_MuSF_ID = (TH2D*) getTH2HistogramFromFile( input , std::string ("MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio"));
  }
  { 
    std::string input = SFfileDir +"/" + "MuonISO_Z_2016runB_2p6fb.root";
    h_MuSF_Iso = (TH2D*) getTH2HistogramFromFile( input , std::string ("MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio") );
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

    weight *= GetBinValueFromXYValues( h_MuSF_ID , abs_eta , pt );
    
  }
  return weight ;

}

double ttHYggdrasilScaleFactors::getTightMuon_IsoSF( ttHYggdrasilEventSelection * event ){

  double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 1 ) continue ; 

    const double abs_eta = std::fabs( event->leptons().at( i )->Eta() ) ; 
    const double pt      =            event->leptons().at( i )->Pt()  ; 

    weight *= GetBinValueFromXYValues( h_MuSF_Iso  , abs_eta , pt );
    
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


  std::string inputFileHF = SFfileDir +"/" + "csv_rwt_fit_hf_v2_final_2016_06_30test.root";
  std::string inputFileLF = SFfileDir +"/" + "csv_rwt_fit_lf_v2_final_2016_06_30test.root";

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
  double PU_MC[50] = {  0.000829312873542,
		      0.00124276120498,
		      0.00339329181587,
		      0.00408224735376,
		      0.00383036590008,
		      0.00659159288946,
		      0.00816022734493,
		      0.00943640833116,
		      0.0137777376066,
		      0.017059392038,
		      0.0213193035468,
		      0.0247343174676,
		      0.0280848773878,
		      0.0323308476564,
		      0.0370394341409,
		      0.0456917721191,
		      0.0558762890594,
		      0.0576956187107,
		      0.0625325287017,
		      0.0591603758776,
		      0.0656650815128,
		      0.0678329011676,
		      0.0625142146389,
		      0.0548068448797,
		      0.0503893295063,
		      0.040209818868,
		      0.0374446988111,
		      0.0299661572042,
		      0.0272024759921,
		      0.0219328403791,
		      0.0179586571619,
		      0.0142926728247,
		      0.00839941654725,
		      0.00522366397213,
		      0.00224457976761,
		      0.000779274977993,
		      0.000197066585944,
		      7.16031761328e-05,
		      0.0,
		      0.0,
		      0.0,
		      0.0,
		      0.0,
		      0.0,
		      0.0,
		      0.0,
		      0.0,
		      0.0,
		      0.0,
		      0.0 };


    { 
      double total = 0 ;
      for( int i = 0 ; i < 50 ; i ++ ){
	total += PU_MC[ i ];
      }
      for( int i = 0 ; i < 50 ; i ++ ){
	PU_MC[ i ] /= total;
      }
    }

    double PU_DATA[50] ;
    {
      TH1D * h;
      TFile * f = TFile::Open( (SFfileDir + "/" + PileupHistogram ).c_str() );
      std::cout << "DEBUG : DATA_PU file path = "  << (SFfileDir + "/" + PileupHistogram) << std::endl ; 
      // --> The file has been produced with following command : 
      //  pileupCalc.py -i ./../data/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt  --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 71300 --maxPileupBin 1000 --numPileupBins 1000    MyDataPileupHistogram.root
      
      f -> GetObject( "pileup" , h ) ;
      for( int i = 0 ; i < 50 ; i++ ){
	PU_DATA[i] = h->GetBinContent( i + 1 );
      } 
    }



    { 
      double total = 0 ;
      for( int i = 0 ; i < 50 ; i ++ ){
	total += PU_DATA[ i ];
      }
      // normalization.
      for( int i = 0 ; i < 50 ; i ++ ){
	PU_DATA[ i ] /= total;
      }
    }

    
    for(int i = 0 ; i < 50 ; i ++){
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
    
    double w_p2 = GetBinValueFromXYValues( h_MuSF_Trig_HLTv4p2 , abs_eta , pt );
    double w_p3 = GetBinValueFromXYValues( h_MuSF_Trig_HLTv4p3 , abs_eta , pt );

    double ratio_p2 =  754.394 / ( 754.394 + 1908.010 ) ;
    double ratio_p3 = 1908.010 / ( 754.394 + 1908.010 ) ;
    
    weight *= 
      w_p2 * ratio_p2 
      +
      w_p3 * ratio_p3 ;

  }
  return weight ;
}


double ttHYggdrasilScaleFactors::get_TrigElSF( ttHYggdrasilEventSelection * event ){

 double weight = 1 ; 

  for( unsigned int i = 0 ; i < event->leptonsIsMuon().size() ; i++ ){
    if( event->leptonsIsMuon().at( i ) != 0 ) continue ; 
    
    const double sc_eta =  event->leptonsSCEta().at(i); 
    const double pt     =  event->leptons().at( i )->Pt() ; 
    
    weight *= GetBinValueFromXYValues( h_EleSF_Trig , pt ,  sc_eta ); // <- Unlike Reco/Iso SF, x=PT and y=SC_Eta.
    
  }
  return weight ;

  
}

void ttHYggdrasilScaleFactors::SetupDataPileupFile( std::string filename ){
  PileupHistogram = filename ; 
}
