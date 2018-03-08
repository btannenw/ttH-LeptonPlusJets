
#ifndef TTHYGGDRASILSCALEFACTORS_H_
#define TTHYGGDRASILSCALEFACTORS_H_

#include <vector>
#include <TLorentzVector.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>

#include <TGraphAsymmErrors.h>

#define NBINS_PU_REWEIGHTING 99


#ifdef STANDALONECOMPILE
#include "ttHYggdrasilEventSelection.h"
#else
#include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilEventSelection.h"
#endif


class ttHYggdrasilScaleFactors{

 public :

  ttHYggdrasilScaleFactors();
  ttHYggdrasilScaleFactors( char * sf_file_directory );
  ~ttHYggdrasilScaleFactors();

  
  double get_csv_wgt( ttHYggdrasilEventSelection * event,
		      int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );

  double get_csv_wgt( std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs, std::vector<int> jetFlavors, 
		      int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF );


  double get_pu_wgt( int mc_pu );

  double getTightMuonSF( ttHYggdrasilEventSelection * event );
  double getTightElectronSF( ttHYggdrasilEventSelection * event );

  double getTightMuon_IDSF( ttHYggdrasilEventSelection * event , int syst = 0 );
  double getTightElectron_IDSF( ttHYggdrasilEventSelection * event , int syst = 0 );

  double getTightMuon_IsoSF( ttHYggdrasilEventSelection * event , int syst = 0 );
  double getTightElectron_RecoSF( ttHYggdrasilEventSelection * event , int syst = 0 );

  // - - -

  double get_TrigMuSF( ttHYggdrasilEventSelection * event , int syst = 0 );
  double get_TrigElSF( ttHYggdrasilEventSelection * event , int syst = 0 );

  double get_TrigElEfficiency( ttHYggdrasilEventSelection * event );
  double get_TrigMuEfficiency( ttHYggdrasilEventSelection * event );

  double get_ElectronZVertexSF();
  
  // Can replace PU file if you want.
  void SetupDataPileupFile( std::string filename ); // just filename. the file should be in the common data directory.

  void SetMCPileupChannel( std::string name );

  void init_all();

  enum analysis_period {
    period_all ,
    period_B ,
    period_C ,
    period_D ,
    period_E ,
    period_F  
  };
  void setAnalysisPeriod( analysis_period p ); // If you use this, use before init_all();
  
 private :

  bool initialized ;

  analysis_period period ;
  
  std::string SFfileDir ;
  std::string PileupHistogram;

  void init_Pileup();
  void init_btagSF();
  void init_ElectronSF();
  void init_MuonSF();
  void init_TrigMuSF();
  void init_TrigElSF();
  TH2 * getTH2HistogramFromFile( std::string input , std::string histoname );
  double GetBinValueFromXYValues( TH2 * h , double xVal , double yVal , int syst = 0
				  , bool useOveflowBinForX = false , bool useOveflowBinForY = false ) ;

  TGraphAsymmErrors* getTGraphFromFile( std::string input , std::string histoname );

  TH2D * ConvertIlldefinedTGraphToTH2D( TGraphAsymmErrors * g , int syst = 0 );

  // CSV reweighting
  TH1D* h_csv_wgt_hf[9][5];
  TH1D* c_csv_wgt_hf[9][5];
  TH1D* h_csv_wgt_lf[9][4][3];

  // PU weighting
  double PU_weight[ NBINS_PU_REWEIGHTING ];


  // Lepton SF
  TH2F * h_EleSF_ID ;
  TH2F * h_EleSF_Reco;
  std::vector< TH2D *> h_MuSF_ID;
  std::vector< double> h_MuSF_ID_Lumi;
  double h_MuSF_ID_LumiTotal;
  
  std::vector< TH2D *> h_MuSF_Iso ;
  std::vector< double > h_MuSF_Iso_Lumi ;
  double h_MuSF_Iso_LumiTotal ;

  // Trif SF
  TH2D * h_MuSF_Trig_SF;
  //  TH2D * h_MuSF_TrigEff_MC;
  TH2D * h_EleSF_Trig_SF;
  //  TH2F * h_EleSF_TrigEff_MC;

  // Trig Efficiency
  // TH2D * h_MUEff_SingleMuonTrig;

  std::vector < TH2D * > h_muTrack ; 
  std::vector < TH2D * > h_muTrack_down ; 

  enum _MC_PU_DIST_CH_NAME{
    MC_PU_DEFAULT  = 0
    , MC_PU_TT_2L    = 1
    , MC_PU_TT_1L    = 2
    , MC_PU_Z        = 3
    , MC_PU_ZLOWMASS = 4
    , MC_PU_W        = 5
    , MC_PU_WW       = 6
    , MC_PU_WZ       = 7
    , MC_PU_ZZ       = 8
    , MC_PU_SINGLETOP_SCH 
    , MC_PU_SINGLETOP_TW     
    , MC_PU_SINGLETOP_TBARW  
    , MC_PU_SINGLETOP_TCH_TBAR
    , MC_PU_SINGLETOP_TCA_T   
    , MC_PU_TTH
    , MC_PU_MuEnrichQCDPt1000toInf 
    , MC_PU_MuEnrichQCDPt120to170  
    , MC_PU_MuEnrichQCDPt15to20    
    , MC_PU_MuEnrichQCDPt170to300  
    , MC_PU_MuEnrichQCDPt20to30    
    , MC_PU_MuEnrichQCDPt300to470  
    , MC_PU_MuEnrichQCDPt30to50    
    , MC_PU_MuEnrichQCDPt470to600  
    , MC_PU_MuEnrichQCDPt50to80    
    , MC_PU_MuEnrichQCDPt600to800  
    , MC_PU_MuEnrichQCDPt800to1000 
    , MC_PU_MuEnrichQCDPt80to120   
  };

  _MC_PU_DIST_CH_NAME   MC_PU_DISTRIBUTION_CHANNEL;

  std::string get_MCPUDistributionFileName( _MC_PU_DIST_CH_NAME ch );

};

#endif
