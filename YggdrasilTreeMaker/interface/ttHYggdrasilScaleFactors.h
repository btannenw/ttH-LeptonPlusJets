
#ifndef TTHYGGDRASILSCALEFACTORS_H_
#define TTHYGGDRASILSCALEFACTORS_H_

#include <vector>
#include <TLorentzVector.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>

#include <TGraphAsymmErrors.h>

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

  double getTightMuon_IDSF( ttHYggdrasilEventSelection * event );
  double getTightElectron_IDSF( ttHYggdrasilEventSelection * event );

  double getTightMuon_IsoSF( ttHYggdrasilEventSelection * event );
  double getTightElectron_RecoSF( ttHYggdrasilEventSelection * event );

  // - - -

  double get_TrigMuSF( ttHYggdrasilEventSelection * event );
  double get_TrigElSF( ttHYggdrasilEventSelection * event );

  double get_TrigElEfficiency( ttHYggdrasilEventSelection * event );
  double get_TrigMuEfficiency( ttHYggdrasilEventSelection * event );


  // Can replace PU file if you want.
  void SetupDataPileupFile( std::string filename ); // just filename. the file should be in the common data directory.

  void init_Pileup();

 private :

  std::string SFfileDir ;
  std::string PileupHistogram;

  void init_all();
  void init_btagSF();
  void init_ElectronSF();
  void init_MuonSF();
  void init_TrigMuSF();
  void init_TrigElSF();
  TH2 * getTH2HistogramFromFile( std::string input , std::string histoname );
  double GetBinValueFromXYValues( TH2 * h , double xVal , double yVal 
				  , bool useOveflowBinForX = false , bool useOveflowBinForY = false ) ;

  TGraphAsymmErrors* getTGraphFromFile( std::string input , std::string histoname );

  TH2D * ConvertIlldefinedTGraphToTH2D( TGraphAsymmErrors * g );

  // CSV reweighting
  TH1D* h_csv_wgt_hf[9][5];
  TH1D* c_csv_wgt_hf[9][5];
  TH1D* h_csv_wgt_lf[9][4][3];

  // PU weighting
  double PU_weight[75];


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

};

#endif
