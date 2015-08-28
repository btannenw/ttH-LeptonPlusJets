// -*- C++ -*-
//
// Package:    ttH-LeptonPlusJets/TTHSyncExercise
// Class:      TTHSyncExercise
// 
/**\class TTHSyncExercise TTHSyncExercise.cc ttH-LeptonPlusJets/AnalysisCode/plugins/TTHSyncExercise.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Darren Puigh
//         Created:  Sun, 15 Feb 2015 20:01:00 GMT
//
//


// system include files
#include <memory>

#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/strbitset.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1HFRings.h"
#include "DataFormats/L1Trigger/interface/L1HFRingsFwd.h"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidateFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "BoostedTTH/BoostedAnalyzer/interface/BoostedUtils.hpp"
#include "MiniAOD/BoostedObjects/interface/SubFilterJet.h"
#include "MiniAOD/BoostedObjects/interface/HTTTopJet.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"
#include "MiniAOD/MiniAODHelper/interface/TopTagger.h"
#include "MiniAOD/MiniAODHelper/interface/HiggsTagger.h"
//
// class declaration
//

class TTHSyncExercise : public edm::EDAnalyzer {
   public:
      explicit TTHSyncExercise(const edm::ParameterSet&);
      ~TTHSyncExercise();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const& iRun,edm::EventSetup const& iSetup) override;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  //--------tt+X categorization
  const edm::EDGetTokenT<int> genTtbarIdToken_;

  //----------

  edm::EDGetTokenT <edm::TriggerResults> triggerResultsToken;
  edm::EDGetTokenT <edm::TriggerResults> filterResultsToken;

  edm::EDGetTokenT <reco::VertexCollection> vertexToken;
  edm::EDGetTokenT <reco::VertexCompositePtrCandidateCollection> secondaryVertexToken;
  edm::EDGetTokenT <pat::ElectronCollection> electronToken;
  edm::EDGetTokenT <pat::MuonCollection> muonToken;
  edm::EDGetTokenT <pat::JetCollection> jetToken;
  edm::EDGetTokenT <pat::METCollection> metToken;

  edm::EDGetTokenT< boosted::HTTTopJetCollection > topJetsToken;
  edm::EDGetTokenT< boosted::SubFilterJetCollection > subFilterJetsToken;

  edm::EDGetTokenT <pat::PackedCandidateCollection> packedpfToken;

  edm::EDGetTokenT <reco::BeamSpot> beamspotToken;
  edm::EDGetTokenT <double> rhoToken;
  edm::EDGetTokenT <reco::GenParticleCollection> mcparicleToken;
  edm::EDGetTokenT <std::vector< PileupSummaryInfo > > puInfoToken;

  edm::EDGetTokenT <GenEventInfoProduct> genInfoProductToken;


  HLTConfigProvider hlt_config_;
  HLTConfigProvider filter_config_;

  std::string hltTag;
  std::string filterTag;

  std::map<std::string, int> hlt_cppath_;
  std::map<std::string, int> flt_cppath_;

  bool verbose_;
  bool dumpHLT_;

  int numEvents_;
  int numPos_;
  int numNeg_;

  int numHbb;
  int numHdecayToOneParticle;
  int numHdecayToOneGluon;
  int numHdecayToOnePhoton;
  int numHdecayToOneZBoson;
  int numHdecayToOneWBoson;


  int insample_;

  std::vector<std::string> hlt_triggerNames_;
  std::vector<std::string> flt_filterNames_;

  edm::Service<TFileService> fs_;


  double mySample_xSec_;
  double mySample_nGen_;
  double intLumi_;

  double ptmax;
  int NptBins;

  double jetptmax;
  int NjetptBins;

  int NjetMin;
  int NjetMax;
  int NjetBins;
  int NtagMin;
  int NtagMax;
  int NtagBins;
  int NpuMin;
  int NpuMax;
  int NpuBins;


  // Declare histograms
  TH1D *h_tth_syncex1_ele;
  TH1D *h_tth_syncex1_mu;
  TH1D *h_tth_syncex1_diele;
  TH1D *h_tth_syncex1_dimuon;
  TH1D *h_tth_syncex1_elemuon;

  TH1D *h_hlt;
  TH1D *h_flt;

  TH1D* h_electron_selection;
  TH1D* h_muon_selection;

  MiniAODHelper miniAODhelper;
  TopTagger toptagger;

  // ofstream syncOutputFile;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TTHSyncExercise::TTHSyncExercise(const edm::ParameterSet& iConfig):
  genTtbarIdToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("genTtbarId")))
{
   //now do what ever initialization is needed
  verbose_ = false;
  dumpHLT_ = false;

  hltTag = "HLT";
  filterTag = "PAT";

  triggerResultsToken = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), hltTag));
  filterResultsToken = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), filterTag));

  vertexToken = consumes <reco::VertexCollection> (edm::InputTag(std::string("offlineSlimmedPrimaryVertices")));
  secondaryVertexToken = consumes <reco::VertexCompositePtrCandidateCollection> (edm::InputTag(std::string("slimmedSecondaryVertices")));
  electronToken = consumes <pat::ElectronCollection> (edm::InputTag(std::string("slimmedElectrons")));
  muonToken = consumes <pat::MuonCollection> (edm::InputTag(std::string("slimmedMuons")));
  jetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("slimmedJets")));
  metToken = consumes <pat::METCollection> (edm::InputTag(std::string("slimmedMETs")));

  topJetsToken = consumes< boosted::HTTTopJetCollection >(edm::InputTag("HTTTopJetMatcher","htttopjets","p"));
  subFilterJetsToken = consumes< boosted::SubFilterJetCollection >(edm::InputTag("SFJetMatcher","subfilterjets","p"));

  packedpfToken = consumes <pat::PackedCandidateCollection> (edm::InputTag(std::string("packedPFCandidates")));

  beamspotToken = consumes <reco::BeamSpot> (edm::InputTag(std::string("offlineBeamSpot")));
  rhoToken = consumes <double> (edm::InputTag(std::string("fixedGridRhoFastjetAll")));
  mcparicleToken = consumes <reco::GenParticleCollection> (edm::InputTag(std::string("prunedGenParticles")));
  puInfoToken = consumes <std::vector< PileupSummaryInfo > > (edm::InputTag(std::string("addPileupInfo")));

  genInfoProductToken = consumes <GenEventInfoProduct> (edm::InputTag(std::string("generator")));


  h_electron_selection = fs_->make<TH1D>("h_electron_selection",";electron cut", 12, 0 , 12 );
  h_muon_selection = fs_->make<TH1D>("h_muon_selection",";muon cut", 12, 0 , 12 );

  h_electron_selection->GetXaxis()->SetBinLabel(1,"All");
  h_electron_selection->GetXaxis()->SetBinLabel(2,"p_{T}>20, |#eta|<2.4, !inCrack");
  h_electron_selection->GetXaxis()->SetBinLabel(3,"full5x5 #sigma_{i#eta i#eta}");
  h_electron_selection->GetXaxis()->SetBinLabel(4,"|#Delta #eta_{in}|");
  h_electron_selection->GetXaxis()->SetBinLabel(5,"|#Delta #phi_{in}|");
  h_electron_selection->GetXaxis()->SetBinLabel(6,"hOverE");
  h_electron_selection->GetXaxis()->SetBinLabel(7,"1/E - 1/p");
  h_electron_selection->GetXaxis()->SetBinLabel(8,"d0");
  h_electron_selection->GetXaxis()->SetBinLabel(9,"dZ");
  h_electron_selection->GetXaxis()->SetBinLabel(10,"expectedMissingInnerHits");
  h_electron_selection->GetXaxis()->SetBinLabel(11,"passConversionVeto");
  h_electron_selection->GetXaxis()->SetBinLabel(12,"relIso (#Delta Beta, 0.3)");

  h_muon_selection->GetXaxis()->SetBinLabel(1,"All");
  h_muon_selection->GetXaxis()->SetBinLabel(2,"p_{T}>20, |#eta|<2.4");
  h_muon_selection->GetXaxis()->SetBinLabel(3,"GlobalMuon || TrackerMuon");
  h_muon_selection->GetXaxis()->SetBinLabel(4,"PFMuon");
  h_muon_selection->GetXaxis()->SetBinLabel(5,"#Chi^{2}");
  h_muon_selection->GetXaxis()->SetBinLabel(6,"validMuonHit");
  h_muon_selection->GetXaxis()->SetBinLabel(7,"validPixelHit");
  h_muon_selection->GetXaxis()->SetBinLabel(8,"trk layers w/meas");
  h_muon_selection->GetXaxis()->SetBinLabel(9,"matched stations");
  h_muon_selection->GetXaxis()->SetBinLabel(10,"d0");
  h_muon_selection->GetXaxis()->SetBinLabel(11,"dZ");
  h_muon_selection->GetXaxis()->SetBinLabel(12,"relIso < 0.1");


  h_tth_syncex1_ele = fs_->make<TH1D>("h_tth_syncex1_ele",";cut", 8, 0 , 8 );
  h_tth_syncex1_ele->GetXaxis()->SetBinLabel(1,"All events");
  h_tth_syncex1_ele->GetXaxis()->SetBinLabel(2,"Single ele trig");
  h_tth_syncex1_ele->GetXaxis()->SetBinLabel(3,"==1 electron");
  h_tth_syncex1_ele->GetXaxis()->SetBinLabel(4,"==0 muons");
  h_tth_syncex1_ele->GetXaxis()->SetBinLabel(5,">=4 jets");
  h_tth_syncex1_ele->GetXaxis()->SetBinLabel(6,">=2 b-tags");
  h_tth_syncex1_ele->GetXaxis()->SetBinLabel(7,">=1 top-tags");
  h_tth_syncex1_ele->GetXaxis()->SetBinLabel(8,">=1 Higgs-tags");

  h_tth_syncex1_mu = fs_->make<TH1D>("h_tth_syncex1_mu",";cut", 8, 0 , 8 );
  h_tth_syncex1_mu->GetXaxis()->SetBinLabel(1,"All events");
  h_tth_syncex1_mu->GetXaxis()->SetBinLabel(2,"Single mu trig");
  h_tth_syncex1_mu->GetXaxis()->SetBinLabel(3,"==1 muon");
  h_tth_syncex1_mu->GetXaxis()->SetBinLabel(4,"==0 electrons");
  h_tth_syncex1_mu->GetXaxis()->SetBinLabel(5,">=4 jets");
  h_tth_syncex1_mu->GetXaxis()->SetBinLabel(6,">=2 b-tags");
  h_tth_syncex1_mu->GetXaxis()->SetBinLabel(7,">=1 top-tags");
  h_tth_syncex1_mu->GetXaxis()->SetBinLabel(8,">=1 Higgs-tags");

  h_tth_syncex1_dimuon = fs_->make<TH1D>("h_tth_syncex1_dimuon",";cut", 8, 0 , 8 );
  h_tth_syncex1_dimuon->GetXaxis()->SetBinLabel(1,"All events");
  h_tth_syncex1_dimuon->GetXaxis()->SetBinLabel(2,"dimuon trig");
  h_tth_syncex1_dimuon->GetXaxis()->SetBinLabel(3,">=2 OS muons");
  h_tth_syncex1_dimuon->GetXaxis()->SetBinLabel(4,"mll>20 GeV");
  h_tth_syncex1_dimuon->GetXaxis()->SetBinLabel(5,"Z veto");
  h_tth_syncex1_dimuon->GetXaxis()->SetBinLabel(6,">=2 jets");
  h_tth_syncex1_dimuon->GetXaxis()->SetBinLabel(7,"MET>40GeV");
  h_tth_syncex1_dimuon->GetXaxis()->SetBinLabel(8,">=1 b-tags");

  h_tth_syncex1_diele = fs_->make<TH1D>("h_tth_syncex1_diele",";cut", 8, 0 , 8 );
  h_tth_syncex1_diele->GetXaxis()->SetBinLabel(1,"All events");
  h_tth_syncex1_diele->GetXaxis()->SetBinLabel(2,"diele trig");
  h_tth_syncex1_diele->GetXaxis()->SetBinLabel(3,">=2 OS electrons");
  h_tth_syncex1_diele->GetXaxis()->SetBinLabel(4,"mll>20 GeV");
  h_tth_syncex1_diele->GetXaxis()->SetBinLabel(5,"Z veto");
  h_tth_syncex1_diele->GetXaxis()->SetBinLabel(6,">=2 jets");
  h_tth_syncex1_diele->GetXaxis()->SetBinLabel(7,"MET>40GeV");
  h_tth_syncex1_diele->GetXaxis()->SetBinLabel(8,">=1 b-tags");

  h_tth_syncex1_elemuon = fs_->make<TH1D>("h_tth_syncex1_elemuon",";cut", 6, 0 , 6 );
  h_tth_syncex1_elemuon->GetXaxis()->SetBinLabel(1,"All events");
  h_tth_syncex1_elemuon->GetXaxis()->SetBinLabel(2,"elemuon trig");
  h_tth_syncex1_elemuon->GetXaxis()->SetBinLabel(3,">=2 OS emu");
  h_tth_syncex1_elemuon->GetXaxis()->SetBinLabel(4,"mll>20 GeV");
  // h_tth_syncex1_elemuon->GetXaxis()->SetBinLabel(5,"Z veto");
  h_tth_syncex1_elemuon->GetXaxis()->SetBinLabel(5,">=2 jets");
  // h_tth_syncex1_elemuon->GetXaxis()->SetBinLabel(7,"MET>40GeV");
  h_tth_syncex1_elemuon->GetXaxis()->SetBinLabel(6,">=1 b-tags");



  std::string era = "2012_53x";
  insample_ = 2500;

  mySample_xSec_ = 831.76;
  mySample_nGen_ = 25446993;
  intLumi_ = 10000;

  analysisType::analysisType iAnalysisType = analysisType::LJ;
  bool isData = !true;

  miniAODhelper.SetUp(era, insample_, iAnalysisType, isData);
  toptagger = TopTagger(TopTag::Likelihood, TopTag::CSV, "toplikelihoodtaggerhistos.root");

  // syncOutputFile.open("ttHSyncEx1_LJ.txt");
}


TTHSyncExercise::~TTHSyncExercise()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TTHSyncExercise::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  numEvents_++;

  int run  = iEvent.id().run();
  int lumi = iEvent.luminosityBlock();
  int event = iEvent.id().event();

  //---------
  //---- comment out checks
  // if(! ( (lumi ==  535 && event ==  53417 ) || (lumi ==  664 && event ==  66386 ) || (lumi == 1061 && event == 106001 ) || (lumi ==   46 && event ==   4599 ) || (lumi == 1001 && event == 100060 ) ) ) return;

  // if(! ( (lumi ==  83  && event ==  8206 ) || (lumi ==  320 && event ==  31947 ) ) )  return;

  // printf("looking at event lumi:event %d, %d\n", lumi, event);  
  //---------

  double wgt = 1;

  double xSec = mySample_xSec_;
  double nGen = mySample_nGen_;

  wgt = intLumi_ * xSec / nGen;

  /// weight for NLO samples
  edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
  iEvent.getByToken(genInfoProductToken,GenEventInfoHandle);

  double GenEventInfoWeight = GenEventInfoHandle.product()->weight();

  if (GenEventInfoWeight > 0) numPos_++;
  else numNeg_++;

  wgt *= GenEventInfoWeight;

  double minTightLeptonPt = 30.;


  h_hlt->Fill(0.,1);
  h_flt->Fill(0.,1);

  bool passSingleElectronTrigger = false;
  bool passSingleMuonTrigger = false;
  bool passDiElectronTrigger = false;
  bool passDiMuonTrigger = false;
  bool passEleMuonTrigger = false;

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken, triggerResults);

  if( triggerResults.isValid() ){
    std::vector<std::string> triggerNames = hlt_config_.triggerNames();

    for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){
      std::string pathName = triggerNames[iPath];
      unsigned int hltIndex = hlt_config_.triggerIndex(pathName);

      if( hltIndex >= triggerResults->size() ) continue;

      int accept = triggerResults->accept(hltIndex);
      int prescale = -1;//hlt_config_.prescaleValue(iEvent, iSetup, pathName);

     if( accept ){
	if( pathName=="HLT_Ele27_eta2p1_WP85_Gsf_HT200_v1" )    passSingleElectronTrigger = true;
	if( pathName=="HLT_IsoMu24_eta2p1_v1" ) passSingleMuonTrigger = true;
	if( pathName=="HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1" ) passDiElectronTrigger = true;
	if( pathName=="HLT_Mu30_TkMu11_v1" || pathName=="HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1" || pathName=="HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1" ) passDiMuonTrigger = true;
	if( pathName=="HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1" || pathName=="HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1" ) passEleMuonTrigger = true;
      }

      if( verbose_ && dumpHLT_ ) std::cout << " =====>  HLT: path name = " << pathName << ",\t prescale = " << prescale << ",\t pass = " << accept << std::endl; 

      std::string pathNameNoVer = hlt_config_.removeVersion(pathName);

      if( accept ) hlt_cppath_[pathNameNoVer]+=1;

      if( accept ){
	TAxis * axis = h_hlt->GetXaxis();
	if( !axis ) continue;
	int bin_num = axis->FindBin(pathNameNoVer.c_str());
	int bn = bin_num - 1;
	h_hlt->Fill(bn, 1);
      }
    }
  }
  else{
    std::cout << "Trigger results not valid for tag " << hltTag << std::endl;
  }
  // std::cout << passSingleMuonTrigger << passSingleElectronTrigger << passDiMuonTrigger << passDiElectronTrigger << passEleMuonTrigger << std::endl;
  if (passSingleMuonTrigger == passSingleElectronTrigger || passDiMuonTrigger || passDiElectronTrigger || passEleMuonTrigger) ;

  edm::Handle<edm::TriggerResults> filterResults;
  iEvent.getByToken(filterResultsToken, filterResults);

  if( filterResults.isValid() ){
    std::vector<std::string> triggerNames = filter_config_.triggerNames();

    for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){
      std::string pathName = triggerNames[iPath];
      unsigned int hltIndex = filter_config_.triggerIndex(pathName);

      if( hltIndex >= filterResults->size() ) continue;

      int accept = filterResults->accept(hltIndex);
      int prescale = -1;//filter_config_.prescaleValue(iEvent, iSetup, pathName);

      if( verbose_ && dumpHLT_ ) std::cout << " =====>  Filter: path name = " << pathName << ",\t prescale = " << prescale << ",\t pass = " << accept << std::endl; 

      std::string pathNameNoVer = filter_config_.removeVersion(pathName);

      if( accept ) flt_cppath_[pathNameNoVer]+=1;

      if( accept ){
	TAxis * axis = h_flt->GetXaxis();
	if( !axis ) continue;
	int bin_num = axis->FindBin(pathNameNoVer.c_str());
	int bn = bin_num - 1;
	h_flt->Fill(bn, 1);
      }
    }
  }
  else{
    std::cout << "Trigger results not valid for tag " << filterTag << std::endl;
  }


  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(vertexToken,vtxHandle);
  reco::VertexCollection vtxs = *vtxHandle;

  edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVtxs;
  iEvent.getByToken(secondaryVertexToken,secondaryVtxs);

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken,electrons);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken,muons);

  edm::Handle<pat::JetCollection> pfjets;
  iEvent.getByToken(jetToken,pfjets);

  edm::Handle<pat::METCollection> pfmet;
  iEvent.getByToken(metToken,pfmet);

  edm::Handle<pat::PackedCandidateCollection> packedPFcands;
  iEvent.getByToken(packedpfToken,packedPFcands);

  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(beamspotToken,bsHandle);

  edm::Handle<reco::GenParticleCollection> mcparticles;
  iEvent.getByToken(mcparicleToken,mcparticles);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken,rhoHandle);

  edm::Handle<std::vector< PileupSummaryInfo > > PupInfo;
  iEvent.getByToken(puInfoToken,PupInfo);

  edm::Handle<boosted::HTTTopJetCollection> h_htttopjet;
  iEvent.getByToken( topJetsToken,h_htttopjet);

  edm::Handle<boosted::SubFilterJetCollection> h_subfilterjet;
  iEvent.getByToken( subFilterJetsToken,h_subfilterjet );

  ////----------------------
  ////---- tt+X Categorization
  ////----------------------
  edm::Handle<int> genTtbarId;
  iEvent.getByToken(genTtbarIdToken_, genTtbarId);
  int additionalJetEventId = *genTtbarId%100;


    ////----------- set up rho for lepton effArea Isolation correction
  double rho_event = ( (rhoHandle.isValid()) ) ? *rhoHandle : -99;
  miniAODhelper.SetRho(rho_event);
  // std::cout << "------>> rho value is "<< rho_event << std::endl;
  ///
  math::XYZPoint beamSpotPosition;
  beamSpotPosition.SetCoordinates(0,0,0);
  double BSx=0,BSy=0,BSz=0;
  if( (bsHandle.isValid()) ){
    reco::BeamSpot bs = *bsHandle;
    BSx = bs.x0();
    BSy = bs.y0();
    BSz = bs.z0();
    beamSpotPosition = bsHandle->position();
  }


  if( verbose_ ) printf("\t BeamSpot: x = %.2f,\t y = %.2f,\t z = %.2f \n", BSx, BSy, BSz );

  int numpv=0;
  bool firstGoodPV = false;
  reco::Vertex vertex;
  if( vtxHandle.isValid() ){
    reco::VertexCollection::const_iterator vtx = vtxs.begin();
    bool isGood = ( !(vtx->isFake()) &&
		    (vtx->ndof() >= 4.0) &&
		    (abs(vtx->z()) <= 24.0) &&
		    (abs(vtx->position().Rho()) <= 2.0) 
		    );

    if( isGood ){
      numpv++;
      firstGoodPV = true;
      vertex = (*vtx);
    }
    // for( reco::VertexCollection::const_iterator vtx = vtxs.begin(); vtx!=vtxs.end(); ++vtx ){
    //   bool isGood = ( !(vtx->isFake()) &&
    // 		      (vtx->ndof() >= 4.0) &&
    // 		      (abs(vtx->z()) <= 24.0) &&
    // 		      (abs(vtx->position().Rho()) <= 2.0) 
    // 		      );

    //   if( !isGood ) continue;

    //   if( numpv==0 ) vertex = (*vtx);

    //   numpv++;
    // }
  }

  if( !firstGoodPV ){
    std::cout << "ERROR!! First PV does not pass requirements! Skipping event" << std::endl;
    return;
  }

  if( verbose_ ) printf("\t Event PV: x = %.3f,\t y = %.3f,\t z = %.3f \n", vertex.x(), vertex.y(), vertex.z());

  if( numpv>0 ) miniAODhelper.SetVertex(vertex);

  const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );   //Get the jet corrector from the event setup

  miniAODhelper.SetJetCorrector(corrector);


  /////////
  ///
  /// Electrons
  ///
  ////////
  std::vector<pat::Electron> selectedElectrons = miniAODhelper.GetSelectedElectrons( *electrons, 10.0, electronID::electronPhys14L, 2.4 );
  std::vector<pat::Electron> selectedElectrons_Tight = miniAODhelper.GetSelectedElectrons( *electrons, minTightLeptonPt, electronID::electronPhys14M, 2.1 );

  int numElectrons = int(selectedElectrons.size());
  int numElectrons_Tight = int(selectedElectrons_Tight.size());


  /////////
  ///
  /// Muons
  ///
  ////////
  std::vector<pat::Muon> selectedMuons = miniAODhelper.GetSelectedMuons( *muons, 10.0, muonID::muonLoose );
  std::vector<pat::Muon> selectedMuons_Tight = miniAODhelper.GetSelectedMuons( *muons, minTightLeptonPt, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.1 );

  int numMuons = int(selectedMuons.size());
  int numMuons_Tight = int(selectedMuons_Tight.size());



  int numMu=0;
  for( std::vector<pat::Muon>::const_iterator pfmu = muons->begin(); pfmu!=muons->end(); ++pfmu ){
    int ncut = 0;
    h_muon_selection->Fill(0.5+ncut++, 1);

    if( pfmu->pt()>20 && fabs(pfmu->eta())<2.4 ) h_muon_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( (pfmu->isGlobalMuon()) ) h_muon_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( (pfmu->isPFMuon()) ) h_muon_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pfmu->globalTrack().isAvailable() ){
      if( (pfmu->globalTrack()->normalizedChi2() < 10.) ) h_muon_selection->Fill(0.5+ncut++, 1);
      else continue;

      if( (pfmu->globalTrack()->hitPattern().numberOfValidMuonHits() > 0) ) h_muon_selection->Fill(0.5+ncut++, 1);
      else continue;
    }
    else continue;

    if( pfmu->innerTrack().isAvailable() ){
      if( (pfmu->innerTrack()->hitPattern().numberOfValidPixelHits() > 0) ) h_muon_selection->Fill(0.5+ncut++, 1);
      else continue;
    }
    else continue;

    if( pfmu->track().isAvailable() ){
      if( (pfmu->track()->hitPattern().trackerLayersWithMeasurement() > 5) ) h_muon_selection->Fill(0.5+ncut++, 1);
      else continue;
    }
    else continue;

    if( (pfmu->numberOfMatchedStations() > 1) ) h_muon_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pfmu->muonBestTrack().isAvailable() ){
      if( (fabs(pfmu->muonBestTrack()->dxy(vertex.position())) < 0.2) ) h_muon_selection->Fill(0.5+ncut++, 1);
      else continue;

      if( (fabs(pfmu->muonBestTrack()->dz(vertex.position())) < 0.5) ) h_muon_selection->Fill(0.5+ncut++, 1);
      else continue;
    }
    else continue;

    if( miniAODhelper.GetMuonRelIso(*pfmu,coneSize::R04,corrType::deltaBeta) < 0.12 ) h_muon_selection->Fill(0.5+ncut++, 1);
    else continue;

    numMu++;
  }


  int numEle=0;
  for( std::vector<pat::Electron>::const_iterator pfele = electrons->begin(); pfele!=electrons->end(); ++pfele ){
    int ncut = 0;
    h_electron_selection->Fill(0.5+ncut++, 1);

    bool inCrack = false;
    if( pfele->superCluster().isAvailable() )
      inCrack = ( fabs(pfele->superCluster()->position().eta())>1.4442 && fabs(pfele->superCluster()->position().eta())<1.5660 );

    if( pfele->pt()>20 && fabs(pfele->eta())<2.4 && !inCrack ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    double SCeta = (pfele->superCluster().isAvailable()) ? pfele->superCluster()->position().eta() : -99;
    double absSCeta = fabs(SCeta);

    bool isEB = ( absSCeta < 1.479 );

    double pfIsoCharged = pfele->pfIsolationVariables().sumChargedHadronPt;
    double pfIsoNeutralHadron = pfele->pfIsolationVariables().sumNeutralHadronEt;
    double pfIsoNeutralPhoton = pfele->pfIsolationVariables().sumPhotonEt;
    double pfIsoSumPUPt = pfele->pfIsolationVariables().sumPUPt;

    double relIso = (pfIsoCharged + std::max( pfIsoNeutralHadron + pfIsoNeutralPhoton - 0.5*pfIsoSumPUPt, 0.0 ))/pfele->pt();

    double full5x5_sigmaIetaIeta = pfele->full5x5_sigmaIetaIeta();
    double dEtaIn = fabs( pfele->deltaEtaSuperClusterTrackAtVtx() );
    double dPhiIn = fabs( pfele->deltaPhiSuperClusterTrackAtVtx() );
    double hOverE = pfele->hcalOverEcal();

    double ooEmooP = -999;
    if( pfele->ecalEnergy() == 0 ) ooEmooP = 1e30;
    else if( !std::isfinite(pfele->ecalEnergy()) ) ooEmooP = 1e30;
    else ooEmooP = fabs(1.0/pfele->ecalEnergy() - pfele->eSuperClusterOverP()/pfele->ecalEnergy() );

    double d0 = -999;
    double dZ = -999;
    double expectedMissingInnerHits = -999;
    if( pfele->gsfTrack().isAvailable() ){
      d0 = fabs(pfele->gsfTrack()->dxy(vertex.position()));
      dZ = fabs(pfele->gsfTrack()->dz(vertex.position()));
      expectedMissingInnerHits = pfele->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    }

    bool passConversionVeto = ( pfele->passConversionVeto() );

    bool pass_full5x5_sigmaIetaIeta = ( isEB ) ? (full5x5_sigmaIetaIeta < 0.010399) : (full5x5_sigmaIetaIeta < 0.029524);
    bool pass_dEtaIn = ( isEB ) ? (dEtaIn < 0.007641) : (dEtaIn < 0.009285);
    bool pass_dPhiIn = ( isEB ) ? (dPhiIn < 0.032643) : (dPhiIn < 0.042447);
    bool pass_hOverE = ( isEB ) ? (hOverE < 0.060662) : (hOverE < 0.104263);
    bool pass_ooEmooP = ( isEB ) ? (ooEmooP < 0.153897) : (ooEmooP < 0.137468);
    bool pass_d0 = ( isEB ) ? (d0 < 0.011811) : (d0 < 0.051682);
    bool pass_dZ = ( isEB ) ? (dZ < 0.070775) : (dZ < 0.180720);
    bool pass_expectedMissingInnerHits = ( isEB ) ? (expectedMissingInnerHits <= 1) : (expectedMissingInnerHits <= 1);
    bool pass_passConversionVeto = ( isEB ) ? (passConversionVeto) : (passConversionVeto);
    bool pass_relIso = ( isEB ) ? (relIso < 0.097213) : (relIso < 0.116708);

    if( pass_full5x5_sigmaIetaIeta ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_dEtaIn ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_dPhiIn ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_hOverE ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_ooEmooP ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_d0 ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_dZ ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_expectedMissingInnerHits ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_passConversionVeto ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    if( pass_relIso ) h_electron_selection->Fill(0.5+ncut++, 1);
    else continue;

    numEle++;
  }



  // Jet selection
  std::vector<pat::Jet> Jets_ID = miniAODhelper.GetSelectedJets(*pfjets, 0., 999, jetID::jetLoose, '-' );
  std::vector<pat::Jet> rawJets = miniAODhelper.GetUncorrectedJets(Jets_ID);
  // std::vector<pat::Jet> jetsNoMu = miniAODhelper.RemoveOverlaps(selectedMuons, rawJets_ID);
  // std::vector<pat::Jet> jetsNoEle = miniAODhelper.RemoveOverlaps(selectedElectrons, jetsNoMu);
  std::vector<pat::Jet> correctedJets = miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup);
  // std::vector<pat::Jet> selectedJets_unsorted = miniAODhelper.GetSelectedJets(correctedJets, 30., 2.4, jetID::jetLoose, '-' );
  // std::vector<pat::Jet> selectedJets_tag_unsorted = miniAODhelper.GetSelectedJets( correctedJets, 30., 2.4, jetID::jetLoose, 'M' );

  std::vector<pat::Jet> selectedJets_unsorted = miniAODhelper.GetSelectedJets(correctedJets, 30., 2.4, jetID::none, '-' );
  std::vector<pat::Jet> selectedJets_tag_unsorted = miniAODhelper.GetSelectedJets( correctedJets, 30., 2.4, jetID::none, 'M' );

  // std::vector<pat::Jet> selectedJets_unsorted = miniAODhelper.GetSelectedJets(*pfjets, 25., 2.4, jetID::none, '-' );
  // std::vector<pat::Jet> selectedJets_tag_unsorted = miniAODhelper.GetSelectedJets( *pfjets, 25., 2.4, jetID::none, 'M' );


  // Sort jet collections by pT
  std::vector<pat::Jet> selectedJets_uncleaned     = miniAODhelper.GetSortedByPt( selectedJets_unsorted );
  std::vector<pat::Jet> selectedJets_tag_uncleaned = miniAODhelper.GetSortedByPt( selectedJets_tag_unsorted );

  std::vector<pat::Jet> selectedJets     = miniAODhelper.GetDeltaRCleanedJets( selectedJets_uncleaned, selectedMuons, selectedElectrons, 0.4);
  std::vector<pat::Jet> selectedJets_tag = miniAODhelper.GetDeltaRCleanedJets( selectedJets_tag_uncleaned, selectedMuons, selectedElectrons, 0.4);

  // std::vector<pat::Jet> selectedJets_tmp     = miniAODhelper.GetSortedByPt( selectedJets_unsorted );
  // std::vector<pat::Jet> selectedJets;

  // for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_tmp.begin(); iJet!=selectedJets_tmp.end(); ++iJet ){
  //   // double dR_lep_jet = 99.;
  //   bool isOverlap = false;

  //   TLorentzVector jet_p4;
  //   jet_p4.SetPxPyPzE(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());

  //   for( std::vector<pat::Electron>::const_iterator iEle = selectedElectrons.begin(); iEle != selectedElectrons.end(); iEle++ ){ 
  //     TLorentzVector ele_p4;
  //     ele_p4.SetPxPyPzE(iEle->px(),iEle->py(),iEle->pz(),iEle->energy());
  //     double delta_tmp = jet_p4.DeltaR(ele_p4);
  //     if(delta_tmp < 0.4){
  // 	isOverlap = true;
  // 	break;
  //     } 
  //   }

  //   if( isOverlap ) continue;

  //   for( std::vector<pat::Muon>::const_iterator iMuon = selectedMuons.begin(); iMuon != selectedMuons.end(); iMuon++ ){ 
  //     TLorentzVector muon_p4;
  //     muon_p4.SetPxPyPzE(iMuon->px(),iMuon->py(),iMuon->pz(),iMuon->energy());
  //     double delta_tmp = jet_p4.DeltaR(muon_p4);
  //     if(delta_tmp < 0.4){
  // 	isOverlap = true;
  // 	break;
  //     } 
  //   }

  //   if( isOverlap ) continue;

  //   selectedJets.push_back(*iJet);

  // }

  // std::vector<pat::Jet> selectedJets_tag = miniAODhelper.GetSelectedJets( selectedJets, 25., 2.4, jetID::none, 'M' );;
  ///// Boosted jet information

  ///// HEP top tagged jet
  int numTopTags = 0;
  vector<boosted::HTTTopJet> syncTopJets;
  if( h_htttopjet.isValid() ){
    boosted::HTTTopJetCollection const &htttopjets_unsorted = *h_htttopjet;
    boosted::HTTTopJetCollection htttopjets = BoostedUtils::GetSortedByPt(htttopjets_unsorted);

    for( boosted::HTTTopJetCollection::iterator topJet = htttopjets.begin() ; topJet != htttopjets.end(); ++topJet ){
      // pt and eta requirements on top jet
      if( !(topJet->fatjet.pt() > 200. && abs(topJet->fatjet.eta()) < 2) ) continue;

      // pt and eta requirements on subjets
      if( !( (topJet->nonW.pt()>20 && abs(topJet->nonW.eta())<2.4 ) &&
	     (topJet->W1.pt()>20 && abs(topJet->W1.eta())<2.4 ) &&
	     (topJet->W2.pt()>20 && abs(topJet->W2.eta())<2.4 ) ) ) continue;

      //////
      if(toptagger.GetTopTaggerOutput(*topJet)>-1){
	numTopTags++;
	syncTopJets.push_back(*topJet);
      }

      // bool toptag = BoostedUtils::GetTopTag(*topJet);
      // // must be top-tagged
      // if( !toptag ) continue;
      // numTopTags++;

    }

  }

  ///// Higgs tagged jet
  int numHiggsTags = 0;
  if( h_subfilterjet.isValid() ){
    boosted::SubFilterJetCollection const &subfilterjets_unsorted = *h_subfilterjet;
    boosted::SubFilterJetCollection subfilterjets = BoostedUtils::GetSortedByPt(subfilterjets_unsorted);

    for( boosted::SubFilterJetCollection::iterator higgsJet = subfilterjets.begin() ; higgsJet != subfilterjets.end(); ++higgsJet ){
      // pt and eta requirements on top jet
      if( !(higgsJet->fatjet.pt() > 200. && abs(higgsJet->fatjet.eta()) < 2) ) continue;

      //remove overlap with topjets
      bool overlapping=false;
      for(auto tj=syncTopJets.begin(); tj!=syncTopJets.end(); tj++){
	if(BoostedUtils::DeltaR(tj->fatjet,higgsJet->fatjet)<1.5){
	  overlapping=true;
	  break;
	}
      }
      if(overlapping) continue;

      int numBtagFiltJets=0;
      std::vector<pat::Jet> filterjets = higgsJet->filterjets;
      int numFiltJets = filterjets.size();
      for( int ijet=0; ijet<numFiltJets; ijet++ ){
	// if( verbose_ ){
	//   printf("\t\t filt jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = %.2f,\t CSVv2 = %+5.3f,\t CSVv1 = %+5.3f \n",
	// 	 ijet, filterjets[ijet].pt(), filterjets[ijet].eta(), filterjets[ijet].phi(), 
	// 	 filterjets[ijet].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"),
	// 	 filterjets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"));
	// }

	if( !(filterjets[ijet].pt()>20. && abs(filterjets[ijet].eta()) < 2.4) ) continue;
	if( !(filterjets[ijet].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89) ) continue;
	numBtagFiltJets++;
      }

      // if( verbose_ ){
      // 	printf("\t Higgs jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = %.2f,\t numFiltJets = %2d,\t numBtagFiltJets = %2d\n",
      // 	       int(higgsJet - subfilterjets.begin()), higgsJet->fatjet.pt(), higgsJet->fatjet.eta(), higgsJet->fatjet.phi(), numFiltJets, numBtagFiltJets );
      // }

      if( numBtagFiltJets>=2 ) numHiggsTags++;
    }
  }


  // Get numJets, numTags
  int numJets = int( selectedJets.size() );
  int numTags = int( selectedJets_tag.size() );

  // Get Corrected MET
  pat::MET correctedMET = pfmet->front();//miniAODhelper.GetCorrectedMET( pfmets.at(0), pfJets_forMET, iSysType );
  double MET_pt = correctedMET.pt();
  double MET_phi = correctedMET.phi();

  //---- comment out checks
  // for(int ijet=0; ijet<numJets; ijet++ ){
  //   double pt = selectedJets.at(ijet).pt();
  //   double eta = selectedJets.at(ijet).eta();
  //   double phi = selectedJets.at(ijet).phi();

  //   std::cout << "---->jet "<< ijet << " pt:eta:phi is" << pt << ":" << eta << ":" << phi << std::endl;
  // }

  // //-------
  // if(lumi == 885278 && event == 88527714){
  //   for( std::vector<pat::Jet>::const_iterator iJet = correctedJets.begin(); iJet!=correctedJets.end(); ++iJet ){
  //     double iPt = iJet->pt();
  //     // double iCSV = iJet->bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");
  //     double iCSV = iJet->chargedEmEnergyFraction();
  //     printf(" pt = %.3f,\t emFrac = %.3f \n", iPt, iCSV );

  //     double dR_lep_jet = 99.;
  //     bool isOverlap = false;
      
  //     TLorentzVector jet_p4;
  //     jet_p4.SetPxPyPzE(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      
  //     for( std::vector<pat::Electron>::const_iterator iEle = selectedElectrons.begin(); iEle != selectedElectrons.end(); iEle++ ){ 
  // 	TLorentzVector ele_p4;
  // 	ele_p4.SetPxPyPzE(iEle->px(),iEle->py(),iEle->pz(),iEle->energy());
  // 	double delta_tmp = jet_p4.DeltaR(ele_p4);
  // 	if(delta_tmp < 0.4){
  // 	  isOverlap = true;
  // 	  dR_lep_jet = delta_tmp;
  // 	  break;
  // 	} 
  //     }

  //     if ( isOverlap ) std::cout<<"is overlap and dR is: "<< dR_lep_jet  <<std::endl;
  //   }
  // }

  double lep1_pt  = -99;
  double lep1_eta = -99;
  double lep1_phi = -99;
  double lep1_iso = -99;
  int lep1_pdgId = -99;

  double lep2_pt  = 0;//-99;
  double lep2_eta = 0;//-99;
  double lep2_phi = 0;//-99;
  // double lep2_iso = 0;
  // double lep2_pdgId = 0;

  double m_ll = -99;
  bool passZVeto = false;
  bool oppositeCharge = false;
  bool twoEle = false;
  bool twoMu = false;
  bool EleMu = false;

  int lep1_charge = -99;
  int lep2_charge = -99;
  TLorentzVector sum_lepton_vect;
  sum_lepton_vect.SetPxPyPzE(0., 0., 0., 0.);

  if( numElectrons==1 && numElectrons_Tight==1 && numMuons==0){
    lep1_pt  = selectedElectrons.at(0).pt();
    lep1_eta = selectedElectrons.at(0).eta();
    lep1_phi = selectedElectrons.at(0).phi();
    lep1_iso = miniAODhelper.GetElectronRelIso(selectedElectrons.at(0), coneSize::R03, corrType::rhoEA);
    lep1_pdgId = selectedElectrons.at(0).pdgId();//11;
  }
  else if( numMuons==1 && numMuons_Tight==1 && numElectrons==0){
    lep1_pt  = selectedMuons.at(0).pt();
    lep1_eta = selectedMuons.at(0).eta();
    lep1_phi = selectedMuons.at(0).phi();
    lep1_iso = miniAODhelper.GetMuonRelIso(selectedMuons.at(0), coneSize::R04, corrType::deltaBeta);
    lep1_pdgId = selectedMuons.at(0).pdgId();//13;
  }

  else if( numElectrons==2 && numElectrons_Tight>0 && numMuons==0 ){
    twoEle = true;
    lep1_pt  = selectedElectrons.at(0).pt();
    lep1_eta = selectedElectrons.at(0).eta();
    lep1_phi = selectedElectrons.at(0).phi();
    lep1_charge = selectedElectrons.at(0).gsfTrack()->charge();

    TLorentzVector lep1P4;	  
    lep1P4.SetPxPyPzE(selectedElectrons.at(0).px(),selectedElectrons.at(0).py(),selectedElectrons.at(0).pz(),selectedElectrons.at(0).energy());
    sum_lepton_vect += lep1P4;

    lep2_pt  = selectedElectrons.at(1).pt();
    lep2_eta = selectedElectrons.at(1).eta();
    lep2_phi = selectedElectrons.at(1).phi();
    lep2_charge = selectedElectrons.at(1).gsfTrack()->charge();

    TLorentzVector lep2P4;	  
    lep2P4.SetPxPyPzE(selectedElectrons.at(1).px(),selectedElectrons.at(1).py(),selectedElectrons.at(1).pz(),selectedElectrons.at(1).energy());
    sum_lepton_vect += lep2P4;

  }
  else if (numElectrons==0 && numMuons==2 && numMuons_Tight>0){
    twoMu = true;
    lep1_pt  = selectedMuons.at(0).pt();
    lep1_eta = selectedMuons.at(0).eta();
    lep1_phi = selectedMuons.at(0).phi();
    lep1_charge = selectedMuons.at(0).muonBestTrack()->charge();

    TLorentzVector lep1P4;	  
    lep1P4.SetPxPyPzE(selectedMuons.at(0).px(),selectedMuons.at(0).py(),selectedMuons.at(0).pz(),selectedMuons.at(0).energy());
    sum_lepton_vect += lep1P4;

    lep2_pt  = selectedMuons.at(1).pt();
    lep2_eta = selectedMuons.at(1).eta();
    lep2_phi = selectedMuons.at(1).phi();
    lep2_charge = selectedMuons.at(1).muonBestTrack()->charge();

    TLorentzVector lep2P4;	  
    lep2P4.SetPxPyPzE(selectedMuons.at(1).px(),selectedMuons.at(1).py(),selectedMuons.at(1).pz(),selectedMuons.at(1).energy());
    sum_lepton_vect += lep2P4;

  }
  else if (numElectrons==1 && numMuons==1 && (numElectrons_Tight+numMuons_Tight)>0 ){
    EleMu = true;
    if(selectedElectrons.at(0).pt() > selectedMuons.at(0).pt()){
      lep1_pt  = selectedElectrons.at(0).pt();
      lep1_eta = selectedElectrons.at(0).eta();
      lep1_phi = selectedElectrons.at(0).phi();
      lep1_charge = selectedElectrons.at(0).gsfTrack()->charge();

      lep2_pt  = selectedMuons.at(0).pt();
      lep2_eta = selectedMuons.at(0).eta();
      lep2_phi = selectedMuons.at(0).phi();
      lep2_charge = selectedMuons.at(0).muonBestTrack()->charge();
    }
    else {
      lep2_pt  = selectedElectrons.at(0).pt();
      lep2_eta = selectedElectrons.at(0).eta();
      lep2_phi = selectedElectrons.at(0).phi();
      lep2_charge = selectedElectrons.at(0).gsfTrack()->charge();

      lep1_pt  = selectedMuons.at(0).pt();
      lep1_eta = selectedMuons.at(0).eta();
      lep1_phi = selectedMuons.at(0).phi();
      lep1_charge = selectedMuons.at(0).muonBestTrack()->charge();
    }

    TLorentzVector lep1P4;	  
    lep1P4.SetPxPyPzE(selectedElectrons.at(0).px(),selectedElectrons.at(0).py(),selectedElectrons.at(0).pz(),selectedElectrons.at(0).energy());
    sum_lepton_vect += lep1P4;

    TLorentzVector lep2P4;	  
    lep2P4.SetPxPyPzE(selectedMuons.at(0).px(),selectedMuons.at(0).py(),selectedMuons.at(0).pz(),selectedMuons.at(0).energy());
    sum_lepton_vect += lep2P4;

  }
  m_ll = sum_lepton_vect.M();
  passZVeto = (m_ll>76 && m_ll<106)? 0:1;
  if( (lep1_charge*lep2_charge)==-1 ) oppositeCharge = true; 
  twoEle = twoEle && oppositeCharge;
  twoMu = twoMu && oppositeCharge;
  EleMu = EleMu && oppositeCharge;

  double jet1_pt = ( numJets>=1 ) ? selectedJets.at(0).pt() : -99;
  double jet2_pt = ( numJets>=2 ) ? selectedJets.at(1).pt() : -99;
  double jet3_pt = ( numJets>=3 ) ? selectedJets.at(2).pt() : -99;
  double jet4_pt = ( numJets>=4 ) ? selectedJets.at(3).pt() : -99;

  // double jet1_CSVv2 = ( numJets>=1 ) ? selectedJets.at(0).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  // double jet2_CSVv2 = ( numJets>=2 ) ? selectedJets.at(1).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  // double jet3_CSVv2 = ( numJets>=3 ) ? selectedJets.at(2).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  // double jet4_CSVv2 = ( numJets>=4 ) ? selectedJets.at(3).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;

  double jet1_CSVv2 = ( numJets>=1 ) ? miniAODhelper.GetJetCSV(selectedJets.at(0), "pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  double jet2_CSVv2 = ( numJets>=2 ) ? miniAODhelper.GetJetCSV(selectedJets.at(1), "pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  double jet3_CSVv2 = ( numJets>=3 ) ? miniAODhelper.GetJetCSV(selectedJets.at(2), "pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  double jet4_CSVv2 = ( numJets>=4 ) ? miniAODhelper.GetJetCSV(selectedJets.at(3), "pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;

  int n_jets = numJets;
  int n_btags = numTags;

  int n_toptags = numTopTags;
  int n_higgstags = numHiggsTags;
  int ttHFCategory = additionalJetEventId;

  //---- comment out checks
	    // printf("===ref %d,%d,%d,1,0,%.2f,%+.2f,%+.2f,%.3f,%d,0,0,0,0,0,%.2f,%.2f,%.2f,%.2f,%+.3f,%+.3f,%+.3f,%+.3f,%.2f,%+.2f,%d,%d,0,%d\n",
	    // 	   run, lumi, event,
	    // 	   // is_SL, is_DL,
	    // 	   lep1_pt, lep1_eta, lep1_phi, lep1_iso, lep1_pdgId,
	    // 	   // lep2_pt, lep2_eta, lep2_phi, lep2_iso, lep2_pdgId,
	    // 	   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	    // 	   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	    // 	   // jet1_corr, jet2_corr, jet3_corr, jet4_corr,
	    // 	   // jet1_corrUp, jet2_corrUp, jet3_corrUp, jet4_corrUp,
	    // 	   // jet1_corrDown, jet2_corrDown, jet3_corrDown, jet4_corrDown,
	    // 	   MET_pt, MET_phi,
	    // 	   n_jets, n_btags,
	    // 	   // bWeight, bWeightLFUp, bWeightLFDown,
	    // 	   ttHFCategory);
	    // printf("numMuons_Tight:numMuons:numElectrons_Tight:numElectrons is %d : %d : %d : %d\n", numMuons_Tight, numMuons, numElectrons_Tight, numElectrons);
  //------------------


  int n_tth_syncex1_mu=0;
  h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);


  ////////////////////
  ///// synchronization with MEM
  ////////////////////
  bool synch_MEM = false;
  if( synch_MEM ){ 
    if (GenEventInfoWeight > 0) h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, 1);
    else h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, -1);

    if( numMuons_Tight==1 ){
      if (GenEventInfoWeight > 0) h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, 1);
      else h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, -1);

      if( numMuons==1 && numElectrons==0 ){
	if (GenEventInfoWeight > 0) h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, 1);
	else h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, -1);

	if( numJets>=4 ){
	  if (GenEventInfoWeight > 0) h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, 1);
	  else h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, -1);

	  if( numTags>=2 ){
	    if (GenEventInfoWeight > 0) h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, 1);
	    else h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++, -1);
	    // printf("===ref %d,%d,%d,1,0,%.2f,%+.2f,%+.2f,%.3f,%d,0,0,0,0,0,%.2f,%.2f,%.2f,%.2f,%+.3f,%+.3f,%+.3f,%+.3f,0,0,0,0,0,0,0,0,0,0,0,0,%.2f,%+.2f,%d,%d,0,0,0,%d\n",
	    printf("===ref %d,%d,%d,1,0,%.3f,%+.3f,%+.3f,%.3f,%d,0,0,0,0,0,%.3f,%.3f,%.3f,%.3f,%+.3f,%+.3f,%+.3f,%+.3f,%.3f,%+.3f,%d,%d,0,%d\n",
		   run, lumi, event,
		   // is_SL, is_DL,
		   lep1_pt, lep1_eta, lep1_phi, lep1_iso, lep1_pdgId,
		   // lep2_pt, lep2_eta, lep2_phi, lep2_iso, lep2_pdgId,
		   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		   // jet1_corr, jet2_corr, jet3_corr, jet4_corr,
		   // jet1_corrUp, jet2_corrUp, jet3_corrUp, jet4_corrUp,
		   // jet1_corrDown, jet2_corrDown, jet3_corrDown, jet4_corrDown,
		   MET_pt, MET_phi,
		   n_jets, n_btags,
		   // bWeight, bWeightLFUp, bWeightLFDown,
		   ttHFCategory);
	  }
	}
      } 
    }
  }


  int n_tth_syncex1_ele=0;
  h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);

  if( synch_MEM ){ 
    if (GenEventInfoWeight > 0) h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, 1);
    else h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, -1);

    if( numElectrons_Tight==1 ){
      if (GenEventInfoWeight > 0) h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, 1);
      else h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, -1);

      if( numElectrons==1 && numMuons==0 ){
	if (GenEventInfoWeight > 0) h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, 1);
	else h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, -1);

	if( numJets>=4 ){
	  if (GenEventInfoWeight > 0) h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, 1);
	  else h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, -1);

	  if( numTags>=2 ){
	    if (GenEventInfoWeight > 0) h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, 1);
	    else h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++, -1);
	    // printf("===ref %d,%d,%d,1,0,%.2f,%+.2f,%+.2f,%.3f,%d,0,0,0,0,0,%.2f,%.2f,%.2f,%.2f,%+.3f,%+.3f,%+.3f,%+.3f,0,0,0,0,0,0,0,0,0,0,0,0,%.2f,%+.2f,%d,%d,0,0,0,%d\n",
	    printf("===ref %d,%d,%d,1,0,%.3f,%+.3f,%+.3f,%.3f,%d,0,0,0,0,0,%.3f,%.3f,%.3f,%.3f,%+.3f,%+.3f,%+.3f,%+.3f,%.3f,%+.3f,%d,%d,0,%d\n",
		   run, lumi, event,
		   // is_SL, is_DL,
		   lep1_pt, lep1_eta, lep1_phi, lep1_iso, lep1_pdgId,
		   // lep2_pt, lep2_eta, lep2_phi, lep2_iso, lep2_pdgId,
		   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		   // jet1_corr, jet2_corr, jet3_corr, jet4_corr,
		   // jet1_corrUp, jet2_corrUp, jet3_corrUp, jet4_corrUp,
		   // jet1_corrDown, jet2_corrDown, jet3_corrDown, jet4_corrDown,
		   MET_pt, MET_phi,
		   n_jets, n_btags,
		   // bWeight, bWeightLFUp, bWeightLFDown,
		   ttHFCategory);

	  }
	}
      } 
    }
  }


  ////-======
  TString printOut = Form("%6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %+7.3f   %2d  %2d   %2d   %2d  %2d",
	   run, lumi, event,
	   lep1_pt, lep1_eta, lep1_phi,
	   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	   MET_pt,
	   n_jets, n_btags,
	   ttHFCategory,
	   n_toptags, n_higgstags);


  //------------------
  //// synchronization with LJ folks
  //////////
  // int n_tth_syncex1_mu=0;
  // h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);
  if( ! synch_MEM ){
  if( passSingleMuonTrigger ){
    std::cout << Form("muon %2d ", n_tth_syncex1_mu) << printOut << std::endl;
    // printf("muon %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %+7.3f   %2d  %2d   %2d   %2d  %2d\n",
    // 	   n_tth_syncex1_mu,
    // 	   run, lumi, event,
    // 	   lep1_pt, lep1_eta, lep1_phi,
    // 	   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
    // 	   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
    // 	   MET_pt,
    // 	   n_jets, n_btags,
    // 	   ttHFCategory,
    // 	   n_toptags, n_higgstags);
    h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);
    if( numMuons_Tight==1 ){
      std::cout << Form("muon %2d ", n_tth_syncex1_mu) << printOut << std::endl;
      h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);

      if( numMuons==1 && numElectrons==0 ){
	std::cout << Form("muon %2d ", n_tth_syncex1_mu) << printOut << std::endl;
	h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);

	if( numJets>=4 ){
	  std::cout << Form("muon %2d ", n_tth_syncex1_mu) << printOut << std::endl;
	  h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);

	  if( numTags>=2 ){
	    std::cout << Form("muon %2d ", n_tth_syncex1_mu) << printOut << std::endl;
	    h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);

	    if( numTopTags>=1 ){
	      std::cout << Form("muon %2d ", n_tth_syncex1_mu) << printOut << std::endl;
	      h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);

	      if( numHiggsTags>=1 ){
		std::cout << Form("muon %2d ", n_tth_syncex1_mu) << printOut << std::endl;
		h_tth_syncex1_mu->Fill(0.5+n_tth_syncex1_mu++);
	      }
	    }
	  }
	}
      }
    }
  }

  // int n_tth_syncex1_ele=0;
  // h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);
  if( passSingleElectronTrigger ){
    std::cout << Form("ele %2d ", n_tth_syncex1_ele) << printOut << std::endl;
    h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);

    if( numElectrons_Tight==1 ){
      std::cout << Form("ele %2d ", n_tth_syncex1_ele) << printOut << std::endl;
      h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);

      if( numElectrons==1 && numMuons==0 ){
	std::cout << Form("ele %2d ", n_tth_syncex1_ele) << printOut << std::endl;
	h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);

	if( numJets>=4 ){
	  std::cout << Form("ele %2d ", n_tth_syncex1_ele) << printOut << std::endl;
	  h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);

	  if( numTags>=2 ){
	    std::cout << Form("ele %2d ", n_tth_syncex1_ele) << printOut << std::endl;
	    h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);

	    if( numTopTags>=1 ){
	      std::cout << Form("ele %2d ", n_tth_syncex1_ele) << printOut << std::endl;
	      h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);

	      if( numHiggsTags>=1 ){
		std::cout << Form("ele %2d ", n_tth_syncex1_ele) << printOut << std::endl;
		h_tth_syncex1_ele->Fill(0.5+n_tth_syncex1_ele++);
	      }
	    }
	  }
	}
      }
    }
  }

  } /// LJ synch

  //-----------double electrons
  int n_tth_syncex1_diele=0;
  h_tth_syncex1_diele->Fill(0.5+n_tth_syncex1_diele++);
  if( false ){ //passDiElectronTrigger ){
    printf("ee %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	   n_tth_syncex1_diele,
	   run, lumi, event,
	   lep1_pt, lep1_eta, lep1_phi,
	   lep2_pt, lep2_eta, lep2_phi,
	   m_ll, MET_pt,
	   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	   n_jets, n_btags);
    h_tth_syncex1_diele->Fill(0.5+n_tth_syncex1_diele++);
    if( twoEle ){
      printf("ee %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	     n_tth_syncex1_diele,
	     run, lumi, event,
	     lep1_pt, lep1_eta, lep1_phi,
	     lep2_pt, lep2_eta, lep2_phi,
	     m_ll, MET_pt,
	     jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	     jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	     n_jets, n_btags);
      h_tth_syncex1_diele->Fill(0.5+n_tth_syncex1_diele++);
      if( m_ll>20 ){
	printf("ee %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	       n_tth_syncex1_diele,
	       run, lumi, event,
	       lep1_pt, lep1_eta, lep1_phi,
	       lep2_pt, lep2_eta, lep2_phi,
	       m_ll, MET_pt,
	       jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	       jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	       n_jets, n_btags);
	h_tth_syncex1_diele->Fill(0.5+n_tth_syncex1_diele++);
	if( passZVeto ){
	  printf("ee %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		 n_tth_syncex1_diele,
		 run, lumi, event,
		 lep1_pt, lep1_eta, lep1_phi,
		 lep2_pt, lep2_eta, lep2_phi,
		 m_ll, MET_pt,
		 jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		 jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		 n_jets, n_btags);
	  h_tth_syncex1_diele->Fill(0.5+n_tth_syncex1_diele++);
	  if( numJets>=2 ){
	    printf("ee %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		   n_tth_syncex1_diele,
		   run, lumi, event,
		   lep1_pt, lep1_eta, lep1_phi,
		   lep2_pt, lep2_eta, lep2_phi,
		   m_ll, MET_pt,
		   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		   n_jets, n_btags);
	    h_tth_syncex1_diele->Fill(0.5+n_tth_syncex1_diele++);
	    if( MET_pt > 40 ){
	      printf("ee %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		     n_tth_syncex1_diele,
		     run, lumi, event,
		     lep1_pt, lep1_eta, lep1_phi,
		     lep2_pt, lep2_eta, lep2_phi,
		     m_ll, MET_pt,
		     jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		     jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		     n_jets, n_btags);
	      h_tth_syncex1_diele->Fill(0.5+n_tth_syncex1_diele++);
	      if( numTags>=1 ){
		printf("ee %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		       n_tth_syncex1_diele,
		       run, lumi, event,
		       lep1_pt, lep1_eta, lep1_phi,
		       lep2_pt, lep2_eta, lep2_phi,
		       m_ll, MET_pt,
		       jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		       jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		       n_jets, n_btags);
		h_tth_syncex1_diele->Fill(0.5+n_tth_syncex1_diele++);
	      }
	    }
	  }
	}
      }
    }
  }

  //-----------double muons
  int n_tth_syncex1_dimuon=0;
  h_tth_syncex1_dimuon->Fill(0.5+n_tth_syncex1_dimuon++);
  if( false ){ //passDiMuonTrigger ){
    printf("mumu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	   n_tth_syncex1_dimuon,
	   run, lumi, event,
	   lep1_pt, lep1_eta, lep1_phi,
	   lep2_pt, lep2_eta, lep2_phi,
	   m_ll, MET_pt,
	   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	   n_jets, n_btags);
    h_tth_syncex1_dimuon->Fill(0.5+n_tth_syncex1_dimuon++);
    if( twoMu ){
      printf("mumu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	     n_tth_syncex1_dimuon,
	     run, lumi, event,
	     lep1_pt, lep1_eta, lep1_phi,
	     lep2_pt, lep2_eta, lep2_phi,
	     m_ll, MET_pt,
	     jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	     jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	     n_jets, n_btags);
      h_tth_syncex1_dimuon->Fill(0.5+n_tth_syncex1_dimuon++);
      if( m_ll>20 ){
	printf("mumu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	       n_tth_syncex1_dimuon,
	       run, lumi, event,
	       lep1_pt, lep1_eta, lep1_phi,
	       lep2_pt, lep2_eta, lep2_phi,
	       m_ll, MET_pt,
	       jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	       jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	       n_jets, n_btags);
	h_tth_syncex1_dimuon->Fill(0.5+n_tth_syncex1_dimuon++);
	if( passZVeto ){
	  printf("mumu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		 n_tth_syncex1_dimuon,
		 run, lumi, event,
		 lep1_pt, lep1_eta, lep1_phi,
		 lep2_pt, lep2_eta, lep2_phi,
		 m_ll, MET_pt,
		 jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		 jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		 n_jets, n_btags);
	  h_tth_syncex1_dimuon->Fill(0.5+n_tth_syncex1_dimuon++);
	  if( numJets>=2 ){
	    printf("mumu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		   n_tth_syncex1_dimuon,
		   run, lumi, event,
		   lep1_pt, lep1_eta, lep1_phi,
		   lep2_pt, lep2_eta, lep2_phi,
		   m_ll, MET_pt,
		   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		   n_jets, n_btags);
	    h_tth_syncex1_dimuon->Fill(0.5+n_tth_syncex1_dimuon++);
	    if( MET_pt > 40 ){
	      printf("mumu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		     n_tth_syncex1_dimuon,
		     run, lumi, event,
		     lep1_pt, lep1_eta, lep1_phi,
		     lep2_pt, lep2_eta, lep2_phi,
		     m_ll, MET_pt,
		     jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		     jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		     n_jets, n_btags);
	      h_tth_syncex1_dimuon->Fill(0.5+n_tth_syncex1_dimuon++);
	      if( numTags>=1 ){
		printf("mumu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		       n_tth_syncex1_dimuon,
		       run, lumi, event,
		       lep1_pt, lep1_eta, lep1_phi,
		       lep2_pt, lep2_eta, lep2_phi,
		       m_ll, MET_pt,
		       jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		       jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		       n_jets, n_btags);
		h_tth_syncex1_dimuon->Fill(0.5+n_tth_syncex1_dimuon++);
	      }
	    }
	  }
	}
      }
    }
  }


  //-----------ele mus
  int n_tth_syncex1_elemuon=0;
  h_tth_syncex1_elemuon->Fill(0.5+n_tth_syncex1_elemuon++);
  if( false ){ //passEleMuonTrigger ){
    printf("emu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	   n_tth_syncex1_elemuon,
	   run, lumi, event,
	   lep1_pt, lep1_eta, lep1_phi,
	   lep2_pt, lep2_eta, lep2_phi,
	   m_ll, MET_pt,
	   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	   n_jets, n_btags);
    h_tth_syncex1_elemuon->Fill(0.5+n_tth_syncex1_elemuon++);
    if( EleMu ){
      printf("emu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	     n_tth_syncex1_elemuon,
	     run, lumi, event,
	     lep1_pt, lep1_eta, lep1_phi,
	     lep2_pt, lep2_eta, lep2_phi,
	     m_ll, MET_pt,
	     jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	     jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	     n_jets, n_btags);
      h_tth_syncex1_elemuon->Fill(0.5+n_tth_syncex1_elemuon++);
      if( m_ll>20 ){
	printf("emu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	       n_tth_syncex1_elemuon,
	       run, lumi, event,
	       lep1_pt, lep1_eta, lep1_phi,
	       lep2_pt, lep2_eta, lep2_phi,
	       m_ll, MET_pt,
	       jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	       jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	       n_jets, n_btags);
	h_tth_syncex1_elemuon->Fill(0.5+n_tth_syncex1_elemuon++);
	// if( passZVeto ){
	//   printf("emu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	// 	 n_tth_syncex1_elemuon,
	// 	 run, lumi, event,
	// 	 lep1_pt, lep1_eta, lep1_phi,
	// 	 lep2_pt, lep2_eta, lep2_phi,
	// 	 m_ll, MET_pt,
	// 	 jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	// 	 jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	// 	 n_jets, n_btags);
	//   h_tth_syncex1_elemuon->Fill(0.5+n_tth_syncex1_elemuon++);
	  if( numJets>=2 ){
	    printf("emu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		   n_tth_syncex1_elemuon,
		   run, lumi, event,
		   lep1_pt, lep1_eta, lep1_phi,
		   lep2_pt, lep2_eta, lep2_phi,
		   m_ll, MET_pt,
		   jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		   jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		   n_jets, n_btags);
	    h_tth_syncex1_elemuon->Fill(0.5+n_tth_syncex1_elemuon++);
	    // if( MET_pt > 40 ){
	    //   printf("emu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
	    // 	     n_tth_syncex1_elemuon,
	    // 	     run, lumi, event,
	    // 	     lep1_pt, lep1_eta, lep1_phi,
	    // 	     lep2_pt, lep2_eta, lep2_phi,
	    // 	     m_ll, MET_pt,
	    // 	     jet1_pt, jet2_pt, jet3_pt, jet4_pt,
	    // 	     jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
	    // 	     n_jets, n_btags);
	    //   h_tth_syncex1_elemuon->Fill(0.5+n_tth_syncex1_elemuon++);
	      if( numTags>=1 ){
		printf("emu %2d %6d %8d %10d   %6.2f %+4.2f %+4.2f   %6.2f %+4.2f %+4.2f   %6.2f  %6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %2d  %2d\n",
		       n_tth_syncex1_elemuon,
		       run, lumi, event,
		       lep1_pt, lep1_eta, lep1_phi,
		       lep2_pt, lep2_eta, lep2_phi,
		       m_ll, MET_pt,
		       jet1_pt, jet2_pt, jet3_pt, jet4_pt,
		       jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
		       n_jets, n_btags);
		h_tth_syncex1_elemuon->Fill(0.5+n_tth_syncex1_elemuon++);
	      }
	    }
	  }
	}
      }
  //   }
  // }

  // if( numElectrons==1 ){

  // printf("%6d %8d %10d   %6.2f %+6.2f %+6.2f   %6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   %+7.3f   %2d  %2d   %2d   %2d  %2d\n",
  // 	 run, lumi, event,
  // 	 lep1_pt, lep1_eta, lep1_phi,
  // 	 jet1_pt, jet2_pt, jet3_pt, jet4_pt,
  // 	 jet1_CSVv2, jet2_CSVv2, jet3_CSVv2, jet4_CSVv2,
  // 	 n_jets, n_btags,
  // 	 n_toptags, n_higgstags);

  // }
}


// ------------ method called once each job just before starting event loop  ------------
void 
TTHSyncExercise::beginJob()
{

  TH1::SetDefaultSumw2(true);

  numEvents_ = 0;
  numPos_  = 0;
  numNeg_  = 0;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
TTHSyncExercise::endJob() 
{
  std::cout << "----------------" << std::endl;
  std::cout << "total = " << numEvents_ << std::endl;
  std::cout << "pos = " << numPos_ << std::endl;
  std::cout << "neg = " << numNeg_ << std::endl;
  // report results of sync exercises
  int nbins_tth_syncex1_mu = h_tth_syncex1_mu->GetNbinsX();
  std::cout << "***************************************************************" << std::endl;
  std::cout << "\t Synchronization exercise 1" << std::endl;
  printf("Selection \t Number of events\n");
  for( int iBin=0; iBin<nbins_tth_syncex1_mu; iBin++ )
    printf("%s\t %.0f\n", h_tth_syncex1_mu->GetXaxis()->GetBinLabel(iBin+1), h_tth_syncex1_mu->GetBinContent(iBin+1));

  int nbins_tth_syncex1_ele = h_tth_syncex1_ele->GetNbinsX();
  std::cout << "***************************************************************" << std::endl;
  std::cout << "\t Synchronization exercise 2" << std::endl;
  printf("Selection \t Number of events\n");
  for( int iBin=0; iBin<nbins_tth_syncex1_ele; iBin++ )
    printf("%s\t %.0f\n", h_tth_syncex1_ele->GetXaxis()->GetBinLabel(iBin+1), h_tth_syncex1_ele->GetBinContent(iBin+1));
  std::cout << "***************************************************************" << std::endl;

  int nbins_tth_syncex1_elemuon = h_tth_syncex1_elemuon->GetNbinsX();
  std::cout << "***************************************************************" << std::endl;
  std::cout << "\t Synchronization exercise 3" << std::endl;
  printf("Selection \t Number of events\n");
  for( int iBin=0; iBin<nbins_tth_syncex1_elemuon; iBin++ )
    printf("%s\t %.0f\n", h_tth_syncex1_elemuon->GetXaxis()->GetBinLabel(iBin+1), h_tth_syncex1_elemuon->GetBinContent(iBin+1));
  std::cout << "***************************************************************" << std::endl;

  int nbins_tth_syncex1_dimuon = h_tth_syncex1_dimuon->GetNbinsX();
  std::cout << "***************************************************************" << std::endl;
  std::cout << "\t Synchronization exercise 4" << std::endl;
  printf("Selection \t Number of events\n");
  for( int iBin=0; iBin<nbins_tth_syncex1_dimuon; iBin++ )
    printf("%s\t %.0f\n", h_tth_syncex1_dimuon->GetXaxis()->GetBinLabel(iBin+1), h_tth_syncex1_dimuon->GetBinContent(iBin+1));
  std::cout << "***************************************************************" << std::endl;

  int nbins_tth_syncex1_diele = h_tth_syncex1_diele->GetNbinsX();
  std::cout << "***************************************************************" << std::endl;
  std::cout << "\t Synchronization exercise 5" << std::endl;
  printf("Selection \t Number of events\n");
  for( int iBin=0; iBin<nbins_tth_syncex1_diele; iBin++ )
    printf("%s\t %.0f\n", h_tth_syncex1_diele->GetXaxis()->GetBinLabel(iBin+1), h_tth_syncex1_diele->GetBinContent(iBin+1));
  std::cout << "***************************************************************" << std::endl;




  // report on triggers fired
  if( dumpHLT_ ){
    std::cout << "***************************************************************" << std::endl;
    std::cout << "  Summary for HLT: Total number of events = " << numEvents_ << std::endl;
    for( std::map<std::string, int>::const_iterator iter = hlt_cppath_.begin(); iter != hlt_cppath_.end(); iter++){
      std::string name = iter->first;
      double eff = 100*double(iter->second)/double(numEvents_);
      printf("\t %s,\t %d,\t %.1f \n",name.c_str(), iter->second, eff);
    }
    std::cout << "***************************************************************" << std::endl;
    std::cout << "  Summary for Filters: Total number of events = " << numEvents_ << std::endl;
    for( std::map<std::string, int>::const_iterator iter = flt_cppath_.begin(); iter != flt_cppath_.end(); iter++){
      std::string name = iter->first;
      double eff = 100*double(iter->second)/double(numEvents_);
      printf("\t %s,\t %d,\t %.1f \n",name.c_str(), iter->second, eff);
    }
    std::cout << "***************************************************************" << std::endl;
  }

  std::cout << "***************************************************************" << std::endl;
  std::cout << "  Total number of events = " << numEvents_ << std::endl;
  std::cout << "***************************************************************" << std::endl;

  // syncOutputFile.close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
TTHSyncExercise::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
void 
TTHSyncExercise::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

  bool hltchanged = true;
  if (!hlt_config_.init(iRun, iSetup, hltTag, hltchanged)) {
    std::cout << "Warning, didn't find trigger process HLT,\t" << hltTag << std::endl;
    return;
  }
  bool filterchanged = true;
  if (!filter_config_.init(iRun, iSetup, filterTag, filterchanged)) {
    std::cout << "Warning, didn't find filter process HLT,\t" << filterTag << std::endl;
    return;
  }

  // Zero out map
  std::vector<std::string> triggerNames = hlt_config_.triggerNames();
  std::vector<std::string> filterNames  = filter_config_.triggerNames();

  hlt_triggerNames_.clear();
  flt_filterNames_.clear();

  hlt_triggerNames_.push_back("All");
  std::string prefix = "HLT_";
  for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){
    std::string name = triggerNames[iPath];
    std::string pathNameNoVer = hlt_config_.removeVersion(name);
    hlt_cppath_[pathNameNoVer] = 0;
    if( name.compare(0, prefix.length(), prefix) == 0 ) hlt_triggerNames_.push_back(pathNameNoVer);
  }

  flt_filterNames_.push_back("All");
  for( unsigned int iPath=0; iPath<filterNames.size(); iPath++ ){
    std::string name = filterNames[iPath];
    std::string pathNameNoVer = filter_config_.removeVersion(name);
    flt_cppath_[pathNameNoVer] = 0;
    flt_filterNames_.push_back(pathNameNoVer);
  }

  unsigned int numHLT = hlt_triggerNames_.size();
  unsigned int numFLT = flt_filterNames_.size();

  h_hlt = fs_->make<TH1D>("h_hlt",";HLT path", numHLT , 0 , numHLT );
  h_flt = fs_->make<TH1D>("h_flt",";Filter path", numFLT , 0 , numFLT );

  for( unsigned int iPath=0; iPath<numHLT; iPath++ ){
    std::string pathNameNoVer = hlt_triggerNames_[iPath];
    int jPath = iPath+1;
    if( h_hlt ){
      TAxis * axis = h_hlt->GetXaxis();
      if( axis ) axis->SetBinLabel(jPath, pathNameNoVer.c_str());
    }
  }

  for( unsigned int iPath=0; iPath<numFLT; iPath++ ){
    std::string pathNameNoVer = flt_filterNames_[iPath];
    int jPath = iPath+1;
    if( h_flt ){
      TAxis * axis = h_flt->GetXaxis();
      if( axis ) axis->SetBinLabel(jPath, pathNameNoVer.c_str());
    }
  }

}

// ------------ method called when ending the processing of a run  ------------
/*
void 
TTHSyncExercise::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
TTHSyncExercise::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
TTHSyncExercise::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTHSyncExercise::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTHSyncExercise);
