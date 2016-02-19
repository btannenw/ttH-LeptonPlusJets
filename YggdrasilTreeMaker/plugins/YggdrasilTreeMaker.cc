// -*- C++ -*-
//
// Package:    ttH-LeptonPlusJets/YggdrasilTreeMaker
// Class:      YggdrasilTreeMaker
// 
/**\class YggdrasilTreeMaker YggdrasilTreeMaker.cc ttH-LeptonPlusJets/YggdrasilTreeMaker/plugins/YggdrasilTreeMaker.cc
 Description: [one line class summary]
 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Darren Puigh
//         Created:  Fri, 12 Sep 2014 16:58:04 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"

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

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"
#include "MiniAOD/MiniAODHelper/interface/TopTagger.h"
#include "BoostedTTH/BoostedAnalyzer/interface/BoostedUtils.hpp"
#include "MiniAOD/MiniAODHelper/interface/HiggsTagger.h"
#include "ttH-LeptonPlusJets/AnalysisCode/interface/YggdrasilEventVars.h"

// #include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

#include "MiniAOD/MiniAODHelper/interface/BDTvars.h"

//
// class declaration
//

class YggdrasilTreeMaker : public edm::EDAnalyzer {
   public:
      explicit YggdrasilTreeMaker(const edm::ParameterSet&);
      ~YggdrasilTreeMaker();

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
  edm::EDGetTokenT<int> genTtbarIdToken_;
  // Histogram for ttbar event categorization ID including information about b jets from top in acceptance
  TH1* h_ttbarId_;
        
  // Histogram for ttbar event categorization ID based on additional jets only
  TH1* h_ttbarAdditionalJetId_;

  
  // Input tags
  const edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  
  const edm::EDGetTokenT<std::vector<int> > genBHadJetIndexToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadFlavourToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadFromTopWeakDecayToken_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> > genBHadPlusMothersToken_;
  const edm::EDGetTokenT<std::vector<std::vector<int> > > genBHadPlusMothersIndicesToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadIndexToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadLeptonHadronIndexToken_;
  const edm::EDGetTokenT<std::vector<int> > genBHadLeptonViaTauToken_;
  
  const edm::EDGetTokenT<std::vector<int> > genCHadJetIndexToken_;
  const edm::EDGetTokenT<std::vector<int> > genCHadFlavourToken_;
  const edm::EDGetTokenT<std::vector<int> > genCHadFromTopWeakDecayToken_;
  const edm::EDGetTokenT<std::vector<int> > genCHadBHadronIdToken_;

  //----------

  edm::EDGetTokenT <edm::TriggerResults> triggerResultsToken;
  edm::EDGetTokenT <edm::TriggerResults> filterResultsToken;

  // new MVAelectron
  edm::EDGetTokenT< edm::View<pat::Electron> > EDMElectronsToken;
  // MVA values and categories
  edm::EDGetTokenT<edm::ValueMap<float> > EDMeleMVAvaluesToken;
  edm::EDGetTokenT<edm::ValueMap<int> > EDMeleMVAcategoriesToken;

  edm::EDGetTokenT <reco::VertexCollection> vertexToken;
  edm::EDGetTokenT <pat::ElectronCollection> electronToken;
  edm::EDGetTokenT <pat::MuonCollection> muonToken;
  edm::EDGetTokenT <pat::JetCollection> jetToken;
  edm::EDGetTokenT <pat::METCollection> metToken;

  // edm::EDGetTokenT< boosted::HEPTopJetCollection > topJetsToken;
  // edm::EDGetTokenT< boosted::SubFilterJetCollection > subFilterJetsToken;

  edm::EDGetTokenT <pat::PackedCandidateCollection> packedpfToken;

  edm::EDGetTokenT <reco::BeamSpot> beamspotToken;
  edm::EDGetTokenT <double> rhoToken;
  edm::EDGetTokenT <reco::GenParticleCollection> mcparicleToken;
  edm::EDGetTokenT <std::vector< PileupSummaryInfo > > puInfoToken;

  edm::EDGetTokenT <GenEventInfoProduct> genInfoProductToken;

  edm::EDGetTokenT <pat::JetCollection> tempjetToken;
  
  // edm::EDGetTokenT <reco::ConversionCollection> EDMConversionCollectionToken;
  edm::EDGetTokenT< boosted::BoostedJetCollection > EDMBoostedJetsToken;
  

  HLTConfigProvider hlt_config_;

  bool verbose_;

  bool isLJ_;

  int insample_;

  std::string mySample_sampleName_;
  double mySample_xSec_;
  double mySample_nGen_;
  double intLumi_;

  std::string hltTag;
  std::string filterTag;

  int neventsFillTree[rNumSys];

  int nevents;
  double nevents_wgt;

  int nevents_clean;
  double nevents_clean_wgt;

  int nevents_1l, nevents_1l_ele, nevents_1l_mu;
  double nevents_1l_wgt, nevents_1l_ele_wgt, nevents_1l_mu_wgt;

  int nevents_1l_4j, nevents_1l_4j_ele, nevents_1l_4j_mu;
  double nevents_1l_4j_wgt, nevents_1l_4j_ele_wgt, nevents_1l_4j_mu_wgt;

  int nevents_1l_4j_2t, nevents_1l_4j_2t_ele, nevents_1l_4j_2t_mu;
  double nevents_1l_4j_2t_wgt, nevents_1l_4j_2t_ele_wgt, nevents_1l_4j_2t_mu_wgt;

  int nevents_2l, nevents_2l_2e, nevents_2l_2m, nevents_2l_em;
  double nevents_2l_wgt, nevents_2l_2e_wgt, nevents_2l_2m_wgt, nevents_2l_em_wgt;

  int nevents_2l_2j, nevents_2l_2j_2e, nevents_2l_2j_2m, nevents_2l_2j_em;
  double nevents_2l_2j_wgt, nevents_2l_2j_2e_wgt, nevents_2l_2j_2m_wgt, nevents_2l_2j_em_wgt;

  int nevents_2l_2j_2t, nevents_2l_2j_2t_2e, nevents_2l_2j_2t_2m, nevents_2l_2j_2t_em;
  double nevents_2l_2j_2t_wgt, nevents_2l_2j_2t_2e_wgt, nevents_2l_2j_2t_2m_wgt, nevents_2l_2j_2t_em_wgt;

  TTree *worldTree;
  yggdrasilEventVars *eve; 

  // EGammaMvaEleEstimatorCSA14* myMVATrig;
 
  MiniAODHelper miniAODhelper;

  BDTvars bdtVARS;
  
  TopTagger toptagger;

  const bool isMC ;

};

//
// constants, enums and typedefs
//
typedef std::vector<std::vector<double> >      vvdouble;
typedef std::vector<std::vector<std::string> > vvstring;
typedef std::vector<double>                    vdouble;
typedef std::vector<string>                    vstring;
typedef std::vector<bool>                      vbool;
typedef std::vector<int>                       vint;
typedef std::vector< TLorentzVector >          vecTLorentzVector;

//
// static data member definitions
//

//
// constructors and destructor
//
YggdrasilTreeMaker::YggdrasilTreeMaker(const edm::ParameterSet& iConfig):
  isMC(iConfig.getParameter<std::string>("isMC") == "MC" )
{
   //now do what ever initialization is needed
  verbose_ = false;
  isLJ_ = true;

  hltTag = "HLT";
  if( isMC ){
    filterTag = "PAT";
  }else{
    filterTag = "HLT";
  }
  triggerResultsToken = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), hltTag));
  filterResultsToken = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), filterTag));

  // new MVAelectron
  EDMElectronsToken = consumes< edm::View<pat::Electron> >(edm::InputTag("slimmedElectrons","",""));
  EDMeleMVAvaluesToken           = consumes<edm::ValueMap<float> >(edm::InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values",""));
  EDMeleMVAcategoriesToken       = consumes<edm::ValueMap<int> >(edm::InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories",""));

  vertexToken = consumes <reco::VertexCollection> (edm::InputTag(std::string("offlineSlimmedPrimaryVertices")));
  electronToken = consumes <pat::ElectronCollection> (edm::InputTag(std::string("slimmedElectrons")));
  muonToken = consumes <pat::MuonCollection> (edm::InputTag(std::string("slimmedMuons")));
  jetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("slimmedJets")));
  metToken = consumes <pat::METCollection> (edm::InputTag(std::string("slimmedMETs")));

  // topJetsToken    = consumes< boosted::HEPTopJetCollection >(edm::InputTag("HEPTopJetsPFMatcher","heptopjets","p"));
  // subFilterJetsToken = consumes< boosted::SubFilterJetCollection >(edm::InputTag("CA12JetsCA3FilterjetsPFMatcher","subfilterjets","p"));

  packedpfToken = consumes <pat::PackedCandidateCollection> (edm::InputTag(std::string("packedPFCandidates")));

  beamspotToken = consumes <reco::BeamSpot> (edm::InputTag(std::string("offlineBeamSpot")));
  rhoToken = consumes <double> (edm::InputTag(std::string("fixedGridRhoFastjetAll")));
  mcparicleToken = consumes <reco::GenParticleCollection> (edm::InputTag(std::string("prunedGenParticles")));
  puInfoToken = consumes <std::vector< PileupSummaryInfo > > (edm::InputTag(std::string("slimmedAddPileupInfo")));

  if( isMC ){
    genInfoProductToken = consumes <GenEventInfoProduct> (edm::InputTag(std::string("generator")));
  }

  tempjetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("slimmedJets")));
  
  // EDMConversionCollectionToken = consumes <reco::ConversionCollection > (edm::InputTag("reducedEgamma","reducedConversions",""));
  EDMBoostedJetsToken     = consumes< boosted::BoostedJetCollection >(edm::InputTag("BoostedJetMatcher","boostedjets","p"));

  if( isMC ){
  genTtbarIdToken_ = consumes<int>( edm::InputTag( "categorizeGenTtbar", "genTtbarId","" ) )  ;
  }

  edm::Service<TFileService> fs_;
  worldTree = fs_->make<TTree>("worldTree", "worldTree");
  eve=0; 
  worldTree->Branch("eve.", "yggdrasilEventVars", &eve, 8000, 1);


  for( int i=0; i<rNumSys; i++ ) neventsFillTree[i] = 0;

  nevents=0;
  nevents_wgt=0;

  nevents_clean=0;
  nevents_clean_wgt=0;

  nevents_1l=0, nevents_1l_ele=0, nevents_1l_mu=0;
  nevents_1l_wgt=0, nevents_1l_ele_wgt=0, nevents_1l_mu_wgt=0;

  nevents_1l_4j=0, nevents_1l_4j_ele=0, nevents_1l_4j_mu=0;
  nevents_1l_4j_wgt=0, nevents_1l_4j_ele_wgt=0, nevents_1l_4j_mu_wgt=0;

  nevents_1l_4j_2t=0, nevents_1l_4j_2t_ele=0, nevents_1l_4j_2t_mu=0;
  nevents_1l_4j_2t_wgt=0, nevents_1l_4j_2t_ele_wgt=0, nevents_1l_4j_2t_mu_wgt=0;

  nevents_2l=0, nevents_2l_2e=0, nevents_2l_2m=0, nevents_2l_em=0;
  nevents_2l_wgt=0, nevents_2l_2e_wgt=0, nevents_2l_2m_wgt=0, nevents_2l_em_wgt=0;

  nevents_2l_2j=0, nevents_2l_2j_2e=0, nevents_2l_2j_2m=0, nevents_2l_2j_em=0;
  nevents_2l_2j_wgt=0, nevents_2l_2j_2e_wgt=0, nevents_2l_2j_2m_wgt=0, nevents_2l_2j_em_wgt=0;

  nevents_2l_2j_2t=0, nevents_2l_2j_2t_2e=0, nevents_2l_2j_2t_2m=0, nevents_2l_2j_2t_em=0;
  nevents_2l_2j_2t_wgt=0, nevents_2l_2j_2t_2e_wgt=0, nevents_2l_2j_2t_2m_wgt=0, nevents_2l_2j_2t_em_wgt=0;


  std::string era = "2012_53x";
  insample_ = 2500;

  mySample_sampleName_ = "TTJets_MSDecaysCKM_central_Tune4C_13TeV_madgraph_PU20bx25_POSTLS170_V5_v1";
  mySample_xSec_ = 689.1;
  mySample_nGen_ = 25474122;
  intLumi_ = 20000;

  analysisType::analysisType iAnalysisType = analysisType::LJ;
  bool isData = ( insample_<0 );

  miniAODhelper.SetUp(era, insample_, iAnalysisType, isData);

  miniAODhelper.SetJetCorrectorUncertainty();

   // miniAODhelper.SetUpElectronMVA("MiniAOD/MiniAODHelper/data/ElectronMVA/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml","MiniAOD/MiniAODHelper/data/ElectronMVA/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml","MiniAOD/MiniAODHelper/data/ElectronMVA/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml");
  

  toptagger = TopTagger(TopTag::Likelihood, TopTag::CSV, "toplikelihoodtaggerhistos.root");



  //std::vector<std::string> myManualCatWeigths;
  //myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_50ns_EB_BDT.weights.xml");
  //myManualCatWeigths.push_back("EgammaAnalysis/ElectronTools/data/CSA14/TrigIDMVA_50ns_EE_BDT.weights.xml");

  //std::vector<std::string> myManualCatWeigthsTrig;
  //std::string the_path;
  //for (unsigned i  = 0 ; i < myManualCatWeigths.size() ; i++){
   // the_path = edm::FileInPath ( myManualCatWeigths[i] ).fullPath();
    //myManualCatWeigthsTrig.push_back(the_path);
  //}
    
  // myMVATrig = new EGammaMvaEleEstimatorCSA14();
  // myMVATrig->initialize("BDT",
  // 			EGammaMvaEleEstimatorCSA14::kTrig,
  // 			true,
  // 			myManualCatWeigthsTrig);
}


YggdrasilTreeMaker::~YggdrasilTreeMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

/*
void getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, float &aplanarity, float &sphericity);
void getFox(vecTLorentzVector jets, float &h0, float &h1, float &h2, float &h3, float &h4);
double getBestHiggsMass(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, vdouble btag, double &minChi, double &dRbb, TLorentzVector &bjet1, TLorentzVector &bjet2, vecTLorentzVector loose_jets, vdouble loose_btag);
*/


// ------------ method called for each event  ------------
void
YggdrasilTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   
   eve->initialize();


  double minTightLeptonPt = ( isLJ_ ) ? 30. : 20.;
  double looseLeptonPt = 10.;

  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(vertexToken,vtxHandle);
  reco::VertexCollection vtxs = *vtxHandle;

  /// old way of getting electrons
  // edm::Handle<pat::ElectronCollection> electrons;
  // iEvent.getByToken(electronToken,electrons);
  //// MVAelectrons
  edm::Handle< edm::View<pat::Electron> > h_electrons;
  iEvent.getByToken( EDMElectronsToken,h_electrons );
  // add electron mva info to electrons
  edm::Handle<edm::ValueMap<float> > h_mvaValues; 
  edm::Handle<edm::ValueMap<int> > h_mvaCategories;
  iEvent.getByToken(EDMeleMVAvaluesToken,h_mvaValues);
  iEvent.getByToken(EDMeleMVAcategoriesToken,h_mvaCategories);  
  std::vector<pat::Electron> electrons = miniAODhelper.GetElectronsWithMVAid(h_electrons,h_mvaValues,h_mvaCategories);

  ////
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
  if( isMC ){
    iEvent.getByToken(mcparicleToken,mcparticles);
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken,rhoHandle);
  ////------- set up rho for lepton effArea Isolation correction
  double rho_event = ( (rhoHandle.isValid()) ) ? *rhoHandle : -99;
  miniAODhelper.SetRho(rho_event);

  edm::Handle<std::vector< PileupSummaryInfo > > PupInfo;
  iEvent.getByToken(puInfoToken,PupInfo);

  // edm::Handle<boosted::HEPTopJetCollection> h_heptopjet;
  // iEvent.getByToken( topJetsToken,h_heptopjet);

  // edm::Handle<boosted::SubFilterJetCollection> h_subfilterjet;
  // iEvent.getByToken( subFilterJetsToken,h_subfilterjet );

  edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
  if( isMC ){
  iEvent.getByToken(genInfoProductToken,GenEventInfoHandle);
  }

  edm::Handle<boosted::BoostedJetCollection> h_boostedjet;
  iEvent.getByToken( EDMBoostedJetsToken,h_boostedjet);
  
  //  edm::Handle<reco::ConversionCollection> h_conversioncollection;
  // iEvent.getByToken( EDMConversionCollectionToken,h_conversioncollection );

  double GenEventInfoWeight = isMC ? GenEventInfoHandle.product()->weight() : 1.0 ;


  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken, triggerResults);
  
//  bool passDiElectronTrigger = false;
//  bool passEleMuonTrigger = false;
//  bool passDiMuonTrigger =false;

  bool passSingleElectronTrigger = false;
  bool passSingleMuonTrigger = false;
  bool passDoubleElectronTrigger = false;
  bool passDoubleMuonTrigger = false;
  bool passElectronMuonTrigger = false;
  bool passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1 = false;
  
  bool passHLT_IsoMu20_v = false;
  bool passHLT_IsoMu20_eta2p1_v = false;
  bool passHLT_IsoMu24_eta2p1_v = false;

  bool passHLT_Ele27_WP85_Gsf_v = false;
  bool passHLT_Ele27_eta2p1_WPLoose_Gsf_v = false;
  bool passHLT_Ele27_eta2p1_WP75_Gsf_v = false;

  bool passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v = false;
  bool passHLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v = false;

  bool passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = false;

  bool passHLT_Mu30_TkMu11_v = false;
  bool passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v = false;
  bool passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v = false;
  bool passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = false;

  bool passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = false;

  bool passHLT_Ele25WP60_SC4_Mass55_v = false;

  if( triggerResults.isValid() ){
    std::vector<std::string> triggerNames = hlt_config_.triggerNames();

    for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){
      std::string pathName = triggerNames[iPath];
      unsigned int hltIndex = hlt_config_.triggerIndex(pathName);

      if( hltIndex >= triggerResults->size() ) continue;

      int accept = triggerResults->accept(hltIndex);

      if( accept ){
      
	const unsigned long MatchedAtTheBegining = 0 ; 

	if( pathName.find( "HLT_Ele27_eta2p1_WPLoose_Gsf_v" ,0) == MatchedAtTheBegining ){    passSingleElectronTrigger = true;}
	if( pathName.find( "HLT_IsoMu20_v" ,0) == MatchedAtTheBegining ){    passSingleMuonTrigger = true; }
//	if( pathName.find( "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"       ,0) == MatchedAtTheBegining ){ passDiElectronTrigger = true;}
//	if( pathName.find( "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v" ,0) == MatchedAtTheBegining ){ passEleMuonTrigger = true;}
//	if( pathName.find( "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v"  ,0) == MatchedAtTheBegining ){ passEleMuonTrigger = true;}
//	if( pathName.find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"             ,0) == MatchedAtTheBegining ){ passDiMuonTrigger = true;}
//	if( pathName.find( "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"           ,0) == MatchedAtTheBegining ){ passDiMuonTrigger = true;}

	if( pathName.find( "HLT_IsoMu20_v"        ,0) == MatchedAtTheBegining ){ passHLT_IsoMu20_v = true;}
	if( pathName.find( "HLT_IsoMu20_eta2p1_v" ,0) == MatchedAtTheBegining ){ passHLT_IsoMu20_eta2p1_v = true;}
	if( pathName.find( "HLT_IsoMu24_eta2p1_v" ,0) == MatchedAtTheBegining ){ passHLT_IsoMu24_eta2p1_v = true;}

	if( pathName.find( "HLT_Ele27_WP85_Gsf_v"          ,0) == MatchedAtTheBegining ){ passHLT_Ele27_WP85_Gsf_v = true;}
	if( pathName.find( "HLT_Ele27_eta2p1_WPLoose_Gsf_v",0) == MatchedAtTheBegining ){ passHLT_Ele27_eta2p1_WPLoose_Gsf_v = true;}
	if( pathName.find( "HLT_Ele27_eta2p1_WP75_Gsf_v"   ,0) == MatchedAtTheBegining ){ passHLT_Ele27_eta2p1_WP75_Gsf_v = true;}

	if( pathName.find( "HLT_Ele27_eta2p1_WP85_Gsf_HT200_v"   ,0) == MatchedAtTheBegining ){ passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v = true;}
	if( pathName.find( "HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v",0) == MatchedAtTheBegining ){ passHLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v = true;}

	if( pathName.find( "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v" ,0) == MatchedAtTheBegining ){passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = true;}

	if( pathName.find( "HLT_Mu30_TkMu11_v" ,0) == MatchedAtTheBegining ){passHLT_Mu30_TkMu11_v = true;}
	if( pathName.find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"               ,0) == MatchedAtTheBegining ){ passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v = true;}
	if( pathName.find( "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"             ,0) == MatchedAtTheBegining ){ passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v = true;}
	if( pathName.find( "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",0) == MatchedAtTheBegining ){ passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = true;}

	if( pathName.find( "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v" ,0) == MatchedAtTheBegining ){ passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = true;}

	if( pathName.find( "HLT_Ele25WP60_SC4_Mass55_v"                       ,0) == MatchedAtTheBegining ){ passHLT_Ele25WP60_SC4_Mass55_v = true;}

	
      }
    }
  }


  edm::Handle<pat::JetCollection> pftempjets;
  iEvent.getByToken(tempjetToken,pftempjets);

  ////----------------------
  ////---- tt+X Categorization
  ////----------------------
  if( isMC ){
  edm::Handle<int> genTtbarId;
  iEvent.getByToken(genTtbarIdToken_, genTtbarId);
  // eve->additionalJetEventId_ = *genTtbarId%100;

  // Fill ID including information about b jets from top in acceptance
  h_ttbarId_->Fill(*genTtbarId);
    
  // Fill ID based only on additional b/c jets
  h_ttbarAdditionalJetId_->Fill(*genTtbarId%100);
  }

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

  int numpv=0; int iPV=0;
  bool firstGoodPV = false;
  reco::Vertex vertex;
  if( vtxHandle.isValid() ){
    for( reco::VertexCollection::const_iterator vtx = vtxs.begin(); vtx!=vtxs.end(); ++vtx ){
      
      iPV++;
      bool isGood = ( !(vtx->isFake()) &&
		      (vtx->ndof() >= 4.0) &&
		      (abs(vtx->z()) <= 24.0) &&
		      (abs(vtx->position().Rho()) <= 2.0) 
		      );
		      
      if( !isGood ) continue;

      if( iPV==1 ){
	firstGoodPV = true;
	vertex = (*vtx);
      }

      numpv++;
    }
  }

  eve->GoodFirstPV_=firstGoodPV;

  // if( !firstGoodPV ){
  //   std::cout << "ERROR!! First PV does not pass requirements! Skipping event" << std::endl;
  //   return;
  // }

  if( numpv>0 ) miniAODhelper.SetVertex(vertex);

  double numTruePV = -1;
  double numGenPV = -1;
  if( (PupInfo.isValid()) ){
    for( std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI ) {
      int BX = PVI->getBunchCrossing();
      if( BX==0 ){
	numTruePV = PVI->getTrueNumInteractions();
	numGenPV  = PVI->getPU_NumInteractions();
      }
    }
  }


  double wgt_lumi = 1.;
  if( insample_>=0 ) wgt_lumi = mySample_xSec_ * intLumi_ *1./ mySample_nGen_;

  // Increment event counter
  nevents++;
  nevents_wgt+=wgt_lumi;

  // Number of events after filtering, triggers, etc. (now nothing)
  nevents_clean++;
  nevents_clean_wgt+=wgt_lumi;

  int run  = iEvent.id().run();
  int lumi = iEvent.luminosityBlock();
  long evt = iEvent.id().event();


  const JetCorrector* corrector = JetCorrector::getJetCorrector( "ak4PFchsL1L2L3", iSetup );   //Get the jet corrector from the event setup

  miniAODhelper.SetJetCorrector(corrector);
  
  int mHdecay = -1;
  mHdecay = isMC ? miniAODhelper.GetHiggsDecay(mcparticles) : -1 ;
  eve->higgsDecayType_=mHdecay;

  eve->ttbarDecayType_ = isMC ? miniAODhelper.GetTTbarDecay(mcparticles) : -10 ;


  /////////
  ///
  /// Electrons
  ///
  ////////
  
  // miniAODhelper.SetElectronMVAinfo(h_conversioncollection, bsHandle);
  
  minTightLeptonPt = 30.0;
  looseLeptonPt = 15.0;
  
  //std::vector<pat::Electron> selectedElectrons_tight = miniAODhelper.GetSelectedElectrons( *electrons, minTightLeptonPt, electronID::electronEndOf15MVAmedium, 2.1 );
  std::vector<pat::Electron> selectedElectrons_tight = miniAODhelper.GetSelectedElectrons( electrons, minTightLeptonPt, electronID::electronEndOf15MVA80iso0p15, 2.1 );
  std::vector<pat::Electron> selectedElectrons_loose = miniAODhelper.GetSelectedElectrons( electrons, looseLeptonPt, electronID::electronEndOf15MVA80iso0p15, 2.4 );

  int numTightElectrons = int(selectedElectrons_tight.size());
  int numLooseElectrons = int(selectedElectrons_loose.size());// - numTightElectrons;


  /////////
  ///
  /// Muons
  ///
  ////////
  //std::vector<pat::Muon> selectedMuons_tight = miniAODhelper.GetSelectedMuons( *muons, 25, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.1 );
  std::vector<pat::Muon> selectedMuons_tight = miniAODhelper.GetSelectedMuons( *muons, 15, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.4);
  std::vector<pat::Muon> selectedMuons_loose = miniAODhelper.GetSelectedMuons( *muons, 15, muonID::muonTightDL, coneSize::R04, corrType::deltaBeta,2.4 );

  int numTightMuons = int(selectedMuons_tight.size());
  int numLooseMuons = int(selectedMuons_loose.size());// - numTightMuons;
  
  eve->numTightMuons_ = numTightMuons;
  eve->numLooseMuons_ = numLooseMuons;
  eve->numTightElectrons_ = numTightElectrons;
  eve->numLooseElectrons_ = numLooseElectrons;


  // Pass one tight lepton, zero loose lepton cuts
  bool passLJ = ( (numTightMuons+numTightElectrons)==1 &&
		  (numLooseMuons+numLooseElectrons)==1 );

  // Pass one tight lepton, exactly two leptons cuts
  bool passDIL = ( (numTightMuons+numTightElectrons)>=1 &&
		   (numLooseMuons+numLooseElectrons)==2 );


  // Event must either be LJ or DIL, as requested
  // remove for now
  //if( !((isLJ_ && passLJ) || (!isLJ_ && passDIL)) ) return;

  eve->PassLJ_ = ( passLJ ) ? 1 : 0;
  eve->PassDIL_ = ( passDIL ) ? 1 : 0;

  //// LJ subcategories
  int OneMuon = 0, OneElectron = 0;
  if( numTightMuons==1 && numLooseMuons==1 && numLooseElectrons == 0 )     OneMuon = 1;
  if( numTightElectrons==1 && numLooseElectrons==1 && numLooseMuons == 0 ) OneElectron = 1;

  eve->OneMuon_ = OneMuon;
  eve->OneElectron_ = OneElectron;

  //// DIL subcategories
  int TwoMuon = 0, TwoElectron = 0, MuonElectron = 0; 
  if ( numTightMuons>=1 && numLooseMuons==2 && numLooseElectrons==0 ) TwoMuon = 1;
  if ( numTightElectrons>=1 && numLooseElectrons==2 && numLooseMuons==0 ) TwoElectron = 1;
  if ( (numTightMuons + numTightElectrons)>=1 && numLooseMuons==1 && numLooseElectrons==1 ) MuonElectron = 1;

  eve->TwoMuon_ = TwoMuon;
  eve->TwoElectron_ = TwoElectron;
  eve->MuonElectron_ = MuonElectron;

   

  // Do jets stuff
  std::vector<pat::Jet> pfJets_ID = miniAODhelper.GetSelectedJets(*pfjets,0.,999,jetID::jetLoose,'-');
  std::vector<pat::Jet> pfJets_ID_clean = miniAODhelper.GetDeltaRCleanedJets( pfJets_ID, selectedMuons_loose, selectedElectrons_loose, 0.4);
  std::vector<pat::Jet> rawJets = miniAODhelper.GetUncorrectedJets(pfJets_ID_clean);
 // std::vector<pat::Jet> jetsNoMu = miniAODhelper.RemoveOverlaps(selectedMuons_loose, rawJets_ID);
 // std::vector<pat::Jet> jetsNoEle = miniAODhelper.RemoveOverlaps(selectedElectrons_loose, jetsNoMu);
  std::vector<pat::Jet> correctedJets_noSys = miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup);
  std::vector<pat::Jet> selectedJets_noSys_unsorted = miniAODhelper.GetSelectedJets(correctedJets_noSys, 30., 2.4, jetID::none, '-' );
  std::vector<pat::Jet> selectedJets_tag_noSys_unsorted = miniAODhelper.GetSelectedJets( correctedJets_noSys, 30., 2.4, jetID::none, 'M' );

  std::vector<pat::Jet> selectedJets_loose_noSys_unsorted = miniAODhelper.GetSelectedJets(correctedJets_noSys, 20., 3.0, jetID::none, '-' );
  std::vector<pat::Jet> selectedJets_loose_tag_noSys_unsorted = miniAODhelper.GetSelectedJets( correctedJets_noSys, 20., 3.0, jetID::none, 'M' );




  // Increment counter for passing cuts
  if( passLJ ){
    nevents_1l++; nevents_1l_wgt += wgt_lumi;
    if( OneMuon ){ nevents_1l_mu++; nevents_1l_mu_wgt += wgt_lumi; }
    else if( OneElectron ){ nevents_1l_ele++; nevents_1l_ele_wgt += wgt_lumi; }

    if( selectedJets_noSys_unsorted.size()>=4 ){
      nevents_1l_4j++; nevents_1l_4j_wgt += wgt_lumi;
      if( OneMuon ){ nevents_1l_4j_mu++; nevents_1l_4j_mu_wgt += wgt_lumi; }
      else if( OneElectron ){ nevents_1l_4j_ele++; nevents_1l_4j_ele_wgt += wgt_lumi; }

      if( selectedJets_tag_noSys_unsorted.size()>=2 ){
	nevents_1l_4j_2t++; nevents_1l_4j_2t_wgt += wgt_lumi;
	if( OneMuon ){ nevents_1l_4j_2t_mu++; nevents_1l_4j_2t_mu_wgt += wgt_lumi; }
	else if( OneElectron ){ nevents_1l_4j_2t_ele++; nevents_1l_4j_2t_ele_wgt += wgt_lumi; }
      }
    }
  }
  else if( passDIL ){
    nevents_2l++; nevents_2l_wgt += wgt_lumi;
    if( TwoMuon ){ nevents_2l_2m++; nevents_2l_2m_wgt += wgt_lumi; }
    else if( TwoElectron ){ nevents_2l_2e++; nevents_2l_2e_wgt += wgt_lumi; }
    else if( MuonElectron ){ nevents_2l_em++; nevents_2l_em_wgt += wgt_lumi; }

    if( selectedJets_noSys_unsorted.size()>=2 ){
      nevents_2l_2j++; nevents_2l_2j_wgt += wgt_lumi;
      if( TwoMuon ){ nevents_2l_2j_2m++; nevents_2l_2j_2m_wgt += wgt_lumi; }
      else if( TwoElectron ){ nevents_2l_2j_2e++; nevents_2l_2j_2e_wgt += wgt_lumi; }
      else if( MuonElectron ){ nevents_2l_2j_em++; nevents_2l_2j_em_wgt += wgt_lumi; }

      if( selectedJets_tag_noSys_unsorted.size()>=2 ){
	nevents_2l_2j_2t++; nevents_2l_2j_2t_wgt += wgt_lumi;
	if( TwoMuon ){ nevents_2l_2j_2t_2m++; nevents_2l_2j_2t_2m_wgt += wgt_lumi; }
	else if( TwoElectron ){ nevents_2l_2j_2t_2e++; nevents_2l_2j_2t_2e_wgt += wgt_lumi; }
	else if( MuonElectron ){ nevents_2l_2j_2t_em++; nevents_2l_2j_2t_em_wgt += wgt_lumi; }
      }
    }
  }


  eve->passElectronTrigger_ = 1;//( electron_triggerPass ) ? 1 : 0;
  eve->passMuonTrigger_     = 1;//( muon_triggerPass     ) ? 1 : 0;

  eve->passSingleElectronTrigger_ = ( passSingleElectronTrigger ) ? 1 : 0;
  eve->passSingleMuonTrigger_     = ( passSingleMuonTrigger ) ? 1 : 0;
  eve->passDoubleElectronTrigger_ = ( passDoubleElectronTrigger ) ? 1 : 0;
  eve->passDoubleMuonTrigger_     = ( passDoubleMuonTrigger ) ? 1 : 0;
  eve->passElectronMuonTrigger_   = ( passElectronMuonTrigger ) ? 1 : 0;
  eve->passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1_ = ( passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1 ) ? 1 : 0;
  
  eve->passHLT_IsoMu20_v_ =  ( passHLT_IsoMu20_v) ? 1 : 0;
  eve->passHLT_IsoMu20_eta2p1_v_ = ( passHLT_IsoMu20_eta2p1_v ) ? 1 : 0;
  eve->passHLT_IsoMu24_eta2p1_v_ = ( passHLT_IsoMu24_eta2p1_v ) ? 1 : 0;
  
  eve->passHLT_Ele27_WP85_Gsf_v_ = ( passHLT_Ele27_WP85_Gsf_v ) ? 1 : 0;
  eve->passHLT_Ele27_eta2p1_WPLoose_Gsf_v_ = ( passHLT_Ele27_eta2p1_WPLoose_Gsf_v ) ? 1 : 0;
  eve->passHLT_Ele27_eta2p1_WP75_Gsf_v_ = ( passHLT_Ele27_eta2p1_WP75_Gsf_v ) ? 1 : 0;
  
  eve->passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v_ = ( passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v ) ? 1 : 0;
  eve->passHLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v_ = ( passHLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v ) ? 1 : 0;
  
  eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = ( passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v ) ? 1 : 0;
  
  eve->passHLT_Mu30_TkMu11_v_ = ( passHLT_Mu30_TkMu11_v ) ? 1 : 0;
  eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v_ = ( passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v ) ? 1 : 0;
  eve->passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v_ = ( passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v ) ? 1 : 0;
  eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = ( passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v ) ? 1 : 0;
  
  eve->passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_ = ( passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v ) ? 1 : 0;
  
  eve->passHLT_Ele25WP60_SC4_Mass55_v_ = ( passHLT_Ele25WP60_SC4_Mass55_v ) ? 1 : 0;


  bool matchSingleMuTrigger = false;
  bool matchSingleElectronTrigger = false;

  // if( isMuon ) matchSingleMuTrigger = miniAODhelper.MuonMatchesSingleMuTrigger( selectedMuons_tight[0], hlt, hltobjs );
  // else         matchSingleElectronTrigger = miniAODhelper.ElectronMatchesSingleEleTrigger(selectedElectrons_tight[0], hlt, hltobjs );

  eve->matchSingleMuTrigger_ = ( matchSingleMuTrigger ) ? 1 : 0;
  eve->matchSingleElectronTrigger_ = ( matchSingleElectronTrigger ) ? 1 : 0;



  eve->run_ = run;
  eve->lumi_ = lumi;
  eve->evt_ = evt;



  eve->Q2ScaleUpWgt_ = 1.;
  eve->Q2ScaleDownWgt_ = 1.;
  // if( mySample.isTTJets ){
  // 	eve->Q2ScaleUpWgt_   = miniAODhelper.GetQ2ScaleUp(events.at(0));//1.402 * event->Q2ScaleUpWgt;
  // 	eve->Q2ScaleDownWgt_ = miniAODhelper.GetQ2ScaleDown(events.at(0));//0.683 * event->Q2ScaleDownWgt;
  // }

  eve->numTruePV_ = numTruePV;
  eve->numGenPV_ = numGenPV;

  eve->numPVs_ = numpv;
  eve->numSys_ = rNumSys;



  int lep_genId=-99;
  int lep_genParentId=-99;
  int lep_genGrandParentId=-99;

  TLorentzVector leptonV;
  double eps = 0.0001;
  leptonV.SetPxPyPzE(eps, eps, eps, eps);

  vdouble vlepton;
      
  double tight_lepton_relIso = -99;
  
  bool outputwords = false;
  if(evt==19226727 || evt==19035242 || evt==19353162 || evt==19681455 || evt==18985031 || evt==19014641)outputwords=false;
if(outputwords)cout<<OneMuon<<" "<<OneElectron<<endl;
  vvdouble vvleptons;
  
  //i want to save everything! set OneMuon and One Electron equal to true!
  
  
  if( OneMuon ) {
    leptonV.SetPxPyPzE( selectedMuons_tight.at(0).px(), selectedMuons_tight.at(0).py(), selectedMuons_tight.at(0).pz(), selectedMuons_tight.at(0).energy());
    vlepton.push_back(selectedMuons_tight.at(0).px());
    vlepton.push_back(selectedMuons_tight.at(0).py());
    vlepton.push_back(selectedMuons_tight.at(0).pz());
    vlepton.push_back(selectedMuons_tight.at(0).energy());

    if( (selectedMuons_tight.at(0).genLepton()) ){
      lep_genId = selectedMuons_tight.at(0).genLepton()->pdgId();

      if( selectedMuons_tight.at(0).genLepton()->numberOfMothers()>=1 ){
	lep_genParentId = selectedMuons_tight.at(0).genLepton()->mother(0)->pdgId();
	if( selectedMuons_tight.at(0).genLepton()->mother(0)->numberOfMothers()>=1 ) 
	  lep_genGrandParentId = selectedMuons_tight.at(0).genLepton()->mother(0)->mother(0)->pdgId();
      }
    }

    tight_lepton_relIso = miniAODhelper.GetMuonRelIso(selectedMuons_tight.at(0), coneSize::R04, corrType::deltaBeta);
  }
  else if( OneElectron){
  if(outputwords)cout<<"OneELECTRON!!"<<endl;
    leptonV.SetPxPyPzE( selectedElectrons_tight.at(0).px(), selectedElectrons_tight.at(0).py(), selectedElectrons_tight.at(0).pz(), selectedElectrons_tight.at(0).energy());
    vlepton.push_back(selectedElectrons_tight.at(0).px());
    vlepton.push_back(selectedElectrons_tight.at(0).py());
    vlepton.push_back(selectedElectrons_tight.at(0).pz());
    vlepton.push_back(selectedElectrons_tight.at(0).energy());

    if( (selectedElectrons_tight.at(0).genLepton()) ){
     if(outputwords)cout<<"why didnt this work? ";
      lep_genId = selectedElectrons_tight.at(0).genLepton()->pdgId();
if(outputwords)cout<<selectedElectrons_tight.at(0).genLepton()->pdgId();
      if( selectedElectrons_tight.at(0).genLepton()->numberOfMothers()>=1 ){
	lep_genParentId = selectedElectrons_tight.at(0).genLepton()->mother(0)->pdgId();
	if( selectedElectrons_tight.at(0).genLepton()->mother(0)->numberOfMothers()>=1 ) 
	  lep_genGrandParentId = selectedElectrons_tight.at(0).genLepton()->mother(0)->mother(0)->pdgId();
      }
    }

    tight_lepton_relIso = miniAODhelper.GetElectronRelIso(selectedElectrons_tight.at(0),coneSize::R03, corrType::rhoEA,effAreaType::spring15);
  }

  if( OneMuon || OneElectron ){
    eve->tight_lepton_pt_  = leptonV.Pt();
    eve->tight_lepton_eta_ = leptonV.Eta();
    eve->tight_lepton_phi_ = leptonV.Phi();
    eve->tight_lepton_relIso_ = tight_lepton_relIso;
    eve->tight_lepton_isMuon_ = ( OneMuon ) ? 1 : 0;

    eve->tight_lepton_genId_ = lep_genId;
    eve->tight_lepton_genParentId_ = lep_genParentId;
    eve->tight_lepton_genGrandParentId_ = lep_genGrandParentId;
  }


  vint lepton_genId, lepton_genParentId, lepton_genGrandParentId, lepton_trkCharge, lepton_isMuon, lepton_isTight, lepton_isLoose;
  vint lepton_isPhys14L, lepton_isPhys14M, lepton_isPhys14T;
  vdouble lepton_pt;
  vdouble lepton_eta;
  vdouble lepton_phi;
  vdouble lepton_relIso;
  vdouble lepton_relIsoR04;
  vdouble lepton_iso_sumChargedHadronPt;
  vdouble lepton_iso_sumNeutralHadronEt;
  vdouble lepton_iso_sumPhotonEt;
  vdouble lepton_iso_sumPUPt;
  vdouble lepton_mvaTrigValue;
  vint    lepton_numMissingHits;
  vdouble lepton_scSigmaIEtaIEta;
  vdouble lepton_full5x5_scSigmaIEtaIEta;
  vdouble lepton_hadronicOverEm;
  vdouble lepton_relEcalIso;
  vdouble lepton_relHcalIso;
  vdouble lepton_relTrackIso;
  vdouble lepton_OneOESuperMinusOneOP;
  vint    lepton_isEB;
  vint    lepton_passHLTId;
  vint    lepton_passConversionVeto;
  vint    lepton_inCrack;
  vdouble lepton_scEta;
  vdouble lepton_dEtaSCTrackAtVtx;
  vdouble lepton_dPhiSCTrackAtVtx;
  vdouble lepton_d0;
  vdouble lepton_dZ;
  vint lepton_isGlobalMuon;
  vint lepton_isTrackerMuon;
  vint lepton_isPFMuon;
  vdouble lepton_normalizedChi2;
  vint lepton_numberOfValidMuonHits;
  vint lepton_numberOfValidPixelHits;
  vint lepton_trackerLayersWithMeasurement;
  vint lepton_numberOfMatchedStations;

  std::vector<TLorentzVector> vec_TLV_lep;
  TLorentzVector sum_lepton_vect;
  sum_lepton_vect.SetPxPyPzE(0., 0., 0., 0.);

  // Loop over muons
  for( std::vector<pat::Muon>::const_iterator iMu = muons->begin(); iMu != muons->end(); iMu++ ){ 
 
    int genId=-99, genParentId=-99, genGrandParentId=-99;
    if( (iMu->genLepton()) ){
      genId = iMu->genLepton()->pdgId();
      if( iMu->genLepton()->numberOfMothers()>=1 ){
	genParentId = iMu->genLepton()->mother(0)->pdgId();
	if( iMu->genLepton()->mother(0)->numberOfMothers()>=1 ) genGrandParentId = iMu->genLepton()->mother(0)->mother(0)->pdgId();
      }
    }

    int trkCharge = -99;
    if( iMu->muonBestTrack().isAvailable() ) trkCharge = iMu->muonBestTrack()->charge();

    int isTight = ( miniAODhelper.isGoodMuon(*iMu, minTightLeptonPt, 2.1, muonID::muonTight, coneSize::R04, corrType::deltaBeta) ) ? 1 : 0;
    int isLoose = ( miniAODhelper.isGoodMuon(*iMu, looseLeptonPt, 2.4, muonID::muonTightDL, coneSize::R04, corrType::deltaBeta) ) ? 1 : 0;

    int isPhys14L = false;
    int isPhys14M = false;
    int isPhys14T = false;

    double d0 = -999;
    double dZ = -999;
    if( iMu->muonBestTrack().isAvailable() ){
      d0 = iMu->muonBestTrack()->dxy(vertex.position());
      dZ = iMu->muonBestTrack()->dz(vertex.position());
    }

    double normalizedChi2 = -999;
    int numberOfValidMuonHits = -999;
    if( iMu->globalTrack().isAvailable() ){
      normalizedChi2 = iMu->globalTrack()->normalizedChi2();
      numberOfValidMuonHits = iMu->globalTrack()->hitPattern().numberOfValidMuonHits();
    }

    int numberOfValidPixelHits = -999;
    if( iMu->innerTrack().isAvailable() ){
      numberOfValidPixelHits = iMu->innerTrack()->hitPattern().numberOfValidPixelHits();
    }

    int trackerLayersWithMeasurement = -999;
    if( iMu->track().isAvailable() ){
      trackerLayersWithMeasurement = iMu->track()->hitPattern().trackerLayersWithMeasurement();
    }

    int numberOfMatchedStations = iMu->numberOfMatchedStations();


    lepton_trkCharge.push_back(trkCharge);
    lepton_isMuon.push_back(1);
    lepton_isTight.push_back(isTight);
    lepton_isLoose.push_back(isLoose);
    lepton_isPhys14L.push_back(isPhys14L);
    lepton_isPhys14M.push_back(isPhys14M);
    lepton_isPhys14T.push_back(isPhys14T);
    lepton_genId.push_back(genId);
    lepton_genParentId.push_back(genParentId);
    lepton_genGrandParentId.push_back(genGrandParentId);
    lepton_pt.push_back(iMu->pt());
    lepton_eta.push_back(iMu->eta());
    lepton_phi.push_back(iMu->phi());
    lepton_relIso.push_back(miniAODhelper.GetMuonRelIso(*iMu));
    lepton_relIsoR04.push_back(miniAODhelper.GetMuonRelIso(*iMu, coneSize::R04, corrType::deltaBeta));
    lepton_iso_sumChargedHadronPt.push_back(iMu->pfIsolationR03().sumChargedHadronPt);
    lepton_iso_sumNeutralHadronEt.push_back(iMu->pfIsolationR03().sumNeutralHadronEt);
    lepton_iso_sumPhotonEt.push_back(iMu->pfIsolationR03().sumPhotonEt);
    lepton_iso_sumPUPt.push_back(iMu->pfIsolationR03().sumPUPt);
    lepton_mvaTrigValue.push_back(-99);
    lepton_scSigmaIEtaIEta.push_back(-99);
    lepton_full5x5_scSigmaIEtaIEta.push_back(-99);
    lepton_hadronicOverEm.push_back(-99);
    lepton_relEcalIso.push_back(-99);
    lepton_relHcalIso.push_back(-99);
    lepton_relTrackIso.push_back(-99);
    lepton_OneOESuperMinusOneOP.push_back(-99);
    lepton_numMissingHits.push_back(-99);
    lepton_isEB.push_back(-99);
    lepton_passHLTId.push_back(-99);
    lepton_passConversionVeto.push_back(-99);
    lepton_inCrack.push_back(-99);
    lepton_scEta.push_back(-99);
    lepton_dEtaSCTrackAtVtx.push_back(-99);
    lepton_dPhiSCTrackAtVtx.push_back(-99);

    lepton_d0.push_back(d0);
    lepton_dZ.push_back(dZ);
    lepton_isGlobalMuon.push_back(iMu->isGlobalMuon());
    lepton_isTrackerMuon.push_back(iMu->isTrackerMuon());
    lepton_isPFMuon.push_back(iMu->isPFMuon());
    lepton_normalizedChi2.push_back(normalizedChi2);
    lepton_numberOfValidMuonHits.push_back(numberOfValidMuonHits);
    lepton_numberOfValidPixelHits.push_back(numberOfValidPixelHits);
    lepton_trackerLayersWithMeasurement.push_back(trackerLayersWithMeasurement);
    lepton_numberOfMatchedStations.push_back(numberOfMatchedStations);


    // Get muon 4Vector and add to vecTLorentzVector for muons
    TLorentzVector leptonP4;	  
    leptonP4.SetPxPyPzE(iMu->px(),iMu->py(),iMu->pz(),iMu->energy());
    vec_TLV_lep.push_back(leptonP4);

    sum_lepton_vect += leptonP4;

    // make vvdouble version of vecTLorentzVector
    vdouble vleptons;
    vleptons.push_back(iMu->px());
    vleptons.push_back(iMu->py());
    vleptons.push_back(iMu->pz());
    vleptons.push_back(iMu->energy());
    vvleptons.push_back(vleptons);
  }

  // Loop over electrons
  for( std::vector<pat::Electron>::const_iterator iEle = electrons.begin(); iEle != electrons.end(); iEle++ ){ 

    int genId=-99, genParentId=-99, genGrandParentId=-99;
    if( (iEle->genLepton()) ){
      genId = iEle->genLepton()->pdgId();
      if( iEle->genLepton()->numberOfMothers()>=1 ){
	genParentId = iEle->genLepton()->mother(0)->pdgId();
	if( iEle->genLepton()->mother(0)->numberOfMothers()>=1 ) genGrandParentId = iEle->genLepton()->mother(0)->mother(0)->pdgId();
      }
    }

    int trkCharge = -99;
    if( iEle->gsfTrack().isAvailable() ) trkCharge = iEle->gsfTrack()->charge();

    int isTight = ( miniAODhelper.isGoodElectron(*iEle, minTightLeptonPt, 2.1, electronID::electronTight) ) ? 1 : 0;
    int isLoose = ( miniAODhelper.isGoodElectron(*iEle, looseLeptonPt, 2.4, electronID::electronLoose) ) ? 1 : 0;

    int isPhys14L = ( miniAODhelper.isGoodElectron(*iEle, looseLeptonPt, 2.4, electronID::electronPhys14L) ) ? 1 : 0;
    int isPhys14M = ( miniAODhelper.isGoodElectron(*iEle, looseLeptonPt, 2.4, electronID::electronPhys14M) ) ? 1 : 0;
    int isPhys14T = ( miniAODhelper.isGoodElectron(*iEle, looseLeptonPt, 2.4, electronID::electronPhys14T) ) ? 1 : 0;


    double mvaTrigValue = -99;//myMVATrig->mvaValue(*iEle,false);

    bool inCrack = false;
    double scEta = -99;
    if( iEle->superCluster().isAvailable() ){
      inCrack = ( fabs(iEle->superCluster()->position().eta())>1.4442 && fabs(iEle->superCluster()->position().eta())<1.5660 );
      scEta = iEle->superCluster()->position().eta();
    }

    double ooEmooP = -999;
    if( iEle->ecalEnergy() == 0 ) ooEmooP = 1e30;
    else if( !std::isfinite(iEle->ecalEnergy()) ) ooEmooP = 1e30;
    else ooEmooP = fabs(1.0/iEle->ecalEnergy() - iEle->eSuperClusterOverP()/iEle->ecalEnergy() );

    bool passHLTId = false;
    double OneOESuperMinusOneOP = ooEmooP;
    int numMissingHits = 99;
    int isEB = -1;
    double relEcalIso = 99, relHcalIso = 99, relTrackIso = 99;
    double d0 = -999;
    double dZ = -999;
    if( (iEle->superCluster().isAvailable()) && (iEle->gsfTrack().isAvailable()) ){
      double SCenergy = iEle->superCluster()->energy();
      double SCeta = iEle->superCluster()->position().eta();
      double absSCeta = fabs( SCeta );
      double SCet = SCenergy * sin (2*atan(exp(-SCeta))); 

      numMissingHits = iEle->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

      d0 = fabs(iEle->gsfTrack()->dxy(vertex.position()));
      dZ = fabs(iEle->gsfTrack()->dz(vertex.position()));

      relEcalIso = iEle->ecalIso()/SCet;
      relHcalIso = iEle->hcalIso()/SCet;
      relTrackIso = iEle->trackIso()/SCet;

      if( absSCeta < 1.479 ){
	isEB = 1;
	passHLTId = ( iEle->scSigmaIEtaIEta() <= 0.011 &&
		      iEle->hadronicOverEm() <= 0.15 &&
		      iEle->ecalIso()/SCet <= 0.16 &&
		      iEle->hcalIso()/SCet <= 0.20 &&
		      iEle->trackIso()/SCet <= 0.05 &&
		      OneOESuperMinusOneOP <= 0.012 &&
		      numMissingHits <= 999
		      );
      }
      else{
	isEB = 0;
	passHLTId = ( iEle->scSigmaIEtaIEta() <= 0.033 &&
		      iEle->hadronicOverEm() <= 0.20 &&
		      iEle->ecalIso()/SCet <= 0.12 &&
		      iEle->hcalIso()/SCet <= 0.30 &&
		      iEle->trackIso()/SCet <= 0.05 &&
		      OneOESuperMinusOneOP <= 0.009 &&
		      numMissingHits <= 1
		      );
      }
    }

    lepton_trkCharge.push_back(trkCharge);
    lepton_isMuon.push_back(0);
    lepton_isTight.push_back(isTight);
    lepton_isLoose.push_back(isLoose);
    lepton_isPhys14L.push_back(isPhys14L);
    lepton_isPhys14M.push_back(isPhys14M);
    lepton_isPhys14T.push_back(isPhys14T);
    lepton_genId.push_back(genId);
    lepton_genParentId.push_back(genParentId);
    lepton_genGrandParentId.push_back(genGrandParentId);
    lepton_pt.push_back(iEle->pt());
    lepton_eta.push_back(iEle->eta());
    lepton_phi.push_back(iEle->phi());
    lepton_relIso.push_back(miniAODhelper.GetElectronRelIso(*iEle));
    lepton_relIsoR04.push_back(miniAODhelper.GetElectronRelIso(*iEle));
    lepton_iso_sumChargedHadronPt.push_back(iEle->pfIsolationVariables().sumChargedHadronPt);
    lepton_iso_sumNeutralHadronEt.push_back(iEle->pfIsolationVariables().sumNeutralHadronEt);
    lepton_iso_sumPhotonEt.push_back(iEle->pfIsolationVariables().sumPhotonEt);
    lepton_iso_sumPUPt.push_back(iEle->pfIsolationVariables().sumPUPt);
    lepton_mvaTrigValue.push_back(mvaTrigValue);
    lepton_scSigmaIEtaIEta.push_back(iEle->scSigmaIEtaIEta());
    lepton_full5x5_scSigmaIEtaIEta.push_back(iEle->full5x5_sigmaIetaIeta());
    lepton_hadronicOverEm.push_back(iEle->hadronicOverEm());
    lepton_relEcalIso.push_back(relEcalIso);
    lepton_relHcalIso.push_back(relHcalIso);
    lepton_relTrackIso.push_back(relTrackIso);
    lepton_OneOESuperMinusOneOP.push_back(OneOESuperMinusOneOP);
    lepton_numMissingHits.push_back(numMissingHits);
    lepton_isEB.push_back(isEB);
    lepton_passHLTId.push_back(passHLTId);
    lepton_passConversionVeto.push_back(iEle->passConversionVeto());
    lepton_inCrack.push_back(inCrack);
    lepton_scEta.push_back(scEta);
    lepton_dEtaSCTrackAtVtx.push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
    lepton_dPhiSCTrackAtVtx.push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
    lepton_d0.push_back(d0);
    lepton_dZ.push_back(dZ);
    lepton_isGlobalMuon.push_back(-99);
    lepton_isTrackerMuon.push_back(-99);
    lepton_isPFMuon.push_back(-99);
    lepton_normalizedChi2.push_back(-99);
    lepton_numberOfValidMuonHits.push_back(-99);
    lepton_numberOfValidPixelHits.push_back(-99);
    lepton_trackerLayersWithMeasurement.push_back(-99);
    lepton_numberOfMatchedStations.push_back(-99);



    // Get electron 4Vector and add to vecTLorentzVector for electrons
    TLorentzVector leptonP4;	  
    leptonP4.SetPxPyPzE(iEle->px(),iEle->py(),iEle->pz(),iEle->energy());
    vec_TLV_lep.push_back(leptonP4);

    sum_lepton_vect += leptonP4;

    // make vvdouble version of vecTLorentzVector
    vdouble vleptons;
    vleptons.push_back(iEle->px());
    vleptons.push_back(iEle->py());
    vleptons.push_back(iEle->pz());
    vleptons.push_back(iEle->energy());
    vvleptons.push_back(vleptons);
  }

  eve->lepton_vect_TLV_         = vvleptons;
  eve->lepton_trkCharge_        = lepton_trkCharge;
  eve->lepton_isMuon_           = lepton_isMuon;
  eve->lepton_isTight_          = lepton_isTight;
  eve->lepton_isLoose_          = lepton_isLoose;
  eve->lepton_isPhys14L_        = lepton_isPhys14L;
  eve->lepton_isPhys14M_        = lepton_isPhys14M;
  eve->lepton_isPhys14T_        = lepton_isPhys14T;
  eve->lepton_genId_            = lepton_genId;
  eve->lepton_genParentId_      = lepton_genParentId;
  eve->lepton_genGrandParentId_ = lepton_genGrandParentId;
  eve->lepton_pt_               = lepton_pt;
  eve->lepton_eta_              = lepton_eta;
  eve->lepton_phi_              = lepton_phi;
  eve->lepton_relIso_           = lepton_relIso;
  eve->lepton_relIsoR04_           = lepton_relIsoR04;
  eve->lepton_iso_sumChargedHadronPt_ = lepton_iso_sumChargedHadronPt;
  eve->lepton_iso_sumNeutralHadronEt_ = lepton_iso_sumNeutralHadronEt;
  eve->lepton_iso_sumPhotonEt_        = lepton_iso_sumPhotonEt;
  eve->lepton_iso_sumPUPt_            = lepton_iso_sumPUPt;
  eve->lepton_mvaTrigValue_     = lepton_mvaTrigValue;
  eve->lepton_scSigmaIEtaIEta_  = lepton_scSigmaIEtaIEta;
  eve->lepton_full5x5_scSigmaIEtaIEta_ = lepton_full5x5_scSigmaIEtaIEta;
  eve->lepton_hadronicOverEm_   = lepton_hadronicOverEm;
  eve->lepton_relEcalIso_       = lepton_relEcalIso;
  eve->lepton_relHcalIso_       = lepton_relHcalIso;
  eve->lepton_relTrackIso_      = lepton_relTrackIso;
  eve->lepton_OneOESuperMinusOneOP_ = lepton_OneOESuperMinusOneOP;
  eve->lepton_numMissingHits_   = lepton_numMissingHits;
  eve->lepton_isEB_             = lepton_isEB;
  eve->lepton_passHLTId_        = lepton_passHLTId;
  eve->lepton_passConversionVeto_ = lepton_passConversionVeto;
  eve->lepton_inCrack_          = lepton_inCrack;
  eve->lepton_scEta_            = lepton_scEta;
  eve->lepton_dEtaSCTrackAtVtx_ = lepton_dEtaSCTrackAtVtx;
  eve->lepton_dPhiSCTrackAtVtx_ = lepton_dPhiSCTrackAtVtx;
  eve->lepton_d0_ = lepton_d0;
  eve->lepton_dZ_ = lepton_dZ;
  eve->lepton_isGlobalMuon_ = lepton_isGlobalMuon;
  eve->lepton_isTrackerMuon_ = lepton_isTrackerMuon;
  eve->lepton_isPFMuon_ = lepton_isPFMuon;
  eve->lepton_normalizedChi2_ = lepton_normalizedChi2;
  eve->lepton_numberOfValidMuonHits_ = lepton_numberOfValidMuonHits;
  eve->lepton_numberOfValidPixelHits_ = lepton_numberOfValidPixelHits;
  eve->lepton_trackerLayersWithMeasurement_ = lepton_trackerLayersWithMeasurement;
  eve->lepton_numberOfMatchedStations_ = lepton_numberOfMatchedStations;


  /// DIL specific, doesn't make sense in current scope
  int oppositeLepCharge = -9;
  if( lepton_trkCharge.size()==2 ){
    int chg0 = lepton_trkCharge[0];
    int chg1 = lepton_trkCharge[1];

    if( (chg0 * chg1)==-1 )     oppositeLepCharge = 1;
    else if( (chg0 * chg1)==1 ) oppositeLepCharge = 0;
    else if( chg0==-99 )        oppositeLepCharge = -1;
    else if( chg1==-99 )        oppositeLepCharge = -2;
    else                        oppositeLepCharge = -3;
  }
  eve->oppositeLepCharge_ = oppositeLepCharge;


  /// DIL specific, doesn't make sense in current scope
  double mass_leplep = -99;
  if( vec_TLV_lep.size()==2 ){
    mass_leplep = sum_lepton_vect.M();
    eve->mass_leplep_ = mass_leplep;
    eve->dR_leplep_ = vec_TLV_lep[0].DeltaR(vec_TLV_lep[1]);
  }

  // PU scale factor
  double wgtPU = 1, wgtPUup = 1, wgtPUdown = 1;

  // wgtPU = miniAODhelper.GetPUweight( numTruePV );
  // wgtPUup = miniAODhelper.GetPUweightUp( numTruePV );
  // wgtPUdown = miniAODhelper.GetPUweightDown( numTruePV );

  eve->wgt_pu_     = wgtPU;
  eve->wgt_puUp_   = wgtPUup;
  eve->wgt_puDown_ = wgtPUdown;


  // Lepton scale factor
  float leptonSF = 1.0;
  // if( insample_>=0 ){
  // 	if( isMuon ) leptonSF = miniAODhelper.GetMuonSF(selectedMuons_tight.at(0));
  // 	else         leptonSF = miniAODhelper.GetElectronSF(selectedElectrons_tight.at(0));
  // }

  eve->wgt_lepSF_  = leptonSF;

  eve->wgt_lumi_  = intLumi_;
  eve->wgt_xs_    = mySample_xSec_;//mySample.xSec;
  eve->wgt_nGen_  = mySample_nGen_;//mySample.nGen;

  eve->wgt_generator_ = GenEventInfoWeight;

  double wgt_topPtSF     = 1;//miniAODhelper.GetTopPtweight( mcparticles );
  double wgt_topPtSFUp   = 1;//miniAODhelper.GetTopPtweightUp( mcparticles );
  double wgt_topPtSFDown = 1;//miniAODhelper.GetTopPtweightDown( mcparticles );

  eve->wgt_topPt_ = wgt_topPtSF;
  eve->wgt_topPtUp_ = wgt_topPtSFUp;
  eve->wgt_topPtDown_ = wgt_topPtSFDown;



  double top_pt = -1;
  double antitop_pt = -1;

  // Get top pTs for reweighting downstream
  // if( mySample.isTTJets ){
  // 	for( int iPar=0; iPar<int(mcparticles.size()); iPar++ ){
  // 	  int id = mcparticles.at(iPar).id;
  // 	  if( id==6 )       top_pt     = mcparticles.at(iPar).pt;
  // 	  else if( id==-6 ) antitop_pt = mcparticles.at(iPar).pt;
  // 	}
  // }
  eve->top_pt_ = top_pt;
  eve->antitop_pt_ = antitop_pt;


  ///////////
  ///// Boosted jet information
  ///////////

  // ///// HEP top tagged jet
  // int numTopTags = 0;
  // vvdouble vvjets_topfatJet;
  // vvdouble vvjets_nonW;
  // vvdouble vvjets_W1;
  // vvdouble vvjets_W2;

  // std::vector<bool> isTopTag;

  // if( h_heptopjet.isValid() ){
  //   boosted::HEPTopJetCollection const &heptopjets_unsorted = *h_heptopjet;
  //   boosted::HEPTopJetCollection heptopjets = BoostedUtils::GetSortedByPt(heptopjets_unsorted);

  //   for( boosted::HEPTopJetCollection::iterator topJet = heptopjets.begin() ; topJet != heptopjets.end(); ++topJet ){

  //     bool toptag = BoostedUtils::GetTopTag(*topJet);
  //     isTopTag.push_back(toptag);

  //     pat::Jet ifatJet = topJet->fatjet;
  //     vdouble vjets_fatJet;
  //     vjets_fatJet.push_back(ifatJet.px());
  //     vjets_fatJet.push_back(ifatJet.py());
  //     vjets_fatJet.push_back(ifatJet.pz());
  //     vjets_fatJet.push_back(ifatJet.energy());
  //     vvjets_topfatJet.push_back(vjets_fatJet);

  //     pat::Jet inonW = topJet->nonW;
  //     vdouble vjets_nonW;
  //     vjets_nonW.push_back(inonW.px());
  //     vjets_nonW.push_back(inonW.py());
  //     vjets_nonW.push_back(inonW.pz());
  //     vjets_nonW.push_back(inonW.energy());
  //     vvjets_nonW.push_back(vjets_nonW);

  //     pat::Jet iW1 = topJet->W1;
  //     vdouble vjets_W1;
  //     vjets_W1.push_back(iW1.px());
  //     vjets_W1.push_back(iW1.py());
  //     vjets_W1.push_back(iW1.pz());
  //     vjets_W1.push_back(iW1.energy());
  //     vvjets_W1.push_back(vjets_W1);

  //     pat::Jet iW2 = topJet->W2;
  //     vdouble vjets_W2;
  //     vjets_W2.push_back(iW2.px());
  //     vjets_W2.push_back(iW2.py());
  //     vjets_W2.push_back(iW2.pz());
  //     vjets_W2.push_back(iW2.energy());
  //     vvjets_W2.push_back(vjets_W2);

  //     // pt and eta requirements on top jet
  //     if( !(topJet->fatjet.pt() > 250. && abs(topJet->fatjet.eta()) < 1.8) ) continue;

  //     // pt and eta requirements on subjets
  //     if( !( (topJet->nonW.pt()>20 && abs(topJet->nonW.eta())<2.5 ) &&
  // 	     (topJet->W1.pt()>20 && abs(topJet->W1.eta())<2.5 ) &&
  // 	     (topJet->W2.pt()>20 && abs(topJet->W2.eta())<2.5 ) ) ) continue;

  //     // must be top-tagged
  //     if( !toptag ) continue;

  //     numTopTags++;
  //   }

  // }

  // eve->topfatJet_vect_TLV_ = vvjets_topfatJet;
  // eve->nonW_vect_TLV_ = vvjets_nonW;
  // eve->W1_vect_TLV_ = vvjets_W1;
  // eve->W2_vect_TLV_ = vvjets_W2;
  // eve->numTopTags_ = numTopTags;

  // ///// Higgs tagged jet
  // int numHiggsTags = 0;
  // vvdouble vvjets_higgsfatJet;
  // std::vector<vvdouble> vvjets_higgsfilterjet_all;
  // vvdouble csv_filterjet_all;
  // if( h_subfilterjet.isValid() ){
  //   boosted::SubFilterJetCollection const &subfilterjets_unsorted = *h_subfilterjet;
  //   boosted::SubFilterJetCollection subfilterjets = BoostedUtils::GetSortedByPt(subfilterjets_unsorted);

  //   for( boosted::SubFilterJetCollection::iterator higgsJet = subfilterjets.begin() ; higgsJet != subfilterjets.end(); ++higgsJet ){

  //     pat::Jet ifatJet = higgsJet->fatjet;
  //     vdouble vjets_fatJet;
  //     vjets_fatJet.push_back(ifatJet.px());
  //     vjets_fatJet.push_back(ifatJet.py());
  //     vjets_fatJet.push_back(ifatJet.pz());
  //     vjets_fatJet.push_back(ifatJet.energy());
  //     vvjets_higgsfatJet.push_back(vjets_fatJet);

  //     vdouble csv_filterjet;
  //     vvdouble vvjets_higgsfilterjet;
  //     std::vector<pat::Jet> filterjets_tmp = higgsJet->filterjets;
  //     int numFiltJets_tmp = filterjets_tmp.size();
  //     for( int ijet=0; ijet<numFiltJets_tmp; ijet++ ){
  // 	pat::Jet ifilterjet = filterjets_tmp[ijet];
  // 	double csv_ifilterjet = ifilterjet.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags");

  // 	vdouble vjets_filterjet;
  // 	vjets_filterjet.push_back(ifilterjet.px());
  // 	vjets_filterjet.push_back(ifilterjet.py());
  // 	vjets_filterjet.push_back(ifilterjet.pz());
  // 	vjets_filterjet.push_back(ifilterjet.energy());
  // 	vvjets_higgsfilterjet.push_back(vjets_filterjet);
	
  // 	csv_filterjet.push_back(csv_ifilterjet);

  //     }

  //     vvjets_higgsfilterjet_all.push_back(vvjets_higgsfilterjet);
  //     csv_filterjet_all.push_back(csv_filterjet);

  //     // pt and eta requirements on top jet
  //     if( !(higgsJet->fatjet.pt() > 250. && abs(higgsJet->fatjet.eta()) < 1.8) ) continue;

  //     int numBtagFiltJets=0;
  //     std::vector<pat::Jet> filterjets = higgsJet->filterjets;
  //     int numFiltJets = filterjets.size();
  //     for( int ijet=0; ijet<numFiltJets; ijet++ ){
  // 	if( verbose_ ){
  // 	  printf("\t\t filt jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = %.2f,\t CSVv2 = %+5.3f,\t CSVv1 = %+5.3f \n",
  // 		 ijet, filterjets[ijet].pt(), filterjets[ijet].eta(), filterjets[ijet].phi(), 
  // 		 filterjets[ijet].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"),
  // 		 filterjets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"));
  // 	}

  // 	if( !(filterjets[ijet].pt()>20. && abs(filterjets[ijet].eta()) < 2.5) ) continue;
  // 	if( !(filterjets[ijet].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags") > 0.814) ) continue;
  // 	numBtagFiltJets++;
  //     }

  //     if( verbose_ ){
  // 	printf("\t Higgs jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = %.2f,\t numFiltJets = %2d,\t numBtagFiltJets = %2d\n",
  // 	       int(higgsJet - subfilterjets.begin()), higgsJet->fatjet.pt(), higgsJet->fatjet.eta(), higgsJet->fatjet.phi(), numFiltJets, numBtagFiltJets );
  //     }

  //     if( numBtagFiltJets>=2 ) numHiggsTags++;
  //   }
  // }

  // eve->numHiggsTags_ = numHiggsTags;
  // eve->higgsfatJet_vect_TLV_ = vvjets_higgsfatJet;
  // eve->higgsfilterjet_all_vect_TLV_ = vvjets_higgsfilterjet_all;
  // eve->csv_filterjet_all_ = csv_filterjet_all;

  //-----------

  // Special collection used for correcting MET collection
  // BNjetCollection pfJets_forMET  = miniAODhelper.GetSelectedJets( pfjets, 12., 4.9, jetID::jetLoose, '-' ); 
  
 

  // bool for passing jet/tag cuts, first is for tree filling, second is for event counting
  bool hasNumJetNumTag = false;


  double totalWgt = wgt_lumi;

  // Loop over systematics
  for( int iSys=0; iSys<rNumSys; iSys++ ){

    eve->wgt_[iSys] = wgt_lumi;

    // no systematic up/down for lepton SF for now
    if( iSys==1 )      eve->wgt_[iSys] *= leptonSF;
    else if( iSys==2 ) eve->wgt_[iSys] *= leptonSF;
    else               eve->wgt_[iSys] *= leptonSF;


    // PU weight up/down
    if( iSys==3 )      eve->wgt_[iSys] *= wgtPUup;
    else if( iSys==4 ) eve->wgt_[iSys] *= wgtPUdown;
    else               eve->wgt_[iSys] *= wgtPU;

    if( iSys==0 ) totalWgt *= leptonSF*wgtPU;

    sysType::sysType iSysType = sysType::NA;
    switch( iSys ){
    case 5 : iSysType = sysType::JERup;    break;
    case 6 : iSysType = sysType::JERdown;  break;
    case 7 : iSysType = sysType::JESup;    break;
    case 8 : iSysType = sysType::JESdown;  break;
    case 9 : iSysType = sysType::CSVLFup;         break;
    case 10: iSysType = sysType::CSVLFdown;       break;
    case 11: iSysType = sysType::CSVHFup;         break;
    case 12: iSysType = sysType::CSVHFdown;       break;
    case 13: iSysType = sysType::CSVHFStats1up;   break;
    case 14: iSysType = sysType::CSVHFStats1down; break;
    case 15: iSysType = sysType::CSVHFStats2up;   break;
    case 16: iSysType = sysType::CSVHFStats2down; break;
    case 17: iSysType = sysType::CSVLFStats1up;   break;
    case 18: iSysType = sysType::CSVLFStats1down; break;
    case 19: iSysType = sysType::CSVLFStats2up;   break;
    case 20: iSysType = sysType::CSVLFStats2down; break;
    default: iSysType = sysType::NA;       break;
    }


    bool IsTauTauLeptonEvent = false;//miniAODhelper.IsTauEvent(taus, pfjets, electrons, muons, iSysType);
    eve->IsTauTauLeptonEvent_[iSys] = ( IsTauTauLeptonEvent ) ? 1 : 0;

    // Add Q2scale systematic
    if( iSys==21 )      eve->wgt_[iSys] *= eve->Q2ScaleUpWgt_;
    else if( iSys==22 ) eve->wgt_[iSys] *= eve->Q2ScaleDownWgt_;


    // top pT reweight
    if( iSys==23 )      eve->wgt_[iSys] *= wgt_topPtSFUp;
    else if( iSys==24 ) eve->wgt_[iSys] *= wgt_topPtSFDown;
    else                eve->wgt_[iSys] *= wgt_topPtSF;

    if( iSys==0 ) totalWgt *= wgt_topPtSF;


    // nJets/Tags
    // THESE MUST BE INITIALIZED TO -1
    eve->numJets_float_[iSys] = -1;
    eve->numTags_float_[iSys] = -1;


    /////////
    ///
    /// Pfjets
    ///
    ////////

    std::vector<pat::Jet> correctedJets = ( !(iSys>=5 && iSys<=8) ) ? correctedJets_noSys : miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, iSysType);
    std::vector<pat::Jet> selectedJets_unsorted = ( !(iSys>=5 && iSys<=8) ) ? selectedJets_noSys_unsorted : miniAODhelper.GetSelectedJets(correctedJets, 30., 2.4, jetID::none, '-' );


    // Get CSVM tagged jet collection
    std::vector<pat::Jet> selectedJets_tag_unsorted = ( !(iSys>=5 && iSys<=8) ) ? selectedJets_tag_noSys_unsorted : miniAODhelper.GetSelectedJets( correctedJets, 30., 2.4, jetID::none, 'M' );

    // Get nontagged jet collection
    std::vector<pat::Jet> selectedJets_untag_unsorted = selectedJets_tag_unsorted;//miniAODhelper.GetSelectedJets( correctedJets, 30., 2.4, jetID::none, 'M' ); 
    //std::vector<pat::Jet> selectedJets_untag_unsorted = miniAODhelper.GetSymmetricDifference( selectedJets_unsorted, selectedJets_tag_unsorted );



    // // Require at least 2 tags
    // if( !(selectedJets_tag_unsorted.size()>=2) ) continue;

  
    // Sort jet collections by pT
    std::vector<pat::Jet> selectedJets_uncleaned       = miniAODhelper.GetSortedByPt( selectedJets_unsorted );
    std::vector<pat::Jet> selectedJets_tag_uncleaned   = miniAODhelper.GetSortedByPt( selectedJets_tag_unsorted );
    std::vector<pat::Jet> selectedJets_untag_uncleaned = miniAODhelper.GetSortedByPt( selectedJets_untag_unsorted );
    
    std::vector<pat::Jet> selectedJets        = miniAODhelper.GetDeltaRCleanedJets( selectedJets_uncleaned,selectedMuons_loose,selectedElectrons_loose,0.4);
    std::vector<pat::Jet> selectedJets_tag    = miniAODhelper.GetDeltaRCleanedJets( selectedJets_tag_uncleaned,selectedMuons_loose,selectedElectrons_loose,0.4);
    std::vector<pat::Jet> selectedJets_untag = miniAODhelper.GetDeltaRCleanedJets( selectedJets_untag_uncleaned,selectedMuons_loose,selectedElectrons_loose,0.4);

    //if( mySample.isTTJets ) splitEvent = miniAODhelper.ttPlusHFKeepEvent( mcparticles, selectedJets );
    //if( !splitEvent ) continue;

    // if( mySample.isTTJets ){
    //   eve->ttbb_algo_result_[iSys] = miniAODhelper.ttPlusBBClassifyEvent( mcparticles, selectedJets );
    //   eve->ttcc_algo_result_[iSys] = miniAODhelper.ttPlusCCClassifyEvent( mcparticles, selectedJets );
    // }


    // Get numJets, numTags
    int numJet = int( selectedJets.size() );
    int numTag = int( selectedJets_tag.size() );
    
    double jet1_pt = ( numJet>=1 ) ? selectedJets.at(0).pt() : -99;
  double jet2_pt = ( numJet>=2 ) ? selectedJets.at(1).pt() : -99;
  double jet3_pt = ( numJet>=3 ) ? selectedJets.at(2).pt() : -99;
  double jet4_pt = ( numJet>=4 ) ? selectedJets.at(3).pt() : -99;


  double jet1_CSVv2 = ( numJet>=1 ) ? miniAODhelper.GetJetCSV(selectedJets.at(0), "pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  double jet2_CSVv2 = ( numJet>=2 ) ? miniAODhelper.GetJetCSV(selectedJets.at(1), "pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  double jet3_CSVv2 = ( numJet>=3 ) ? miniAODhelper.GetJetCSV(selectedJets.at(2), "pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;
  double jet4_CSVv2 = ( numJet>=4 ) ? miniAODhelper.GetJetCSV(selectedJets.at(3), "pfCombinedInclusiveSecondaryVertexV2BJetTags") : -99;

    
   

  ///// Boosted jet information

  ///// HEP top tagged jet
  int numTopTags = 0;
  int  n_fatjets=0;
  double  pt_fatjet_1=0, pt_fatjet_2=0;
  double  pt_nonW_1=0, pt_nonW_2=0;
  double pt_W1_1=0, pt_W1_2=0;
  double pt_W2_1=0, pt_W2_2=0;
  double m_top_1=0, m_top_2=0;
  double pt_top_1=0, pt_top_2=0;
  ///// Higgs tagged jet
  int numHiggsTags = 0;
  double higgstag_fatjet_1=0;
  double higgstag_fatjet_2=0;
  double csv2_fatjet_1=0, csv2_fatjet_2=0;
  
  boosted::BoostedJetCollection selectedBoostedJets;
  if(h_boostedjet.isValid()){
  boosted::BoostedJetCollection const &boostedjets_unsorted = *h_boostedjet;
    boosted::BoostedJetCollection boostedjets = BoostedUtils::GetSortedByPt(boostedjets_unsorted);
    selectedBoostedJets = miniAODhelper.GetSelectedBoostedJets(boostedjets, 200., 2.0, 20., 2.4, jetID::jetLoose);
    
  
  vector<boosted::BoostedJet> syncTopJets;
 // if( h_htttopjet.isValid() ){
    //boosted::HTTTopJetCollection const &htttopjets_unsorted = *h_htttopjet;
    //boosted::HTTTopJetCollection htttopjets = BoostedUtils::GetSortedByPt(htttopjets_unsorted);

    // int itop = 0;
    for( boosted::BoostedJetCollection::iterator topJet = boostedjets.begin() ; topJet != boostedjets.end(); ++topJet ){
      // itop++;
    // n_fatjets++;
      // pt and eta requirements on top jet
      if( !(topJet->fatjet.pt() > 200. && abs(topJet->fatjet.eta()) < 2) ) continue;
n_fatjets++;
if(n_fatjets==1)pt_fatjet_1=topJet->fatjet.pt();
if(n_fatjets==2)pt_fatjet_2=topJet->fatjet.pt();
      //if( !(topJet->topjet.pt() > 200. && abs(topJet->topjet.eta()) < 2) ) continue;

      // pt and eta requirements on subjets
      if( !( (topJet->nonW.pt()>20 && abs(topJet->nonW.eta())<2.4 ) &&
	     (topJet->W1.pt()>20 && abs(topJet->W1.eta())<2.4 ) &&
	     (topJet->W2.pt()>20 && abs(topJet->W2.eta())<2.4 ) ) ) continue;

     // n_fatjets++;
      if (n_fatjets == 1){
	//pt_fatjet_1 = topJet->fatjet.pt();
	//cout<<n_fatjets<<" "<<pt_fatjet_1<<" "<<endl;
	pt_nonW_1 = topJet->nonW.pt();
	pt_W1_1 = topJet->W1.pt() ;
	pt_W2_1 = topJet->W2.pt() ;
	m_top_1 = topJet->topMass ;
	pt_top_1 = topJet->topMass;
      }
      if (n_fatjets ==2 ){
	pt_fatjet_2 = topJet->fatjet.pt();
	pt_nonW_2 = topJet->nonW.pt();
	pt_W1_2 = topJet->W1.pt() ;
	pt_W2_2 = topJet->W2.pt() ;
	m_top_2 = topJet->topMass ;
	pt_top_2 = topJet->topMass ;
      }
 


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

  

  
  //if( h_subfilterjet.isValid() ){
   // boosted::SubFilterJetCollection const &subfilterjets_unsorted = *h_subfilterjet;
    //boosted::SubFilterJetCollection subfilterjets = BoostedUtils::GetSortedByPt(subfilterjets_unsorted);

    for( boosted::BoostedJetCollection::iterator higgsJet = boostedjets.begin() ; higgsJet != boostedjets.end(); ++higgsJet ){
      // pt and eta requirements on top jet
      if( !(higgsJet->fatjet.pt() > 200. && abs(higgsJet->fatjet.eta()) < 2) ) continue;
      numHiggsTags++;
      if(numHiggsTags==1) higgstag_fatjet_1 = higgsJet->fatjet.pt(); 
      if(numHiggsTags==2) higgstag_fatjet_2 = higgsJet->fatjet.pt(); 

      //remove overlap with topjets
       bool overlapping=false;
       for(auto tj=syncTopJets.begin(); tj!=syncTopJets.end(); tj++){
       	if(BoostedUtils::DeltaR(tj->fatjet,higgsJet->fatjet)<1.5){
       	  overlapping=true;
       	  break;
       	}
       }
       if(overlapping) continue;
       if(overlapping) continue;
    std::vector<pat::Jet> filterjets = higgsJet->filterjets;
    int subjettags=0;
    for(auto j=filterjets.begin(); j!=filterjets.end(); j++ ){
      if(j->pt()<20 || fabs(j->eta())>2.4) continue;
      //if(BoostedUtils::PassesCSV(*j)){
	//subjettags++;
      //}	
    }
    if(subjettags>=2) numHiggsTags++;

      // int numBtagFiltJets=0;
      // std::vector<pat::Jet> filterjets = higgsJet->filterjets;
      // int numFiltJets = filterjets.size();
      // for( int ijet=0; ijet<numFiltJets; ijet++ ){
      // 	// if( verbose_ ){
      // 	//   printf("\t\t filt jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = %.2f,\t CSVv2 = %+5.3f,\t CSVv1 = %+5.3f \n",
      // 	// 	 ijet, filterjets[ijet].pt(), filterjets[ijet].eta(), filterjets[ijet].phi(), 
      // 	// 	 filterjets[ijet].bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"),
      // 	// 	 filterjets[ijet].bDiscriminator("combinedSecondaryVertexBJetTags"));
      // 	// }

      // 	if( !(filterjets[ijet].pt()>20. && abs(filterjets[ijet].eta()) < 2.4) ) continue;
      // 	if( !(filterjets[ijet].bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.800) ) continue;
      // 	numBtagFiltJets++;
      // }

      // if( verbose_ ){
      // 	printf("\t Higgs jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = %.2f,\t numFiltJets = %2d,\t numBtagFiltJets = %2d\n",
      // 	       int(higgsJet - subfilterjets.begin()), higgsJet->fatjet.pt(), higgsJet->fatjet.eta(), higgsJet->fatjet.phi(), numFiltJets, numBtagFiltJets );
      // }

      // if( numBtagFiltJets>=2 ) numHiggsTags++;

    }
    
    
    if(n_fatjets>0){
    
    higgstag_fatjet_1 = selectedBoostedJets.at(0).fatjet.bDiscriminator("pfBoostedDoubleSecondaryVertexCA15BJetTags");
    
    if(selectedBoostedJets.at(0).filterjets.size()>1){
      std::vector<pat::Jet> filterjets = BoostedUtils::GetHiggsFilterJets(selectedBoostedJets.at(0));
      csv2_fatjet_1 = MiniAODHelper::GetJetCSV(filterjets.at(1));
    }
  }
  
  if(n_fatjets>1){
  
    higgstag_fatjet_2 = selectedBoostedJets.at(1).fatjet.bDiscriminator("pfBoostedDoubleSecondaryVertexCA15BJetTags");
    
    if(selectedBoostedJets.at(1).filterjets.size()>1){
      std::vector<pat::Jet> filterjets = BoostedUtils::GetHiggsFilterJets(selectedBoostedJets.at(1));
      csv2_fatjet_2 = MiniAODHelper::GetJetCSV(filterjets.at(1));
    }
  }
    
  }




         
    // Get Corrected MET (propagating JEC and JER)
    // pat::MET correctedMET = pfmet->front();//miniAODhelper.GetCorrectedMET( pfmets.at(0), pfJets_forMET, iSysType );
    std::vector<pat::Jet> oldJetsForMET = miniAODhelper.GetSelectedJets(*pfjets, 0., 999, jetID::jetMETcorrection, '-' );
    std::vector<pat::Jet> oldJetsForMET_uncorr = miniAODhelper.GetUncorrectedJets(oldJetsForMET);
    std::vector<pat::Jet> newJetsForMET = miniAODhelper.GetCorrectedJets(oldJetsForMET_uncorr, iEvent, iSetup, iSysType);
    std::vector<pat::MET> newMETs = miniAODhelper.CorrectMET(oldJetsForMET, newJetsForMET, *pfmet);

    pat::MET correctedMET = newMETs.at(0); 
    TLorentzVector metV(correctedMET.px(),correctedMET.py(),0.0,correctedMET.pt());

       
    // Initialize event vars, selected jets
    double mht_px = 0.;
    double mht_py = 0.;
    eve->HT_[iSys] = 0.;


    double min_dR_lep_jet = 99.;
    vecTLorentzVector jetV;
    std::vector<double> csvV;
    std::vector<double> jet_combinedMVABJetTags;
    std::vector<double> jet_combinedInclusiveSecondaryVertexV2BJetTags;
    std::vector<double> jet_vtxMass;
    std::vector<double> jet_vtxNtracks;
    std::vector<double> jet_vtx3DVal;
    std::vector<double> jet_vtx3DSig;

    double numtag = 0, numuntag = 0;
    double sum_btag_disc_btags = 0;
    double sum_btag_disc_non_btags = 0;

    double max_m3_pt = -1;
    double max_m3 = -1;
    double max_m3_1tag_pt = -1;
    double max_m3_1tag = -1;

    vvdouble vvjets;
   // vint jetFlavor;	
	
    vdouble dR2Mean_vect;
    vdouble dRMean_vect;
    vdouble frac01_vect;
    vdouble frac02_vect;
    vdouble frac03_vect;
    vdouble beta_vect;
    vdouble betaStar_vect;
    vdouble leadCandDistFromPV_vect;
    vdouble minPVDist_vect;


    vint jet_genId_vect;
    vint jet_partonflavour_vect;
    vint jet_flavour_vect;
    vint jet_genParentId_vect;
    vint jet_genGrandParentId_vect;

int jcntn=0;
    // Loop over selected jets
    for( std::vector<pat::Jet>::const_iterator iJet = selectedJets.begin(); iJet != selectedJets.end(); iJet++ ){ 
jcntn++;
      jet_partonflavour_vect.push_back(iJet->partonFlavour());
      jet_flavour_vect.push_back(iJet->hadronFlavour());

      int genPartonId=-99, genPartonMotherId=-99, genPartonGrandMotherId=-99;
      if( (iJet->genParton()) ){ // if there is a matched parton, fill variables
	genPartonId = iJet->genParton()->pdgId();

	if( iJet->genParton()->numberOfMothers()>=1 ){
	  genPartonMotherId = iJet->genParton()->mother(0)->pdgId();

	  if( iJet->genParton()->mother(0)->numberOfMothers()>=1 ) genPartonMotherId = iJet->genParton()->mother(0)->mother(0)->pdgId();
	}
      }

      jet_genId_vect.push_back(genPartonId);
      jet_genParentId_vect.push_back(genPartonMotherId);
      jet_genGrandParentId_vect.push_back(genPartonGrandMotherId);
      
      

      // DR(lepton system, jet)
      double dR_lep_jet = reco::deltaR(sum_lepton_vect.Eta(),sum_lepton_vect.Phi(),iJet->eta(),iJet->phi());
      if( dR_lep_jet<min_dR_lep_jet ) min_dR_lep_jet = dR_lep_jet;

      // Get jet 4Vector and add to vecTLorentzVector for jets
      TLorentzVector jet0p4;	  
      jet0p4.SetPxPyPzE(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      jetV.push_back(jet0p4);
	  
      // make vvdouble version of vecTLorentzVector
      vdouble vjets;
      vjets.push_back(iJet->px());
      vjets.push_back(iJet->py());
      vjets.push_back(iJet->pz());
      vjets.push_back(iJet->energy());
      vvjets.push_back(vjets);

      jet_vtxMass.push_back(iJet->userFloat("vtxMass"));
      jet_vtxNtracks.push_back(iJet->userFloat("vtxNtracks"));
      jet_vtx3DVal.push_back(iJet->userFloat("vtx3DVal"));
      jet_vtx3DSig.push_back(iJet->userFloat("vtx3DSig"));

      // Get CSV discriminant, check if passes Med WP 
      double myCSV = iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      csvV.push_back(myCSV);
    
      int csvM0 = ( myCSV > 0.8000 ) ? 1 : 0;
      if( myCSV>0.8000 ){
	numtag += 1;
	sum_btag_disc_btags += myCSV;
      }
      else{
	numuntag += 1;
	sum_btag_disc_non_btags += myCSV;
      }
      
      

      jet_combinedMVABJetTags.push_back( iJet->bDiscriminator("combinedMVABJetTags") );
      jet_combinedInclusiveSecondaryVertexV2BJetTags.push_back( iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );

      // Second Loop over selected jets
      for( std::vector<pat::Jet>::const_iterator jJet = iJet; jJet != selectedJets.end(); jJet++ ){ 

	// Continue if same jet as above loop
	if( iJet==jJet ) continue;

	// Get second jet 4Vector and check bTag discriminant
	TLorentzVector jet1p4;
	jet1p4.SetPxPyPzE(jJet->px(),jJet->py(),jJet->pz(),jJet->energy());
	int csvM1 = ( jJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.8000 ) ? 1 : 0;

	// Third loop over selected jets
	for( std::vector<pat::Jet>::const_iterator kJet = jJet; kJet != selectedJets.end(); kJet++ ){ 
     
	  // Continue is third jet is same as the above two jets
	  if( iJet==kJet || jJet==kJet ) continue;

	  // Get third jet 4Vector and chekc bTag discriminant
	  TLorentzVector jet2p4;
	  jet2p4.SetPxPyPzE(kJet->px(),kJet->py(),kJet->pz(),kJet->energy());
	  int csvM2 = ( kJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.8000 ) ? 1 : 0;

	  // Get sum of three jet 4Vectors
	  TLorentzVector sum_jet = jet0p4 + jet1p4 + jet2p4;
	  double temp_pt = sum_jet.Pt();
	  double temp_mass = sum_jet.M();
             
	  // Check if this combination has the highest pT, for M3
	  if( temp_pt>max_m3_pt ){
	    max_m3_pt = temp_pt;
	    max_m3 = temp_mass;
	  }

	  // Check if this combination has highest pT with only 1 bTag
	  if( (csvM0 + csvM1 + csvM2)==1 && temp_pt>max_m3_1tag_pt ){
	    max_m3_1tag_pt = temp_pt;
	    max_m3_1tag = temp_mass;
	  }

	}// end loop over kJet
      }// end loop over jJet
    }// end loop over iJet


    // Compute Avg bTag Discriminant for event
    double ave_btag_disc_non_btags = ( numuntag > 0. ) ? sum_btag_disc_non_btags/numuntag : 0.;
    double ave_btag_disc_btags = ( numtag > 0. ) ? sum_btag_disc_btags/numtag : 0.;


    // Add loose jet container
    std::vector<pat::Jet> selectedJets_loose_unsorted = ( !(iSys>=5 && iSys<=8) ) ? selectedJets_loose_noSys_unsorted : miniAODhelper.GetSelectedJets( correctedJets, 20., 3.0, jetID::none, '-' );
    std::vector<pat::Jet> selectedJets_loose = miniAODhelper.GetSortedByPt( selectedJets_loose_unsorted );

    std::vector<pat::Jet> selectedJets_loose_tag_unsorted = ( !(iSys>=5 && iSys<=8) ) ? selectedJets_loose_tag_noSys_unsorted : miniAODhelper.GetSelectedJets( correctedJets, 20., 3.0, jetID::none, 'M' );

    vvdouble vvjets_loose;
    std::vector<double> csvV_loose;
    std::vector<double> jet_loose_combinedMVABJetTags;
    std::vector<double> jet_loose_combinedInclusiveSecondaryVertexV2BJetTags;
    std::vector<double> jet_loose_vtxMass;
    std::vector<double> jet_loose_vtxNtracks;
    std::vector<double> jet_loose_vtx3DVal;
    std::vector<double> jet_loose_vtx3DSig;
    std::vector<double> jet_loose_pileupJetId_fullDiscriminant;
    vint jet_partonflavour_vect_loose;
    vint jet_flavour_vect_loose;
    vecTLorentzVector jetV_loose;

    int numJet_loose = 0;
    int numTag_loose = 0;

    vvdouble jet_all_vect_TLV;
    vdouble  jet_all_CSV;
    vint     jet_all_flavour;

    // Loop over selected jets
    for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_loose.begin(); iJet != selectedJets_loose.end(); iJet++ ){ 

      vdouble vjets_loose;
      vjets_loose.push_back(iJet->px());
      vjets_loose.push_back(iJet->py());
      vjets_loose.push_back(iJet->pz());
      vjets_loose.push_back(iJet->energy());
      jet_all_vect_TLV.push_back(vjets_loose);

      double myCSV = iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      jet_all_CSV.push_back(myCSV);

      jet_all_flavour.push_back(iJet->partonFlavour());

      jet_loose_combinedMVABJetTags.push_back( iJet->bDiscriminator("combinedMVABJetTags") );
      jet_loose_combinedInclusiveSecondaryVertexV2BJetTags.push_back( iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );

      // MHT
      mht_px += - iJet->px();
      mht_py += - iJet->py();
      eve->HT_[iSys] += iJet->pt();


      // loose represents all jets with pT > 20
      //if( iJet->pt()>=30. ) continue;

      numJet_loose++;
      jet_partonflavour_vect_loose.push_back(iJet->partonFlavour());
      jet_flavour_vect_loose.push_back(iJet->hadronFlavour());

      jet_loose_vtxMass.push_back(iJet->userFloat("vtxMass"));
      jet_loose_vtxNtracks.push_back(iJet->userFloat("vtxNtracks"));
      jet_loose_vtx3DVal.push_back(iJet->userFloat("vtx3DVal"));
      jet_loose_vtx3DSig.push_back(iJet->userFloat("vtx3DSig"));

      jet_loose_pileupJetId_fullDiscriminant.push_back(iJet->userFloat("pileupJetId:fullDiscriminant"));

      TLorentzVector jet0p4;	  
      jet0p4.SetPxPyPzE(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
      jetV_loose.push_back(jet0p4);

      // make vvdouble version of vecTLorentzVector
      vvjets_loose.push_back(vjets_loose);

      // Get CSV discriminant, check if passes Med WP 
      csvV_loose.push_back(myCSV);
      if( myCSV>0.8000 ) numTag_loose++;
    }




 


  /// DIL specific, doesn't make sense in current scope
    // Add lepton 4Vector quantities to MHT
    mht_px += - sum_lepton_vect.Px();
    mht_py += - sum_lepton_vect.Py();
    double MHT = sqrt( mht_px*mht_px + mht_py*mht_py );

    bool PassBigDiamondZmask = ( MuonElectron || 
				 (mass_leplep < (65.5 + 3*MHT/8)) || 
				 (mass_leplep > (108 - MHT/4)) || 
				 (mass_leplep < (79 - 3*MHT/4)) || 
				 (mass_leplep > (99 + MHT/2)) 
				 );
    eve->PassZmask_ = ( PassBigDiamondZmask ) ? 1 : 0;


    double MET = correctedMET.pt();
    bool PassBigDiamondZmaskMET = ( MuonElectron || 
				    (mass_leplep < (65.5 + 3*MET/8)) || 
				    (mass_leplep > (108 - MET/4)) || 
				    (mass_leplep < (79 - 3*MET/4)) || 
				    (mass_leplep > (99 + MET/2)) 
				    );
    eve->PassZmaskMET_ = ( PassBigDiamondZmaskMET ) ? 1 : 0;


    // Get sorted vector of bTags
    std::vector<double> csvV_temp;
    for( int j=0; j<int(csvV.size()); j++ ) csvV_temp.push_back(csvV[j]);

    std::sort(csvV_temp.begin(),csvV_temp.end());
    std::reverse(csvV_temp.begin(),csvV_temp.end());



    // Intialize Event Variables for Tagged Jets
    double min_btag = 999;
    double min_dR_tag_tag = 99.;
    double sum_dR_tag_tag = 0.;
    double sum_mass_tag_tag = 0.;
    double num_tag_tag = 0;
    double sum_dev_from_ave_disc_btags = 0;
    std::vector<double> good_pfjet_tag_pt;
    double tag_tag_dijet_mass_min_dR = 0;
    double tag_tag_dijet_mass_closest_to_125 = 0;
    double min_dR_tag_lep = 99.;
    double min_tag_tag_dijet_mass_closest_to_125 = 99.;
    double mlb_temp = -1;

    // Loop over selected, tagged jets
    for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_tag.begin(); iJet != selectedJets_tag.end(); iJet++ ){ 

      // Get bTag Value
      double myCSV = iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

      // Compute Deviation from Avg bTag
      double dev = myCSV - ave_btag_disc_btags;
      sum_dev_from_ave_disc_btags += dev*dev;

      // Get lowest bTag
      min_btag = std::min( min_btag, myCSV );

      // Get this jet's 4Vector
      TLorentzVector jet0p4;
      jet0p4.SetPxPyPzE(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());

      // Min dR(lep, tag), M(lep, tag)
      double dR_lep_tag = jet0p4.DeltaR(leptonV);
      if( dR_lep_tag<min_dR_tag_lep ){
	min_dR_tag_lep = dR_lep_tag;
	TLorentzVector sum_lep_b = leptonV + jet0p4;//leptonV + jet0p4;
	mlb_temp = sum_lep_b.M();
      }

      // Second loop over selected, tagged jets
      for( std::vector<pat::Jet>::const_iterator jJet = iJet; jJet != selectedJets_tag.end(); jJet++ ){ 

	// Continue if jet is same as above
	if( iJet==jJet ) continue;

	// Get second jet's 4Vector
	TLorentzVector jet1p4;
	jet1p4.SetPxPyPzE(jJet->px(),jJet->py(),jJet->pz(),jJet->energy());

	// Compute this combination's diJet mass
	TLorentzVector diJet = jet0p4 + jet1p4;
	double diJetMass = diJet.M();

	// Compute min dR(tag, tag)
	double dR_tag_tag = jet0p4.DeltaR(jet1p4);
	if( dR_tag_tag<min_dR_tag_tag ){
	  min_dR_tag_tag = dR_tag_tag;
	  tag_tag_dijet_mass_min_dR = diJetMass;
	}

	// compute min Mass(tag,tag, closest to M=125)
	if( fabs(diJetMass-125)<min_tag_tag_dijet_mass_closest_to_125 ){
	  min_tag_tag_dijet_mass_closest_to_125 = fabs(diJetMass-125);
	  tag_tag_dijet_mass_closest_to_125 = diJetMass;
	}

	// Sum tag, tag quantities
	sum_dR_tag_tag += dR_tag_tag;
	num_tag_tag++;

	sum_mass_tag_tag += diJetMass;

      }// end jJet loop
    }// end iJet loop


    // Compute deviation from avg discriminant bTags
    double ave_dev_from_ave_disc_btags = ( numtag > 0. ) ? sum_dev_from_ave_disc_btags/numtag : 0.;

    // Compute avg dR(tag,tag)
    double ave_dR_tag_tag = ( num_tag_tag>0 ) ? sum_dR_tag_tag/num_tag_tag : -1;
    // Compute avg M(tag,tag)
    double ave_mass_tag_tag = ( num_tag_tag>0 ) ? sum_mass_tag_tag/num_tag_tag : -1;



    // Intialize Event Variables for Untagged Jets
    double min_dR_untag_untag = 99.;
    double sum_dR_untag_untag = 0.;
    double sum_mass_untag_untag = 0.;
    double num_untag_untag = 0;
    double sum_dev_from_ave_disc_non_btags = 0;

    // Loop over selected, untagged jets
    for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_untag.begin(); iJet != selectedJets_untag.end(); iJet++ ){ 

      // Get CSV discriminant
      double myCSV = iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");

      // Compute deviation form Avg bTag Disc
      double dev = myCSV - ave_btag_disc_non_btags;
      sum_dev_from_ave_disc_non_btags += dev*dev;

      // Get this jet's 4Vect
      TLorentzVector jet0p4;
      jet0p4.SetPxPyPzE(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());

      // Second loop over selected, untagged jets
      for( std::vector<pat::Jet>::const_iterator jJet = iJet; jJet != selectedJets_untag.end(); jJet++ ){ 

	// Continue if this jet matches the jet above
	if( iJet==jJet ) continue;

	// Get this jet's 4Vector
	TLorentzVector jet1p4;
	jet1p4.SetPxPyPzE(jJet->px(),jJet->py(),jJet->pz(),jJet->energy());

	// Compute min dR(untag, untag)
	double dR_untag_untag = jet0p4.DeltaR(jet1p4);
	if( dR_untag_untag<min_dR_untag_untag ) min_dR_untag_untag = dR_untag_untag;
	// Compute summed untag variables
	sum_dR_untag_untag += dR_untag_untag;
	num_untag_untag++;

	// Compute untagged diJet mass
	TLorentzVector diJet = jet0p4 + jet1p4;
	double diJetMass = diJet.M();

	sum_mass_untag_untag += diJetMass;

      }// end jJet loop
    }// end iJet loop


    // Compute avg M(untag, untag)
    double ave_mass_untag_untag = ( num_untag_untag>0 ) ? sum_mass_untag_untag/num_untag_untag : -1;


   
    // Sum 4Vectors of all objects in event
   // TLorentzVector all_objects = sum_lepton_vect + metV;//leptonV + metV;
	TLorentzVector all_objects = leptonV + metV;
	
    for( int j=0; j<int(jetV.size()); j++ ) all_objects += jetV[j];
	  

    // Compute MHT_
    eve->MHT_[iSys] = sqrt( mht_px*mht_px + mht_py*mht_py );
    eve->MHT_phi_[iSys] = std::atan2(mht_py,mht_px);


    // Comput Ht+Lepton
    //double ht_lep = HT_[iSys] + leptonV.Pt();


    // Compute angular quantities
    // getSp(leptonV,metV,jetV,eve->aplanarity_[iSys],eve->sphericity_[iSys]);
    // getFox(jetV,eve->h0_[iSys],eve->h1_[iSys],eve->h2_[iSys],eve->h3_[iSys],eve->h4_[iSys]);

    // Intialize 4Vectors, to use in Best Higgs Mass function
    TLorentzVector mydummyguy;
    mydummyguy.SetPxPyPzE(0,0,0,0);
    TLorentzVector bjet1_tmp = mydummyguy;
    TLorentzVector bjet2_tmp = mydummyguy;

    // Get Best Higgs Mass (X^2 method)
    // eve->best_higgs_mass_[iSys] = getBestHiggsMass(leptonV, metV, jetV, csvV, eve->minChi2_[iSys], eve->dRbb_[iSys], bjet1_tmp, bjet2_tmp, jetV_loose, csvV_loose);

    //double chiSq2=10000;

    // eve->minChi2_bjet1_px_[iSys] = bjet1_tmp.Px();
    // eve->minChi2_bjet1_py_[iSys] = bjet1_tmp.Py();
    // eve->minChi2_bjet1_pz_[iSys] = bjet1_tmp.Pz();
    // eve->minChi2_bjet1_E_[iSys]  = bjet1_tmp.E();

    // eve->minChi2_bjet2_px_[iSys] = bjet2_tmp.Px();
    // eve->minChi2_bjet2_py_[iSys] = bjet2_tmp.Py();
    // eve->minChi2_bjet2_pz_[iSys] = bjet2_tmp.Pz();
    // eve->minChi2_bjet2_E_[iSys]  = bjet2_tmp.E();


    // nJets/Tags
    eve->numJets_float_[iSys] = numJet;
    eve->numTags_float_[iSys] = numTag;


    // jet1-4 pT
    if( selectedJets.size()>0 ) eve->first_jet_pt_[iSys]  = selectedJets.at(0).pt();
    if( selectedJets.size()>1 ) eve->second_jet_pt_[iSys] = selectedJets.at(1).pt();
    if( selectedJets.size()>2 ) eve->third_jet_pt_[iSys]  = selectedJets.at(2).pt();
    if( selectedJets.size()>3 ) eve->fourth_jet_pt_[iSys] = selectedJets.at(3).pt();


    // MET
    eve->MET_[iSys]      = correctedMET.pt();
    eve->uMET_[iSys]     = correctedMET.pt();//pfmet->front().shiftedPt(pat::MET::NoShift, pat::MET::Raw);
    eve->MET_phi_[iSys]  = correctedMET.phi();
    eve->uMET_phi_[iSys] = correctedMET.phi();//pfmet->front().shiftedPhi(pat::MET::NoShift, pat::MET::Raw);


    // Multi-object Quantities
    eve->all_sum_pt_[iSys]                   = eve->HT_[iSys]+sum_lepton_vect.Pt();//leptonV.Pt();
    eve->all_sum_pt_with_met_[iSys]          = eve->HT_[iSys]+sum_lepton_vect.Pt()+correctedMET.pt();;//leptonV.Pt()+correctedMET.pt();
    eve->M3_[iSys]                           = max_m3;
    eve->M3_1tag_[iSys]                      = max_m3_1tag;
    eve->Mlb_[iSys]     = mlb_temp;
    eve->invariant_mass_of_everything_[iSys] = all_objects.M();


    // dR Quantities
    eve->min_dr_tagged_jets_[iSys]             = min_dR_tag_tag;
    eve->avg_dr_tagged_jets_[iSys]             = ave_dR_tag_tag;
    eve->dr_between_lep_and_closest_jet_[iSys] = min_dR_lep_jet;


    // diJet Quantities 
    eve->avg_tagged_dijet_mass_[iSys]     = ave_mass_tag_tag;
    eve->avg_untagged_dijet_mass_[iSys]   = ave_mass_untag_untag;
    eve->closest_tagged_dijet_mass_[iSys] = tag_tag_dijet_mass_min_dR;
    eve->tagged_dijet_mass_closest_to_125_[iSys] = tag_tag_dijet_mass_closest_to_125;
       
	
    // bTag Quantities
    eve->avg_btag_disc_btags_[iSys]     = ave_btag_disc_btags; 
    eve->avg_btag_disc_non_btags_[iSys] = ave_btag_disc_non_btags; 
    eve->dev_from_avg_disc_btags_[iSys] = ave_dev_from_ave_disc_btags;
    if( csvV_temp.size()>0 ) eve->first_highest_btag_[iSys]  = csvV_temp[0];
    if( csvV_temp.size()>1 ) eve->second_highest_btag_[iSys] = csvV_temp[1];
    if( csvV_temp.size()>2 ) eve->third_highest_btag_[iSys]  = csvV_temp[2];
    if( csvV_temp.size()>3 ) eve->fourth_highest_btag_[iSys] = csvV_temp[3];
    if( csvV_temp.size()>4 ) eve->fifth_highest_CSV_[iSys]   = csvV_temp[4];
    if( csvV_temp.size()>5 ) eve->sixth_highest_CSV_[iSys]   = csvV_temp[5];
    eve->lowest_btag_[iSys]             = min_btag;


    //getBDTvars
    bdtVARS.getSp(leptonV,metV,jetV,eve->aplanarity_[iSys],eve->sphericity_[iSys]);
    bdtVARS.getFox(jetV,eve->h0_[iSys],eve->h1_[iSys],eve->h2_[iSys],eve->h3_[iSys],eve->h4_[iSys]);
    eve->best_higgs_mass_[iSys] = bdtVARS.getBestHiggsMass(leptonV, metV, jetV, csvV, eve->minChi2_[iSys], eve->dRbb_[iSys], bjet1_tmp, bjet2_tmp, jetV_loose, csvV_loose);

    eve->maxeta_jet_jet_[iSys] = bdtVARS.get_jet_jet_etamax(vvjets);
    eve->maxeta_jet_tag_[iSys] = bdtVARS.get_jet_tag_etamax(vvjets,csvV);
    eve->maxeta_tag_tag_[iSys] = bdtVARS.get_tag_tag_etamax(vvjets,csvV);
    
    //eve->maxdeta_JetAvgJet_[iSys] = bdtVARS.JetDelta_EtaAvgEta(vvjets,csvV,"Jet","Jet");
    //eve->maxdeta_TagAvgJet_[iSys] = bdtVARS.JetDelta_EtaAvgEta(vvjets,csvV,"Tag","Jet");
    //eve->maxdeta_TagAvgTag_[iSys] = bdtVARS.JetDelta_EtaAvgEta(vvjets,csvV,"Tag","Tag");
    //eve->maxdeta_JetAvgTag_[iSys] = bdtVARS.JetDelta_EtaAvgEta(vvjets,csvV,"Jet","Tag");
    

    
    
    eve->median_bb_mass_[iSys]  = bdtVARS.get_median_bb_mass(vvjets,csvV);
    eve->pt_all_jets_over_E_all_jets_[iSys]  = bdtVARS.pt_E_ratio_jets(vvjets);
    
    
    TLorentzVector newmet, b1, b2;
    double chi2lepW, chi2leptop, chi2hadW, chi2hadtop;
    double mass_lepW, mass_leptop, mass_hadW, mass_hadtop, drbb_dummy, minChi;
    double bbleptopeta, bbhadtopeta, dphifn, bbeta, avgetatops, detafn, anglbbtops;	
    double MET_phi = correctedMET.phi();
    double MET_pt = correctedMET.pt();
    
    bdtVARS.study_tops_bb_syst(MET_pt, MET_phi, newmet, leptonV, vvjets, 
			       csvV, minChi, chi2lepW, chi2leptop, chi2hadW, chi2hadtop, mass_lepW, mass_leptop, mass_hadW, mass_hadtop, 
			       drbb_dummy, bbleptopeta, bbhadtopeta, dphifn, bbeta, avgetatops, detafn, anglbbtops, b1, b2);	
    
    eve->dEta_leptop_bb_[iSys] = bbhadtopeta;
    eve->dEta_hadtop_bb_[iSys] = bbleptopeta;
    eve->dEta_f_[iSys] = detafn;
    ////

    eve->minChi2_bjet1_px_[iSys] = bjet1_tmp.Px();
    eve->minChi2_bjet1_py_[iSys] = bjet1_tmp.Py();
    eve->minChi2_bjet1_pz_[iSys] = bjet1_tmp.Pz();
    eve->minChi2_bjet1_E_[iSys]  = bjet1_tmp.E();

    eve->minChi2_bjet2_px_[iSys] = bjet2_tmp.Px();
    eve->minChi2_bjet2_py_[iSys] = bjet2_tmp.Py();
    eve->minChi2_bjet2_pz_[iSys] = bjet2_tmp.Pz();
    eve->minChi2_bjet2_E_[iSys]  = bjet2_tmp.E();

    eve->jet_vect_TLV_[iSys]   = vvjets;
    eve->jet_CSV_[iSys]        = csvV;
    eve->jet_combinedMVABJetTags_[iSys] = jet_combinedMVABJetTags;
    eve->jet_combinedInclusiveSecondaryVertexV2BJetTags_[iSys] = jet_combinedInclusiveSecondaryVertexV2BJetTags;

    eve->jet_vtxMass_[iSys]    = jet_vtxMass;
    eve->jet_vtxNtracks_[iSys] = jet_vtxNtracks;
    eve->jet_vtx3DVal_[iSys]   = jet_vtx3DVal;
    eve->jet_vtx3DSig_[iSys]   = jet_vtx3DSig;
    eve->jet_partonflavour_[iSys]          = jet_partonflavour_vect;
    eve->jet_flavour_[iSys]          = jet_flavour_vect;
    eve->jet_genId_[iSys]            = jet_genId_vect;
    eve->jet_genParentId_[iSys]      = jet_genParentId_vect;
    eve->jet_genGrandParentId_[iSys] = jet_genGrandParentId_vect;

    // loose jets

    // nJets/Tags
    eve->numJets_loose_float_[iSys] = numJet_loose;
    eve->numTags_loose_float_[iSys] = numTag_loose;

    eve->jet_loose_vect_TLV_[iSys] = vvjets_loose;
    eve->jet_loose_CSV_[iSys]      = csvV_loose;
    eve->jet_loose_combinedMVABJetTags_[iSys] = jet_loose_combinedMVABJetTags;
    eve->jet_loose_combinedInclusiveSecondaryVertexV2BJetTags_[iSys] = jet_loose_combinedInclusiveSecondaryVertexV2BJetTags;
    eve->jet_loose_vtxMass_[iSys]    = jet_loose_vtxMass;
    eve->jet_loose_vtxNtracks_[iSys] = jet_loose_vtxNtracks;
    eve->jet_loose_vtx3DVal_[iSys]   = jet_loose_vtx3DVal;
    eve->jet_loose_vtx3DSig_[iSys]   = jet_loose_vtx3DSig;
    eve->jet_loose_pileupJetId_fullDiscriminant_[iSys] = jet_loose_pileupJetId_fullDiscriminant;
    eve->jet_loose_partonflavour_[iSys]  = jet_partonflavour_vect_loose;
    eve->jet_loose_flavour_[iSys]  = jet_flavour_vect_loose;


   



    double csvWgtHF=1, csvWgtLF=1, csvWgtCF=1;
    double csvWgt = 1;//get_csv_wgt(jet_all_vect_TLV, jet_all_CSV, jet_all_flavour, iSys, csvWgtHF, csvWgtLF, csvWgtCF);
   double bWeight = 1;
   // double bWeight = get_csv_wgt( vvjets, csvV, jet_flavour_vect, 0, csvWgtHF, csvWgtLF, csvWgtCF);
    if( insample_<0 ){
      csvWgt = 1.; csvWgtHF = 1.; csvWgtLF = 1.; csvWgtCF = 1.;
    }

    eve->wgt_csvSF_[iSys]    = csvWgt;
    eve->wgt_csvSF_HF_[iSys] = csvWgtHF;
    eve->wgt_csvSF_LF_[iSys] = csvWgtLF;
    eve->wgt_csvSF_CF_[iSys] = csvWgtCF;


    eve->wgt_[iSys] *= csvWgt;

    if( iSys==0 ) totalWgt *= csvWgt;


    unsigned int minNumLooseJet = ( isLJ_ ) ? 4 : 1;
    unsigned int minNumLooseTag = ( isLJ_ ) ? 1 : 0;
    
     bool passedSyncSelections = false;
  
   if(passSingleMuonTrigger && numTightMuons==1 && numLooseMuons==1 && numLooseElectrons==0 && numJet>=4 && numTag>=2)passedSyncSelections=true;
   if(passSingleElectronTrigger && numTightElectrons==1 && numLooseElectrons==1 && numLooseMuons==0 && numJet>=4 && numTag>=2)passedSyncSelections=true;
 if(passedSyncSelections && false && iSys==0){
         int met_passed=0;
	//	if(m_ll>20)mll_pass=1;
		if(MET_pt>40)met_passed=1; //printf("%d,%d,%d,1,0,%.3f,%+.3f,%+.3f,%.3f,%d,0,0,0,0,0,%.3f,%.3f,%.3f,%.3f,%+.3f,%+.3f,%+.3f,%+.3f,%.3f,%+.3f,%d,%d,%.3f,%d,0,0,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
	     cout<<std::setprecision(6)<<run<<",";
	     cout<<std::setprecision(6)<<lumi<<",";
	     cout<<std::setprecision(6)<<evt<<",1,0,";
	     cout<<std::setprecision(6)<<eve->tight_lepton_pt_<<",";
	     cout<<std::setprecision(6)<<eve->tight_lepton_eta_<<",";
	     cout<<std::setprecision(6)<<eve->tight_lepton_phi_<<",";
	     cout<<std::setprecision(6)<<tight_lepton_relIso<<",";
	     cout<<std::setprecision(6)<<lep_genId<<",0,0,0,0,0,0,0,";
	     cout<<std::setprecision(6)<<jet1_pt<<",";
	     cout<<std::setprecision(6)<<jet2_pt<<",";
	     cout<<std::setprecision(6)<<jet3_pt<<",";
	     cout<<std::setprecision(6)<<jet4_pt<<",";
	     cout<<std::setprecision(6)<<jet1_CSVv2<<",";
	     cout<<std::setprecision(6)<<jet2_CSVv2<<",";
	     cout<<std::setprecision(6)<<jet3_CSVv2<<",";
	     cout<<std::setprecision(6)<<jet4_CSVv2<<",";
	     cout<<std::setprecision(6)<<MET_pt<<",";
	     cout<<std::setprecision(6)<<MET_phi<<",";
	     cout<<std::setprecision(6)<<met_passed<<",";
	     cout<<std::setprecision(6)<<numJet<<",";
	     cout<<std::setprecision(6)<<numTag<<",";
	     cout<<std::setprecision(6)<<bWeight<<",";
	     //	     cout<<std::setprecision(6)<<additionalJetEventId<<",0,0,";
	     cout<<std::setprecision(6)<<n_fatjets<<",";
	     cout<<std::setprecision(6)<<pt_fatjet_1<<",";
	     cout<<std::setprecision(6)<<pt_fatjet_2<<",";
	     cout<<std::setprecision(6)<<pt_nonW_1<<",";
	     cout<<std::setprecision(6)<<pt_nonW_2<<",";
	     cout<<std::setprecision(6)<<pt_W1_1<<",";
	     cout<<std::setprecision(6)<<pt_W1_2<<",";
	     cout<<std::setprecision(6)<<pt_W2_1<<",";
	     cout<<std::setprecision(6)<<pt_W2_2<<",";
	     cout<<std::setprecision(6)<<pt_top_1<<",";
	     cout<<std::setprecision(6)<<pt_top_2<<",";
	     cout<<std::setprecision(6)<<m_top_1<<",";
	     cout<<std::setprecision(6)<<m_top_2<<",";
	     cout<<std::setprecision(6)<<higgstag_fatjet_1<<",";
	     cout<<std::setprecision(6)<<higgstag_fatjet_2<<",";
	     cout<<std::setprecision(6)<<csv2_fatjet_1<<",";
	     cout<<std::setprecision(6)<<csv2_fatjet_2<<endl;
  }

    // if( iSys>=7 && iSys<=8 ){
    //   for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_loose.begin(); iJet != selectedJets_loose.end(); iJet++ ){ 
    // 	printf("\t selectedJets_loose iSys = %d,\t jet %d: pt = %.2f,\t raw pt = %.2f,\t eta = %.2f,\t phi = %.2f \n", 
    // 	       iSys, int(iJet - selectedJets_loose.begin()), iJet->pt(), iJet->correctedJet(0).pt(), iJet->eta(), iJet->phi() );
    //   }
    //   for( std::vector<pat::Jet>::const_iterator iJet = jetsNoEle.begin(); iJet != jetsNoEle.end(); iJet++ ){ 
    // 	if( !( iJet->pt()>20. && fabs(iJet->eta())<2.4 ) ) continue;
    // 	printf("\t jetsNoEle iSys = %d,\t jet %d: pt = %.2f,\t raw pt = %.2f,\t eta = %.2f,\t phi = %.2f \n", 
    // 	       iSys, int(iJet - jetsNoEle.begin()), iJet->pt(), iJet->correctedJet(0).pt(), iJet->eta(), iJet->phi() );
    //   }
    //   for( std::vector<pat::Jet>::const_iterator iJet = correctedJets_noSys.begin(); iJet != correctedJets_noSys.end(); iJet++ ){ 
    // 	if( !( iJet->pt()>20. && fabs(iJet->eta())<2.4 ) ) continue;
    // 	printf("\t correctedJets_noSy iSys = %d,\t jet %d: pt = %.2f,\t raw pt = %.2f,\t eta = %.2f,\t phi = %.2f \n", 
    // 	       iSys, int(iJet - correctedJets_noSys.begin()), iJet->pt(), iJet->correctedJet(0).pt(), iJet->eta(), iJet->phi() );
    //   }
    //   for( std::vector<pat::Jet>::const_iterator iJet = pfjets->begin(); iJet != pfjets->end(); iJet++ ){ 
    // 	if( !( iJet->pt()>20. && fabs(iJet->eta())<2.4 ) ) continue;
    // 	printf("\t pfjets iSys = %d,\t jet %d: pt = %.2f,\t raw pt = %.2f,\t eta = %.2f,\t phi = %.2f \n", 
    // 	       iSys, int(iJet - pfjets->begin()), iJet->pt(), iJet->correctedJet(0).pt(), iJet->eta(), iJet->phi() );
    //   }
    // }
    if( selectedJets_loose_unsorted.size()>=minNumLooseJet && selectedJets_loose_tag_unsorted.size()>=minNumLooseTag ){
      hasNumJetNumTag = true;
      neventsFillTree[iSys]++;
    }

  } // end loop over systematics


  //
  // Fill tree if pass full selection
  //
  if( hasNumJetNumTag || true ){
    worldTree->Fill();
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
YggdrasilTreeMaker::beginJob()
{
  edm::Service<TFileService> fileService;
  if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");
    
  h_ttbarId_ = fileService->make<TH1F>("ttbarId", "ttbarId", 260, 0, 260);
  h_ttbarAdditionalJetId_ = fileService->make<TH1F>("ttbarAdditionalJetId", "ttbarAdditionalJetId", 60, 0, 60);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
YggdrasilTreeMaker::endJob() 
{

  std::cout << " *********************************************************** " << std::endl;
  if( isLJ_ ) std::cout << " ** Lepton + jets event selection ** " << std::endl;
  else       std::cout << " ** Opposite sign dilepton event selection ** " << std::endl;
  std::cout << " " << std::endl;
  std::cout << "\t Event information: " << std::endl;
  std::cout << "\t\t Sample = " << mySample_sampleName_ << std::endl;
  std::cout << "\t\t Integrated lumi (pb^-1) = " << intLumi_ << std::endl;
  std::cout << "\t\t Cross section (pb)      = " << mySample_xSec_ << std::endl;
  std::cout << "\t\t Number generated        = " << mySample_nGen_ << std::endl;
  std::cout << " " << std::endl;
  std::cout << "   Number of Events Processed  = " << nevents << std::endl;
  std::cout << "   Number of Events Processed (wgt) = " << nevents_wgt << std::endl;

  std::cout << "   Number of Events Pass Trigger  = " << nevents_clean << std::endl;
  std::cout << "   Number of Events Pass Trigger (wgt) = " << nevents_clean_wgt << std::endl;

  if( isLJ_ ){
    std::cout << "   Number of Events with single lepton = " << nevents_1l << " (ele = " << nevents_1l_ele << ", mu = " << nevents_1l_mu << ")" << std::endl;
    std::cout << "   Number of Events with single lepton (wgt) = " << nevents_1l_wgt << " (ele = " << nevents_1l_ele_wgt << ", mu = " << nevents_1l_mu_wgt << ")" << std::endl;
    std::cout << "   Number of Events with single lepton, >=4 jets = " << nevents_1l_4j << " (ele = " << nevents_1l_4j_ele << ", mu = " << nevents_1l_4j_mu << ")" << std::endl;
    std::cout << "   Number of Events with single lepton, >=4 jets (wgt) = " << nevents_1l_4j_wgt << " (ele = " << nevents_1l_4j_ele_wgt << ", mu = " << nevents_1l_4j_mu_wgt << ")" << std::endl;
    std::cout << "   Number of Events with single lepton, >=4 jets, >=2 tags = " << nevents_1l_4j_2t << " (ele = " << nevents_1l_4j_2t_ele << ", mu = " << nevents_1l_4j_2t_mu << ")" << std::endl;
    std::cout << "   Number of Events with single lepton, >=4 jets, >=2 tags (wgt) = " << nevents_1l_4j_2t_wgt << " (ele = " << nevents_1l_4j_2t_ele_wgt << ", mu = " << nevents_1l_4j_2t_mu_wgt << ")" << std::endl;
  }
  else {
    std::cout << "   Number of Events with two leptons = " << nevents_2l << " (2e = " << nevents_2l_2e << ", 2m = " << nevents_2l_2m << ", em = " << nevents_2l_em << ")" << std::endl;
    std::cout << "   Number of Events with two leptons (wgt) = " << nevents_2l_wgt << " (2e = " << nevents_2l_2e_wgt << ", 2m = " << nevents_2l_2m_wgt << ", em = " << nevents_2l_em_wgt << ")" << std::endl;
    std::cout << "   Number of Events with two leptons, >=2 jets = " << nevents_2l_2j << " (2e = " << nevents_2l_2j_2e << ", 2m = " << nevents_2l_2j_2m << ", em = " << nevents_2l_2j_em << ")" << std::endl;
    std::cout << "   Number of Events with two leptons, >=2 jets (wgt) = " << nevents_2l_2j_wgt << " (2e = " << nevents_2l_2j_2e_wgt << ", 2m = " << nevents_2l_2j_2m_wgt << ", em = " << nevents_2l_2j_em_wgt << ")" << std::endl;
    std::cout << "   Number of Events with two leptons, >=2 jets, >=2 tags = " << nevents_2l_2j_2t << " (2e = " << nevents_2l_2j_2t_2e << ", 2m = " << nevents_2l_2j_2t_2m << ", em = " << nevents_2l_2j_2t_em << ")" << std::endl;
    std::cout << "   Number of Events with two leptons, >=2 jets, >=2 tags (wgt) = " << nevents_2l_2j_2t_wgt << " (2e = " << nevents_2l_2j_2t_2e_wgt << ", 2m = " << nevents_2l_2j_2t_2m_wgt << ", em = " << nevents_2l_2j_2t_em_wgt << ")" << std::endl;
  }

  //for( int i=0; i<rNumSys; i++ ) printf(" i = %d,\t Num Pass = %d \n", i, neventsFillTree[i] );
  std::cout << " *********************************************************** " << std::endl;

}

// ------------ method called when starting to processes a run  ------------
/*
void 
YggdrasilTreeMaker::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
void 
YggdrasilTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  bool hltchanged = true;
  if (!hlt_config_.init(iRun, iSetup, hltTag, hltchanged)) {
    std::cout << "Warning, didn't find trigger process HLT,\t" << hltTag << std::endl;
    return;
  }

}


// ------------ method called when ending the processing of a run  ------------
/*
void 
YggdrasilTreeMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
YggdrasilTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
YggdrasilTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
YggdrasilTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

/*
void getSp(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, float &aplanarity, float &sphericity) {
	//
	// Aplanarity and sphericity
	//

	int nJets = int(jets.size());

	float mxx = lepton.Px()*lepton.Px() + met.Px()*met.Px();
	float myy = lepton.Py()*lepton.Py() + met.Py()*met.Py();
	float mzz = lepton.Pz()*lepton.Pz() + met.Pz()*met.Pz();
	float mxy = lepton.Px()*lepton.Py() + met.Px()*met.Py();
	float mxz = lepton.Px()*lepton.Pz() + met.Px()*met.Pz();
	float myz = lepton.Py()*lepton.Pz() + met.Px()*met.Pz();

	for (int i=0; i<nJets; i++) {
		mxx += jets[i].Px()*jets[i].Px();
		myy += jets[i].Py()*jets[i].Py();
		mzz += jets[i].Pz()*jets[i].Pz();
		mxy += jets[i].Px()*jets[i].Py();
		mxz += jets[i].Px()*jets[i].Pz();
		myz += jets[i].Py()*jets[i].Pz();		
	}
	float sum = mxx + myy + mzz;
	mxx /= sum;
	myy /= sum;
	mzz /= sum;
	mxy /= sum;
	mxz /= sum;
	myz /= sum;

	TMatrix tensor(3,3);
	tensor(0,0) = mxx;
	tensor(1,1) = myy;
	tensor(2,2) = mzz;
	tensor(0,1) = mxy;
	tensor(1,0) = mxy;
	tensor(0,2) = mxz;
	tensor(2,0) = mxz;
	tensor(1,2) = myz;
	tensor(2,1) = myz;
	TVector eigenval(3);
	tensor.EigenVectors(eigenval);

	sphericity = 3.0*(eigenval(1)+eigenval(2))/2.0;
	aplanarity = 3.0*eigenval(2)/2.0;

	return;
}

void getFox(vecTLorentzVector jets, float &h0, float &h1, float &h2, float &h3, float &h4) {
	

	int visObjects = int(jets.size());

	float eVis = 0.0;
	for (int i=0; i<visObjects; i++) {
		eVis += jets[i].E();
	}

	h0 = 0.0;
	h1 = 0.0;
	h2 = 0.0;
	h3 = 0.0;
	h4 = 0.0;
	for (int i=0; i<visObjects-1; i++) {
		for (int j=i+1; j<visObjects; j++) {
			float costh = cos(jets[i].Angle(jets[j].Vect()));
			float p0 = 1.0;
			float p1 = costh;
			float p2 = 0.5*(3.0*costh*costh - 1.0);
			float p3 = 0.5*(5.0*costh*costh - 3.0*costh);
			float p4 = 0.125*(35.0*costh*costh*costh*costh - 30.0*costh*costh + 3.0);
			float pipj = jets[i].P()*jets[j].P();
			h0 += (pipj/(eVis*eVis))*p0;
			h1 += (pipj/(eVis*eVis))*p1;
			h2 += (pipj/(eVis*eVis))*p2;
			h3 += (pipj/(eVis*eVis))*p3;
			h4 += (pipj/(eVis*eVis))*p4;
		}
	}

	return;
}


double getBestHiggsMass(TLorentzVector lepton, TLorentzVector met, vecTLorentzVector jets, vdouble btag, double &minChi, double &dRbb, TLorentzVector &bjet1, TLorentzVector &bjet2, vecTLorentzVector loose_jets, vdouble loose_btag)
{

  if( jets.size()<6 && loose_jets.size()>0 ){
    jets.push_back( loose_jets[0] );
    btag.push_back( loose_btag[0] );
  }

  int nJets = int(jets.size());

  double chi_top_lep=10000;
  double chi_top_had=10000;
  //double chi_W_lep=10000; //isn't really used
  double chi_W_had=10000;

  minChi = 1000000;
  dRbb = 1000000;
  double btagCut = 0.814;
  double W_mass = 80.0;
  double top_mass = 172.5;
  //double H_mass=120.0;

  // updated 8/22/2012 from J. Timcheck
  //sigma's from >=6j >=4t, muon, no imaginary neutrino pz ttH
  double sigma_hadW   = 12.77;
  double sigma_hadTop = 18.9;
  double sigma_lepTop = 32.91;

  // //sigma's from >=6j >=4t, muon, no imaginary neutrino pz ttH
  // double sigma_hadW   = 12.59;
  // double sigma_hadTop = 19.9;
  // double sigma_lepTop = 39.05;

  //sigma's from >=6j >=4t, muon, no imaginary neutrino pz ttJets
  //double sigma_hadW		= 12.72,
    // sigma_hadTop	= 18.12,
    // sigma_lepTop	= 38.72;


  double metPz[2];
  double chi=999999;

  //stuff to find:
  double higgs_mass_high_energy=0;

  int nBtags = 0;
  for(int i=0;i<nJets;i++){
    if(btag[i]>btagCut) nBtags++;
  }

  int nUntags = nJets-nBtags;

  double lowest_btag = 99.;
  double second_lowest_btag = 999.;
  int ind_lowest_btag = 999;
  int ind_second_lowest_btag = 999;

  if( nJets>=6 && nBtags>=4 ){
    if( nUntags<2 ){
      for(int i=0;i<nJets;i++){
	if( btag[i]<lowest_btag ){
	  second_lowest_btag = lowest_btag;
	  ind_second_lowest_btag = ind_lowest_btag;

	  lowest_btag = btag[i];
	  ind_lowest_btag = i;
	}
	else if( btag[i]<second_lowest_btag ){
	  second_lowest_btag = btag[i];
	  ind_second_lowest_btag = i;
	}
      }
    }
  }


  //Handle 6j3t.
  int ind_promoted_btag = 999;

  if( nJets>=6 && nBtags==3 ){
    for(int i=0;i<nJets;i++){
      int rank = 0;
      for(int j=0;j<nJets;j++){
	if( btag[j] > btag[i] ){
	  rank++;
	}
      }
      if( rank == 3 ) ind_promoted_btag = i;
    }
  }

  // First get the neutrino z
  double energyLep = lepton.E();
  double a = (W_mass*W_mass/(2.0*energyLep)) + (lepton.Px()*met.Px() + lepton.Py()*met.Py())/energyLep;
  double radical = (2.0*lepton.Pz()*a/energyLep)*(2.0*lepton.Pz()*a/energyLep);
  radical = radical - 4.0*(1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep))*(met.Px()*met.Px() + met.Py()*met.Py()- a*a);
  if (radical < 0.0) radical = 0.0;
  metPz[0] = (lepton.Pz()*a/energyLep) + 0.5*sqrt(radical);
  metPz[0] = metPz[0] / (1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep));
  metPz[1] = (lepton.Pz()*a/energyLep) - 0.5*sqrt(radical);
  metPz[1] = metPz[1] / (1.0 - (lepton.Pz()/energyLep)*(lepton.Pz()/energyLep));


  // Loop over all jets, both Pz, calcaulte chi-square
  TLorentzVector metNew;
  for( int ipznu=0; ipznu<2; ipznu++ ){
    metNew.SetXYZM(met.Px(),met.Py(),metPz[ipznu],0.0); //neutrino has mass 0
    //with b-tag info
    if( (nJets>=6 && nBtags>=4) || (nJets>=6 && nBtags==3) ){
      vecTLorentzVector not_b_tagged,b_tagged;
      //fill not_b_tagged and b_tagged
      for( int i=0;i<nJets;i++ ){
	if( (btag[i]>btagCut && i!=ind_second_lowest_btag && i!=ind_lowest_btag) || (i==ind_promoted_btag) ) b_tagged.push_back(jets[i]);
	else not_b_tagged.push_back(jets[i]);
      }
      //first make possible t_lep's with b-tagged jets (includes making W_lep)
      for( int i=0; i<int(b_tagged.size()); i++ ){
	TLorentzVector W_lep=metNew+lepton; //used for histogram drawing only
	TLorentzVector top_lep=metNew+lepton+b_tagged.at(i);
	chi_top_lep=pow((top_lep.M()-top_mass)/sigma_lepTop,2);
	//next make possible W_had's with not b-tagged jets
	for( int j=0; j<int(not_b_tagged.size()); j++ ){
	  for( int k=0; k<int(not_b_tagged.size()); k++ ){
	    if( j!=k ){
	      TLorentzVector W_had=not_b_tagged.at(j)+not_b_tagged.at(k);
	      chi_W_had=pow((W_had.M()-W_mass)/sigma_hadW,2);
	      //now make possible top_had's (using the W_had + some b-tagged jet)
	      for( int l=0; l<int(b_tagged.size()); l++ ){
		if( l!=i ){
		  TLorentzVector top_had=W_had+b_tagged.at(l);
		  chi_top_had=pow((top_had.M()-top_mass)/sigma_hadTop,2);
		  chi=chi_top_lep+chi_W_had+chi_top_had;
		  //accept the lowest chi
		  if( chi<minChi ){
		    minChi=chi;
		    //pick the other two b's that have the highest et (energy in transverse plane) as higgs mass constituents
		    TLorentzVector H2;
		    int numH2Constituents=0;
		    TLorentzVector bBest[2];
		    for( int m=0; m<int(b_tagged.size()); m++ ){
		      if( m!=i && m!=l && numH2Constituents<2 ){
			bBest[numH2Constituents] = b_tagged.at(m);
			numH2Constituents++;
			H2+=b_tagged.at(m);
		      }
		    }
		    dRbb = bBest[0].DeltaR( bBest[1] );
		    higgs_mass_high_energy=H2.M();
		    bjet1 = bBest[0];
		    bjet2 = bBest[1];
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return higgs_mass_high_energy;
}



*/

//define this as a plug-in
DEFINE_FWK_MODULE(YggdrasilTreeMaker);
