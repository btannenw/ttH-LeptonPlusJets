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


#define EVENTSYNCMODE false
// In case of not EVENTSYNCMODE, event (no lepton, less jets events) can be not recorded in the output tree.

#define DISABLE_LHE_FOR_DIBOSON_SAMPLES false

// system include files
#include <memory>
#include <algorithm>

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
#include "FWCore/Common/interface/TriggerNames.h"


#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

#include "MiniAOD/MiniAODHelper/interface/MiniAODHelper.h"
#include "MiniAOD/MiniAODHelper/interface/TopTagger.h"
#include "MiniAOD/MiniAODHelper/interface/HiggsTagger.h"
#include "ttH-LeptonPlusJets/AnalysisCode/interface/YggdrasilEventVars.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

// #include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"

#include "MiniAOD/MiniAODHelper/interface/BDTvars.h"

#include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilEventSelection.h"
#include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/ttHYggdrasilScaleFactors.h"
#include "ttH-LeptonPlusJets/YggdrasilTreeMaker/interface/RoccoR.h"

#include "LHAPDF/LHAPDF.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

//#include "TTH/CommonClassifier/interface/BDTClassifier.h"


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

  TH1D * h_event;
        
  // Histogram for ttbar event categorization ID based on additional jets only
  TH1* h_ttbarAdditionalJetId_;

  
  // Input tags
  edm::EDGetTokenT<reco::GenJetCollection> genJetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection> fatgenJetsToken_;
  edm::EDGetTokenT<reco::GenJetCollection> reclusteredfatgenJetsToken_;
  
  edm::EDGetTokenT<std::vector<int> > genBHadJetIndexToken_;
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
  edm::EDGetTokenT <pat::TriggerObjectStandAloneCollection> TriggerObjectStandAloneToken ;

  // // new MVAelectron
  // edm::EDGetTokenT< edm::View<pat::Electron> > EDMElectronsToken;
  // // MVA values and categories
  // edm::EDGetTokenT<edm::ValueMap<float> > EDMeleMVAvaluesToken;
  // edm::EDGetTokenT<edm::ValueMap<int> > EDMeleMVAcategoriesToken;

  edm::EDGetTokenT <reco::VertexCollection> vertexToken;
  edm::EDGetTokenT <pat::ElectronCollection> electronToken;
  edm::EDGetTokenT <pat::MuonCollection> muonToken;
  edm::EDGetTokenT < edm::View<pat::Muon> > muonview_Token;

  edm::EDGetTokenT <pat::JetCollection> jetToken;
  edm::EDGetTokenT <pat::JetCollection> jetTokenWithSeeds; // BBT 10-12-18
  edm::EDGetTokenT <pat::JetCollection> jetTokenWithSeeds_v2; // BBT 10-15-18
  edm::EDGetTokenT <pat::JetCollection> puppijetToken;
  edm::EDGetTokenT <pat::JetCollection> fatjetToken;
  edm::EDGetTokenT <pat::JetCollection> rerun_fatjetToken;
  edm::EDGetTokenT <pat::METCollection> metToken;
  edm::EDGetTokenT <pat::METCollection> puppimetToken;

  // edm::EDGetTokenT< boosted::HEPTopJetCollection > topJetsToken;
  // edm::EDGetTokenT< boosted::SubFilterJetCollection > subFilterJetsToken;

  edm::EDGetTokenT <pat::PackedCandidateCollection> packedpfToken;

  edm::EDGetTokenT <reco::BeamSpot> beamspotToken;
  edm::EDGetTokenT <double> rhoToken;
  edm::EDGetTokenT <reco::GenParticleCollection> mcparicleToken;
  edm::EDGetTokenT <std::vector< PileupSummaryInfo > > puInfoToken;

  edm::EDGetTokenT <GenEventInfoProduct> genInfoProductToken;

  edm::EDGetTokenT <LHEEventProduct> LHEEventProductToken;

  // edm::EDGetTokenT <reco::ConversionCollection> EDMConversionCollectionToken;
  edm::EDGetTokenT< boosted::BoostedJetCollection > EDMBoostedJetsToken;

  edm::EDGetTokenT< TtGenEvent >   TtGenEventToken ;
  
  edm::EDGetTokenT<reco::JetCorrector> jetCorrectorToken_;
  edm::EDGetTokenT<reco::JetCorrector> puppijetCorrectorToken_;
  edm::EDGetTokenT<reco::JetCorrector> fatjetCorrectorToken_;

  // BBT, 10-08-18: add token for electron ID
  edm::EDGetTokenT<edm::ValueMap<bool> > eleIdMapToken_;


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

  int nevents;
  double nevents_wgt;

  int nevents_clean;
  double nevents_clean_wgt;

  TTree *worldTree;
  yggdrasilEventVars *eve; 

  // EGammaMvaEleEstimatorCSA14* myMVATrig;
 
  MiniAODHelper miniAODhelper;
  MiniAODHelper miniAODhelper_Puppi;
  MiniAODHelper miniAODhelper_fatjet ; 

  RoccoR * muon_roc ;

  BDTvars bdtVARS;
  //BDTClassifier bdt;
  
  TopTagger toptagger;

  const bool isMC ;
  const bool isTTbarMC ;
  const bool usePUPPI ;

  bool checkIfRegisterd( const reco::Candidate * candidate , std::vector< const reco::Candidate * > list ) ;
  const reco::Candidate * TraceBackToJustAfterBirth( const reco::Candidate * particle );
  void DumpDecay (  const reco::Candidate * c ,  int depth = 0 );



  // - - - - - - PDF uncertainty - - - - - - - - -
  LHAPDF::PDFSet * CT14nlo_PDFSet;
  std::vector<LHAPDF::PDF*> CT14nlo_systPDFs ;         

  LHAPDF::PDFSet * NNPDF30_nlo_as_0118_PDFSet;
  std::vector<LHAPDF::PDF*> NNPDF30_nlo_as_0118_systPDFs ;         

  ttHYggdrasilEventSelection selection;
  ttHYggdrasilScaleFactors   scalefactors;

  bool b_EventSyncFirstEvent;

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
  genJetsToken_ ( consumes <reco::GenJetCollection> ( iConfig.getParameter<edm::InputTag>("genjet") ) )
  ,  isMC(iConfig.getParameter<std::string>("inputfiletype") != "data" )
  ,  isTTbarMC(iConfig.getParameter<std::string>("inputfiletype") == "TTbarMC" )
  , usePUPPI(iConfig.getParameter<std::string>("jetPU") == "PUPPI" )
  , b_EventSyncFirstEvent( true )
{
   //now do what ever initialization is needed
  verbose_ = false;
  isLJ_ = true;

  if( isMC ){
    filterTag = "PAT";
    hltTag    = "HLT";
    filterTag = "RECO"; // BBT 10-19-18
  }else{
    filterTag = "RECO"; // BBT 02-20-19
    hltTag = "HLT";
  }
  triggerResultsToken = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), hltTag));
  filterResultsToken = consumes <edm::TriggerResults> (edm::InputTag(std::string("TriggerResults"), std::string(""), filterTag));

  fatgenJetsToken_ =  consumes <reco::GenJetCollection> (edm::InputTag("slimmedGenJetsAK8","","")) ;
  reclusteredfatgenJetsToken_ =  consumes <reco::GenJetCollection> (edm::InputTag("ak8ReclusteredGenJets","","")) ;

  TriggerObjectStandAloneToken = consumes <pat::TriggerObjectStandAloneCollection>
    //( edm::InputTag( std::string ( "slimmedPatTrigger" ), std::string("") , std::string(isMC ? "PAT" : "RECO") )) ; // ygg core, but no RECO in 2017 data files
    ( edm::InputTag( std::string ( "slimmedPatTrigger" ), std::string("") , std::string(isMC ? "PAT" : "PAT") )) ;  // BBT, 10-29-18
  //    ( edm::InputTag( std::string ( "selectedPatTrigger" ), std::string("") , std::string(isMC ? "PAT" : "RECO") )) ; // ygg core, but commented out by default

  if( isMC ){
    jetCorrectorToken_ = consumes< reco::JetCorrector > (edm::InputTag("ak4PFCHSL1FastL2L3Corrector","","")) ;
  }else{
    jetCorrectorToken_ = consumes< reco::JetCorrector > (edm::InputTag("ak4PFCHSL1FastL2L3ResidualCorrector","","")) ;
  }

  if( isMC ){
    puppijetCorrectorToken_ = consumes< reco::JetCorrector > (edm::InputTag("ak4PFPuppiL1FastL2L3Corrector","","")) ;
  }else{
    puppijetCorrectorToken_ = consumes< reco::JetCorrector > (edm::InputTag("ak4PFPuppiL1FastL2L3ResidualCorrector","","")) ;
  }

  if( isMC ){
    fatjetCorrectorToken_ = consumes< reco::JetCorrector > (edm::InputTag("ak8PFPuppiL1FastL2L3Corrector","","")) ;
  }else{
    fatjetCorrectorToken_ = consumes< reco::JetCorrector > (edm::InputTag("ak8PFPuppiL1FastL2L3ResidualCorrector","","")) ;
  }


  // // new MVAelectron
  // EDMElectronsToken = consumes< edm::View<pat::Electron> >(edm::InputTag("slimmedElectrons","",""));
  // EDMeleMVAvaluesToken           = consumes<edm::ValueMap<float> >(edm::InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values",""));
  // EDMeleMVAcategoriesToken       = consumes<edm::ValueMap<int> >(edm::InputTag("electronMVAValueMapProducer",  "ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories",""));

  vertexToken = consumes <reco::VertexCollection> (edm::InputTag(std::string("offlineSlimmedPrimaryVertices")));
  electronToken = consumes <pat::ElectronCollection> (edm::InputTag(std::string("slimmedElectrons")));
  // BBT, 10-08-18: token for electron ID
  eleIdMapToken_ = consumes<edm::ValueMap<bool> > (edm::InputTag(std::string("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight")));

   muonToken = consumes <pat::MuonCollection> (edm::InputTag(std::string("slimmedMuons"))); // uncomment BBT, 10-12-18
  // muonToken = consumes <pat::MuonCollection> (edm::InputTag(std::string("MuonWithPuppiIsolation")));
  //muonToken = consumes <pat::MuonCollection> (edm::InputTag("PUPPIMuonRelIso", "MuonWithPuppiIsolation","" ) ) ; // comment out ygg core

  //muonToken = consumes <pat::MuonCollection> (edm::InputTag("deterministicSeeds", "muonsWithSeed",""));

  //  muonview_Token = consumes < edm::View<pat::Muon> > (edm::InputTag("PUPPIMuonRelIso", "MuonWithPuppiIsolation","" ) ) ;
  muonview_Token = consumes < edm::View<pat::Muon> > (edm::InputTag(std::string("slimmedMuons")));

  jetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("slimmedJets")));
  //jetTokenWithSeeds = consumes <pat::JetCollection> (edm::InputTag("deterministicSeeds", "jetsWithSeed","")); // BBT 10-12-18
  jetTokenWithSeeds = consumes <pat::JetCollection> (edm::InputTag("deterministicSeeds", "jetsWithSeed","MAOD")); // BBT 10-29-18
  //edm::InputTag jetCollection_v2 = iConfig.getParameter<edm::InputTag>("jetCollection"); // BBT 10-29-18 
  //jetTokenWithSeeds_v2 = consumes <pat::JetCollection> (jetCollection_v2); // BBT 10-29-18 

  //edm::InputTag jetCollection_config; // BBT 10-15-18
  //jetCollection_config = iConfig.getParameter<edm::InputTag>("jetCollection"); // BBT 10-15-18
  //jetTokenWithSeeds_v2 = consumes <pat::JetCollection> (jetCollection_config); // BBT 10-15-18 iConfig.getParameter
  puppijetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("slimmedJetsPuppi")));
  fatjetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("slimmedJetsAK8")));

  //  rerun_fatjetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("selectedPatJetsAK8PFPuppi")));// Rerun of PUPPI ak8
  //rerun_fatjetToken = consumes <pat::JetCollection> (edm::InputTag(std::string("selectedPatJetsCA15PFPuppi")));// Rerun of PUPPI ca15

  //(In moriond17 analysis, met needs to be recalculated.) 
  //   consumes <pat::METCollection> (edm::InputTag("slimmedMETs","","PAT") );

  if( isMC ){
    metToken = consumes <pat::METCollection> (edm::InputTag("slimmedMETs","","PAT") );
  }else{
    //    metToken = consumes <pat::METCollection> (edm::InputTag("slimmedMETsMuEGClean","","") );
    //metToken = consumes <pat::METCollection> (edm::InputTag("slimmedMETs","","RECO") ); // ygg core, but not present in data files?
    metToken = consumes <pat::METCollection> (edm::InputTag("slimmedMETs","","PAT") ); // BBT, 10-24-18
  }

  puppimetToken = consumes <pat::METCollection> (edm::InputTag("slimmedMETsPuppi","","") );

  // topJetsToken    = consumes< boosted::HEPTopJetCollection >(edm::InputTag("HEPTopJetsPFMatcher","heptopjets","p"));
  // subFilterJetsToken = consumes< boosted::SubFilterJetCollection >(edm::InputTag("CA12JetsCA3FilterjetsPFMatcher","subfilterjets","p"));

  packedpfToken = consumes <pat::PackedCandidateCollection> (edm::InputTag(std::string("packedPFCandidates")));

  beamspotToken = consumes <reco::BeamSpot> (edm::InputTag(std::string("offlineBeamSpot")));
  rhoToken = consumes <double> (edm::InputTag(std::string("fixedGridRhoFastjetAll")));
  puInfoToken = consumes <std::vector< PileupSummaryInfo > > (edm::InputTag(std::string("slimmedAddPileupInfo")));

  if( isMC ){
    mcparicleToken = consumes <reco::GenParticleCollection> (edm::InputTag(std::string("prunedGenParticles")));
    genInfoProductToken = consumes <GenEventInfoProduct> (edm::InputTag(std::string("generator")));
    if( ! DISABLE_LHE_FOR_DIBOSON_SAMPLES ){ LHEEventProductToken = consumes<LHEEventProduct> ( edm::InputTag(std::string("externalLHEProducer") )  ); }
    TtGenEventToken = consumes< TtGenEvent >( edm::InputTag("genEvt") );
  }

  // EDMConversionCollectionToken = consumes <reco::ConversionCollection > (edm::InputTag("reducedEgamma","reducedConversions",""));
  EDMBoostedJetsToken     = consumes< boosted::BoostedJetCollection >(edm::InputTag("BoostedJetMatcher","boostedjets","p"));

  if( isMC ){
  genTtbarIdToken_ = consumes<int>( edm::InputTag( "categorizeGenTtbar", "genTtbarId","" ) )  ;
  }

  edm::Service<TFileService> fs_;
  worldTree = fs_->make<TTree>("worldTree", "worldTree");
  eve=0; 
  worldTree->Branch("eve.", "yggdrasilEventVars", &eve, 8000, 1);

  genBHadJetIndexToken_ = consumes< std::vector<int> >( edm::InputTag("matchGenBHadron", "genBHadJetIndex", "") );

  nevents=0;

  std::string era = "2012_53x";
  insample_ = 2500;

  mySample_sampleName_ = "TTJets_MSDecaysCKM_central_Tune4C_13TeV_madgraph_PU20bx25_POSTLS170_V5_v1";
  mySample_xSec_ = 689.1;
  mySample_nGen_ = 25474122;
  intLumi_ = 20000;

  analysisType::analysisType iAnalysisType = analysisType::LJ;

  miniAODhelper.SetUp(era, insample_, iAnalysisType, ! isMC );
  miniAODhelper_Puppi.SetUp(era, insample_, iAnalysisType, ! isMC );
  miniAODhelper_fatjet.SetUp(era, insample_, iAnalysisType, ! isMC );

  muon_roc = new RoccoR ( std::string(  getenv("CMSSW_BASE") ) + "/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/rcdata.2016.v3" );

   // miniAODhelper.SetUpElectronMVA("MiniAOD/MiniAODHelper/data/ElectronMVA/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml","MiniAOD/MiniAODHelper/data/ElectronMVA/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml","MiniAOD/MiniAODHelper/data/ElectronMVA/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml");
  

  toptagger = TopTagger(TopTag::Likelihood, TopTag::CSV, "toplikelihoodtaggerhistos.root");


  // - - - - - - PDF uncertainty - - - - - - - - -
  CT14nlo_PDFSet = new   LHAPDF::PDFSet("CT14nlo");
  CT14nlo_systPDFs = CT14nlo_PDFSet->mkPDFs();


  NNPDF30_nlo_as_0118_PDFSet = new   LHAPDF::PDFSet("NNPDF30_nlo_as_0118");
  NNPDF30_nlo_as_0118_systPDFs = NNPDF30_nlo_as_0118_PDFSet->mkPDFs();


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


  scalefactors.init_all();

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


  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(vertexToken,vtxHandle);
  reco::VertexCollection vtxs = *vtxHandle;

  /// old way of getting electrons
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken,electrons);

  // BBT, 10-08-18: electron ID
  edm::Handle<edm::ValueMap<bool> > ele_id_decisions;
  iEvent.getByToken(eleIdMapToken_ ,ele_id_decisions);

  // //// MVAelectrons
  // edm::Handle< edm::View<pat::Electron> > h_electrons;
  // iEvent.getByToken( EDMElectronsToken,h_electrons );
  // // add electron mva info to electrons
  // edm::Handle<edm::ValueMap<float> > h_mvaValues; 
  // edm::Handle<edm::ValueMap<int> > h_mvaCategories;
  // iEvent.getByToken(EDMeleMVAvaluesToken,h_mvaValues);
  // iEvent.getByToken(EDMeleMVAcategoriesToken,h_mvaCategories);  
  // std::vector<pat::Electron> electrons = miniAODhelper.GetElectronsWithMVAid(h_electrons,h_mvaValues,h_mvaCategories);


  ////
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken,muons);

  edm::Handle<pat::JetCollection> pfjets;
  //iEvent.getByToken(jetTokenWithSeeds_v2,pfjets); // BBT 10-15-18
  if (isMC)
    iEvent.getByToken(jetTokenWithSeeds,pfjets); // BBT 10-12-18
  else
    iEvent.getByToken(jetToken,pfjets); // ygg core

  edm::Handle<pat::JetCollection> pfpuppijets;
  iEvent.getByToken(puppijetToken,pfpuppijets);

  edm::Handle<pat::JetCollection> fatjets;
  iEvent.getByToken(fatjetToken,fatjets);

  //edm::Handle<pat::JetCollection> rerun_fatjets;
  //iEvent.getByToken(rerun_fatjetToken,rerun_fatjets);

  edm::Handle<pat::METCollection> pfmet;
  iEvent.getByToken(metToken,pfmet);

  edm::Handle<pat::METCollection> puppimet;
  iEvent.getByToken(puppimetToken,puppimet);

  edm::Handle<pat::PackedCandidateCollection> packedPFcands;
  iEvent.getByToken(packedpfToken,packedPFcands);


  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(beamspotToken,bsHandle);

  edm::Handle<reco::GenParticleCollection> mcparticles;
  edm::Handle<reco::GenJetCollection> genjetCollection;
  edm::Handle<reco::GenJetCollection> fatgenjetCollection;
  edm::Handle<reco::GenJetCollection> reclustered_fatgenjetCollection ; 
  if( isMC ){
    iEvent.getByToken(mcparicleToken,mcparticles);
    iEvent.getByToken( genJetsToken_ , genjetCollection );
    iEvent.getByToken( fatgenJetsToken_ , fatgenjetCollection );
    iEvent.getByToken( reclusteredfatgenJetsToken_ , reclustered_fatgenjetCollection );
  }

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken,rhoHandle);
  ////------- set up rho for lepton effArea Isolation correction
  double rho_event = ( (rhoHandle.isValid()) ) ? *rhoHandle : -99;
  miniAODhelper.SetRho(rho_event);
  miniAODhelper_Puppi.SetRho(rho_event);
  miniAODhelper_fatjet.SetRho(rho_event);

  edm::Handle<std::vector< PileupSummaryInfo > > PupInfo;
  iEvent.getByToken(puInfoToken,PupInfo);

  // edm::Handle<boosted::HEPTopJetCollection> h_heptopjet;
  // iEvent.getByToken( topJetsToken,h_heptopjet);

  // edm::Handle<boosted::SubFilterJetCollection> h_subfilterjet;
  // iEvent.getByToken( subFilterJetsToken,h_subfilterjet );

  edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
  edm::Handle<LHEEventProduct> LHEEventProductHandle;
  if( isMC ){
  iEvent.getByToken(genInfoProductToken,GenEventInfoHandle);
  if( ! DISABLE_LHE_FOR_DIBOSON_SAMPLES ){ iEvent.getByToken(LHEEventProductToken,  LHEEventProductHandle) ; }
  }

  edm::Handle<boosted::BoostedJetCollection> h_boostedjet;
  iEvent.getByToken( EDMBoostedJetsToken,h_boostedjet);
  
  //  edm::Handle<reco::ConversionCollection> h_conversioncollection;
  // iEvent.getByToken( EDMConversionCollectionToken,h_conversioncollection );

  double GenEventInfoWeight = isMC ? GenEventInfoHandle.product()->weight() : 1.0 ;


  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken, triggerResults);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken( TriggerObjectStandAloneToken , triggerObjects ) ; 
  edm::Handle<edm::TriggerResults> filterResults; // MET Filter, 10-19-18
  iEvent.getByToken(filterResultsToken, filterResults); // MET Filter, 10-19-18

  bool passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1 = false;
  bool passHLT_Ele27_eta2p1_WPTight_Gsf_v = false;

  bool passHLT_Ele27_WPTight_Gsf_v = false;
  

  bool passHLT_IsoMu24_v = false;
  bool passHLT_IsoMu22_v = false;
  bool passHLT_IsoTkMu22_v = false;
  bool passHLT_IsoTkMu24_v = false;
  bool passHLT_IsoMu20_eta2p1_v = false;
  bool passHLT_IsoMu24_eta2p1_v = false;

  bool passHLT_Ele27_WP85_Gsf_v = false;
  bool passHLT_Ele27_eta2p1_WPLoose_Gsf_v = false;
  bool passHLT_Ele27_eta2p1_WP75_Gsf_v = false;

  bool passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v = false;
  bool passHLT_Ele27_eta2p1_WPLoose_Gsf_HT200_v = false;

  bool passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = false;
  bool passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = false;

  bool passHLT_Mu30_TkMu11_v = false;
  bool passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v = false;
  bool passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v = false;
  bool passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = false;
  bool passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = false;
  bool passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = false;
  bool passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = false;

  bool passHLT_Ele25WP60_SC4_Mass55_v = false;

  // bool passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = false ;
  // bool passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = false ;
  // bool passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v = false ;
  bool passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v = false ;
  bool passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v = false ;

  bool passHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v  = false ; 
  bool passHLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v  = false ; 
  bool passHLT_Ele32_WPTight_Gsf_L1DoubleEG_v  = false ; 
  bool passHLT_Ele32_WPTight_Gsf_v  = false ; 
  bool passHLT_Ele35_WPTight_Gsf_v  = false ; 
  bool passHLT_Ele38_WPTight_Gsf_v  = false ; 
  bool passHLT_Ele40_WPTight_Gsf_v  = false ; 
  bool passHLT_IsoMu24_2p1_v  = false ; 
  bool passHLT_IsoMu27_v  = false ; 
  bool passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v  = false ; 
  bool passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v  = false ; 
  bool passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v  = false ; 
  bool passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v  = false ; 
  bool passHLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v  = false ; 
  bool passHLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v  = false ; 
  bool passHLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v  = false ; 
  bool passHLT_PFHT430_SixJet40_BTagCSV_p080_v  = false ; 
  bool passHLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v  = false ; 

  // MET Filters, BBT 02-20-19
  bool passMETFilter_globalTightHalo2016Filter_v = false;
  bool passMETFilter_globalSuperTightHalo2016Filter_v = false;    
  bool passMETFilter_HBHENoiseFilter_v = false;                   
  bool passMETFilter_HBHENoiseIsoFilter_v = false;                
  bool passMETFilter_EcalDeadCellTriggerPrimitiveFilter_v = false;
  bool passMETFilter_BadPFMuonFilter_v = false;                   
  bool passMETFilter_BadChargedCandidateFilter_v = false;               
  bool passMETFilter_ecalBadCalibFilter_v = false;                
  bool passMETFilter_goodVertices_v = false;

  // 2017 MET triggers
  bool passHLT_PFHT500_PFMET100_PFMHT100_IDTight_v = false ;
  bool passHLT_PFHT500_PFMET110_PFMHT110_IDTight_v = false ;
  bool passHLT_PFHT700_PFMET85_PFMHT85_IDTight_v = false ;
  bool passHLT_PFHT700_PFMET95_PFMHT95_IDTight_v = false ;
  bool passHLT_PFHT800_PFMET75_PFMHT75_IDTight_v = false ;
  bool passHLT_PFHT800_PFMET85_PFMHT85_IDTight_v = false ;
  bool passHLT_PFMET110_PFMHT110_IDTight_v = false ;
  bool passHLT_PFMET120_PFMHT120_IDTight_v = false ;
  bool passHLT_PFMET130_PFMHT130_IDTight_v = false ;
  bool passHLT_PFMET140_PFMHT140_IDTight_v = false ;
  bool passHLT_PFMETTypeOne110_PFMHT110_IDTight_v = false ;
  bool passHLT_PFMETTypeOne120_PFMHT120_IDTight_v = false ;
  bool passHLT_PFMETTypeOne130_PFMHT130_IDTight_v = false ;
  bool passHLT_PFMETTypeOne140_PFMHT140_IDTight_v = false ;
  bool passHLT_DiJet110_35_Mjj650_PFMET110_v = false ;
  bool passHLT_DiJet110_35_Mjj650_PFMET120_v = false ;
  bool passHLT_DiJet110_35_Mjj650_PFMET130_v = false ;
  bool passHLT_TripleJet110_35_35_Mjj650_PFMET110_v = false ;
  bool passHLT_TripleJet110_35_35_Mjj650_PFMET120_v = false ;
  bool passHLT_TripleJet110_35_35_Mjj650_PFMET130_v = false ;
  bool passHLT_MET105_IsoTrk50_v = false ;
  bool passHLT_MET120_IsoTrk50_v = false ;
  bool passHLT_PFMET120_PFMHT120_IDTight_PFHT60_v = false ;
  bool passHLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v = false ;
  bool passHLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v = false ;
  bool passHLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v = false ;
  bool passHLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v = false ;
  bool passHLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v = false ;
  bool passHLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v = false ;
  bool passHLT_CaloMET100_HBHECleaned_v = false ;
  bool passHLT_CaloMET250_HBHECleaned_v = false ;
  bool passHLT_CaloMET300_HBHECleaned_v = false ;
  bool passHLT_CaloMET350_HBHECleaned_v = false ;
  bool passHLT_CaloMET70_HBHECleaned_v = false ;
  bool passHLT_CaloMET80_HBHECleaned_v = false ;
  bool passHLT_CaloMET90_HBHECleaned_v = false ;
  bool passHLT_PFMET200_HBHE_BeamHaloCleaned_v = false ;
  bool passHLT_PFMET200_HBHECleaned_v = false ;
  bool passHLT_PFMET250_HBHECleaned_v = false ;
  bool passHLT_PFMET300_HBHECleaned_v = false ;
  bool passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v = false ;
  bool passHLT_PFMET100_PFMHT100_IDTight_PFHT60_v = false ;
  bool passHLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v = false ;
  bool passHLT_PFHT250_v = false ;


  if( triggerResults.isValid() ){
    std::vector<std::string> triggerNames = hlt_config_.triggerNames();

    for( unsigned int iPath=0; iPath<triggerNames.size(); iPath++ ){
      std::string pathName = triggerNames[iPath];
      unsigned int hltIndex = hlt_config_.triggerIndex(pathName);
      
      if( hltIndex >= triggerResults->size() ) continue;
      int accept = triggerResults->accept(hltIndex);

      if( accept ){
	const unsigned long MatchedAtTheBegining = 0 ; 

	if( pathName.find( "HLT_IsoMu24_v"        ,0) == MatchedAtTheBegining ){ passHLT_IsoMu24_v = true;}
	if( pathName.find( "HLT_IsoTkMu24_v"      ,0) == MatchedAtTheBegining ){ passHLT_IsoTkMu24_v = true;}
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
	if( pathName.find( "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",0) == MatchedAtTheBegining ){ passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = true;}
	if( pathName.find( "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v" ,0) == MatchedAtTheBegining ){ passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = true;}

	if( pathName.find( "HLT_Ele25WP60_SC4_Mass55_v"                       ,0) == MatchedAtTheBegining ){ passHLT_Ele25WP60_SC4_Mass55_v = true;}

	// if( pathName.find( "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"      ,0) == MatchedAtTheBegining ){ passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v        = true ; }
	// if( pathName.find( "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",0) == MatchedAtTheBegining ){ passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v  = true ; }
	// if( pathName.find( "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v" ,0) == MatchedAtTheBegining ){ passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v   = true ; }
	if( pathName.find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"            ,0) == MatchedAtTheBegining ){ passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v	        = true ; }
	if( pathName.find( "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"          ,0) == MatchedAtTheBegining ){ passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v            = true ; }
	

	if( pathName.find( "HLT_Ele27_eta2p1_WPTight_Gsf_v"        ,0) == MatchedAtTheBegining ){ passHLT_Ele27_eta2p1_WPTight_Gsf_v = true;}
	if( pathName.find( "HLT_Ele27_WPTight_Gsf_v"        ,0) == MatchedAtTheBegining ){ passHLT_Ele27_WPTight_Gsf_v = true;}
	if( pathName.find( "HLT_IsoMu22_v"        ,0) == MatchedAtTheBegining ){ passHLT_IsoMu22_v = true;}
	if( pathName.find( "HLT_IsoTkMu22_v"        ,0) == MatchedAtTheBegining ){ passHLT_IsoTkMu22_v = true;}
	if( pathName.find( "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"        ,0) == MatchedAtTheBegining ){ passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = true;}


	if( pathName.find( "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v = true ;} 
	if( pathName.find( "HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v = true ;} 
	if( pathName.find( "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Ele32_WPTight_Gsf_L1DoubleEG_v = true ;} 
	if( pathName.find( "HLT_Ele32_WPTight_Gsf_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Ele32_WPTight_Gsf_v = true ;} 
	if( pathName.find( "HLT_Ele35_WPTight_Gsf_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Ele35_WPTight_Gsf_v = true ;} 
	if( pathName.find( "HLT_Ele38_WPTight_Gsf_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Ele38_WPTight_Gsf_v = true ;} 
	if( pathName.find( "HLT_Ele40_WPTight_Gsf_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Ele40_WPTight_Gsf_v = true ;} 
	if( pathName.find( "HLT_IsoMu24_2p1_v" , 0 ) == MatchedAtTheBegining ){ passHLT_IsoMu24_2p1_v = true ;} 
	if( pathName.find( "HLT_IsoMu27_v" , 0 ) == MatchedAtTheBegining ){ passHLT_IsoMu27_v = true ;} 
	if( pathName.find( "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = true ;} 
	if( pathName.find( "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = true ;} 
	if( pathName.find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v = true ;} 
	if( pathName.find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v" , 0 ) == MatchedAtTheBegining ){ passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v = true ;} 
	if( pathName.find( "HLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v = true ;} 
	if( pathName.find( "HLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v = true ;} 
	if( pathName.find( "HLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v = true ;} 
	if( pathName.find( "HLT_PFHT430_SixJet40_BTagCSV_p080_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT430_SixJet40_BTagCSV_p080_v = true ;} 
	if( pathName.find( "HLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v = true ;}

	// MET Filters, BBT 11-06-18
	// 2017 MET triggers
	if( pathName.find( "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT500_PFMET100_PFMHT100_IDTight_v = true ;}
	if( pathName.find( "HLT_PFHT500_PFMET110_PFMHT110_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT500_PFMET110_PFMHT110_IDTight_v = true ;}
	if( pathName.find( "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT700_PFMET85_PFMHT85_IDTight_v = true ;}
	if( pathName.find( "HLT_PFHT700_PFMET95_PFMHT95_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT700_PFMET95_PFMHT95_IDTight_v = true ;}
	if( pathName.find( "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT800_PFMET75_PFMHT75_IDTight_v = true ;}
	if( pathName.find( "HLT_PFHT800_PFMET85_PFMHT85_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT800_PFMET85_PFMHT85_IDTight_v = true ;}
	if( pathName.find( "HLT_PFMET110_PFMHT110_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET110_PFMHT110_IDTight_v = true ;}
	if( pathName.find( "HLT_PFMET120_PFMHT120_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET120_PFMHT120_IDTight_v = true ;}
	if( pathName.find( "HLT_PFMET130_PFMHT130_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET130_PFMHT130_IDTight_v = true ;}
	if( pathName.find( "HLT_PFMET140_PFMHT140_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET140_PFMHT140_IDTight_v = true ;}
	if( pathName.find( "HLT_PFMETTypeOne110_PFMHT110_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMETTypeOne110_PFMHT110_IDTight_v = true ;}
	if( pathName.find( "HLT_PFMETTypeOne120_PFMHT120_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMETTypeOne120_PFMHT120_IDTight_v = true ;}
	if( pathName.find( "HLT_PFMETTypeOne130_PFMHT130_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMETTypeOne130_PFMHT130_IDTight_v = true ;}
	if( pathName.find( "HLT_PFMETTypeOne140_PFMHT140_IDTight_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMETTypeOne140_PFMHT140_IDTight_v = true ;}
	if( pathName.find( "HLT_DiJet110_35_Mjj650_PFMET110_v" , 0 ) == MatchedAtTheBegining ){ passHLT_DiJet110_35_Mjj650_PFMET110_v = true ;}
	if( pathName.find( "HLT_DiJet110_35_Mjj650_PFMET120_v" , 0 ) == MatchedAtTheBegining ){ passHLT_DiJet110_35_Mjj650_PFMET120_v = true ;}
	if( pathName.find( "HLT_DiJet110_35_Mjj650_PFMET130_v" , 0 ) == MatchedAtTheBegining ){ passHLT_DiJet110_35_Mjj650_PFMET130_v = true ;}
	if( pathName.find( "HLT_TripleJet110_35_35_Mjj650_PFMET110_v" , 0 ) == MatchedAtTheBegining ){ passHLT_TripleJet110_35_35_Mjj650_PFMET110_v = true ;}
	if( pathName.find( "HLT_TripleJet110_35_35_Mjj650_PFMET120_v" , 0 ) == MatchedAtTheBegining ){ passHLT_TripleJet110_35_35_Mjj650_PFMET120_v = true ;}
	if( pathName.find( "HLT_TripleJet110_35_35_Mjj650_PFMET130_v" , 0 ) == MatchedAtTheBegining ){ passHLT_TripleJet110_35_35_Mjj650_PFMET130_v = true ;}
	if( pathName.find( "HLT_MET105_IsoTrk50_v" , 0 ) == MatchedAtTheBegining ){ passHLT_MET105_IsoTrk50_v = true ;}
	if( pathName.find( "HLT_MET120_IsoTrk50_v" , 0 ) == MatchedAtTheBegining ){ passHLT_MET120_IsoTrk50_v = true ;}
	if( pathName.find( "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET120_PFMHT120_IDTight_PFHT60_v = true ;}
	if( pathName.find( "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v = true ;}
	if( pathName.find( "HLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v = true ;}
	if( pathName.find( "HLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v = true ;}
	if( pathName.find( "HLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v = true ;}
	if( pathName.find( "HLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v = true ;}
	if( pathName.find( "HLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v = true ;}
	if( pathName.find( "HLT_CaloMET100_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_CaloMET100_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_CaloMET250_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_CaloMET250_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_CaloMET300_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_CaloMET300_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_CaloMET350_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_CaloMET350_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_CaloMET70_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_CaloMET70_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_CaloMET80_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_CaloMET80_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_CaloMET90_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_CaloMET90_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_PFMET200_HBHE_BeamHaloCleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET200_HBHE_BeamHaloCleaned_v = true ;}
	if( pathName.find( "HLT_PFMET200_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET200_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_PFMET250_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET250_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_PFMET300_HBHECleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET300_HBHECleaned_v = true ;}
	if( pathName.find( "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v = true ;}
	if( pathName.find( "HLT_PFMET100_PFMHT100_IDTight_PFHT60_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMET100_PFMHT100_IDTight_PFHT60_v = true ;}
	if( pathName.find( "HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v = true ;}
	if( pathName.find( "HLT_PFHT250_v" , 0 ) == MatchedAtTheBegining ){ passHLT_PFHT250_v = true ;}

      }
    }
  }
  
  // MET Filters, BBT 02-20-19
  if ( filterResults.isValid() ){
    const edm::TriggerNames& filterNames = iEvent.triggerNames(*filterResults);
    unsigned int nFilters = filterNames.triggerNames().size();

    unsigned int passMETFilter_globalTightHalo2016FilterIndex          = filterNames.triggerIndex("Flag_globalTightHalo2016Filter");
    unsigned int passMETFilter_globalSuperTightHalo2016FilterIndex     = filterNames.triggerIndex("Flag_globalSuperTightHalo2016Filter");
    unsigned int passMETFilter_HBHENoiseFilterIndex                    = filterNames.triggerIndex("Flag_HBHENoiseFilter");
    unsigned int passMETFilter_HBHENoiseIsoFilterIndex                 = filterNames.triggerIndex("Flag_HBHENoiseIsoFilter");
    unsigned int passMETFilter_EcalDeadCellTriggerPrimitiveFilterIndex = filterNames.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter");
    unsigned int passMETFilter_BadPFMuonFilterIndex                    = filterNames.triggerIndex("Flag_BadPFMuonFilter"); 
    unsigned int passMETFilter_BadChargedCandidateFilterIndex          = filterNames.triggerIndex("Flag_BadChargedCandidateFilter"); 
    unsigned int passMETFilter_ecalBadCalibFilterIndex                 = filterNames.triggerIndex("Flag_ecalBadCalibFilter");
    unsigned int passMETFilter_goodVerticesIndex                       = filterNames.triggerIndex("Flag_goodVertices");


    // MET Filters, BBT 02-20-19
    passMETFilter_globalTightHalo2016Filter_v = passMETFilter_globalTightHalo2016FilterIndex < nFilters ? filterResults->accept(passMETFilter_globalTightHalo2016FilterIndex) : false;
    passMETFilter_globalSuperTightHalo2016Filter_v = passMETFilter_globalSuperTightHalo2016FilterIndex < nFilters ? filterResults->accept(passMETFilter_globalSuperTightHalo2016FilterIndex) : false;
    passMETFilter_HBHENoiseFilter_v = passMETFilter_HBHENoiseFilterIndex < nFilters ? filterResults->accept(passMETFilter_HBHENoiseFilterIndex) : false;     
    passMETFilter_HBHENoiseIsoFilter_v = passMETFilter_HBHENoiseIsoFilterIndex < nFilters ? filterResults->accept(passMETFilter_HBHENoiseIsoFilterIndex) : false;
    passMETFilter_EcalDeadCellTriggerPrimitiveFilter_v = passMETFilter_EcalDeadCellTriggerPrimitiveFilterIndex < nFilters ? filterResults->accept(passMETFilter_EcalDeadCellTriggerPrimitiveFilterIndex) : false;
    passMETFilter_BadPFMuonFilter_v = passMETFilter_BadPFMuonFilterIndex < nFilters ? filterResults->accept(passMETFilter_BadPFMuonFilterIndex) : false;
    passMETFilter_BadChargedCandidateFilter_v = passMETFilter_BadChargedCandidateFilterIndex < nFilters ? filterResults->accept(passMETFilter_BadChargedCandidateFilterIndex) : false;
    passMETFilter_ecalBadCalibFilter_v = passMETFilter_ecalBadCalibFilterIndex < nFilters ? filterResults->accept(passMETFilter_ecalBadCalibFilterIndex) : false;
    passMETFilter_goodVertices_v = passMETFilter_goodVerticesIndex < nFilters ? filterResults->accept(passMETFilter_goodVerticesIndex) : false;
  }



  ///-- Genjet Information 
  vdouble genjet_pt;
  vdouble genjet_eta;
  vdouble genjet_phi;
  vdouble genjet_m;
  vint   genjet_BhadronMatch;

  vdouble fatgenjet_pt;
  vdouble fatgenjet_eta;
  vdouble fatgenjet_phi;
  vdouble fatgenjet_m;

  if( isMC  ){

    //    genjets
    
    edm::Handle<std::vector<int> > genBHadJetIndex;
    iEvent.getByToken(genBHadJetIndexToken_, genBHadJetIndex);
    
    const std::vector<reco::GenJet> * genjets  = genjetCollection.product();

    for( unsigned int iGen = 0 ; iGen < genjets->size() ; iGen ++){
	reco::GenJet jet = genjets->at( iGen );

	if( jet . pt() < 20.0 ) continue ;
	if( fabs( jet . eta() ) > 5 ) continue ;

	genjet_pt  . push_back( jet . pt() );
	genjet_phi . push_back( jet . phi() );
	genjet_eta . push_back( jet . eta() );
	genjet_m   . push_back( jet . mass () );
	
	bool b_associateWithBHadron = false;
	for( unsigned int i = 0 ; i < genBHadJetIndex->size() && ! b_associateWithBHadron ; i++ ){
	  if( iGen == (unsigned int) genBHadJetIndex->at(i) ){ b_associateWithBHadron = true ; }
	}
	genjet_BhadronMatch . push_back( b_associateWithBHadron ? 1 : 0 );

      }
	   

    // memo : https://github.com/hsatoshi/MiniAODPrivateTuple/blob/master/TupleMaker/plugins/TupleMaker.cc#L3303
    
    
    
    {// ---------- AK 8 gen jets
      
      const std::vector<reco::GenJet> * genjets  = fatgenjetCollection.product();
      
      for( unsigned int iGen = 0 ; iGen < genjets->size() ; iGen ++){
	reco::GenJet jet = genjets->at( iGen );
	
	if( jet . pt() < 150.0 ) continue ;
	if( fabs( jet . eta() ) > 5 ) continue ;
	
	fatgenjet_pt  . push_back( jet . pt() );
	fatgenjet_phi . push_back( jet . phi() );
	fatgenjet_eta . push_back( jet . eta() );
	fatgenjet_m   . push_back( jet . mass () );
	
      }
    }


  }


  ////----------------------
  ////---- tt+X Categorization
  ////----------------------
  if( isMC ){
  edm::Handle<int> genTtbarId;
  iEvent.getByToken(genTtbarIdToken_, genTtbarId);

  // Fill ID including information about b jets from top in acceptance
  h_ttbarId_->Fill(*genTtbarId);
   
  // Fill ID based only on additional b/c jets
  h_ttbarAdditionalJetId_->Fill(*genTtbarId%100);
  
  eve->additionalJetEventId_ = *genTtbarId;
  
  if( ! DISABLE_LHE_FOR_DIBOSON_SAMPLES ){
    const int idx_Q2_upup     = 1005;
    eve->weight_q2_upup_     = LHEEventProductHandle -> weights()[idx_Q2_upup]    .wgt / LHEEventProductHandle -> originalXWGTUP(); 
    const int idx_Q2_downdown = 1009;
    eve->weight_q2_downdown_ = LHEEventProductHandle -> weights()[idx_Q2_downdown].wgt / LHEEventProductHandle -> originalXWGTUP(); 
  }
    //

//(debug)    std::cout <<"Satoshi debug : wgt size " << LHEEventProductHandle -> weights().size() << std::endl ; 
//(debug)    std::cout <<"Satoshi debug : original wgt = " << LHEEventProductHandle -> originalXWGTUP()  << std::endl ; 
//(debug)    std::cout <<"Satoshi debug : wgt  [zero]  = " << LHEEventProductHandle -> weights()[0].wgt  << std::endl ; 

    std::vector<int> mcWeight_key ; 


    if( ! DISABLE_LHE_FOR_DIBOSON_SAMPLES && LHEEventProductHandle -> weights().size() > 10 ){

    // Memo : The size of "weights()" was 1080 as of 94x, Fall17 MiniAID.
    // Memo : It means that the index of the array is ID-1000.

     mcWeight_key . push_back ( 0 ); 
     mcWeight_key . push_back ( 1 );
     mcWeight_key . push_back ( 2 );
     mcWeight_key . push_back ( 3 );
     mcWeight_key . push_back ( 4 );
     mcWeight_key . push_back ( 5 );
     mcWeight_key . push_back ( 6 );
     mcWeight_key . push_back ( 7 );
     mcWeight_key . push_back ( 8 );
     mcWeight_key . push_back ( 9 );
    for( std::vector<int>::iterator index =  mcWeight_key . begin(); 
	 index !=  mcWeight_key . end(); 
	 index ++ ){
      eve -> mcWeight_value. push_back(  LHEEventProductHandle -> weights()[ *index ].wgt  );
    }

     mcWeight_key . push_back ( 2000 ); // my ID. use to store Original XWGTUP();
    eve -> mcWeight_value. push_back( LHEEventProductHandle -> originalXWGTUP() );

    }

    if( GenEventInfoHandle -> weights().size() >=  14 ){

     mcWeight_key . push_back ( 3000 ); 
     mcWeight_key . push_back ( 3001 ); 
     mcWeight_key . push_back ( 3002 ); 
     mcWeight_key . push_back ( 3003 ); 
     mcWeight_key . push_back ( 3004 ); 
     mcWeight_key . push_back ( 3005 ); 
     mcWeight_key . push_back ( 3006 ); 
     mcWeight_key . push_back ( 3007 );
     mcWeight_key . push_back ( 3008 ); 
     mcWeight_key . push_back ( 3009 ); 
     mcWeight_key . push_back ( 3010 ); 
     mcWeight_key . push_back ( 3011 ); 
     mcWeight_key . push_back ( 3012 ); 
     mcWeight_key . push_back ( 3013 ); 

    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[0]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[1]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[2]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[3]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[4]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[5]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[6]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[7]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[8]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[9]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[10]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[11]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[12]  );
    eve -> mcWeight_value. push_back(  GenEventInfoHandle -> weights()[13]  );

    }

    auto pdfInfos = GenEventInfoHandle -> pdf();
    double pdfNominal = pdfInfos->xPDF.first * pdfInfos->xPDF.second;


    { /// PDF uncertainty for CT14
    std::vector<double> pdfs;
    for (size_t j = 0; j < CT14nlo_systPDFs.size(); ++j) {
      double xpdf1 = CT14nlo_systPDFs[j]->xfxQ(pdfInfos->id.first,  pdfInfos->x.first,  pdfInfos->scalePDF);
      double xpdf2 = CT14nlo_systPDFs[j]->xfxQ(pdfInfos->id.second, pdfInfos->x.second, pdfInfos->scalePDF);
      pdfs.push_back(xpdf1 * xpdf2);
    }
    //  Combine the products and compute the 1 sigma shift 
    const LHAPDF::PDFUncertainty pdfUnc = CT14nlo_PDFSet -> uncertainty(pdfs, 68.);

    //  Calculate the up/down weights
    if (std::isfinite(1./pdfNominal)) {
      eve-> weight_PDF_CT14nlo_up_   = (pdfUnc.central + pdfUnc.errplus)  / pdfNominal;
      eve-> weight_PDF_CT14nlo_down_ = (pdfUnc.central - pdfUnc.errminus) / pdfNominal;
    }else{
      eve-> weight_PDF_CT14nlo_up_   = 1.0;
      eve-> weight_PDF_CT14nlo_down_ = 1.0;
    }
    }

    { // Syst for NNPDF 
    std::vector<double> pdfs;
    for (size_t j = 0; j < NNPDF30_nlo_as_0118_systPDFs.size(); ++j) {
      double xpdf1 = NNPDF30_nlo_as_0118_systPDFs[j]->xfxQ(pdfInfos->id.first,  pdfInfos->x.first,  pdfInfos->scalePDF);
      double xpdf2 = NNPDF30_nlo_as_0118_systPDFs[j]->xfxQ(pdfInfos->id.second, pdfInfos->x.second, pdfInfos->scalePDF);
      pdfs.push_back(xpdf1 * xpdf2);
    }

    //  Combine the products and compute the 1 sigma shift 
    const LHAPDF::PDFUncertainty pdfUnc = NNPDF30_nlo_as_0118_PDFSet -> uncertainty(pdfs, 68.);

    //  Calculate the up/down weights
    if (std::isfinite(1./pdfNominal)) {
      eve-> weight_PDF_NNPDF30NLO_up_   = (pdfUnc.central + pdfUnc.errplus)  / pdfNominal;
      eve-> weight_PDF_NNPDF30NLO_down_ = (pdfUnc.central - pdfUnc.errminus) / pdfNominal;
    }else{
      eve-> weight_PDF_NNPDF30NLO_up_   = 1.0;
      eve-> weight_PDF_NNPDF30NLO_down_ = 1.0;
    }
    }

  }else{
    eve->additionalJetEventId_ = -10 ; 
    eve->weight_q2_upup_     = 1.0 ;
    eve->weight_q2_downdown_ = 1.0 ;

    eve-> weight_PDF_CT14nlo_up_   = 1.0;
    eve-> weight_PDF_CT14nlo_down_ = 1.0;
    
    eve-> weight_PDF_NNPDF30NLO_up_   = 1.0;
    eve-> weight_PDF_NNPDF30NLO_down_ = 1.0;
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
		      (fabs(vtx->z()) <= 24.0) &&
		      (fabs(vtx->position().Rho()) <= 2.0) 
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

  if( numpv>0 ){
    miniAODhelper.SetVertex(vertex);
    miniAODhelper_Puppi.SetVertex(vertex);
    miniAODhelper_fatjet.SetVertex(vertex);
  }

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


  edm::Handle<reco::JetCorrector> corrector ; 
  iEvent.getByToken(jetCorrectorToken_, corrector );
  miniAODhelper.SetJetCorrector( &(*corrector) );

  edm::Handle<reco::JetCorrector> puppi_corrector ; 
  iEvent.getByToken(puppijetCorrectorToken_, puppi_corrector );
  miniAODhelper_Puppi.SetJetCorrector( &(*puppi_corrector) );

  edm::Handle<reco::JetCorrector> fatjet_corrector ; 
  iEvent.getByToken( fatjetCorrectorToken_, fatjet_corrector );
  miniAODhelper_fatjet.SetBoostedJetCorrector( &(*fatjet_corrector) ); // set ak8 corrector.


  int mHdecay = -1;
  mHdecay = isMC ? miniAODhelper.GetHiggsDecay(mcparticles) : -1 ;
  eve->higgsDecayType_=mHdecay;

  TLorentzVector GenTopQuark, GenAntitopQuark;
  eve->ttbarDecayType_ = isMC ? miniAODhelper.GetTTbarDecay(mcparticles , & GenTopQuark , & GenAntitopQuark ) : -10 ;

  if( ! isTTbarMC ){
    eve -> weight_topPt_ = 1.0;
  }else{
    
    edm::Handle<TtGenEvent> genEvt ;
    iEvent.getByToken( TtGenEventToken, genEvt );

    // **** Ignore the genEvt since I can not pass the packed (in config.py).
//    eve -> weight_topPt_ = sqrt( 
//				exp( 0.156 -0.00137 * genEvt->leptonicDecayTop()->pt() )
//				*
//				exp( 0.156 -0.00137 * genEvt->hadronicDecayTop()->pt() )
//				 );


    eve -> weight_topPt_ = sqrt( 
				exp( 0.156 -0.00137 * GenTopQuark . Pt() )
				*
				exp( 0.156 -0.00137 * GenAntitopQuark . Pt() )
				 );

    // parameters taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting?rev=19
  }


  if( isMC ){

      MiniAODHelper::_topquarkdecayobjects topPosDecay = { }; 
      MiniAODHelper::_topquarkdecayobjects topNegDecay = { }; 

      std::vector< const reco::Candidate * > registeredTop ;
      std::vector< const reco::Candidate * > list ;
      vint parentIdx ; 

//      if( isTTbarMC ){
//      // gather information of top quark decay in ttbar.
//
//      for(size_t i=0; i<mcparticles->size();i++){
//    
//	if( abs( (*mcparticles)[i].pdgId()  ) == 6 ){
//
//	  const reco::Candidate * cand =  & (*mcparticles)[i] ; 
//	  cand = miniAODhelper.GetObjectJustBeforeDecay ( cand );
//
//	  if(std::find( registeredTop.begin(), registeredTop.end(), cand ) != registeredTop.end()) {
//	    // already registered. skip.
//	    continue ; 
//	  }
//
//	  std::cout <<"top quark (pdg id = " <<  cand -> pdgId() << ") is being registerd." << std::endl ; 
//      
//	  registeredTop.push_back(cand);
//	  miniAODhelper.FillTopQuarkDecayInfomration ( cand ,
//						       ( cand -> pdgId() == 6 ) ? ( & topPosDecay ) : ( & topNegDecay ) ) ; 
//	  
//      
//	} // end if : |PDGID|==6
//      }// end MC particle loop 
//
////  struct _topquarkdecayobjects {
////    const reco::Candidate * top ;
////    const reco::Candidate * bottom ;
////    const reco::Candidate * W ;
////    const reco::Candidate * WChild_up;
////    const reco::Candidate * WChild_down;
////    bool isWChild_tau ;
////    const reco::Candidate * Tau_Neu ;
////    std::vector< const reco::Candidate *> TauChildren ;
////
////    bool isLeptonicDecay(){
////      return
////	abs( WChild_down->pdgId() ) == 11
////	||
////	abs( WChild_down->pdgId() ) == 13
////	||
////	isWChild_tau ;
////    }
//
//      int idx = 0 ; 
//      list.push_back(  topPosDecay.top );         parentIdx.push_back( -1 ) ; idx = parentIdx.size()-1;
//      list.push_back(  topPosDecay.bottom );      parentIdx.push_back( idx );
//      list.push_back(  topPosDecay.W );           parentIdx.push_back( idx ); idx = parentIdx.size()-1;
//      list.push_back(  topPosDecay.WChild_down ); parentIdx.push_back( idx ); idx = parentIdx.size()-1; // this particle can be tau.
//      if( topPosDecay.isWChild_tau ){
//	list .push_back( topPosDecay.Tau_Neu ) ;  parentIdx.push_back( idx );      
//	for( std::vector< const reco::Candidate *>::iterator it = topPosDecay.TauChildren.begin() ;
//	     it != topPosDecay.TauChildren.end() ;
//	     it ++
//	     ){
//	  list.push_back( (*it) );  parentIdx.push_back( idx ); 
//	}
//      }
//      // anti-top quark
//      list.push_back(  topNegDecay.top );         parentIdx.push_back( -1 ) ; idx = parentIdx.size()-1;
//      list.push_back(  topNegDecay.bottom );      parentIdx.push_back( idx );
//      list.push_back(  topNegDecay.W );           parentIdx.push_back( idx ); idx = parentIdx.size()-1;
//      list.push_back(  topNegDecay.WChild_up  );  parentIdx.push_back( idx );
//      list.push_back(  topNegDecay.WChild_down ); parentIdx.push_back( idx ); idx = parentIdx.size()-1; // this particle can be tau.
//      if( topNegDecay.isWChild_tau ){
//	list .push_back( topNegDecay.Tau_Neu ) ;  parentIdx.push_back( idx );      
//	for( std::vector< const reco::Candidate *>::iterator it = topNegDecay.TauChildren.begin() ;
//	     it != topNegDecay.TauChildren.end() ;
//	     it ++
//	     ){
//	  list.push_back( (*it) );  parentIdx.push_back( idx ); 
//	}
//      }
//
//      } // end if "isTTbar MC"


      std::vector<const reco::Candidate * > idx_higgs ; 
      std::vector<const reco::Candidate * > idx_Z ; 
      std::vector<const reco::Candidate * > idx_W ; 
      std::vector<const reco::Candidate * > idx_singletop ; 
      std::vector<const reco::Candidate * > idx_rescured_zw ;  
      for(size_t i=0; i < mcparticles ->size();i++){ 
	
	if ( abs( (*mcparticles)[i].pdgId() )  == 6 ){
	  const reco::Candidate * cand =  & (*mcparticles)[i] ; 
	  cand = miniAODhelper.GetObjectJustBeforeDecay ( cand );
	  if ( ! checkIfRegisterd( cand , idx_singletop ) ){
	    idx_singletop . push_back( cand );
	  }
	}

	// - - - Check Higgs - - - - 
	if ( (*mcparticles)[i].pdgId()  == 25 ){
	  const reco::Candidate * cand =  & (*mcparticles)[i] ; 
	  cand = miniAODhelper.GetObjectJustBeforeDecay ( cand );
	  if ( ! checkIfRegisterd( cand , idx_higgs ) ){
	    idx_higgs . push_back( cand );
	  }
	}
       
	// - - - Check Z boson - - - - 
	if ( (*mcparticles)[i].pdgId()  == 23 ){
	  const reco::Candidate * cand =  & (*mcparticles)[i] ; 
	  cand = miniAODhelper.GetObjectJustBeforeDecay ( cand );
	  if ( ! checkIfRegisterd( cand , idx_Z ) ){

	    const reco::Candidate * cand_afterbirth = TraceBackToJustAfterBirth( & (*mcparticles)[i] ) ;
	    bool isZbosonFromHiggs = false ;
	    for ( unsigned int i = 0 ; i <  cand_afterbirth -> numberOfMothers(); i++ ){
	      if(       cand_afterbirth -> mother( i ) -> pdgId()  ==  25 ){
		isZbosonFromHiggs = true ; 
	      }
	    }
	    if( ! isZbosonFromHiggs ){
	      if( ! checkIfRegisterd( cand , idx_Z ) ){ idx_Z . push_back( cand ) ; }
	    }
	  }
	}//// end Z boson

	///  =----- check W
	if ( abs( (*mcparticles)[i].pdgId() )  == 24 ){
	  const reco::Candidate * cand =  & (*mcparticles)[i] ; 
	  cand = miniAODhelper.GetObjectJustBeforeDecay ( cand );
	  if ( ! checkIfRegisterd( cand , idx_W ) ){

	    // - check its parent, and requore not top/H;.
	    const reco::Candidate * cand_afterbirth = TraceBackToJustAfterBirth( cand ) ;

	    bool isWbosonFromTopOrHiggs = false ; 
	    for ( unsigned int i = 0 ; i <  cand_afterbirth -> numberOfMothers(); i++ ){
	      if(       cand_afterbirth -> mother( i ) -> pdgId()  ==  25 ||
			fabs( cand_afterbirth -> mother( i ) -> pdgId() ) ==  6 ||
			fabs( cand_afterbirth -> mother( i ) -> pdgId() ) == 15 ){
		isWbosonFromTopOrHiggs = true ; 
	      }
	    }
	    if( ! isWbosonFromTopOrHiggs ){
	      if( ! checkIfRegisterd( cand , idx_W ) ){ idx_W . push_back( cand ); }
	    }
	  }
	}//// end W boson

	// Tricky Z boson (and also W.)
	//  - I found non negligible events where chaged leptons emrge skipping Zboson(pdgid 23) but from partons.
	//  - to catch such leptons as prompt, check if the gen particle has two leptons in pair (consistent with W/Z assumption))
	//    Require no hadron particles in that decay, allow photons and quarks.
	{
	  if( abs( (*mcparticles)[i].pdgId() )  != 24 && // not w
	      abs( (*mcparticles)[i].pdgId() )  != 23 && // not z
	      abs( (*mcparticles)[i].pdgId() )  != 25 && // not H
	      abs( (*mcparticles)[i].pdgId() )  != 15 && // not tau
	      abs( (*mcparticles)[i].pdgId() )  < 100 && // not hadron
	      (*mcparticles)[i].numberOfDaughters() >=2 ){
	    const reco::Candidate * posPDGID_lepton = 0 ;
	    const reco::Candidate * negPDGID_lepton = 0 ;
	    bool noExtraHandrons = true; 
	    for ( unsigned int iDaug = 0 ; iDaug < (*mcparticles)[i].numberOfDaughters() ; iDaug++ ){

	      int pdgid = (*mcparticles)[i].daughter( iDaug ) -> pdgId() ;

	      if(  abs( pdgid ) >= 100 ){ noExtraHandrons = false ; break ; }

	      if(  11 <= pdgid && pdgid <=  16 ) posPDGID_lepton = TraceBackToJustAfterBirth( ( (*mcparticles)[i].daughter( iDaug ) ) );
	      if( -11 >= pdgid && pdgid >= -16 ) negPDGID_lepton = TraceBackToJustAfterBirth( ( (*mcparticles)[i].daughter( iDaug ) ) );
	    }
	    
	    if( noExtraHandrons && posPDGID_lepton != 0 && negPDGID_lepton != 0 
		&& ! checkIfRegisterd( posPDGID_lepton , idx_rescured_zw )
		&& ! checkIfRegisterd( negPDGID_lepton , idx_rescured_zw )
		){
	      // require minimum energy to cut off noise...
	      if( posPDGID_lepton -> pt() > 10.0 && fabs( posPDGID_lepton -> eta() ) < 3.0 ){ 
		idx_rescured_zw . push_back( posPDGID_lepton ); 
	      }
	      if( negPDGID_lepton -> pt() > 10.0 && fabs( negPDGID_lepton -> eta() ) < 3.0 ){ 
		idx_rescured_zw . push_back( negPDGID_lepton );
	      }
	    }
	  }
	}


      }// end of W/H/Z search.

      // debug 
      if( false ){
	for(size_t i=0; i < mcparticles ->size();i++){ 
	  
	  if( (*mcparticles)[i]. numberOfMothers() == 0 ){
	    DumpDecay( &((*mcparticles)[i]) , 0 );
	  }
	}//debug
      }


      // - - - 
      // Store the information of gen particles in interest.
      // - - - -


      // idx_singletop
      for( unsigned int iTop  = 0 ; iTop < idx_singletop . size ()  ; iTop ++ ){

	MiniAODHelper::_topquarkdecayobjects top = { };
	miniAODhelper.FillTopQuarkDecayInfomration (  idx_singletop[iTop] , & top ) ;

	int idx = -1 ; 
	list.push_back(  top.top );         parentIdx.push_back( -1 ) ; idx = parentIdx.size()-1;
	list.push_back(  top.bottom );      parentIdx.push_back( idx );
	list.push_back(  top.W );           parentIdx.push_back( idx ); idx = parentIdx.size()-1;
	list.push_back(  top.WChild_up  );  parentIdx.push_back( idx );
	list.push_back(  top.WChild_down ); parentIdx.push_back( idx ); idx = parentIdx.size()-1; // this particle can be tau.
	if( top.isWChild_tau ){
	  list .push_back( top.Tau_Neu ) ;  parentIdx.push_back( idx );      
	  for( std::vector< const reco::Candidate *>::iterator it = top.TauChildren.begin() ;
	       it != top.TauChildren.end() ;
	       it ++
	       ){
	    list.push_back( (*it) );  parentIdx.push_back( idx ); 
	  }
	} // tau decay

      } // top loop 




      std::vector<const reco::Candidate * > idx_sum;
      for( std::vector<const reco::Candidate * >::iterator i = idx_higgs.begin() ; i != idx_higgs.end(); i++ ){ idx_sum . push_back( * i );}
      for( std::vector<const reco::Candidate * >::iterator i = idx_W    .begin() ; i != idx_W    .end(); i++ ){ idx_sum . push_back( * i );}
      for( std::vector<const reco::Candidate * >::iterator i = idx_Z    .begin() ; i != idx_Z    .end(); i++ ){ idx_sum . push_back( * i );}


      std::vector<const reco::Candidate * > idx_tau ;  // keep tau later.
      std::vector< int  > idx_tausParent;  

      for( unsigned int iParent  = 0 ; iParent < idx_sum . size ()  ; iParent ++ ){

	list.push_back( idx_sum[iParent] );
	parentIdx.push_back( -1 ) ; // -1 = no parent is assigned
	int ThisParentIdx = parentIdx.size()-1;

	for( unsigned int iChild  = 0 ; iChild < idx_sum[iParent]-> numberOfDaughters()  ; iChild ++ ){
	  if( abs( idx_sum[iParent]->daughter( iChild ) ->pdgId() ) == 15 ){
	    idx_tau . push_back( miniAODhelper.GetObjectJustBeforeDecay(  idx_sum[iParent]->daughter( iChild ) ) );
	    idx_tausParent.push_back( ThisParentIdx );
	  }else{
	    list.push_back( miniAODhelper.GetObjectJustBeforeDecay(  idx_sum[iParent]->daughter( iChild ) ) );
	    parentIdx.push_back( ThisParentIdx  ) ;
	  }
	}

      }

      // idx_rescured_zw ; 
      for( unsigned int iLep  = 0 ; iLep < idx_rescured_zw . size ()  ; iLep ++ ){

	if( fabs( idx_rescured_zw[iLep] -> pdgId() ) == 15 ){
	  idx_tau . push_back( miniAODhelper.GetObjectJustBeforeDecay(  idx_rescured_zw[iLep])  );
	  idx_tausParent . push_back(-1);
	}else{
	  list.push_back( idx_rescured_zw[iLep] );
	  parentIdx.push_back( -1 ) ; // -1 = no parent is assigned
	}

      }

      // fill tau;
      for( unsigned int iTau = 0 ; iTau  < idx_tau . size() ;  iTau ++){
	
	list.push_back( idx_tau[iTau] );
	parentIdx.push_back( idx_tausParent[iTau] );
	const int idx_ThisTau = parentIdx.size()-1;
	
	for( unsigned int iChild = 0 ; iChild < idx_tau[iTau] -> numberOfDaughters()  ; iChild ++ ){ 
	  list.push_back( miniAODhelper.GetObjectJustBeforeDecay(  idx_tau[iTau]->daughter( iChild ) ) );
	  parentIdx.push_back( idx_ThisTau  ) ;
	}

      }



      vfloat pt, eta, phi, m ;
      vint pdgid; 
      for( std::vector< const reco::Candidate *>::iterator p = list.begin();
	   p != list.end();
	   p ++ ){
	pt.push_back((*p) -> pt());
	eta.push_back((*p) -> eta() );
	phi.push_back( (*p) -> phi() );
	m.push_back( (*p) -> mass() );
	pdgid . push_back( (*p) -> pdgId() );
      }

      eve -> truth_pt_ = pt;
      eve -> truth_eta_ = eta;
      eve -> truth_phi_ = phi;
      eve -> truth_m_ = m;
      eve -> truth_pdgid_ = pdgid;
      eve -> truth_parentIdx_ = parentIdx ;

  } // end of "isMC".




  /////////
  ///
  /// Electrons
  ///
  ////////
  
  // miniAODhelper.SetElectronMVAinfo(h_conversioncollection, bsHandle);
  
  /////////
  ///
  /// Muons
  ///
  ////////
  //std::vector<pat::Muon> selectedMuons_tight = miniAODhelper.GetSelectedMuons( *muons, 25, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.1 );
  std::vector<pat::Muon> selectedMuons_tight = miniAODhelper.GetSelectedMuons( *muons, 15, muonID::muonTight, coneSize::R04, corrType::deltaBeta, 2.4);
  std::vector<pat::Muon> selectedMuons_loose = miniAODhelper.GetSelectedMuons( *muons, 15, muonID::muonTightDL, coneSize::R04, corrType::deltaBeta,2.4 );

  



  eve->passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1_ = ( passHLT_Ele27_eta2p1_WP85_Gsf_HT200_v1 ) ? 1 : 0;
  
  eve->passHLT_IsoMu24_v_ =  ( passHLT_IsoMu24_v) ? 1 : 0;
  eve->passHLT_IsoTkMu24_v_ =  ( passHLT_IsoTkMu24_v) ? 1 : 0;
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
  eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ = ( passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v ) ? 1 : 0;
  eve->passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = ( passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v ) ? 1 : 0;
  
  eve->passHLT_Ele25WP60_SC4_Mass55_v_ = ( passHLT_Ele25WP60_SC4_Mass55_v ) ? 1 : 0;

  // eve->passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_       = ( passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v       ) ? 1 : 0 ;
  // eve->passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_ = ( passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v ) ? 1 : 0 ;
  // eve->passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v_  = ( passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v  ) ? 1 : 0 ;
  eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_		    = ( passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v	      ) ? 1 : 0 ;
  eve->passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v_           = ( passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v           ) ? 1 : 0 ;

  eve->passHLT_Ele27_eta2p1_WPTight_Gsf_v_ =  ( passHLT_Ele27_eta2p1_WPTight_Gsf_v) ? 1 : 0;
  eve->passHLT_Ele27_WPTight_Gsf_v =  ( passHLT_Ele27_WPTight_Gsf_v) ? 1 : 0;
  eve->passHLT_IsoMu22_v_ =  ( passHLT_IsoMu22_v) ? 1 : 0;
  eve->passHLT_IsoTkMu22_v_ =  ( passHLT_IsoTkMu22_v) ? 1 : 0;
  eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_ =  ( passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v ) ? 1 : 0;


  eve->passHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_ = passHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v ? 1 : 0 ; 
  eve->passHLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v_ = passHLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v ? 1 : 0 ; 
  eve->passHLT_Ele32_WPTight_Gsf_L1DoubleEG_v_ = passHLT_Ele32_WPTight_Gsf_L1DoubleEG_v ? 1 : 0 ; 
  eve->passHLT_Ele32_WPTight_Gsf_v_ = passHLT_Ele32_WPTight_Gsf_v ? 1 : 0 ; 
  eve->passHLT_Ele35_WPTight_Gsf_v_ = passHLT_Ele35_WPTight_Gsf_v ? 1 : 0 ; 
  eve->passHLT_Ele38_WPTight_Gsf_v_ = passHLT_Ele38_WPTight_Gsf_v ? 1 : 0 ; 
  eve->passHLT_Ele40_WPTight_Gsf_v_ = passHLT_Ele40_WPTight_Gsf_v ? 1 : 0 ; 
  eve->passHLT_IsoMu24_2p1_v_ = passHLT_IsoMu24_2p1_v ? 1 : 0 ; 
  eve->passHLT_IsoMu27_v_ = passHLT_IsoMu27_v ? 1 : 0 ; 
  eve->passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v_ = passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v ? 1 : 0 ; 
  eve->passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ = passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v ? 1 : 0 ; 
  eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ = passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v ? 1 : 0 ; 
  eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_ = passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v ? 1 : 0 ; 
  eve->passHLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v_ = passHLT_PFHT380_SixJet32_DoubleBTagCSV_p075_v ? 1 : 0 ; 
  eve->passHLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v_ = passHLT_PFHT380_SixPFJet32_DoublePFBTagCSV_2p2_v ? 1 : 0 ; 
  eve->passHLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v_ = passHLT_PFHT380_SixPFJet32_DoublePFBTagDeepCSV_2p2_v ? 1 : 0 ; 
  eve->passHLT_PFHT430_SixJet40_BTagCSV_p080_v_ = passHLT_PFHT430_SixJet40_BTagCSV_p080_v ? 1 : 0 ; 
  eve->passHLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v_ = passHLT_PFHT430_SixPFJet40_PFBTagCSV_1p5_v ? 1 : 0 ; 

  // MET Filters, BBT 02-20-19
  eve->passMETFilter_Flag_goodVertices_v_ = passMETFilter_goodVertices_v; 
  eve->passMETFilter_Flag_globalTightHalo2016Filter_v_ = passMETFilter_globalTightHalo2016Filter_v; 
  eve->passMETFilter_Flag_globalSuperTightHalo2016Filter_v_ = passMETFilter_globalSuperTightHalo2016Filter_v; 
  eve->passMETFilter_Flag_HBHENoiseFilter_v_ = passMETFilter_HBHENoiseFilter_v;
  eve->passMETFilter_Flag_HBHENoiseIsoFilter_v_ = passMETFilter_HBHENoiseIsoFilter_v; 
  eve->passMETFilter_Flag_EcalDeadCellTriggerPrimitiveFilter_v_ = passMETFilter_EcalDeadCellTriggerPrimitiveFilter_v;
  eve->passMETFilter_Flag_BadPFMuonFilter_v_ = passMETFilter_BadPFMuonFilter_v; 
  eve->passMETFilter_Flag_BadChargedCandidateFilter_v_ = passMETFilter_BadChargedCandidateFilter_v;  
  eve->passMETFilter_Flag_ecalBadCalibFilter_v_ = passMETFilter_ecalBadCalibFilter_v;

  // 2017 MET triggers
  eve->passHLT_PFHT500_PFMET100_PFMHT100_IDTight_v_ = passHLT_PFHT500_PFMET100_PFMHT100_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFHT500_PFMET110_PFMHT110_IDTight_v_ = passHLT_PFHT500_PFMET110_PFMHT110_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFHT700_PFMET85_PFMHT85_IDTight_v_ = passHLT_PFHT700_PFMET85_PFMHT85_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFHT700_PFMET95_PFMHT95_IDTight_v_ = passHLT_PFHT700_PFMET95_PFMHT95_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFHT800_PFMET75_PFMHT75_IDTight_v_ = passHLT_PFHT800_PFMET75_PFMHT75_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFHT800_PFMET85_PFMHT85_IDTight_v_ = passHLT_PFHT800_PFMET85_PFMHT85_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFMET110_PFMHT110_IDTight_v_ = passHLT_PFMET110_PFMHT110_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFMET120_PFMHT120_IDTight_v_ = passHLT_PFMET120_PFMHT120_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFMET130_PFMHT130_IDTight_v_ = passHLT_PFMET130_PFMHT130_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFMET140_PFMHT140_IDTight_v_ = passHLT_PFMET140_PFMHT140_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFMETTypeOne110_PFMHT110_IDTight_v_ = passHLT_PFMETTypeOne110_PFMHT110_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFMETTypeOne120_PFMHT120_IDTight_v_ = passHLT_PFMETTypeOne120_PFMHT120_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFMETTypeOne130_PFMHT130_IDTight_v_ = passHLT_PFMETTypeOne130_PFMHT130_IDTight_v ? 1 : 0 ; 
  eve->passHLT_PFMETTypeOne140_PFMHT140_IDTight_v_ = passHLT_PFMETTypeOne140_PFMHT140_IDTight_v ? 1 : 0 ; 
  eve->passHLT_DiJet110_35_Mjj650_PFMET110_v_ = passHLT_DiJet110_35_Mjj650_PFMET110_v ? 1 : 0 ; 
  eve->passHLT_DiJet110_35_Mjj650_PFMET120_v_ = passHLT_DiJet110_35_Mjj650_PFMET120_v ? 1 : 0 ; 
  eve->passHLT_DiJet110_35_Mjj650_PFMET130_v_ = passHLT_DiJet110_35_Mjj650_PFMET130_v ? 1 : 0 ; 
  eve->passHLT_TripleJet110_35_35_Mjj650_PFMET110_v_ = passHLT_TripleJet110_35_35_Mjj650_PFMET110_v ? 1 : 0 ; 
  eve->passHLT_TripleJet110_35_35_Mjj650_PFMET120_v_ = passHLT_TripleJet110_35_35_Mjj650_PFMET120_v ? 1 : 0 ; 
  eve->passHLT_TripleJet110_35_35_Mjj650_PFMET130_v_ = passHLT_TripleJet110_35_35_Mjj650_PFMET130_v ? 1 : 0 ; 
  eve->passHLT_MET105_IsoTrk50_v_ = passHLT_MET105_IsoTrk50_v ? 1 : 0 ; 
  eve->passHLT_MET120_IsoTrk50_v_ = passHLT_MET120_IsoTrk50_v ? 1 : 0 ; 
  eve->passHLT_PFMET120_PFMHT120_IDTight_PFHT60_v_ = passHLT_PFMET120_PFMHT120_IDTight_PFHT60_v ? 1 : 0 ; 
  eve->passHLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v_ = passHLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v ? 1 : 0 ; 
  eve->passHLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v_ = passHLT_PFMET100_PFMHT100_IDTight_CaloBTagCSV_3p1_v ? 1 : 0 ; 
  eve->passHLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v_ = passHLT_PFMET110_PFMHT110_IDTight_CaloBTagCSV_3p1_v ? 1 : 0 ; 
  eve->passHLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v_ = passHLT_PFMET120_PFMHT120_IDTight_CaloBTagCSV_3p1_v ? 1 : 0 ; 
  eve->passHLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v_ = passHLT_PFMET130_PFMHT130_IDTight_CaloBTagCSV_3p1_v ? 1 : 0 ; 
  eve->passHLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v_ = passHLT_PFMET140_PFMHT140_IDTight_CaloBTagCSV_3p1_v ? 1 : 0 ; 
  eve->passHLT_CaloMET100_HBHECleaned_v_ = passHLT_CaloMET100_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_CaloMET250_HBHECleaned_v_ = passHLT_CaloMET250_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_CaloMET300_HBHECleaned_v_ = passHLT_CaloMET300_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_CaloMET350_HBHECleaned_v_ = passHLT_CaloMET350_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_CaloMET70_HBHECleaned_v_ = passHLT_CaloMET70_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_CaloMET80_HBHECleaned_v_ = passHLT_CaloMET80_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_CaloMET90_HBHECleaned_v_ = passHLT_CaloMET90_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_PFMET200_HBHE_BeamHaloCleaned_v_ = passHLT_PFMET200_HBHE_BeamHaloCleaned_v ? 1 : 0 ; 
  eve->passHLT_PFMET200_HBHECleaned_v_ = passHLT_PFMET200_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_PFMET250_HBHECleaned_v_ = passHLT_PFMET250_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_PFMET300_HBHECleaned_v_ = passHLT_PFMET300_HBHECleaned_v ? 1 : 0 ; 
  eve->passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v_ = passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v ? 1 : 0 ; 
  eve->passHLT_PFMET100_PFMHT100_IDTight_PFHT60_v_ = passHLT_PFMET100_PFMHT100_IDTight_PFHT60_v ? 1 : 0 ; 
  eve->passHLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v_ = passHLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v ? 1 : 0 ; 
  eve->passHLT_PFHT250_v_ = passHLT_PFHT250_v ? 1 : 0 ; 
  

  std::vector< std::pair<float , float> > SingleMuonTriggerDirection;
  std::vector< std::pair<float , float> > SingleElTriggerDirection;

  std::vector< std::pair<float , float> > DiMuonTriggerDirection;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);
  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
       
    unsigned long location1  = names.triggerName(i).find( "HLT_IsoMu27_v" ,0);
    unsigned long location3  = names.triggerName(i).find( "HLT_Ele35_WPTight_Gsf_v" ,0);
    
    if( ( location1 == 0 || location3 == 0 ) && triggerResults->accept(i) ){

      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
	obj.unpackPathNames(names);

	std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	std::vector<std::string> pathNamesLast = obj.pathNames(true);
	for (unsigned int iPathName = 0; iPathName < pathNamesAll.size(); iPathName ++ ) {

	  if( ( location1 == 0 && pathNamesAll[iPathName] . find( "HLT_IsoMu27_v" , 0 ) == 0 )
	      ||
	      ( location3 == 0 && pathNamesAll[iPathName] . find( "HLT_Ele35_WPTight_Gsf_v" , 0 ) == 0 )
	      ){

	    bool isBoth = obj.hasPathName( pathNamesAll[iPathName], true, true ); 
	    if (isBoth){

	      if( location3 == 0 ) {
		// electron trigger 
		SingleElTriggerDirection   . push_back( std::pair<float, float>(  obj.eta(), obj.phi() ) );
	      }else{
		// = muon trigger
		SingleMuonTriggerDirection . push_back( std::pair<float, float>(  obj.eta(), obj.phi() ) );
	      }


	    }// end if-"the trigger-path fired due to this object."


	  } // end if the trigger object matches the trigger.
	           
	} // end loop all trigger path
	       
	       
      } // end of trigger object loop
    } // end of if the trigger fired.
  } // end of trigger-bit loop (look into all HLT path)



  for (unsigned int i = 0; i < triggerResults->size(); ++i) {
       
    unsigned long loc1 = names.triggerName(i).find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v" ,0);
    unsigned long loc2 = names.triggerName(i).find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v" ,0);
    unsigned long loc3 = names.triggerName(i).find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v" ,0);


    if( ( loc1 == 0 || loc2 == 0 || loc3 == 0 ) && triggerResults->accept(i) ){

      for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
	obj.unpackPathNames(names);

	std::vector<std::string> pathNamesAll  = obj.pathNames(false);
	std::vector<std::string> pathNamesLast = obj.pathNames(true);
	for (unsigned int iPathName = 0; iPathName < pathNamesAll.size(); iPathName ++ ) {

	  if( ( loc1 == 0 && pathNamesAll[iPathName] . find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v" , 0 ) == 0 )
	      ||
	      ( loc2 == 0 && pathNamesAll[iPathName] . find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v" , 0 ) == 0 )
	      || 
	      ( loc3 == 0 && pathNamesAll[iPathName] . find( "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v" , 0 ) == 0 )
	      ){

	    bool isBoth = obj.hasPathName( pathNamesAll[iPathName], true, true ); 
	    if (isBoth){

	      DiMuonTriggerDirection . push_back( std::pair<float, float>(  obj.eta(), obj.phi() ) );

	    }// end if-"the trigger-path fired due to this object."


	  } // end if the trigger object matches the trigger.
	           
	} // end loop all trigger path
	       
	       
      } // end of trigger object loop
    } // end of if the trigger fired.
  } // end of trigger-bit loop (look into all HLT path)




  eve->run_ = run;
  eve->lumi_ = lumi;
  eve->evt_ = evt;

  eve->numTruePV_ = numTruePV;
  eve->numGenPV_ = numGenPV;

  eve->numPVs_ = numpv;
  eve->numSys_ = rNumSys;


  TLorentzVector leptonV;
  double eps = 0.0001;
  leptonV.SetPxPyPzE(eps, eps, eps, eps);

  vdouble vlepton;
      
  vvdouble vvleptons;

  
  vint lepton_genId, lepton_genParentId, lepton_genGrandParentId, lepton_trkCharge, lepton_charge, lepton_isMuon, lepton_isTight, lepton_isLoose, lepton_isLooseAlt, lepton_isTightV1;//, lepton_isTightV2;
  vdouble lepton_pt;
  vdouble lepton_eta;
  vdouble lepton_phi;
  vdouble lepton_e;
  vdouble lepton_seed;  // BBT 10-12-18
  vdouble lepton_IDSF;  // BBT 11-02-18
  vdouble lepton_recoIsoSF;  // BBT 11-02-18
  vdouble lepton_relIso;
  vdouble lepton_puppirelIso;
  vdouble lepton_dbiso_CH ;
  vdouble lepton_dbiso_NH ;
  vdouble lepton_dbiso_PH ;
  vdouble lepton_dbiso_PU ;
  vdouble lepton_puppiIsoWithLep_CH   ;
  vdouble lepton_puppiIsoWithLep_NH   ;
  vdouble lepton_puppiIsoWithLep_PH   ;
  vdouble lepton_puppiIsoWithoutLep_CH;
  vdouble lepton_puppiIsoWithoutLep_NH;
  vdouble lepton_puppiIsoWithoutLep_PH;
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
  vdouble lepton_dRSingleLepTrig ; 
  vdouble lepton_dRDiLepTrig ; 
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
    int charge = iMu->charge();
   
    int isPOGTight = miniAODhelper.passesMuonPOGIdTight(*iMu) ? 1 : 0;
    int isPOGLoose = 1 ; //this is needed for the consistency of variables with electron.


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

    // our pre-selections 
    if( iMu->pt() < 15 ){ continue;}
    if( fabs( iMu->eta() ) > 2.4 ){ continue;}

    //


    lepton_trkCharge.push_back(trkCharge);
    lepton_charge.push_back(charge);
    lepton_isMuon.push_back(1);
    lepton_isTight.push_back(isPOGTight);
    lepton_isTightV1.push_back(isPOGTight);
    //lepton_isTightV2.push_back(isPOGTight);
    lepton_isLoose.push_back(isPOGLoose);
    lepton_isLooseAlt.push_back(isPOGLoose);
    lepton_genId.push_back(genId);
    lepton_genParentId.push_back(genParentId);
    lepton_genGrandParentId.push_back(genGrandParentId);


    {
      double _mu_pt = iMu->pt() ;
      double sf = 1.0 ; 
      if( isMC ){
//	TRandom3 rnd ;
//        rnd.SetSeed((uint32_t)( iMu -> userInt("deterministicSeed")));
//	if( (iMu->genLepton()) ){// todo
//	  sf = muon_roc->kScaleFromGenMC ( trkCharge , _mu_pt, iMu->eta(), iMu->phi(), trackerLayersWithMeasurement, iMu->genLepton()->pt() , rnd.Rndm() ) ;
//	}else{
//	  sf = muon_roc->kScaleAndSmearMC( trkCharge , _mu_pt, iMu->eta(), iMu->phi() , trackerLayersWithMeasurement, rnd.Rndm(), rnd.Rndm() ) ;
//	}
	
      }else{
	// = data.
	// sf = muon_roc->kScaleDT( trkCharge , _mu_pt, iMu->eta(), iMu->phi());
      }
      
      lepton_pt.push_back(iMu->pt() * sf  );
    }


    // muon SF stuff, BBT 11-02-18
    double lepIDSF  = ( ! isMC ? 1 : scalefactors.getTightMuon_IDSF_single( iMu->pt(), iMu->eta() ) );
    double lepIsoSF = ( ! isMC ? 1 : scalefactors.getTightMuon_IsoSF_single( iMu->pt(), iMu->eta() ) );
    //std::cout << "muon IDSF: " << lepIDSF <<" , muon IsoSF: "<< lepIsoSF << std::endl;

    lepton_eta.push_back(iMu->eta());
    lepton_phi.push_back(iMu->phi());
    lepton_e.push_back(iMu->energy());
    //lepton_seed.push_back( iMu -> userInt("deterministicSeed") );  // BBT 10-12-18
    lepton_relIso.push_back( miniAODhelper.GetMuonRelIso(*iMu, coneSize::R04, corrType::deltaBeta) ) ;
    lepton_IDSF.push_back(lepIDSF);
    lepton_recoIsoSF.push_back(lepIsoSF);

    //lepton_puppirelIso.push_back( iMu -> userFloat( "reliso_puppicombined" ) ) ;
    lepton_dbiso_CH . push_back( iMu -> pfIsolationR04().sumChargedHadronPt );
    lepton_dbiso_NH . push_back( iMu -> pfIsolationR04().sumNeutralHadronEt );
    lepton_dbiso_PH . push_back( iMu -> pfIsolationR04().sumPhotonEt );
    lepton_dbiso_PU . push_back( iMu -> pfIsolationR04().sumPUPt );
    //lepton_puppiIsoWithLep_CH    . push_back( iMu -> userFloat( "reliso_puppiwithlepton_ch" )   );
    //lepton_puppiIsoWithLep_NH    . push_back( iMu -> userFloat( "reliso_puppiwithlepton_nh" )   );
    //lepton_puppiIsoWithLep_PH    . push_back( iMu -> userFloat( "reliso_puppiwithlepton_ph" )   );
    //lepton_puppiIsoWithoutLep_CH . push_back( iMu -> userFloat( "reliso_puppinolepton_ch" )   );
    //lepton_puppiIsoWithoutLep_NH . push_back( iMu -> userFloat( "reliso_puppinolepton_nh" )   );
    //lepton_puppiIsoWithoutLep_PH . push_back( iMu -> userFloat( "reliso_puppinolepton_ph" )   );
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

    float lep_trig_dr = 10 ; 
    for( std::vector< std::pair<float,float>>::iterator trig = SingleMuonTriggerDirection.begin ();
	 trig != SingleMuonTriggerDirection.end() ; 
	 trig ++ ){
      float d_eta = iMu->eta() - trig->first;
      float d_phi = iMu->phi() - trig->second;
      d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ; 
      double dr =  sqrt( d_eta*d_eta + d_phi*d_phi );
      lep_trig_dr = ( lep_trig_dr < dr ) ? lep_trig_dr : dr ;
    }
    lepton_dRSingleLepTrig.push_back( lep_trig_dr );

    {
      float lep_trig_dr = 10 ; 
      for( std::vector< std::pair<float,float>>::iterator trig = DiMuonTriggerDirection.begin ();
	   trig != DiMuonTriggerDirection.end() ; 
	   trig ++ ){
	float d_eta = iMu->eta() - trig->first;
	float d_phi = iMu->phi() - trig->second;
	d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ; 
	double dr =  sqrt( d_eta*d_eta + d_phi*d_phi );
	lep_trig_dr = ( lep_trig_dr < dr ) ? lep_trig_dr : dr ;
      }
      lepton_dRDiLepTrig.push_back( lep_trig_dr );
    }
  }

  // Loop over electrons
  // for( std::vector<pat::Electron>::const_iterator iEle = electrons.begin(); iEle != electrons.end(); iEle++ ){ 
  for( std::vector<pat::Electron>::const_iterator iEle = electrons->begin(); iEle != electrons->end(); iEle++ ){ 

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
    int charge = iEle->charge();

    // BBT, 10-10-18
    //auto corrP4  = pat::Electron::p4() * pat::Electron::userFloat("ecalTrkEnergyPostCorr") / pat::Electron::energy();
    auto corrP4  = iEle->p4() * iEle->userFloat("ecalTrkEnergyPostCorr") / iEle->energy();
    //std::cout << "[Electron loop] corrected pt: " << corrP4.pt() << ", nominal pt: " << iEle->pt() << std::endl;

    // revert to usual vector if running on data, comment out for test 11-05-18
    //if (!isMC)
    //  corrP4 = iEle->p4();

    bool inCrack = false;
    double scEta = -99;
    if( iEle->superCluster().isAvailable() ){
      inCrack = ( fabs(iEle->superCluster()->position().eta())>1.4442 
		  && 
		  fabs(iEle->superCluster()->position().eta())<1.5660 );
      scEta = iEle->superCluster()->position().eta();
    }

    //int isPOGTight = iEle->electronID("cutBasedElectronID-Fall17-94X-V1-tight") && ? 1 : 0  ; // BBT 10-17-18
    int isPOGTight = ! inCrack && miniAODhelper.PassElectron94XId(*iEle ,electronID::electron94XCutBasedT ) ? 1 : 0  ; // ygg core
    int isPOGLoose = ! inCrack && miniAODhelper.PassElectron94XId(*iEle ,electronID::electron94XCutBasedV ) ? 1 : 0  ;
    int isPOGLooseAlt = ! inCrack && miniAODhelper.PassElectron94XId(*iEle ,electronID::electron94XCutBasedL ) ? 1 : 0  ;

    //int isPOGTightV2 = !inCrack && miniAODhelper.PassElectron94XId(*iEle ,electronID::electron94XCutBasedT_V2 ) ? 1 : 0; // new V2 definition
    //int isPOGTightV2 = isPOGTight;

    //bool isPassEleId = (*ele_id_decisions)[*iEle];
    //int isVIDTight = ! inCrack && isPassEleId ? 1 : 0  ;

    // our pre-selections 
    if( iEle->pt() < 15 ){ continue;}
    if( fabs( iEle->eta() ) > 2.4 ){ continue;}
    //


    double mvaTrigValue = -99;//myMVATrig->mvaValue(*iEle,false);

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

    // electron SF stuff, BBT 11-02-18
    double lepIDSF =  ( ! isMC ? 1 : scalefactors.getTightElectron_IDSF_single( corrP4.pt(), scEta ) );
    double lepRecoSF =  ( ! isMC ? 1 : scalefactors.getTightElectron_RecoSF_single( corrP4.pt(), scEta ) ) ; 
    //std::cout << "electron IDSF: " << lepIDSF <<" , electron RecoSF: "<< lepRecoSF << std::endl;


    lepton_trkCharge.push_back(trkCharge);
    lepton_charge.push_back(charge);
    lepton_isMuon.push_back(0);
    //lepton_isTight.push_back( iEle->electronID("cutBasedElectronID-Fall17-94X-V1-tight") ); // BBT, 10-10-18
    lepton_isTight.push_back(isPOGTight);
    lepton_isTightV1.push_back(isPOGTight);
    //lepton_isTightV2.push_back(isPOGTightV2);
    lepton_isLoose.push_back(isPOGLoose);
    lepton_isLooseAlt.push_back(isPOGLooseAlt);
    lepton_genId.push_back(genId);
    lepton_genParentId.push_back(genParentId);
    lepton_genGrandParentId.push_back(genGrandParentId);
    //lepton_pt.push_back(iEle->pt()); // ygg core
    lepton_pt.push_back( corrP4.pt() ); // BBT, 10-10-18
    lepton_eta.push_back(iEle->eta());
    lepton_phi.push_back(iEle->phi());
    //lepton_e.push_back(iEle->energy()); // ygg core
    lepton_e.push_back( corrP4.energy() ); // BBT, 10-10-18
    //lepton_seed.push_back( iEle -> userInt("deterministicSeed") ); // BBT, 10-12-18
    lepton_IDSF.push_back(lepIDSF);
    lepton_recoIsoSF.push_back(lepRecoSF);
    lepton_relIso.push_back(     miniAODhelper.GetElectronRelIso(*iEle, coneSize::R03, corrType::rhoEA,effAreaType::fall17) );
    lepton_puppirelIso.push_back(miniAODhelper.GetElectronRelIso(*iEle, coneSize::R03, corrType::rhoEA,effAreaType::fall17) );
    lepton_dbiso_CH . push_back(0);
    lepton_dbiso_NH . push_back(0);
    lepton_dbiso_PH . push_back(0);
    lepton_dbiso_PU . push_back(0);
    lepton_puppiIsoWithLep_CH   .push_back(0);
    lepton_puppiIsoWithLep_NH   .push_back(0);
    lepton_puppiIsoWithLep_PH   .push_back(0);
    lepton_puppiIsoWithoutLep_CH.push_back(0);
    lepton_puppiIsoWithoutLep_NH.push_back(0);
    lepton_puppiIsoWithoutLep_PH.push_back(0);
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


    float lep_trig_dr = 10 ; 
    for( std::vector< std::pair<float,float>>::iterator trig = SingleElTriggerDirection.begin ();
	 trig != SingleElTriggerDirection.end() ; 
	 trig ++ ){
      float d_eta = scEta       - trig->first;
      float d_phi = iEle->phi() - trig->second;
      d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ; 
      double dr =  sqrt( d_eta*d_eta + d_phi*d_phi );
      lep_trig_dr = ( lep_trig_dr < dr ) ? lep_trig_dr : dr ;
    }
    lepton_dRSingleLepTrig.push_back( lep_trig_dr );

    {
      lepton_dRDiLepTrig.push_back( 10 ); // <- not prepared for di-electron yet.
    }
  }

  eve->wgt_lumi_  = intLumi_;
  eve->wgt_xs_    = mySample_xSec_;//mySample.xSec;
  eve->wgt_nGen_  = mySample_nGen_;//mySample.nGen;

  eve->wgt_generator_ = GenEventInfoWeight;

  h_event -> Fill( 0.0 , GenEventInfoWeight > 0 ? +1.0 : -1.0 ) ; // with weight 
  h_event -> Fill( 1.0 ) ; // no weight

  if( evt % 2 == 0 ){
    h_event -> Fill( 2.0 , GenEventInfoWeight > 0 ? +1.0 : -1.0 ) ; // with weight 
  }else{
    h_event -> Fill( 3.0 , GenEventInfoWeight > 0 ? +1.0 : -1.0 ) ; // with weight 
  }
  
  if( ! EVENTSYNCMODE && lepton_pt . size() == 0 ){
    return ;  // No data recording 
  }
  

  eve ->  genjet_pt_  = genjet_pt ;
  eve ->  genjet_eta_ = genjet_eta ; 
  eve ->  genjet_phi_ = genjet_phi ;  
  eve ->  genjet_m_   = genjet_m ; 
  eve ->  genjet_BhadronMatch_ = genjet_BhadronMatch ; 
  
  eve ->  fatgenjet_pt_  = fatgenjet_pt ;
  eve ->  fatgenjet_eta_ = fatgenjet_eta ; 
  eve ->  fatgenjet_phi_ = fatgenjet_phi ;  
  eve ->  fatgenjet_m_   = fatgenjet_m ; 



  eve->lepton_trkCharge_        = lepton_trkCharge;
  eve->lepton_charge_           = lepton_charge;
  eve->lepton_isMuon_           = lepton_isMuon;
  eve->lepton_isTight_          = lepton_isTight;
  eve->lepton_isTightV1_        = lepton_isTightV1;
  //eve->lepton_isTightV2_        = lepton_isTightV2;
  eve->lepton_isLoose_          = lepton_isLoose;
  eve->lepton_isLooseAlt_       = lepton_isLooseAlt ;
  eve->lepton_pt_               = lepton_pt;
  eve->lepton_eta_              = lepton_eta;
  eve->lepton_phi_              = lepton_phi;
  eve->lepton_e_                = lepton_e;
  eve->lepton_seed_             = lepton_seed;  // BBT 10-12-18
  eve->lepton_IDSF_             = lepton_IDSF;  // BBT 11-02-18
  eve->lepton_recoIsoSF_        = lepton_recoIsoSF;  // BBT 11-02-18
  eve->lepton_relIso_           = lepton_relIso;
  eve->lepton_puppirelIso_      = lepton_puppirelIso;

  eve -> lepton_dbiso_CH_    = lepton_dbiso_CH    ;
  eve -> lepton_dbiso_NH_    = lepton_dbiso_NH    ;
  eve -> lepton_dbiso_PH_    = lepton_dbiso_PH    ;
  eve -> lepton_dbiso_PU_    = lepton_dbiso_PU    ;

  eve -> lepton_puppiIsoWithLep_CH_    = lepton_puppiIsoWithLep_CH    ;
  eve -> lepton_puppiIsoWithLep_NH_    = lepton_puppiIsoWithLep_NH    ;
  eve -> lepton_puppiIsoWithLep_PH_    = lepton_puppiIsoWithLep_PH    ;
  eve -> lepton_puppiIsoWithoutLep_CH_ = lepton_puppiIsoWithoutLep_CH ;
  eve -> lepton_puppiIsoWithoutLep_NH_ = lepton_puppiIsoWithoutLep_NH ;
  eve -> lepton_puppiIsoWithoutLep_PH_ = lepton_puppiIsoWithoutLep_PH ;

  eve->lepton_scEta_           = lepton_scEta;

  eve->lepton_dRSingleLepTrig_        = lepton_dRSingleLepTrig ;
  eve->lepton_dRDiLepTrig_        = lepton_dRDiLepTrig ;



  bool b_TheEventHasFourJets_ForAtLeastOneSystematicVariation = false ; 

  // Loop over systematics
  for( int iSys=0; iSys<rNumSys; iSys++ ){

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
      
      // JES up 
    case 25 : iSysType =  sysType::JESPileUpDataMCup; break ;
    case 26 : iSysType =  sysType::JESPileUpPtRefup; break ;
    case 27 : iSysType =  sysType::JESPileUpPtBBup; break ;
    case 28 : iSysType =  sysType::JESPileUpPtEC1up; break ;
    case 29 : iSysType =  sysType::JESPileUpPtEC2up; break ;
    case 30 : iSysType =  sysType::JESPileUpPtHFup; break ;
    case 31 : iSysType =  sysType::JESRelativeJEREC1up; break ;
    case 32 : iSysType =  sysType::JESRelativeJEREC2up; break ;
    case 33 : iSysType =  sysType::JESRelativeJERHFup; break ;
    case 34 : iSysType =  sysType::JESRelativeFSRup; break ;
    case 35 : iSysType =  sysType::JESRelativeStatFSRup; break ;
      //    case 36 : iSysType =  sysType::JESRelativeStatEC2up; break ; ) RelativeStatEC2 is not supported._
    case 37 : iSysType =  sysType::JESRelativeStatECup; break ;
    case 38 : iSysType =  sysType::JESRelativeStatHFup; break ;
    case 39 : iSysType =  sysType::JESRelativePtBBup; break ;
    case 40 : iSysType =  sysType::JESRelativePtEC1up; break ;
    case 41 : iSysType =  sysType::JESRelativePtEC2up; break ;
    case 42 : iSysType =  sysType::JESRelativePtHFup; break ;
    case 43 : iSysType =  sysType::JESTimeEtaup; break ;
    case 44 : iSysType =  sysType::JESAbsoluteScaleup; break ;
    case 45 : iSysType =  sysType::JESAbsoluteMPFBiasup; break ;
    case 46 : iSysType =  sysType::JESAbsoluteStatup; break ;
    case 47 : iSysType =  sysType::JESSinglePionECALup; break ;
    case 48 : iSysType =  sysType::JESSinglePionHCALup; break ;
    case 49 : iSysType =  sysType::JESFragmentationup; break ;
    case 50 : iSysType =  sysType::JESTimePtup; break ;
    case 51 : iSysType =  sysType::JESFlavorQCDup; break ;

    case 52 : iSysType =  sysType::JESPileUpDataMCdown; break ;
    case 53 : iSysType =  sysType::JESPileUpPtRefdown; break ;
    case 54 : iSysType =  sysType::JESPileUpPtBBdown; break ;
    case 55 : iSysType =  sysType::JESPileUpPtEC1down; break ;
    case 56 : iSysType =  sysType::JESPileUpPtEC2down; break ;
    case 57 : iSysType =  sysType::JESPileUpPtHFdown; break ;
    case 58 : iSysType =  sysType::JESRelativeJEREC1down; break ;
    case 59 : iSysType =  sysType::JESRelativeJEREC2down; break ;
    case 60 : iSysType =  sysType::JESRelativeJERHFdown; break ;
    case 61 : iSysType =  sysType::JESRelativeFSRdown; break ;
    case 62 : iSysType =  sysType::JESRelativeStatFSRdown; break ;
      //    case 63 : iSysType =  sysType::JESRelativeStatEC2down; break ; (  RelativeStatEC2 is not supported. )
    case 64 : iSysType =  sysType::JESRelativeStatECdown; break ;
    case 65 : iSysType =  sysType::JESRelativeStatHFdown; break ;
    case 66 : iSysType =  sysType::JESRelativePtBBdown; break ;
    case 67 : iSysType =  sysType::JESRelativePtEC1down; break ;
    case 68 : iSysType =  sysType::JESRelativePtEC2down; break ;
    case 69 : iSysType =  sysType::JESRelativePtHFdown; break ;
    case 70 : iSysType =  sysType::JESTimeEtadown; break ;
    case 71 : iSysType =  sysType::JESAbsoluteScaledown; break ;
    case 72 : iSysType =  sysType::JESAbsoluteMPFBiasdown; break ;
    case 73 : iSysType =  sysType::JESAbsoluteStatdown; break ;
    case 74 : iSysType =  sysType::JESSinglePionECALdown; break ;
    case 75 : iSysType =  sysType::JESSinglePionHCALdown; break ;
    case 76 : iSysType =  sysType::JESFragmentationdown; break ;
    case 77 : iSysType =  sysType::JESTimePtdown; break ;
    case 78 : iSysType =  sysType::JESFlavorQCDdown; break ;

    default: iSysType = sysType::NA;       break;
    }

    /////////
    ///
    /// Pfjets
    ///
    ////////
    //std::cout << "START AK4 jets" << std::endl;

    std::vector<pat::Jet> rawJets = miniAODhelper.GetUncorrectedJets( *pfjets );
    //cout << "!!! before ak4pfchs correction" << endl;
    //std::vector<pat::Jet> correctedJets =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection , iSysType ); // ygg core
    //std::vector<pat::Jet> correctedJets =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection , iSysType, true, true, true ); // BBT, 10-29-18
    //                                                                                                                         doJES, doJER, isAK4
    std::vector<pat::Jet> correctedJets =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection , iSysType, true,  true, true ); // BBT, 10-29-18
    //cout << "!!! after ak4pfchs correction" << endl;
    std::vector<pat::Jet> selectedJets_unsorted =  miniAODhelper.GetSelectedJets(correctedJets, 20., 5.0 ,
										 ( jetID::jetTight ) // <- For 2017, no LooseID, only tight.
										 , '-' );
    //std::cout << "END AK4 jets" << std::endl;
    
    double JecUpdatePropagationToMET_x = 0 ;
    double JecUpdatePropagationToMET_y = 0 ;
    for( unsigned int i = 0 ; i < pfjets -> size() ; i++ ){
      /* // BBT 10-05-18
      if(iSys == 0 ){
	std::cout <<" JEC diff " << i << " " << pfjets->at(i).px() - correctedJets.at(i).px() << " " << pfjets->at(i).py() - correctedJets.at(i).py()  << std::endl ; 
	std::cout <<" JEC  orig " << i <<" "<< pfjets->at(i).px() << " " <<  pfjets->at(i).py() << std::endl ; 
      }
      */
      JecUpdatePropagationToMET_x +=  pfjets->at(i).px() - correctedJets.at(i).px() ; 
      JecUpdatePropagationToMET_y +=  pfjets->at(i).py() - correctedJets.at(i).py() ; 
    }

    ///// HEP top tagged jet
    int numTopTags = 0;
    int  n_fatjets=0;
    ///// Higgs tagged jet
    int numHiggsTags = 0;
    
    boosted::BoostedJetCollection selectedBoostedJets;
    if(h_boostedjet.isValid()){
      boosted::BoostedJetCollection const &boostedjets_unsorted = *h_boostedjet;
      
      //    boosted::BoostedJetCollection boostedjets = BoostedUtils::GetSortedByPt(boostedjets_unsorted);
      boosted::BoostedJetCollection boostedjets = boostedjets_unsorted;
      
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
      //if( !(topJet->topjet.pt() > 200. && abs(topJet->topjet.eta()) < 2) ) continue;

      // pt and eta requirements on subjets
      if( !( (topJet->nonW.pt()>20 && abs(topJet->nonW.eta())<2.4 ) &&
	     (topJet->W1.pt()>20 && abs(topJet->W1.eta())<2.4 ) &&
	     (topJet->W2.pt()>20 && abs(topJet->W2.eta())<2.4 ) ) ) continue;


       if(toptagger.GetTopTaggerOutput(*topJet)>-1){
       	numTopTags++;
       	syncTopJets.push_back(*topJet);
       }

    }

    
    for( boosted::BoostedJetCollection::iterator higgsJet = boostedjets.begin() ; higgsJet != boostedjets.end(); ++higgsJet ){
      // pt and eta requirements on top jet
      if( !(higgsJet->fatjet.pt() > 200. && abs(higgsJet->fatjet.eta()) < 2) ) continue;
      numHiggsTags++;

      //remove overlap with topjets
       bool overlapping=false;
       for(auto tj=syncTopJets.begin(); tj!=syncTopJets.end(); tj++){

//       	if(BoostedUtils::DeltaR(tj->fatjet,higgsJet->fatjet)<1.5){
//       	  overlapping=true;
//       	  break;
//       	}

       }
       if(overlapping) continue;
       if(overlapping) continue;
    std::vector<pat::Jet> filterjets = higgsJet->filterjets;
       
    }

    }
         
    pat::MET correctedMET = pfmet->at(0); 

    std::vector<double> csvV;
    std::vector<double> jet_combinedMVABJetTags;
    std::vector<double> jet_combinedInclusiveSecondaryVertexV2BJetTags;
    std::vector<double> jet_combinedMVABJetTags_HIP;
    std::vector<double> jet_combinedInclusiveSecondaryVertexV2BJetTags_HIP;

    std::vector<double> jet_DeepCSV_b;
    std::vector<double> jet_DeepCSV_bb;

    std::vector<double> jet_vtxMass;
    std::vector<double> jet_vtxNtracks;
    std::vector<double> jet_vtx3DVal;
    std::vector<double> jet_vtx3DSig;

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

    vdouble jet_pt;
    vdouble jet_eta;
    vdouble jet_phi;
    vdouble jet_m;

    vint    jet_puid;
    vint    jet_seed;  // BBT 10-12-18
    vdouble jet_DeepCSV_SF;  // BBT 10-12-18

    vdouble jet_AssociatedGenJet_pt;
    vdouble jet_AssociatedGenJet_eta;
    vdouble jet_AssociatedGenJet_phi;
    vdouble jet_AssociatedGenJet_m;

    vdouble jet_precore_phi ; 
    vdouble jet_precore_pt ; 

    vint jet_genId_vect;
    vint jet_partonflavour_vect;
    vint jet_flavour_vect;
    vint jet_genParentId_vect;
    vint jet_genGrandParentId_vect;


    // Loop over selected jets

    int ijet = 0 ;
    for( std::vector<pat::Jet>::const_iterator iJet = selectedJets_unsorted.begin(); iJet != selectedJets_unsorted.end(); iJet++ ){ 
      jet_pt  .push_back( iJet -> pt()  );
      jet_phi .push_back( iJet -> phi() );
      jet_eta .push_back( iJet -> eta() );
      jet_m   .push_back( iJet -> mass()   );

      jet_puid . push_back( iJet -> userInt("pileupJetId:fullId") ) ;
      //std::cout << " %%% jet puid: " << iJet -> userInt("pileupJetId:fullId") << std::endl;     
      if (isMC)
	jet_seed . push_back( iJet -> userInt("deterministicSeed") ) ; // BBT 10-12-18
      //std::cout << " %%% jet pt: " << iJet->pt() << " , jet eta: " << iJet->eta() << " , jet seed: " << iJet -> userInt("deterministicSeed") << std::endl;

      // DeepCSV SFs, BBT 11-02-18
      double deepCSV_SF  = ( ! isMC ? 1 : scalefactors.get_csv_wgt_single( iJet->hadronFlavour(), iJet -> eta(), iJet -> pt(), iJet->bDiscriminator("pfDeepCSVJetTags:probb")+iJet->bDiscriminator("pfDeepCSVJetTags:probbb")) );
      jet_DeepCSV_SF.push_back( deepCSV_SF );

      jet_precore_pt . push_back( iJet->userFloat( "OrigPt"  ) );
      jet_precore_phi. push_back( iJet->userFloat( "OrigPhi" ) );
   
      const reco::GenJet* ref = iJet -> genJet();
      if (ref) {
	jet_AssociatedGenJet_pt  .push_back( ref -> pt() );
	jet_AssociatedGenJet_eta .push_back( ref -> eta() );
	jet_AssociatedGenJet_phi .push_back( ref -> phi() );
	jet_AssociatedGenJet_m   .push_back( ref -> mass() );
      } else {
	jet_AssociatedGenJet_pt  .push_back( -999 ) ;
	jet_AssociatedGenJet_eta .push_back( -999 ) ;
	jet_AssociatedGenJet_phi .push_back( -999 ) ;
	jet_AssociatedGenJet_m   .push_back( -999 );
      }

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
      
	  
      // make vvdouble version of vecTLorentzVector
      vdouble vjets;
      vjets.push_back(iJet->px());
      vjets.push_back(iJet->py());
      vjets.push_back(iJet->pz());
      vjets.push_back(iJet->energy());
      vvjets.push_back(vjets);

      // Get CSV discriminant, check if passes Med WP 
      double myCSV = iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
      csvV.push_back(myCSV);

      jet_combinedMVABJetTags.push_back( iJet->bDiscriminator("pfCombinedMVAV2BJetTags") );
      jet_combinedInclusiveSecondaryVertexV2BJetTags.push_back( iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );

      jet_DeepCSV_b .push_back( iJet->bDiscriminator("pfDeepCSVJetTags:probb") );
      jet_DeepCSV_bb.push_back( iJet->bDiscriminator("pfDeepCSVJetTags:probbb") );

      jet_combinedMVABJetTags_HIP.push_back( iJet->bDiscriminator("newpfCombinedMVAV2BJetTags") );
      jet_combinedInclusiveSecondaryVertexV2BJetTags_HIP.push_back( iJet->bDiscriminator("newpfCombinedInclusiveSecondaryVertexV2BJetTags") );


    }// end loop over iJet


    std::vector<pat::Jet> fatrawJets = miniAODhelper_fatjet.GetUncorrectedJets( *fatjets );

    /* 02-05-19
    std::vector<double>  fatjet_pt            ;
    std::vector<double>  fatjet_eta	      ;
    std::vector<double>  fatjet_phi	      ;
    std::vector<double>  fatjet_m  	      ;
    std::vector<double>  fatjet_doublebtagging  	      ;
    std::vector<int>     fatjet_nSubjet 	      ;
    std::vector<double>  fatjet_sdmass_miniaod ;
    std::vector<double>  fatjet_sdmass_uncorr  ;
    std::vector<double>  fatjet_tau1	      ;
    std::vector<double>  fatjet_tau2	      ;
    std::vector<double>  fatjet_tau3	      ;
    std::vector<double>  fatjet_tau4	      ;
    std::vector<double>  fatjet_chstau1	      ;
    std::vector<double>  fatjet_chstau2	      ;
    std::vector<double>  fatjet_chstau3	      ;
    std::vector<double>  fatjet_nb1N2 	      ;
    std::vector<double>  fatjet_nb1N3 	      ;
    std::vector<double>  fatjet_nb2N2 	      ;
    std::vector<double>  fatjet_nb2N3 	      ;
    std::vector<double>  fatjet_chsprunedmass  ;

    std::vector<std::vector<double> >fatjet_subjet_pt ;
    std::vector<std::vector<double> >fatjet_subjet_eta;
    std::vector<std::vector<double> >fatjet_subjet_phi;
    std::vector<std::vector<double> >fatjet_subjet_m  ;
    std::vector<std::vector<double> >fatjet_subjet_beepcsv;
    std::vector<std::vector<double> >fatjet_subjet_csvv2  ;

    std::vector<double>  re_fatjet_pt            ;
    std::vector<double>  re_fatjet_eta	      ;
    std::vector<double>  re_fatjet_phi	      ;
    std::vector<double>  re_fatjet_tau21      ;
    std::vector<double>  re_fatjet_tau32      ;
    std::vector<double>  re_fatjet_sdmass_miniaod ;
    std::vector<double>  re_fatjet_sdmass_uncorr  ;

    for( unsigned int i = 0 ; i < fatjets -> size() ; i++  ){
      pat::Jet correctedJet =  miniAODhelper_fatjet.GetCorrectedAK8Jet( fatrawJets.at(i) , iEvent, iSetup, fatgenjetCollection , iSysType );
      pat::Jet originalJet  = fatjets->at(i); 


    std::vector<double> subjet_pt ;    
    std::vector<double> subjet_eta;    
    std::vector<double> subjet_phi;    
    std::vector<double> subjet_m  ;    
    std::vector<double> subjet_beepcsv;
    std::vector<double> subjet_csvv2  ;


      TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
      auto const & sdSubjetsPuppi = originalJet.subjets("SoftDropPuppi");
      long n = 0 ;
      for ( auto const & it : sdSubjetsPuppi ) {
	n ++ ; 
	puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
	puppi_softdrop+=puppi_softdrop_subjet;

	subjet_pt      . push_back( it -> pt()  );    
	subjet_eta     . push_back( it -> eta() );    
	subjet_phi     . push_back( it -> phi() );    
	subjet_m       . push_back( it -> mass());    
	subjet_beepcsv . push_back( it -> bDiscriminator("pfDeepCSVJetTags:probb")  +  it -> bDiscriminator("pfDeepCSVJetTags:probbb")  );
	subjet_csvv2   . push_back( it -> bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")  );
      }

      fatjet_subjet_pt       .push_back( subjet_pt      );
      fatjet_subjet_eta      .push_back( subjet_eta     );
      fatjet_subjet_phi      .push_back( subjet_phi     );
      fatjet_subjet_m        .push_back( subjet_m       );
      fatjet_subjet_beepcsv  .push_back( subjet_beepcsv );
      fatjet_subjet_csvv2    .push_back( subjet_csvv2   );

      fatjet_pt . push_back( correctedJet .pt() ) ;  
      fatjet_eta. push_back( correctedJet .eta() ) ;  
      fatjet_phi. push_back( correctedJet .phi() ) ;  
      fatjet_m  . push_back( correctedJet .mass() ) ;  
      fatjet_doublebtagging . push_back( correctedJet . bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags") ) ;  

      fatjet_nSubjet  . push_back( n ) ;  

      fatjet_sdmass_miniaod . push_back(  originalJet . userFloat("ak8PFJetsPuppiSoftDropMass") ) ;
      fatjet_sdmass_uncorr  . push_back( puppi_softdrop . M()  ) ; 

      fatjet_tau1 . push_back( originalJet . userFloat("NjettinessAK8Puppi:tau1") ) ; 
      fatjet_tau2 . push_back( originalJet . userFloat("NjettinessAK8Puppi:tau2") ) ; 
      fatjet_tau3 . push_back( originalJet . userFloat("NjettinessAK8Puppi:tau3") ) ; 
      fatjet_tau4 . push_back( originalJet . userFloat("NjettinessAK8Puppi:tau4") ) ; 

      fatjet_chstau1 . push_back( originalJet . userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1") );
      fatjet_chstau2 . push_back( originalJet . userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2") );
      fatjet_chstau3 . push_back( originalJet . userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3") );

      fatjet_nb1N2 . push_back( originalJet . userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2") );
      fatjet_nb1N3 . push_back( originalJet . userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3") );
      fatjet_nb2N2 . push_back( originalJet . userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN2") );
      fatjet_nb2N3 . push_back( originalJet . userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN3") );

      fatjet_chsprunedmass. push_back( originalJet . userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass") );

    }
    */

    /*
    { /// --------- Add information of reclusted (by myself) jets.
      std::vector<pat::Jet> rerun_fatrawJets = miniAODhelper_fatjet.GetUncorrectedJets( *rerun_fatjets );

      for( unsigned int i = 0 ; i < rerun_fatrawJets . size() ; i++  ){
	pat::Jet correctedJet =  miniAODhelper_fatjet.GetCorrectedAK8Jet( rerun_fatrawJets.at(i) , iEvent, iSetup, reclustered_fatgenjetCollection , iSysType ); //todo

	if(        correctedJet .pt()   < 80 ) continue ; 
	if( fabs( correctedJet .eta() ) > 2.4 ) continue ; 


        re_fatjet_pt . push_back( correctedJet .pt() ) ;  
	re_fatjet_eta. push_back( correctedJet .eta() ) ;  
	re_fatjet_phi. push_back( correctedJet .phi() ) ;  

	re_fatjet_tau21 . push_back( correctedJet.userFloat("NjettinessCA15Puppi:tau2")
				     / 
				     correctedJet.userFloat("NjettinessCA15Puppi:tau1") );

	re_fatjet_tau32 . push_back( correctedJet.userFloat("NjettinessCA15Puppi:tau3")
				     / 
				     correctedJet.userFloat("NjettinessCA15Puppi:tau2") );

        re_fatjet_sdmass_miniaod .push_back( rerun_fatjets -> at(i). userFloat("ca15PFJetsPuppiSoftDropMass")  );


// Not Yet implemented 	pat::Jet originalJet = rerun_fatjets -> at(i) ; 
// Not Yet implemented 	TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
// Not Yet implemented 	auto const & sdSubjetsPuppi = originalJet.subjets("SoftDrop");
// Not Yet implemented 	long n = 0 ;
// Not Yet implemented 	for ( auto const & it : sdSubjetsPuppi ) {
// Not Yet implemented 	  n ++ ;
// Not Yet implemented 	  puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
// Not Yet implemented 	  puppi_softdrop+=puppi_softdrop_subjet;
// Not Yet implemented 	}
// Not Yet implemented 
//
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetToolbox#addSubjets_from_prunning_or_soft : 
	// addSubjets (from prunning or softdrop)
	// 
	//        re_fatjet_sdmass_uncorr  .push_back( puppi_softdrop . M()   );

        re_fatjet_sdmass_uncorr  .push_back( 0 ) ;

      }
    }
    */
    // - - - - MET works - - - - - 

    {
      double met_x = correctedMET.corPx(pat::MET::Type1) + JecUpdatePropagationToMET_x ; 
      double met_y = correctedMET.corPy(pat::MET::Type1) + JecUpdatePropagationToMET_y ; 
      
      eve->MET_[iSys]      = sqrt( met_x * met_x + met_y * met_y );
      eve->MET_phi_[iSys]  = atan2( met_y , met_x   );
    }


    {

      double met_x = correctedMET.corPx(pat::MET::Type1XY) + JecUpdatePropagationToMET_x ; 
      double met_y = correctedMET.corPy(pat::MET::Type1XY) + JecUpdatePropagationToMET_y ; 

      eve->MET_Type1xy_[iSys]      = sqrt( met_x * met_x + met_y * met_y );
      eve->MET_Type1xy_phi_[iSys]  = atan2( met_y , met_x   );		   

    }

    {
      // BBT, 10-04-18
      //double met_x = correctedMET.shiftedP4(pat::MET::NoShift, pat::MET::Type1XY) ; 
      //double met_x = correctedMET.shiftedP4(pat::MET::NoShift, pat::MET::Type1XY) ; 
      //double met_y = pfmet[0].shiftedP4(pat::MET::NoShift, pat::MET::Type1XY) ; 
      //double met_y = pfmet[0].shiftedP4(pat::MET::NoShift, pat::MET::Type1XY) ; 
      double met_x = correctedMET.corPx(pat::MET::Type1XY);
      double met_y = correctedMET.corPy(pat::MET::Type1XY);
      eve->MET_Type1xy_sync_[iSys]      = sqrt( met_x * met_x + met_y * met_y );
      eve->MET_Type1xy_phi_sync_[iSys]  = atan2( met_y , met_x   );		   

    }


    if( false ){ // For test of met correction type : 
      //  std::cout <<"Test satoshi_et: " << correctedMET.pt() << " " << correctedMET.corPt() << " " << correctedMET.corPt(pat::MET::Type1) << " " << correctedMET.corPt(pat::MET::Type1XY) << std::endl ; 
      // std::cout <<"Test satoshi_phi : " << correctedMET.phi() << " " << correctedMET.corPhi() << " " << correctedMET.corPhi(pat::MET::Type1) << " " << correctedMET.corPhi(pat::MET::Type1XY) << std::endl ; 
      //// list of correction types are here : https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/MET.h#L151-L168
    }

    
    if( ! EVENTSYNCMODE ){

    if( jet_pt . size() < 4 ){
      jet_pt. clear() ;
      jet_eta. clear() ;
      jet_phi. clear() ;
      jet_m. clear() ;
      jet_puid. clear() ;
      jet_seed. clear() ; // BBT 10-12-18
      jet_DeepCSV_SF. clear() ; // BBT 11-02-18
      jet_AssociatedGenJet_pt. clear() ;
      jet_AssociatedGenJet_eta. clear() ;
      jet_AssociatedGenJet_phi. clear() ;
      jet_AssociatedGenJet_m. clear() ;
      jet_precore_phi . clear() ; 
      jet_precore_pt . clear() ; 
      jet_genId_vect. clear() ;
      jet_partonflavour_vect. clear() ;
      jet_flavour_vect. clear() ;
      jet_genParentId_vect. clear() ;
      jet_genGrandParentId_vect. clear() ;
    }

    if( jet_pt . size() ){
      b_TheEventHasFourJets_ForAtLeastOneSystematicVariation = true ; 
    }

    }
    
    
    eve->jet_combinedMVABJetTags_[iSys] = jet_combinedMVABJetTags;
    eve->jet_combinedInclusiveSecondaryVertexV2BJetTags_[iSys] = jet_combinedInclusiveSecondaryVertexV2BJetTags;

    eve->jet_combinedMVABJetTags_HIP_[iSys] = jet_combinedMVABJetTags_HIP;
    eve->jet_combinedInclusiveSecondaryVertexV2BJetTags_HIP_[iSys] = jet_combinedInclusiveSecondaryVertexV2BJetTags_HIP;

    eve->jet_DeepCSV_b_  [iSys] = jet_DeepCSV_b;
    eve->jet_DeepCSV_bb_ [iSys] = jet_DeepCSV_bb;

    eve->jet_vtxMass_[iSys]    = jet_vtxMass;
    eve->jet_vtxNtracks_[iSys] = jet_vtxNtracks;
    eve->jet_vtx3DVal_[iSys]   = jet_vtx3DVal;
    eve->jet_vtx3DSig_[iSys]   = jet_vtx3DSig;
    eve->jet_partonflavour_[iSys]          = jet_partonflavour_vect;
    eve->jet_flavour_[iSys]          = jet_flavour_vect;
    eve->jet_genId_[iSys]            = jet_genId_vect;
    eve->jet_genParentId_[iSys]      = jet_genParentId_vect;
    eve->jet_genGrandParentId_[iSys] = jet_genGrandParentId_vect;

    eve->jet_pt_  [iSys]= jet_pt  ;
    eve->jet_phi_ [iSys]= jet_phi ;
    eve->jet_eta_ [iSys]= jet_eta ;
    eve->jet_m_   [iSys]= jet_m   ;

    eve->jet_puid_  [iSys]= jet_puid  ;
    eve->jet_seed_  [iSys]= jet_seed  ; // BBT 10-12-18
    eve->jet_DeepCSV_SF_  [iSys]= jet_DeepCSV_SF  ; // BBT 11-02-18

    eve->jet_precorr_pt_  [iSys]= jet_precore_pt  ;
    eve->jet_precorr_phi_ [iSys]= jet_precore_phi ;

    eve->jet_AssociatedGenJet_pt_[iSys] = jet_AssociatedGenJet_pt;
    eve->jet_AssociatedGenJet_eta_[iSys]= jet_AssociatedGenJet_eta;
    eve->jet_AssociatedGenJet_phi_[iSys]= jet_AssociatedGenJet_phi;
    eve->jet_AssociatedGenJet_m_[iSys]  = jet_AssociatedGenJet_m;

    /* 02-05-19
    eve ->  fatjet_pt         [iSys] =          fatjet_pt            ;     
    eve ->  fatjet_eta	      [iSys] = 	        fatjet_eta	      ;	       
    eve ->  fatjet_phi	      [iSys] = 	        fatjet_phi	      ;	       
    eve ->  fatjet_m  	      [iSys] = 	        fatjet_m  	      ;	       
    eve ->  fatjet_doublebtagging  [iSys] = 	        fatjet_doublebtagging   ;	       
    eve ->  fatjet_nSubjet    [iSys] =          fatjet_nSubjet 	      ;
    eve ->  fatjet_sdmass_miniaod [iSys] =      fatjet_sdmass_miniaod ;    
    eve ->  fatjet_sdmass_uncorr  [iSys] =      fatjet_sdmass_uncorr  ;    
    eve ->  fatjet_tau1	      [iSys] = 	        fatjet_tau1	      ;	       
    eve ->  fatjet_tau2	      [iSys] = 	        fatjet_tau2	      ;	       
    eve ->  fatjet_tau3	      [iSys] = 	        fatjet_tau3	      ;	       
    eve ->  fatjet_tau4	      [iSys] = 	        fatjet_tau4	      ;	       
    eve ->  fatjet_chstau1	      [iSys] =  fatjet_chstau1	      ;
    eve ->  fatjet_chstau2	      [iSys] =  fatjet_chstau2	      ;
    eve ->  fatjet_chstau3	      [iSys] =  fatjet_chstau3	      ;
    eve ->  fatjet_nb1N2 	      [iSys] =  fatjet_nb1N2 	      ;
    eve ->  fatjet_nb1N3 	      [iSys] =  fatjet_nb1N3 	      ;
    eve ->  fatjet_nb2N2 	      [iSys] =  fatjet_nb2N2 	      ;
    eve ->  fatjet_nb2N3 	      [iSys] =  fatjet_nb2N3 	      ;
    eve ->  fatjet_chsprunedmass  [iSys] =      fatjet_chsprunedmass  ;    

    eve -> fatjet_subjet_pt      [iSys] = fatjet_subjet_pt       ;
    eve -> fatjet_subjet_eta     [iSys] = fatjet_subjet_eta      ;
    eve -> fatjet_subjet_phi     [iSys] = fatjet_subjet_phi      ;
    eve -> fatjet_subjet_m       [iSys] = fatjet_subjet_m        ;
    eve -> fatjet_subjet_beepcsv [iSys] = fatjet_subjet_beepcsv  ;
    eve -> fatjet_subjet_csvv2   [iSys] = fatjet_subjet_csvv2    ;


    eve ->  re_fatjet_pt             [iSys] =  re_fatjet_pt            ;     
    eve ->  re_fatjet_eta	     [iSys] =  re_fatjet_eta	      ;	       
    eve ->  re_fatjet_phi	     [iSys] =  re_fatjet_phi	      ;	       
    eve ->  re_fatjet_tau21             [iSys] =  re_fatjet_tau21            ;     
    eve ->  re_fatjet_tau32             [iSys] =  re_fatjet_tau32            ;     
    eve ->  re_fatjet_sdmass_miniaod [iSys] =  re_fatjet_sdmass_miniaod ;    
    eve ->  re_fatjet_sdmass_uncorr  [iSys] =  re_fatjet_sdmass_uncorr  ;    
    */
  } // end loop over systematics


  if( ! EVENTSYNCMODE && ! b_TheEventHasFourJets_ForAtLeastOneSystematicVariation ){
    return ; // Do not record the data.
  }

  //
  // Fill tree if pass full selection
  //
  worldTree->Fill();

  if( EVENTSYNCMODE ){

    selection . EnableInfoDumpForDebug();

    // -----------------------
    // start setting variables --> 


    selection . SetEl_ORTrigger( 2 ,
				 & eve->passHLT_Ele35_WPTight_Gsf_v_ ,
				 & eve->passHLT_Ele28_eta2p1_WPTight_Gsf_HT150_v_) ;

    selection . SetMu_ORTrigger( 2 , 
			      & eve->passHLT_IsoMu24_2p1_v_ ,
			      & eve->passHLT_IsoMu27_v_ 
			      );


    selection . SetElEl_ORTrigger( 2 , 
				   &  eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v_ ,
				   &  eve->passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_  );

    selection . SetElMu_ORTrigger( 4 , 
				   & eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v_,
				   & eve->passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v_,
				   & eve->passHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_,
				   & eve->passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v_ );

    selection . SetMuMu_ORTrigger( 3 , 
				   & eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v_,
				   & eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v_,
				   & eve->passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v_ ) ; 

    selection . SetGoodVtx( & ( eve->GoodFirstPV_ ) );

    selection . SetLepTrigDR( & eve->lepton_dRSingleLepTrig_ );

    selection . SetLeptons( & eve->lepton_pt_, 
			    & eve->lepton_eta_, 
			    & eve->lepton_scEta_, 
			    & eve->lepton_phi_,
			    & eve->lepton_e_,
			    & eve->lepton_charge_, 
			    & eve->lepton_isMuon_, 
			    & eve->lepton_relIso_,
			    & eve->lepton_isLoose_,
			    & eve->lepton_isTight_ );
    selection . SetJets( & eve->jet_pt_  [0] , 
			 & eve->jet_eta_ [0] , 
			 & eve->jet_phi_ [0] , 
			 & eve->jet_m_   [0] , 
			 & eve->jet_combinedInclusiveSecondaryVertexV2BJetTags_[0]  ,
			 & eve->jet_flavour_[0]  );

    selection . SetJetsPTBeforRecorrection( & eve->jet_precorr_pt_  [0] , 
					    & eve->jet_precorr_phi_ [0] ); 

    selection . SetMet( & ( eve->MET_Type1xy_[ 0 ] ) , &( eve->MET_Type1xy_phi_[ 0 ] ) );

    // ---> end setting variables.
    // -----------------------

    
    selection . doEventSelection();

    std::cout.setf(std::ios::fixed);

   // "cout" for event sync. //todo
    if( b_EventSyncFirstEvent ){
      std::cout << "run,lumi,event,is_e,is_mu,is_ee,is_emu,is_mumu,n_jets,n_btags,lep1_pt,lep1_eta,lep1_iso,lep1_pdgId,lep1_idSF,lep1_isoSF,lep2_pt,lep2_eta,lep2_iso,lep2_idSF,lep2_isoSF,lep2_pdgId,jet1_pt,jet1_eta,jet1_phi,jet1_jesSF,jet1_jesSF_up,jet1_jesSF_down,jet1_jesSF_PileUpDataMC_down,jet1_jesSF_RelativeFSR_up,jet1_jerSF_nominal,jet1_csv,jet1_PUJetId,jet1_PUJetDiscriminant,jet2_pt,jet2_eta,jet2_phi, jet2_jesSF,jet2_jesSF_up,jet2_jesSF_down,jet2_jesSF_PileUpDataMC_down,jet2_jesSF_RelativeFSR_up,jet2_jerSF_nominal,jet2_csv,jet2_PUJetId,jet2_PUJetDiscriminant,MET_pt,MET_phi,mll,ttHFCategory,ttHFGenFilterTag,n_interactions,puWeight,csvSF,csvSF_lf_up,csvSF_hf_down,csvSF_cErr1_down" << std::endl;
      b_EventSyncFirstEvent = false ;
    }

    std::cout << eve->run_ << "," ;
    std::cout <<eve->lumi_ << "," ;
    std::cout <<eve->evt_  << "," ;
    
    //      is_e, is_mu, is_ee, is_emu, is_mumu,
    std::cout <<(selection . PassSingleElCh() ?1:0) << "," ;
    std::cout <<(selection . PassSingleMuCh() ?1:0) << "," ;
    std::cout <<(selection . PassElEl()       ?1:0) << "," ;
    std::cout <<(selection . PassElMu()       ?1:0) << "," ;
    std::cout <<(selection . PassMuMu()       ?1:0) << "," ;


    std::cout << selection.jets().size() << "," ;
    std::cout << selection.bjets().size() << ",";

    if( selection.looseLeptons().size() >=1 ){
      long pdgid = 
	( selection.looseLeptonsIsMuon().at(0) == 1 ? 13 : 11 )
	*
	( selection.looseLeptonsCharge().at(0) > 0 ? -1 : +1 ) ;
      std::cout<< std::setprecision(4) << selection.looseLeptons().at(0)->Pt()<< "," ;
      std::cout<< std::setprecision(4) << selection.looseLeptons().at(0)->Eta()<< "," ;
      std::cout<< std::setprecision(4) << selection.looseLeptonsRelIso().at(0)<< "," ;
      std::cout << pdgid << "," ;

      float lep1_ID_SF = 0 ;  //todo. to be implemented.
      std::cout <<  lep1_ID_SF << "," ; 

      float lep1_Iso_SF = 0 ;  //todo. to be implemented.
      std::cout <<  lep1_Iso_SF << "," ; 
    }else{
      std::cout << "-1,-1,-1,-1,-1,-1," ; // pt/eta/iso/pdgif/idSF/isoSF.
    }

    if( selection.looseLeptons().size() >=2 ){
      long pdgid = 
	( selection.looseLeptonsIsMuon().at(1) == 1 ? 13 : 11 )
	*
	( selection.looseLeptonsCharge().at(1) > 0 ? -1 : +1 ) ;
      std::cout<< std::setprecision(4) << selection.looseLeptons().at(1)->Pt()<< "," ;
      std::cout<< std::setprecision(4) << selection.looseLeptons().at(1)->Eta()<< "," ;
      std::cout<< std::setprecision(4) << selection.looseLeptonsRelIso().at(1)<< "," ;

      float lep2_ID_SF = 0 ;  //todo. to be implemented.
      std::cout <<  lep2_ID_SF << "," ; 

      float lep2_Iso_SF = 0 ;  //todo. to be implemented.
      std::cout <<  lep2_Iso_SF << "," ; 

      std::cout << pdgid << "," ; // Note. For Lepton-2, PDFID comes after ID_SF and iso_SF.

    }else{
      std::cout << "-1,-1,-1,-1,-1,-1," ; // pt/eta/iso/pdgif/idSF/isoSF.
    }
    
    {
      bool nJet_ge_one = selection.jets().size() >=1;
      bool nJet_ge_two = selection.jets().size() >=2;

      double JEC[2]     ={ -1,  -1 };
      double JECup[2]   ={ -1,  -1 };      
      double JECdown[2] ={ -1,  -1 };
      double JEC_PileupData_down[2] ={ -1,  -1 };
      double JEC_RelativeSFR_up[2]   ={ -1,  -1 };      
      double JER[2]     ={ -1,  -1 };
      int    PileupID[2] ={ -1,  -1 };

      const bool  doJES = true;
      const bool  dontJER = false;

      std::vector<pat::Jet> rawJets = miniAODhelper.GetUncorrectedJets( *pfjets );
      std::vector<pat::Jet> jet_JESNOMI =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection ,sysType::NA     , doJES, dontJER );
      std::vector<pat::Jet> jet_JESUP   =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection ,sysType::JESup  , doJES, dontJER );
      std::vector<pat::Jet> jet_JESDOWN =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection ,sysType::JESdown, doJES, dontJER );

      std::vector<pat::Jet> jet_JES_PileupData_DOWN =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection , sysType::JESPileUpDataMCdown, doJES, dontJER );
      std::vector<pat::Jet> jet_JES_RelativeSFR_UP  =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection ,sysType::JESRelativeFSRup    , doJES, dontJER );

      std::vector<pat::Jet> jet_JER  =  miniAODhelper.GetCorrectedJets(rawJets, iEvent, iSetup, genjetCollection ,sysType::NA   , false, true ); // <- Only JER

      for( unsigned int iJet = 0 ; iJet < 2 && iJet < selection.jets().size() ; iJet ++ ){

	const double eta = selection.jets().at(iJet)->Eta();
	const double phi = selection.jets().at(iJet)->Phi();

	// Loop for JEC_nominal
	for( unsigned int idxJet = 0 ; idxJet < rawJets.size(); idxJet ++ ){

	  pat::Jet * iRawJet = & rawJets.at( idxJet );
	  double d_eta =       eta -  iRawJet->eta();
	  double d_phi = fabs( phi -  iRawJet->phi() ) ; 
	  d_phi = ( d_phi < M_PI ) ? d_phi : 2 * M_PI - d_phi ; 

	  if(  d_eta*d_eta + d_phi*d_phi < 0.01 * 0.01 ){ // matching btw Raw and Corrected (physics) jet.
	    JEC    [iJet] = jet_JESNOMI.at( idxJet ).pt() / iRawJet->pt();
	    JECup  [iJet] = jet_JESUP  .at( idxJet ).pt() / iRawJet->pt();
	    JECdown[iJet] = jet_JESDOWN.at( idxJet ).pt() / iRawJet->pt();
	    JEC_PileupData_down[iJet] = jet_JES_PileupData_DOWN.at(idxJet).pt() / iRawJet->pt() ; // with respect to raw ? 
	    JEC_RelativeSFR_up [iJet] = jet_JES_RelativeSFR_UP .at(idxJet).pt() / iRawJet->pt() ; // with respect to raw ? 
	    JER                [iJet] = jet_JER                .at(idxJet).pt() / iRawJet->pt() ; // with respect to raw ? 

	    PileupID[iJet] = iRawJet  -> userInt("pileupJetId:fullId");

	  }
	}

      }




      std::cout<< std::setprecision(4) << ( nJet_ge_one ? selection.jets().at(0)->Pt() : -1 )<< "," ;
      std::cout<< std::setprecision(4) << ( nJet_ge_one ? selection.jets().at(0)->Eta() : -1 )<< "," ;
      std::cout<< std::setprecision(4) << ( nJet_ge_one ? selection.jets().at(0)->Phi() : -1 )<< "," ;

      std::cout << JEC[0] << ","
		<< JECup[0] << ","
		<< JECdown[0] << "," 
		<< JEC_PileupData_down[0] << "," 
		<< JEC_RelativeSFR_up[0] << "," 
		<< JER[0] << "," ;

      std::cout<< std::setprecision(4) << ( nJet_ge_one ? selection.jetsBdiscriminant().at(0) : -1 )<< "," ;

      std::cout << PileupID[0] << ","
		<< PileupID[0] << "," ;

      // - - - - Second Jet - - - 

      std::cout<< std::setprecision(4) << ( nJet_ge_two ? selection.jets().at(1)->Pt() : -1 )<< "," ;
      std::cout<< std::setprecision(4) << ( nJet_ge_two ? selection.jets().at(1)->Eta() : -1 )<< "," ;
      std::cout<< std::setprecision(4) << ( nJet_ge_two ? selection.jets().at(1)->Phi() : -1 )<< "," ;

      std::cout << JEC[1] << ","
		<< JECup[1] << ","
		<< JECdown[1] << "," 
		<< JEC_PileupData_down[1] << "," 
		<< JEC_RelativeSFR_up[1] << "," 
		<< JER[1] << "," ;

      std::cout<< std::setprecision(4) << ( nJet_ge_two ? selection.jetsBdiscriminant().at(1) : -1 )<< "," ;

      std::cout << PileupID[1] << ","
		<< PileupID[1] << "," ;

    }


    
    std::cout<< std::setprecision(4) << selection . metAbs() << "," ;
    std::cout<< std::setprecision(4) << selection . metPhi() << "," ;

    double mll = -1 ; 
    if(   selection.looseLeptons().size() >=2 ){
      mll = 
	( * (selection.looseLeptons().at(0)) 
	+
	  * (selection.looseLeptons().at(1)) ) . M();
    }
    std::cout<< std::setprecision(4) << mll << "," ;

    double ttHFenFilterTag = -1 ; 

    if( isMC ){
    std::cout << eve->additionalJetEventId_ <<",";
    std::cout << ttHFenFilterTag << ",";
    std::cout << eve->numTruePV_ << ",";
    std::cout << scalefactors.get_pu_wgt( eve -> numTruePV_ ) << "," ;    // PUWeight,
    }else{
    std::cout << -1<<",";
    std::cout << -1<<",";
    std::cout << -1 << ",";
    std::cout << -1<< "," ;    // PUWeight,
    }
    /*
    double bWeight = -1 ;
    double bWeight_lf_up   = -1 ;
    double bWeight_hf_down = -1 ;
    double bWeight_cErr1   = -1 ;
    
    if( isMC ){
      int iSYS = 0 ; 
      double dummy = - 1 ;
      bWeight = scalefactors.get_csv_wgt( & selection , iSYS,  dummy , dummy , dummy );

      bWeight_lf_up   = scalefactors.get_csv_wgt( & selection ,  9 ,  dummy , dummy , dummy ); //   case 9 : iSysType = sysType::CSVLFup;      
      bWeight_hf_down = scalefactors.get_csv_wgt( & selection , 12 ,  dummy , dummy , dummy ); //   case 12: iSysType = sysType::CSVHFdown;  

    }

    std::cout << bWeight <<",";
    std::cout << bWeight_lf_up   <<",";
    std::cout << bWeight_hf_down <<",";
    std::cout << bWeight_cErr1   <<",";
    */


// 
//     double triggerSF =( ! isMC ?
// 			1 :
// 			scalefactors.get_TrigMuSF( & selection )
// 			*
// 			scalefactors.get_TrigElSF( & selection )
// 			) ;
//     std::cout << triggerSF <<",";
// 
//     double lepIDSF =  ( ! isMC ? 1 : 
// 			scalefactors.getTightMuon_IDSF( & selection )
// 			* 
// 			scalefactors.getTightElectron_IDSF( & selection )
// 			);
//     double lepISOSF =  ( ! isMC ? 1 : 
// 			 scalefactors.getTightElectron_RecoSF( & selection ) 
// 			*
// 			 scalefactors.getTightMuon_IsoSF( & selection ) 
// 			) ; 
//     std::cout << lepIDSF <<","<< lepISOSF <<",";



//    if( isMC ){
//    std::cout << eve->weight_q2_upup_ <<",";
//    std::cout << eve->weight_q2_downdown_ <<",";
//    std::cout << eve-> weight_PDF_NNPDF30NLO_up_ <<",";
//    std::cout << eve-> weight_PDF_NNPDF30NLO_down_ ;
//    }else{
//    std::cout << 1 <<",";
//    std::cout << 1 <<",";
//    std::cout << 1 <<",";
//    std::cout << 1  ;
//    }
    
    std::cout << std::endl ;
  }

}


void YggdrasilTreeMaker::DumpDecay (  const reco::Candidate * c ,
			      int depth ){

  for( unsigned int i = 0 ; i < c->numberOfDaughters(); i ++ ){
    
    for( int _k = 0 ; _k < depth ; _k ++ ){
      std::cout <<"\t";
    }
    std::cout << c->daughter(i) ->pdgId() << " ("
	      << c->daughter(i) ->pt()    << ", "
	      << c->daughter(i) ->eta()   <<  ")" << std::endl ; 
    DumpDecay( c->daughter(i) , depth + 1  );

  }

}


const reco::Candidate * YggdrasilTreeMaker::TraceBackToJustAfterBirth( const reco::Candidate * particle ){

  for ( unsigned int i = 0 ; i <  particle -> numberOfMothers(); i++ ){
    if( particle -> mother( i ) -> pdgId()  ==  particle -> pdgId() ){
      return TraceBackToJustAfterBirth( particle -> mother (i) ) ; 
    }
  }

  return particle ; 
}


bool YggdrasilTreeMaker::checkIfRegisterd( const reco::Candidate * candidate , std::vector< const reco::Candidate * > list ){

  for ( std::vector< const reco::Candidate * >::iterator it = list.begin() ; 
	it != list.end() ; 
	it ++ ){
    if( candidate == * it  ) return true  ;
  }
  
  return false ;

} 



// ------------ method called once each job just before starting event loop  ------------
void 
YggdrasilTreeMaker::beginJob()
{
  edm::Service<TFileService> fileService;
  if(!fileService) throw edm::Exception(edm::errors::Configuration, "TFileService is not registered in cfg file");
    
  h_ttbarId_ = fileService->make<TH1F>("ttbarId", "ttbarId", 260, 0, 260);
  h_ttbarAdditionalJetId_ = fileService->make<TH1F>("ttbarAdditionalJetId", "ttbarAdditionalJetId", 60, 0, 60);

  h_event = fileService->make<TH1D>("event", "event", 4, 0 , 4 ) ;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
YggdrasilTreeMaker::endJob() 
{

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

  std::cout <<"[debug message] YggdrasilTreeMaker::beginRun() is called." << std::endl ; 
  miniAODhelper.UpdateJetCorrectorUncertainties( iSetup );
  miniAODhelper_Puppi.UpdateJetCorrectorUncertainties( iSetup );
  miniAODhelper_fatjet.SetAK8JetCorrectorUncertainty( iSetup ); // Ak8

  std::cout << "[debug message] YggdrasilTreeMaker::beginRun() has finished JetCorrectorWork for all jet collections" << std::endl;


  miniAODhelper.SetJER_SF_Tool( iSetup );
  std::cout <<"[debug message] YggdrasilTreeMaker::beginRun() has finished JER SF work for resolved jets" << std::endl ; 
  miniAODhelper_Puppi.SetJER_SF_Tool( iSetup );
  std::cout <<"[debug message] YggdrasilTreeMaker::beginRun() has finished JER SF work for PUPPI jets" << std::endl ; 
  miniAODhelper_fatjet.SetJER_SF_Tool( iSetup );
  std::cout <<"[debug message] YggdrasilTreeMaker::beginRun() has finished JER SF work for all jets" << std::endl ; 
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
