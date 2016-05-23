import FWCore.ParameterSet.Config as cms

#isMC=True
isMC=True


isPUPPI=False
#--> if isPUPPI=true, change the code too.

enableJECFromLocalDB=False

process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#### caution: use the correct global tag for MC or Data 
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

if isMC:
    process.GlobalTag.globaltag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1'
else :
    process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'




# Load the producer for MVA IDs. Make sure it is also added to the sequence!
process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )



if enableJECFromLocalDB :
    print "JEC is applied with LOCAL DB."

    import sys
    import os.path
   

    #
    # taken from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecSqliteFile
    #
    process.load("CondCore.DBCommon.CondDBCommon_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
      cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_DATA_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            ),
      ), 
#      connect = cms.string( 'sqlite:' + os.environ['CMSSW_BASE'] + '/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/Fall15_25nsV2_DATA.db' )
      connect = cms.string( 'sqlite:Fall15_25nsV2_DATA.db' )
                            )
      ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')



from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

if isPUPPI :
    process.ak4PFCHSL1Fastjet = cms.ESProducer(
        'L1FastjetCorrectionESProducer',
        level       = cms.string('L1FastJet'),
        algorithm   = cms.string('AK4PFPuppi'),
        srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
    )
    process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFPuppi' )
    process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFPuppi' )

else :
    process.ak4PFCHSL1Fastjet = cms.ESProducer(
        'L1FastjetCorrectionESProducer',
        level       = cms.string('L1FastJet'),
        algorithm   = cms.string('AK4PFchs'),
        srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
    )
    process.ak4PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK4PFchs' )
    process.ak4PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK4PFchs' )


process.ak4PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring(
	'ak4PFCHSL1Fastjet', 
        'ak4PFchsL2Relative', 
        'ak4PFchsL3Absolute')
)

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/00DF0A73-17C2-E511-B086-E41D2D08DE30.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/020D0AF1-4BC2-E511-BDFC-0026B95ADB18.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/02837459-03C2-E511-8EA2-002590A887AC.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/0400D1B1-4CC2-E511-AD14-AC853DA06A1A.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/0634456A-08C2-E511-A0C1-001E6739722E.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/0A5A5364-03C2-E511-9352-00259048AE00.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/0A5C7429-0CC2-E511-BC6B-002590D9D8B8.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/0AE230FC-4BC2-E511-B35C-000F530E4644.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/0C5BB11A-E2C1-E511-8C02-002590A831B6.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/10D0AF2E-33C2-E511-BCE5-001E67398683.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/14C8A87C-DCC1-E511-BE67-0CC47A0AD704.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/14F4047E-0CC2-E511-B1CE-002590A3C970.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/16467C82-D8C1-E511-9546-001E67E71417.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/167A53AB-08C2-E511-A58B-001E67397F2B.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/18897870-17C2-E511-9DC5-E41D2D08DCA0.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/18CCDFE4-02C2-E511-B0FA-002590A37114.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/20A0646D-D8C1-E511-9587-001E67396A09.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/2235F83F-DBC1-E511-A83F-001E673986B5.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/22761DF1-E1C1-E511-8F86-003048CB8584.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/244D8186-D8C1-E511-B684-001E67E71BF5.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/26CAB9BA-CEC1-E511-8B0B-001E67398E49.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/2A18493A-04C2-E511-935C-002590D9D8C2.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/2C7F8762-0CC2-E511-9F77-001E673981C4.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/2CA9FB72-17C2-E511-A3E3-0025907DCA7E.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/2CCB8F4A-07C2-E511-A3C6-002590A3C970.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/2EF6150C-1AC2-E511-A67B-002590A83192.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/324D4BDC-DCC1-E511-925D-00259029E922.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/32D38D09-E2C1-E511-986A-001E67397BC5.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/32E942A0-4DC2-E511-AD4A-001E674FCBF7.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/343CDCDC-E1C1-E511-80C4-00259029E91C.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/36171880-02C2-E511-9B6E-0CC47A0AD630.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/363584E8-4BC2-E511-9DE3-008CFAF060D8.root',
'/store/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/3660CB69-35C2-E511-86B8-002590A83354.root'

#        '/store/user/puigh/TTHSync/ttjets_phys14_20bx25_withfatjets_v2.root'
            #'/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/08B36E8F-5E7F-E411-9D5A-002590200AE4.root'
            #'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root'
            #'/store/mc/Spring14miniaod/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1E4F9BDC-3E1E-E411-A56C-001E67396EAA.root'
            #'/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/004C6DA7-FB03-E411-96BD-0025905A497A.root'
            )
)


########## PUPPI Muon Isolation 
process.PUPPIMuonRelIso = cms.EDProducer('PuppiLeptonIsolation'
                                         , srcLepton = cms.string( "slimmedMuons" )
                                         , dR = cms.double( 0.4 )
                                         , mixFraction = cms.double( 0.5 )
                                         , configuration = cms.string( "##" )
)



###############
#### tt+X
###############
# Setting input particle collections to be used by the tools
genJetCollection = 'ak4GenJetsCustom'
genParticleCollection = 'prunedGenParticles'
genJetInputParticleCollection = 'packedGenParticles'

## producing a subset of particles to be used for jet clustering
from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
process.genParticlesForJetsNoNu = genParticlesForJetsNoNu.clone(
	src = genJetInputParticleCollection
)

# Supplies PDG ID to real name resolution of MC particles
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

# Producing own jets for testing purposes
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsCustom = ak4GenJets.clone(
    src = 'genParticlesForJetsNoNu',
#    src = genJetInputParticleCollection,
    rParam = cms.double(0.4),
    jetAlgorithm = cms.string("AntiKt")
)

# Ghost particle collection used for Hadron-Jet association 
# MUST use proper input particle collection
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
    particles = genParticleCollection
)

# Input particle collection for matching to gen jets (partons + leptons) 
# MUST use use proper input jet collection: the jets to which hadrons should be associated
# rParam and jetAlgorithm MUST match those used for jets to be associated with hadrons
# More details on the tool: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools#New_jet_flavour_definition
from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
    rParam = cms.double(0.4),
    jetAlgorithm = cms.string("AntiKt")
)


# Plugin for analysing B hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
    genParticles = genParticleCollection
)

# Plugin for analysing C hadrons
# MUST use the same particle collection as in selectedHadronsAndPartons
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenCHadron
process.matchGenCHadron = matchGenCHadron.clone(
    genParticles = genParticleCollection
)

## Producer for ttbar categorisation ID
# MUST use same genJetCollection as used for tools above
from TopQuarkAnalysis.TopTools.GenTtbarCategorizer_cfi import categorizeGenTtbar
process.categorizeGenTtbar = categorizeGenTtbar.clone(
    genJetPtMin = 20.,
    genJetAbsEtaMax = 2.4,
    genJets = genJetCollection,
)


# , MuonIso = cms.string( "PUPPIMUONISO" )

if isMC :
    if isPUPPI :
        process.ttHTreeMaker = cms.EDAnalyzer('YggdrasilTreeMaker',
                                          isMC    =  cms.string("MC"),
                                          jetPU = cms.string( "PUPPI" )
                                         , MuonIso = cms.string( "PUPPIMUONISO" )
                                          )
    else:
        process.ttHTreeMaker = cms.EDAnalyzer('YggdrasilTreeMaker',
                                          isMC    =  cms.string("MC"),
                                          jetPU = cms.string( "CHS" )
                                         , MuonIso = cms.string( "PUPPIMUONISO" )
                                          )

else :
    if isPUPPI :
        process.ttHTreeMaker = cms.EDAnalyzer('YggdrasilTreeMaker',
                                          isMC    =  cms.string("data"),
                                          jetPU = cms.string( "PUPPI" )
                                         , MuonIso = cms.string( "PUPPIMUONISO" )
                                          )
    else:
        process.ttHTreeMaker = cms.EDAnalyzer('YggdrasilTreeMaker',
                                          isMC    =  cms.string("data"),
                                          jetPU = cms.string( "CHS" )
                                         , MuonIso = cms.string( "PUPPIMUONISO" )
                                          )
        
    
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('yggdrasil_treeMaker.root')
)

process.p = cms.Path(process.electronMVAValueMapProducer * process.PUPPIMuonRelIso * process.ttHTreeMaker)
