import FWCore.ParameterSet.Config as cms


isPUPPI=False
#--> if isPUPPI=true, change the code too.


#isMC=True
isMC=True

# 
isTTBARMC=False

isGridJob=False
# isGridJob=True

genjetInputTag = cms.InputTag("slimmedGenJets","","")
#genjetInputTag = cms.InputTag("ak4GenJetsReproduced","","")
#genjetInputTag = cms.InputTag("ak4GenJetsWithChargedLepFromTop","","")


enableJECFromLocalDB=True

# - - - - - - - - - - - - - - - - - - - - 
# Special option for Morind17 analysis
#  This is a flag used to apply dedicated JEC for each data set.
#  The placeholder will be replaced by crab job make script.
isPeriodBCD=False
isPeriodEF1=False
isPeriodF2G=False
isPeriodH=False
# - - - - - - - - - - - - - - - - - - - -

# Switch to perform lumi-mask inside this python script.
# It is preferable to set this option False as default,
#   because, when we submit job to grid, we also set the file too.
ForDebugAndEventSync_EnableLumiMaskByHand=False


process = cms.Process("MAOD")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#### caution: use the correct global tag for MC or Data 
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')


# Update global tag based on : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions?rev=568
if isMC:
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8'
else :
    process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'



# Load the producer for MVA IDs. Make sure it is also added to the sequence!
##process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    )



import sys
import os.path

JecLocalDataBaseName = \
    'Summer16_23Sep2016BCDV3_DATA' if isPeriodBCD else \
    'Summer16_23Sep2016EFV3_DATA'  if isPeriodEF1 else \
    'Summer16_23Sep2016GV3_DATA'   if isPeriodF2G else \
    'Summer16_23Sep2016HV3_DATA'   if isPeriodH   else 'Summer16_23Sep2016V3_MC'

JecDBPathPrefix = 'sqlite://.' if isGridJob else 'sqlite:///'+os.environ.get('CMSSW_BASE') 
# This switch is needed because the variable CMSSW_BASE remains the same as local job (directory where you do "cmsenv") when the job runs on the grid.



if enableJECFromLocalDB :


    print "JEC is applied with LOCAL DB. -- " + JecLocalDataBaseName
                           
    process.GlobalTag.toGet.append(
        cms.PSet(
            connect = cms.string( JecDBPathPrefix +'/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/' + JecLocalDataBaseName +'.db' ),
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_'+JecLocalDataBaseName+'_AK4PFchs'),
            label  = cms.untracked.string('AK4PFchs')
            )
        )
    process.GlobalTag.toGet.append(
        cms.PSet(
            connect = cms.string( JecDBPathPrefix + '/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/' + JecLocalDataBaseName +'.db' ),
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_'+JecLocalDataBaseName+'_AK8PFchs'),
            label  = cms.untracked.string('AK8PFchs')
            )
        )
    
    #  line for PUPPI for temporal.
    process.GlobalTag.toGet.append(
        cms.PSet(
            connect = cms.string( JecDBPathPrefix +'/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/' + 'Summer16_25nsV5_MC'+'.db' ),
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_'+'Summer16_25nsV5_MC'+'_AK4PFPuppi'),
            label  = cms.untracked.string('AK4PFPuppi')
            )
        )



# Set up JetCorrections chain to be used in MiniAODHelper
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
# process.ak4PFCHSL1Fastjet = cms.ESProducer(
#  'L1FastjetCorrectionESProducer',
#  level = cms.string('L1FastJet'),
#  algorithm = cms.string('AK4PFchs'),
#  srcRho = cms.InputTag( 'fixedGridRhoFastjetAll' )
#  )
# process.ak4PFchsL2RelativeCorrector = ak4CaloL2RelativeCorrector.clone( algorithm = 'AK4PFchs' )
# process.ak4PFchsL3AbsoluteCorrector = ak4CaloL3AbsoluteCorrector.clone( algorithm = 'AK4PFchs' )
# process.ak4PFchsResidualCorrector   = ak4CaloResidualCorrector.clone( algorithm = 'AK4PFchs' )
# process.ak4PFchsL1L2L3 = cms.ESProducer("ChainedJetCorrectorProducer",
#  correctors = cms.VInputTag(
#    'ak4PFCHSL1FastjetCorrector',
#    'ak4PFchsL2RelativeCorrector',
#    'ak4PFchsL3AbsoluteCorrector')
# )
# if not isMC :
#     process.ak4PFchsL1L2L3.correctors.append('ak4PFchsResidual') # add residual JEC for data
# 
# 
# process.ak8PFCHSL1Fastjet = cms.ESProducer(
#  'L1FastjetCorrectionESProducer',
#  level = cms.string('L1FastJet'),
#  algorithm = cms.string('AK8PFchs'),
#  srcRho = cms.InputTag( 'fixedGridRhoFastjetAll' )
#  )
# process.ak8PFchsL2Relative = ak4CaloL2Relative.clone( algorithm = 'AK8PFchs' )
# process.ak8PFchsL3Absolute = ak4CaloL3Absolute.clone( algorithm = 'AK8PFchs' )
# process.ak8PFchsResidual = ak4CaloResidual.clone( algorithm = 'AK8PFchs' )
# process.ak8PFchsL1L2L3 = cms.ESProducer("JetCorrectionESChain",
#  correctors = cms.vstring(
#    'ak8PFCHSL1Fastjet',
#    'ak8PFchsL2Relative',
#    'ak8PFchsL3Absolute')
# )
# 
# if not isMC :
#  process.ak8PFchsL1L2L3.correctors.append('ak8PFchsResidual') # add residual JEC for data
# 


process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        'file:///uscms/home/satoshi/temporal_strage/0CB1DE33-60BF-E611-B520-0025905A4964.root'
#        ' /store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-600to800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/120000/0CB1DE33-60BF-E611-B520-0025905A4964.root'
#        '/store/user/puigh/TTHSync/ttjets_phys14_20bx25_withfatjets_v2.root'
            #'/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/08B36E8F-5E7F-E411-9D5A-002590200AE4.root'
            #'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root'
            #'/store/mc/Spring14miniaod/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/1E4F9BDC-3E1E-E411-A56C-001E67396EAA.root'
            #'/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v2/00000/004C6DA7-FB03-E411-96BD-0025905A497A.root'
            )
)





# Egamma energy smearing for mc
if isMC :
    process.load("Configuration.StandardSequences.Services_cff")
    process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                       calibratedPatElectrons = cms.PSet(
            initialSeed = cms.untracked.uint32(81),
            engineName  = cms.untracked.string("TRandom3")
        )
    )
    process.load("EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi")
    from EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi import files
    process.calibratedPatElectrons.isMC = cms.bool(True)
    process.calibratedPatElectrons.correctionFile = cms.string(files["Moriond2017_JEC"])
  


###############
### GenJet production from ChargedLeptonVetoedGenParticles

if not isMC :
    if ForDebugAndEventSync_EnableLumiMaskByHand :
        import sys
        import os.path
        import FWCore.PythonUtilities.LumiList as Lumilist
        process.source.lumisToProcess = LumiList.LumiList(filename =
                                                          os.environ.get('CMSSW_BASE')+'/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt'
                                                          ).getVLuminosityBlockRange()



if isMC :
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets

    setattr( process , 'myGenParticlesForJets' ,
             cms.EDProducer("InputGenJetsParticleSelector",
                            src = cms.InputTag("packedGenParticles"),
                            ignoreParticleIDs = cms.vuint32(
                1000022,
                1000012, 1000014, 1000016,
                2000012, 2000014, 2000016,
                1000039, 5100039,
                4000012, 4000014, 4000016,
                9900012, 9900014, 9900016,
                39, 12,14,16),
                            partonicFinalState = cms.bool(False),
                            excludeResonances = cms.bool(False),
                            excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
                            tausAsJets = cms.bool(False)
                            ) )
    
    genJetInputParticleCollection = 'myGenParticlesForJets' 
    
    from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
    process.ak4GenJetsReproduced = ak4GenJets.clone(
        src = genJetInputParticleCollection,
        rParam = cms.double(0.4),
        jetAlgorithm = cms.string("AntiKt")
        )
    
    
    process.GenParticleWithoutChargedLeptonFropTop = cms.EDProducer('GenParticleTopOriginChargedleptonFilter',
                                                                    PackedGenParticle = cms.InputTag( "packedGenParticles" )
                                                                    , PrunedGenParticle = cms.InputTag( "prunedGenParticles" )
                                                                    )
    
    process.myGenParticlesWithChargedLeptonFromTopForJet = process.myGenParticlesForJets . clone(
        src = cms.InputTag("GenParticleWithoutChargedLeptonFropTop","TopOriginChargedleptonFilteredGenParticle","")
        )
    process.ak4GenJetsWithChargedLepFromTop = ak4GenJets.clone(
        src = cms.InputTag( 'myGenParticlesWithChargedLeptonFromTopForJet' ),
        rParam = cms.double(0.4),
        jetAlgorithm = cms.string("AntiKt")
        )

    if isTTBARMC :
        process.load("TopQuarkAnalysis.TopEventProducers.sequences.ttGenEvent_cff")
        process.initSubset.src=cms.InputTag("prunedGenParticles")
        process.decaySubset.src=cms.InputTag("prunedGenParticles")
        process.decaySubset.fillMode = cms.string("kStable")
          ## define fill mode. The following modes are available:
          ## 'kStable' : status 2 equivalents (after parton shower) are
          ##             calculated and saved (as status 2 particles)
          ## 'kME'     : status 3 particles (from matrix element, before
          ##             parton shower) are saved (as status 3 particles)


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




#######
## it's included in CMSSW(>80X) by default now
################# hip mitigation
##process.load("Configuration.StandardSequences.MagneticField_cff")
##process.load("Configuration.Geometry.GeometryRecoDB_cff")
##from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
##if isMC :
##    updateJetCollection(
##        process,
##        jetSource = cms.InputTag('slimmedJets'),
##        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute' ]),'None' ),
##        btagDiscriminators = ['pfCombinedMVAV2BJetTags' , 'pfCombinedInclusiveSecondaryVertexV2BJetTags'],
##        runIVF=True,
##        btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
##        )
##else :
##    updateJetCollection(
##        process,
##        jetSource = cms.InputTag('slimmedJets'),
##        jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']),'None' ),
##        btagDiscriminators = ['pfCombinedMVAV2BJetTags' , 'pfCombinedInclusiveSecondaryVertexV2BJetTags'],
##        runIVF=True,
##        btagPrefix = 'new' # optional, in case interested in accessing both the old and new discriminator values
##        )
##


################################################
# Instruction from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription?rev=59
#   to make met(puppi met) with correct calibration
# 
# Those create met object which can be obtained by 
#     cms.InputTag("slimmedMETs","","YourProcessName")
# and 
# cms.InputTag("slimmedMETsPuppi","","YourProcessName")


from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData= not isMC , 
#(relative path from src/)  jecUncFile = ( 'filepath' ),
                           )

process.load("Configuration.StandardSequences.GeometryDB_cff")
#escape 
#escape from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
#escape makePuppiesFromMiniAOD( process, True );
#escape runMetCorAndUncFromMiniAOD(process,
#escape                            isData = not isMC ,
#escape                            metType="Puppi",
#escape                            pfCandColl=cms.InputTag("puppiForMET"),
#escape                            recoMetFromPFCs=True,
#escape                            jetFlavor="AK4PFPuppi",
#escape                            postfix="Puppi",
#escape #                           jecUncFile = ( '' ),
#escape                            )
#escape process.puppiNoLep.useExistingWeights = False
#escape process.puppi.useExistingWeights = False
#escape 


if isMC :
    if isPUPPI :
        process.ttHTreeMaker = cms.EDAnalyzer('YggdrasilTreeMaker',
                                              genjet =  genjetInputTag,
                                              inputfiletype   =  cms.string("TTbarMC") if isTTBARMC else cms.string("MC") ,
                                              jetPU = cms.string( "PUPPI" )
                                          )
    else:
        process.ttHTreeMaker = cms.EDAnalyzer('YggdrasilTreeMaker',
                                              genjet =  genjetInputTag,
                                              inputfiletype   =  cms.string("TTbarMC") if isTTBARMC else cms.string("MC") ,
                                          jetPU = cms.string( "CHS" )
                                          )

else :
    if isPUPPI :
        process.ttHTreeMaker = cms.EDAnalyzer('YggdrasilTreeMaker',
                                              genjet =  genjetInputTag,
                                          inputfiletype    =  cms.string("data"),
                                          jetPU = cms.string( "PUPPI" )
                                          )
    else:
        process.ttHTreeMaker = cms.EDAnalyzer('YggdrasilTreeMaker',
                                              genjet =  genjetInputTag,
                                          inputfiletype    =  cms.string("data"),
                                          jetPU = cms.string( "CHS" )
                                          )
        
    
process.TFileService = cms.Service("TFileService",
	fileName = cms.string('yggdrasil_treeMaker.root')
)

process.PUPPIMuonRelIso = cms.EDProducer('PuppiLeptonIsolation'
                                         , srcLepton =cms.string('slimmedMuons')
                                         , dR = cms.double( 0.4 ) 
                                         , mixFraction = cms.double( 0.5 ) 
                                         , configuration = cms.string( "##" )
                                         )




# Electron Energy Regression
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)
process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

process.load('RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi')
process.load('RecoEgamma.PhotonIdentification.PhotonMVAValueMapProducer_cfi')
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
process.photonMVAValueMapProducer.srcMiniAOD   = cms.InputTag('slimmedPhotons')


if isMC : 
    print "[debuyg message ]  process process.calibratedPatElectrons is added."
    process.p = cms.Path(
        process.regressionApplication *   #<- Electron energy regression 
        process.calibratedPatElectrons * # <- Electron Smearing. MC only.
#        process.egmPhotonIDSequence * process.puppiMETSequence * process.fullPatMetSequencePuppi *
        process.GenParticleWithoutChargedLeptonFropTop * process.myGenParticlesWithChargedLeptonFromTopForJet * process.ak4GenJetsWithChargedLepFromTop *  
        process.PUPPIMuonRelIso * process.ttHTreeMaker)
#        process.PUPPIMuonRelIso * process.electronMVAValueMapProducer * process.ttHTreeMaker)
else :
    process.p = cms.Path(
        process.regressionApplication *   #<- Electron energy regression 
#        process.egmPhotonIDSequence * process.puppiMETSequence * process.fullPatMetSequencePuppi *
        process.PUPPIMuonRelIso * process.ttHTreeMaker)
#        process.PUPPIMuonRelIso * process.electronMVAValueMapProducer * process.ttHTreeMaker)
