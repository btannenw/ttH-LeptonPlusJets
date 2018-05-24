import FWCore.ParameterSet.Config as cms


isPUPPI=False
#--> if isPUPPI=true, change the code too.


#isMC=True
isMC=True

# 
isTTBARMC=False

# isGridJob=False
isGridJob=True

genjetInputTag = cms.InputTag("slimmedGenJets","","")
#genjetInputTag = cms.InputTag("ak4GenJetsReproduced","","")
#genjetInputTag = cms.InputTag("ak4GenJetsWithChargedLepFromTop","","")


enableJECFromLocalDB=True

# - - - - - - - - - - - - - - - - - - - - 
# Special option for Morind17 analysis
#  This is a flag used to apply dedicated JEC for each data set.
#  The placeholder will be replaced by crab job make script.

period="XXXPERIODXXX"
# "2016B" , "2016C", "2016D"
# "2016E" , "2016F1"
# "2016F2", "2016G" 
# "2016H1", "2016H2"


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


process.load("Configuration.StandardSequences.GeometryDB_cff")

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
    input = cms.untracked.int32(-1)
    )




import sys
import os.path

JecLocalDataBaseName = \
    'Summer16_23Sep2016BCDV4_DATA' if period in ("2016B" , "2016C", "2016D") else \
    'Summer16_23Sep2016EFV4_DATA'  if period in ("2016E" , "2016F1") else \
    'Summer16_23Sep2016GV4_DATA'   if period in ("2016F2", "2016G" ) else \
    'Summer16_23Sep2016HV4_DATA'   if period in ("2016H1", "2016H2") else 'Summer16_23Sep2016V4_MC'

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
            connect = cms.string( JecDBPathPrefix +'/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/' + JecLocalDataBaseName +'.db' ),
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_'+JecLocalDataBaseName+'_AK4PFPuppi'),
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
# ttH
#         'file:///uscms/home/satoshi/temporal_strage/44949CF4-96C6-E611-B9A0-0025905A6122.root'

# 
#  xrdcp root://cmsxrootd.fnal.gov//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/60000/CEAB3688-1CC7-E611-8BC3-C4346BBCB6A8.root ~/temporal_strage/DY_80x_MINIAODSIM.root 
 'file:///uscms/home/satoshi/temporal_strage/DY_80x_MINIAODSIM.root'
            )
)



# collection placeholders
#
electronCollection = cms.InputTag("slimmedElectrons", "", "PAT")
muonCollection     = cms.InputTag("slimmedMuons", "", "PAT")
tauCollection      = cms.InputTag("slimmedTaus", "", "PAT")
photonCollection   = cms.InputTag("slimmedPhotons", "", "PAT")
METCollection      = cms.InputTag("slimmedMETs", "", "PAT")

if isPUPPI :
    jetCollection      = cms.InputTag("slimmedJetsPuppi", "", "PAT")
else : 
    jetCollection      = cms.InputTag("slimmedJets", "", "PAT")


# deterministic seed producer
#
process.load("PhysicsTools.PatUtils.deterministicSeeds_cfi")
process.deterministicSeeds.produceCollections = cms.bool(True)
process.deterministicSeeds.produceValueMaps   = cms.bool(False)
process.deterministicSeeds.electronCollection = electronCollection
process.deterministicSeeds.muonCollection     = muonCollection
process.deterministicSeeds.tauCollection      = tauCollection
process.deterministicSeeds.photonCollection   = photonCollection
process.deterministicSeeds.jetCollection      = jetCollection
process.deterministicSeeds.METCollection      = METCollection
# overwrite output collections
electronCollection = cms.InputTag("deterministicSeeds", "electronsWithSeed", process.name_())
muonCollection     = cms.InputTag("deterministicSeeds", "muonsWithSeed", process.name_())
tauCollection      = cms.InputTag("deterministicSeeds", "tausWithSeed", process.name_())
photonCollection   = cms.InputTag("deterministicSeeds", "photonsWithSeed", process.name_())
jetCollection      = cms.InputTag("deterministicSeeds", "jetsWithSeed", process.name_())
METCollection      = cms.InputTag("deterministicSeeds", "METsWithSeed", process.name_())





###################################


from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)
process.load("EgammaAnalysis.ElectronTools.regressionApplication_cff")
# User add "process.regressionApplication" to the sequence.

# set the electron and photon sources
process.slimmedElectrons.src = electronCollection
process.slimmedPhotons.src = photonCollection
# overwrite output collections
electronCollection = cms.InputTag("slimmedElectrons", "", process.name_())
photonCollection = cms.InputTag("slimmedPhotons", "", process.name_())


###################################
# Egamma energy smearing for data and mc

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
    src = electronCollection,
    cut = cms.string("pt>5 && abs(superCluster.eta)<2.5")
)
electronCollection = cms.InputTag("selectedElectrons", "", process.name_())
# setup the smearing
process.load("EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi")
from EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi import files
process.calibratedPatElectrons.isMC           = cms.bool( isMC )
process.calibratedPatElectrons.correctionFile = cms.string(files[ "Moriond17_23Jan" ])
process.calibratedPatElectrons.electrons      = electronCollection
#  seq += process.calibratedPatElectrons -> Satoshi added to the seqyence by hand. See the bottom part of this script.

# use our deterministic seeds :
process.calibratedPatElectrons.seedUserInt = process.deterministicSeeds.seedUserInt
# overwrite output collections
electronCollection = cms.InputTag("calibratedPatElectrons", "", process.name_())
  

###################################


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



#
# electron VIDs
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat, \
        switchOnVIDElectronIdProducer, setupAllVIDIdsInModule, setupVIDElectronSelection

eleVIDModules = [
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff",
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff",
    "RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff"
]

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)




for mod in eleVIDModules:
    setupAllVIDIdsInModule(process, mod, setupVIDElectronSelection)

# update some VID modules to work with potentially changed electron collections
process.egmGsfElectronIDs.physicsObjectSrc = electronCollection
process.electronRegressionValueMapProducer.srcMiniAOD = electronCollection
process.electronMVAValueMapProducer.srcMiniAOD = electronCollection







from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
# do not use a postfix here!
runMetCorAndUncFromMiniAOD(process,
    isData           = not isMC,
    electronColl     = electronCollection.value(),
    muonColl         = muonCollection.value(),
    tauColl          = tauCollection.value(),
    photonColl       = photonCollection.value(),
    jetCollUnskimmed = jetCollection.value(),
    recoMetFromPFCs  = True
)
METCollection = cms.InputTag("slimmedMETs", "", process.name_())

# also add MET corrections due to e/g corrections, such as the slew rate fix in reMiniAOD
if not isMC :
    from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
    corMETFromMuonAndEG(process,
        pfCandCollection      = "",
        electronCollection    = "slimmedElectronsBeforeGSFix",
        photonCollection      = "slimmedPhotonsBeforeGSFix",
        corElectronCollection = electronCollection.value(),
        corPhotonCollection   = photonCollection.value(),
        allMETEGCorrected     = True,
        muCorrection          = False,
        eGCorrection          = True,
        runOnMiniAOD          = True,
        postfix               = "MuEGClean"
    )
    process.slimmedMETsMuEGClean = getattr(process, METCollection.getModuleLabel()).clone(
        src             = cms.InputTag("patPFMetT1MuEGClean"),
        rawVariation    = cms.InputTag("patPFMetRawMuEGClean"),
        t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
    )
    del process.slimmedMETsMuEGClean.caloMET

    # overwrite output collections
    METCollection = cms.InputTag("slimmedMETsMuEGClean", "", process.name_())


# patch the phi correction parameter sets that are used in runMetCorAndUncFromMiniAOD,
# we only need to overwrite patMultPhiCorrParams_T1Txy_25ns with the new one
if not isMC :
    if period in ("2016B", "2016C", "2016D", "2016E", "2016F1" , "2016F2"):
        from MetTools.MetPhiCorrections.tools.multPhiCorr_ReMiniAOD_Data_BCDEF_80X_sumPt_cfi \
                import multPhiCorr_Data_BCDEF_80X as metPhiCorrParams
    else: 
        from MetTools.MetPhiCorrections.tools.multPhiCorr_ReMiniAOD_Data_GH_80X_sumPt_cfi \
                import multPhiCorr_Data_GH_80X as metPhiCorrParams
else:
    from MetTools.MetPhiCorrections.tools.multPhiCorr_Summer16_MC_DY_80X_sumPt_cfi \
            import multPhiCorr_MC_DY_sumPT_80X as metPhiCorrParams
# actual patch
getattr(process, "patPFMetTxyCorr").parameters = cms.VPSet(pset for pset in metPhiCorrParams)



#
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#runMetCorAndUncFromMiniAOD(process,
#                           isData= not isMC , 
##(relative path from src/)  jecUncFile = ( 'filepath' ),
#                           )
#
#process.load("Configuration.StandardSequences.GeometryDB_cff")
#
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


process.load("RecoMET.METFilters.BadPFMuonFilter_cfi")
process.BadPFMuonFilter.muons        = muonCollection
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load("RecoMET.METFilters.BadChargedCandidateFilter_cfi")
process.BadChargedCandidateFilter.muons        = muonCollection
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load("RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff")
process.badGlobalMuonTaggerMAOD.muons         = muonCollection
process.badGlobalMuonTaggerMAOD.taggingMode   = cms.bool(True)
process.cloneGlobalMuonTaggerMAOD.muons       = muonCollection
process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)



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
                                         , configuration = cms.string( "#detail#" )
                                         )



# process.load('RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi')
# process.load('RecoEgamma.PhotonIdentification.PhotonMVAValueMapProducer_cfi')
# process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
# process.photonMVAValueMapProducer.srcMiniAOD   = cms.InputTag('slimmedPhotons')


if isMC : 
    print "[debuyg message ]  process process.calibratedPatElectrons is added."
    process.p = cms.Path(
        process.regressionApplication *   #<- Electron energy regression 
        process.calibratedPatElectrons * 
#        process.egmPhotonIDSequence * process.puppiMETSequence * process.fullPatMetSequencePuppi *
        process.GenParticleWithoutChargedLeptonFropTop * process.myGenParticlesWithChargedLeptonFromTopForJet * process.ak4GenJetsWithChargedLepFromTop *  
        process.PUPPIMuonRelIso * process.ttHTreeMaker)
#        process.PUPPIMuonRelIso * process.electronMVAValueMapProducer * process.ttHTreeMaker)
else :
    process.p = cms.Path(
        process.regressionApplication *   #<- Electron energy regression 
        process.calibratedPatElectrons * 
#        process.egmPhotonIDSequence * process.puppiMETSequence * process.fullPatMetSequencePuppi *
        process.PUPPIMuonRelIso * process.ttHTreeMaker)
#        process.PUPPIMuonRelIso * process.electronMVAValueMapProducer * process.ttHTreeMaker)
