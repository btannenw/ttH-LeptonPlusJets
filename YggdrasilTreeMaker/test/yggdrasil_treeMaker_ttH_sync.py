import FWCore.ParameterSet.Config as cms


isPUPPI=False
#--> if isPUPPI=true, change the code too.


#isMC=True
isMC=True

# 
isTTBARMC=False

# isGridJob=False
isGridJob=False

genjetInputTag = cms.InputTag("slimmedGenJets","","")



enableJECFromLocalDB = True 
# 2017_V6 JECs for Data/MC are ready.

# - - - - - - - - - - - - - - - - - - - - 
# Special option for Morind17 analysis
#  This is a flag used to apply dedicated JEC for each data set.
#  The placeholder will be replaced by crab job make script.

period="XXXPERIODXXX"
# e.g "2017B"

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
    process.GlobalTag.globaltag = '94X_mc2017_realistic_v12'
else :
    process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
    #input = cms.untracked.int32( -1 )
    )


import sys
import os.path

JecLocalDataBaseName = \
    'Fall17_17Nov2017B_V6_DATA' if period in ("2017B") else \
    'Fall17_17Nov2017C_V6_DATA' if period in ("2017C") else \
    'Fall17_17Nov2017D_V6_DATA' if period in ("2017D") else \
    'Fall17_17Nov2017E_V6_DATA' if period in ("2017E") else \
    'Fall17_17Nov2017F_V6_DATA' if period in ("2017F") else \
    'Fall17_17Nov2017_V6_MC'


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
            tag    = cms.string('JetCorrectorParametersCollection_'+ JecLocalDataBaseName +'_AK8PFPuppi'),
            label  = cms.untracked.string('AK8PFPuppi')
            )
        )

    process.GlobalTag.toGet.append(
        cms.PSet(
            connect = cms.string( JecDBPathPrefix + '/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/' + JecLocalDataBaseName +'.db' ),
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_'+ JecLocalDataBaseName +'_AK8PFPuppi'),
            label  = cms.untracked.string('AK15PFPuppi')
            )
        )

    #  line for PUPPI for temporal.
    process.GlobalTag.toGet.append(
        cms.PSet(
            connect = cms.string( JecDBPathPrefix +'/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/data/' + JecLocalDataBaseName +'.db' ),
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_'+ JecLocalDataBaseName +'_AK4PFPuppi'),
            label  = cms.untracked.string('AK4PFPuppi')
            )
        )


process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
# where 
#  -   ak4PFCHSL1FastL2L3CorrectorChain 
#  - ak4PFPuppiL1FastL2L3CorrectorChain
#  -   ak4PFCHSL1FastL2L3ResidualCorrectorChain
#  - ak4PFPuppiL1FastL2L3ResidualCorrectorChain
#  are defined.
# https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/JetMETCorrections/Configuration/python/JetCorrectors_cff.py


process.ak8PFPuppiL1FastjetCorrector = cms.EDProducer(
    'L1FastjetCorrectorProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK8PFPuppi'),
    srcRho      = cms.InputTag( 'fixedGridRhoFastjetAll' )
    )
process.ak8PFPuppiL2RelativeCorrector = process.ak4CaloL2RelativeCorrector.clone( algorithm = 'AK8PFPuppi' )
process.ak8PFPuppiL3AbsoluteCorrector = process.ak4CaloL3AbsoluteCorrector.clone( algorithm = 'AK8PFPuppi' )
process.ak8PFPuppiL1FastL2L3Corrector = process.ak4PFPuppiL2L3Corrector.clone()
process.ak8PFPuppiL1FastL2L3Corrector.correctors.insert(0,'ak8PFPuppiL1FastjetCorrector')

process.ak8PFPuppiL1FastL2L3CorrectorChain = cms.Sequence(
    process.ak8PFPuppiL1FastjetCorrector * process.ak8PFPuppiL2RelativeCorrector * process.ak8PFPuppiL3AbsoluteCorrector * process.ak8PFPuppiL1FastL2L3Corrector
)


process.ak8PFPuppiResidualCorrector  = process.ak4CaloResidualCorrector.clone( algorithm = 'AK8PFPuppi' )
process.ak8PFPuppiL1FastL2L3ResidualCorrector = cms.EDProducer(
    'ChainedJetCorrectorProducer',
    correctors = cms.VInputTag('ak8PFPuppiL1FastjetCorrector','ak8PFPuppiL2RelativeCorrector','ak8PFPuppiL3AbsoluteCorrector','ak8PFPuppiResidualCorrector')
    )
process.ak8PFPuppiL1FastL2L3ResidualCorrectorTask = cms.Task(
    process.ak8PFPuppiL1FastjetCorrector, process.ak8PFPuppiL2RelativeCorrector, process.ak8PFPuppiL3AbsoluteCorrector, process.ak8PFPuppiResidualCorrector, process.ak8PFPuppiL1FastL2L3ResidualCorrector
)
process.ak8PFPuppiL1FastL2L3ResidualCorrectorChain = cms.Sequence( process.ak8PFPuppiL1FastL2L3ResidualCorrectorTask)




process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(

# xrdcp root://cms-xrd-global.cern.ch//store/data/Run2017E/SingleMuon/MINIAOD/17Nov2017-v1/50000/000DCB8B-2ADD-E711-9100-008CFAF35AC0.root ~/temporal_strage/Run2017E_SingleMuon_MINIAOD.root 
# 
# xrdcp root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAOD/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/0044D634-7FED-E711-B0EF-0242AC130002.root ~/temporal_strage/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_MINIAODSIM.root

# 2017-data analysis
#'file:///uscmst1b_scratch/lpc1/3DayLifetime/satoshi/Run2017E_SingleMuon_MINIAOD.root'
# 'file:///uscmst1b_scratch/lpc1/3DayLifetime/satoshi/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_MINIAODSIM.root'
            
# 2018 ttH sync exercise 
  #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAOD/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/CC0FCC49-B50A-E811-9694-02163E0144C8.root' # v10
  #'root://cmsxrootd.fnal.gov//store/mc/RunIIFall17MiniAODv2/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/60000/AC628CE7-0169-E811-9C5E-00010100096B.root' # v14
  'file:///afs/cern.ch/work/b/btannenw/ttH/yggdrasil2018/triggerEff/localUpdates/CMSSW_9_4_9/src/ttH-LeptonPlusJets/YggdrasilTreeMaker/test/inputFiles/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8__MINIAODSIM__PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1__00096B.root' # v-14 local
        )
)



# collection placeholders
#triggerCollection  = cms.InputTag("triggerResults", "", "RECO")
electronCollection = cms.InputTag("slimmedElectrons", "", "PAT")
tauCollection      = cms.InputTag("slimmedTaus", "", "PAT")
photonCollection   = cms.InputTag("slimmedPhotons", "", "PAT")
METCollection      = cms.InputTag("slimmedMETs", "", "PAT")
jetCollection      = cms.InputTag("slimmedJets", "", "PAT")

muonCollection     = cms.InputTag("slimmedMuons", "", "PAT")
if not isMC : 
    muonCollection     = cms.InputTag("slimmedMuons", "", "RECO")


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

    process.ak8ReclusteredGenJets = ak4GenJets.clone(
        src = 'genParticlesForJetsNoNu',
        rParam = cms.double(0.8),
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



from PhysicsTools.SelectorUtils.tools.vid_id_tools import DataFormat, \
        switchOnVIDElectronIdProducer, setupAllVIDIdsInModule, setupVIDElectronSelection

eleVIDModules = [
    # 10-08-18, BBT: comment out old IDs and add 94x
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff",
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff",
    "RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff"
    #"RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff"
]

switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

for mod in eleVIDModules:
    setupAllVIDIdsInModule(process, mod, setupVIDElectronSelection)

# update some VID modules to work with potentially changed electron collections
process.egmGsfElectronIDs.physicsObjectSrc = electronCollection
process.electronRegressionValueMapProducer.srcMiniAOD = electronCollection
process.electronMVAValueMapProducer.srcMiniAOD = electronCollection



### EGamma routine to apply IDs + Scale & Smearing Corrections
### REF: https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#How_to_run_the_Scale_Smearing_co
# BBT, 10-10-18
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

#setupEgammaPostRecoSeq(process, isMiniAOD=True,
#  applyEnergyCorrections    = False,
#  applyVIDOnCorrectedEgamma = False,
#)

# BBT, 10-10-18
#from RecoEgamma.EgammaTools.calibratedEgammas_cff import calibratedPatElectrons as calibratedElectrons
#proc.correctedElectrons = calibratedElectrons.clone(src = collection)
#proc.correctedElectrons.produceCalibratedObjs = cms.bool(True)
#proc.correctedElectrons.semiDeterministic = cms.bool(True)


# BBT, 10-11-18
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=False, #saves CPU time by not needlessly re-running VID
                       #era='Run2017_17Nov2017_v1')  
                       era='2017-Nov17ReReco')  
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)






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
	fileName = cms.string('yggdrasil_treeMaker_ttH_sync_11-03-18_v22_addSFs.root')
)


process.load( "PuppiLeptonIsolationhelper.PuppiLeptonIsolation.PuppiMuonIsolationProducer_cfi" )
process.PUPPIMuonRelIso.addIsolationComponentInfo = cms.bool( True )
if not isMC :
    process.PUPPIMuonRelIso.pfCandidates  = cms.InputTag("packedPFCandidates", "", "RECO")
    process.PUPPIMuonRelIso.muonCollection  = muonCollection
    

#BBT, 10-12-18, add deterministic seeds
process.load("PhysicsTools.PatUtils.deterministicSeeds_cfi")
process.deterministicSeeds.produceCollections = cms.bool(True)
process.deterministicSeeds.produceValueMaps   = cms.bool(False)
#process.deterministicSeeds.seedUserInt        = cms.string('deterministicSeed')
process.deterministicSeeds.electronCollection = electronCollection
process.deterministicSeeds.muonCollection     = muonCollection
process.deterministicSeeds.tauCollection      = tauCollection
process.deterministicSeeds.photonCollection   = photonCollection
process.deterministicSeeds.jetCollection      = jetCollection
process.deterministicSeeds.METCollection      = METCollection

#BBT, 10-12-18, add deterministic seeds
# overwrite output collections (skip electrons as the EGamma tools handle the deterministic seed creation internally)
#print "@@@@@@@@@@@ ", process.name_()

Muoncollection     = cms.InputTag("deterministicSeeds", "muonsWithSeed", process.name_())
tauCollection      = cms.InputTag("deterministicSeeds", "tausWithSeed", process.name_())
photonCollection   = cms.InputTag("deterministicSeeds", "photonsWithSeed", process.name_())
jetCollection      = cms.InputTag("deterministicSeeds", "jetsWithSeed", process.name_())
METCollection      = cms.InputTag("deterministicSeeds", "METsWithSeed", process.name_())


## The line below should always be included, BBT, 10-09-18
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
## Example 1: If you only want to re-correct MET and get the proper uncertainties [e.g. when updating JEC]
runMetCorAndUncFromMiniAOD(process,
                           isData           = not isMC
#            jecUncFile       = os.path.basename(options.JESUncFiles[0]),
#            electronColl     = electronCollection.value(),
#            muonColl         = muonCollection.value(),
#            tauColl          = tauCollection.value(),
#            photonColl       = photonCollection.value(),
#            jetCollUnskimmed = jetCollection.value(),
#            #pfCandColl=cms.InputTag("packedPFCandidates"),
#            #recoMetFromPFCs  = True
                           )


#process.load("Configuration.StandardSequences.MagneticField_cff")
#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox

process.load('CommonTools/PileupAlgos/Puppi_cff')
## e.g. to run on miniAOD                                                                                                                                                     
process.puppi.candName = cms.InputTag('packedPFCandidates')
process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')

process.puppiOnTheFly = process.puppi.clone()
#process.puppiOnTheFly.useExistingWeights = False
process.puppiOnTheFly.useExistingWeights = True
# --- To be consistent with MiniAOD jets, keep using the same PUPPI weight (Is this correct approach ?)
#jetToolbox( process, 'ca15', 'jetSequence', 'out', 
#            runOnMC= isMC , addNsub=True, 
#            PUMethod = 'Puppi', newPFCollection=True, nameNewPFCollection='puppiOnTheFly' , addSoftDrop=True , addSoftDropSubjets=True )



if isMC : 
    process.p = cms.Path(
        process.genParticlesForJetsNoNu * process.ak4GenJetsCustom *
        process.selectedHadronsAndPartons * 
        process.genJetFlavourInfos *
        process.matchGenCHadron *
        process.matchGenBHadron * 
        process.categorizeGenTtbar *
        process.ak4PFCHSL1FastL2L3CorrectorChain *
        process.ak4PFPuppiL1FastL2L3CorrectorChain *
        process.ak8PFPuppiL1FastL2L3CorrectorChain *
        process.GenParticleWithoutChargedLeptonFropTop * process.myGenParticlesWithChargedLeptonFromTopForJet * process.ak4GenJetsWithChargedLepFromTop *  
        process.ak8ReclusteredGenJets *
        process.PUPPIMuonRelIso * 
        process.fullPatMetSequence *  # BBT, 10-04-18
        process.egmGsfElectronIDSequence * # BBT, 10-08-18
        process.egammaPostRecoSeq * # BBT, 10-11-18
        process.deterministicSeeds * # BBT, 10-12-18
        process.ttHTreeMaker)
else :
    process.p = cms.Path(
        process.ak4PFCHSL1FastL2L3ResidualCorrectorChain *
        process.ak4PFPuppiL1FastL2L3ResidualCorrectorChain *
        process.ak8PFPuppiL1FastL2L3ResidualCorrectorChain *
        process.PUPPIMuonRelIso * process.ttHTreeMaker)
