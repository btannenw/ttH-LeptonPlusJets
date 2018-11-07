from WMCore.Configuration import Configuration 
config = Configuration() 
config.section_("General") 
config.General.requestName = 'ICHEP18_postMCsync_v0_WJetsToLNu' ## change//Name of the directory that crab will create

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'yggdrasil_treeMaker_ttH_sync.py' ###Config file that crab will use for the job
#config.JobType.psetName = 'yggdrasil_treeMaker_cfg_addLeptonSFs_ttbarMC.py' ###Config file that crab will use for the job

config.JobType.allowUndistributedCMSSW = False

config.section_("Data") ###Might need to chance to MC
config.Data.inputDataset = '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM' # W+jets (W->lnu)
#'/ttHTobb_ttTo2L2Nu_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttH, H->bb, tt DL
#'/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttH, H-> bb
#'/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM' # ttbar SL
#'/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' ## ttbar DL
#'/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttbar, fully hadronic

###Add the full Dataset you find from CMSDAS site here.
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/' ###Probably wont need to change
config.Data.splitting = 'Automatic' ###Can be changed to Automatic, FileBased, LumiBased
#config.Data.unitsPerJob = 10
#Commented OUT Feb13  
	#config.Data.publication = True 
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' ###More a factor if generating MC samples, usually just fine.
#Commented OUT Feb13
	#config.Data.publishDataName = 'Validation_v1'
### change user Space 
config.Data.outLFNDirBase = '/store/group/lpctthrun2/UVA/ICHEP2018/MC/'
#config.Data.outLFNDirBase = '/store/user/btannenw/'

#'/store/group/lpcmj/data/UVa/2018Data'###Wherever we can find space and start after '/eos/uscms' in the directory path.
#'/store/user/lwming/'



config.Data.ignoreLocality = False #Changed from True Feb 13

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC'
