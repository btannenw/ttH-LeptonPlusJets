from WMCore.Configuration import Configuration 
config = Configuration() 
config.section_("General") 
#config.General.requestName = 'preApproval2017_v2_TTTo2L2Nu' ## change//Name of the directory that crab will create
#config.General.requestName = 'preApproval2017_v2_TTToHadronic' 
#config.General.requestName = 'preApproval2017_v2_TTToSemiLeptonic'
#config.General.requestName = 'preApproval2017_v2_ttHtobb' 
#config.General.requestName = 'preApproval2017_v2_ttHTobb_ttTo2L2Nu' 
#config.General.requestName = 'preApproval2017_v2_TTGJets' 
#config.General.requestName = 'preApproval2017_v2_TTZToQQ' 
#config.General.requestName = 'preApproval2017_v2_TTZToLLNuNu' 
#config.General.requestName = 'preApproval2017_v2_TTWJetstoQQ' 
#config.General.requestName = 'preApproval2017_v2_TTWJetsToLNu' 
#config.General.requestName = 'preApproval2017_v2_WW' 
#config.General.requestName = 'preApproval2017_v2_ZZ' 
#config.General.requestName = 'preApproval2017_v2_WZ' 
config.General.requestName = 'preApproval2017_v2_WJetsToLNu' 
#config.General.requestName = 'preApproval2017_v2_ST_s-channel' 
#config.General.requestName = 'preApproval2017_v2_ST_t-channel_top' 
#config.General.requestName = 'preApproval2017_v2_ST_t-channel_antitop' 
#config.General.requestName = 'preApproval2017_v2_ST_tW_top' 
#config.General.requestName = 'preApproval2017_v2_ST_tW_antitop'
#config.General.requestName = 'preApproval2017_v2_DYJetsToLL_M-10to50'
#config.General.requestName = 'preApproval2017_v2_DYJetsToLL_M-50'

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'yggdrasil_treeMaker_ttH_sync_combined.py' ###Config file that crab will use for the job
#config.JobType.psetName = 'yggdrasil_treeMaker_cfg_addLeptonSFs_ttbarMC.py' ###Config file that crab will use for the job

config.JobType.allowUndistributedCMSSW = False

config.section_("Data") ###Might need to chance to MC
#config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' ## ttbar DL
#config.Data.inputDataset = '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttbar, fully hadronic
#config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM' # ttbar SL
#config.Data.inputDataset = '/ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttH, H-> bb
#config.Data.inputDataset = '/ttHTobb_ttTo2L2Nu_M125_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttH, H->bb, tt DL
#config.Data.inputDataset = '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # tt+gjets
#config.Data.inputDataset = '/TTZToQQ_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttZ, Z->qq
#config.Data.inputDataset = '/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttZ, Z->ll/nunu
#config.Data.inputDataset = '/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttW, W->qq
#config.Data.inputDataset = '/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ttW, W->lnu
#config.Data.inputDataset = '/WW_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # WW+jets
#config.Data.inputDataset = '/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # ZZ+jets
#config.Data.inputDataset = '/WZ_TuneCP5_13TeV-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # WZ+jets
config.Data.inputDataset = '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM' # W+jets (W->lnu)
#config.Data.inputDataset = '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM' # single top s-channel, W->lnu
#config.Data.inputDataset = '/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM' # single top t-channel (top)
#config.Data.inputDataset = '/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM' # single top t-channel, (anti-top)
#config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM' # single top tW, (top)
#config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM' # single top tW, (anti-top)
config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # Z+jets, mll10to50
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' # Z+jets, mll50toInf




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
config.Data.outLFNDirBase = '/store/group/lpctthrun2/preApproval2017/MC/'
#config.Data.outLFNDirBase = '/store/user/btannenw/'

#'/store/group/lpcmj/data/UVa/2018Data'###Wherever we can find space and start after '/eos/uscms' in the directory path.
#'/store/user/lwming/'



config.Data.ignoreLocality = False #Changed from True Feb 13

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC'
