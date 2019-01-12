from WMCore.Configuration import Configuration 
config = Configuration() 
config.section_("General") 
config.General.requestName = 'triggerSF_newJECJER_v0_PeriodD-31Mar2018-v1' ## change//Name of the directory that crab will create

config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'yggdrasil_treeMaker_ttH_sync_data.py' ###Config file that crab will use for the job
#config.JobType.psetName = 'yggdrasil_treeMaker_cfg_addLeptonSFs_ttbarMC.py' ###Config file that crab will use for the job

config.JobType.allowUndistributedCMSSW = False

config.section_("Data") ###Might need to chance to MC
config.Data.inputDataset = '/MET/Run2017D-31Mar2018-v1/MINIAOD' #MET, period D
#'/MET/Run2017E-31Mar2018-v1/MINIAOD' #MET, period E
#'/MET/Run2017F-31Mar2018-v1/MINIAOD' #MET, period F
#'/MET/Run2017C-31Mar2018-v1/MINIAOD' #MET, period C
#'/MET/Run2017B-31Mar2018-v1/MINIAOD' #MET, period B


#'/MuonEG/Run2017F-31Mar2018-v1/MINIAOD' # e+mu, period F
#'/MuonEG/Run2017E-31Mar2018-v1/MINIAOD' # e+mu, period E
#'/MuonEG/Run2017D-31Mar2018-v1/MINIAOD' # e+mu, period D
#'/MuonEG/Run2017C-31Mar2018-v1/MINIAOD' # e+mu, period C
#'/MuonEG/Run2017B-31Mar2018-v1/MINIAOD' # e+mu, period B

#'/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD' # DoubleEl, period F
#'/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD' # DoubleEl, period E
#'/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD' # DoubleEl, period D
#'/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD' # DoubleEl, period C
#'/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD' # DoubleEl, period B

#'/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD' # DoubleMu, period F
#'/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD' # DoubleMu, period E
#'/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD' # DoubleMu, period D
#'/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD' # DoubleMu, period C
#'/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD' # DoubleMu, period B


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
config.Data.outLFNDirBase = '/store/group/lpctthrun2/analysis2017Data/triggerSF/MET/'
#config.Data.outLFNDirBase = '/store/user/btannenw/'

#'/store/group/lpcmj/data/UVa/2018Data'###Wherever we can find space and start after '/eos/uscms' in the directory path.
#'/store/user/lwming/'



config.Data.ignoreLocality = False #Changed from True Feb 13

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC'
