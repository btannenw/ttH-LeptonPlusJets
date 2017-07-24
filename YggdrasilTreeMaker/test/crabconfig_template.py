
# information : https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

##########################################
# General : In this section, the user specifies generic parameters about the request (e.g. request name).

config.General.requestName  = 'XXXXX_YYYYY'
config.Data.outputDatasetTag= 'XXXXX_YYYYY'

config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

##########################################
# JobType : This section aims to contain all the parameters of the user job type and related configurables (e.g. CMSSW parameter-set configuration file, additional input files, etc.).

# 'Analusis' means that this job uses an existing file(s). If create private smaple, set it to'PrivateMC'.
config.JobType.pluginName = 'Analysis'

# :The name of the CMSSW parameter-set configuration file that should be run via cmsRun. 
config.JobType.psetName = '__yggdrasil_treeMaker_QQQQQ_cfg.py'

config.JobType.allowUndistributedCMSSW = True
config.JobType.maxJobRuntimeMin = 2000 # default 1315 = 21 hours.
# config.JobType.inputFiles = ['']

config.JobType.sendExternalFolder = True

##########################################
# Information of data to be analyzed, including the splitting parameters
                            
#config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/MINIAODSIM'
config.Data.inputDataset = 'ZZZZZ'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
# config.Data.unitsPerJob = 10  # <<<<<<<- normal.
config.Data.outLFNDirBase = '/store/group/lpctthrun2/'
config.Data.publication = True

#PPPPP#


config.Site.storageSite = 'T3_US_FNALLPC'

