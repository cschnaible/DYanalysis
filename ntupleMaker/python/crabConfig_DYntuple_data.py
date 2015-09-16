from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'singleMuon_data_50ns'
config.General.workArea = 'crab'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/DYntupleMaker_cfg.py'

config.Data.inputDataset = '/SingleMuon/Run2015B-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/DoubleMuon/Run2015B-PromptReco-v1/MINIAOD'
config.Data.inputDBS = 'global'

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON.txt'


config.Data.outLFNDirBase = '/store/user/cschnaib/'
config.Data.publication = False

config.Site.storageSite = 'T2_CH_CERN'
