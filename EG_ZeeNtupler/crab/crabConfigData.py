from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName  = 'makeZeeTreeData.py'

config.Data.inputDBS = 'global'
config.Data.splitting    = 'LumiBased'
config.Data.lumiMask     = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt'
config.Data.unitsPerJob = 100
config.Data.publication = False
config.Data.outLFNDirBase   = '/store/user/%s/EGamma/ZeeNtupless/miniAOD-16Dec2015/DataCorrectRho_bugFix/' % (getUsernameFromSiteDB())
#config.Data.publishDataName = 

config.Site.storageSite = 'T2_FR_GRIF_IRFU'



if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects_DATA_v1_CorrectedRho_bugFix'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    

    config.General.requestName  = 'RunD_25ns'
    config.Data.inputDataset    = '/DoubleEG/Run2015D-16Dec2015-v2/MINIAOD'
    config.Data.allowNonValidInputDataset = True
    submit(config)

    config.General.requestName  = 'RunC_25ns'
    config.Data.inputDataset    = '/DoubleEG/Run2015C_25ns-16Dec2015-v1/MINIAOD'
    config.Data.allowNonValidInputDataset = False
    submit(config)


