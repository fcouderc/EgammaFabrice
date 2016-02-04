from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.transferLogs = False

config.JobType.pluginName  = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName  = 'makeZeeTreeMC.py'

config.Data.inputDBS     = 'global'
config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
config.Data.outLFNDirBase   = '/store/user/%s/EGamma/ZeeNtupless/miniAOD-16Dec2015/MC/' % (getUsernameFromSiteDB())
#config.Data.publishDataName = 

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_FR_GRIF_IRFU'



if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'crab_projects_MC_v1'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)


    

    config.General.requestName  = 'DY_powheg'
    config.Data.inputDataset    = '/DYToEE_NNPDF30_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM'
    submit(config)

#    config.General.requestName  = 'RunC_25ns'
#    config.Data.inputDataset    = '/DoubleEG/Run2015C_25ns-16Dec2015-v1/MINIAOD'
#    config.Data.allowNonValidInputDataset = False
#    submit(config)


