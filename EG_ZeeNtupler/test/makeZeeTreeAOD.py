import FWCore.ParameterSet.Config as cms

process = cms.Process("ZeeTree")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10


inputFilesMC   = cms.untracked.vstring(
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00C0BECF-6F14-E511-96F8-0025904B739A.root',
    '/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0260F225-7614-E511-A79F-00A0D1EE8EB4.root',
    )

inputFilesData = cms.untracked.vstring(
    '/store/data/Run2015D/DoubleEG/AOD/16Dec2015-v2/00000/0033E989-FEA5-E511-A55E-0CC47A78A41C.root',
    '/store/data/Run2015D/DoubleEG/AOD/16Dec2015-v2/00000/003C08D3-BAA5-E511-92D9-002618943935.root',
    '/store/data/Run2015D/DoubleEG/AOD/16Dec2015-v2/00000/0050F783-D9A5-E511-98FC-0CC47A4C8F0A.root',
    '/store/data/Run2015D/DoubleEG/AOD/16Dec2015-v2/00000/00914683-DBA5-E511-8A50-0CC47A4D7694.root'
    )




isData = True
outputFile = 'unknown'

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50'
#    for Spring15 25ns MC: global tag is 'auto:run2_mc'
#    for Run 2 data: global tag is 'auto:run2_data'
#  as a rule, find the "auto" global tag in $CMSSW_RELEASE_BASE/src/Configuration/AlCa/python/autoCond.py
from Configuration.AlCa.GlobalTag import GlobalTag
if isData:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data' , '')
    inputFiles = inputFilesData
    outputFile = "ZeeTreeFab_Data.root"
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc'   , '')
    inputFiles = inputFilesMC
    outputFile = "ZeeTreeFab_MC.root"



#
# Define input data to read
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

#
# Set up VID framework
#

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
switchOnVIDElectronIdProducer( process, DataFormat.AOD)
switchOnVIDPhotonIdProducer(   process, DataFormat.AOD)

# define which IDs we want to produce
eleid_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']

#add them to the VID producer for electrons and photons
for idmod in eleid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

for idmod in phoid_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection  )

    

#
# Configure an example module for user analysis with electrons
#

process.ntupler = cms.EDAnalyzer(
    'EG_ZeeNtupler',
    # The module automatically detects AOD vs miniAOD, so we configure both
    #
    # Common to all formats objects
    #
    #    beamSpot = cms.InputTag('offlineBeamSpot'),
    electrons = cms.InputTag("gedGsfElectrons"),
    photons   = cms.InputTag("photons"),
    vertices  = cms.InputTag("offlinePrimaryVertices"),
    rho       = cms.InputTag("fixedGridRhoAll"),
    genParticles =  cms.InputTag("prunedGenParticles"),
#
    # ID decisions (common to all formats)
    #
    eleId = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    phoId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose")
#    phoId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
#    eleId = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
#    eleId = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),

    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

# Make sure to add the ID sequence upstream from the user analysis module
process.p = cms.Path( process.egmGsfElectronIDSequence* process.egmPhotonIDSequence * process.ntupler)
