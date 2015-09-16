import FWCore.ParameterSet.Config as cms

process = cms.Process("DYntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.MessageLogger.cerr.threshold = 'INFO'
#process.MessageLogger.categories.append('DYntuple')
#process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#        limit = cms.untracked.int32(-1)
#        )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
#)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00C0BECF-6F14-E511-96F8-0025904B739A.root',
'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0260F225-7614-E511-A79F-00A0D1EE8EB4.root',
'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02B810EA-7214-E511-BDAB-0025905964C2.root',
#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/02CEA7DD-7714-E511-A93E-00266CFAEA68.root',
#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0453351C-7014-E511-A296-0025905B85AA.root',
#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0679BC6F-7714-E511-945E-0025905B8562.root',
#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0823BF6F-7814-E511-8E48-00A0D1EE8B08.root',
#'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/08271551-9714-E511-B209-0025907FD2DA.root'
    )
)

process.accepted = cms.EDAnalyzer("DYntupleMaker",
    #srcRho          = cms.InputTag('kt6PFJets','rho'),
    #srcRho25        = cms.InputTag('kt6PFJetsCentralNeutral','rho'),
    #srcRhoQG        = cms.InputTag('kt6PFJetsIsoQG','rho'),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    genJets = cms.InputTag("slimmedGenJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    minLLMass = cms.double(10),
    minLepPt = cms.double(10),
    maxLepEta = cms.double(2.4),
    minLeadMuoPt = cms.double(20),
    minNjets        = cms.int32(1),
    #jetLepIsoRadius = cms.double(0.3),
    minJetPt        = cms.double(25),
    maxJetEta       = cms.double(2.5)
)

#process.TFileService = cms.Service("TFileService", fileName = cms.string('/afs/cern.ch/work/c/cschnaib/DYanalysis/outputRoot/DYJetsToLL_events_crabTest.root'))
process.TFileService = cms.Service("TFileService", fileName = cms.string('/afs/cern.ch/user/c/cschnaib/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/DYJetsToLL_events.root'))
#process.TFileService = cms.Service("TFileService", fileName = cms.string('/afs/cern.ch/user/c/cschnaib/CMSSW_7_4_6_patch6/src/DYanalysis/outputRoot/DYJetsToLL_Test.root'))

process.p = cms.Path(process.accepted)

