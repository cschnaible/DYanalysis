import FWCore.ParameterSet.Config as cms

process = cms.Process("DYntuple")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root',
'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/00C0BECF-6F14-E511-96F8-0025904B739A.root',
'/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/0260F225-7614-E511-A79F-00A0D1EE8EB4.root')
)


from Configuration.StandardSequences.GeometryRecoDB_cff import *
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# for PF weighted candidates (muon isolation)
process.load("CommonTools.ParticleFlow.deltaBetaWeights_cff")

#for data in 720pre7
process.GlobalTag.globaltag ='MCRUN2_74_V9'

process.analyzer = cms.EDAnalyzer("DYntupleMaker",
    is25ns          = cms.bool(False), # true for 25ns / false for 50ns
    srcRho          = cms.InputTag('fixedGridRhoFastjetAll'),
    #srcRho25        = cms.InputTag('kt6PFJetsCentralNeutral','rho'),
    #srcRhoQG        = cms.InputTag('kt6PFJetsIsoQG','rho'),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    gsfelectrons = cms.InputTag("reducedEGamma"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    genJets = cms.InputTag("slimmedGenJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversions = cms.InputTag("reducedEgamma:reducedConversions"),
    #prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),
    packedGenParticles = cms.InputTag("packedGenParticles"),
    prunedGenParticles = cms.InputTag("prunedGenParticles"),
    minLLMass = cms.double(10),
    minLepPt = cms.double(10),
    maxLepEta = cms.double(2.4),
    minLeadMuoPt = cms.double(20),
    maxCombRelIso03 = cms.double(0.15),
    maxCombRelIso04 = cms.double(0.12),
    minNjets        = cms.int32(1),
    #jetLepIsoRadius = cms.double(0.3),
    minJetPt        = cms.double(25),
    maxJetEta       = cms.double(2.5),
    processName     = cms.string('HLT'), 
    OnlyMC          = cms.bool(False),
    triggerName = cms.vstring("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
                              "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
                              "HLT_Mu27_TkMu8_v",
                              "HLT_Mu30_TkMu11_v",
                              "HLT_Mu40_TkMu11_v",
                              "HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
                              "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
                              "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
                              "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
                              "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",
                              "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v",
                              "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v",
                              "HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v",
                              "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v",
                              "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
                              "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
                              "HLT_IsoMu17_eta2p1",
                              "HLT_IsoMu20",
                              "HLT_IsoMu20_eta2p1",
                              "HLT_IsoMu24_eta2p1",
                              "HLT_IsoTkMu20",
                              "HLT_IsoTkMu20_eta2p1",
                              "HLT_IsoTkMu24_eta2p1"),
    triggerMuMu = cms.vstring("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
                              "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",
                              "HLT_Mu27_TkMu8_v",
                              "HLT_Mu30_TkMu11_v",
                              "HLT_Mu40_TkMu11_v"),
    triggerMuEle = cms.vstring("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
                               "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",
                               "HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
                               "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",
                               "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",
                               "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v",
                               "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v"),
    triggerEleEle = cms.vstring("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v",
                                "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v",
                                "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",
                                "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
                                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
    triggerMu = cms.vstring("HLT_IsoMu17_eta2p1",
                            "HLT_IsoMu20",
                            "HLT_IsoMu20_eta2p1",
                            "HLT_IsoMu24_eta2p1",
                            "HLT_IsoTkMu20",
                            "HLT_IsoTkMu20_eta2p1",
                            "HLT_IsoTkMu24_eta2p1")
)
'''
# Trigger Matching
pathTriggerMuMu = '(path("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",1,1) || path("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",1,1) || path("HLT_Mu27_TkMu8_v*",1,1) || path("HLT_Mu30_TkMu11_v*",1,1) || path("HLT_Mu40_TkMu11_v*",1,1))'

pathTriggerMuEle = '(path("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",1,1) || path("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",1,1) || path("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v*",1,1) || path("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",1,1) || path("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v*",1,1) || path("HLT_Mu8_DiEle12_CaloIdL_TrkIdL_v*",1,1) || path("HLT_DiMu9_Ele9_CaloIdL_TrkIdL_v*",1,1))'

pathTriggerEleEle = '(path("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v",1,1) || path("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v",1,1) || path("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",1,1) || path("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",1,1) || path("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",1,1))'

process.MuonTriggerMatchHLT = cms.EDProducer(
    "PATTriggerMatcherDRLessByR",
    src  = cms.InputTag('slimmedMuons'),
    matched = cms.InputTag('selectedPatTrigger'),
    matchedCuts = cms.InputTag(pathTriggerMuMu),
    maxDPtRel = cms.InputTag(0.5),
    maxDeltaR = cms.InputTag(0.2),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)
process.ElectronTriggerMatchHLT = cms.EDProducer(
    "PATTriggerMatcherDRLessByR",
    src  = cms.InputTag('slimmedElectrons'),
    matched = cms.InputTag('selectedPatTrigger'),
    matchedCuts = cms.InputTag(pathTriggerEleEle),
    maxDPtRel = cms.InputTag(0.5),
    maxDeltaR = cms.InputTag(0.2),
    resolveAmbiguities = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)

from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
removeCleaning(process)

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process , triggerMatchers = ['MuonTriggerMatchHLT','ElectronTriggerMatchHLT'])

# Muons from MuMu Tiggers
process.patMuonTriggerMatchHLTMuMu = process.patMuonTriggerMatch.clone()
process.patDefaultSequence.replace(process.patMuonTriggerMatch, process.patMuonTriggerMatchHLTMuMu)
# Electrons from EleEle Triggers
process.patElectronTriggerMatchHLTEleEle = process.patElectronTriggerMatch.clone()
process.patDefaultSequence.replace(process.patElectronTriggerMatch, process.patElectronTriggerMatchHLTEleEle)
# Muons from MuEle Triggers
process.MuonTriggerMatchHLTMuEle = process.MuonTriggerMatchHLTMuMu.clone(matchedCuts = pathTriggerMuEle)
process.patDefaultSequence.replace(process.MuonTriggerMatchHLTMuMu,
                                   process.MuonTriggerMatchHLTMuMu *
                                   process.MuonTriggerMatchHLTMuEle)
process.patMuonTriggerMatchHLTMuEle = process.patMuonTriggerMatchHLTMuMu.clone(matches = cms.VInputTag('MuonTriggerMatchHLTMuEle'))
process.patDefaultSequence.replace(process.patMuonTriggerMatchHLTMuMu,
                                   process.patMuonTriggerMatchHLTMuMu *
                                   process.patMuonTriggerMatchHLTMuEle) 
# Electrons from MuEle Triggers
process.MuonElectronTriggerMatchHLTMuEle = process.ElectronTriggerMatchHLTEleEle.clone(matchedCuts = pathTriggerMuEle)
process.patDefaultSequence.replace(process.ElectronTriggerMatchHLTEleEle, 
                                   process.ElectronTriggerMatchHLTEleEle * 
                                   process.ElectronTriggerMatchHLTMuEle)
process.patElectronTriggerMatchHLTMuEle = process.patElectronTriggerMatchHLTEleEle.clone(matches = cms.VInputTag('ElectronTriggerMatchHLTMuEle'))
process.patDefaultSequence.replace(process.patElectronTriggerMatchHLTEleEle,
                                   process.patElectronTriggerMatchHLTEleEle *
                                   process.patElectronTriggerMatchHLTMuEle) 
'''


#process.TFileService = cms.Service("TFileService", fileName = cms.string('/afs/cern.ch/work/c/cschnaib/DYanalysis/outputRoot/DYJetsToLL_events_Test1.root'))
process.TFileService = cms.Service("TFileService", fileName = cms.string('DYanalysis_50ns_ntuple.root'))

process.p = cms.Path(process.analyzer)







