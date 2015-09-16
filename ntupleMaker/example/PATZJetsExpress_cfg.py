from PhysicsTools.PatAlgos.patTemplate_cfg import *

isMC=True

process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
if(isMC):
	process.load("RecoJets.Configuration.GenJetParticles_cff")
	process.load("RecoJets.Configuration.GenJetParticles_cff")
	process.load('RecoJets.Configuration.RecoGenJets_cff')

from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.metTools import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector

# ---- access the global tag (needed for the JEC) -----------------------
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
if(isMC):
	#process.GlobalTag.globaltag = 'START53_V27::All'    # MC 53Y release > CMSSW_5_3_8_patch3
	process.GlobalTag.globaltag = 'START53_V26::All'    # MC 53Y release < CMSSW_5_3_8_patch3
else:
	process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'  # Winter13 2012 A, B, C, D datasets re-reco with CMSSW_5_3_7_patch6 
	#-----------------2012 A, B, C, D datasets re-reco + prompt with CMSSW > 5_3_2 (official recommendation) 
	#process.GlobalTag.globaltag = 'FT_53_V6C_AN3::All'  # 2012AB - July13 2012 - re-reco of 2012AB in 53X
	#process.GlobalTag.globaltag = 'FT_53_V6C_AN3::All' # 2012A - Aug06 2012 - re-reco of run-range 190782-190949 
	#process.GlobalTag.globaltag = 'FT53_V10A_AN3::All' # 2012C-v1 - Aug24 2012 - re-reco of 2012C (v1)
	#process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'  # 2012C-v2 - prompt reco for 2012C_v2 - prompt reco
	#process.GlobalTag.globaltag = 'FT_P_V42C_AN3::All' # 2012C (only 201191) - 11Dec 2012 - ecal recovery of run 201191
	#process.GlobalTag.globaltag = 'GR_P_V42_AN3::All'  # 2012D - prompt reco for 2012D - prompt reco	
	#process.GlobalTag.globaltag = 'FT_P_V43E_AN3::All' # 2012D (only range: 207883-208307) - 16Jan 2013 - recovery of 2012D run range 207883-208307

##--------- good primary vertices ---------------
process.goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    src          = cms.InputTag('offlinePrimaryVertices'),
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) )
)

from amarini.VPlusJets.hggPhotonIDCuts_cfi import *

##--------- remove cleaning --------------------
removeCleaning(process)
##--------- jets -------------------------------
process.patJets.embedPFCandidates = False
process.patJets.embedCaloTowers = False
process.patJets.addTagInfos = True

# ---- format the message service ---------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# ---- load geometry package --------------------------------------------
#process.load("Configuration.StandardSequences.Geometry_cff")
# ---- maximum number of events to run over -----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
#process.maxLuminosityBlocks = cms.untracked.PSet(input = cms.untracked.int32(1))

# ---- define the source ------------------------------------------------
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/0AA8984E-B1DB-E111-99C5-E0CB4E29C4D9.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/20A78F18-C7DB-E111-B3A2-485B39800BB4.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/2A37EB73-C1DB-E111-8370-0030487CDB24.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/2E6B20C2-E2DB-E111-9024-90E6BA19A242.root',
	#part2
	#'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/4CBED865-BBDB-E111-AFBA-BCAEC50971D0.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/50D343DC-DEDB-E111-9185-00259073E32A.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/52FC4CF0-B9DB-E111-BD2E-20CF305B0509.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/6887CC13-E9DB-E111-995A-20CF3027A5AF.root',
       #part3
       '/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/6C692933-B9DB-E111-8DAC-0030487C6A1E.root',
       '/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/706380E5-C0DB-E111-96B9-485B39897219.root',
       '/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/765D42AC-BFDB-E111-8EF7-90E6BA19A257.root',
       '/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/768EC66F-DDDB-E111-AC5C-00261834B5BE.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/7A155E7E-B6DB-E111-8915-00259073E382.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/7CF16BEB-E4DB-E111-90EF-20CF3027A574.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/848ED1CB-B8DB-E111-A29C-001EC9D825CD.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/88C0A26D-A6DB-E111-89F4-E0CB4EA0A937.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/9C335DA6-53DC-E111-9487-00259073E4D4.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/A0BD0DEA-B3DB-E111-AA32-00259073E3A6.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/A6D85956-B7DB-E111-A405-20CF305B0509.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/A8BECC30-BEDB-E111-BCE7-001E4F3F28D4.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/B269C124-57DC-E111-AD52-001EC9D7F217.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/B283CDB4-ECDB-E111-A431-485B39800C03.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/C82F09CC-B4DB-E111-9B2D-E0CB4E553658.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/CAE19663-E0DB-E111-A384-E0CB4E19F99B.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/D00569E0-BCDB-E111-8682-20CF305B0519.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/D45504AF-DBDB-E111-9F00-00259073E3F2.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/E2419980-B7DB-E111-9889-00259073E33A.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/E485CF03-86DB-E111-9000-90E6BA0D09AD.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/F614CE49-D7DB-E111-BF55-001EC9D8B52E.root',
       #'/store/mc/Summer12_DR53X/TTZJets_8TeV-madgraph_v2/AODSIM/PU_S10_START53_V7A-v1/0000/FAA1A368-D9DB-E111-93B1-E0CB4E29C4F1.root' 
 ] );

# load the PU JetID sequence
process.load("CMGTools.External.pujetidsequence_cff")
#change jet type to our
process.puJetId.jets = cms.InputTag("jetExtender",'extendedPatJets')
process.puJetMva.jets = cms.InputTag("jetExtender",'extendedPatJets')
###
#### keep the PU JetID products
###process.out.extend(["keep *_puJetId_*_*", # input variables
###		    "keep *_puJetMva_*_*" # final MVAs and working point flags
###		])

from CMGTools.External.pujetidsequence_cff import loadPujetId
#(PUIdSequence,inputsTag,mvaTags,idTags,outputCommands)=loadPujetId(process, 'selectedPatJets',mvaOnly=False,isChs=False,release="53X")
#print PUIdSequence,inputsTag,mvaTags,idTags,outputCommands

# load QGL Tagger
process.load('QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff') 
process.QGTagger.srcJets = cms.InputTag('jetExtender','extendedPatJets')
#process.QGTagger.jecService     = cms.string('') #useless for pat or corrected jets
process.QGTagger.dataDir        = cms.untracked.string("QuarkGluonTagger/EightTeV/data/")
process.QGTagger.isPatJet = cms.untracked.bool(True)


##--------- remove MC matching -----------------
if not isMC:
	removeMCMatching(process)

addPfMET(process, 'PF')

if(isMC):
	Corrections=cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
else:
	Corrections=cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'])

#BTagDiscriminators= ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags','combinedSecondaryVertexMVABJetTags','softMuonBJetTags','softMuonByPtBJetTags','softMuonByIP3dBJetTags','simpleSecondaryVertexNegativeHighEffBJetTags','simpleSecondaryVertexNegativeHighPurBJetTags']
BTagDiscriminators= ['jetBProbabilityBJetTags','jetProbabilityBJetTags','trackCountingHighPurBJetTags','trackCountingHighEffBJetTags','simpleSecondaryVertexHighEffBJetTags','simpleSecondaryVertexHighPurBJetTags','combinedSecondaryVertexBJetTags','combinedSecondaryVertexMVABJetTags']

switchJetCollection(process,cms.InputTag('ak5PFJets'),
                 doJTA              = True,
                 doBTagging         = True,
                 jetCorrLabel       = ('AK5PF', Corrections),
                 genJetCollection   = cms.InputTag('ak5GenJets'),
                 doType1MET         = False,
                 doJetID            = True,
                 btagdiscriminators = BTagDiscriminators 
                 )

process.selectedPatJets.cut = "pt > 10 && abs(eta) < 4.7"

##--------- keep only jet and MET PAT objects ---
#removeAllPATObjectsBut(process,["Jets","METs"])
if(isMC):
	process.patJets.addGenPartonMatch               = cms.bool(True)
	process.patJets.embedGenPartonMatch             = cms.bool(True)
#process.patJets.genPartonMatch  = cms.InputTag('patJetPartonMatch');

process.jetExtender = cms.EDProducer("JetExtendedProducer",
    jets    = cms.InputTag('selectedPatJets'),
    result  = cms.string('extendedPatJets'),
    payload = cms.string('AK5PF')
)

# ---- define the output file -------------------------------------------
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM_pt2.root"),
    closeFileFast = cms.untracked.bool(True)
)


# ---- Gen-Jet flavour matching -----------------------------
process.load("PhysicsTools.JetMCAlgos.SelectPartons_cff")

process.jetPartonAssociationAK5PF = cms.EDProducer("JetPartonMatcher",
               jets                = cms.InputTag("ak5PFJets"),
               coneSizeToAssociate = cms.double(0.3),
               partons             = cms.InputTag("myPartons")
               )

process.jetFlavourAssociationAK5PF = cms.EDProducer("JetFlavourIdentifier",
                srcByReference    = cms.InputTag("patJetPartonAssociationAK5PF"),
                physicsDefinition = cms.bool(False)
                )

process.GenJetFlavourMatching = cms.Sequence(process.myPartons*process.jetPartonAssociationAK5PF*process.jetFlavourAssociationAK5PF)


# ---- Recompute electron's PF iso deposits -----------------------------
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(False)
  )

from JetMETCorrections.Type1MET.pfMETCorrections_cff import *
if not isMC:
	process.pfJetMETcorr.jetCorrLabel = "ak5PFL1FastL2L3Residual"

# ---- ZJetsExpress analyzer --------------------------------------------
process.accepted = cms.EDAnalyzer('PATZJetsExpress',
    jets            = cms.InputTag('jetExtender','extendedPatJets'),
    srcRho          = cms.InputTag('kt6PFJets','rho'),
    srcRho25        = cms.InputTag('kt6PFJetsCentralNeutral','rho'),
    srcRhoQG        = cms.InputTag('kt6PFJetsIsoQG','rho'),
    pfIsoValEleCH03 = cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
    pfIsoValEleNH03 = cms.InputTag('elPFIsoValueNeutral03PFIdPFIso'),
    pfIsoValEleG03  = cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
    PfMETType1      = cms.InputTag('pfType1CorrectedMet'),
    #PfMETType01     = cms.InputTag('pfMetT0pcT1'),
    GenMETTrue      = cms.InputTag('genMetTrue'),                         
    minNjets        = cms.int32(1),
    jetLepIsoRadius = cms.double(0.3),
    jetLepPhoRadius = cms.double(0.3),
    minJetPt        = cms.double(25),
    maxJetEta       = cms.double(2.5),
    #photons are saved from minPhoPtId if they pass
    #CaloId in addition, minPhoPt all photons are saved
    minPhoPt        = cms.double(50),
    minPhoPtId      = cms.double(30),
    maxPhoEta       = cms.double(3.0),
    minLepPt        = cms.double(18),
    minNuPt         = cms.double(10),			  
    maxLepEta       = cms.double(2.4),
    maxCombRelIso03 = cms.double(0.15),
    maxCombRelIso04 = cms.double(0.12),
    minLLMass       = cms.double(15),
    OnlyMC          = cms.bool(False),
    ReducedPh       = cms.bool(False), #if TRUE: reduce photon information to 4-Vector  
    dressedRadius   = cms.double(0.1),
    #GENCrossCleaning= cms.int32(1), #OBSOLETE
    GENType         = cms.int32(0), #1 for dressed level, 0 for gen default
    processName     = cms.string('HLT'),
    ##SuperCluster Foot Print Removal
    tag_jets        = cms.InputTag('jetExtender','extendedPatJets'),
    #isolation_cone_size_forSCremoval = cms.double(0.4),#default=0.4 

    triggerName     = cms.vstring('HLT_DoubleMu7_v',
                                  'HLT_Mu13_Mu8_v',
                                  'HLT_Mu17_Mu8_v',
                                  'HLT_Mu17_TkMu8_v', #end 2011
                                  'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v', # end 2011
                                  'HLT_Mu17_Ele8_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
                                  'HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v',
                                  'HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                  'HLT_Mu30_Ele30_CaloIdL_v',
                                  'HLT_Mu7_Ele7_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                  'HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v',
                                  'HLT_Photon20_CaloIdVL_v',
                                  'HLT_Photon20_CaloIdVL_IsoL_v',
                                  'HLT_Photon30_v',
                                  'HLT_Photon30_CaloIdVL_v',
                                  'HLT_Photon30_CaloIdVL_IsoL_v',
                                  'HLT_Photon50_CaloIdVL_v',
                                  'HLT_Photon50_CaloIdVL_IsoL_v',
                                  'HLT_Photon75_CaloIdVL_v',
                                  'HLT_Photon75_CaloIdVL_IsoL_v',
                                  'HLT_Photon90_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_IsoL_v',
                                  'HLT_Photon135_v',
                                  'HLT_Photon150_v',
                                  'HLT_Mu15_v',
                                  'HLT_Mu24_v',
                                  'HLT_Mu30_v',
                                  'HLT_Mu40_v',
                                  'HLT_Mu40_eta2p1_v',
                                  'HLT_IsoMu17_v',
                                  'HLT_IsoMu20_v',
                                  'HLT_IsoMu24_v',
                                  'HLT_IsoMu24_eta2p1_v', # end 2011
                                  'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v',
                                  'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v',
                                  'HLT_Ele52_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele65_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele80_CaloIdVT_TrkIdT_v', # end of 2011
                                  'HLT_Ele27_WP80_v',
				  'HLT_SingleE_v',
				  'HLT_HT200_v',
				  'HLT_HT250_v',
                                  'HLT_HT300_v',
                                  'HLT_HT350_v',
                                  'HLT_HT400_v',
                                  'HLT_HT450_v',
                                  'HLT_HT500_v',
                                  'HLT_HT550_v',
                                  'HLT_HT600_v',
                                  'HLT_HT650_v',
                                  'HLT_HT700_v',
                                  'HLT_HT750_v'
                            ),                     
    triggerResults  = cms.InputTag("TriggerResults","","HLT"),
    triggerEvent    = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    triggerFamily1  = cms.vstring('HLT_DoubleMu7_v',
                                  'HLT_Mu13_Mu8_v',
                                  'HLT_Mu17_Mu8_v',
                                  'HLT_Mu17_TkMu8_v'),
    triggerFamily2  = cms.vstring('HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v',
                                  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v'),
    triggerFamily3  = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
                                  'HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v',
                                  'HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                  'HLT_Mu30_Ele30_CaloIdL_v',
                                  'HLT_Mu7_Ele7_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                  'HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v'),
    triggerFamily4  = cms.vstring('HLT_Photon20_CaloIdVL_v',
                                  'HLT_Photon20_CaloIdVL_IsoL_v',
                                  'HLT_Photon30_v',
                                  'HLT_Photon30_CaloIdVL_v',
                                  'HLT_Photon30_CaloIdVL_IsoL_v',
                                  'HLT_Photon50_CaloIdVL_v',
                                  'HLT_Photon50_CaloIdVL_IsoL_v',
                                  'HLT_Photon75_CaloIdVL_v',
                                  'HLT_Photon75_CaloIdVL_IsoL_v',
                                  'HLT_Photon90_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_IsoL_v',
                                  'HLT_Photon135_v',
                                  'HLT_Photon150_v'),
    triggerFamily5  = cms.vstring('HLT_Mu15_v',
                                  'HLT_Mu24_v',
                                  'HLT_Mu30_v',
                                  'HLT_Mu40_v',
                                  'HLT_Mu40_eta2p1_v',
                                  'HLT_IsoMu17_v',
                                  'HLT_IsoMu20_v',
                                  'HLT_IsoMu24_v',
                                  'HLT_IsoMu24_eta2p1_v'), # end 2011
    triggerFamily6  = cms.vstring('HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v',
                                  'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v',
                                  'HLT_Ele52_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele65_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele80_CaloIdVT_TrkIdT_v', # end of 2011
                                  'HLT_Ele27_WP80_v'), # end 2011
    triggerFamily7  = cms.vstring('HLT_SingleE_v'), ##TO BE CHECKED!
    triggerFamily8  = cms.vstring('HLT_HT200_v',
				  'HLT_HT250_v',
                                  'HLT_HT300_v',
                                  'HLT_HT350_v',
                                  'HLT_HT400_v',
                                  'HLT_HT450_v',
                                  'HLT_HT500_v',
                                  'HLT_HT550_v',
                                  'HLT_HT600_v',
                                  'HLT_HT650_v',
                                  'HLT_HT700_v',
                                  'HLT_HT750_v'), ##TO BE CHECKED!

    hggPhotonIDConfiguration = cms.PSet(hggPhotonIDCuts),
    prescaleDontAsk = cms.vstring('HLT_Mu17_Ele8_CaloIdL_v', # don't ask for L1 prescales for these bits
                                  'HLT_Mu8_Ele17_CaloIdL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v',
                                  'HLT_Mu17_TkMu8_v',
                                  'HLT_Mu13_Mu8_v',
                                  'HLT_Mu17_Mu8_v',
                                  'HLT_Ele27_WP80_v',
				  'HLT_Photon30_v',
                                  'HLT_Photon30_CaloIdVL_v',
                                  'HLT_Photon30_CaloIdVL_IsoL_v',
                                  'HLT_Photon50_CaloIdVL_v',
                                  'HLT_Photon50_CaloIdVL_IsoL_v',
                                  'HLT_Photon75_CaloIdVL_v',
                                  'HLT_Photon75_CaloIdVL_IsoL_v',
                                  'HLT_Photon90_CaloIdVL_v',
                                  'HLT_Photon90_CaloIdVL_IsoL_v',
                                  'HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v',
                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
                                  'HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v',
                                  'HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v',
                                  'HLT_Ele65_CaloIdVT_TrkIdT_v',
                                  'HLT_Ele80_CaloIdVT_TrkIdT_v',
				  'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
				  'HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v',
				  'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v',
				  'HLT_Mu30_Ele30_CaloIdL_v',
				  'HLT_Mu7_Ele7_CaloIdT_CaloIsoVL_v',
				  'HLT_Photon135_v',
				  'HLT_Photon150_v',
				  'HLT_Mu24_v',
				  'HLT_Mu30_v',
				  'HLT_Mu40_v',
				  'HLT_Mu40_eta2p1_v',
				  'HLT_IsoMu24_v',
				  'HLT_IsoMu24_eta2p1_v',
				  'HLT_HT200_v',
				  'HLT_HT250_v',
                                  'HLT_HT300_v',
                                  'HLT_HT350_v',
                                  'HLT_HT400_v',
                                  'HLT_HT450_v',
                                  'HLT_HT500_v',
                                  'HLT_HT550_v',
                                  'HLT_HT600_v',
                                  'HLT_HT650_v',
                                  'HLT_HT700_v',
                                  'HLT_HT750_v'),
)
# ---- duplicate the analyzer with different name -----------------------
process.rejected = process.accepted.clone()
# ---- filter the required HLT bits -------------------------------------
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag  = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths           = cms.vstring('HLT_DoubleMu7_v*',
                                     'HLT_Mu13_Mu8_v*',
                                     'HLT_Mu17_Mu8_v*',
                                     'HLT_Mu17_TkMu8_v*', #end 2011
                                     'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
                                     'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*', # end 2011
                                     'HLT_Mu17_Ele8_CaloIdL_v*',
                                     'HLT_Mu8_Ele17_CaloIdL_v*',
                                     'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*',
                                     'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*',
                                     'HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v*',
                                     'HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v*',
                                     'HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
                                     'HLT_Mu30_Ele30_CaloIdL_v*',
                                     'HLT_Mu7_Ele7_CaloIdT_CaloIsoVL_v*',
                                     'HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v*',
                                     'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
                                     'HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v*',
                                     'HLT_Photon20_CaloIdVL_v*',
                                     'HLT_Photon20_CaloIdVL_IsoL_v*',
                                     'HLT_Photon30_v*',
                                     'HLT_Photon30_CaloIdVL_v*',
                                     'HLT_Photon30_CaloIdVL_IsoL_v*',
                                     'HLT_Photon50_CaloIdVL_v*',
                                     'HLT_Photon50_CaloIdVL_IsoL_v*',
                                     'HLT_Photon75_CaloIdVL_v*',
                                     'HLT_Photon75_CaloIdVL_IsoL_v*',
                                     'HLT_Photon90_CaloIdVL_v*',
                                     'HLT_Photon90_CaloIdVL_IsoL_v*',
                                     'HLT_Photon135_v*',
                                     'HLT_Photon150_v*', 
                                     'HLT_Mu15_v2*',
                                     'HLT_Mu24_v*',
                                     'HLT_Mu30_v*',
                                     'HLT_Mu40_v*',
                                     'HLT_Mu40_eta2p1_v*',
                                     'HLT_IsoMu17_v*',
                                     'HLT_IsoMu20_v*',
                                     'HLT_IsoMu24_v*',
                                     'HLT_IsoMu24_eta2p1_v*', # end 2011
                                     'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                     'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
                                     'HLT_Ele52_CaloIdVT_TrkIdT_v*',
                                     'HLT_Ele65_CaloIdVT_TrkIdT_v*',
                                     'HLT_Ele80_CaloIdVT_TrkIdT_v*', # end of 2011
                                     'HLT_Ele27_WP80_v*'),
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #----- True = OR, False = AND between the HLTPaths
    throw              = cms.bool(False)
)
############# turn-on the fastjet area calculation needed for the L1Fastjet ##############
process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True
############# turn-on the fastjet area in |eta|<2.5 needed for the photonISO #############
process.kt6PFJets25 = process.kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)

############# Trigger Matching #############
pathTriggerMuonsFamily1 = '(path("HLT_DoubleMu7_v*", 1, 1) || path("HLT_Mu13_Mu8_v*", 1, 1) || path("HLT_Mu17_Mu8_v*", 1, 1) || path("HLT_Mu17_TkMu8_v*", 1, 1))' # selecting the trigger objects

pathTriggerElectronsFamily2 = '(path("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*", 1, 1) || path("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", 1, 1))' 

pathTriggerMuonElectronFamily3 = '(path("HLT_Mu17_Ele8_CaloIdL_v*", 1, 1) || path("HLT_Mu8_Ele17_CaloIdL_v*", 1, 1) || path("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*", 1, 1) || path("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v*", 1, 1) || path("HLT_DoubleMu5_Ele8_CaloIdT_TrkIdVL_v*", 1, 1) || path("HLT_DoubleMu8_Ele8_CaloIdT_TrkIdVL_v*", 1, 1) || path("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", 1, 1) || path("HLT_Mu30_Ele30_CaloIdL_v*", 1, 1) || path("HLT_Mu7_Ele7_CaloIdT_CaloIsoVL_v*", 1, 1) || path("HLT_Mu8_DoubleEle8_CaloIdT_TrkIdVL_v*", 1, 1) || path("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", 1, 1) || path("HLT_Mu8_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v*", 1, 1))' 

pathTriggerPhotonsFamily4 = '(path("HLT_Photon20_CaloIdVL_v*", 1, 1) || path("HLT_Photon20_CaloIdVL_IsoL_v*", 1, 1) || path("HLT_Photon30_v*", 1, 1) || path("HLT_Photon30_CaloIdVL_v*", 1, 1) || path("HLT_Photon30_CaloIdVL_IsoL_v*", 1, 1)|| path("HLT_Photon50_CaloIdVL_v*", 1, 1) || path("HLT_Photon50_CaloIdVL_IsoL_v*", 1, 1) || path("HLT_Photon75_CaloIdVL_v*", 1, 1) || path("HLT_Photon75_CaloIdVL_IsoL_v*", 1, 1) || path("HLT_Photon90_CaloIdVL_v*", 1, 1) || path("HLT_Photon90_CaloIdVL_IsoL_v*", 1, 1) || path("HLT_Photon135_v*", 1, 1)|| path("HLT_Photon150_v*", 1, 1))'

pathTriggerMuonsFamily5 = '(path("HLT_Mu15_v2*", 1, 1) || path("HLT_Mu24_v*", 1, 1) || path("HLT_Mu30_v*", 1, 1) || path("HLT_Mu40_v*", 1, 1) || path("HLT_Mu40_eta2p1_v*", 1, 1) || path("HLT_IsoMu17_v*", 1, 1) || path("HLT_IsoMu20_v*", 1, 1) || path("HLT_IsoMu24_v*", 1, 1) || path("HLT_IsoMu24_eta2p1_v*", 1, 1))' # selecting the trigger objects

pathTriggerElectronsFamily6 = '(path("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*", 1, 1) || path("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*", 1, 1) || path("HLT_Ele52_CaloIdVT_TrkIdT_v*", 1, 1) || path("HLT_Ele65_CaloIdVT_TrkIdT_v*", 1, 1) || path("HLT_Ele80_CaloIdVT_TrkIdT_v*", 1, 1) || path("HLT_Ele27_WP80_v*", 1, 1))' 

process.MuonsTriggerMatchHLTFamily1 = cms.EDProducer(
    "PATTriggerMatcherDRLessByR", 
    src                   = cms.InputTag('patMuons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string(pathTriggerMuonsFamily1),
    maxDPtRel             = cms.double(0.5),
    maxDeltaR             = cms.double(0.5),
    # maxDeltaEta         = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)

process.ElectronsTriggerMatchHLTFamily2 = cms.EDProducer(
    "PATTriggerMatcherDRLessByR",
    src                   = cms.InputTag('patElectrons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string(pathTriggerElectronsFamily2),
    maxDPtRel             = cms.double(0.5),
    maxDeltaR             = cms.double(0.5),
    # maxDeltaEta         = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)

process.PhotonsTriggerMatchHLTFamily4 = cms.EDProducer(
    "PATTriggerMatcherDRLessByR", 
    src                   = cms.InputTag('patPhotons'),
    matched               = cms.InputTag('patTrigger'),
    matchedCuts           = cms.string(pathTriggerPhotonsFamily4),
    maxDPtRel             = cms.double(1.0),
    maxDeltaR             = cms.double(0.2),
    # maxDeltaEta         = cms.double(0.2),
    resolveAmbiguities    = cms.bool(True),
    resolveByMatchQuality = cms.bool(True)
)

from PhysicsTools.PatAlgos.tools.coreTools import removeCleaning
removeCleaning(process)

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process , triggerMatchers = ['MuonsTriggerMatchHLTFamily1', 'ElectronsTriggerMatchHLTFamily2', 'PhotonsTriggerMatchHLTFamily4'])

# Rename existing:
process.patMuonsTriggerMatchHLTFamily1 = process.patMuonsTriggerMatch.clone()
process.patDefaultSequence.replace(process.patMuonsTriggerMatch, process.patMuonsTriggerMatchHLTFamily1)
#process.out.outputCommands.append('keep *_patMuonsTriggerMatchHLTFamily1_*_*')

process.patElectronsTriggerMatchHLTFamily2 = process.patElectronsTriggerMatch.clone()
process.patDefaultSequence.replace(process.patElectronsTriggerMatch, process.patElectronsTriggerMatchHLTFamily2)
#process.out.outputCommands.append('keep *_patElectronsTriggerMatchHLTFamily2_*_*')

process.patPhotonsTriggerMatchHLTFamily4 = process.patPhotonsTriggerMatch.clone()
process.patDefaultSequence.replace(process.patPhotonsTriggerMatch, process.patPhotonsTriggerMatchHLTFamily4)
#process.out.outputCommands.append('keep *_patPhotonsTriggerMatchHLTFamily4_*_*')

# Clone, what is needed for second matcher:
#Muons
process.MuonElectronTriggerMatchHLTFamily3mu = process.MuonsTriggerMatchHLTFamily1.clone(matchedCuts =
                                                                                        pathTriggerMuonElectronFamily3)
process.MuonsTriggerMatchHLTFamily5 = process.MuonsTriggerMatchHLTFamily1.clone(matchedCuts = pathTriggerMuonsFamily5)
process.patDefaultSequence.replace(process.MuonsTriggerMatchHLTFamily1, 
                                   process.MuonsTriggerMatchHLTFamily1 * 
                                   process.MuonElectronTriggerMatchHLTFamily3mu * 
                                   process.MuonsTriggerMatchHLTFamily5)
                                   
process.patMuonElectronTriggerMatchHLTFamily3mu = process.patMuonsTriggerMatchHLTFamily1.clone(matches = 
                                                                           cms.VInputTag('MuonElectronTriggerMatchHLTFamily3mu'))
process.patMuonsTriggerMatchHLTFamily5 = process.patMuonsTriggerMatchHLTFamily1.clone(matches = 
                                                                           cms.VInputTag('MuonsTriggerMatchHLTFamily5'))
process.patDefaultSequence.replace(process.patMuonsTriggerMatchHLTFamily1, 
                                   process.patMuonsTriggerMatchHLTFamily1 * 
                                   process.patMuonElectronTriggerMatchHLTFamily3mu *
                                   process.patMuonsTriggerMatchHLTFamily5)
#process.out.outputCommands.append('keep *_patMuonElectronTriggerMatchHLTFamily3mu_*_*')
#process.out.outputCommands.append('keep *_patMuonsTriggerMatchHLTFamily5_*_*')

#Electrons
process.MuonElectronTriggerMatchHLTFamily3ele = process.ElectronsTriggerMatchHLTFamily2.clone(matchedCuts =
                                                                                        pathTriggerMuonElectronFamily3)
process.ElectronsTriggerMatchHLTFamily6 = process.ElectronsTriggerMatchHLTFamily2.clone(matchedCuts = pathTriggerElectronsFamily6)
process.patDefaultSequence.replace(process.ElectronsTriggerMatchHLTFamily2, 
                                   process.ElectronsTriggerMatchHLTFamily2 * 
                                   process.MuonElectronTriggerMatchHLTFamily3ele * 
                                   process.ElectronsTriggerMatchHLTFamily6)
                                   
process.patMuonElectronTriggerMatchHLTFamily3ele = process.patElectronsTriggerMatchHLTFamily2.clone(matches = 
                                                                           cms.VInputTag('MuonElectronTriggerMatchHLTFamily3ele'))
process.patElectronsTriggerMatchHLTFamily6 = process.patElectronsTriggerMatchHLTFamily2.clone(matches = 
                                                                           cms.VInputTag('ElectronsTriggerMatchHLTFamily6'))
process.patDefaultSequence.replace(process.patElectronsTriggerMatchHLTFamily2, 
                                   process.patElectronsTriggerMatchHLTFamily2 * 
                                   process.patMuonElectronTriggerMatchHLTFamily3ele *
                                   process.patElectronsTriggerMatchHLTFamily6)
#process.out.outputCommands.append('keep *_patMuonElectronTriggerMatchHLTFamily3ele_*_*')
#process.out.outputCommands.append('keep *_patElectronsTriggerMatchHLTFamily6_*_*')

process.dump = cms.EDAnalyzer("EventContentAnalyzer")

del process.outpath


# ---- save all events for any trigger ---------------------
process.p = cms.Path(process.pfParticleSelectionSequence
                     + process.eleIsoSequence
                     + process.phoIsoSequence
                     + process.kt6PFJets
                     + process.ak5PFJets
                     + process.kt6PFJets25
                     + process.goodOfflinePrimaryVertices)

if(isMC):
	process.p += process.genParticlesForJets
process.tail = cms.Sequence(process.patDefaultSequence + process.jetExtender + 
		process.puJetIdSqeuence + 
		#PUIdSequence + 
		process.QuarkGluonTagger +
		process.producePFMETCorrections+
		#process.dump*
		process.accepted)

process.p += process.tail



