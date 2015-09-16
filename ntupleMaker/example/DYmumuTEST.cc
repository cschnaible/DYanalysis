// -*- C++ -*-
//
// Package:    DYanalysis/ntupleMaker
// Class:      DYntupleMaker
// 
/**\class DYntupleMaker DYntupleMaker.cc DYanalysis/ntupleMaker/plugins/DYntupleMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christian John Schnaible
//         Created:  Mon, 13 Jul 2015 15:15:29 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>

// Root Header Files
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TSystem.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/GenMET.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"

#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

//Added photon stuff
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
//#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
// Energy Regression
//#include "HiggsAnalysis/GBRLikelihoodEGTools/interface/EGEnergyCorrectorSemiParm.h"
//#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForest.h"

//Marco Isolation
//#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"
//#include "PFIsolation/SCFootprintRemoval/interface/SuperClusterFootprintRemoval.h"

//#include "amarini/VPlusJets/interface/CiCPhotonID.h"
//JER
#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include <sys/stat.h>
#include "CLHEP/Random/RandGauss.h"

//#include "CMGTools/External/interface/PileupJetIdentifier.h"
//STD
#include <exception>

//
// class declaration
//


using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;

class DYntupleMaker : public edm::EDAnalyzer {
   public:
      explicit DYntupleMaker(const edm::ParameterSet&);
      ~DYntupleMaker();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      inline std::vector<float> getDataMCResFactor(const double& eta) const;
      bool checkTriggerName(string,std::vector<string>); //checks if string belongs to any of the vector<string>
      inline double getEffectiveAreaForElectrons(const double& ) const;//eta


      // Method that builds the tree
      void buildTree();
      // Method that re-initializes the tree branches
      void clearTree();

      // ----------member data ---------------------------
      //
      // Tokens
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<reco::GsfElectronCoreCollection> gsfelectronToken_;
      edm::EDGetTokenT<pat::TauCollection> tauToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection> genPackedCollToken_;
      edm::EDGetTokenT<reco::GenParticleCollection> genPrunedCollToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
      edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;


      
      edm::InputTag mSrcRho;//,mSrcRho25,mSrcRhoQG; //mJetsName,
      // Tree Declarations
      edm::Service<TFileService> fTFileService;
      TTree *myTree_;
      TTree *processedDataTree_;
      
      // define some event cuts
      double mMinLLMass;
      double mMaxLepEta;
      double mMinLepPt;
      double mMinLeadMuoPt;
      double mMinNjets;
      //double mJetLepIsoR;
      double mMinJetPt;
      double mMaxJetEta;
      double mMaxCombRelIso03;
      double mMaxCombRelIso04; 
      edm::InputTag pfIsoValEleCH03Name,pfIsoValEleNH03Name,pfIsoValEleG03Name;
 

      int runNum_, eventNum_, nVtx_, nLeptons_, nMuons_, nElectrons_, lumi_;
      int isRealData_; // isMC_;

      vector<float> *vtxZ_, *vtxNdof_;

      struct GENPARTICLE {
         // ---- momentum 4-vector ---------------------------------------
         TLorentzVector p4;
         // ---- pdgid ---------------------------------------------------
         int pdgId;
         // ---- motherid ---------------------------------------------------
         int motherId;
         int mmotherId;
         // ---- Charge ----------------------------
         int charge;
      };

      struct MUON {
         // Momentum 4-Vector
         TLorentzVector p4;
         float pt,px,py,pz;
         // Charge ID
         int chid;
         // tight ID
         int id;
         int muType; //muon type
         // ---- standard PF isolation uncorrected for PU ----------------
         float isoPFUnc;
         // ---- modified isolation Delta beta corrected -----------------
         float isoPFDb;
         // ---- modified isolation rho corrected ------------------------
         float isoPFRho;
         // --- flag bit
         int bit;
         // global track
         float chi2ndof;
         // muon global track
         int muonHits;
         int globalMuonHits;
         int globalHits;
         int hits;
         int pixelHits;
         int trackerLayers;

         float qoverp;
         float theta;
         float lambda;
         float dxy;
         float d0;
         float dsz;
         float dz;
         float dxyPV;
         float dzPV;
         int matchedMuStations;
         float dxyBS;
         float dzBS;
         //float dszBS;
         float vx;
         float vy;
         float vz;
         int MuonInnerTrackHits;
         int MuonLostInnerTrackHits;
         // --- Trigger Matching
         //int TriMatchF1Path, TriMatchF2Path, TriMatchF3muPath, TriMatchF3elePath, TriMatchF4Path, TriMatchF5Path,TriMatchF6Path;


      };

      struct ELECTRON {
         // Momentum 4-Vector
         TLorentzVector p4;
         float pt,px,py,pz;
         // Charge ID
         int chid;
         // tight ID
         int id;
         int elType; //muon type
         // ---- standard PF isolation uncorrected for PU ----------------
         float isoPFUnc;
         // ---- modified isolation Delta beta corrected -----------------
         float isoPFDb;
         // ---- modified isolation rho corrected ------------------------
         float isoPFRho;
         // --- flag bit
         int bit;
         // global track
         float r9;
         // --- sigma IetaIeta
         float sigmaIEtaIEta;
         float hadronicOverEm;
         float gsfTrackPt;
         float SCRawE;
         float eSCOverP;
         float dPhiSCTrackAtVtx;
         float dEtaSCTrackAtVtx;
         float epDiff;
         int MissedInnerTrackHits;
         int hasMatchedConversion;
         int EleInnerTrackHits;
         int EleLostInnerTrackHits;
         //float RegressionCorr;
         //float RegressionCorrErr;

         // --- Trigger Matching
         //int TriMatchF1Path, TriMatchF2Path, TriMatchF3muPath, TriMatchF3elePath, TriMatchF4Path, TriMatchF5Path,TriMatchF6Path;


      };

      struct JET {
         // Momentum 4 vector
         TLorentzVector p4;
         // tight id
         //int id;
         // btag info
         //float btag;
         //float taginfoNtracks;
         //float taginfoNvtx;
         //float taginfoVtxMass;
         // ---- MC flavour (0 in case of data) --------------------------
         int mcflavour;
         // ---- jet energy correction factor ----------------------------
         float jec;
         // ---- jet energy uncertainty ----------------------------------
         float unc;
         // ---- charged hadron energy fraction --------------------------
         float chf;
         // ---- neutral hadron energy fraction --------------------------
         float nhf;
         // ---- photon energy fraction ---------------------------------- 
         float phf;
         // ---- muon energy fraction ------------------------------------ 
         float muf;
         // ---- electron energy fraction --------------------------------
         float elf;
         // ---- JER corrected pt ---------------------------------------
         float rms;
      };
      
      struct GENJET {
          // pdg id
          int pdgId;
          // 4 vector
          TLorentzVector p4;
          // nparton
          //int npartons;
          // veto
          //int veto;
      };

      static bool muonSortingRule(MUON x, MUON y) {return x.pt > y.pt;}
      static bool muonSortingRuleGen(GENPARTICLE x, GENPARTICLE y) {return x.p4.Pt() > y.p4.Pt();}
      static bool p4SortingRule(TLorentzVector x, TLorentzVector y) {return x.Pt() > y.Pt();}

      TH1D  *hWEvents_; 
      /*std::string processName_;
      edm::InputTag triggerBits_;
      // ---- trigger decisions -----------------------------------------
      std::vector<int> *fired_;
      int isTriggered_;
      // ---- L1 prescale -----------------------------------------------
      std::vector<int> *prescaleL1_;
      // ---- HLT prescale -----------------------------------------------
      std::vector<int> *prescaleHLT_;*/



      // ------- Tree Variables ------------
      //
      // ------- Event Variables -----------
      float mcWeight_;
      float qScale_;
      float alphaQED_;
      float alphaQCD_;
      float x1_;
      float x2_;
      int pdf1Id_;
      int pdf2Id_;
      float scalePDF_;
      bool mOnlyMC;
      bool is25ns;
      // ---- pf pt density ---------------------------------------------
      float rho_;
      //float rho25_;

      // Trigger variables
      std::string   processName_;
      vector<int> *fired_;//, *prescaleL1_, *prescaleHLT_;
      int isTriggered_;
      int singleMu_;
      HLTConfigProvider hltConfig_;
      std::vector<std::string> triggerNames_, triggerNamesFull_;
      std::vector<std::string> triggerMuMu_;//, triggerMuEle_, triggerEleEle_;
      std::vector<std::string> triggerMu_;
      std::vector<unsigned int> triggerIndex_;
      
      // Muon Variables
      vector<int>   *muonType_,*muonId_;
      vector<float> *muonpT_,*muonEta_,*muonPhi_,*muonE_;
      vector<float> *muonIsoDBCor_,*muonIsoUnc_;
      vector<float> *muonChi2ndof_;
      vector<float> *muonQoverP_,*muonTheta_,*muonLambda_;
      vector<float> *muonDxyPV_,*muonDzPV_,*muonD0PV_;
      vector<float> *muonDxyBS_,*muonDzBS_,*muonD0BS_;
      vector<float> *muonVx_,*muonVx_,*muonVy_,*muonVz;
      vector<int>   *MuonInnerTrackHits_,*MuonLostInnerTrackHits_,*lepMissedInnerTrackHits_;
      vector<int>   *muonGlobalMuonHits_, *muonGlobalHits_, *muonHits_, *muonPixelHits_, *muonTrackerLayers_, *muonMatchedStations_;
      // DiMuon Variables
      float mmM_,mmY_,mmEta_,mmPt_,mmPhi_,mmDPhi_,mmChId_;


      // Generator Variables
      int nGenLeptons_;
      vector<float> *muonPtGEN_, *muonEtaGEN_, *muonPhiGEN_, *muonEGEN_;
      vector<int> *muonPdgIdGEN_;
      
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DYntupleMaker::DYntupleMaker(const edm::ParameterSet& iConfig):
   vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
   muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
   tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
   triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
   triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
   //triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
   genPackedCollToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
   genPrunedCollToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
   //conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
   //beamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter <edm::InputTag>("beamSpot")))
{
     //now do what ever initialization is needed
   
     edm::Service<TFileService> fs;

     mMinLLMass         = iConfig.getParameter<double>                    ("minLLMass");
     mMaxLepEta         = iConfig.getParameter<double>                    ("maxLepEta"); 
     mMinLepPt          = iConfig.getParameter<double>                    ("minLepPt");
     mMinLeadMuoPt      = iConfig.getParameter<double>                    ("minLeadMuoPt");
     mMinLepPt          = iConfig.getParameter<double>                    ("minLepPt");
     mMaxCombRelIso03   = iConfig.getParameter<double>                    ("maxCombRelIso03");
     mMaxCombRelIso04   = iConfig.getParameter<double>                    ("maxCombRelIso04");
     processName_       = iConfig.getParameter<std::string>               ("processName");
     mSrcRho            = iConfig.getParameter<edm::InputTag>             ("srcRho");
     //mSrcRho25          = iConfig.getParameter<edm::InputTag>             ("srcRho25");
     //mSrcRhoQG          = iConfig.getParameter<edm::InputTag>             ("srcRhoQG");
     triggerNames_      = iConfig.getParameter<std::vector<std::string>>  ("triggerName");
     triggerMuMu_       = iConfig.getParameter<std::vector<std::string>>  ("triggerMuMu");
     triggerMu_         = iConfig.getParameter<std::vector<std::string>>  ("triggerMu");
     mOnlyMC            = iConfig.getParameter<bool>                      ("OnlyMC"); 
     is25ns             = iConfig.getParameter<bool>                      ("is25ns");

     conversionsToken_ = mayConsume< reco::ConversionCollection >(iConfig.getParameter<edm::InputTag>("conversions"));
     beamSpotToken_    = consumes<reco::BeamSpot> (iConfig.getParameter <edm::InputTag>("beamSpot"));



  
     //ptResol.initialize(filePtResol, doGaussian);

}

// Destructor
DYntupleMaker::~DYntupleMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions for DYntupleMaker class
//


// ------------ method called for each event  ------------
void
DYntupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    clearTree(); 
    isRealData_ = iEvent.isRealData() ? 1:0;

    vector <GENJET> myGenJets;
    vector<GENPARTICLE> myGenLeptons;
    if (!isRealData_) {

        Handle<GenEventInfoProduct> geninfo;
        if(!isRealData_) {
            iEvent.getByLabel("generator",geninfo);
            mcWeight_ = geninfo->weight();
             qScale_   = geninfo->qScale();
             alphaQED_ = geninfo->alphaQED();
             alphaQCD_ = geninfo->alphaQCD();
             x1_       = geninfo->pdf()->x.first;
             x2_       = geninfo->pdf()->x.second;
             pdf1Id_   = geninfo->pdf()->id.first;
             pdf2Id_   = geninfo->pdf()->id.second;
             scalePDF_ = geninfo->pdf()->scalePDF;
             if (mcWeight_>0){
                 hWEvents_->Fill(0.5,mcWeight_);
                 hWEvents_->Fill(1.5,mcWeight_*mcWeight_);
             }
             if (mcWeight_<0) {
                 hWEvents_->Fill(2.5,-1*mcWeight_);
                 hWEvents_->Fill(3.5,mcWeight_*mcWeight_);
             }
             hWEvents_->Fill(4.5,mcWeight_);
             hWEvents_->Fill(5.5,1);
         }
         
     
         // Loop over generator particles 
         // Store 4-vector, PDG ID#
         edm::Handle<pat::PackedGenParticleCollection> genPackedColl;
         edm::Handle<reco::GenParticleCollection> genPrunedColl;
         if (!isRealData_){
             iEvent.getByToken(genPackedCollToken_, genPackedColl);
             iEvent.getByToken(genPrunedCollToken_, genPrunedColl);
         }
          
         //for (const pat::PackedGenParticle& i_gen : *genPackedColl) {
         for (const reco::GenParticle &i_gen : *genPrunedColl) {

             // Find the generated Z
              
             // ----------- Muons and Electrons -------------
             if ( abs(i_gen.pdgId()) != 13 ) continue; 
             // require stable particles  
             if ( i_gen.status() != 1 ) continue; 
             // set lepton 4-vector and ID
             GENPARTICLE aGenLepton;
             TLorentzVector lepP4GEN(i_gen.px(),i_gen.py(),i_gen.pz(),i_gen.energy());
             aGenLepton.p4 = lepP4GEN;
             aGenLepton.pdgId = i_gen.pdgId();
             // Mother particles are of type reco::GenParticle
             const reco::Candidate * motherCand = i_gen.mother(0) ;
             aGenLepton.motherId = motherCand->pdgId();
             if (fabs(aGenLepton.pdgId)==13){
                 if (aGenLepton.motherId == aGenLepton.pdgId) {
	             aGenLepton.motherId = motherCand->mother()->pdgId();
	             if (aGenLepton.motherId == aGenLepton.pdgId) {
                         aGenLepton.motherId = motherCand->mother()->mother()->pdgId();
	                 if (aGenLepton.motherId == aGenLepton.pdgId) {
	     	        aGenLepton.motherId = motherCand->mother()->mother()->mother()->pdgId();
	     	     }
	             }
	         }
             }
             aGenLepton.charge = i_gen.charge();
             // Keep leptons only within kinematic and geometric acceptance
             // if ( (aGenLepton.p4.Pt() < mMinLepPt) && (fabs(aGenLepton.p4.Eta()) > mMaxLepEta) ) continue; 
             myGenLeptons.push_back(aGenLepton);
             // Acceptance Counters (placed in example/extra.txt)
     
         } // end PackedGenParticles loop
         
    } // end if (!isRealData_)
      // Fetch pf iso maps
    Handle<ValueMap<double> > pfIsoValEleCH03;
    iEvent.getByLabel(pfIsoValEleCH03Name, pfIsoValEleCH03);
    Handle<ValueMap<double> > pfIsoValEleNH03;
    iEvent.getByLabel(pfIsoValEleNH03Name, pfIsoValEleNH03);
    Handle<ValueMap<double> > pfIsoValEleG03;
    iEvent.getByLabel(pfIsoValEleG03Name, pfIsoValEleG03);
  
    
    //---- Rho ------------------------------------------------------------
    Handle<double> rho;
    iEvent.getByLabel(mSrcRho,rho);
    rho_        = *rho;
    //---- Rho25 ------------------------------------------------------------
    //Handle<double> rho25;
    //iEvent.getByLabel(mSrcRho25,rho25);
    //rho25_      = *rho25;
    

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV is found
    const reco::Vertex &PV = vertices->front();
    for(reco::VertexCollection::const_iterator i_vtx = vertices->begin(); i_vtx != vertices->end(); ++i_vtx) {  
      // find out what this if statement is for
      if (!i_vtx->isFake() && (fabs(i_vtx->z()) < 24) && (i_vtx->ndof() >= 4)) {
	vtxZ_   ->push_back(i_vtx->z());
        vtxNdof_ ->push_back(i_vtx->ndof());
      } // end if 
    } // end vertices loop

    edm::Handle<reco::BeamSpot> beamspot;
    iEvent.getByToken(beamSpotToken_,beamspot);
    const reco::BeamSpot beamSpot = (*beamspot);
    edm::Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(conversionsToken_, conversions);

    vector<MUON> myMuons;
   
    // Loop over muons 
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    for (const pat::Muon &i_mu : *muons) {
        // Apply Kinematic and geometric acceptance for muons
        //if ( (i_mu.pt() < mMinLepPt) || fabs(i_mu.eta()) > mMaxLepEta ) continue;
        
        // pf isolation
        float muonIsoPFUnc = (i_mu.pfIsolationR03().sumChargedHadronPt + i_mu.pfIsolationR03().sumNeutralHadronEt + i_mu.pfIsolationR03().sumPhotonEt)/i_mu.pt();
        
        // Muon pf isolation rho corrected
        //float muonIsoPFRho = (i_mu.pfIsolationR03().sumChargedHadronPt
        //                  + max(0.,i_mu.pfIsolationR03().sumNeutralHadronEt + i_mu.pfIsolationR03().sumPhotonEt - (*rho25)*getEffectiveAreaForMuons(i_mu.eta())))/i_mu.pt();

        // Muon pf isolation DB corrected
        float muonIsoPFdb  = (i_mu.pfIsolationR03().sumChargedHadronPt
                          + max(0., i_mu.pfIsolationR03().sumNeutralHadronEt + i_mu.pfIsolationR03().sumPhotonEt - 0.5*i_mu.pfIsolationR03().sumPUPt))/i_mu.pt();

        //muon isolation
        //if (muonIsoPFdb > mMaxCombRelIso03) continue;

        MUON aMuon;
        TLorentzVector lepP4(i_mu.px(), i_mu.py(), i_mu.pz(), i_mu.energy() );
        aMuon.p4 = lepP4;
        aMuon.px = i_mu.px();
        aMuon.py = i_mu.py();
        aMuon.pz = i_mu.pz();
        aMuon.pt = i_mu.pt();
        aMuon.energy = i_mu.energy();
        aMuon.chid = 2*i_mu.charge();
        aMuon.id = 0;
        if (i_mu.isTightMuon(PV)){ 
            aMuon.id = 1; // =1 for tight muons
        }
        // add here for global, tracker, stand alone
        aMuon.muType = 0;
        aMuon.muType |= i_mu.isGlobalMuon() << 0;
        aMuon.muType |= i_mu.isTrackerMuon() << 1;
        aMuon.muType |= i_mu.isStandAloneMuon() << 2;
        aMuon.muType |= i_mu.isPFMuon() << 3;
        
        aMuon.isoPFUnc = muonIsoPFUnc;//(i_mu.pfIsolationR04().sumChargedHadronPt + i_mu.pfIsolationR04().sumNeutralHadronEt + i_mu.pfIsolationR04().sumPhotonEt)/i_mu.pt();
        aMuon.isoPFDb  = muonIsoPFdb;
        //aMuon.isoPFRho = muonIsoPFRho;
        // Muon Inner Track information

        aMuon.dxyPV = i_mu.muonBestTrack()->dxy(PV.position());
        aMuon.dzPV  = i_mu.muonBestTrack()->dz(PV.position());
        aMuon.dxyBS = i_mu.muonBestTrack()->dxy(beamSpot.position());
        aMuon.dzBS  = i_mu.muonBestTrack()->dz(beamSpot.position());
        //aMuon.matchedMuStations = 0;

        aMuon.matchedMuStations = i_mu.numberOfMatchedStations();
        // Muon Global Track information
        reco::TrackRef trackerTrack = i_mu.innerTrack();
        reco::TrackRef muonTrack    = i_mu.outerTrack();
        reco::TrackRef glbTrack     = i_mu.globalTrack();
        if (i_mu.globalTrack().isNonnull()) {
            const reco::HitPattern & glbhit = glbTrack->hitPattern();
            aMuon.muonHits        = glbhit.numberOfValidMuonHits();
            aMuon.globalMuonHits  = glbhit.numberOfValidMuonHits();
            aMuon.globalHits      = glbTrack->numberOfValidHits();
            aMuon.hits            = glbTrack->numberOfValidHits();
            aMuon.qoverp          = glbTrack->qoverp();
            aMuon.theta           = glbTrack->theta();
            aMuon.lambda          = glbTrack->lambda();
            aMuon.d0              = glbTrack->d0();
            aMuon.dsz             = glbTrack->dsz();
            //aMuon.dz              = glbTrack->dz(beamSpot.position());
            //aMuon.dxy             = glbTrack->dxy(beamSpot.position());
            //aMuon.dxyBS           = glbTrack->dxy(beamspot.position());
            //aMuon.dszBS           = glbTrack->dsz(beamspot.position());
            //aMuon.dzBS            = glbTrack->dz(beamspot.position());
            aMuon.vx           = glbTrack->vx();
            aMuon.vy           = glbTrack->vy();
            aMuon.vz           = glbTrack->vz();
            aMuon.chi2ndof  = i_mu.globalTrack()->normalizedChi2();

            if (trackerTrack.isNonnull()){
                const reco::HitPattern & inhit = trackerTrack->hitPattern();
                aMuon.pixelHits    = inhit.numberOfValidPixelHits();
                aMuon.trackerLayers = inhit.trackerLayersWithMeasurement();
            }

        }
        else  {
            if (trackerTrack.isNonnull()){
                aMuon.chi2ndof = trackerTrack->normalizedChi2();
                aMuon.hits           = trackerTrack->numberOfValidHits();
                if (i_mu.isTrackerMuon()){
                    //aMuon.matchedMuStations = trackerTrack.numberOfMatchedStations();
                }

                const reco::HitPattern & inhit = trackerTrack->hitPattern();
                aMuon.pixelHits    = inhit.numberOfValidPixelHits();
                aMuon.trackerLayers = inhit.trackerLayersWithMeasurement();

                if (muonTrack.isNonnull()) {
                    const reco::HitPattern & muonhit = muonTrack->hitPattern();
                    aMuon.muonHits = muonhit.numberOfValidMuonHits();
                }
                else {
                    aMuon.muonHits = 0;
                }
                
                aMuon.qoverp          = trackerTrack->qoverp();
                aMuon.theta           = trackerTrack->theta();
                aMuon.lambda          = trackerTrack->lambda();
                aMuon.d0              = trackerTrack->d0();
                aMuon.dsz             = trackerTrack->dsz();
                //aMuon.dz              = trackerTrack->dz(beamSpot.position());
                //aMuon.dxy             = trackerTrack->dxy(beamSpot.position());
                //aMuon.dxyBS           = trackerTrack->dxy(beamspot.position());
                //aMuon.dszBS           = trackerTrack->dsz(beamspot.position());
                //aMuon.dzBS            = trackerTrack->dz(beamspot.position());
                aMuon.vx           = trackerTrack->vx();
                aMuon.vy           = trackerTrack->vy();
                aMuon.vz           = trackerTrack->vz();
                aMuon.MuonInnerTrackHits = trackerTrack->numberOfValidHits();
                aMuon.MuonLostInnerTrackHits  = trackerTrack->numberOfLostHits();
            } //end if tracker track
        } // end else
    
        // Muon Trigger Matching
        // Muons from MuMu Triggers
        // Muons from MuEle Triggers

        // Muon Trigger Match information
        //aMuon.TriMatchF1Path    = MuisTriMatchF1Path;
        //aMuon.TriMatchF2Path    = 0;
        //aMuon.TriMatchF3muPath  = MuisTriMatchF3Path;
        //aMuon.TriMatchF3elePath = 0;
        //aMuon.TriMatchF5Path    = MuisTriMatchF5Path;
        //aMuon.TriMatchF6Path    = 0;

        myMuons.push_back(aMuon);
    } // end loop over muons

    // Sort leptons by Pt
    sort(myMuons.begin(), myMuons.end(), muonSortingRule);
    sort(myGenMuons.begin(), myGenMuons.end(), lepSortingRuleGen);

    // Define di-lepton pair
    TLorentzVector mmP4(0,0,0,0);
    if (myMuons.size() > 1) mmP4 = myMuons[0].p4 + myMuons[1].p4;
    TLorentzVector mmP4GEN(0,0,0,0);
    if (myGenMuon.size() > 1) mmP4GEN = myGenMuon[0].p4 + myGenMuon[1].p4;
    
    nVtx_         = int(vtxZ_->size());
    nMuons_     = int(myMuons.size()); 
    nGenMuons_  = int(myGenLeptons.size()); 
    
    bool selectionRECO(false);
     if (!mOnlyMC) {
        selectionRECO = ( (nVtx_>0) && (nLeptons_>1) && llP4.M()>mMinLLMass );
    }
    bool selection(selectionRECO);
    bool selectionGEN(false);
     
    if (!isRealData_) {
        selectionGEN = ( (nGenLeptons_>1)  && llP4GEN.M()>mMinLLMass );
        selection    += selectionGEN;
    }

    //if (selection) {
    if (selectionRECO||selectionGEN) {

        eventNum_   = iEvent.id().event();
        runNum_     = iEvent.id().run();
        lumi_       = iEvent.luminosityBlock();

        // -------------- Fill Generator Variables ----------------- //
        if (selectionGEN) {
            // Gen DiMuons
            if(myGenLeptons.size()>1) {
                llMGEN_                      = llP4GEN.M();
                llPtGEN_                     = llP4GEN.Pt();
                llPhiGEN_                    = llP4GEN.Phi();
                if(llPtGEN_>0)llYGEN_        = llP4GEN.Rapidity();
                if(llPtGEN_>0)llEtaGEN_      = llP4GEN.Eta();
                llDPhiGEN_                   = fabs(myGenLeptons[0].p4.DeltaPhi(myGenLeptons[1].p4));
            }
            // Gen Muon
            TLorentzVector lepP4GEN(0,0,0,0);
            for(unsigned l = 0; l < myGenLeptons.size(); l++) {
                lepP4GEN += myGenLeptons[l].p4;
                lepPtGEN_     ->push_back(myGenLeptons[l].p4.Pt());
                lepEtaGEN_    ->push_back(myGenLeptons[l].p4.Eta());
                lepPhiGEN_    ->push_back(myGenLeptons[l].p4.Phi());
                lepEGEN_      ->push_back(myGenLeptons[l].p4.Energy());
                lepPdgIdGEN_   ->push_back(myGenLeptons[l].pdgId);
            }
            
        } // end if (selectionGEN)

        // -------------- Fill Reco Variables ----------------- //
        if (selectionRECO) {
            // DiMuons
            if (myMuons.size() > 1) {
                mmM_                = mmP4.M();
                mmPt_               = mmP4.Pt();
                mmPhi_              = mmP4.Phi();
                if (mmPt_>0) mmY_   = mmP4.Rapidity();
                if (mmPt_>0) mmEta_ = mmP4.Eta();
                mmDPhi_             = fabs( myMuons[0].p4.DeltaPhi(myMuons[1].p4) );
                mmchid_             = myMuons[0].chid*myMuons[1].chid;
            }
            // Electron + Muon
            TLorentzVector lepP4(0,0,0,0);
            for (unsigned l = 0; l < myLeptons.size(); l++) {
                lepP4 += myLeptons[l].p4;
                lepPt_          ->push_back(myLeptons[l].p4.Pt());
                lepEta_         ->push_back(myLeptons[l].p4.Eta());
                lepPhi_         ->push_back(myLeptons[l].p4.Phi());
                lepE_           ->push_back(myLeptons[l].p4.Energy());
                lepId_          ->push_back(myLeptons[l].id);
                lepMuType_      ->push_back(myLeptons[l].muType);
                lepChId_        ->push_back(myLeptons[l].chid);
                lepPFIsoUnc_    ->push_back(myLeptons[l].isoPFUnc);
                lepPFIsoDBCor_ ->push_back(myLeptons[l].isoPFDb);
                lepPFIsoRhoCor_ ->push_back(myLeptons[l].isoPFRho);
                lepR9orChi2ndof_ ->push_back(myLeptons[l].r9_or_chi2ndof);
                lepgsfTrackPt_  ->push_back(myLeptons[l].gsfTrackPt);
                lepSCRawE_      ->push_back(myLeptons[l].SCRawE);
                lepeSCOverP_    ->push_back(myLeptons[l].eSCOverP);
                lepdPhiSCTrackAtVtx_->push_back(myLeptons[l].dPhiSCTrackAtVtx);
                lepdEtaSCTrackAtVtx_->push_back(myLeptons[l].dEtaSCTrackAtVtx);
                lepepDiff_      ->push_back(myLeptons[l].epDiff);
                lepMissedInnerTrackHits_->push_back(myLeptons[l].MissedInnerTrackHits);
                lepMuonInnerTrackHits_->push_back(myLeptons[l].MuonInnerTrackHits);
                lepEleInnerTrackHits_->push_back(myLeptons[l].EleInnerTrackHits);
                lepMuonLostInnerTrackHits_->push_back(myLeptons[l].MuonLostInnerTrackHits);
                lepEleLostInnerTrackHits_->push_back(myLeptons[l].EleLostInnerTrackHits);
                lephasMatchedConversion_->push_back(myLeptons[l].hasMatchedConversion);
                lepHadronicOverEm_     ->push_back(myLeptons[l].hadronicOverEm);
                lepSigmaIEtaIEta_      ->push_back(myLeptons[l].sigmaIEtaIEta);
                lepQoverP_             ->push_back(myLeptons[l].qoverp);
                lepTheta_              ->push_back(myLeptons[l].theta);
                lepLambda_             ->push_back(myLeptons[l].lambda);
                lepDxy_                ->push_back(myLeptons[l].dxy);
                lepDz_                 ->push_back(myLeptons[l].dz);
                lepDsz_                ->push_back(myLeptons[l].dsz);
                lepD0_                 ->push_back(myLeptons[l].d0);
                lepDxyBS_              ->push_back(myLeptons[l].dxyBS);
                lepDzBS_               ->push_back(myLeptons[l].dzBS);
                //lepDszBS_              ->push_back(myLeptons[l].dszBS);
                lepVx_                 ->push_back(myLeptons[l].vx);
                lepVy_                 ->push_back(myLeptons[l].vy);
                lepVz_                 ->push_back(myLeptons[l].vz);
                lepMuonHits_           ->push_back(myLeptons[l].muonHits);
                lepGlobalMuonHits_     ->push_back(myLeptons[l].globalMuonHits);//muon detectors hits from global muon track
                lepMuonGlobalHits_     ->push_back(myLeptons[l].globalHits);//number of hits from global muon track (unfortunate naming scheme)
                lepHits_               ->push_back(myLeptons[l].hits);
                lepMuonPixelHits_          ->push_back(myLeptons[l].pixelHits);
                lepMuonTrackerLayers_  ->push_back(myLeptons[l].trackerLayers);
                lepMuonDxyPV_          ->push_back(myLeptons[l].dxyPV);
                lepMuonDzPV_           ->push_back(myLeptons[l].dzPV);
                lepMatchedMuStations_->push_back(myLeptons[l].matchedMuStations);
                //TriMatchF1Path_doubleMu_->push_back(myLeptons[l].TriMatchF1Path);
                //TriMatchF2Path_doubleEle_->push_back(myLeptons[l].TriMatchF2Path);
                //TriMatchF3Path_MuEle_muon_->push_back(myLeptons[l].TriMatchF3muPath);
                //TriMatchF3Path_MuEle_electron_->push_back(myLeptons[l].TriMatchF3elePath);
                //TriMatchF5Path_singleMu_->push_back(myLeptons[l].TriMatchF5Path);
                //TriMatchF6Path_singleEle_->push_back(myLeptons[l].TriMatchF6Path);
            }

            // Check to see if triggers fired
            // ----  Trigger block: Bother for trigger info only if event is selected (saves time)-------------
            edm::Handle<edm::TriggerResults> triggerBits;
            iEvent.getByToken(triggerBits_,triggerBits);
            if (!triggerBits.isValid()) {
                cout << "ProcessedTreeProducer::analyze: Error in getting TriggerResults product from Event!" << endl;
                return;
            }
            // sanity check
            assert(triggerBits->size() == hltConfig_.size());
            //------ loop over all trigger names ---------
            for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
                //cout<<" trigger "<<itrig<<"/"<<triggerNames_[itrig] <<endl;
                bool accept(false);
                int tmpFired(-1);
                if (triggerIndex_[itrig] < hltConfig_.size()) {
                    accept = triggerBits->accept(triggerIndex_[itrig]);
                    string reducedTriggerName = "";
                    string reducedTriggerName2 = "";
                    int arraySize = int(triggerNamesFull_[itrig].size());
                    if(arraySize-1>0) {
                        reducedTriggerName=triggerNamesFull_[itrig].substr(0,triggerNamesFull_[itrig].size()-1); // remove last char from the str
                        reducedTriggerName2=triggerNamesFull_[itrig].substr(0,triggerNamesFull_[itrig].size()-2); // remove 2 last chars from the str
                    }
                } 
                if (!accept) {
                    tmpFired = 0;
                }
                else {
                    std::string ss(triggerNames_[itrig]);
                    //hTriggerPass_->Fill((ss.erase(ss.find("v")-1,ss.find("v"))).c_str(),1);
                    tmpFired = 1;
                    // save trigger bit (0001) if family1 has fired, (0100) if family 3 has triggered
                    isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerMuMu_) << 0; // if true 0001
                    isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerMuEle_) << 1; // if true 0010
                    isTriggered_ |= checkTriggerName(triggerNames_[itrig],triggerEleEle_) << 2; // if true 0100
                    singleMu_    |= checkTriggerName(triggerNames_[itrig],triggerMu_) << 0; // if true 0001
                }
                fired_->push_back(tmpFired);
            } // and loop over trigger names
        } //end selectionRECO
        myTree_->Fill();
    } // end selection


} // End DYntupleMaker::analyze()


// ------------ method called once each job just before starting event loop  ------------
void 
DYntupleMaker::beginJob()
{

  myTree_                = fTFileService->make<TTree>("events", "events");
  hWEvents_              = fTFileService->make<TH1D>("WEvents", "Weighted Events",6,0,6);hWEvents_->Sumw2();


  // ---- build the tree ------------------------------------------------
  buildTree();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DYntupleMaker::endJob()
{
// cout statements for acceptance placed in example/extras.txt    

}

// ------------ method called when starting to processes a run  ------------

void 
DYntupleMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
 
// Insert Trigger Things Here?
    if(!mOnlyMC){
        if (triggerNames_.size() > 0) {
            bool changed(true);
            if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
                if (changed) {
                    triggerNamesFull_.clear();
                    // check if trigger names in (new) config
                    cout<<"New trigger menu found !!!"<<endl;
                    triggerIndex_.clear();
                    const unsigned int n(hltConfig_.size());
                    for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
                        bool found(false);
                        for(unsigned iName=0;iName<n;iName++) {
                            std::string ss(hltConfig_.triggerName(iName));
                            if (int(ss.find(triggerNames_[itrig])) > -1) {
                                triggerNamesFull_.push_back(ss);
                                found = true;
                                continue;
                            } // if trigger is found in HLT config
                        } // for loop over trigger names in HLT config
                        if (!found) {
                            triggerNamesFull_.push_back("");
                        } // if trigger is not found
                        triggerIndex_.push_back(hltConfig_.triggerIndex(triggerNamesFull_[itrig]));
                        cout<<triggerNames_[itrig]<<" "<<triggerNamesFull_[itrig]<<" "<<triggerIndex_[itrig]<<" ";
                        if (triggerIndex_[itrig] >= n)
                            cout<<"does not exist in the current menu"<<endl;
                        else
                            cout<<"exists"<<endl;
                    }// trigger names loop
                } // end if (changed)
            } // end if (hltConfig_init.(iRun,iSetup,processName_,changed))
            else {
                cout << "ProcessedTreeProducer::analyze:"
                << " config extraction failure with process name "
                << processName_ << endl;
            }
        } // end if (triggernames_size>0)
    } // end if (!mOnlyMC)
}


// ------------ method called when ending the processing of a run  ------------

void 
DYntupleMaker::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------

void 
DYntupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


// ------------ method called when ending the processing of a luminosity block  ------------

void 
DYntupleMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void DYntupleMaker::buildTree() {

   // define some of the tree variables

   fired_               = new std::vector<int>();
   //prescaleL1_          = new std::vector<int>();
   //prescaleHLT_         = new std::vector<int>();

   // lepton variables
   lepPt_                    = new std::vector<float>();
   lepEta_                   = new std::vector<float>();
   lepPhi_                   = new std::vector<float>();
   lepE_                     = new std::vector<float>(); 
   lepChId_                  = new std::vector<int>();
   lepId_                    = new std::vector<int>();
   lepMuType_                = new std::vector<int>();
   lepPFIsoUnc_              = new std::vector<float>();
   lepPFIsoDBCor_          = new std::vector<float>();
   lepPFIsoRhoCor_         = new std::vector<float>();
   lepR9orChi2ndof_          = new std::vector<float>();
   //TriMatch_doubleMu_      = new std::vector<float>();
   //TriMatch_singleMu_      = new std::vector<float>();
   lepgsfTrackPt_            = new std::vector<float>();
   lepSCRawE_                = new std::vector<float>();
   lepeSCOverP_              = new std::vector<float>(); 
   lepdPhiSCTrackAtVtx_      = new std::vector<float>();
   lepdEtaSCTrackAtVtx_      = new std::vector<float>();
   lepepDiff_                = new std::vector<float>();              
   lepMissedInnerTrackHits_  = new std::vector<int>();
   lepMuonInnerTrackHits_   = new std::vector<int>(); 
   lepEleInnerTrackHits_   = new std::vector<int>(); 
   lepMuonLostInnerTrackHits_    = new std::vector<int>();  
   lepEleLostInnerTrackHits_    = new std::vector<int>();  
   lephasMatchedConversion_  = new std::vector<int>();
   lepHadronicOverEm_        = new std::vector<float>();      
   lepSigmaIEtaIEta_         = new std::vector<float>();       
   lepQoverP_                  = new std::vector<float>();
   lepTheta_                   = new std::vector<float>();
   lepLambda_                  = new std::vector<float>();
   lepDxy_                     = new std::vector<float>();
   lepDz_                      = new std::vector<float>();
   lepDsz_                     = new std::vector<float>();
   lepD0_                      = new std::vector<float>();
   lepDxyBS_                   = new std::vector<float>();
   lepDzBS_                    = new std::vector<float>();
   //lepDszBS_                   = new std::vector<float>();
   lepVx_                      = new std::vector<float>();
   lepVy_                      = new std::vector<float>();
   lepVz_                      = new std::vector<float>();
   lepMuonHits_               = new std::vector<int>();
   lepGlobalMuonHits_               = new std::vector<int>();
   lepMuonGlobalHits_         = new std::vector<int>();
   lepHits_                   = new std::vector<int>();
   lepMuonPixelHits_              = new std::vector<int>();
   lepMuonTrackerLayers_       = new std::vector<int>();
   lepMuonDxyPV_                   = new std::vector<float>();
   lepMuonDzPV_                    = new std::vector<float>();
   lepMatchedMuStations_         = new std::vector<int>();

   vtxZ_                = new std::vector<float>();
   vtxNdof_            = new std::vector<float>();

   lepPtGEN_          = new std::vector<float>();
   lepEtaGEN_         = new std::vector<float>();
   lepPhiGEN_         = new std::vector<float>();
   lepEGEN_           = new std::vector<float>();
   lepPdgIdGEN_        = new std::vector<int>();
   
   // Define myTree branches
   // Event properties
   myTree_->Branch("isRealData"       ,&isRealData_        ,"isRealData/I");
   myTree_->Branch("runNum"           ,&runNum_            ,"runNum/I");
   myTree_->Branch("eventNum"         ,&eventNum_          ,"eventNum/I");
   myTree_->Branch("lumi"             ,&lumi_              ,"lumi/I");
   myTree_->Branch("vtxZ"             ,"vector<float>"     ,&vtxZ_);
   myTree_->Branch("vtxNdof"          ,"vector<float>"     ,&vtxNdof_);
   myTree_->Branch("nVtx"             ,&nVtx_              ,"nVtx/I");
   myTree_->Branch("nMuons"           ,&nMuons_            ,"nMuons/I");
   myTree_->Branch("mcWeight"         ,&mcWeight_          ,"mcWeight/F");
   myTree_->Branch("qScale"           ,&qScale_            ,"qScale/F");
   myTree_->Branch("alphaQED"         ,&alphaQED_          ,"alphaQED/F");
   myTree_->Branch("alphaQCD"         ,&alphaQCD_          ,"alphaQCD/F");
   myTree_->Branch("x1"               ,&x1_                ,"x1/F");
   myTree_->Branch("x2"               ,&x2_                ,"x2/F");
   myTree_->Branch("pdf1Id"           ,&pdf1Id_            ,"pdf1Id/I");
   myTree_->Branch("pdf2Id"           ,&pdf2Id_            ,"pdf2Id/I");
   myTree_->Branch("scalePDF"         ,&scalePDF_          ,"scalePDF/F");

   // Reco/Trigger

   // Trigger Variables
   myTree_->Branch("fired",           "vector<int>"        ,&fired_);
   //myTree_->Branch("prescaleL1",      "vector<int>"        ,&prescaleL1_);
   //myTree_->Branch("prescaleHLT",     "vector<int>"        ,&prescaleHLT_);
   myTree_->Branch("isTriggered"      ,&isTriggered_       ,"isTriggered/I");
   myTree_->Branch("singleMuTrig"      ,&singleMu_       ,"singleMu/I");
   // DiMuon
   myTree_->Branch("mmM"              ,&mmM_               ,"mmM/F");
   myTree_->Branch("mmPt"             ,&mmPt_              ,"mmPt/F");
   myTree_->Branch("mmPhi"            ,&mmPhi_             ,"mmPhi/F");
   myTree_->Branch("mmDPhi"           ,&mmDPhi_            ,"mmDPhi/F");
   myTree_->Branch("mmY"              ,&mmY_               ,"mmY/F");
   myTree_->Branch("mmEta"            ,&mmEta_             ,"mmEta/F");
   myTree_->Branch("mmChId"           ,&mmchid_            ,"mmchid/I");
   // Lepton variables
   myTree_->Branch("lepPt"             ,"vector<float>"     ,&lepPt_);
   myTree_->Branch("lepEta"            ,"vector<float>"     ,&lepEta_);
   myTree_->Branch("lepPhi"            ,"vector<float>"     ,&lepPhi_);
   myTree_->Branch("lepE"              ,"vector<float>"     ,&lepE_);
   myTree_->Branch("lepChId"           ,"vector<int>"       ,&lepChId_);
   myTree_->Branch("lepId"             ,"vector<int>"       ,&lepId_);
   myTree_->Branch("lepMuType"         ,"vector<int>"       ,&lepMuType_);
   myTree_->Branch("lepPFIsoUnc"       ,"vector<float>"     ,&lepPFIsoUnc_);
   myTree_->Branch("lepPFIsoDBCor"    ,"vector<float>"     ,&lepPFIsoDBCor_);
   myTree_->Branch("lepPFIsoRhoCor"   ,"vector<float>"     ,&lepPFIsoRhoCor_);
   myTree_->Branch("lepR9orChi2ndof"   ,"vector<float>"     ,&lepR9orChi2ndof_);
   myTree_->Branch("lepGsfTrackPt"    ,"vector<float>"     ,&lepgsfTrackPt_);
   myTree_->Branch("lepSCRawE"        ,"vector<float>"     ,&lepSCRawE_);
   myTree_->Branch("lepSCEOverP"      ,"vector<float>"     ,&lepeSCOverP_);
   myTree_->Branch("lepDPhiSCTrackAtVtx","vector<float>"   ,&lepdPhiSCTrackAtVtx_);
   myTree_->Branch("lepDEtaSCTrackAtVtx","vector<float>"   ,&lepdEtaSCTrackAtVtx_);
   myTree_->Branch("lepEPDiff"        ,"vector<float>"     ,&lepepDiff_ );
   myTree_->Branch("lepMissedInnerTrackHits","vector<int>" ,&lepMissedInnerTrackHits_);
   myTree_->Branch("lepHasMatchedConv","vector<int>"       ,&lephasMatchedConversion_);
   myTree_->Branch("lepMuonInnerTrackHits","vector<int>"   ,&lepMuonInnerTrackHits_);
   myTree_->Branch("lepEleInnerTrackHits","vector<int>"    ,&lepEleInnerTrackHits_);
   myTree_->Branch("lepMuonLostInnerTrackHits","vector<int>"   ,&lepMuonLostInnerTrackHits_);
   myTree_->Branch("lepEleLostInnerTrackHits","vector<int>"   ,&lepEleLostInnerTrackHits_);
   myTree_->Branch("lepSigmaIEtaIEta"     ,"vector<float>"     ,&lepSigmaIEtaIEta_);
   myTree_->Branch("lepHadronicOverEm"    ,"vector<float>"     ,&lepHadronicOverEm_);
   myTree_->Branch("lepQoverP"            ,"vector<float>"     ,&lepQoverP_);
   myTree_->Branch("lepTheta"             ,"vector<float>"     ,&lepTheta_);
   myTree_->Branch("lepLambda"            ,"vector<float>"     ,&lepLambda_);
   myTree_->Branch("lepDxy"               ,"vector<float>"     ,&lepDxy_);
   myTree_->Branch("lepD0"                ,"vector<float>"     ,&lepD0_);
   myTree_->Branch("lepDz"                ,"vector<float>"     ,&lepDz_);
   myTree_->Branch("lepDsz"               ,"vector<float>"     ,&lepDsz_);
   myTree_->Branch("lepVx"                ,"vector<float>"     ,&lepVx_);
   myTree_->Branch("lepVy"                ,"vector<float>"     ,&lepVy_);
   myTree_->Branch("lepVz"                ,"vector<float>"     ,&lepVz_);
   myTree_->Branch("lepMuonHits"          ,"vector<int>"       ,&lepMuonHits_);
   myTree_->Branch("globalMuonHits"          ,"vector<int>"       ,&lepGlobalMuonHits_);
   myTree_->Branch("totalGlobalMuonHits"    ,"vector<int>"       ,&lepMuonGlobalHits_);
   myTree_->Branch("lepHits"              ,"vector<int>"       ,&lepHits_);
   myTree_->Branch("lepMuonPixelHits"     ,"vector<int>"       ,&lepMuonPixelHits_);
   myTree_->Branch("lepMuonTrackerLayers" ,"vector<int>"       ,&lepMuonTrackerLayers_);
   myTree_->Branch("lepMuonDxyPV"         ,"vector<float>"     ,&lepMuonDxyPV_);
   myTree_->Branch("lepMuonDzPV"          ,"vector<float>"     ,&lepMuonDzPV_);
   myTree_->Branch("matchedMuStations" ,"vector<int>"       ,&lepMatchedMuStations_);
   myTree_->Branch("lepdxyBS"         ,"vector<float>"     ,&lepDxyBS_);
   myTree_->Branch("lepdzBS"          ,"vector<float>"     ,&lepDzBS_);
   //myTree_->Branch("lepdszBS"         ,"vector<float>"     ,&lepDszBS_);
   //myTree_->Branch("TriMatchF1Path_doubleMu"      ,"vector<int>"        ,&TriMatchF1Path_doubleMu_);
   //myTree_->Branch("TriMatchF2Path_doubleEle"     ,"vector<int>"        ,&TriMatchF2Path_doubleEle_);
   //myTree_->Branch("TriMatchF3Path_MuEle_muon"    ,"vector<int>"        ,&TriMatchF3Path_MuEle_muon_);
   //myTree_->Branch("TriMatchF3Path_MuEle_electron","vector<int>"        ,&TriMatchF3Path_MuEle_electron_);
   //myTree_->Branch("TriMatchF5Path_singleMu"      ,"vector<int>"        ,&TriMatchF5Path_singleMu_);
   //myTree_->Branch("TriMatchF6Path_singleEle"     ,"vector<int>"        ,&TriMatchF6Path_singleEle_);

   // Jet Variables
   myTree_->Branch("jetPt"            ,"vector<float>"     ,&jetPt_);
   myTree_->Branch("jetEta"           ,"vector<float>"     ,&jetEta_);
   myTree_->Branch("jetPhi"           ,"vector<float>"     ,&jetPhi_);
   myTree_->Branch("jetE"             ,"vector<float>"     ,&jetE_);
   myTree_->Branch("jetRMS"           ,"vector<float>"     ,&jetRMS_);
   myTree_->Branch("jetMCFlavour"     ,"vector<int>"       ,&jetMCFlavour_);
   myTree_->Branch("jetJEC"           ,"vector<float>"     ,&jetJEC_);
   myTree_->Branch("jetUNC"           ,"vector<float>"     ,&jetUNC_);
   myTree_->Branch("jetllDPhi"        ,"vector<float>"     ,&jetllDPhi_);

   // Generator Variables

   // Generator Lepton Variables
   myTree_->Branch("nGenLeptons"      ,&nGenLeptons_       ,"nGenLeptons/I");
   myTree_->Branch("lepPtGEN"         ,"vector<float>"     ,&lepPtGEN_);
   myTree_->Branch("lepEtaGEN"        ,"vector<float>"     ,&lepEtaGEN_);
   myTree_->Branch("lepPhiGEN"        ,"vector<float>"     ,&lepPhiGEN_);
   myTree_->Branch("lepPdgIdGEN"      ,"vector<int>"       ,&lepPdgIdGEN_);
   myTree_->Branch("lepEGEN"          ,"vector<float>"     ,&lepEGEN_);
   myTree_->Branch("llMGEN"           ,&llMGEN_            ,"llMGEN/F");
   myTree_->Branch("llPtGEN"          ,&llPtGEN_           ,"llPtGEN/F");
   myTree_->Branch("llEtaGEN"         ,&llEtaGEN_          ,"llEtaGEN/F");
   myTree_->Branch("llPhiGEN"         ,&llPhiGEN_          ,"llPhiGEN/F");
   myTree_->Branch("llYGEN"           ,&llYGEN_            ,"llYGEN/F");
   myTree_->Branch("llDPhiGEN"        ,&llDPhiGEN_         ,"llDPhiGEN/F");
   // Generator Jet Variables
   myTree_->Branch("nJetsGEN"            ,&nJetsGEN_          ,"nJetsGEN/I");
   myTree_->Branch("jetPtGEN"            ,"vector<float>"     ,&jetPtGEN_);
   myTree_->Branch("jetEtaGEN"           ,"vector<float>"     ,&jetEtaGEN_);
   myTree_->Branch("jetPhiGEN"           ,"vector<float>"     ,&jetPhiGEN_);
   myTree_->Branch("jetEGEN"             ,"vector<float>"     ,&jetEGEN_);
   myTree_->Branch("jetllDPhiGEN"        ,"vector<float>"     ,&jetllDPhiGEN_);
   //myTree_->Branch("jetIdGEN"            ,"vector<int>"        ,&jetIdGEN_);
   //myTree_->Branch("jetVetoGEN"          ,"vector<int>"        ,&jetVetoGEN_);  
   //myTree_->Branch("jetNpartonsGEN"      ,"vector<int>"        ,&jetNpartonsGEN_);
   
}

void DYntupleMaker::clearTree() {
   
   // if vector<type> var then var->clear();
   // if float/int/double var then var = -999; 
   isRealData_ = -999;   
   runNum_     = -999;
   eventNum_   = -999;
   nVtx_       = -999;
   lumi_       = -999;
   nLeptons_   = -999;
   nMuons_     = -999;
   nElectrons_ = -999;
   nJets_      = -999;
   nJetsGEN_      = -999;
   mcWeight_          = -999;
   qScale_            = -999;
   alphaQED_          = -999;
   alphaQCD_          = -999;
   x1_                = -999;
   x2_                = -999;
   pdf1Id_            = -999;
   pdf2Id_            = -999;
   scalePDF_          = -999;
   isTriggered_       = 0;
   singleMu_       = 0;

   llM_               = -999;
   llPt_              = -999; 
   llPhi_             = -999;
   llDPhi_            = -999;
   llY_               = -999;
   llEta_             = -999;
   llchid_            = -999;
   mmM_               = -999;
   mmPt_              = -999; 
   mmPhi_             = -999;
   mmDPhi_            = -999;
   mmY_               = -999;
   mmEta_             = -999;
   mmchid_            = -999;
   eeM_               = -999;
   eePt_              = -999; 
   eePhi_             = -999;
   eeDPhi_            = -999;
   eeY_               = -999;
   eeEta_             = -999;
   eechid_            = -999;
   vtxZ_              ->clear();
   vtxNdof_              ->clear();
   //prescaleL1_        ->clear();
   //prescaleHLT_       ->clear();
   fired_             ->clear();

   lepPt_             ->clear();
   lepEta_            ->clear();
   lepPhi_            ->clear();
   lepE_              ->clear();
   lepId_             ->clear();
   lepMuonHits_       ->clear();
   lepGlobalMuonHits_       ->clear();
   lepMuonGlobalHits_       ->clear();
   lepHits_       ->clear();
   lepMuonPixelHits_  ->clear();
   lepMuonTrackerLayers_->clear();
   lepMuType_         ->clear();
   lepChId_           ->clear();
   lepPFIsoUnc_       ->clear();
   lepPFIsoDBCor_     ->clear();
   lepPFIsoRhoCor_    ->clear();
   lepR9orChi2ndof_           ->clear();
   //TriMatch_doubleMu_ ->clear();
   //TriMatch_singleMu_ ->clear();
   lepSigmaIEtaIEta_->clear();
   lepHadronicOverEm_->clear();
   lepgsfTrackPt_       -> clear();
   lepSCRawE_           -> clear();
   lepeSCOverP_         -> clear();
   lepdPhiSCTrackAtVtx_    -> clear();
   lepdEtaSCTrackAtVtx_ -> clear();
   lepepDiff_           -> clear();
   lepMissedInnerTrackHits_   -> clear();
   lephasMatchedConversion_-> clear();
   lepMuonInnerTrackHits_-> clear();
   lepEleInnerTrackHits_-> clear();
   lepMuonLostInnerTrackHits_-> clear();
   lepEleLostInnerTrackHits_-> clear();
   lepQoverP_            ->clear();
   lepTheta_             ->clear();
   lepLambda_            ->clear();
   lepDxy_               ->clear();
   lepDz_                ->clear();
   lepDsz_               ->clear();
   lepD0_                ->clear();
   lepMuonDxyPV_         ->clear();
   lepMuonDzPV_          ->clear();
   lepMatchedMuStations_ ->clear();
   lepDxyBS_             ->clear();
   lepDzBS_              ->clear();
   //lepDszBS_             ->clear();
   lepVx_                ->clear();
   lepVy_                ->clear();
   lepVz_                ->clear();
   
   nGenMuons_          = -999;
   muonPtGEN_             ->clear();
   muonEtaGEN_            ->clear();
   muonPhiGEN_            ->clear();
   muonEGEN_              ->clear();
   muonPdgIdGEN_          ->clear();
   mmMGEN_               = -999;
   mmPtGEN_              = -999; 
   mmPhiGEN_             = -999;
   mmYGEN_               = -999;
   mmDPhiGEN_            = -999;
   mmEtaGEN_             = -999;
   
} // end clearTree()

// Check Trigger Names
bool DYntupleMaker::checkTriggerName(string aString,std::vector<string> aFamily)
{
  bool result(false);
  for(unsigned int i=0;i<aFamily.size();i++) // checks if any of the aFamily strings contains aString
  {
    size_t found = aFamily[i].find(aString);
    if(found!=string::npos)result=true;
  }
  return result;
}

double DYntupleMaker::getEffectiveAreaForElectrons(const double& eta) const 
{
  double abseta = abs(eta);
    // values for 0.3 cone
    // values for 25ns runs
  if(abseta <= 1.)                            return 0.1752;
  else if(abseta > 1.    && abseta <= 1.479)  return 0.1862;
  else if(abseta > 1.479 && abseta <= 2.0)    return 0.1411;
  else if(abseta > 2.    && abseta <= 2.2)    return 0.1534;
  else if(abseta > 2.2   && abseta <=  2.3)   return 0.1903;
  else if(abseta > 2.3   && abseta <= 2.4)    return 0.2243;
  else if(abseta > 2.4   && abseta <= 2.5)    return 0.2687;
  else                                        return 9999;
}


/*
bool DYntupleMaker::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
//particle is already the ancestor
    if( ancestor == particle ) return true;

    //otherwise loop on mothers, if any and return true if the ancestor is found
    for(size_t i=0;i< particle->numberOfMothers();i++) {
        if(isAncestor(ancestor,particle->mother(i))) return true;
    }
    //if we did not return yet, then particle and ancestor are not relatives
    return false;
}
*/
//define this as a plug-in
DEFINE_FWK_MODULE(DYntupleMaker);
