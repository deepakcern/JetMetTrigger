// -*- C++ -*-
//
// Package:    HLTJetMETNtupleProducer
// Class:      HLTJetMETNtupleProducer
// 
/**\class HLTJetMETNtupleProducer HLTJetMETNtupleProducer.cc HLTStudies/HLTFilter/src/HLTJetMETNtupleProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arun Nayak
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
//
// class declaration
//
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h" 
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h" 
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include <iostream>
#include <vector>
#include "TTree.h"
#include <string>


class HLTJetMETNtupleProducer : public edm::EDAnalyzer {
   public:
      explicit HLTJetMETNtupleProducer(const edm::ParameterSet&);
      ~HLTJetMETNtupleProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      bool isICHEPMuon(const reco::Muon&);

      // ----------member data ---------------------------
  //edm::InputTag hlt_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection> genParticleToken_;
  edm::EDGetTokenT<pat::METCollection> MetCollectionToken_;
  edm::EDGetTokenT<reco::GenMETCollection> GenMetCollectionToken_;
  edm::EDGetTokenT<pat::MuonCollection> MuonCollectionToken_;
  edm::EDGetTokenT<edm::View<pat::Electron> > ElectronCollectionToken_;
  edm::EDGetTokenT<reco::VertexCollection> PVToken_;
  edm::EDGetTokenT<pat::JetCollection> PFJetCollectionToken_;
  edm::EDGetTokenT<reco::CaloJetCollection> CaloJetCollectionToken_;
  edm::EDGetTokenT<reco::GenJetCollection> GenJetCollectionToken_;
  edm::EDGetTokenT<reco::PFJetCollection> HLTPFJetCollectionToken_;
  edm::EDGetTokenT<reco::CaloJetCollection> HLTCaloJetCollectionToken_;
  edm::EDGetTokenT<reco::PFMETCollection> HLTPFMetCollectionToken_;
  edm::EDGetTokenT<reco::PFMETCollection> HLTPFMetType1CollectionToken_;
  edm::EDGetTokenT<reco::CaloMETCollection> HLTCaloMetCollectionToken_;
  edm::EDGetTokenT<edm::TriggerResults> hlt_;
  edm::EDGetTokenT<edm::TriggerResults> hlt_rerun;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  std::vector<std::string> triggerPaths_;
  //// e-ID decisions objects
  /* edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaNonTrigWP80MapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaNonTrigWP90MapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaTrigWP80MapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaTrigWP90MapToken_;
  //// MVA values and categories
  edm::EDGetTokenT<edm::ValueMap<float> > mvaNonTrigValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> >   mvaNonTrigCategoriesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaTrigValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> >   mvaTrigCategoriesMapToken_;
  //Spring16 MVA
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaSpring16WPMediumMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaSpring16WPTightMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaSpring16ValuesMapToken_;
  edm::EDGetTokenT<edm::ValueMap<int> > mvaSpring16CategoriesMapToken_;
  //Summer16 cut-based*/
  edm::EDGetTokenT<edm::ValueMap<bool> > eleFall17VetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleFall17LooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleFall17MediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleFall17TightIdMapToken_;


  /* edm::EDGetTokenT<edm::ValueMap<bool> > eleMvanoIsoIdWP80MapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvanoIsoIdWP90MapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvanoIsoIdwpLooseMapToken_;*/
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaIsoIdWP80MapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaIsoIdWP90MapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > eleMvaIsoIdwpLooseMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > mvaFall17IsoValuesMapToken_;
  // edm::EDGetTokenT<edm::ValueMap<float> > mvaFall17noIsoValuesMapToken_;
  bool runJets_, runMets_, isData_;

  TTree* tree_;
  unsigned long run_,event_,lumi_;
  float metPx_, metPy_, metPt_, metPhi_, caloMetPt_, caloMetPhi_;
  float genMetPt_, genMetPhi_;
  float hltPFMetPt_, hltPFMetPhi_, hltPFMetType1Pt_, hltPFMetType1Phi_;
  float hltCaloMetPt_, hltCaloMetPhi_;
  bool passMETFilter_;
  std::vector<std::string> triggerResults_;
  std::vector<std::string> triggerResults_rerun;
  UInt_t  nPV_;
  Float_t PVx_;
  Float_t PVy_;
  Float_t PVz_;

  UInt_t passedHLTCaloMETClean_;
  UInt_t passedHLTCaloMET_;
  UInt_t passedL1MET_;

  //gen Jets
  std::vector<float> genjetPt_;
  std::vector<float> genjetEta_;
  std::vector<float> genjetPhi_;
  //pat jets
  std::vector<float> pfjetPt_;
  std::vector<float> pfjetEta_;
  std::vector<float> pfjetPhi_;
  std::vector<float> calojetPt_;
  std::vector<float> calojetEta_;
  std::vector<float> calojetPhi_;
  //hlt jets
  std::vector<float> hltpfjetPt_;
  std::vector<float> hltpfjetEta_;
  std::vector<float> hltpfjetPhi_;
  std::vector<float> hltcalojetPt_;
  std::vector<float> hltcalojetEta_;
  std::vector<float> hltcalojetPhi_;

  // pat muons
  std::vector<float> muonPx_;
  std::vector<float> muonPy_;
  std::vector<float> muonPz_;
  std::vector<float> muonPt_;
  std::vector<float> muonEta_;
  std::vector<float> muonPhi_;
  std::vector<float> muonCharge_;
  std::vector<float> muonR04SumChargedHadronPt_;
  std::vector<float> muonR04SumChargedParticlePt_;
  std::vector<float> muonR04SumNeutralHadronEt_;
  std::vector<float> muonR04SumPhotonEt_;
  std::vector<float> muonR04SumPUPt_;
  std::vector<bool> muonIsPF_;
  std::vector<bool> muonIsGlobal_;
  std::vector<bool> muonIsTracker_;
  std::vector<bool> muonIsICHEPMedium_;
  std::vector<float> muonDz_;
  std::vector<float> muonDxy_;
  //std::vector<float> muonNormChi2_;

  //pat electrons
  std::vector<float> elecPx_;
  std::vector<float> elecPy_;
  std::vector<float> elecPz_;
  std::vector<float> elecPt_;
  std::vector<float> elecEta_;
  std::vector<float> elecPhi_;
  std::vector<float> elecCharge_;
  std::vector<float> elecDz_;
  std::vector<float> elecDxy_;
  std::vector<float> elecR03SumChargedHadronPt_;
  std::vector<float> elecR03SumChargedParticlePt_;
  std::vector<float> elecR03SumNeutralHadronEt_;
  std::vector<float> elecR03SumPhotonEt_;
  std::vector<float> elecR03SumPUPt_;
  //std::vector<float> elec_mva_value_nontrig_Spring15_v1_;
  //std::vector<float> elec_mva_value_trig_Spring15_v1_;
  //std::vector<int> elec_mva_category_nontrig_Spring15_v1_;
  //std::vector<int> elec_mva_category_trig_Spring15_v1_;
  //std::vector<bool> elec_mva_wp80_nontrig_Spring15_v1_;
  //std::vector<bool> elec_mva_wp90_nontrig_Spring15_v1_;
  //std::vector<bool> elec_mva_wp80_trig_Spring15_v1_;
  //std::vector<bool> elec_mva_wp90_trig_Spring15_v1_;
  //std::vector<bool> elec_cutId_veto_Spring15_;
  //std::vector<bool> elec_cutId_loose_Spring15_;
  //std::vector<bool> elec_cutId_medium_Spring15_;
  //std::vector<bool> elec_cutId_tight_Spring15_;
  //std::vector<bool> elec_mva_medium_Spring16_v1_;
  //std::vector<bool> elec_mva_tight_Spring16_v1_;
  //std::vector<float> elec_mva_value_Spring16_v1_;
  //std::vector<int> elec_mva_category_Spring16_v1_;
  std::vector<bool> elec_cutId_veto_Fall17_;
  std::vector<bool> elec_cutId_loose_Fall17_;
  std::vector<bool> elec_cutId_medium_Fall17_;
  std::vector<bool> elec_cutId_tight_Fall17_;
  /* std::vector<bool> elec_mva_wp80_noIso_Fall17_;
  std::vector<bool> elec_mva_wp90_noIso_Fall17_; 
  std::vector<bool> elec_mva_wpLoose_noIso_Fall17_;*/
  std::vector<bool> elec_mva_wp80_Iso_Fall17_;
  std::vector<bool> elec_mva_wp90_Iso_Fall17_;
  std::vector<bool> elec_mva_wpLoose_Iso_Fall17_;
  std::vector<float> elec_mva_Iso_value_Fall17_;
  // std::vector<float> elec_mva_noIso_value_Fall17_;
  std::vector<bool> elec_pass_conversion_;
  std::vector<unsigned int> elec_nmissinginnerhits_;
};

