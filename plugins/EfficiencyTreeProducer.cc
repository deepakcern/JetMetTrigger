#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TString.h"

#include "HLTrigger/Configuration/plugins/EfficiencyTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h" 
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

using namespace std;
using namespace reco;

namespace pat {
  class Jet;
  typedef std::vector<Jet>              JetCollection;
  typedef edm::Ref<JetCollection>       JetRef;
  typedef edm::RefVector<JetCollection> JetRefVector;
}


// here goes member data:

  edm::InputTag srcOfflinePFJets_; 
  edm::InputTag srcHLTPFJets_;
  edm::InputTag srcHLTCaloJets_;     
  edm::InputTag srcRho_;  
  edm::InputTag srcVtx_; 
  edm::InputTag HBHEdepthEnergy_;
  edm::InputTag hltpfCondidates_;
  edm::InputTag pfCondidates_ ; 
  double ptMin_;  
  double weight_;

  std::vector<float> DepthFractions_;
  std::vector<float> offlinefractions_;
//  edm::InputTag srcMuons_;
   
  edm::EDGetTokenT<reco::PFJetCollection> T_OfflinePFJets_; 
  edm::EDGetTokenT<reco::PFJetCollection> T_HLTPFJets_;
  edm::EDGetTokenT<edm::View<reco::PFJet>> Test_HLTPFJets_; 
  edm::EDGetTokenT<double> T_Rho_;  
  edm::EDGetTokenT<reco::VertexCollection> T_Vtx_;  
  edm::EDGetTokenT<edm::ValueMap<float> > HBHEdepthperEnergy_;
  edm::EDGetTokenT<reco::PFCandidateCollection> hltpfConds_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfConds_;

EfficiencyTreeProducer::EfficiencyTreeProducer(edm::ParameterSet const& cfg) 
{
  srcOfflinePFJets_   = cfg.getParameter<edm::InputTag>             ("offlinepfjets");
  srcHLTPFJets_       = cfg.getParameter<edm::InputTag>             ("hltpfjets");
  HBHEdepthEnergy_    = cfg.getParameter<edm::InputTag>             ("depthEnergyTag");
  hltpfCondidates_    = cfg.getParameter<edm::InputTag>             ("hltpfCondidates");
  pfCondidates_       = cfg.getParameter<edm::InputTag>             ("pfCondidates");
  srcRho_             = cfg.getParameter<edm::InputTag>             ("rho");
  srcVtx_             = cfg.getParameter<edm::InputTag>             ("vertices");
  ptMin_              = cfg.getParameter<double>                    ("ptMin");
  etaMin_             = cfg.getParameter<double>		    ("etaMin");
  etaMax_             = cfg.getParameter<double>		    ("etaMax");

// here goes consume:
  
  T_OfflinePFJets_ = consumes<reco::PFJetCollection> (srcOfflinePFJets_);
  T_HLTPFJets_ = consumes<reco::PFJetCollection> (srcHLTPFJets_);
  Test_HLTPFJets_ = consumes<edm::View<reco::PFJet>> (srcHLTPFJets_);
  HBHEdepthperEnergy_ = consumes<edm::ValueMap<float>> (HBHEdepthEnergy_);
  hltpfConds_         = consumes<reco::PFCandidateCollection>(hltpfCondidates_);
  pfConds_            = consumes<reco::PFCandidateCollection>(pfCondidates_);
  T_Rho_ = consumes<double> (srcRho_);
  T_Vtx_ = consumes<reco::VertexCollection> (srcVtx_);
  
}
//////////////////////////////////////////////////////////////////////////////////////////
void EfficiencyTreeProducer::beginJob() 
{
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo",&run_    ,"run_/I");
  outTree_->Branch("nvtx" ,&nvtx_   ,"nvtx_/I");
  outTree_->Branch("wt"   ,&weight_ ,"wt_/F");
  outTree_->Branch("rho"  ,&rho_    ,"rho_/F");
  outTree_->Branch("pass" ,&pass_   ,"rho_/O");
  offlinepfPt_      = new std::vector<float>;
  offlinepfEta_     = new std::vector<float>;
  offlinepfPhi_     = new std::vector<float>;
  offlinepfMatchDR_ = new std::vector<float>;  
  hltPt_      = new std::vector<float>;
  hltEta_     = new std::vector<float>;
  hltPhi_     = new std::vector<float>;
  hltMatchDR_ = new std::vector<float>;
  hltCandEta_ = new std::vector<float>; 

  depthEnergyFraction = new std::vector<std::vector<float>>;
  depthEnergyFraction_offline = new std::vector<std::vector<float>>;
  outTree_->Branch("offlinepfPt"     ,"vector<float>",&offlinepfPt_);
  outTree_->Branch("offlinepfEta"    ,"vector<float>",&offlinepfEta_); 
  outTree_->Branch("offlinepfPhi"    ,"vector<float>",&offlinepfPhi_);
  outTree_->Branch("offlinepfMatchDR","vector<float>",&offlinepfMatchDR_);
  outTree_->Branch("hltPt"     ,"vector<float>",&hltPt_);
  outTree_->Branch("hltEta"    ,"vector<float>",&hltEta_);
  outTree_->Branch("hltPhi"    ,"vector<float>",&hltPhi_);
  outTree_->Branch("hltMatchDR","vector<float>",&hltMatchDR_);
  outTree_->Branch("hltCandEta", "vector<float>",&hltCandEta_);
  outTree_->Branch("depthFractions","vector<vector<float>>",&depthEnergyFraction);
  outTree_->Branch("depthFractions_offline","vector<vector<float>>",&depthEnergyFraction_offline);

  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void EfficiencyTreeProducer::endJob() 
{  
  delete offlinepfPt_;
  delete offlinepfEta_;
  delete offlinepfPhi_;
  delete offlinepfMatchDR_;
  delete hltPt_;
  delete hltEta_;
  delete hltPhi_;
  delete hltMatchDR_;
  delete hltCandEta_;
  delete depthEnergyFraction;
  delete depthEnergyFraction_offline;
}
//////////////////////////////////////////////////////////////////////////////////////////
void EfficiencyTreeProducer::initialize()
{
  pass_ = false;
  run_  = -999;
  nvtx_ = -999;
  rho_  = -999;
  offlinepfPt_ ->clear();
  offlinepfEta_->clear();
  offlinepfPhi_->clear();
  offlinepfMatchDR_->clear();
  hltPt_ ->clear();
  hltEta_->clear();
  hltPhi_->clear();
  hltMatchDR_->clear();
  hltCandEta_->clear();
  depthEnergyFraction->clear();
  depthEnergyFraction_offline->clear();
}
//////////////////////////////////////////////////////////////////////////////////////////
void EfficiencyTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);
  
  edm::Handle<reco::PFJetCollection> offlinepfjets;
  iEvent.getByLabel(srcOfflinePFJets_,offlinepfjets);
  
  edm::Handle<reco::PFJetCollection> pfjets;
   
  edm::Handle<reco::VertexCollection> vtxs;
  iEvent.getByLabel(srcVtx_,vtxs);

  //std::cout << "Entered inside analyser step1"  << std::endl;

  //edm::Handle<edm::ValueMap<float> > depth_handle_;
  //std::cout << "Entered inside analyser step2 "  << std::endl;

  //if (not iEvent.getByToken(HBHEdepthperEnergy_, depth_handle_)){
  //   std::cout << ">>> depthEnrgy does not exit !!!\n";
  //   return;
  // }

  edm::Handle<reco::PFCandidateCollection> cands_test;
  iEvent.getByToken(pfConds_, cands_test);
  for (unsigned int ij = 0, cn = cands_test->size(); ij < cn; ++ij){
     const reco::PFCandidate &cand_test = (*cands_test)[ij];
     auto fraction = cand_test.hcalDepthEnergyFractions();
     if (fabs(cand_test.eta()) < 1.3) continue;
     if (cand_test.pt() < 1) continue;
     //std::cout << "offl condidates pdgId  "  << fabs(cand_test.pdgId()) << std::endl;
     offlinefractions_.clear();
     for (unsigned int i=0; i < fraction.size(); i++){
        offlinefractions_.push_back(fraction[i]);
        //std::cout << "offl frac"  << fraction[i] << std::endl;
   }
  depthEnergyFraction_offline->push_back(offlinefractions_);
  }


  edm::Handle<reco::PFCandidateCollection> cands;
  iEvent.getByToken(hltpfConds_, cands);
  for (unsigned int ic = 0, nc = cands->size(); ic < nc; ++ic) {
    const reco::PFCandidate &cand = (*cands)[ic];
    auto DepthFractions = cand.hcalDepthEnergyFractions();
    //std::cout << "fraction"  << DepthFractions[0]  << std::endl;
    //std::cout << "cand eta"  << cand.eta() << endl;  
    //std::cout << "condidates pdgId  "  << fabs(cand.pdgId() << std::endl;
    if (fabs(cand.eta()) < 1.3 ) continue;
    if (cand.pt() < 1) continue;
    //std::cout << "onl condidates pdgId  "  << fabs(cand.pdgId()) << std::endl;
    hltCandEta_->push_back(cand.eta()); 
    //std::cout << "condidates pdgId  "  << fabs(cand.pdgId()) << "condidate pt " << cand.pt() << std::endl;
    DepthFractions_.clear(); 
    for (unsigned int i=0; i < DepthFractions.size(); i++){
        DepthFractions_.push_back(DepthFractions[i]);
        //std::cout << "onl frac "  << DepthFractions[i] << std::endl;
         }
    depthEnergyFraction->push_back(DepthFractions_);       
        
    
   
  }

  //std::cout << "Entered inside analyser step2 "  << std::endl;
  edm::Handle<edm::View<reco::PFJet>>  hltpfjets;
  //std::cout << "Entered inside analyser step3 "  << std::endl;
  iEvent.getByToken(Test_HLTPFJets_,hltpfjets);
  //std::cout << "Entered inside analyser"  << std::endl;
  //std::cout << "size of hltpfjets " << hltpfjets->size() << std::endl;
  //for (unsigned int i = 0; i < hltpfjets->size(); ++i){
  //    float value = 1.0;
  //    value = (*depth_handle_)[hltpfjets->refAt(i)];
  //    std::cout << "value:  "  << value << std::endl;
  // }

    iEvent.getByLabel(srcHLTPFJets_,pfjets);
    for(reco::PFJetCollection::const_iterator ijet = pfjets->begin();ijet != pfjets->end(); ++ijet) {  
      if (ijet->pt() > ptMin_ && (ijet->eta()) >= etaMin_ && (ijet->eta()) < etaMax_) {
        hltPt_ ->push_back(ijet->pt());
        hltEta_->push_back(ijet->eta()); 
        hltPhi_->push_back(ijet->phi()); 

        //std::cout << "jet fraction"  << ijet->hcalDepthEnergyFractions() << std::endl;
 
        pass_ = true;
        // ---- loop over the genjets and find the closest match -----
        float dRmin(1000);
        for(reco::PFJetCollection::const_iterator ipf = offlinepfjets->begin();ipf != offlinepfjets->end(); ++ipf) { 
          float dR = deltaR(ipf->eta(),ipf->phi(),ijet->eta(),ijet->phi());
          if (dR < dRmin) {
            dRmin = dR;
          }
        }
        hltMatchDR_->push_back(dRmin);
      }
    }
   
  


  // ---- loop over the offlinepfjets and check the matching radius with any of the firing jets ------
   for(reco::PFJetCollection::const_iterator ipf = offlinepfjets->begin();ipf != offlinepfjets->end(); ++ipf) {  
    if (ipf->pt() > 20 && (ipf->eta()) >= etaMin_ && (ipf->eta()) < etaMax_) {
      offlinepfPt_ ->push_back(ipf->pt());
      offlinepfEta_->push_back(ipf->eta()); 
      offlinepfPhi_->push_back(ipf->phi()); 
      float dRmin(1000);
      for(unsigned k=0;k<hltEta_->size();k++) {
        float dR = deltaR(ipf->eta(),ipf->phi(),(*hltEta_)[k],(*hltPhi_)[k]);
        if (dR < dRmin) {
          dRmin = dR;
        }
      }
      offlinepfMatchDR_->push_back(dRmin);
    }
   }
  
  rho_  = *rho;
  nvtx_ = vtxs->size();
  run_  = iEvent.id().run();
  //weight_ = eventInfo->weight();
  weight_ = 1;

  outTree_->Fill();
  
   
}
//////////////////////////////////////////////////////////////////////////////////////////
EfficiencyTreeProducer::~EfficiencyTreeProducer() 
{
}

DEFINE_FWK_MODULE(EfficiencyTreeProducer);
