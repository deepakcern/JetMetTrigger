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
  std::vector<float> onljetfrac_tmp;
  std::vector<float> offljetfrac_tmp;
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

  OfflJet_DepthFractions  = new std::vector<std::vector<float>>;
  OnlJet_DepthFractions   = new std::vector<std::vector<float>>;

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
  outTree_->Branch("OfflJet_DepthFractions","vector<vector<float>>",&OfflJet_DepthFractions);
  outTree_->Branch("OnlJet_DepthFractions","vector<vector<float>>",&OnlJet_DepthFractions);

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
  delete OfflJet_DepthFractions;
  delete OnlJet_DepthFractions;


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

  OfflJet_DepthFractions->clear();
  OnlJet_DepthFractions->clear();


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
  edm::Handle<edm::View<reco::PFJet>>  hltpfjets_;
  //std::cout << "Entered inside analyser step3 "  << std::endl;
  iEvent.getByToken(Test_HLTPFJets_,hltpfjets_);

  iEvent.getByLabel(srcHLTPFJets_,pfjets);


    for(reco::PFJetCollection::const_iterator ijet = pfjets->begin();ijet != pfjets->end(); ++ijet) { 
      //std::vector<edm::Ptr<reco::Candidate>> constituentPtrs = ijet->getJetConstituents(); 
      std::vector< reco::PFCandidatePtr > iparticles = ijet->getPFConstituents();

      if (fabs(ijet->eta()) < 1.3 ) continue;

      float hcalEnergy = 0.0;
      float hcalEnergy_ = 0.0;
      float hcalEnergyDepth01 = 0.0;
      float hcalEnergyDepth02 = 0.0;
      float hcalEnergyDepth03 = 0.0;
      float hcalEnergyDepth04 = 0.0;
      float hcalEnergyDepth05 = 0.0;
      float hcalEnergyDepth06 = 0.0;
      float hcalEnergyDepth07 = 0.0;


      for ( std::vector<reco::PFCandidatePtr>::const_iterator partBegin = iparticles.begin(), 
           partEnd = iparticles.end(), ipart = partBegin; ipart != partEnd; ++ipart ) {
           hcalEnergy_ = ipart->get()->hcalEnergy();
           auto jet_DepthFractions = ipart->get()->hcalDepthEnergyFractions();
           //std::cout << "DepthFractions" << jet_DepthFractions[0]  << std::endl;

          hcalEnergyDepth01 +=  hcalEnergy_ * jet_DepthFractions[0];
          hcalEnergyDepth02 +=  hcalEnergy_ * jet_DepthFractions[1];
          hcalEnergyDepth03 +=  hcalEnergy_ * jet_DepthFractions[2];
          hcalEnergyDepth04 +=  hcalEnergy_ * jet_DepthFractions[3];
          hcalEnergyDepth05 +=  hcalEnergy_ * jet_DepthFractions[4];
          hcalEnergyDepth06 +=  hcalEnergy_ * jet_DepthFractions[5];
          hcalEnergyDepth07 +=  hcalEnergy_ * jet_DepthFractions[6];
          hcalEnergy += hcalEnergy_;

     }
      onljetfrac_tmp.clear();
      onljetfrac_tmp.push_back(hcalEnergyDepth01/hcalEnergy);
      onljetfrac_tmp.push_back(hcalEnergyDepth02/hcalEnergy);
      onljetfrac_tmp.push_back(hcalEnergyDepth03/hcalEnergy);
      onljetfrac_tmp.push_back(hcalEnergyDepth04/hcalEnergy);
      onljetfrac_tmp.push_back(hcalEnergyDepth05/hcalEnergy);
      onljetfrac_tmp.push_back(hcalEnergyDepth06/hcalEnergy);
      onljetfrac_tmp.push_back(hcalEnergyDepth07/hcalEnergy);
      
      OnlJet_DepthFractions->push_back(onljetfrac_tmp);

      //std::cout << "online jet eta " << ijet->eta() << "    depth1:  "  << hcalEnergyDepth01 << "  depth2:  "  << hcalEnergyDepth02 << "   depth3:   "  << hcalEnergyDepth03 << " depth4:   "  << hcalEnergyDepth04 << "  depth5:  "  << hcalEnergyDepth05  << "   depth6:   "  << hcalEnergyDepth06  << "  depth7:   "  << hcalEnergyDepth07 << std::endl;

      if (ijet->pt() > ptMin_ && (ijet->eta()) >= etaMin_ && (ijet->eta()) < etaMax_) {
        hltPt_ ->push_back(ijet->pt());
        hltEta_->push_back(ijet->eta()); 
        hltPhi_->push_back(ijet->phi()); 

        //std::cout << "jet fraction"  << ijet->hcalDepthEnergyFractions() << std::endl;
 
        pass_ = true;
        // ---- loop over the genjets and find the closest match -----
        float dRmin(1000);
        for(reco::PFJetCollection::const_iterator ipf = offlinepfjets->begin();ipf != offlinepfjets->end(); ++ipf) { 

          //std::vector< reco::PFCandidatePtr > offiparticles = ipf->getPFConstituents();


        //for ( std::vector<reco::PFCandidatePtr>::const_iterator offpartBegin = offiparticles.begin(),          
          //  offpartEnd = offiparticles.end(), offipart = offpartBegin; offipart != offpartEnd; ++offipart ) {
            //auto offjet_DepthFractions = offipart->get()->hcalDepthEnergyFractions();
              //   if (0) {
                //        std::cout << "off DepthFractions"  << offjet_DepthFractions[0]  << std::endl;
                  //      }
          //}

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
    std::vector< reco::PFCandidatePtr > offiparticles = ipf->getPFConstituents();

    if (fabs(ipf->eta()) < 1.3 ) continue;

    float offhcalEnergy = 0.0;
    float offhcalEnergy_ = 0.0;
    float offhcalEnergyDepth01 = 0.0;
    float offhcalEnergyDepth02 = 0.0;
    float offhcalEnergyDepth03 = 0.0;
    float offhcalEnergyDepth04 = 0.0;
    float offhcalEnergyDepth05 = 0.0;
    float offhcalEnergyDepth06 = 0.0;
    float offhcalEnergyDepth07 = 0.0;

    for ( std::vector<reco::PFCandidatePtr>::const_iterator offpartBegin = offiparticles.begin(),
            offpartEnd = offiparticles.end(), offipart = offpartBegin; offipart != offpartEnd; ++offipart ) {

            offhcalEnergy_ = offipart->get()->hcalEnergy();
            auto offjet_DepthFractions = offipart->get()->hcalDepthEnergyFractions();
            //std::cout << "off DepthFractions"  << offjet_DepthFractions[0]  << std::endl;

            offhcalEnergyDepth01 +=  offhcalEnergy_ * offjet_DepthFractions[0];
            offhcalEnergyDepth02 +=  offhcalEnergy_ * offjet_DepthFractions[1];
            offhcalEnergyDepth03 +=  offhcalEnergy_ * offjet_DepthFractions[2];
            offhcalEnergyDepth04 +=  offhcalEnergy_ * offjet_DepthFractions[3];
            offhcalEnergyDepth05 +=  offhcalEnergy_ * offjet_DepthFractions[4];
            offhcalEnergyDepth06 +=  offhcalEnergy_ * offjet_DepthFractions[5];
            offhcalEnergyDepth07 +=  offhcalEnergy_ * offjet_DepthFractions[6];

            offhcalEnergy += offhcalEnergy_;

          }
         offljetfrac_tmp.clear();
         offljetfrac_tmp.push_back(offhcalEnergyDepth01/offhcalEnergy);
         offljetfrac_tmp.push_back(offhcalEnergyDepth02/offhcalEnergy);
         offljetfrac_tmp.push_back(offhcalEnergyDepth03/offhcalEnergy);
         offljetfrac_tmp.push_back(offhcalEnergyDepth04/offhcalEnergy);
         offljetfrac_tmp.push_back(offhcalEnergyDepth05/offhcalEnergy);
         offljetfrac_tmp.push_back(offhcalEnergyDepth06/offhcalEnergy);
         offljetfrac_tmp.push_back(offhcalEnergyDepth07/offhcalEnergy);

    OfflJet_DepthFractions->push_back(offljetfrac_tmp);
 
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
