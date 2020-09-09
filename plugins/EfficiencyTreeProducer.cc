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

//#include "HLTrigger/Configuration/plugins/EfficiencyTreeProducer.h"
#include "JMETriggerAnalysis/Common/plugins/EfficiencyTreeProducer.h"
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
  edm::InputTag srcGenJets_;
  edm::InputTag inputTagJets_;
  // edm::InputTag hltpfCondidates_;
  // edm::InputTag pfCondidates_ ;
  double ptMin_;
  double weight_;

  std::vector<float> DepthFractions_;
  std::vector<float> hltDepthFractions_;
  std::vector<float> offlinefractions_;
  std::vector<float> hltjetEnergyfrac;
  std::vector<float> jetEnergyfrac;
  // std::vector<float> offljetfrac_tmp;
//  edm::InputTag srcMuons_;

  edm::EDGetTokenT<pat::JetCollection> T_OfflinePFJets_;
  edm::EDGetTokenT<reco::PFJetCollection> T_HLTPFJets_;
  edm::EDGetTokenT<reco::GenJetCollection> T_GenJets_;
  edm::EDGetTokenT<edm::View<pat::Jet>> jetToken_;
  // edm::EDGetTokenT<edm::View<reco::PFJet>> Test_HLTPFJets_;
  edm::EDGetTokenT<double> T_Rho_;
  edm::EDGetTokenT<reco::VertexCollection> T_Vtx_;
  // edm::EDGetTokenT<edm::ValueMap<float> > HBHEdepthperEnergy_;
  // edm::EDGetTokenT<reco::PFCandidateCollection> hltpfConds_;
  // edm::EDGetTokenT<reco::PFCandidateCollection> pfConds_;

EfficiencyTreeProducer::EfficiencyTreeProducer(edm::ParameterSet const& cfg)
{
  srcOfflinePFJets_   = cfg.getParameter<edm::InputTag>             ("offlinepfjets");
  srcHLTPFJets_       = cfg.getParameter<edm::InputTag>             ("hltpfjets");
  srcGenJets_         = cfg.getParameter<edm::InputTag>             ("genJets");
  inputTagJets_       = cfg.getParameter<edm::InputTag>             ( "JetTag" );

  // HBHEdepthEnergy_    = cfg.getParameter<edm::InputTag>             ("depthEnergyTag");
  // hltpfCondidates_    = cfg.getParameter<edm::InputTag>             ("hltpfCondidates");
  // pfCondidates_       = cfg.getParameter<edm::InputTag>             ("pfCondidates");
  srcRho_             = cfg.getParameter<edm::InputTag>             ("rho");
  srcVtx_             = cfg.getParameter<edm::InputTag>             ("vertices");
  ptMin_              = cfg.getParameter<double>                    ("ptMin");
  etaMin_             = cfg.getParameter<double>		                ("etaMin");
  etaMax_             = cfg.getParameter<double>		                ("etaMax");

// here goes consume:

  T_OfflinePFJets_ = consumes<pat::JetCollection> (srcOfflinePFJets_);
  T_HLTPFJets_     = consumes<reco::PFJetCollection> (srcHLTPFJets_);
  T_GenJets_       = consumes<reco::GenJetCollection> (srcGenJets_);
  jetToken_        = consumes<edm::View<pat::Jet> >(inputTagJets_);
  // jetsToken_       = consumes<edm::View<reco::Jet> >(cfg.getParameter<edm::InputTag>("jetSource"));
  // Test_HLTPFJets_ = consumes<edm::View<reco::PFJet>> (srcHLTPFJets_);
  // HBHEdepthperEnergy_ = consumes<edm::ValueMap<float>> (HBHEdepthEnergy_);
  // hltpfConds_         = consumes<reco::PFCandidateCollection>(hltpfCondidates_);
  // pfConds_            = consumes<reco::PFCandidateCollection>(pfCondidates_);
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


  JetPt_                = new std::vector<float>;
  JetEta_               = new std::vector<float>;
  JetPhi_               = new std::vector<float>;
  OfflinehltJetMatchDR_ = new std::vector<float>;
  offlineGenJetMatchDR_ = new std::vector<float>;
  hltofflineJetMatchDR_ = new std::vector<float>;
  hltGenJetMatchDR_     = new std::vector<float>;
  hltJetPt_             = new std::vector<float>;
  hltJetEta_            = new std::vector<float>;
  hltJetPhi_            = new std::vector<float>;

  DepthFractions        = new std::vector<std::vector<float>>;
  hltDepthFractions     = new std::vector<std::vector<float>>;
  // hltMatchDR_           = new std::vector<float>;
  // hltCandEta_       = new std::vector<float>;

  // depthEnergyFraction         = new std::vector<std::vector<float>>;
  // depthEnergyFraction_offline = new std::vector<std::vector<float>>;

  outTree_->Branch("JetPt"     ,"vector<float>",&JetPt_);
  outTree_->Branch("JetEta"    ,"vector<float>",&JetEta_);
  outTree_->Branch("JetPhi"    ,"vector<float>",&JetPhi_);
  outTree_->Branch("OfflinehltJetMatchDR","vector<float>",&OfflinehltJetMatchDR_);
  outTree_->Branch("hltofflineJetMatchDR","vector<float>",&hltofflineJetMatchDR_);
  outTree_->Branch("hltJetPt"     ,"vector<float>",&hltJetPt_);
  outTree_->Branch("hltJetEta"    ,"vector<float>",&hltJetEta_);
  outTree_->Branch("hltJetPhi"    ,"vector<float>",&hltJetPhi_);
  // outTree_->Branch("hltMatchDR","vector<float>",&hltMatchDR_);
  // outTree_->Branch("hltCandEta", "vector<float>",&hltCandEta_);
  // outTree_->Branch("depthFractions","vector<vector<float>>",&depthEnergyFraction);
  // outTree_->Branch("depthFractions_offline","vector<vector<float>>",&depthEnergyFraction_offline);
  outTree_->Branch("DepthFractions","vector<vector<float>>",&DepthFractions);
  outTree_->Branch("hltDepthFractions","vector<vector<float>>",&hltDepthFractions);

  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void EfficiencyTreeProducer::endJob()
{
  delete JetPt_;
  delete JetEta_;
  delete JetPhi_;
  delete hltJetPt_;
  delete hltJetEta_;
  delete hltJetPhi_;
  // delete hltMatchDR_;
  delete offlineGenJetMatchDR_;
  delete hltGenJetMatchDR_;
  delete OfflinehltJetMatchDR_;
  delete hltofflineJetMatchDR_;
  // delete hltCandEta_;
  // delete depthEnergyFraction;
  // delete hltdepthEnergyFraction;
  delete DepthFractions;
  delete hltDepthFractions;


}
//////////////////////////////////////////////////////////////////////////////////////////
void EfficiencyTreeProducer::initialize()
{
  pass_ = false;
  run_  = -999;
  nvtx_ = -999;
  rho_  = -999;
  JetPt_ ->clear();
  JetEta_->clear();
  JetPhi_->clear();
  OfflinehltJetMatchDR_->clear();
  hltofflineJetMatchDR_->clear();
  hltJetPt_ ->clear();
  hltJetEta_->clear();
  hltJetPhi_->clear();
  DepthFractions->clear();
  hltDepthFractions->clear();
  hltGenJetMatchDR_->clear();
  offlineGenJetMatchDR_->clear();
  // hltMatchDR_->clear();
  // hltCandEta_->clear();
  // depthEnergyFraction->clear();
  // depthEnergyFraction_offline->clear();




}
//////////////////////////////////////////////////////////////////////////////////////////
void EfficiencyTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
{
  initialize();

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<pat::JetCollection> offlinepfjets;
  iEvent.getByLabel(srcOfflinePFJets_,offlinepfjets);

  edm::Handle<reco::PFJetCollection> hltpfjets;
  iEvent.getByLabel(srcHLTPFJets_,hltpfjets);

  edm::Handle<reco::GenJetCollection> genjets;
  iEvent.getByLabel(srcGenJets_,genjets);

  edm::Handle<edm::View<pat::Jet> > testjets;
  iEvent.getByToken(jetToken_, testjets);

  for (edm::View<pat::Jet>::const_iterator itJet = testjets->begin(); itJet != testjets->end(); itJet++) {
    std::cout << itJet->getPFConstituents().size() << std::endl;

 }


  edm::Handle<reco::VertexCollection> vtxs;
  iEvent.getByLabel(srcVtx_,vtxs);


    for(reco::PFJetCollection::const_iterator ijet = hltpfjets->begin();ijet != hltpfjets->end(); ++ijet) {
      std::vector< reco::PFCandidatePtr > iparticles = ijet->getPFConstituents();

      if (fabs(ijet->eta()) < 1.3 or ijet->pt() < 20.0) continue;

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

          hcalEnergyDepth01 +=  hcalEnergy_ * jet_DepthFractions[0];
          hcalEnergyDepth02 +=  hcalEnergy_ * jet_DepthFractions[1];
          hcalEnergyDepth03 +=  hcalEnergy_ * jet_DepthFractions[2];
          hcalEnergyDepth04 +=  hcalEnergy_ * jet_DepthFractions[3];
          hcalEnergyDepth05 +=  hcalEnergy_ * jet_DepthFractions[4];
          hcalEnergyDepth06 +=  hcalEnergy_ * jet_DepthFractions[5];
          hcalEnergyDepth07 +=  hcalEnergy_ * jet_DepthFractions[6];
          hcalEnergy += hcalEnergy_;

     }
      hltjetEnergyfrac.clear();
      hltjetEnergyfrac.push_back(hcalEnergyDepth01/hcalEnergy);
      hltjetEnergyfrac.push_back(hcalEnergyDepth02/hcalEnergy);
      hltjetEnergyfrac.push_back(hcalEnergyDepth03/hcalEnergy);
      hltjetEnergyfrac.push_back(hcalEnergyDepth04/hcalEnergy);
      hltjetEnergyfrac.push_back(hcalEnergyDepth05/hcalEnergy);
      hltjetEnergyfrac.push_back(hcalEnergyDepth06/hcalEnergy);
      hltjetEnergyfrac.push_back(hcalEnergyDepth07/hcalEnergy);

      hltDepthFractions->push_back(hltjetEnergyfrac);


      // if (1){//ijet->pt() > ptMin_ && (ijet->eta()) >= etaMin_ && (ijet->eta()) < etaMax_)
        //std::cout << "pt" << ijet->pt() << "eta"  << ijet->eta() << "phi"  << ijet->phi() << std::endl;
        hltJetPt_->push_back(ijet->pt());
        hltJetEta_->push_back(ijet->eta());
        hltJetPhi_->push_back(ijet->phi());
      //
      //
        pass_ = true;
        // ---- loop over the genjets and find the closest match -----
        float dRmin(1000);
        for(pat::JetCollection::const_iterator ipf = offlinepfjets->begin();ipf != offlinepfjets->end(); ++ipf) {
          if (fabs(ipf->eta()) < 1.3 or ipf->pt() < 20.0) continue;

          float dR = deltaR(ipf->eta(),ipf->phi(),ijet->eta(),ijet->phi());
          if (dR < dRmin) {
            dRmin = dR;
          }
        }//end of offline jets

      hltofflineJetMatchDR_->push_back(dRmin);
      // }
      // loop over genjets
      float dRminjen(1000);
      for (reco::GenJetCollection::const_iterator igenjet = genjets->begin(); igenjet != genjets->end(); ++igenjet ){
        if (fabs(igenjet->eta()) < 1.3 or igenjet->pt() < 20.0) continue;

        float dR = deltaR(igenjet->eta(),igenjet->phi(),ijet->eta(),ijet->phi());
        if (dR < dRmin){
          dRminjen =dR;
        }
        hltGenJetMatchDR_->push_back(dRminjen);
      }//end of gen jets loop


    }// end of loop over reco/hlt jets




  // ---- loop over the offlinepfjets and check the matching radius with any of the firing jets ------

   for(pat::JetCollection::const_iterator ipf = offlinepfjets->begin();ipf != offlinepfjets->end(); ++ipf) {
    // std::vector< reco::PFCandidatePtr > offiparticles = ipf->getPFConstituents();

    if (fabs(ipf->eta()) < 1.3 or ipf->pt() < 20.0) continue;
    // std::cout << ipf->userFloat("hcalFracDepth01") << std::endl;

    float hcalEnergy = 0.0;
    float hcalEnergy_ = 0.0;
    float hcalEnergyDepth01 = 0.0;
    float hcalEnergyDepth02 = 0.0;
    float hcalEnergyDepth03 = 0.0;
    float hcalEnergyDepth04 = 0.0;
    float hcalEnergyDepth05 = 0.0;
    float hcalEnergyDepth06 = 0.0;
    float hcalEnergyDepth07 = 0.0;
    //
    // for ( std::vector<reco::PFCandidatePtr>::const_iterator partBegin = offiparticles.begin(),
    //         partEnd = offiparticles.end(), ipart = partBegin; ipart != partEnd; ++ipart ) {
    //
    //         hcalEnergy_ = ipart->get()->hcalEnergy();
    //         auto jet_DepthFractions = ipart->get()->hcalDepthEnergyFractions();
    //         //std::cout << "off DepthFractions"  << offjet_DepthFractions[0]  << std::endl;
    //
    //         hcalEnergyDepth01 +=  hcalEnergy_ * jet_DepthFractions[0];
    //         hcalEnergyDepth02 +=  hcalEnergy_ * jet_DepthFractions[1];
    //         hcalEnergyDepth03 +=  hcalEnergy_ * jet_DepthFractions[2];
    //         hcalEnergyDepth04 +=  hcalEnergy_ * jet_DepthFractions[3];
    //         hcalEnergyDepth05 +=  hcalEnergy_ * jet_DepthFractions[4];
    //         hcalEnergyDepth06 +=  hcalEnergy_ * jet_DepthFractions[5];
    //         hcalEnergyDepth07 +=  hcalEnergy_ * jet_DepthFractions[6];
    //
    //         hcalEnergy += hcalEnergy_;
    //
    //       }
    //      jetEnergyfrac.clear();
    //      jetEnergyfrac.push_back(hcalEnergyDepth01/hcalEnergy);
    //      jetEnergyfrac.push_back(hcalEnergyDepth02/hcalEnergy);
    //      jetEnergyfrac.push_back(hcalEnergyDepth03/hcalEnergy);
    //      jetEnergyfrac.push_back(hcalEnergyDepth04/hcalEnergy);
    //      jetEnergyfrac.push_back(hcalEnergyDepth05/hcalEnergy);
    //      jetEnergyfrac.push_back(hcalEnergyDepth06/hcalEnergy);
    //      jetEnergyfrac.push_back(hcalEnergyDepth07/hcalEnergy);
    //
    // DepthFractions->push_back(jetEnergyfrac);
    //
    // if (1){//ipf->pt() > 20 && (ipf->eta()) >= etaMin_ && (ipf->eta()) < etaMax_) {
    //   offlinepfPt_ ->push_back(ipf->pt());
    //   offlinepfEta_->push_back(ipf->eta());
    //   offlinepfPhi_->push_back(ipf->phi());
    //   float dRmin(1000);
    //   for(unsigned k=0;k<hltEta_->size();k++) {
    //     float dR = deltaR(ipf->eta(),ipf->phi(),(*hltEta_)[k],(*hltPhi_)[k]);
    //     if (dR < dRmin) {
    //       dRmin = dR;
    //     }
    //   }
    //   offlinepfMatchDR_->push_back(dRmin);
    // }

  //loop over gen jets
  // float dRmin(1000);
  // for (reco::GenJetCollection::const_iterator igenjet = genjets->begin(); igenjet != genjets->end(); ++igenjet ){
  //   if (fabs(igenjet->eta()) < 1.3 or igenjet->pt() < 20.0) continue;
  //
  //   float dR = deltaR(igenjet->eta(),igenjet->phi(),ipf->eta(),ipf->phi());
  //   if (dR < dRmin){
  //     dRmin =dR;
  //   }
  //   offlineGenJetMatchDR_->push_back(dRmin);
  // }

}// end of jet loop


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
