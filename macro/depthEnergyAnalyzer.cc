// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
//
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Deepak Kumar
//         Created:  Wed, 10 Apr 2019 11:55:44 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "Math/VectorUtil.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Common/interface/TriggerNames.h"

// pf condidate library


#include "RecoParticleFlow/PFProducer/interface/PFAlgo.h"
#include "RecoParticleFlow/PFProducer/interface/PFMuonAlgo.h"
#include "RecoParticleFlow/PFProducer/interface/PFElectronAlgo.h"
#include "RecoParticleFlow/PFProducer/interface/PFPhotonAlgo.h"
#include "RecoParticleFlow/PFProducer/interface/PFElectronExtraEqual.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackAlgoTools.h"

#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibrationHF.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFSCEnergyCalibration.h"

#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"


//
// class declaration


class depthEnergyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
   public:
      explicit depthEnergyAnalyzer(const edm::ParameterSet&);
      ~depthEnergyAnalyzer();
      //void SetBranches();
    //  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);



   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      TTree* tree_;
      //edm::EDGetTokenT<reco::PFMETCollection> pfMETRawToken_;
      edm::EDGetTokenT<edm::View<reco::PFCluster>> hcaldepth_ENToken_;
      edm::EDGetTokenT<edm::ValueMap<float> > depth_ENToken_;
      float ENGdepth;

      // ----------member data ---------------------------
};

depthEnergyAnalyzer::depthEnergyAnalyzer(const edm::ParameterSet& iConfig):
hcaldepth_ENToken_(consumes<edm::View<reco::PFCluster>>(iConfig.getParameter<edm::InputTag>("hltdepthLabel"))),
depth_ENToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("depthenergyLabel")))

//jetToken_( iConfig.getParameter<InputTag> ( "JetTag" ) )
{
usesResource("TFileService");

}


depthEnergyAnalyzer::~depthEnergyAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
depthEnergyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//  Clear();
   using namespace edm;
   using namespace std;

 //Handle<reco::PFClusterCollection>  pfclusterHcal;
 Handle<edm::View<reco::PFCluster>>  pfclusterHcal;
 if (not iEvent.getByToken(hcaldepth_ENToken_, pfclusterHcal)){
      std::cout << ">>>> pfCluster does not exit !!!\n";
    }


 edm::Handle<edm::ValueMap<float> > depth_handle_;
 //Handle<reco::PFClusterCollection>  pfclusterHcal;
 if (not iEvent.getByToken(depth_ENToken_, depth_handle_)){
    std::cout << ">>>> pfCluster does not exit !!!\n";
  }
  //const std::vector<edm::ValueMap<float>>& depthEn = *depth_handle_;
  //reco::PFClusterCollection cls(*pfclusterHcal.product());
  //std::vector<reco::PFCluster>::const_iterator cl = cls.begin();
  //const std::vector<reco::PFCluster>& pfClusters = *pfclusterHcal;
  for (unsigned int i = 0; i < pfclusterHcal->size(); ++i){
  
  float value = 1.0;
  value = (*depth_handle_)[pfclusterHcal->refAt(i)];
  std::cout << "value:    " << value << std::endl;
  }

  //auto testv=depth_handle_.product()->begin();
  //std::cout << "testv" << testv.values[0] << std::endl;
 //
 // unsigned int nclust =  pfClusters.size();




 


tree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
depthEnergyAnalyzer::beginJob()
{
	edm::Service<TFileService> fs;
	tree_ = fs->make<TTree>("tree_","tree");

	// tree_->Branch("run", &run_, "run/l");
	// tree_->Branch("event", &event_, "event/l");
	// tree_->Branch("lumi", &lumi_, "lumi/l");


}

// ------------ method called once each job just after ending the event loop  ------------
void
depthEnergyAnalyzer::endJob()
{/*
	float dummy = -99999;
	trigName_.clear();
	trigResult_.clear();
	CaloMET=dummy;
*/
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

//define this as a plug-in
DEFINE_FWK_MODULE(depthEnergyAnalyzer);
