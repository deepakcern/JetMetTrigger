// -*- C++ -*-
//
// Package:    depthEnergyProducer
// Class:      depthEnergyProducer
//
/**\class depthEnergyProducer depthEnergyProducer.cc ProdTutorial/TrackAndPointsProducer/src/depthEnergyProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Deepak Kumar
//         Created:  Tue Nov 18 17:42:49 CST 2019
// $Id$
//
//

#include <numeric>

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <vector>
#include <iostream>
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
//
#include "DataFormats/Common/interface/ValueMap.h"
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
// class decleration
//

class depthEnergyProducer : public edm::EDProducer {
   public:
      explicit depthEnergyProducer(const edm::ParameterSet&);
      ~depthEnergyProducer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

    //edm::InputTag src_;
    //edm::InputTag culsterTag;
   // edm::EDGetTokenT<View<reco::PFCandidate> > clusterToken_;
    edm::EDGetTokenT<edm::View<reco::PFCandidate>> pfCond_Token_;
    edm::EDGetTokenT<reco::PFClusterCollection> pfClusterToken_;
    bool debug_;

};


depthEnergyProducer::depthEnergyProducer(const edm::ParameterSet& iConfig):

pfCond_Token_(consumes<edm::View<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfCondtag"))),
pfClusterToken_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("clustertag")))

{


  //clusterToken_ = consumes<View<reco::PFCandidate> >(culsterTag;
  produces<edm::ValueMap<float> >("");
 

}
depthEnergyProducer::~depthEnergyProducer()

{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void depthEnergyProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

////////////////
   using namespace edm;
   using namespace reco;
   using namespace std;

   debug_ = false;

   Handle<edm::View<reco::PFCandidate>>  pfcondidates;
   if (not iEvent.getByToken(pfCond_Token_, pfcondidates)){
    std::cout << ">>>> pfCondidates does not exit !!!\n";
   }
   
   Handle<reco::PFClusterCollection>  pfclusterHcal;
   if (not iEvent.getByToken(pfClusterToken_, pfclusterHcal)){
      std::cout << ">>>> pfCluster does not exit !!!\n";
    }
   const std::vector<reco::PFCluster>& pfClusters = *pfclusterHcal;

   unsigned int nclust =  pfClusters.size();

   for (unsigned int icl = 0; icl < nclust; ++icl) {
 
      Ptr<reco::PFCandidate> pcand = pfcondidates->ptrAt( icl );
      reco::PFCandidate fcand = reco::PFCandidate( *pcand );
      std::array<float,7> energyPerDepth;
      std::fill(energyPerDepth.begin(), energyPerDepth.end(), 0.0);

     for (auto & hitRefAndFrac : pfClusters[icl].recHitFractions()) {
        const auto & hit = *hitRefAndFrac.recHitRef();

          if (DetId(hit.detId()).det() == DetId::Hcal) {
             if (hit.depth() == 0) {
                  edm::LogWarning("setHcalDepthInfo") << "Depth zero found";
                  continue;
              }
            if (hit.depth() < 1 || hit.depth() > 7) {
                  throw cms::Exception("CorruptData") << "Bogus depth " << hit.depth() << " at detid " << hit.detId() << "\n";
              }
          energyPerDepth[hit.depth()-1] += hitRefAndFrac.fraction()*hit.energy();

         }

      
       }

    double sum = std::accumulate(energyPerDepth.begin(), energyPerDepth.end(), 0.);
    std::array<float,7> depthFractions;
    if (sum > 0) {
        for (unsigned int i = 0; i < depthFractions.size(); ++i) {
            depthFractions[i] = energyPerDepth[i]/sum;
            //std::cout << "depth energy" << energyPerDepth[i]/sum << std::endl;
        }
    } else {
        std::fill(depthFractions.begin(), depthFractions.end(), 0.f);
    }
    fcand.setHcalDepthEnergyFractions(depthFractions);
    auto DepthFractions = fcand.hcalDepthEnergyFractions();
    //std::cout << "DepthFractions size:   " << DepthFractions.size() << std::endl;
    if(debug_){
     for (int i=0; i < 7; ++i){
      std::cout << "count:   "  << i << std::endl;   
      std::cout << "DepthFractions:   "  << DepthFractions[i]  << std::endl;
      }

    }
    }


}

// ------------ method called once each job just before starting event loop  ------------
void
depthEnergyProducer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
depthEnergyProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(depthEnergyProducer);
