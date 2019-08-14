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

    edm::InputTag src_;
    edm::InputTag culsterTag;
   // edm::EDGetTokenT<View<reco::PFCandidate> > clusterToken_;

    edm::EDGetTokenT<reco::PFClusterCollection> depth_ENToken_;
    std::vector<float> values;
    typedef math::XYZPointD Point;
    typedef std::vector<Point> PointCollection;

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
depthEnergyProducer::depthEnergyProducer(const edm::ParameterSet& iConfig):

depth_ENToken_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("hltdepthLabel")))

{


//  src_  = iConfig.getParameter<edm::InputTag>( "src" );
   //culsterTag = iConfig.getParameter<edm::InputTag>( "culsterTag" );
  // produces<PointCollection>( "innerPoint" ).setBranchAlias( "innerPoints");
  // produces<PointCollection>( "outerPoint" ).setBranchAlias( "outerPoints");
  //clusterToken_ = consumes<View<reco::PFCandidate> >(culsterTag;
  //produces<std::vector<reco::PFClusterCollection>>( "depthEnergy" ).setBranchAlias( "depthEnergy");
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


   //std::auto_ptr<float> tr(new float(0));


   //unique_ptr<std::vector<reco::PFCandidate>> depthEnergy( new std::vector<reco::PFCandidate> );
   //std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>() );
   
   
   //edm::ValueMap<float>::Filler filler(*out);

   values.clear();
   //Handle<View<reco::PFCluster>>  pfclusterHcal;
   Handle<reco::PFClusterCollection>  pfclusterHcal;
   if (not iEvent.getByToken(depth_ENToken_, pfclusterHcal)){
      std::cout << ">>>> pfCluster does not exit !!!\n";
    }
   const std::vector<reco::PFCluster>& pfClusters = *pfclusterHcal;

   unsigned int nclust =  pfClusters.size();//pfclusterHcal->size();//pfClusters.size();
   values.reserve(nclust);

   for (unsigned int icl = 0; icl < nclust; ++icl) {

      // Ptr<reco::PFCluster> pcand = pfclusterHcal->ptrAt( icl );//pfClusters.ptrAt( icl );
      //reco::PFCluster fcand = reco::PFCluster( *pcand );
//       std::cout << "depth lenth:   " << pfClusters[icl]->depth()  << std::endl;
       double energyPerDepth = 0.0 ;
       //for (auto & hitRefAndFrac : pcand->recHitFractions()) {
       //   const auto & hit = *hitRefAndFrac.recHitRef();
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
          energyPerDepth += hitRefAndFrac.fraction()*hit.energy() ;

         }


       }
       //std::cout << "energy :   "  << energyPerDepth  << std::endl;
       //depthEnergy->push_back(energyPerDepth)
       values.push_back(energyPerDepth);
        //fcand.addUserFloat("energyPerDepth",energyPerDepth);
    }
    std::cout << "size of cluster:   " << nclust  << "lenth of values vector"  << values.size() << std::endl;
    auto out = std::make_unique<edm::ValueMap<float>>();
    edm::ValueMap<float>::Filler filler(*out);
    filler.insert(pfclusterHcal, values.begin(), values.end() );
    filler.fill();
    iEvent.put(std::move(out),"");
    //iEvent.put(fcand);// depthEnergy, "depthEnergy" );


   // retrieve the tracks
  //  Handle<TrackCollection> tracks;
  //  iEvent.getByLabel( src_, tracks );
	// std::cout << "I am in the producer" << std::endl;
  //  // create the vectors. Use auto_ptr, as these pointers will automatically
  //  // delete when they go out of scope, a very efficient way to reduce memory leaks.
  //  auto_ptr<PointCollection> innerPoints( new PointCollection );
  //  auto_ptr<PointCollection> outerPoints( new PointCollection );
  //  // and already reserve some space for the new data, to control the size
  //  // of your executible's memory use
  //  const int size = tracks->size();
  //  innerPoints->reserve( size );
  //  outerPoints->reserve( size );
  //  // loop over the tracks:
  //  for( TrackCollection::const_iterator track = tracks->begin();
  //      track != tracks->end(); ++ track ) {
  //    // fill the points in the vectors
  //    innerPoints->push_back( track->innerPosition() );
  //    outerPoints->push_back( track->outerPosition() );
	// std::cout << "innerpoint " << track->innerPosition() <<std::endl;
  //  }
   // and save the vectors
   // iEvent.put( innerPoints, "innerPoint" );
   // iEvent.put( outerPoints, "outerPoint" );


///////////////



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
