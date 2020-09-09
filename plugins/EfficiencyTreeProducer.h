#ifndef EfficiencyTreeProducer_h
#define EfficiencyTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1F.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class EfficiencyTreeProducer : public edm::EDAnalyzer
{
  public:
    explicit EfficiencyTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~EfficiencyTreeProducer();

  private:
    void initialize();
    //---- configurable parameters --------
    edm::InputTag srcHLTJets_,srcOfflinePFJets_,srcVtx_,srcRho_;
    edm::Service<TFileService> fs_;
    int run_,nvtx_;
    double ptMin_;
    double etaMax_;
    double etaMin_;
    float weight_,rho_;
    bool pass_;
    std::vector<float> *offlinepfPt_,*offlinepfPhi_,*offlinepfEta_,*offlinepfMatchDR_,*hltPt_,*hltEta_,*hltPhi_,*hltMatchDR_, *hltCandEta_;
    std::vector<std::vector<float>> *DepthFractions;
    std::vector<std::vector<float>> *hltDepthFractions;
    TTree *outTree_;

    std::vector<float> *JetPt_,*JetEta_,*JetPhi_,*OfflinehltJetMatchDR_,*offlineGenJetMatchDR_,*hltofflineJetMatchDR_,*hltGenJetMatchDR_,*hltJetPt_,*hltJetEta_,*hltJetPhi_;
};

#endif
