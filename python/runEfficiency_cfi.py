import FWCore.ParameterSet.Config as cms

def configureJetMetNtuple(process):
    import FWCore.ParameterSet.Config as cms
    process.load('FWCore.MessageService.MessageLogger_cfi')
    process.TFileService=cms.Service("TFileService",fileName=cms.string('tree_hlt.root'))

    ##-------------------- Define the source  ----------------------------
    #process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

    #process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
    #process.source = cms.Source('PoolSource',
    #    fileNames = cms.untracked.vstring( (

     #       'file:hltoutput_hlt.root'

     #   ) )
    #)


    process.pfRaw     = cms.EDAnalyzer('EfficiencyTreeProducer',
      offlinepfjets   = cms.InputTag('ak4PFJets'),
      depthEnergyTag  = cms.InputTag('HBHEdepthEnergy'),
      hltpfCondidates = cms.InputTag('hltParticleFlow'),
      hltpfjets       = cms.InputTag('hltAK4PFJets'),
      vertices        = cms.InputTag('offlinePrimaryVertices'),#'goodOfflinePrimaryVertices'),
      ptMin           = cms.double(100.0),
      etaMin          = cms.double(-5.0),
      etaMax          = cms.double(5.0),
      rho             = cms.InputTag('fixedGridRhoFastjetAll')#,'','TEST'),
    )
    

    #process.out = cms.OutputModule("PoolOutputModule",
    #                           fileName = cms.untracked.string('output.root'),
    #                           outputCommands = cms.untracked.vstring('keep *' )
    #                           )
    print ('before sequence ') 
    process.JetMetNtupleSequence = cms.Sequence(process.pfRaw)
    print ('after sequence')


