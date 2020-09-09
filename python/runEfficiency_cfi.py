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

    '''
    process.HBHEdepthEnergy    = cms.EDProducer("depthEnergyProducer",
        clustertag             = cms.InputTag('hltParticleFlowClusterHBHE'),
        pfCondtag              = cms.InputTag("hltParticleFlow"),
       )


    '''
    process.pfRaw     = cms.EDAnalyzer('EfficiencyTreeProducer',
      offlinepfjets   = cms.InputTag('ak4PFJets'),#'slimmedJets'),
      JetTag         = cms.InputTag('slimmedJets'),
      #depthEnergyTag  = cms.InputTag('HBHEdepthEnergy'),
      #hltpfCondidates = cms.InputTag('hltParticleFlow'),
      #pfCondidates    = cms.InputTag('particleFlow'),
      hltpfjets       = cms.InputTag('hltAK4PFJets'),
      genJets         = cms.InputTag('ak4GenJets'),
      vertices        = cms.InputTag('offlineSlimmedPrimaryVertices'),#'offlinePrimaryVertices'),#'goodOfflinePrimaryVertices'),
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
    #process.p = cms.Path(process.HBHEdepthEnergy)
    #process.outpath = cms.EndPath(process.out) 

    '''
    process.JetMetNtupleSequence = cms.Sequence(process.HBHEdepthEnergy*
                                                process.pfRaw
    )
    '''

