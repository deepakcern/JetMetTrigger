import FWCore.ParameterSet.Config as cms

import FWCore.PythonUtilities.LumiList as LumiList

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

isData = True
#isData = False
runOnData=isData #data/MC switch

def configureJetMetNtuple(process):
  #  Options parsing  #
  ###################
  
  from FWCore.ParameterSet.VarParsing import VarParsing
  
  options = VarParsing('analysis')
  options.register('applyMETFilters',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply MET filters')
  options.register('runJets',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Store Jet Variables')

  process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
  process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
  process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

  process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
  process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
  process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

  # Electron ID ==========================================================================================

  #from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
  # turn on VID producer, indicate data format  to be
  # DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
  useAOD = False
  
  if useAOD == True :
    dataFormat = DataFormat.AOD
  else :
    dataFormat = DataFormat.MiniAOD
    
    switchOnVIDElectronIdProducer(process, dataFormat)
    
  # define which IDs we want to produce
  #my_id_modules = [#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
  #  'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
    #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
    #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
   # 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff'
   # ]


# define which IDs we want to produce
  my_id_modules = [#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
                 ]



  #add them to the VID producer
  for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

  ### END Electron ID ====================================================================================


    process.hltJetMetNtuple = cms.EDAnalyzer('hltAnalyzer',#'HLTJetMETNtupleProducer',
  #                                           runJets = cms.bool(options.runJets),
 #              isData = cms.bool(True),
#                                             PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                             #MetCollectionTag = cms.InputTag('slimmedMETs'),
                                             #GenMetCollectionTag = cms.InputTag('genMetTrue'),
               #applyMETFilters = cms.bool(options.applyMETFilters),
                                            # BadMuonFilter              = cms.InputTag("BadPFMuonFilter",""),
                                             #BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
                                             MuonCollectionTag = cms.InputTag('hltMuons'),
                                             #ElectronCollectionTag = cms.InputTag('slimmedElectrons'),
                                             #PFJetCollectionTag = cms.InputTag('slimmedJets'),
                                             #CaloJetCollectionTag = cms.InputTag('slimmedCaloJets'),
                                             #GenJetCollectionTag = cms.InputTag('slimmedGenJets'),
                                             #HLTPFJetCollectionTag = cms.InputTag('hltAK4PFJetsCorrected'),
                                             # HLTCaloJetCollectionTag = cms.InputTag('hltAK4CaloJetsCorrected'),
                                             hltdepthLabel = cms.InputTag('hltParticleFlowClusterHCAL'),
               #HLTPFMetCollectionTag = cms.InputTag('hltPFMETProducer'),
                                             HLTCaloMetCollectionTag = cms.InputTag('hltMET'),
               #HLTCaloMetCollectionTag = cms.InputTag('hltCaloMETMuon'),
               
                                             #eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
                                             #eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
                                             #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
                                             #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
                                             #eleMvaNonTrigIdWP80Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
                                             #eleMvaNonTrigIdWP90Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
                                             #eleMvaTrigIdWP80Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
                                             #eleMvaTrigIdWP90Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
                                             #mvaNonTrigValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                                             #mvaNonTrigCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
                                             #mvaTrigValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
                                             #mvaTrigCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
                                             #eleMvaSpring16WPMediumMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
                                             #eleMvaSpring16WPTightMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
                                             #mvaSpring16ValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                                             #mvaSpring16CategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                                             #eleSummer16VetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                                             #eleSummer16LooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                                             #eleSummer16MediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                                            # eleSummer16TightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
                                             hltprocess = cms.InputTag("TriggerResults::HLT"),#changed HLT->TEST
                                             hltprocess_rerun = cms.InputTag("TriggerResults::TEST"),
                                            # triggerObjects = cms.InputTag("slimmedPatTrigger"),
					     caloMetRaw=cms.InputTag("hltMet","","TEST"),   #hltMet,hltPFMETTypeOne,hltPFMETProducer
					     pfMetRaw=cms.InputTag("hltPFMETProducer","","TEST"),
                                             triggerLabel=cms.InputTag("TriggerResults","","TEST"),
                                             #triggerPaths = cms.untracked.vstring('HLT_PFMET200_NotCleaned_v',           
                                             #                                     'HLT_PFMET200_HBHECleaned_v',
                                              #                                    'HLT_PFMET250_HBHECleaned_v',
                                               #                                   'HLT_PFMET300_HBHECleaned_v',
                                               #                                   'HLT_PFMET200_HBHE_BeamHaloCleaned_v',
                                               #                                   'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v', 
                                                #                                  'HLT_PFMET110_PFMHT110_IDTight_v', 
                      # 'HLT_CaloMET80_NotCleaned_v',
                                                  #                                'HLT_IsoMu20_v',
                                                  #                                'HLT_IsoMu24_eta2p1_v',
                                                   #                               'HLT_IsoMu24_v',
                                                   #                               'HLT_CaloJet500_NoJetID_v12',
                                                                                  # 'HLT_CaloMET90_NotCleaned_v',
                                                                                  # 'HLT_CaloMET100_NotCleaned_v',
                                                                                  # 'HLT_CaloMET110_NotCleaned_v',
                                                                                  # 'HLT_CaloMET250_NotCleaned_v',
                                                                                  # 'HLT_CaloMET70_HBHECleaned_v',
                                                                                  # 'HLT_CaloMET80_HBHECleaned_v',
                                                                                  # 'HLT_CaloMET90_HBHECleaned_v',
                                                                                  # 'HLT_CaloMET100_HBHECleaned_v',
                                                                                  # 'HLT_CaloMET250_HBHECleaned_v',
                                                                                  # 'HLT_CaloMET300_HBHECleaned_v',
                                                                                  # 'HLT_CaloMET350_HBHECleaned_v',
                                                                                  # 'HLT_Ele25_WPTight_Gsf_v',
                                                                                  # 'HLT_Ele25_eta2p1_WPTight_Gsf_v',
                                                                                  # 'HLT_Ele27_WPTight_Gsf_v',
                                                                                  # 'HLT_Ele27_eta2p1_WPTight_Gsf_v',
                                                                                  # 'HLT_Ele30_WPTight_Gsf_v',
                                                                                  # 'HLT_Ele30_eta2p1_WPTight_Gsf_v',
                                                      #                            'HLT_Ele32_WPTight_Gsf_v',
                                                                                  # 'HLT_Ele32_eta2p1_WPTight_Gsf_v',
                                                                                  # 'HLT_PFJet40_v',
                                                                                  # 'HLT_PFJet60_v',
                                                                                  # 'HLT_PFJet80_v',
                                                                                  # 'HLT_PFJet140_v',
                                                                                  # 'HLT_PFJet200_v',
                                                                                  # 'HLT_PFJet260_v',
                                                                                  # 'HLT_PFJet320_v',
                                                                                  # 'HLT_PFJet400_v',
                                                       #                           'HLT_PFJet450_v',
                                                        #                          'HLT_PFJet500_v',
                                                       #                           'HLT_PFJet550_v'
                                                                                  # 'HLT_PFJetFwd40_v', 
                                                                                  # 'HLT_PFJetFwd60_v',
                                                                                  # 'HLT_PFJetFwd80_v',
                                                                                  # 'HLT_PFJetFwd140_v',
                                                                                  # 'HLT_PFJetFwd200_v',
                                                                                  # 'HLT_PFJetFwd260_v',
                                                                                  # 'HLT_PFJetFwd320_v',
                                                                                  # 'HLT_PFJetFwd400_v',
                                                                                  # 'HLT_PFJetFwd450_v',
                                                                                  # 'HLT_PFJetFwd500_v'

                                                     #                             )


                      
                                             )
    process.JetMetNtupleSequence = cms.Sequence(process.hltJetMetNtuple) 
    #process.JetMetNtupleSequence = cms.Sequence(process.egmGsfElectronIDSequence*
    #                                            process.BadChargedCandidateFilter*
    #                                            process.BadPFMuonFilter* 
    #                                            process.hltJetMetNtuple)
