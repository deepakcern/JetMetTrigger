import FWCore.ParameterSet.Config as cms

process = cms.Process("ANA")

process.load('FWCore/MessageService/MessageLogger_cfi')
'''
process.MessageLogger.categories.append('TriggerSummaryAnalyzerAOD')
process.MessageLogger.categories.append('TriggerSummaryAnalyzerRAW')
process.MessageLogger.categories.append('HLTEventAnalyzerAOD')
process.MessageLogger.categories.append('HLTEventAnalyzerRAW')
process.MessageLogger.categories.append('L1GtTrigReport')
process.MessageLogger.categories.append('L1TGlobalSummary')
process.MessageLogger.categories.append('HLTrigReport')
process.MessageLogger.categories.append('HLTSummaryFilter')
process.MessageLogger.categories.append('HLTConfigProvider')
process.MessageLogger.categories.append('HLTPrescaleProvider')
process.MessageLogger.categories.append('HLTConfigData')
'''
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('Configuration.StandardSequences.CondDBESSource_cff')
from Configuration.AlCa.GlobalTag import GlobalTag as customiseGlobalTag
# process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:run2_hlt_GRun')
process.GlobalTag = customiseGlobalTag(None, globaltag = 'auto:run2_hlt_GRun')

# process.Timing = cms.Service("Timing")

# process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#     ignoreTotal = cms.untracked.int32(-1) ## default is one
# )

# process.Tracer = cms.Service("Tracer")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2)
)


#import FWCore.Utilities.FileUtils as FileUtils
#files=FileUtils.loadListFromFile('/afs/cern.ch/work/d/dekumar/public/JETMET_Work/work3/CMSSW_10_1_11_patch1/src/HLTrigger/Configuration/Filelists/HLT_withoutFilterStudy_v4_v21.txt')#withoutFilter1.txt')

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/d/dekumar/public/JETMET_Work/depthStudeis/CMSSW_10_6_1_patch1/src/HLTrigger/Configuration/test/hltoutput.root')
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False) ## default is false
)

'''
import HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi
process.tsaAOD = HLTrigger.HLTcore.triggerSummaryAnalyzerAOD_cfi.triggerSummaryAnalyzerAOD.clone()
import HLTrigger.HLTcore.triggerSummaryAnalyzerRAW_cfi
process.tsaRAW = HLTrigger.HLTcore.triggerSummaryAnalyzerRAW_cfi.triggerSummaryAnalyzerRAW.clone()
process.tsa = cms.Path(process.tsaAOD)#+process.tsaRAW)


import HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi
process.hltAOD = HLTrigger.HLTcore.hltEventAnalyzerAOD_cfi.hltEventAnalyzerAOD.clone()
process.hltAOD.processName = cms.string("TEST")
process.hltAOD.triggerResults = cms.InputTag("TriggerResults","","TEST")
process.hltAOD.triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","TEST")
#process.hltAOD.pfmet   = cms.InputTag("hltPFMETTypeOne","TEST")
#process.hltAOD.calomet   = cms.InputTag("hltMet","TEST")

import HLTrigger.HLTcore.hltEventAnalyzerRAW_cfi
process.hltRAW = HLTrigger.HLTcore.hltEventAnalyzerRAW_cfi.hltEventAnalyzerRAW.clone()
process.hlt = cms.Path(process.hltAOD)#+process.hltRAW)


import HLTrigger.HLTanalyzers.hltTrigReport_cfi
process.hltReport = HLTrigger.HLTanalyzers.hltTrigReport_cfi.hltTrigReport.clone()
process.hltReport.HLTriggerResults = cms.InputTag("TriggerResults","","TEST")
'''
#process.aom = cms.OutputModule("AsciiOutputModule")
#process.eca = cms.EDAnalyzer("EventContentAnalyzer")
process.depthEnergy = cms.EDProducer("depthEnergyProducer",
	caloMetRaw=cms.InputTag("hltMet","","TEST"),   #hltMet,hltPFMETTypeOne,hltPFMETProducer
	pfMetRaw=cms.InputTag("hltPFMETProducer","","TEST"),
        hltdepthLabel = cms.InputTag('hltParticleFlowClusterHBHE'),
	MuonCollectionTag = cms.InputTag('hltMuons'),
	triggerLabel=cms.InputTag("TriggerResults","","TEST")

 )


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('output.root'),
                               outputCommands = cms.untracked.vstring('keep *' )
                               )


'''
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("test_ntuples.root")
                                   )
'''


process.p = cms.Path(process.depthEnergy)

process.outpath = cms.EndPath(process.out)
#process.final = cms.EndPath(process.hltReport+process.aom)#+process.eca)


