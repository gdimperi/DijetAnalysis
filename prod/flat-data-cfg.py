import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'

process.TFileService=cms.Service("TFileService",fileName=cms.string('dijetTree_data.root'))
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring('/store/cmst3/group/das2014/EXODijetsLE/test_data.root')
)

#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load('cmsdas2014.exo_dijets_exercise.PAT_ca8jets_simple_cff')

process.patJetCorrFactorsCA8CHS.levels.append('L2L3Residual')
process.patJetCorrFactorsCA8CHSpruned.levels.append('L2L3Residual')

##-------------------- User analyzer  --------------------------------
process.dijets     = cms.EDAnalyzer('DijetTreeProducer',
  jets             = cms.InputTag('patJetsCA8CHSwithNsub'),
  jetsPruned       = cms.InputTag('patJetsCA8CHSpruned'),
  met              = cms.InputTag('pfMet'),
  vtx              = cms.InputTag('goodOfflinePrimaryVertices'),
  mjjMin           = cms.double(700),
  ptMin            = cms.double(40),
  dEtaMax          = cms.double(1.3),
  ## trigger ###################################
  triggerAlias     = cms.vstring('Fat','PFHT650','PFNoPUHT650','HT750','HT550'),
  triggerSelection = cms.vstring(
    'HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v*',
    'HLT_PFHT650_v*',
    'HLT_PFNoPUHT650_v*',
    'HLT_HT750_v*',  
    'HLT_HT550_v*'
  ),
  triggerConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','HLT'),
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMask         = cms.bool(False),
    l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  )
)
##----------------- hlt filter -------------------------------------
process.hltFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag = cms.InputTag('TriggerResults','','HLT'),
    HLTPaths          = cms.vstring(
      'HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v*',
      'HLT_PFHT650_v*',
      'HLT_PFNoPUHT650_v*',
      'HLT_HT750_v*',  
      'HLT_HT550_v*'
    ),  
    eventSetupPathsKey = cms.string(''),
    andOr              = cms.bool(True), #---- True = OR, False = AND between the HLT paths
    throw              = cms.bool(False)
)

process.p = cms.Path(process.hltFilter * process.ca8Jets * process.dijets)
