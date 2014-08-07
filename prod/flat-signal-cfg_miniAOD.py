import FWCore.ParameterSet.Config as cms 

process = cms.Process('myprocess')

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'START53_V27::All'
#process.GlobalTag.globaltag = 'START70_V6::All'
#process.GlobalTag.globaltag = 'POSTLS170_V5::All'
#process.GlobalTag.globaltag = 'POSTLS170_V7::All'
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'


process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService=cms.Service("TFileService",fileName=cms.string('dijetTree_signal.root'))

##-------------------- Define the source  ----------------------------

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:RSGravToJJ_kMpl01_M-1000_test.root')
    fileNames = cms.untracked.vstring('file:/cmshome/santanas/CMS/data/Spring14miniaod__RSGravToJJ_kMpl01_M-1000_Tune4C_13TeV-pythia8__MINIAODSIM__PU20bx25_POSTLS170_V5-v1__00000__6AACD832-3707-E411-A167-001E672489D5.root')
    #fileNames = cms.untracked.vstring('file:/cmshome/santanas/CMS/data/Spring14drAODSIM__RSGravToJJ_kMpl01_M-1000_Tune4C_13TeV-pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__00000__0622C950-58E4-E311-A595-0025904B130A.root')
    
)

##-------------------- Output  --------------------------------

# process.OUT = cms.OutputModule("PoolOutputModule",
#     fileName = cms.untracked.string('test.root'),
#     outputCommands = cms.untracked.vstring(['drop *','keep patJets_patJetsAK4PFCHS_*_*']) 
# )
# process.endpath= cms.EndPath(process.OUT)

##-------------------- Create collections  --------------------------------

# Select candidates that would pass CHS requirements
process.chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

##---- AK4 ----
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets

process.ak4PFJetsCHS = ak4PFJets.clone(
    src = 'chs',
    #jetPtMin = cms.double(10.0),   
    doAreaFastjet = cms.bool(True), 
    #doRhoFastjet = cms.bool(True), 
    rParam = cms.double(0.4),
    jetAlgorithm = cms.string("AntiKt"),
    )
#process.ak4PFJets = ak4PFJets.clone(src = 'packedPFCandidates') 
process.ak4GenJetsMy = ak4GenJets.clone(src = 'packedGenParticles')

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(
   process,
   postfix   = "",
   labelName = 'AK4PFCHS',
   jetSource = cms.InputTag('ak4PFJetsCHS'),
   trackSource = cms.InputTag('unpackedTracksAndVertices'), 
   pvSource = cms.InputTag('unpackedTracksAndVertices'), 
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
   btagDiscriminators = [      'combinedSecondaryVertexBJetTags'     ]
   ,algo= 'AK', rParam = 0.4
   )

#adjust MC matching
process.patJetGenJetMatchAK4PFCHS.matched = "ak4GenJetsMy"
process.patJetPartonMatchAK4PFCHS.matched = "prunedGenParticles"
process.patJetPartons.particles = "prunedGenParticles"

#adjust PV used for Jet Corrections
process.patJetCorrFactorsAK4PFCHS.primaryVertices = "offlineSlimmedPrimaryVertices"

##---- AK8 ----
process.ak8PFJetsCHS = ak4PFJets.clone(
    src = 'chs',
    #jetPtMin = cms.double(10.0),   
    doAreaFastjet = cms.bool(True),
    #doRhoFastjet = cms.bool(True),  
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("AntiKt"),
    )
process.ak8GenJetsMy = ak4GenJets.clone(
    src = 'packedGenParticles',
    rParam = cms.double(0.8),
    )

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(
   process,
   postfix   = "",
   labelName = 'AK8PFCHS',
   jetSource = cms.InputTag('ak8PFJetsCHS'),
   trackSource = cms.InputTag('unpackedTracksAndVertices'), 
   pvSource = cms.InputTag('unpackedTracksAndVertices'), 
   jetCorrections = ('AK8PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-2'),
   btagDiscriminators = [      'combinedSecondaryVertexBJetTags'     ]
   ,algo= 'AK', rParam = 0.8
   )

#adjust MC matching
process.patJetGenJetMatchAK8PFCHS.matched = "ak8GenJetsMy"
process.patJetPartonMatchAK8PFCHS.matched = "prunedGenParticles"
process.patJetPartons.particles = "prunedGenParticles"

#adjust PV used for Jet Corrections
process.patJetCorrFactorsAK8PFCHS.primaryVertices = "offlineSlimmedPrimaryVertices"

##-- B-tagging for all jets

# the following part is needed if you want to run b-tagging on the freshly made jets
# CAVEAT: it is not 100% the same b-tagging as in RECO, but performance plots are almost identical
# As tracks are not stored in miniAOD, and b-tag fwk for CMSSW < 72X does not accept candidates
# we need to recreate tracks and pv for btagging in standard reco format:
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')
process.combinedSecondaryVertex.trackMultiplicityMin = 1  #needed for CMSSW < 71X
#new PAT default running is "unscheduled" so we just need to say in the outputCommands what we want to store
process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(True),
)


##-------------------- User analyzer  --------------------------------
process.dijets     = cms.EDAnalyzer('DijetTreeProducer',
  ## JETS/MET ########################################
  jetsAK4             = cms.InputTag('patJetsAK4PFCHS'),
  jetsAK8         = cms.InputTag('patJetsAK8PFCHS'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  #mjjMin           = cms.double(700),
  #dEtaMax          = cms.double(1.3),
  ## MC ########################################
  pu               = cms.untracked.InputTag('addPileupInfo'),
  ## trigger ###################################
  triggerAlias     = cms.vstring('Fat','PFHT650','PFNoPUHT650','HT750','HT550'),
  triggerSelection = cms.vstring(
    'HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v*',
     #'HLT_PFHT650_v*', #giulia : commented because not found in new entuples
     ### giulia
    'HLT_HT650_v*',
     ### end giulia
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

process.p = cms.Path(process.chs + 
                     process.ak4PFJetsCHS +
                     process.ak4GenJetsMy + 
                     process.ak8PFJetsCHS +
                     process.ak8GenJetsMy + 
                     process.dijets )


