import FWCore.ParameterSet.Config as cms

##--------- good primary vertices ---------------
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
goodOfflinePrimaryVertices = cms.EDFilter("PrimaryVertexObjectFilter",
    src          = cms.InputTag('offlinePrimaryVertices'),
    filterParams = pvSelector.clone( minNdof = cms.double(4.0), maxZ = cms.double(24.0) )
)

from CommonTools.ParticleFlow.pfNoPileUp_cff import * 
from CommonTools.ParticleFlow.pfParticleSelection_cff import *

pfPileUp.checkClosestZVertex = False
pfPileUp.Vertices = 'goodOfflinePrimaryVertices'
#giulia debug
#pfPileUp.PFCandidates = 'particleFlow'      
#pfNoPileUp.bottomCollection = 'particleFlow'
pfPileUp.PFCandidates = 'particleFlowPtrs'      
pfNoPileUp.bottomCollection = 'particleFlowPtrs'

from CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi import goodOfflinePrimaryVertices
pfNoPileUpSequence.insert(0, goodOfflinePrimaryVertices)

pileUpSubtractionSequence = cms.Sequence(
    pfNoPileUpSequence +
    pfParticleSelectionSequence
    )

# JETS  CA8 ----------------------------

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
ca8PFJetsCHS = ak5PFJets.clone(
  src = 'pfNoPileUp',
  jetPtMin = cms.double(30.0),
  doAreaFastjet = cms.bool(True),
  rParam = cms.double(0.8),
  jetAlgorithm = cms.string("CambridgeAachen"),
)

jetSource = 'ca8PFJetsCHS'

# corrections 
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
patJetCorrFactorsCA8CHS = patJetCorrFactors.clone()
patJetCorrFactorsCA8CHS.src = jetSource
# will need to add L2L3 corrections in the cfg
patJetCorrFactorsCA8CHS.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
patJetCorrFactorsCA8CHS.payload = 'AK7PFchs'
patJetCorrFactorsCA8CHS.useRho = True

from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import *
patJetsCA8CHS = patJets.clone()
patJetsCA8CHS.jetSource = jetSource
patJetsCA8CHS.addJetCharge = False
#giulia debug
#patJetsCA8CHS.embedCaloTowers = False
patJetsCA8CHS.embedPFCandidates = False
patJetsCA8CHS.addAssociatedTracks = False
patJetsCA8CHS.addBTagInfo = False
patJetsCA8CHS.addDiscriminators = False
patJetsCA8CHS.addJetID = False
patJetsCA8CHS.addGenPartonMatch = False
patJetsCA8CHS.embedGenPartonMatch = False
patJetsCA8CHS.addGenJetMatch = False
patJetsCA8CHS.getJetMCFlavour = False
patJetsCA8CHS.jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsCA8CHS'))

#### Adding Nsubjetiness

patJetsCA8CHSwithNsub = cms.EDProducer("NjettinessAdder",
  src=cms.InputTag("patJetsCA8CHS"),
  cone=cms.double(0.8)
)

# JETS PRUNED CA8 ----------------------------

from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
ca8PFJetsCHSpruned = ak5PFJetsPruned.clone(
  src = 'pfNoPileUp',
  jetPtMin = cms.double(30.0),
  doAreaFastjet = cms.bool(True),
  rParam = cms.double(0.8),
  jetAlgorithm = cms.string("CambridgeAachen"),
)

jetSource = 'ca8PFJetsCHSpruned'

# corrections 
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
patJetCorrFactorsCA8CHSpruned = patJetCorrFactors.clone()
patJetCorrFactorsCA8CHSpruned.src = jetSource
# will need to add L2L3 corrections in the cfg
patJetCorrFactorsCA8CHSpruned.levels = ['L1FastJet', 'L2Relative', 'L3Absolute']
patJetCorrFactorsCA8CHSpruned.payload = 'AK7PFchs'
patJetCorrFactorsCA8CHSpruned.useRho = True


patJetsCA8CHSpruned = patJets.clone()
patJetsCA8CHSpruned.jetSource = jetSource
patJetsCA8CHSpruned.addJetCharge = False
#giulia debug
#patJetsCA8CHSpruned.embedCaloTowers = False
patJetsCA8CHSpruned.embedPFCandidates = False
patJetsCA8CHSpruned.addAssociatedTracks = False
patJetsCA8CHSpruned.addBTagInfo = False
patJetsCA8CHSpruned.addDiscriminators = False
patJetsCA8CHSpruned.addJetID = False
patJetsCA8CHSpruned.addGenPartonMatch = False
patJetsCA8CHSpruned.embedGenPartonMatch = False
patJetsCA8CHSpruned.addGenJetMatch = False
patJetsCA8CHSpruned.getJetMCFlavour = False
patJetsCA8CHSpruned.jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsCA8CHSpruned'))


ca8Jets = cms.Sequence(
  goodOfflinePrimaryVertices +
  pfNoPileUpSequence +
  ca8PFJetsCHS + 
  patJetCorrFactorsCA8CHS +
  patJetsCA8CHS +  
  patJetsCA8CHSwithNsub +
  ca8PFJetsCHSpruned +
  patJetCorrFactorsCA8CHSpruned +
  patJetsCA8CHSpruned
)
