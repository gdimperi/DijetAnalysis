import FWCore.ParameterSet.Config as cms 

process = cms.Process('jetToolbox')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

## ----------------- Global Tag ------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.GlobalTag.globaltag = 'START53_V27::All'
#process.GlobalTag.globaltag = 'START70_V6::All'
#process.GlobalTag.globaltag = 'POSTLS170_V5::All'
#process.GlobalTag.globaltag = 'POSTLS170_V7::All'
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'

#--------------------- Report and output ---------------------------

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService=cms.Service("TFileService",
                                 fileName=cms.string('dijetTree_signal.root'),
                                 closeFileFast = cms.untracked.bool(True)
                                 )

## --- suppress long output ---> wantSummary = cms.untracked.bool(False) 

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False),
)

############## output test giulia edm format ###############
process.out = cms.OutputModule('PoolOutputModule',                                                                                                                  
                               fileName = cms.untracked.string('jettoolbox.root'),                                                                              
                               outputCommands = cms.untracked.vstring(['keep *_ak4PFJetsCHS_*_*',                                                                    
                                                                       'keep *_patJetsAK4PFCHS_*_*',                                                                  
                                                                       'keep *_ca8PFJetsCHS_*_*',                                                                     
                                                                       'keep *_patJetsCA8PFCHS_*_*',                                                                  
                                                                       'keep *_ak8PFJetsCHS_*_*',                                                                     
                                                                       'keep *_patJetsAK8PFCHS_*_*',                                                                  
                                                                       ])                                                                                           
                               )                                                                                                                                     
process.endpath = cms.EndPath(process.out)    


# ----------------------- Jet Tool Box  -----------------

##Load the toolBoxMiniHelper
#from RecoJets.JetProducers.jettoolboxMiniHelper_cff import * 
##for some reason doesn't work  ?__?
## just copy&paste the cff here

process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV'))

process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')

process.ak4PFJets.src = 'packedPFCandidates'
process.ak4PFJets.doAreaFastjet = True

process.ak4PFJetsCHS = process.ak4PFJetsCHS.clone(src = 'chs', doAreaFastjet = True) #boh ?
process.ak8PFJetsCHS = process.ak4PFJetsCHS.clone(src = 'chs', doAreaFastjet = True, rParam = 0.8)
process.ca8PFJetsCHS = process.ca4PFJets.clone(src = 'chs', doAreaFastjet = True, rParam = 0.8)

process.ak4GenJets.src = 'packedGenParticles'
process.ak8GenJets = process.ak4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)
process.ca8GenJets = process.ca4GenJets.clone(src = 'packedGenParticles', rParam = 0.8)

process.fixedGridRhoFastjetAll.pfCandidatesTag = 'packedPFCandidates'

# process.ak4PFJetsCHSPruned = ak5PFJetsCHSPruned.clone(
#     src='chs'
#     rParam = 0.4, 
#     jetPtMin = 15.0
#     )
# process.ak4PFJetsCHSFiltered = ak5PFJetsCHSFiltered.clone(
#     src='chs'
#     rParam = 0.4, 
#     jetPtMin = 15.0 
#     )                                                                             
# process.ak4PFJetsCHSTrimmed = ak5PFJetsCHSTrimmed.clone(
#     src='chs'
#     rParam = 0.4,
#     jetPtMin = 15.0  
#     )

process.ak8PFJetsCHSPruned.src = 'chs'
process.ak8PFJetsCHSTrimmed.src = 'chs'
process.ak8PFJetsCHSFiltered.src = 'chs'

process.ca8PFJetsCHSPruned.src = 'chs'
process.ca8PFJetsCHSTrimmed.src = 'chs'
process.ca8PFJetsCHSFiltered.src = 'chs'

process.cmsTopTagPFJetsCHS.src = 'chs'

from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection

addJetCollection(
    process,
    labelName = 'AK4PFCHS',
    jetSource = cms.InputTag('ak4PFJetsCHS'),
    algo = 'ak4',
    rParam = 0.4,
    jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag('unpackedTracksAndVertices'),
    btagDiscriminators = ['combinedSecondaryVertexBJetTags'],
    ) 

addJetCollection(
    process,
    labelName = 'CA8PFCHS',
    jetSource = cms.InputTag('ca8PFJetsCHS'),
    algo = 'ca8',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag('unpackedTracksAndVertices'),
    )                                                                                                                                                                                  
addJetCollection(
    process,
    labelName = 'AK8PFCHS',
    jetSource = cms.InputTag('ak8PFJetsCHS'),
    algo = 'ak8',
    rParam = 0.8,
    jetCorrections = ('AK7PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag('unpackedTracksAndVertices'),
    ) 
"""
switchJetCollection(
    process,
    jetSource = cms.InputTag('ak4PFJets'),
    algo = 'ak4',
    rParam = 0.4,
    jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'Type-1'),
    # btagDiscriminators = ['jetBProbabilityBJetTags',
    #                       'jetProbabilityBJetTags',
    #                       'trackCountingHighPurBJetTags',
    #                       'trackCountingHighEffBJetTags',
    #                       'simpleSecondaryVertexHighEffBJetTags',
    #                       'simpleSecondaryVertexHighPurBJetTags',
    #                       'combinedSecondaryVertexBJetTags'
    #                       ],
    trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pvSource = cms.InputTag('unpackedTracksAndVertices'),
    )
"""

process.patJetsAK4PFCHS.addJetCharge   = False
process.patJetsAK4PFCHS.addBTagInfo    = True
process.patJetsAK4PFCHS.getJetMCFlavour = False
process.patJetsAK4PFCHS.addAssociatedTracks = False
process.patJetPartonMatchAK4PFCHS.matched='prunedGenParticles'
process.patJetCorrFactorsAK4PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

process.patJetsCA8PFCHS.addJetCharge   = False
process.patJetsCA8PFCHS.addBTagInfo    = False   #For some reason this has to be False
process.patJetsCA8PFCHS.getJetMCFlavour = False
process.patJetsCA8PFCHS.addAssociatedTracks = False
process.patJetPartonMatchCA8PFCHS.matched='prunedGenParticles'
process.patJetCorrFactorsCA8PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

process.patJetsAK8PFCHS.addJetCharge   = False
process.patJetsAK8PFCHS.addBTagInfo    = False    #For some reason this has to be False
process.patJetsAK8PFCHS.getJetMCFlavour = False
process.patJetsAK8PFCHS.addAssociatedTracks = False
process.patJetPartonMatchAK8PFCHS.matched='prunedGenParticles'
process.patJetCorrFactorsAK8PFCHS.primaryVertices = 'offlineSlimmedPrimaryVertices'

process.patJets.addJetCharge   = False
process.patJets.addBTagInfo    = True
process.patJets.getJetMCFlavour = False
process.patJets.addAssociatedTracks = False
process.patJetPartonMatch.matched = 'prunedGenParticles'
process.patJetCorrFactors.primaryVertices = 'offlineSlimmedPrimaryVertices'


process.load('RecoBTag.Configuration.RecoBTag_cff')
process.load('RecoJets.Configuration.RecoJetAssociations_cff')
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

process.ak4JetTracksAssociatorAtVertexPF.jets = cms.InputTag('ak4PFJetsCHS')
process.ak4JetTracksAssociatorAtVertexPF.tracks = cms.InputTag('unpackedTracksAndVertices')
process.ak8JetTracksAssociatorAtVertexPF=process.ak4JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ak8PFJetsCHS'),
                                                                                        coneSize = 0.8)
process.ca8JetTracksAssociatorAtVertexPF=process.ak4JetTracksAssociatorAtVertexPF.clone(jets = cms.InputTag('ca8PFJetsCHS'),
                                                                                        coneSize = 0.8)

process.impactParameterTagInfos.primaryVertex = cms.InputTag('unpackedTracksAndVertices')
process.inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag('unpackedTracksAndVertices','secondary','')
process.combinedSecondaryVertex.trackMultiplicityMin = 1 #silly sv, uses un filtered tracks.. i.e. any pt

process.load('FWCore.MessageLogger.MessageLogger_cfi')
#process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter','manystripclus53X','toomanystripclus53X')
#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
#process.options.allowUnscheduled = cms.untracked.bool(True)


#Load the toolbox
process.load('RecoJets.JetProducers.jettoolbox_cff')


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#Test with Edmund Code... not working
#Load the GenEvent and GenPraticles cff files
#process.load('CMSROMA.DijetAnalysis.RootTupleMakerV2_GenEventInfo_cfi')
#process.load('CMSROMA.DijetAnalysis.RootTupleMakerV2_GenParticles_cfi')

#------------- Recluster Gen Jets to access the constituents -------
#already in toolbox, just add keep statements

process.out.outputCommands.append("keep *_ak4GenJets_*_*")
process.out.outputCommands.append("keep *_ak8GenJets_*_*")
process.out.outputCommands.append("keep *_ca8GenJets_*_*")

##-------------------- Define the source  ----------------------------

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:RSGravToJJ_kMpl01_M-1000_test.root')
    #fileNames = cms.untracked.vstring('file:/cmshome/santanas/CMS/data/Spring14miniaod__RSGravToJJ_kMpl01_M-1000_Tune4C_13TeV-pythia8__MINIAODSIM__PU20bx25_POSTLS170_V5-v1__00000__6AACD832-3707-E411-A167-001E672489D5.root')
    #fileNames = cms.untracked.vstring('file:/cmshome/santanas/CMS/data/Spring14drAODSIM__RSGravToJJ_kMpl01_M-1000_Tune4C_13TeV-pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__00000__0622C950-58E4-E311-A595-0025904B130A.root')
    fileNames = cms.untracked.vstring('file:/cmshome/gdimperi/data_test/6AACD832-3707-E411-A167-001E672489D5.root')    

)

# #Keep statements for valueMaps (link Reco::Jets to associated quantities)
# #You don't have to keep them to deswizzle below

#process.out.outputCommands += ['keep *_pileupJetIdEvaluator_*_*',
#                       'keep *_QGTagger_*_*']

process.out.outputCommands += ['keep *_NjettinessCA8_*_*',
                               'keep *_NjettinessAK8_*_*',
                               
#                                'keep *_QJetsAdderCA8_*_*',
                               'keep *_ca8PFJetsCHSPrunedLinks_*_*',
                               'keep *_ca8PFJetsCHSTrimmedLinks_*_*',
                               'keep *_ca8PFJetsCHSFilteredLinks_*_*',
                               'keep *_cmsTopTagPFJetsCHSLinksCA8_*_*',
#                                'keep *_QJetsAdderAK8_*_*',
                               'keep *_ak8PFJetsCHSPrunedLinks_*_*',
                               'keep *_ak8PFJetsCHSTrimmedLinks_*_*',
                               'keep *_ak8PFJetsCHSFilteredLinks_*_*',
                               'keep *_cmsTopTagPFJetsCHSLinksAK8_*_*',
                               ]


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#Deswizzle valueMaps and attach to PAT::Jets as userFloats

#process.patJetsAK4PFCHS.userData.userFloats.src += ['pileupJetIdEvaluator:fullDiscriminant','QGTagger:qgLikelihood']
#process.patJetsAK4PFCHS.userData.userInts.src   += ['pileupJetIdEvaluator:cutbasedId','pileupJetIdEvaluator:fullId']

#run these modules before pat in the sequence

process.patJetsCA8PFCHS.userData.userFloats.src += ['NjettinessCA8:tau1',
                                                    'NjettinessCA8:tau2',
                                                    'NjettinessCA8:tau3',
#                                                    'QJetsAdderCA8:QjetsVolatility',
                                                    'ca8PFJetsCHSPrunedLinks',
                                                    'ca8PFJetsCHSTrimmedLinks',
                                                    'ca8PFJetsCHSFilteredLinks',
                                                    'cmsTopTagPFJetsCHSLinksCA8'
                                                    ]



process.patJetsAK8PFCHS.userData.userFloats.src += ['NjettinessAK8:tau1',
                                                    'NjettinessAK8:tau2',
                                                    'NjettinessAK8:tau3',
 #                                                   'QJetsAdderAK8:QjetsVolatility',
                                                    'ak8PFJetsCHSPrunedLinks',
                                                    'ak8PFJetsCHSTrimmedLinks',
                                                    'ak8PFJetsCHSFilteredLinks',
                                                    'cmsTopTagPFJetsCHSLinksAK8'
                                                    ]





##-------------------- User analyzer  --------------------------------


process.dijets     = cms.EDAnalyzer('DijetTreeProducer',
  ## JETS/MET ########################################
  jetsAK4             = cms.InputTag('patJetsAK4PFCHS'),
  jetsAK8         = cms.InputTag('patJetsAK8PFCHS'),
  jetsCA8         = cms.InputTag('patJetsCA8PFCHS'),
  genJetsAK4             = cms.InputTag('ak4GenJets'),
  genJetsAK8         = cms.InputTag('ak8GenJets'),
  genJetsCA8         = cms.InputTag('ca8GenJets'),
  met              = cms.InputTag('slimmedMETs'),
  vtx              = cms.InputTag('offlineSlimmedPrimaryVertices'),
  ptMinAK4         = cms.double(10),
  ptMinAK8         = cms.double(10),
  ptMinCA8         = cms.double(10),
  #mjjMin           = cms.double(700),
  #dEtaMax          = cms.double(1.3),

  #forse non serve perche' gia' aggiunti al pat::jet -giulia-
  #Cutbasedid        = cms.InputTag("pileupJetIdEvaluator:cutbasedId"),
  #fullDiscriminant  = cms.InputTag("pileupJetIdEvaluator:fullDiscriminant")
  #fullId            = cms.InputTag("pileupJetIdEvaluator:fullId"),
  ## MC ########################################
  pu               = cms.untracked.InputTag('addPileupInfo'),
  ptHat            = cms.untracked.InputTag('generator'),

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


# ------------------ path --------------------------

process.p = cms.Path(process.chs * 
                     process.ak4PFJetsCHS *
                     process.ak4GenJets *
                     process.patJetsAK4PFCHS *
                     #process.patJetPartonMatchAK4PFCHS *
                     #process.patJetCorrFactorsAK4PFCHS* 
                     process.ak8PFJetsCHS *
                     process.ak8GenJets *
                     process.NjettinessAK8 *
                     #process.QJetsAdderAK8 * #questo rallenta - da capire
                     process.ak8PFJetsCHSPrunedLinks *
                     process.ak8PFJetsCHSFilteredLinks *
                     process.ak8PFJetsCHSTrimmedLinks *
                     process.patJetsAK8PFCHS *
                     process.ca8PFJetsCHS *
                     process.ca8GenJets * 
                     process.NjettinessCA8 *
                     #process.QJetsAdderCA8 * #questo rallenta - da capire
                     process.ca8PFJetsCHSPrunedLinks *
                     process.ca8PFJetsCHSFilteredLinks *
                     process.ca8PFJetsCHSTrimmedLinks *
                     process.patJetsCA8PFCHS *    

                     #process.pileupJetId        ##recipe not working for now
                     #process.puJetIdSequence    ##recipe not working for now
                     #process.pileupJetIdCalculator *  ##recipe not working for now
                     #process.pileupJetIdEvaluator     ##recipe not working for now
                     #process.QGTagger *
                     process.dijets 
                     ) 




