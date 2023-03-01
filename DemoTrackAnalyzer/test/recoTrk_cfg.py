import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process("ReReco", Run3)

#from Configuration.Eras.Era_Run3_noMkFit_cff import Run3_noMkFit
#process = cms.Process("ReReco", Run3_noMkFit)

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_13_0_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_130X_mcRun3_2022_realistic_v2-v1/00000/0e8ae23c-7729-4ba7-b3fd-91d3183893d3.root'
        #'/store/relval/CMSSW_7_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3_CondDBv2-v1/00000/7E6493F0-5F21-E411-90CD-0025905A613C.root',
        #'/store/relval/CMSSW_7_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3_CondDBv2-v1/00000/A22ECB6B-5F21-E411-9968-0025905B8576.root',
        #'/store/relval/CMSSW_7_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3_CondDBv2-v1/00000/C0FABB8A-B021-E411-80D0-0025905B85A2.root',
        #'/store/relval/CMSSW_7_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3_CondDBv2-v1/00000/EE538882-B021-E411-B970-0025905B85AE.root'
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# Load Manually
#process.GlobalTag.globaltag = 'POSTLS172_V4::All'
# or get it automatically
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2022_realistic_v2', '')

### standard includes (can't live without)
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### define tracking sequences to run
# need to recreate hits from clusters (which are stored in RECO)
process.clustToHits = cms.Sequence(
    process.siPixelRecHits*process.siStripMatchedRecHits
)

noTrackingAndDependent = list()
noTrackingAndDependent.append(process.siPixelClustersPreSplitting)
noTrackingAndDependent.append(process.siStripZeroSuppression)
noTrackingAndDependent.append(process.siStripClusters)
noTrackingAndDependent.append(process.initialStepSeedLayersPreSplitting)
noTrackingAndDependent.append(process.trackerClusterCheckPreSplitting)
noTrackingAndDependent.append(process.initialStepTrackingRegionsPreSplitting)
noTrackingAndDependent.append(process.initialStepHitDoubletsPreSplitting)
noTrackingAndDependent.append(process.initialStepHitTripletsPreSplitting)
noTrackingAndDependent.append(process.initialStepSeedsPreSplitting)
noTrackingAndDependent.append(process.initialStepTrackCandidatesPreSplitting)
noTrackingAndDependent.append(process.initialStepTracksPreSplitting)
noTrackingAndDependent.append(process.firstStepPrimaryVerticesPreSplitting)
noTrackingAndDependent.append(process.initialStepTrackRefsForJetsPreSplitting)
noTrackingAndDependent.append(process.caloTowerForTrkPreSplitting)
noTrackingAndDependent.append(process.ak4CaloJetsForTrkPreSplitting)
noTrackingAndDependent.append(process.jetsForCoreTrackingPreSplitting)
noTrackingAndDependent.append(process.siPixelClusterShapeCachePreSplitting)
noTrackingAndDependent.append(process.siPixelClusters)
noTrackingAndDependent.append(process.clusterSummaryProducer)
noTrackingAndDependent.append(process.siPixelRecHitsPreSplitting)
noTrackingAndDependent.append(process.MeasurementTrackerEventPreSplitting)
noTrackingAndDependent.append(process.PixelLayerTriplets)
noTrackingAndDependent.append(process.pixelTracks)
noTrackingAndDependent.append(process.pixelVertices)

process.tracking_fromRECO = process.trackingGlobalReco.copyAndExclude(noTrackingAndDependent)

# This is needed for tracking to work properly
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")
process.load("RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi")
# this is the official tracking sequence
process.tracking = cms.Sequence(
    process.MeasurementTrackerEvent
    + process.siPixelClusterShapeCache
    #+ process.trackingGlobalReco
    + process.tracking_fromRECO
)


### define output EDMFile
MyReRecoEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_siPixelClusters_*_*',
    'keep *_siStripClusters_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_pixelVertices_*_*',
    'keep *_*_*_RERECO'),
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize=cms.untracked.int32(5*1024*1024)
)
process.out = cms.OutputModule("PoolOutputModule",
    MyReRecoEventContent,
    fileName = cms.untracked.string('reco_trk.root')
)

### paths
process.p = cms.Path(process.clustToHits * process.tracking)
process.o = cms.EndPath(process.out)

### sequence of paths to run
process.schedule = cms.Schedule(process.p,
                                process.o)

### debug time and memory information (more at Validation/Performance/python//TimeMemory*.py)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.SimpleMemoryCheck=cms.Service("SimpleMemoryCheck",
                                      ignoreTotal=cms.untracked.int32(1),
                                      oncePerEventMode=cms.untracked.bool(False))

### CUSTOMIZATION FUNCTIONS TO BE USED ###

# CUSTOMIZE LAYERS USED TO PRODUCE SEEDS IN THE INITIAL_STEP
def customize_initialStepSeeds(process):
    process.initialStepSeedLayers.layerList.remove('BPix1+BPix2+BPix3')
    process.out.fileName = cms.untracked.string('reco_trk_mod.root')
    return process

### INVOKE APPROPRIATE CUSTOMIZATION FUNCTION(S) ###

# process = customize_initialStepSeeds(process)



# Local Variables:
# truncate-lines: t
# show-trailing-whitespace: t
# End:
