import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# Load Manually
#process.GlobalTag.globaltag = 'POSTLS172_V4::All'
# or get it automatically
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
### standard includes (can't live without)
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_7_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3_CondDBv2-v1/00000/7E6493F0-5F21-E411-90CD-0025905A613C.root',
    '/store/relval/CMSSW_7_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3_CondDBv2-v1/00000/A22ECB6B-5F21-E411-9968-0025905B8576.root',
    '/store/relval/CMSSW_7_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3_CondDBv2-v1/00000/C0FABB8A-B021-E411-80D0-0025905B85A2.root',
    '/store/relval/CMSSW_7_2_0_pre3/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS172_V3_CondDBv2-v1/00000/EE538882-B021-E411-B970-0025905B85AE.root'
    )
)

# need to recreate hits from clusters (which are stored in RECO)
process.clustToHits = cms.Sequence(
    process.siPixelRecHits
    *process.siStripMatchedRecHits
)

# This is needed for tracking to work properly
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

### Track Refitter
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.demo = cms.EDAnalyzer('DemoTrackAnalyzer',
                              tracks = cms.untracked.InputTag('TrackRefitter'),#'TrackRefitter'
                              TTRHBuilder = cms.string("WithAngleAndTemplate")
                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histo.root')
                                   )

process.p = cms.Path(process.clustToHits *
                     process.MeasurementTrackerEvent *
                     process.TrackRefitter *
                     #                     process.trackingGlobalReco *
                     process.demo)

process.MessageLogger.categories.extend(["GetManyWithoutRegistration","GetByLabelWithoutRegistration"])
_messageSettings = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1),
    optionalPSet = cms.untracked.bool(True),
    limit = cms.untracked.int32(10000000)
    )

process.MessageLogger.cerr.GetManyWithoutRegistration = _messageSettings
process.MessageLogger.cerr.GetByLabelWithoutRegistration = _messageSettings

def customizeSelectHPTracks(process):
    process.selectedTracks = cms.EDFilter("TrackSelector",
                                          src = cms.InputTag("TrackRefitter"),
                                          cut = cms.string("quality('highPurity') & (algo = 4) & abs(eta) < 0.9")
                                          )
    process.demo.tracks = cms.untracked.InputTag("selectedTracks")
    # Now run the filter that produces the input collection for our
    # demo analyzer before the demo itself.
    process.p.insert(process.p.index(process.demo), process.selectedTracks)
    return process

def customizeForSeeds(process):
    if not process.p.replace(process.TrackRefitter, process.trackingGlobalReco):
        raise Exception("Cannot customize process to run full reco instead of refit.")
    # Add also siPixelClusterShapeCache that is needed by the full tracking.
    process.load("RecoPixelVertexing.PixelLowPtUtilities.siPixelClusterShapeCache_cfi")
    process.p.insert(process.p.index(process.trackingGlobalReco), process.siPixelClusterShapeCache)

    # Check if we have been called after the customizeSelectHPTracks
    # has been called: if so, change the input collection of the
    # filter, otherwise change the input collection of the
    # DemoTrackAnalyzer module directly.
    if getattr(process, "selectedTracks", None):
        process.selectedTracks.src = cms.untracked.InputTag("generalTracks")
    else:
        process.demo.tracks = cms.untracked.InputTag("generalTracks")

    return process

### CUSTOMIZE PROCESS

#process = customizeSelectHPTracks(process)
process = customizeForSeeds(process)


# Local Variables:
# truncate-lines: t
# show-trailing-whitespace: t
# End:
