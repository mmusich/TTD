import FWCore.ParameterSet.Config as cms
from TTD.DemoTrackAnalyzer.demo_customizations import *

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

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

process.selectedTracksInitialStep = cms.EDFilter("TrackSelector",
                                                 src = cms.InputTag("generalTracks"),
                                                 cut = cms.string("quality('highPurity') & (algo = 4)")
                                                 )

process.selectedTracksPixelLessStep = cms.EDFilter("TrackSelector",
                                                 src = cms.InputTag("generalTracks"),
                                                 cut = cms.string("quality('highPurity') & (algo = 9)")
                                                 )

process.demo = cms.EDAnalyzer('DemoTrackAnalyzer',
                              tracks = cms.untracked.InputTag('TrackRefitter'),#'TrackRefitter'
                              seed = cms.untracked.InputTag('pixelLessStepSeeds'),
                              tracks_for_seed = cms.untracked.InputTag('selectedTracksPixelLessStep'),
                              do_rereco = cms.untracked.bool(False)
                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('trackAnalysis_PixelLessStep_small.root')
                                   )

process.p = cms.Path(process.clustToHits *
                     process.MeasurementTrackerEvent *
                     process.TrackRefitter *
                     process.selectedTracksPixelLessStep *
                     process.demo)

### CUSTOMIZE PROCESS

# Select specific Tracks out of the passed collection
# process = customizeSelectHPTracks(process, 'TrackRefitter')

# Re-run **FULL** tracking
process = customizeForSeeds(process)

# Enable MessageLogger for specific LogDebug messages in various
# tracking modules (require local installation&compilation of those
# modules with -DEDM_ML_DEBUG set.)
#process = customizeMessageLogger(process)

# Local Variables:
# truncate-lines: t
# show-trailing-whitespace: t
# End:
