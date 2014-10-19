import FWCore.ParameterSet.Config as cms

# inspired from RecoTracker/TrackProducer/test/TrackRefit.py

process = cms.Process("Refitting")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### Standard Configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# Load Manually
#process.GlobalTag.globaltag = 'POSTLS172_V4::All'
# or get it automatically
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# This is needed for tracking to work properly
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

### Track Refitter
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:reco_trk.root')
)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *_*_*_*',
                                                                      'keep recoTracks_*_*_*',
                                                                      'keep recoTrackExtras_*_*_*',
                                                                      'keep TrackingRecHitsOwned_*_*_*'),
                               fileName = cms.untracked.string('refit_trk.root')
                               )

process.p = cms.Path(process.MeasurementTrackerEvent
                     * process.TrackRefitter)
process.o = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.p,process.o)

### CUSTOMIZATION FUNCTIONS TO BE USED ###

# CUSTOMIZE LAYERS USED TO PRODUCE SEEDS IN THE INITIAL_STEP
def customize_UseAnalyticalPropagator(process):
    process.TrackRefitter.Propagator = cms.string('AnalyticalPropagator')
    process.RKTrajectoryFitter.Propagator = cms.string('AnalyticalPropagator')
    process.RKTrajectorySmoother.Propagator = cms.string('AnalyticalPropagator')
    process.out.fileName = cms.untracked.string('refit_trk_analytical_propagator.root')
    return process

### INVOKE APPROPRIATE CUSTOMIZATION FUNCTION(S) ###

# process = customize_UseAnalyticalPropagator(process)



