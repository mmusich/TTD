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

### Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

### Track Refitter
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:reco_trk.root')
) 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.out = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('drop *_*_*_*', 
                                                                      'keep recoTracks_*_*_*',
                                                                      'keep recoTrackExtras_*_*_*',
                                                                      'keep TrackingRecHitsOwned_*_*_*'),
                               fileName = cms.untracked.string('refit_trk.root')
                               )

process.p = cms.Path(process.TrackRefitter)
process.o = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.p,process.o)

 
