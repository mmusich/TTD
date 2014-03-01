import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# same as input file
process.GlobalTag.globaltag = 'PRE_ST62_V8::All'
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
    "/store/relval/CMSSW_7_0_0_pre10/RelValTTbar/GEN-SIM-RECO/START70_V3-v1/00000/94F51711-7165-E311-99E4-002618943966.root"
    #        '/store/relval/CMSSW_7_0_0_pre5/RelValSingleMuPt10/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/26B1495D-F42F-E311-A3B5-002618943800.root'
    #        'file:../../refit_trk.root'
    )
)

# need to recreate hits from clusters (which are stored in RECO)
process.clustToHits = cms.Sequence(
    process.siPixelRecHits*process.siStripMatchedRecHits
)

# Create MeasurementTrackerEvent

process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

### Track Refitter
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.demo = cms.EDAnalyzer('DemoTrackAnalyzer',
                              tracks = cms.untracked.InputTag('generalTracks'),#'TrackRefitter'
                              TTRHBuilder = cms.string("WithAngleAndTemplate")
                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histo.root')
                                   )

#process.p = cms.Path(process.clustToHits*process.TrackRefitter*process.demo)
process.p = cms.Path(process.clustToHits *
                     process.MeasurementTrackerEvent *
                     process.trackingGlobalReco *
                     process.demo)

process.MessageLogger.categories.extend(["GetManyWithoutRegistration","GetByLabelWithoutRegistration"])
_messageSettings = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1),
    optionalPSet = cms.untracked.bool(True),
    limit = cms.untracked.int32(10000000)
    )

process.MessageLogger.cerr.GetManyWithoutRegistration = _messageSettings
process.MessageLogger.cerr.GetByLabelWithoutRegistration = _messageSettings

# Local Variables:
# truncate-lines: t
# show-trailing-whitespace: t
# End:
