import FWCore.ParameterSet.Config as cms

def customizeMessageLogger(process):
    ### Easy customisation of MessageLogger ###
    # 1. Extend MessageLogger to monitor all modules: the * means any
    #    label for all defined python modules
    process.MessageLogger.debugModules.extend(['*'])
    # 2. Define destination and its default logging properties
    destination = 'debugTrackingModules'
    how_to_debug = cms.untracked.PSet(threshold = cms.untracked.string("DEBUG"),
                                      DEBUG = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                                      default = cms.untracked.PSet(limit = cms.untracked.int32(0)),
                                      )
    # 3. Attach destination and its logging properties to the main process
    process.MessageLogger.destinations.extend([destination])
    process.MessageLogger._Parameterizable__addParameter(destination, how_to_debug)
    # 4. Define and extend the categories we would like to monitor
    log_debug_categories = ['CkfTrackCandidateMakerBase',
                            'TkStripMeasurementDet',
                            'TkPixelMeasurementDet',
                            'CkfPattern']
    process.MessageLogger.categories.extend(log_debug_categories)

    # 5. Extend the configuration of the configured destination so that it
    #    will trace all messages coming from the list of specified
    #    categories.
    unlimit_debug = cms.untracked.PSet(limit = cms.untracked.int32(-1))
    for val in log_debug_categories:
        process.MessageLogger.debugTrackingModules._Parameterizable__addParameter(val, unlimit_debug)

    return process

def customizeSelectHPTracks(process, track_collection_name):
    process.selectedTracks = cms.EDFilter("TrackSelector",
                                          src = cms.InputTag(track_collection_name),
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
    process.demo.do_rereco = cms.untracked.bool(True)

    # Check if we have been called after the customizeSelectHPTracks
    # has been called: if so, change the input collection of the
    # filter, otherwise change the input collection of the
    # DemoTrackAnalyzer module directly.
    if getattr(process, "selectedTracks", None):
        process.selectedTracks.src = cms.InputTag("generalTracks")
    else:
        process.demo.tracks = cms.untracked.InputTag("generalTracks")
    process.demo.traj_tracks = cms.untracked.InputTag("generalTracks")

    return process

# Local Variables:
# truncate-lines: t
# show-trailing-whitespace: t
# End:
