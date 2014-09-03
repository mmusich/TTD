import FWCore.ParameterSet.Config as cms

process = cms.Process("RERECO")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_7_0_0_pre5/RelValSingleMuPt10/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/26B1495D-F42F-E311-A3B5-002618943800.root'
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# same as input file
process.GlobalTag.globaltag = 'PRE_ST62_V8::All'
# or get it automatically
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

### standard includes (can't live without)
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### define tracking sequences to run
# need to recreate hits from clusters (which are stored in RECO)
process.clustToHits = cms.Sequence(
    process.siPixelRecHits*process.siStripMatchedRecHits
)
# this is the official tracking sequence
process.tracking = cms.Sequence(
    process.trackingGlobalReco
)

process.pixellayertriplets.layerList = cms.vstring('BPix1+BPix2+FPix1_pos',
                                                   'BPix1+BPix2+FPix1_neg',
                                                   'BPix1+FPix1_pos+FPix2_pos',
                                                   'BPix1+FPix1_neg+FPix2_neg')

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
process.schedule = cms.Schedule(process.p, process.o)

### debug time and memory information (more at Validation/Performance/python//TimeMemory*.py)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.SimpleMemoryCheck=cms.Service("SimpleMemoryCheck",
                                      ignoreTotal=cms.untracked.int32(1),
                                      oncePerEventMode=cms.untracked.bool(False))

# Local Variables:
# truncate-lines: t
# show-trailing-whitespace: t
# End:
