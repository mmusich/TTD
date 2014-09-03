import FWCore.ParameterSet.Config as cms

process = cms.Process("MULTITRACKVALIDATOR")

# message logger
process.MessageLogger = cms.Service("MessageLogger",
     default = cms.untracked.PSet( limit = cms.untracked.int32(10) )
)

# source
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/26E2FF89-C22F-E311-8629-002618943901.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/2EE0C319-C82F-E311-B2C7-0026189438A9.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-RECO/PRE_ST62_V8-v1/00000/44313B0F-C72F-E311-8F88-002618943940.root'
] );


secFiles.extend( [
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/128849D8-AC2F-E311-A7E6-003048FFD76E.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/3236B511-AD2F-E311-BFF3-0030486792B4.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/38223CDA-AC2F-E311-AEBF-003048FFD752.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/862D57E2-B22F-E311-B088-003048FFD720.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/9E15651B-AD2F-E311-9956-0025905938A8.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/B46A2B17-AD2F-E311-8416-002590596490.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/C08E3BD7-AC2F-E311-9E1D-00259059642A.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/C2D019D3-AC2F-E311-8057-002618943834.root',
        '/store/relval/CMSSW_7_0_0_pre5/RelValTTbar/GEN-SIM-DIGI-RAW-HLTDEBUG/PRE_ST62_V8-v1/00000/FAA408D8-AC2F-E311-B255-002590593902.root'
] );

process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'PRE_ST62_V8::All'

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### validation-specific includes
#process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("Validation.RecoTrack.cuts_cff")
process.load("Validation.RecoTrack.MultiTrackValidator_cff")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("Validation.Configuration.postValidation_cff")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

########### configuration MultiTrackValidator ########
process.multiTrackValidator.outputFile = 'multitrackvalidator.root'
process.multiTrackValidator.associators = ['quickTrackAssociatorByHits']
process.multiTrackValidator.skipHistoFit=cms.untracked.bool(False)
#process.cutsRecoTracks.quality = cms.vstring('','highPurity')
#process.cutsRecoTracks.quality = cms.vstring('')
process.multiTrackValidator.label = ['cutsRecoTracks']
process.multiTrackValidator.useLogPt=cms.untracked.bool(True)
process.multiTrackValidator.minpT = cms.double(0.1)
process.multiTrackValidator.maxpT = cms.double(3000.0)
process.multiTrackValidator.nintpT = cms.int32(40)
process.multiTrackValidator.UseAssociators = cms.bool(True)
process.multiTrackValidator.runStandalone = cms.bool(True)

#process.load("Validation.RecoTrack.cuts_cff")
#process.cutsRecoTracks.ptMin    = cms.double(0.5)
#process.cutsRecoTracks.minHit   = cms.int32(10)
#process.cutsRecoTracks.minRapidity  = cms.int32(-1.0)
#process.cutsRecoTracks.maxRapidity  = cms.int32(1.0)

process.quickTrackAssociatorByHits.useClusterTPAssociation = cms.bool(True)
process.load("SimTracker.TrackerHitAssociation.clusterTpAssociationProducer_cfi")

process.validation = cms.Sequence(
    process.tpClusterProducer *
    process.multiTrackValidator
)

# paths
process.p = cms.Path(
      process.cutsRecoTracks
    * process.validation
)
process.schedule = cms.Schedule(
      process.p
)



# Local Variables:
# truncate-lines: t
# show-trailing-whitespace: t
# End:
