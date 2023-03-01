// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"

//
// class declaration
//

class DemoTrackAnalyzer : public edm::one::EDAnalyzer<> {
 public:
  explicit DemoTrackAnalyzer(const edm::ParameterSet&);
  ~DemoTrackAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  float computeMinimumTrackDistance(reco::Track const&, reco::Track const&,
                                    const MagneticField*);
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ----------member data ---------------------------
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopologyToken_;

  edm::InputTag trackTags_;  // used to select what tracks to read from
                             // configuration file
  edm::InputTag genericSeedsTag_;
  edm::InputTag genericStepTracksTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackCollection_token_;
  edm::EDGetTokenT<SiPixelRecHitCollection> pixelHits_token_;
  edm::EDGetTokenT<reco::TrackCollection> genericStepTracks_token_;
  edm::EDGetTokenT<TrajectorySeedCollection> genericSeeds_token_;
  edm::EDGetTokenT<TrajTrackAssociationCollection> trajTrackAssociation_token_;

  TH1I* h_crossed_;
  TH1I* h_missed_;
  TH1I* h_hit_pixel_layers_;
  TH1F* h_charge_;
  TH1F* h_track_pt_;
  TH1F* h_seed_pt_;
  TH1F* h_eta_genericstep_seeds_;
  TH1F* h_eta_genericstep_hp_;
  TH1F* h_tob_xpull_;
  TH2F* h_hit_map_;
  TH2F* h_hit_pixelbarrel_map_;
  bool do_rereco_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoTrackAnalyzer::DemoTrackAnalyzer(const edm::ParameterSet& iConfig)
  : magFieldToken_(esConsumes()),
    trackerTopologyToken_(esConsumes()),
    trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
    genericSeedsTag_(iConfig.getUntrackedParameter<edm::InputTag>("seed")),
    genericStepTracksTag_(iConfig.getUntrackedParameter<edm::InputTag>("tracks_for_seed")),
    do_rereco_(iConfig.getUntrackedParameter<bool>("do_rereco")) {
  edm::Service<TFileService> fs;
  h_charge_ = fs->make<TH1F>("charge", "Charges", 200, -2., 2.);
  h_track_pt_ = fs->make<TH1F>("Track_Pt", "Track_Pt", 200, 0., 100.);
  h_seed_pt_ = fs->make<TH1F>("Seed_Pt", "Seed_Pt", 200, 0., 100.);
  h_crossed_ = fs->make<TH1I>("CrossedLayers", "CrossedLayers", 50, 0., 50.);
  h_missed_ = fs->make<TH1I>("MissedLayers", "MissedLayers", 50, 0., 50.);
  h_hit_map_ =
      fs->make<TH2F>("Hit_map", "Hit_map", 1200, -300., 300., 280, 0., 140.);
  h_hit_pixelbarrel_map_ = fs->make<TH2F>("Hit_pixelbarrel", "Hit_pixelbarrel",
                                          240, -12., 12., 240, -12., 12.);
  h_hit_pixel_layers_ =
      fs->make<TH1I>("Hits_pixellayers", "Hits_pixellayers", 5, 1, 6);
  h_eta_genericstep_seeds_ =
      fs->make<TH1F>("GenericSeed_Eta", "GenericSeed_Eta", 100, -2.5, 2.5);
  h_eta_genericstep_hp_ =
      fs->make<TH1F>("GenericHP_Eta", "GenericHP_Eta", 100, -2.5, 2.5);
  h_tob_xpull_ = fs->make<TH1F>("TOB_Pull_x", "TOB_Pull_x", 100, -5., 5.);

  // Declare what we need to consume.
  using edm::InputTag;
  using reco::TrackCollection;

  trackCollection_token_ = consumes<TrackCollection>(trackTags_);
  genericStepTracks_token_ = mayConsume<TrackCollection>(genericStepTracksTag_);
  genericSeeds_token_ =
      mayConsume<TrajectorySeedCollection>(genericSeedsTag_);
  pixelHits_token_ =
      consumes<SiPixelRecHitCollection>(InputTag("siPixelRecHits"));
  trajTrackAssociation_token_ =
      mayConsume<TrajTrackAssociationCollection>(trackTags_);
}

DemoTrackAnalyzer::~DemoTrackAnalyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

float DemoTrackAnalyzer::computeMinimumTrackDistance(
    reco::Track const& trk1, reco::Track const& trk2,
    const MagneticField* theMagField) {
  GlobalPoint position1(trk1.vertex().x(), trk1.vertex().y(),
                        trk1.vertex().z());

  GlobalVector momentum1(trk1.momentum().x(), trk1.momentum().y(),
                         trk1.momentum().z());

  GlobalTrajectoryParameters gtp1(position1, momentum1, trk1.charge(),
                                  theMagField);
  GlobalPoint position2(trk2.vertex().x(), trk2.vertex().y(),
                        trk2.vertex().z());

  GlobalVector momentum2(trk2.momentum().x(), trk2.momentum().y(),
                         trk2.momentum().z());

  GlobalTrajectoryParameters gtp2(position2, momentum2, trk2.charge(),
                                  theMagField);

  TwoTrackMinimumDistance two_track_min_distance;
  if (two_track_min_distance.calculate(gtp1, gtp2)) {
    return two_track_min_distance.distance();
  }
  return -1.0;
}

// ------------ method called for each event  ------------
void DemoTrackAnalyzer::analyze(const edm::Event& iEvent,
                                const edm::EventSetup& iSetup) {
  //  use strictly needed namespaces
  using edm::Handle;
  using edm::ESHandle;
  using edm::Ref;
  using reco::TrackCollection;

  ESHandle<MagneticField> magneticField = iSetup.getHandle(magFieldToken_);

  // Retrieve tracker topology from geometry
  const TrackerTopology* const tTopo = &iSetup.getData(trackerTopologyToken_);

  Handle<TrackCollection> tracks;
  iEvent.getByToken(trackCollection_token_, tracks);

  Handle<SiPixelRecHitCollection> pixelHits;
  iEvent.getByToken(pixelHits_token_, pixelHits);

  Handle<TrajectorySeedCollection> genericSeeds;
  Handle<TrackCollection> genericStepTracks;
  Handle<TrajTrackAssociationCollection> trajTrackAssociation;
  // Do not attempt to get collection from the event unless we know
  // they are there since we are running the full reconstruction
  // sequence
  TrajTrackAssociationCollection::const_iterator tji;
  if (do_rereco_) {
    iEvent.getByToken(genericSeeds_token_, genericSeeds);
    iEvent.getByToken(genericStepTracks_token_, genericStepTracks);
    iEvent.getByToken(trajTrackAssociation_token_, trajTrackAssociation);
    tji = trajTrackAssociation->begin();
  }

  reco::Track const& trk1 = *tracks->begin();
  bool do_two_tracks_min_distance = true;
  // auto --> "reco::Track"
  for (auto const& itTrack : *tracks) {
    // Plot track variables and hit-pattern
    h_charge_->Fill(itTrack.charge());
    h_track_pt_->Fill(itTrack.pt());

    // hit pattern of the track
    const reco::HitPattern& p = itTrack.hitPattern();
    h_crossed_->Fill(p.trackerLayersWithMeasurement());
    h_missed_->Fill(p.numberOfAllHits(reco::HitPattern::MISSING_OUTER_HITS));

    // Access to hit information
    auto bi = itTrack.recHitsBegin();
    auto be = itTrack.recHitsEnd();
    for (; bi != be; ++bi) {
      const TrackingRecHit & hit = **bi;
      if (hit.isValid())
        h_hit_map_->Fill(hit.globalPosition().z(),
                         hit.globalPosition().perp());
    }
    // Two track minimum distance
    if (&itTrack != &trk1 && do_two_tracks_min_distance) {
      do_two_tracks_min_distance = false;
      auto trk2 = itTrack;
      std::cout << "Minimum distance between the first two tracks is: "
                << computeMinimumTrackDistance(
                       trk1, trk2, magneticField.product()) << std::endl;
    }

    if (do_rereco_) {
      // Access to seeds
      // auto --> edm::RefToBase<TrajectorySeed>
      auto tkseed = itTrack.seedRef();
      // To pick-up the last hit on the seed you can either use
      // seed->recHits().first + (num_rechits_in_seed - 1) or use
      // seed->recHits().second - 1.
      TrajectoryStateOnSurface state = trajectoryStateTransform::transientState(
          tkseed->startingState(),
          (tkseed->recHits().end() - 1)->surface(),
          magneticField.product());
      h_seed_pt_->Fill(state.globalMomentum().perp());

      // Access to Trajectory: pull distribution for TOB
      Ref<std::vector<Trajectory> > traj = tji->key;
      for (auto const& measurement : traj->measurements()) {
        TrackingRecHit::ConstRecHitPointer hit = measurement.recHit();
        DetId hitId = hit->geographicalId();
        if (hit->isValid() &&
            hitId.subdetId() == static_cast<int>(StripSubdetector::TOB)) {
          TrajectoryStateOnSurface fwdState = measurement.forwardPredictedState();
          float delta = hit->localPosition().x() - fwdState.localPosition().x();
          float err2 = hit->localPositionError().xx() +
              fwdState.localError().positionError().xx();
          if (err2) h_tob_xpull_->Fill(delta / sqrt(err2));
        }
      }  // End Loop over trajectory measurements
      // Increment the TrajectoryTrackAssociation to be in sync w.r.t
      // the main track collection
      ++tji;
    }  // Endif do_rereco_
  }  //  End loop over tracks (and associated trajectories)

  if (do_rereco_) {
    for (auto const& a_seed : *genericSeeds) {
      TrajectoryStateOnSurface state = trajectoryStateTransform::transientState(
										a_seed.startingState(), (a_seed.recHits().end() - 1)->surface(),
          magneticField.product());
      h_eta_genericstep_seeds_->Fill(state.globalMomentum().eta());
    }
    for (auto const & a_track : *genericStepTracks) {
      h_eta_genericstep_hp_->Fill(a_track.eta());
    }
  }

  // Access to hit information: loop over all pixel hits
  // auto --> edmNew::DetSetVector<SiPixelRecHit>
  for (auto const& pixel_rechit_per_detid : *pixelHits) {
    DetId hitId = pixel_rechit_per_detid.detId();
    for (auto const& a_pixel_rechit : pixel_rechit_per_detid) {
      if (a_pixel_rechit.isValid()) {
        if (hitId.subdetId() ==
            static_cast<int>(PixelSubdetector::PixelBarrel)) {
          h_hit_pixelbarrel_map_->Fill(a_pixel_rechit.globalPosition().x(),
                                       a_pixel_rechit.globalPosition().y());
          h_hit_pixel_layers_->Fill(tTopo->pxbLayer(hitId));
        } else {
          h_hit_pixel_layers_->Fill(3 + tTopo->pxfDisk(hitId));
        }
      }
    }
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------
void DemoTrackAnalyzer::fillDescriptions(
    edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  // Specify that only 'tracks' is allowed
  // To use, remove the default given above and uncomment below
  // ParameterSetDescription desc;
  // desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  // descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(DemoTrackAnalyzer);

// Local Variables:
// truncate-lines: t
// show-trailing-whitespace: t
// End:
