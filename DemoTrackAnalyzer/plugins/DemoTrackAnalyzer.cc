// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"

//
// class declaration
//

class DemoTrackAnalyzer : public edm::EDAnalyzer {
 public:
  explicit DemoTrackAnalyzer(const edm::ParameterSet&);
  ~DemoTrackAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  float computeMinimumTrackDistance(reco::Track const&, reco::Track const&,
                                    const MagneticField*);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::InputTag trackTags_;  // used to select what tracks to read from
                             // configuration file
  edm::EDGetTokenT<reco::TrackCollection> trackCollection_token_;

  TH1I* h_crossed_;
  TH1I* h_missed_;
  TH1F* h_charge_;
  TH1F* h_track_pt_;
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
    : trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")) {
  edm::Service<TFileService> fs;
  h_charge_ = fs->make<TH1F>("charge", "Charges", 200, -2., 2.);
  h_track_pt_ = fs->make<TH1F>("Track_Pt", "Track_Pt", 200, 0., 100.);
  h_crossed_ = fs->make<TH1I>("CrossedLayers", "CrossedLayers", 50, 0., 50.);
  h_missed_ = fs->make<TH1I>("MissedLayers", "MissedLayers", 50, 0., 50.);

  // Declare what we need to consume.
  using edm::InputTag;
  using reco::TrackCollection;

  trackCollection_token_ = consumes<TrackCollection>(trackTags_);
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

  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  Handle<TrackCollection> tracks;
  iEvent.getByToken(trackCollection_token_, tracks);

  // Do not attempt to get collection from the event unless we know
  // they are there since we are running the full reconstruction
  // sequence
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
    h_missed_->Fill(p.numberOfHits(reco::HitPattern::MISSING_OUTER_HITS));

    // Two track minimum distance
    if (&itTrack != &trk1 && do_two_tracks_min_distance) {
      do_two_tracks_min_distance = false;
      auto trk2 = itTrack;
      std::cout << "Minimum distance between the first two tracks is: "
                << computeMinimumTrackDistance(
                       trk1, trk2, magneticField.product()) << std::endl;
    }
  }  //  End loop over tracks (and associated trajectories)
}

// ------------ method called once each job just before starting event loop
// ------------
void DemoTrackAnalyzer::beginJob() {}

// ------------ method called once each job just after ending the event loop
// ------------
void DemoTrackAnalyzer::endJob() {}

// ------------ method called when starting to processes a run  ------------

void DemoTrackAnalyzer::beginRun(edm::Run const& run,
                                 edm::EventSetup const& setup) {
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

