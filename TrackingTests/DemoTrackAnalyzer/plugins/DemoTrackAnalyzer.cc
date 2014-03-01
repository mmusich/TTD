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
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
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

class DemoTrackAnalyzer : public edm::EDAnalyzer {
 public:
  explicit DemoTrackAnalyzer(const edm::ParameterSet&);
  ~DemoTrackAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


 private:
  float computeMinimumTrackDistance(reco::TrackCollection::const_iterator,
                                    reco::TrackCollection::const_iterator,
                                    const MagneticField*);
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::InputTag trackTags_; //used to select what tracks to read from configuration file
  edm::EDGetTokenT<reco::TrackCollection> trackCollection_token_;
  edm::EDGetTokenT<reco::TrackCollection> tracks_token_;
  edm::EDGetTokenT<SiPixelRecHitCollection> pixelHits_token_;
  edm::EDGetTokenT<TrajectorySeedCollection> initialStepSeeds_token_;
  edm::EDGetTokenT<TrajTrackAssociationCollection> trajTrackAssociation_token_;

  TH1I * h_crossed;
  TH1I * h_missed;
  TH1I * h_hit_pixel_layers;
  TH1F * histo;
  TH1F * h_track_pt;
  TH1F * h_seed_pt;
  TH1F * h_eta_initialstep_seeds;
  TH1F * h_tob_xpull;
  TH2F * h_hit_map;
  TH2F * h_hit_pixelbarrel_map;
  std::string builderName;
  const TransientTrackingRecHitBuilder* builder;
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
    :trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
     builderName(iConfig.getParameter<std::string>("TTRHBuilder"))
{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  histo      = fs->make<TH1F>("charge", "Charges", 200, -2., 2.);
  h_track_pt = fs->make<TH1F>("Track_Pt", "Track_Pt", 200, 0., 100.);
  h_seed_pt  = fs->make<TH1F>("Seed_Pt", "Seed_Pt", 200, 0., 100.);
  h_crossed  = fs->make<TH1I>("CrossedLayers", "CrossedLayers", 50, 0., 50.);
  h_missed   = fs->make<TH1I>("MissedLayers", "MissedLayers", 50, 0., 50.);
  h_hit_map  = fs->make<TH2F>("Hit_map", "Hit_map",
                             1200, -300., 300.,
                             280, 0., 140. );
  h_hit_pixelbarrel_map  = fs->make<TH2F>("Hit_pixelbarrel", "Hit_pixelbarrel",
                                          240, -12., 12.,
                                          240, -12., 12. );
  h_hit_pixel_layers = fs->make<TH1I>("Hits_pixellayers", "Hits_pixellayers",
                                      5, 1, 6);
  h_eta_initialstep_seeds = fs->make<TH1F>("InitialSeed_Eta",
                                           "InitialSeed_Eta",
                                           70, -3.5, 3.5);
  h_tob_xpull = fs->make<TH1F>("TOB_Pull_x",
                               "TOB_Pull_x",
                               100, -5., 5.);

  // Declare what we need to consume.
  using namespace edm;
  using reco::TrackCollection;

  trackCollection_token_ = consumes<TrackCollection>(trackTags_);
  tracks_token_ = consumes<TrackCollection>(trackTags_);
  pixelHits_token_ = consumes<SiPixelRecHitCollection>(edm::InputTag("siPixelRecHits"));
  initialStepSeeds_token_ = consumes<TrajectorySeedCollection>(edm::InputTag("initialStepSeeds"));
  trajTrackAssociation_token_ = consumes<TrajTrackAssociationCollection>(trackTags_);
}


DemoTrackAnalyzer::~DemoTrackAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

float DemoTrackAnalyzer::computeMinimumTrackDistance(reco::TrackCollection::const_iterator trk1,
                                                     reco::TrackCollection::const_iterator trk2,
                                                      const MagneticField* theMagField)
{
  GlobalPoint position1(trk1->vertex().x(),
                        trk1->vertex().y(),
                        trk1->vertex().z());

  GlobalVector momentum1(trk1->momentum().x(),
                         trk1->momentum().y(),
                         trk1->momentum().z());

  GlobalTrajectoryParameters gtp1(position1,
                                  momentum1,
                                  trk1->charge(),
                                  theMagField);
  GlobalPoint position2(trk2->vertex().x(),
                        trk2->vertex().y(),
                        trk2->vertex().z());

  GlobalVector momentum2(trk2->momentum().x(),
                         trk2->momentum().y(),
                         trk2->momentum().z());

  GlobalTrajectoryParameters gtp2(position2,
                                  momentum2,
                                  trk2->charge(),
                                  theMagField);

  TwoTrackMinimumDistance two_track_min_distance;
  if (two_track_min_distance.calculate(gtp1, gtp2)) {
    return two_track_min_distance.distance();
  }
  return -1.0;
}



// ------------ method called for each event  ------------
void
DemoTrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using reco::TrackCollection;

  ESHandle< MagneticField > magneticField;
  iSetup.get< IdealMagneticFieldRecord >().get( magneticField );

  Handle<TrackCollection> tracks;
  iEvent.getByToken(tracks_token_, tracks);

  Handle<SiPixelRecHitCollection> pixelHits;
  iEvent.getByToken(pixelHits_token_, pixelHits);

  Handle<TrajectorySeedCollection> initialStepSeeds;
  iEvent.getByToken(initialStepSeeds_token_, initialStepSeeds);

  Handle<TrajTrackAssociationCollection> trajTrackAssociation;
  iEvent.getByToken(trajTrackAssociation_token_, trajTrackAssociation);

  auto tji = trajTrackAssociation->begin();

  TrackCollection::const_iterator trk1 = tracks->begin();
  TrackCollection::const_iterator trk2 = tracks->begin();
  for (TrackCollection::const_iterator itTrack = tracks->begin();
      itTrack != tracks->end();
      ++itTrack, ++tji) {
    int charge = 0;
    charge = itTrack->charge();
    histo->Fill( charge );
    h_track_pt->Fill(itTrack->pt());
    // hit pattern of the track
    const reco::HitPattern& p = itTrack->hitPattern();
    h_crossed->Fill(p.trackerLayersWithMeasurement());
    h_missed->Fill(itTrack->trackerExpectedHitsOuter().numberOfHits());
    auto bi = itTrack->recHitsBegin();
    auto be = itTrack->recHitsEnd();
    for (; bi != be; ++bi) {
      TransientTrackingRecHit::RecHitPointer thit = builder->build(&**bi);
      if (thit->isValid())
        h_hit_map->Fill(thit->globalPosition().z(), thit->globalPosition().perp());
    }
    edm::RefToBase<TrajectorySeed> tkseed = itTrack->seedRef();
    if (tkseed->nHits() == 3) {
      TransientTrackingRecHit::RecHitPointer recHit =
          builder->build(&*(tkseed->recHits().first+2));
      TrajectoryStateOnSurface state = trajectoryStateTransform::transientState(
          tkseed->startingState(), recHit->surface(), magneticField.product());
      h_seed_pt->Fill(state.globalMomentum().perp());
    }
    Ref<std::vector<Trajectory> > traj = tji->key;
    std::vector<TrajectoryMeasurement> trajMeas = traj->measurements();
    auto tjmi = trajMeas.begin();
    auto tjme = trajMeas.end();
    for (; tjmi != tjme; ++tjmi) {
      TransientTrackingRecHit::ConstRecHitPointer hit = tjmi->recHit();
      DetId hitId = hit->geographicalId();
      if (hit->isValid() && hitId.subdetId() == (int) StripSubdetector::TOB) {
        TrajectoryStateOnSurface fwdState = tjmi->forwardPredictedState();
        float delta = hit->localPosition().x() - fwdState.localPosition().x();
        float err2 = hit->localPositionError().xx()
            + fwdState.localError().positionError().xx();
        if (err2)
          h_tob_xpull->Fill(delta/sqrt(err2));
      }
    }
    // Two track minimum distance
    if (itTrack != trk1 && trk2 == trk1) {
      trk2 = itTrack;
      std::cout << "Minimum distance between the first two tracks is: "
          << computeMinimumTrackDistance(trk1, trk2, magneticField.product())
          << std::endl;
    }
  } //  Loop over tracks (and associated trajectories)

  auto si = initialStepSeeds->begin();
  auto se = initialStepSeeds->end();
  for (; si != se; ++si) {
    TransientTrackingRecHit::RecHitPointer recHit =
        builder->build(&*(si->recHits().first+2));
    TrajectoryStateOnSurface state = trajectoryStateTransform::transientState(
        si->startingState(), recHit->surface(), magneticField.product());
    h_eta_initialstep_seeds->Fill(state.globalMomentum().eta());
  }

  // Loop over all pixel hits
  auto pi = pixelHits->begin();
  auto pe = pixelHits->end();
  for (; pi != pe; ++pi) {
    DetId hitId = pi->detId();
    auto ppi = pi->begin();
    auto ppe = pi->end();
    for (; ppi != ppe; ++ppi) {
      TransientTrackingRecHit::RecHitPointer ttrh = builder->build(&*ppi);
      if (ttrh->isValid()) {
        if (hitId.subdetId() == (int) PixelSubdetector::PixelBarrel) {
          h_hit_pixelbarrel_map->Fill(ttrh->globalPosition().x(), ttrh->globalPosition().y());
          h_hit_pixel_layers->Fill(PXBDetId(hitId).layer());
        } else {
          h_hit_pixel_layers->Fill(3 + PXFDetId(hitId).disk());
        }
      }
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void
DemoTrackAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DemoTrackAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------

void
DemoTrackAnalyzer::beginRun(edm::Run const &run, edm::EventSetup const &setup)
{
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  setup.get<TransientRecHitRecord>().get(builderName, theBuilder);
  builder = theBuilder.product();
}


// ------------ method called when ending the processing of a run  ------------
/*
  void
  DemoTrackAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void
  DemoTrackAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void
  DemoTrackAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoTrackAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoTrackAnalyzer);

// Local Variables:
// truncate-lines: t
// show-trailing-whitespace: t
// End:
