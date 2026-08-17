// Pixel hit efficiency from the stored track hit pattern.
//
// Reads reco::Track and reco::Vertex only. Every layer a track should have
// crossed, and whether it found a hit there, was already decided at
// reconstruction time and is kept in reco::HitPattern so this runs on AOD
// and needs no refitting

#include <algorithm>
#include <cmath>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TTree.h"

namespace {

  // Regions the hit patterns store: 0-3 = BPIX L1-L4,
  // 4-6 = FPIX D1-D3, -1 = anything else.
  constexpr int kNRegion = 7;

  int pixelRegion(uint16_t pattern) {
    if (reco::HitPattern::pixelBarrelHitFilter(pattern)) {
      const uint32_t layer = reco::HitPattern::getLayer(pattern);
      return (layer >= 1 && layer <= 4) ? static_cast<int>(layer) - 1 : -1;
    }
    if (reco::HitPattern::pixelEndcapHitFilter(pattern)) {
      const uint32_t disk = reco::HitPattern::getLayer(pattern);
      return (disk >= 1 && disk <= 3) ? 3 + static_cast<int>(disk) : -1;
    }
    return -1;
  }

}  // namespace

class SiPixelTrackHitPatternTree : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SiPixelTrackHitPatternTree(const edm::ParameterSet&);
  ~SiPixelTrackHitPatternTree() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  bool isGoodVertex(const reco::Vertex&) const;

  // ---- configuration ----
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

  const double trackPtCut_; // track quality cuts
  const int trackNStripCut_;
  const bool requireHighPurity_;

  const double pvMinNdof_; // pv quality cuts
  const double pvMaxAbsZ_;
  const double pvMaxRho_;
  const int pvMinNTrack_;

  const double d0Max_; //track max distances to leading PV
  const double dzMax_;

  // ---- counters, written to the jobInfo tree for bookkeepign ----
  ULong64_t nEvents_ = 0;
  ULong64_t nEventsWithPV_ = 0;
  ULong64_t nTracksRead_ = 0;
  ULong64_t nTracksQuality_ = 0;
  ULong64_t nTracksSelected_ = 0;
  ULong64_t nNoL1Entry_ = 0;

  // Which hit pattern category the L1 entries came from. MISSING_OUTER is not
  // a mistake: for an outside-in trajectory setSecondHitPattern labels the
  // inner layers missing_outer, so L1 lands can be MISSING_OUTER
  ULong64_t nL1FromTrack_ = 0;
  ULong64_t nL1FromInner_ = 0;
  ULong64_t nL1FromOuter_ = 0;

  // ---- trees ----
  TTree* tree_ = nullptr;
  TTree* info_ = nullptr;

  UInt_t b_run_, b_ls_;
  ULong64_t b_event_;
  Int_t b_bx_;

  Int_t b_npv_, b_pv_ntrk_;
  Float_t b_pv_z_;

  Float_t b_trk_pt_, b_trk_eta_, b_trk_phi_, b_trk_d0_, b_trk_dz_;
  
  // reco::TrackBase::TrackAlgorithm. Stored because muon-seeded
  // outside-in steps (13, 14) are the ones whose inner layers arrive via
  // MISSING_OUTER_HITS, possibly affected by https://github.com/cms-sw/cmssw/pull/51478
  Int_t b_trk_algo_;

  // Per region, counted rather than reduced to one status code.
  // In unlikely case of track crossing a layer multiple times keeps the info
  // Also makes it possible to check if a layer was not crossed at all (0,0,0)
  UChar_t b_nvalid_[kNRegion], b_nmissing_[kNRegion], b_ninactive_[kNRegion];
};

SiPixelTrackHitPatternTree::SiPixelTrackHitPatternTree(const edm::ParameterSet& iConfig)
    : trackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
      vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
      trackPtCut_(iConfig.getParameter<double>("trackPtCut")),
      trackNStripCut_(iConfig.getParameter<int>("trackNStripCut")),
      requireHighPurity_(iConfig.getParameter<bool>("requireHighPurity")),
      pvMinNdof_(iConfig.getParameter<double>("pvMinNdof")),
      pvMaxAbsZ_(iConfig.getParameter<double>("pvMaxAbsZ")),
      pvMaxRho_(iConfig.getParameter<double>("pvMaxRho")),
      pvMinNTrack_(iConfig.getParameter<int>("pvMinNTrack")),
      d0Max_(iConfig.getParameter<double>("d0Max")),
      dzMax_(iConfig.getParameter<double>("dzMax")) {
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tracks", "pixel hit pattern, one row per selected track");

  tree_->Branch("run", &b_run_, "run/i");
  tree_->Branch("ls", &b_ls_, "ls/i");
  tree_->Branch("event", &b_event_, "event/l");
  tree_->Branch("bx", &b_bx_, "bx/I");

  tree_->Branch("npv", &b_npv_, "npv/I");
  tree_->Branch("pv_ntrk", &b_pv_ntrk_, "pv_ntrk/I");
  tree_->Branch("pv_z", &b_pv_z_, "pv_z/F");

  tree_->Branch("trk_pt", &b_trk_pt_, "trk_pt/F");
  tree_->Branch("trk_eta", &b_trk_eta_, "trk_eta/F");
  tree_->Branch("trk_phi", &b_trk_phi_, "trk_phi/F");
  // Measured to the leading primary vertex, vertices->at(0)
  tree_->Branch("trk_d0", &b_trk_d0_, "trk_d0/F");
  tree_->Branch("trk_dz", &b_trk_dz_, "trk_dz/F");
  tree_->Branch("trk_algo", &b_trk_algo_, "trk_algo/I");

  // Region order: BPIX L1-L4, then FPIX D1-D3.
  tree_->Branch("nvalid", b_nvalid_, "nvalid[7]/b");
  tree_->Branch("nmissing", b_nmissing_, "nmissing[7]/b");
  tree_->Branch("ninactive", b_ninactive_, "ninactive[7]/b");

  // One entry, filled at endJob.
  info_ = fs->make<TTree>("jobInfo", "cuts applied and tracks seen");
  info_->Branch("trackPtCut", const_cast<double*>(&trackPtCut_), "trackPtCut/D");
  info_->Branch("trackNStripCut", const_cast<int*>(&trackNStripCut_), "trackNStripCut/I");
  info_->Branch("requireHighPurity", const_cast<bool*>(&requireHighPurity_), "requireHighPurity/O");
  info_->Branch("pvMinNdof", const_cast<double*>(&pvMinNdof_), "pvMinNdof/D");
  info_->Branch("pvMaxAbsZ", const_cast<double*>(&pvMaxAbsZ_), "pvMaxAbsZ/D");
  info_->Branch("pvMaxRho", const_cast<double*>(&pvMaxRho_), "pvMaxRho/D");
  info_->Branch("pvMinNTrack", const_cast<int*>(&pvMinNTrack_), "pvMinNTrack/I");
  info_->Branch("d0Max", const_cast<double*>(&d0Max_), "d0Max/D");
  info_->Branch("dzMax", const_cast<double*>(&dzMax_), "dzMax/D");
  info_->Branch("nEvents", &nEvents_, "nEvents/l");
  info_->Branch("nEventsWithPV", &nEventsWithPV_, "nEventsWithPV/l");
  info_->Branch("nTracksRead", &nTracksRead_, "nTracksRead/l");
  info_->Branch("nTracksQuality", &nTracksQuality_, "nTracksQuality/l");
  info_->Branch("nTracksSelected", &nTracksSelected_, "nTracksSelected/l");
  info_->Branch("nNoL1Entry", &nNoL1Entry_, "nNoL1Entry/l");
  info_->Branch("nL1FromTrack", &nL1FromTrack_, "nL1FromTrack/l");
  info_->Branch("nL1FromInner", &nL1FromInner_, "nL1FromInner/l");
  info_->Branch("nL1FromOuter", &nL1FromOuter_, "nL1FromOuter/l");
}

void SiPixelTrackHitPatternTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("primaryVertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<double>("trackPtCut", 1.0);
  desc.add<int>("trackNStripCut", 10);
  desc.add<bool>("requireHighPurity", true);

  // To be checked
  desc.add<double>("pvMinNdof", 4.0);
  desc.add<double>("pvMaxAbsZ", 24.0);
  desc.add<double>("pvMaxRho", 2.0);
  desc.add<int>("pvMinNTrack", 10);

  desc.add<double>("d0Max", 0.1);
  desc.add<double>("dzMax", 0.5);
  descriptions.addWithDefaultLabel(desc);
}

bool SiPixelTrackHitPatternTree::isGoodVertex(const reco::Vertex& v) const {
  return !v.isFake() && v.isValid() && v.ndof() >= pvMinNdof_ && std::abs(v.z()) < pvMaxAbsZ_ &&
         v.position().rho() < pvMaxRho_;
}

void SiPixelTrackHitPatternTree::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  ++nEvents_;

  const edm::Handle<reco::TrackCollection> tracks = iEvent.getHandle(trackToken_);
  const edm::Handle<reco::VertexCollection> vertices = iEvent.getHandle(vertexToken_);
  if (!tracks.isValid() || !vertices.isValid() || vertices->empty())
    return;

  const reco::Vertex& pv = vertices->front();
  if (!isGoodVertex(pv) || static_cast<int>(pv.tracksSize()) <= pvMinNTrack_)
    return;
  ++nEventsWithPV_;

  b_npv_ = 0;
  for (const auto& v : *vertices)
    if (isGoodVertex(v))
      ++b_npv_;

  b_run_ = iEvent.id().run();
  b_ls_ = iEvent.luminosityBlock();
  b_event_ = iEvent.id().event();
  b_bx_ = iEvent.bunchCrossing();
  b_pv_ntrk_ = static_cast<int>(pv.tracksSize());
  b_pv_z_ = pv.z();

  for (const reco::Track& track : *tracks) {
    ++nTracksRead_;

    // ---- track-level selection ----
    if (requireHighPurity_ && !track.quality(reco::TrackBase::highPurity))
      continue;
    if (track.pt() <= trackPtCut_)
      continue;
    const int nstrip = track.hitPattern().numberOfValidStripHits();
    if (nstrip <= trackNStripCut_)
      continue;
    ++nTracksQuality_;

    const double d0 = -1.0 * track.dxy(pv.position());
    const double dz = track.dz(pv.position());
    if (d0Max_ > 0 && std::abs(d0) >= d0Max_)
      continue;
    if (dzMax_ > 0 && std::abs(dz) >= dzMax_)
      continue;

    b_trk_pt_ = track.pt();
    b_trk_eta_ = track.eta();
    b_trk_phi_ = track.phi();
    b_trk_d0_ = d0;
    b_trk_dz_ = dz;
    b_trk_algo_ = static_cast<int>(track.algo());

    std::fill(std::begin(b_nvalid_), std::end(b_nvalid_), 0);
    std::fill(std::begin(b_nmissing_), std::end(b_nmissing_), 0);
    std::fill(std::begin(b_ninactive_), std::end(b_ninactive_), 0);

    // Loop over all three categories and classify by decoded pattern:
    // an invalid innermost hit is stripped from the trajectory, so L1 usually
    // arrives via a MISSING_* (secondary) pattern rather than TRACK_HITS.
    const reco::HitPattern& hp = track.hitPattern();
    for (auto category : {reco::HitPattern::TRACK_HITS,
                          reco::HitPattern::MISSING_INNER_HITS,
                          reco::HitPattern::MISSING_OUTER_HITS}) {
      for (int i = 0; i < hp.numberOfAllHits(category); ++i) {
        const uint16_t pattern = hp.getHitPattern(category, i);
        const int region = pixelRegion(pattern);
        if (region < 0)
          continue;

        switch (reco::HitPattern::getHitType(pattern)) {
          case reco::HitPattern::VALID:
            ++b_nvalid_[region];
            break;
          case reco::HitPattern::MISSING:
            ++b_nmissing_[region];
            break;
          case reco::HitPattern::INACTIVE:
            ++b_ninactive_[region];
            break;
          default:
            break;
        }

        if (region == 0) {
          if (category == reco::HitPattern::TRACK_HITS)
            ++nL1FromTrack_;
          else if (category == reco::HitPattern::MISSING_INNER_HITS)
            ++nL1FromInner_;
          else
            ++nL1FromOuter_;
        }
      }
    }

    if (b_nvalid_[0] == 0 && b_nmissing_[0] == 0 && b_ninactive_[0] == 0)
      ++nNoL1Entry_;

    ++nTracksSelected_;
    tree_->Fill();
  }
}

void SiPixelTrackHitPatternTree::endJob() {
  info_->Fill();

  edm::LogPrint("SiPixelTrackHitPatternTree")
      << "SiPixelTrackHitPatternTree summary:"
      << "\n  events                    " << nEvents_
      << "\n  ... with a good PV        " << nEventsWithPV_
      << "\n  tracks read               " << nTracksRead_
      << "\n  ... passing track quality " << nTracksQuality_
      << "\n  ... passing d0/dz         " << nTracksSelected_
      << "\n  ... with no BPIX L1 entry " << nNoL1Entry_
      << "\n  BPIX L1 entries from TRACK_HITS / MISSING_INNER / MISSING_OUTER  " << nL1FromTrack_ << " / "
      << nL1FromInner_ << " / " << nL1FromOuter_;
}

DEFINE_FWK_MODULE(SiPixelTrackHitPatternTree);
