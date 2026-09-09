#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexGNNReco/interface/VertexGNNHostCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "RecoVertex/PrimaryVertexProducer/interface/TrackFilterForPVFinding.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <algorithm>
#include <cmath>

namespace vertexgnn {

  class TrackFeatureProducer : public edm::stream::EDProducer<> {
  public:
    explicit TrackFeatureProducer(const edm::ParameterSet& params)
        : trackToken_(consumes<reco::TrackCollection>(params.getParameter<edm::InputTag>("tracks"))),
          beamSpotToken_(consumes<reco::BeamSpot>(params.getParameter<edm::InputTag>("beamSpot"))),
          mvaToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"))),
          tmtdToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("tmtdSrc"))),
          sigmatmtdToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("sigmatmtdSrc"))),
          pathLengthToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("pathmtd"))),
          tofPiToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("tofPi"))),
          tofKToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("tofK"))),
          tofPToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("tofP"))),
          sigmaTofPiToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("sigmatofpiSrc"))),
          sigmaTofKToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("sigmatofkSrc"))),
          sigmaTofPToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("sigmatofpSrc"))),
          npixBarrelToken_(consumes<edm::ValueMap<int>>(params.getParameter<edm::InputTag>("npixBarrelSrc"))),
          npixEndcapToken_(consumes<edm::ValueMap<int>>(params.getParameter<edm::InputTag>("npixEndcapSrc"))),
          outermostHitPositionToken_(
              consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("outermostHitPositionSrc"))),
          ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
          trackFilter_(params.getParameter<edm::ParameterSet>("TkFilterParameters")) {
      produces<TrackFeaturesHostCollection>();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
      desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
      desc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
      desc.add<edm::InputTag>("tmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
      desc.add<edm::InputTag>("sigmatmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
      desc.add<edm::InputTag>("pathmtd", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
      desc.add<edm::InputTag>("tofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
      desc.add<edm::InputTag>("tofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
      desc.add<edm::InputTag>("tofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
      desc.add<edm::InputTag>("sigmatofpiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
      desc.add<edm::InputTag>("sigmatofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
      desc.add<edm::InputTag>("sigmatofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
      desc.add<edm::InputTag>("npixBarrelSrc", edm::InputTag("trackExtenderWithMTD", "npixBarrel"));
      desc.add<edm::InputTag>("npixEndcapSrc", edm::InputTag("trackExtenderWithMTD", "npixEndcap"));
      desc.add<edm::InputTag>("outermostHitPositionSrc",
                              edm::InputTag("trackExtenderWithMTD:generalTrackOutermostHitPosition"));
      edm::ParameterSetDescription filterDesc;
      TrackFilterForPVFinding::fillPSetDescription(filterDesc);
      desc.add<edm::ParameterSetDescription>("TkFilterParameters", filterDesc)
          ->setComment("must be identical to the TkFilterParameters of the PrimaryVertexProducer using the GNN");
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(edm::Event& event, const edm::EventSetup& eventSetup) override {
      const auto trackHandle = event.getHandle(trackToken_);
      const auto& beamSpot = event.get(beamSpotToken_);
      const auto& mva = event.get(mvaToken_);
      const auto& tmtd = event.get(tmtdToken_);
      const auto& sigmatmtd = event.get(sigmatmtdToken_);
      const auto& pathLength = event.get(pathLengthToken_);
      const auto& tofPi = event.get(tofPiToken_);
      const auto& tofK = event.get(tofKToken_);
      const auto& tofP = event.get(tofPToken_);
      const auto& sigmaTofPi = event.get(sigmaTofPiToken_);
      const auto& sigmaTofK = event.get(sigmaTofKToken_);
      const auto& sigmaTofP = event.get(sigmaTofPToken_);
      const auto& npixBarrel = event.get(npixBarrelToken_);
      const auto& npixEndcap = event.get(npixEndcapToken_);
      const auto& outermostHitPosition = event.get(outermostHitPositionToken_);

      const auto& ttBuilder = eventSetup.getData(ttbToken_);
      const std::vector<reco::TransientTrack> selectedTracks =
          trackFilter_.select(ttBuilder.build(trackHandle, beamSpot));
      const int N = selectedTracks.size();

      auto features = std::make_unique<TrackFeaturesHostCollection>(cms::alpakatools::host(), N);
      auto view = features->view();
      for (int i = 0; i < N; ++i) {
        const reco::Track& track = selectedTracks[i].track();
        const reco::TrackBaseRef ref = selectedTracks[i].trackBaseRef();

        float trackMva = mva[ref];
        float trackPathLength = pathLength[ref];
        if (!std::isfinite(trackMva))
          trackMva = -1.f;
        if (!std::isfinite(trackPathLength))
          trackPathLength = -1.f;

        const float tMtd = tmtd[ref];
        const float sigmaTMtd = sigmatmtd[ref];
        const bool tmtdValid = (sigmaTMtd >= 0.f);
        float t_pi = 0.f, t_k = 0.f, t_p = 0.f;
        float s_pi = 0.2f, s_k = 0.2f, s_p = 0.2f;
        if (tmtdValid && sigmaTofPi[ref] >= 0.f) {
          t_pi = tMtd - tofPi[ref];
          s_pi = combineSigma(sigmaTofPi[ref], sigmaTMtd);
        }
        if (tmtdValid && sigmaTofK[ref] >= 0.f) {
          t_k = tMtd - tofK[ref];
          s_k = combineSigma(sigmaTofK[ref], sigmaTMtd);
        }
        if (tmtdValid && sigmaTofP[ref] >= 0.f) {
          t_p = tMtd - tofP[ref];
          s_p = combineSigma(sigmaTofP[ref], sigmaTMtd);
        }

        view[i].vz() = sanitize(track.vz());
        view[i].dz() = sanitize(track.dzError());
        view[i].pt() = sanitize(track.pt());
        view[i].eta() = sanitize(track.eta());
        view[i].mva() = sanitize(trackMva);
        view[i].pl() = sanitize(trackPathLength);
        view[i].t_pi() = sanitize(t_pi);
        view[i].t_k() = sanitize(t_k);
        view[i].t_p() = sanitize(t_p);
        view[i].s_pi() = sanitize(s_pi);
        view[i].s_k() = sanitize(s_k);
        view[i].s_p() = sanitize(s_p);
        view[i].has_time() = (trackMva >= 0.f) ? 1.f : 0.f;
        view[i].trk_chi2() = sanitize(static_cast<float>(track.chi2()));
        view[i].trk_ndof() = sanitize(static_cast<float>(track.ndof()));
        view[i].trk_validhits() = static_cast<float>(track.numberOfValidHits());
        view[i].npix_barrel() = static_cast<float>(npixBarrel[ref]);
        view[i].npix_endcap() = static_cast<float>(npixEndcap[ref]);
        view[i].outermost_hit_pos() = sanitize(outermostHitPosition[ref]);
      }
      event.put(std::move(features));
    }

  private:
    static float combineSigma(float sigmaTof, float sigmaTMtd) {
      const float a = std::max(sigmaTof, 0.f);
      const float b = std::max(sigmaTMtd, 0.f);
      const float out = std::sqrt(a * a + b * b + 1e-12f);
      return std::isfinite(out) ? out : 1e9f;
    }
    static float sanitize(float v) { return std::isfinite(v) ? v : 1e9f; }

    const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> mvaToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> sigmatmtdToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> tofPiToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> tofKToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> tofPToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> sigmaTofPiToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> sigmaTofKToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> sigmaTofPToken_;
    const edm::EDGetTokenT<edm::ValueMap<int>> npixBarrelToken_;
    const edm::EDGetTokenT<edm::ValueMap<int>> npixEndcapToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> outermostHitPositionToken_;
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
    const TrackFilterForPVFinding trackFilter_;
  };

}  // namespace vertexgnn

DEFINE_FWK_MODULE(vertexgnn::TrackFeatureProducer);
