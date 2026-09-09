#ifndef RecoVertex_PrimaryVertexProducer_GNNClusterizerFromAlpaka_h
#define RecoVertex_PrimaryVertexProducer_GNNClusterizerFromAlpaka_h

/**
 * GNNClusterizerFromAlpaka - builds TransientVertex candidates from the vertex-slot GNN outputs
 *
 * No inference happens here: GNNVertexProducerAlpaka runs the model and PrimaryVertexProducer
 * unpacks its output SoA (per track: highest-probability slot and that probability; per slot:
 * predicted z and existence probability). A slot becomes a vertex candidate if its existence
 * probability exceeds existenceThreshold and at least one track points to it; a track joins its
 * slot if its assignment probability is at least trackAssignmentThreshold. The candidate position
 * is the predicted z (the vertex fit refines it) and the track weights are the probabilities.
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include <vector>

namespace vertexgnn {

  class GNNClusterizerFromAlpaka {
  public:
    GNNClusterizerFromAlpaka(const edm::ParameterSet& conf);

    // tracks[i] points to slot trackSlot[i] with probability trackMaxProb[i]; z_hat and p have one
    // entry per slot
    std::vector<TransientVertex> vertices(const std::vector<reco::TransientTrack>& tracks,
                                          const std::vector<int>& trackSlot,
                                          const std::vector<float>& trackMaxProb,
                                          const std::vector<float>& z_hat,
                                          const std::vector<float>& p) const;

    static void fillPSetDescription(edm::ParameterSetDescription& desc);

  private:
    const double existenceThreshold_;
    const double trackAssignmentThreshold_;
    const bool verbose_;
  };

}  // namespace vertexgnn

#endif
