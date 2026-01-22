/**
 * GNNClusterizer.cc - End-to-end Vertex Slot Model based track clustering
 * 
 * This implementation uses a VertexSlotModel ONNX export that directly outputs
 * track-to-vertex assignments and vertex positions, eliminating the need for
 * external clustering algorithms.
 * 
 * ONNX Model Inputs:
 *   - x: (1, N, 13) track features
 *   - row_splits: (1, 2) event boundaries [0, N]
 * 
 * ONNX Model Outputs:
 *   - A: (N, K) soft assignment probabilities
 *   - z_hat: (1, K) predicted vertex z positions [cm]
 *   - t_hat: (1, K) predicted vertex times [ns]
 *   - p: (1, K) slot existence probabilities
 *   - pi: (N, 3) PID mixture weights
 */

#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GNNClusterizer.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryError.h"
#include "TrackingTools/TrajectoryState/interface/PerigeeConversions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include <iostream>
#include <unordered_map>
#include <array>
#include <cmath>
#include <cassert>
#include <limits>
#include <iomanip>
#include <map>
#include <vector>
#include <algorithm>

#define cputime
#ifdef cputime
#include <chrono>
typedef std::chrono::duration<int, std::micro> microseconds_type;
#endif

using namespace std;

// Speed of light in cm/ns (for time-to-spatial conversion)
constexpr float C_LIGHT = 30.0f;

/**
 * Combine timing uncertainties in quadrature
 */
float combine_sigma(float sig_h, float sig_tmtd, float eps = 1e-12f) {
    float h = std::max(sig_h, 0.0f);
    float t = std::max(sig_tmtd, 0.0f);
    float out = std::sqrt(h*h + t*t + eps);
    if (!std::isfinite(out)) {
        out = 1e9f;
    }
    return out;
}

GNNClusterizer::GNNClusterizer(const edm::ParameterSet& conf, const ONNXRuntime* onnxRuntime)
    : onnxRuntime_(onnxRuntime),
      nnVersion_(conf.getParameter<std::string>("nnVersion")),
      verbose_(conf.getUntrackedParameter<bool>("verbose", true)),
      d0CutOff_(conf.getParameter<double>("d0CutOff")),
      vertexSize_(conf.getParameter<double>("vertexSize")),
      existenceThreshold_(conf.getParameter<double>("existenceThreshold")),
      trackAssignmentThreshold_(conf.getParameter<double>("trackAssignmentThreshold")),
      // Note: minTrackCount removed - rely on ndof for filtering
      numSlots_(conf.getParameter<int>("numSlots")) {}

void GNNClusterizer::globalEndJob(const ONNXRuntime* cache) {}

std::unique_ptr<ONNXRuntime> GNNClusterizer::initializeGlobalCache(const edm::ParameterSet& conf) {
  std::string backend_str = conf.getParameter<std::string>("onnxBackend");
  cms::Ort::Backend backend = cms::Ort::Backend::cpu;
  if (backend_str == "CUDA") {
    backend = cms::Ort::Backend::cuda;
  }

  auto session_options = ONNXRuntime::defaultSessionOptions(backend);
  edm::LogInfo("GNNClusterizer") << "Initializing ONNXRuntime with backend: " << backend_str;
  return std::make_unique<ONNXRuntime>(conf.getParameter<edm::FileInPath>("onnxModelPath").fullPath(), &session_options);
}

std::vector<TransientVertex> GNNClusterizer::vertices(const std::vector<reco::TransientTrack>& tracks) const {
  cms::Ort::FloatArrays input_values;
  std::vector<TransientVertex> clusters;
  int N = tracks.size();
  int K = numSlots_;
  int shapeFeatures = 13;  // Input feature dimension (matches model.py)

  // Reset caches
  last_valid_ = false;
  last_n_tracks_ = N;
  last_num_slots_ = K;
  last_assignments_.clear();
  last_z_hat_.clear();
  last_t_hat_.clear();
  last_p_.clear();
  last_pi_.clear();

  if (N == 0) {
    return clusters;
  }

  // =========================================================================
  // PREPARE INPUT FEATURES (matching model.py input_dim=13)
  // =========================================================================
  // Features: [z, dz, pt, eta, mva, pl, t_pi, t_k, t_p, s_pi, s_k, s_p, has_time]
  
  std::vector<float> features;
  features.reserve(N * shapeFeatures);

  for (size_t i = 0; i < tracks.size(); ++i) {
    const reco::TransientTrack& ttrack = tracks[i];
    
    // Basic track properties
    float track_vz = ttrack.track().vz();
    float track_dz = ttrack.track().dzError();
    float track_pt = ttrack.track().pt();
    float track_eta = ttrack.track().eta();
    float track_mva = ttrack.MVAquality();
    float track_pl = ttrack.pathLength();

    // v8p9p2: Ensure missing MVA/PL are -1.0f to match model's missing feature check
    if (!std::isfinite(track_mva)) track_mva = -1.0f;
    if (!std::isfinite(track_pl))  track_pl  = -1.0f;
    
    // Timing information
    float track_t_pi = ttrack.MTDtime() - ttrack.trackTime_pi();
    float track_t_k = ttrack.MTDtime() - ttrack.trackTime_k();
    float track_t_p = ttrack.MTDtime() - ttrack.trackTime_p();
    float track_s_tmtd = ttrack.MTDtimeErr();
    float track_s_pi = ttrack.sigma_time_pi();
    float track_s_k = ttrack.sigma_time_k();
    float track_s_p = ttrack.sigma_time_p();

    // Build timing masks
    bool tmtd_ok = (track_s_tmtd >= 0.0f);
    bool pi_ok = (track_s_pi >= 0.0f) && tmtd_ok;
    bool k_ok = (track_s_k >= 0.0f) && tmtd_ok;
    bool p_ok = (track_s_p >= 0.0f) && tmtd_ok;

    // Apply masks to timing features (default to 0 time, 0.2 sigma if no timing)
    float feat_t_pi = 0.0f, feat_t_k = 0.0f, feat_t_p = 0.0f;
    float feat_s_pi = 0.2f, feat_s_k = 0.2f, feat_s_p = 0.2f;

    if (pi_ok) {
        feat_t_pi = track_t_pi;
        feat_s_pi = combine_sigma(track_s_pi, track_s_tmtd);
    }
    if (k_ok) {
        feat_t_k = track_t_k;
        feat_s_k = combine_sigma(track_s_k, track_s_tmtd);
    }
    if (p_ok) {
        feat_t_p = track_t_p;
        feat_s_p = combine_sigma(track_s_p, track_s_tmtd);
    }

    // has_time flag: 1.0 if MTD timing available, 0.0 otherwise
    float has_time = tmtd_ok ? 1.0f : 0.0f;

    // Push all 13 features
    features.push_back(track_vz);
    features.push_back(track_dz);
    features.push_back(track_pt);
    features.push_back(track_eta);
    features.push_back(track_mva);
    features.push_back(track_pl);
    features.push_back(feat_t_pi);
    features.push_back(feat_t_k);
    features.push_back(feat_t_p);
    features.push_back(feat_s_pi);
    features.push_back(feat_s_k);
    features.push_back(feat_s_p);
    features.push_back(has_time);
  }

  // Replace NaN/Inf with large number
  for (auto& v : features) {
    if (!std::isfinite(v)) v = 1e9f;
  }

  // =========================================================================
  // RUN ONNX INFERENCE
  // =========================================================================
  // New export: single input "x" only, row_splits are handled internally
  std::vector<std::string> input_names = {"x"};

  input_values.emplace_back(features);

  // Input dimensions: x is (1, N, 13) - batch dim required by CMSSW wrapper
  std::vector<std::vector<long int>> input_dims;
  input_dims.push_back({1, N, shapeFeatures});

#ifdef cputime
  std::chrono::duration<int, std::micro> tcpu_inference(0), tcpu_out(0);
  auto start_inference = std::chrono::high_resolution_clock::now();
#endif

  auto output_values = onnxRuntime_->run(input_names, input_values, input_dims);

#ifdef cputime
  auto stop_inference = std::chrono::high_resolution_clock::now();
  tcpu_inference = std::chrono::duration_cast<std::chrono::microseconds>(stop_inference - start_inference);
  edm::LogInfo("GNNClusterizer") << "###TIME inference " << tcpu_inference.count() << " us";
#endif

  // =========================================================================
  // PARSE ONNX OUTPUTS
  // =========================================================================
  // Expected outputs from VertexSlotModel ONNX export:
  //   0: A      - (N, K) assignment matrix
  //   1: z_hat  - (1, K) or (K,) vertex z positions
  //   2: t_hat  - (1, K) or (K,) vertex times
  //   3: p      - (1, K) or (K,) existence probabilities
  //   4: pi     - (N, 3) PID mixture weights

  std::vector<float>& A_flat = output_values[0];      // N*K
  std::vector<float>& z_hat_flat = output_values[1];  // K or 1*K
  std::vector<float>& t_hat_flat = output_values[2];  // K or 1*K
  std::vector<float>& p_flat = output_values[3];      // K or 1*K
  std::vector<float>& pi_flat = output_values[4];     // N*3

  // Cache raw outputs
  last_assignments_ = A_flat;
  last_z_hat_ = z_hat_flat;
  last_t_hat_ = t_hat_flat;
  last_p_ = p_flat;
  last_pi_ = pi_flat;
  last_valid_ = true;

  // Infer K from z_hat size (handles both (K,) and (1,K) shapes)
  K = static_cast<int>(z_hat_flat.size());
  last_num_slots_ = K;

  if (verbose_) {
    edm::LogInfo("GNNClusterizer") << "ONNX outputs: N=" << N << ", K=" << K 
                                   << ", A_size=" << A_flat.size()
                                   << ", p_size=" << p_flat.size();
  }

#ifdef cputime
  auto start_out = std::chrono::high_resolution_clock::now();
#endif

  // =========================================================================
  // PROCESS ASSIGNMENTS - End-to-end clustering
  // =========================================================================
  
  // Step 1: For each track, find the slot with maximum assignment probability
  // Assign to best slot regardless of probability (no threshold)
  std::vector<int> track_to_slot(N, -1);
  std::vector<float> track_max_prob(N, 0.0f);
  
  for (int i = 0; i < N; ++i) {
    float max_prob = -1.0f;
    int best_slot = -1;
    for (int k = 0; k < K; ++k) {
      float prob = A_flat[i * K + k];
      if (prob > max_prob) {
        max_prob = prob;
        best_slot = k;
      }
    }
    // Assign track to best slot (no threshold - use all assignments)
    track_to_slot[i] = best_slot;
    track_max_prob[i] = max_prob;
  }

  // Step 2: Count tracks per slot
  std::vector<int> slot_track_count(K, 0);
  for (int i = 0; i < N; ++i) {
    int slot = track_to_slot[i];
    if (slot >= 0 && slot < K) {
      slot_track_count[slot]++;
    }
  }

  // Step 3: Identify active slots (existence probability > threshold)
  // Note: minTrackCount filtering removed - rely on ndof for downstream filtering
  std::vector<bool> slot_active(K, false);
  int num_active = 0;
  for (int k = 0; k < K; ++k) {
    if (p_flat[k] > existenceThreshold_ && slot_track_count[k] >= 1) {
      slot_active[k] = true;
      num_active++;
    }
  }

  if (verbose_) {
    edm::LogInfo("GNNClusterizer") << "Active slots: " << num_active 
                                   << " (threshold p>" << existenceThreshold_ << ")";
  }

  // Step 4: Build vertex clusters from active slots (Pure GNN approach)
  // Use 4D error matrix for position (x, y, z) and time (t)
  // Note: sigma_xy and sigma_z now computed dynamically from Hessian inverse  
  
  for (int k = 0; k < K; ++k) {
    if (!slot_active[k]) continue;

    // Collect tracks assigned to this slot with original indices for pi_flat access
    std::vector<reco::TransientTrack> cluster_tracks;
    std::vector<float> cluster_weights;
    std::vector<int> cluster_orig_indices;  // Original track indices for pi_flat access
    for (int i = 0; i < N; ++i) {
      if (track_to_slot[i] == k) {
        cluster_tracks.push_back(tracks[i]);
        cluster_weights.push_back(track_max_prob[i]);
        cluster_orig_indices.push_back(i);  // Store original index
      }
    }

    if (cluster_tracks.empty()) continue;

    // Use model's predicted z and t positions directly
    double z_vertex = static_cast<double>(z_hat_flat[k]);
    double t_vertex = static_cast<double>(t_hat_flat[k]);
    
    // Compute vertex time error from track timing uncertainties using GNN PID weights
    // Similar to VertexTimeAlgorithmFromTracksPID but using GNN's learned PID probabilities
    // sigma_t_vertex = sqrt(sum(w²*σ²)) / sum(w) where w includes PID weights
    double wsum = 0.0;
    double w2sum = 0.0;
    double z_wsum = 0.0;      // For z error
    double z_w2_s2_sum = 0.0; // For z error
    int n_tracks_with_time = 0;
    
    for (size_t ii = 0; ii < cluster_tracks.size(); ++ii) {
      const reco::TransientTrack& ttrack = cluster_tracks[ii];
      int orig_idx = cluster_orig_indices[ii];  // Use stored original index
      
      float assignment_prob = cluster_weights[ii];
      
      // Get PID weights from GNN output (pi_flat is N x 3)
      float pi_pion = (orig_idx >= 0 && orig_idx < N) ? pi_flat[orig_idx * 3 + 0] : 0.7f;
      float pi_kaon = (orig_idx >= 0 && orig_idx < N) ? pi_flat[orig_idx * 3 + 1] : 0.2f;
      float pi_proton = (orig_idx >= 0 && orig_idx < N) ? pi_flat[orig_idx * 3 + 2] : 0.1f;
      
      // Get PID time uncertainties (combined MTD + TOF errors)
      float sigma_t_tmtd = ttrack.MTDtimeErr();
      float sigma_t_pi = ttrack.sigma_time_pi();
      float sigma_t_k = ttrack.sigma_time_k();
      float sigma_t_p = ttrack.sigma_time_p();
      
      // Combine uncertainties (like classical approach lines 118-126)
      // sigma_hyp = sqrt(MTD_err² + TOF_sigma²)
      bool has_time = (sigma_t_tmtd > 0 && std::isfinite(sigma_t_tmtd));
      if (has_time && sigma_t_pi > 0) {
        float s_pi = std::sqrt(sigma_t_tmtd*sigma_t_tmtd + sigma_t_pi*sigma_t_pi);
        float s_k = (sigma_t_k > 0) ? std::sqrt(sigma_t_tmtd*sigma_t_tmtd + sigma_t_k*sigma_t_k) : s_pi;
        float s_p = (sigma_t_p > 0) ? std::sqrt(sigma_t_tmtd*sigma_t_tmtd + sigma_t_p*sigma_t_p) : s_pi;
        
        // Weight by assignment prob, PID weight, and inverse variance
        // Like classical line 176: w = wt * trkWeight / (dt[j] * dt[j])
        float w_pi = assignment_prob * pi_pion / (s_pi * s_pi);
        float w_k = assignment_prob * pi_kaon / (s_k * s_k);
        float w_p = assignment_prob * pi_proton / (s_p * s_p);
        
        // Per-track sums (like classical lines 173-180)
        float wsum_trk = w_pi + w_k + w_p;
        float wsum_sigma_trk = w_pi * s_pi + w_k * s_k + w_p * s_p;
        
        wsum += wsum_trk;
        // Classical w2sum: square of per-track sum (line 183)
        w2sum += wsum_sigma_trk * wsum_sigma_trk;
        n_tracks_with_time++;
      }
      
      // Also compute z error from track dz
      float dz_track = ttrack.track().dzError();
      if (dz_track > 0 && std::isfinite(dz_track)) {
        float w_z = assignment_prob / (dz_track * dz_track);
        z_wsum += w_z;
        z_w2_s2_sum += (w_z * dz_track) * (w_z * dz_track);
      }
    }  // End of track loop for time/z errors
    
    // =========================================================================
    // COMPUTE VERTEX PROPERTIES EXACTLY LIKE CLASSICAL FITTER
    // (AdaptiveChisquarePrimaryVertexFitter)
    // =========================================================================
    
    // Vertex time error: weighted combination (like classical line 197)
    double sigma_t_vertex = (wsum > 0) ? std::sqrt(w2sum) / wsum : 0.030;
    
    // -------------------------------------------------------------------------
    // Step 1: Compute vertex x,y position and Hessian matrix for covariance
    // (Classical uses beamspot + iterative fit; we use weighted track positions)
    // -------------------------------------------------------------------------
    double x_sum = 0.0, y_sum = 0.0, w_pos_sum = 0.0;
    
    // Hessian matrix S (3x3) for position covariance - same as classical lines 271-278
    double S00 = 0, S01 = 0, S02 = 0, S11 = 0, S12_S = 0, S22 = 0;
    
    // Track projection vectors for chi2 calculation - stored for second pass
    std::vector<std::array<double, 10>> track_perigee_info;  // Store perigee data
    track_perigee_info.reserve(cluster_tracks.size());
    
    for (size_t i = 0; i < cluster_tracks.size(); ++i) {
      const reco::TransientTrack& trk = cluster_tracks[i];
      float track_weight = cluster_weights[i];
      
      try {
        const auto tspca = trk.stateAtBeamLine().trackStateAtPCA();
        const auto tspca_pe = PerigeeConversions::ftsToPerigeeError(tspca);
        const auto momentum = tspca.momentum();
        const double cos_phi = momentum.x() / momentum.perp();
        const double sin_phi = momentum.y() / momentum.perp();
        const double tan_lambda = momentum.z() / momentum.perp();
        
        // Perigee covariance matrix (indices 3,4 = d0, z0)
        double cov11 = tspca_pe.covarianceMatrix()(3, 3);
        double cov22 = tspca_pe.covarianceMatrix()(4, 4);
        double cov12 = tspca_pe.covarianceMatrix()(3, 4);
        
        double DetV = cov11 * cov22 - cov12 * cov12;
        if (std::fabs(DetV) > 1.e-16) {
          double invS11 = cov22 / DetV;
          double invS22 = cov11 / DetV;
          double invS12 = -cov12 / DetV;
          
          // H matrices (same as classical lines 70-76)
          double H1[3] = {-sin_phi, cos_phi, 0.0};
          double H2[3] = {cos_phi * tan_lambda, sin_phi * tan_lambda, -1.0};
          
          // Compute C matrix for Hessian: C(l,k) = H1[l]*SH1[k] + H2[l]*SH2[k]
          // (same as classical lines 78-84)
          for (int kk = 0; kk < 3; kk++) {
            double SH1k = invS11 * H1[kk] + invS12 * H2[kk];
            double SH2k = invS12 * H1[kk] + invS22 * H2[kk];
            for (int ll = 0; ll <= kk; ll++) {
              double C_lk = H1[ll] * SH1k + H2[ll] * SH2k;
              // Accumulate weighted Hessian: S += w * C (same as line 278)
              if (ll == 0 && kk == 0) S00 += track_weight * C_lk;
              else if (ll == 0 && kk == 1) S01 += track_weight * C_lk;
              else if (ll == 0 && kk == 2) S02 += track_weight * C_lk;
              else if (ll == 1 && kk == 1) S11 += track_weight * C_lk;
              else if (ll == 1 && kk == 2) S12_S += track_weight * C_lk;
              else if (ll == 2 && kk == 2) S22 += track_weight * C_lk;
            }
          }
          
          // Store perigee info for chi2 calculation
          double b1 = tspca.position().x() * sin_phi - tspca.position().y() * cos_phi;
          double b2 = tspca.position().z() - (tspca.position().x() * cos_phi + tspca.position().y() * sin_phi) * tan_lambda;
          track_perigee_info.push_back({{invS11, invS22, invS12, b1, b2, H1[0], H1[1], H2[0], H2[1], static_cast<double>(track_weight)}});
          
          // Weighted position for x,y
          x_sum += track_weight * tspca.position().x();
          y_sum += track_weight * tspca.position().y();
          w_pos_sum += track_weight;
        } else {
          track_perigee_info.push_back({{0,0,0,0,0,0,0,0,0,0}});  // Invalid track
        }
      } catch (...) {
        track_perigee_info.push_back({{0,0,0,0,0,0,0,0,0,0}});  // Invalid track
      }
    }
    
    // Compute vertex x,y position (weighted average of track positions at beamline)
    double x_vertex = (w_pos_sum > 0) ? x_sum / w_pos_sum : 0.0;
    double y_vertex = (w_pos_sum > 0) ? y_sum / w_pos_sum : 0.0;
    
    // -------------------------------------------------------------------------
    // Step 2: Invert Hessian for covariance (same as classical lines 298-303)
    // S is symmetric 3x3 matrix
    // -------------------------------------------------------------------------
    double cov_xx = 0.0001, cov_yy = 0.0001, cov_zz = 0.0001;  // Defaults
    double cov_xy = 0, cov_xz = 0, cov_yz = 0;
    
    // Invert 3x3 symmetric matrix S -> Covariance
    double det_S = S00*(S11*S22 - S12_S*S12_S) - S01*(S01*S22 - S12_S*S02) + S02*(S01*S12_S - S11*S02);
    if (std::fabs(det_S) > 1.e-20) {
      // Inverse using cofactors
      cov_xx = (S11*S22 - S12_S*S12_S) / det_S;
      cov_yy = (S00*S22 - S02*S02) / det_S;
      cov_zz = (S00*S11 - S01*S01) / det_S;
      cov_xy = -(S01*S22 - S02*S12_S) / det_S;
      cov_xz = (S01*S12_S - S02*S11) / det_S;
      cov_yz = -(S00*S12_S - S01*S02) / det_S;
    }
    
    // 4D error matrix: (xx, xy, yy, xz, yz, zz, xt, yt, zt, tt)
    GlobalError error4D(cov_xx, cov_xy, cov_yy, 
                        cov_xz, cov_yz, cov_zz,
                        0., 0., 0., sigma_t_vertex*sigma_t_vertex);
    
    // -------------------------------------------------------------------------
    // Step 3: Compute chi2 and ndof (same as classical lines 352-375)
    // -------------------------------------------------------------------------
    float vtx_chi2 = 0.0f;
    float vtx_ndof = -3.0f;
    
    for (size_t i = 0; i < cluster_tracks.size(); ++i) {
      float track_weight = cluster_weights[i];
      vtx_ndof += 2.0f * track_weight;
      
      const auto& pi = track_perigee_info[i];
      if (pi[9] > 0) {  // Valid track
        double invS11 = pi[0], invS22 = pi[1], invS12 = pi[2];
        double b1 = pi[3], b2 = pi[4];
        double H1_0 = pi[5], H1_1 = pi[6];
        double H2_0 = pi[7], H2_1 = pi[8];
        
        // Residuals at fitted vertex position (same as lines 25-26)
        double F1 = b1 + x_vertex * H1_0 + y_vertex * H1_1;
        double F2 = b2 + x_vertex * H2_0 + y_vertex * H2_1 - z_vertex;
        
        // Chi2 contribution (same as line 27)
        double track_chi2 = F1 * F1 * invS11 + F2 * F2 * invS22 + 2. * F1 * F2 * invS12;
        vtx_chi2 += track_weight * track_chi2;
      }
    }

    // Create vertex with computed x,y position and Hessian-based covariance
    GlobalPoint pos(x_vertex, y_vertex, z_vertex);
    
    // Use time-aware constructor with classical-like chi2 and ndof
    TransientVertex vertex(pos, t_vertex, error4D, cluster_tracks, vtx_chi2, vtx_ndof);
    
    // Set track weights from GNN assignment probabilities
    TransientVertex::TransientTrackToFloatMap weightMap;
    for (size_t i = 0; i < cluster_tracks.size(); ++i) {
      weightMap[cluster_tracks[i]] = cluster_weights[i];
    }
    vertex.weightMap(weightMap);
    
    clusters.push_back(vertex);

    if (verbose_) {
      edm::LogInfo("GNNClusterizer") << "  Vertex k=" << k 
                                     << ": z=" << z_vertex
                                     << ", t=" << t_vertex
                                     << ", p=" << p_flat[k]
                                     << ", chi2=" << vtx_chi2
                                     << ", ndof=" << vtx_ndof
                                     << ", ntracks=" << cluster_tracks.size();
    }
  }

#ifdef cputime
  auto stop_out = std::chrono::high_resolution_clock::now();
  tcpu_out = std::chrono::duration_cast<std::chrono::microseconds>(stop_out - start_out);
  edm::LogInfo("GNNClusterizer") << "###TIME clustering " << tcpu_out.count() << " us";
#endif

  if (verbose_) {
    edm::LogInfo("GNNClusterizer") << "Final vertex count: " << clusters.size();
  }

  return clusters;
}

std::vector<std::vector<reco::TransientTrack>> GNNClusterizer::clusterize(
    const std::vector<reco::TransientTrack>& tracks) const {
  return std::vector<std::vector<reco::TransientTrack>>();
}

void GNNClusterizer::fillPSetDescription(edm::ParameterSetDescription& desc) {
  DAClusterizerInZ_vect::fillPSetDescription(desc);
  
  // ONNX model configuration
  desc.add<edm::FileInPath>("onnxModelPath", 
      edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/vertex_slot_model.onnx"))
      ->setComment("Path to the VertexSlotModel ONNX export");
  desc.add<std::string>("nnVersion", "vertex_slot_v8p9")
      ->setComment("Model version tag");
  
  // VertexSlotModel-specific parameters
  desc.add<double>("existenceThreshold", 0.5)
      ->setComment("Slot existence probability threshold (0-1)");
  desc.add<double>("trackAssignmentThreshold", 0.4)
      ->setComment("Minimum track assignment probability - tracks below this are not assigned");
  // Note: minTrackCount removed - rely on ndof for downstream filtering
  desc.add<int>("numSlots", 200)
      ->setComment("Number of slots K in the model (should match ONNX export)");
  desc.add<std::string>("onnxBackend", "CPU")
      ->setComment("ONNX execution backend: CPU, CUDA, etc.");

  // General vertex parameters (retained from base)
  desc.addUntracked<bool>("verbose", false);
}
