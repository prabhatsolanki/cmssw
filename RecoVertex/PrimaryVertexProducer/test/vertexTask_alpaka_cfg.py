"""
vertexTask_alpaka_cfg.py - primary-vertex reconstruction with the vertex-slot GNN (Alpaka backend)

Re-runs the vertex reconstruction of a Phase-2 GEN-SIM-RECO file with the vertexSlotGNN process
modifier: 3D (offlinePrimaryVertices), 4D (offlinePrimaryVertices4D) and GNN (offlinePrimaryVerticesGNN)
vertices plus the corresponding tofPID collections, for the 3-way MTD validation
(mtdValidation_3way_cfg.py).

Model, thresholds, input and output names are meant to be overridden by the caller, e.g.
    from vertexTask_alpaka_cfg import *
    process.gnnVertexProducer.model = cms.FileInPath('RecoVertex/PrimaryVertexProducer/data/vertex_slot_v17s29_alpaka.pt')
    process.unsortedOfflinePrimaryVerticesGNN.TkClusParameters.TkDAClusParameters.existenceThreshold = cms.double(0.01)
"""

import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9
from Configuration.ProcessModifiers.vertexSlotGNN_cff import vertexSlotGNN

process = cms.Process('REVTX', Phase2C22I13M9, vertexSlotGNN)

process.load('Configuration.StandardSequences.Services_cff')
process.load('HeterogeneousCore.AlpakaCore.ProcessAcceleratorAlpaka_cfi')
process.load('HeterogeneousCore.CUDACore.ProcessAcceleratorCUDA_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D121Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_20_0_0/RelValTTbar_14TeV/GEN-SIM-RECO/PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU_16Aug26-v6/2590000/02886033-95df-4ef9-9bc6-b148f6aaf487.root',
    ),
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# the vertexSlotGNN modifier adds trackFeatureProducer, gnnVertexProducer, unsortedOfflinePrimaryVerticesGNN,
# offlinePrimaryVerticesGNN and tofPIDGNN to process.vertexreco (RecoVertex.Configuration.RecoVertex_cff)
process.exe = cms.Path(process.firstStepPrimaryVerticesUnsorted * process.vertexreco)

# per-track GNN outputs (slot, probability, PID weights) as histograms
process.TFileService = cms.Service("TFileService", fileName = cms.string("gnn_inspector.root"))
process.gnnInspector = cms.EDAnalyzer("GNNTrackInspector",
    pvModule = cms.InputTag("unsortedOfflinePrimaryVerticesGNN", "", "REVTX"),
    trackSrc = cms.InputTag("generalTracks"),
    printFirstN = cms.uint32(0),
    dropNaNs = cms.bool(True),
)
process.inspect_path = cms.EndPath(process.gnnInspector)

# output: FEVTDEBUGHLT without the RECO-step products superseded in this process
outputCommands = cms.untracked.vstring(list(process.FEVTDEBUGHLTEventContent.outputCommands) + [
    'drop *_offlinePrimaryVertices_*_RECO',
    'drop *_offlinePrimaryVertices4D_*_RECO',
    'drop *_offlinePrimaryVertices4DWithBS_*_RECO',
    'drop *_offlinePrimaryVerticesWithBS_*_RECO',
    'drop *_TriggerResults_*_RECO',
    'drop *_trackTimeValueMapProducer_generalTracksConfigurableFlatResolutionModel_RECO',
    'drop *_trackTimeValueMapProducer_generalTracksConfigurableFlatResolutionModelResolution_RECO',
    'drop *_trackTimeValueMapProducer_generalTracksPerfectResolutionModel_RECO',
    'drop *_trackTimeValueMapProducer_generalTracksPerfectResolutionModelResolution_RECO',
    'drop *_tofPID_probK_RECO',
    'drop *_tofPID_probP_RECO',
    'drop *_tofPID_probPi_RECO',
    'drop *_tofPID_sigmat0_RECO',
    'drop *_tofPID_sigmat0safe_RECO',
    'drop *_tofPID_t0_RECO',
    'drop *_tofPID_t0safe_RECO',
    'drop *_ak4CaloJetsForTrk_*_RECO',
    'drop *_inclusiveSecondaryVertices_*_RECO',
    'drop *_generalV0Candidates_Kshort_RECO',
    'drop *_generalV0Candidates_Lambda_RECO',
])
process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(dataTier = cms.untracked.string('GEN-SIM-RECO'), filterName = cms.untracked.string('')),
    fileName = cms.untracked.string('file:revtx_step3_alpaka.root'),
    outputCommands = outputCommands,
    splitLevel = cms.untracked.int32(0),
)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

from Validation.Performance.TimeMemorySummary import customiseWithTimeMemorySummary
process = customiseWithTimeMemorySummary(process)
