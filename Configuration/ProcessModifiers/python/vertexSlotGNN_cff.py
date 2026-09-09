import FWCore.ParameterSet.Config as cms

# Modifier to run the vertex-slot GNN primary-vertex reconstruction
# (Phase-2 with MTD timing layer): produces offlinePrimaryVerticesGNN and tofPIDGNN
# in addition to the standard 3D and 4D primary vertices

vertexSlotGNN = cms.Modifier()
