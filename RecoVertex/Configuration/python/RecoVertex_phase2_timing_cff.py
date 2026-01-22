import FWCore.ParameterSet.Config as cms
from RecoVertex.Configuration.RecoVertex_cff import unsortedOfflinePrimaryVertices, trackWithVertexRefSelector, trackRefsForJets, sortedPrimaryVertices, offlinePrimaryVertices, offlinePrimaryVerticesWithBS,vertexrecoTask

unsortedOfflinePrimaryVertices4D = unsortedOfflinePrimaryVertices.clone(
    TkClusParameters = cms.PSet(algorithm = cms.string("DA2D_vect"),
        TkDAClusParameters = cms.PSet(
            Tmin = cms.double(4.0),
            Tpurge = cms.double(4.0),
            Tstop = cms.double(2.0),
        )
    ),
    TrackTimesLabel = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModel"),
    TrackTimeResosLabel = cms.InputTag("trackTimeValueMapProducer","generalTracksConfigurableFlatResolutionModelResolution"),
    trackMTDTimeQualityVMapTag = cms.InputTag("mtdTrackQualityMVA:mtdQualMVA"),
    useMVACut = cms.bool(False),
    minTrackTimeQuality = cms.double(0.8),
    vertexCollections = {0: dict(vertexTimeParameters = cms.PSet( algorithm = cms.string('fromTracksPID'))),
                         1: dict(vertexTimeParameters = cms.PSet( algorithm = cms.string('fromTracksPID')))}
    )
trackWithVertexRefSelectorBeforeSorting4D = trackWithVertexRefSelector.clone(
    vertexTag = "unsortedOfflinePrimaryVertices4D",
    ptMax = 9e99,
    ptErrorCut = 9e99
)
trackRefsForJetsBeforeSorting4D = trackRefsForJets.clone(
    src = "trackWithVertexRefSelectorBeforeSorting4D"
)
offlinePrimaryVertices4D = sortedPrimaryVertices.clone(
    vertices = "unsortedOfflinePrimaryVertices4D",
    particles = "trackRefsForJetsBeforeSorting4D",
    trackTimeTag = "trackTimeValueMapProducer:generalTracksConfigurableFlatResolutionModel",
    trackTimeResoTag = "trackTimeValueMapProducer:generalTracksConfigurableFlatResolutionModelResolution",
    assignment = dict(useTiming = True)
)
offlinePrimaryVertices4DWithBS = offlinePrimaryVertices4D.clone(
    vertices = "unsortedOfflinePrimaryVertices4D:WithBS"
)

unsortedOfflinePrimaryVertices4DwithPID = unsortedOfflinePrimaryVertices4D.clone(
    TrackTimesLabel = "tofPID4DnoPID:t0safe",
    TrackTimeResosLabel = "tofPID4DnoPID:sigmat0safe"
)
trackWithVertexRefSelectorBeforeSorting4DwithPID = trackWithVertexRefSelector.clone(
    vertexTag = "unsortedOfflinePrimaryVertices4DwithPID",
    ptMax = 9e99,
    ptErrorCut = 9e99
)
trackRefsForJetsBeforeSorting4DwithPID = trackRefsForJets.clone(
    src = "trackWithVertexRefSelectorBeforeSorting4DwithPID"
)
offlinePrimaryVertices4DwithPID=offlinePrimaryVertices4D.clone(
    vertices = "unsortedOfflinePrimaryVertices4DwithPID",
    particles = "trackRefsForJetsBeforeSorting4DwithPID",
    trackTimeTag = "tofPID4DnoPID:t0safe",
    trackTimeResoTag = "tofPID4DnoPID:sigmat0safe"
)
offlinePrimaryVertices4DwithPIDWithBS = offlinePrimaryVertices4DwithPID.clone(
    vertices = "unsortedOfflinePrimaryVertices4DwithPID:WithBS"
)



unsortedOfflinePrimaryVerticesGNN = unsortedOfflinePrimaryVertices4D.clone(
    TkClusParameters = cms.PSet(algorithm = cms.string("GNN2D_vect"),
        TkDAClusParameters = cms.PSet(
            # VertexSlotModel-specific parameters (Pure GNN mode)
            existenceThreshold = cms.double(0.5),           # Slot existence threshold
            trackAssignmentThreshold = cms.double(0.0),     # No threshold - use all assignments
            # Note: minTrackCount removed - rely on ndof for downstream filtering
            numSlots = cms.int32(200),
            nnVersion = cms.string("vertex_slot_v8p9"),
            onnxBackend = cms.string("CUDA"),
            # Inherited parameters
            Tmin = cms.double(4.0),
            Tpurge = cms.double(4.0),
            Tstop = cms.double(2.0),
            vertexSize = cms.double(0.01),
            d0CutOff = cms.double(3.0),
            onnxModelPath = cms.FileInPath('RecoVertex/PrimaryVertexProducer/data/vertex_slot_model.onnx'),
        )
    ),
    TrackTimesLabel = "tofPID4DnoPID:t0safe",
    TrackTimeResosLabel = "tofPID4DnoPID:sigmat0safe",
    # Pure GNN: Skip fitter, use GNN vertices directly
    vertexCollections = cms.VPSet(
        cms.PSet(
            label = cms.string(""),
            algorithm = cms.string(""),  # Empty = no fitter, use GNN clusters directly
            minNdof = cms.double(-1.0),  # Accept all vertices (GNN provides its own filtering)
            useBeamConstraint = cms.bool(False),
            maxDistanceToBeam = cms.double(1.0),
            # Skip time estimator - GNN provides time via t_hat
            vertexTimeParameters = cms.PSet(algorithm = cms.string(""))
        ),
        cms.PSet(
            label = cms.string("WithBS"),
            algorithm = cms.string(""),
            minNdof = cms.double(-1.0),
            useBeamConstraint = cms.bool(True),
            maxDistanceToBeam = cms.double(1.0),
            vertexTimeParameters = cms.PSet(algorithm = cms.string(""))
        )
    )
)
trackWithVertexRefSelectorBeforeSortingGNN = trackWithVertexRefSelector.clone(
    vertexTag = "unsortedOfflinePrimaryVerticesGNN",
    ptMax = 9e99,
    ptErrorCut = 9e99
)
trackRefsForJetsBeforeSortingGNN = trackRefsForJets.clone(
    src = "trackWithVertexRefSelectorBeforeSortingGNN"
)
offlinePrimaryVerticesGNN = sortedPrimaryVertices.clone(
    vertices = "unsortedOfflinePrimaryVerticesGNN",
    particles = "trackRefsForJetsBeforeSortingGNN",
    trackTimeTag = "tofPID4DnoPID:t0safe",
    trackTimeResoTag = "tofPID4DnoPID:sigmat0safe",
    assignment = dict(useTiming = True)
)


from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import tpClusterProducer
from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import quickTrackAssociatorByHits
from SimTracker.TrackAssociation.trackTimeValueMapProducer_cfi import trackTimeValueMapProducer
from RecoMTD.TimingIDTools.tofPIDProducer_cfi import tofPIDProducer
from RecoMTD.TimingIDTools.tofPIDProducerMCTruth_cfi import tofPIDProducerMCTruth

#tofPIDMC=tofPIDProducerMCTruth.clone()
tofPID4DnoPID=tofPIDProducer.clone(vtxsSrc='unsortedOfflinePrimaryVertices')
tofPID=tofPIDProducer.clone()
tofPID3D=tofPIDProducer.clone(vtxsSrc='unsortedOfflinePrimaryVertices')
tofPIDGNN=tofPIDProducer.clone(vtxsSrc='unsortedOfflinePrimaryVerticesGNN')

from Configuration.Eras.Modifier_phase2_timing_layer_cff import phase2_timing_layer
phase2_timing_layer.toModify(tofPID, vtxsSrc='unsortedOfflinePrimaryVertices4D', vertexReassignment=False)
phase2_timing_layer.toModify(tofPID3D, vertexReassignment=False)
phase2_timing_layer.toModify(tofPIDGNN, vertexReassignment=False)
