import FWCore.ParameterSet.Config as cms

gctSumFileInputProducer = cms.EDProducer(
    "GCTSumFileInputProducer",
    inputFile=cms.string(""),
    posOffset=cms.uint32(0),
    negOffset=cms.uint32(12),
    eventOffset=cms.uint32(0),
    wrapAround=cms.bool(False),
    debug=cms.bool(False),
)
