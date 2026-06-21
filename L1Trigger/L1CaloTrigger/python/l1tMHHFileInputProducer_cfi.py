import FWCore.ParameterSet.Config as cms

mhhFileInputProducer = cms.EDProducer(
    "MHHFileInputProducer",
    inputFile=cms.string(""),
    payloadOffset=cms.uint32(0),
    eventOffset=cms.uint32(0),
    wrapAround=cms.bool(False),
    debug=cms.bool(False),
)
