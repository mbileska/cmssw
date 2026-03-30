import FWCore.ParameterSet.Config as cms

gctSumTestVectorProducer = cms.EDProducer(
    "GCTSumTestVectorProducer",
    patternMode = cms.string("cyclic"),
    debug = cms.bool(False),
)
