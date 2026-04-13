import FWCore.ParameterSet.Config as cms

gctSumTestVectorProducer = cms.EDProducer(
    "GCTSumTestVectorProducer",
    patternMode = cms.string("validation"),
    debug = cms.bool(False),
)
