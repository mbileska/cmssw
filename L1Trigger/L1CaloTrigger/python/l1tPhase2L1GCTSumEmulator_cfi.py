import FWCore.ParameterSet.Config as cms

phase2L1GCTSumEmulator = cms.EDProducer(
    "Phase2L1GCTSumEmulator",
    inputLinks=cms.VInputTag(),
    debug=cms.bool(False),
)
