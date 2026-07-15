import FWCore.ParameterSet.Config as cms

l1tPhase2L1GCTSumEmulator = cms.EDProducer(
    "Phase2L1GCTSumEmulator",
    inputLinks=cms.VInputTag(),
    debug=cms.bool(False),
)
