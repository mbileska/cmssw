import FWCore.ParameterSet.Config as cms

phase2L1MHHEmulator = cms.EDProducer(
    "Phase2L1MHHEmulator",
    inputLinks=cms.VInputTag(),
    debug=cms.bool(False),
)
