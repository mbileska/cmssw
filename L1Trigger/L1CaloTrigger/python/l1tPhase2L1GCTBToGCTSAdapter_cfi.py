import FWCore.ParameterSet.Config as cms
import math

l1tPhase2L1GCTBToGCTSAdapter = cms.EDProducer(
    "Phase2L1GCTBToGCTSAdapter",
    preIP2=cms.InputTag("l1tPhase2GCTEmulatorProducer", "GCT1PreIP2"),
    postIP2=cms.InputTag("l1tPhase2GCTEmulatorProducer", "GCT1PostIP2"),
    cardIndex=cms.uint32(0),
    # Phase of the 24 ST centers relative to absolute phi=0. This is a
    # proof-of-concept convention pending a validated RCT/GCTB geometry map.
    # RCT card pair 0 begins at -pi in RCT_IO.h. Verify this mapping against
    # the intended GCTB card geometry before interpreting global MET/phi.
    phiPhaseRadians=cms.double(-math.pi),
)
