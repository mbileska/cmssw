import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase2L1GCTSumEmulator_cfi import (
    l1tPhase2L1GCTSumEmulator,
)
from L1Trigger.L1CaloTrigger.l1tPhase2L1MHHFromHFAndDMX_cff import (
    l1tPhase2L1DMXEmulatorNegative,
    l1tPhase2L1DMXEmulatorPositive,
    l1tPhase2L1HFMHHLinkProducerNegative,
    l1tPhase2L1HFMHHLinkProducerPositive,
    l1tPhase2L1MHHFromHFAndDMXNegative,
    l1tPhase2L1MHHFromHFAndDMXPositive,
)


def _gctsum_side_inputs(mhh_label, zero_label):
    # Input triplet 0 is MHH.  The three GCTB triplets are deliberately zero
    # in this smoke chain until GCTB provides the per-side sums-v2 contract.
    inputs = [
        cms.InputTag(mhh_label, "LinkOut{}".format(link))
        for link in range(3)
    ]
    inputs.extend(
        cms.InputTag(zero_label, "LinkIn3")
        for _ in range(9)
    )
    return inputs


l1tPhase2L1GCTSumFromMHHSmoke = l1tPhase2L1GCTSumEmulator.clone(
    inputLinks=cms.VInputTag(
        *(
            _gctsum_side_inputs(
                "l1tPhase2L1MHHFromHFAndDMXPositive",
                "l1tPhase2L1DMXEmulatorPositive",
            )
            + _gctsum_side_inputs(
                "l1tPhase2L1MHHFromHFAndDMXNegative",
                "l1tPhase2L1DMXEmulatorNegative",
            )
        )
    )
)

l1tPhase2L1MHHGCTSumSmokeTask = cms.Task(
    l1tPhase2L1HFMHHLinkProducerPositive,
    l1tPhase2L1HFMHHLinkProducerNegative,
    l1tPhase2L1DMXEmulatorPositive,
    l1tPhase2L1DMXEmulatorNegative,
    l1tPhase2L1MHHFromHFAndDMXPositive,
    l1tPhase2L1MHHFromHFAndDMXNegative,
    l1tPhase2L1GCTSumFromMHHSmoke,
)
