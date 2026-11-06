import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase2RCTEmulatorProducer_cfi import (
    l1tPhase2RCTEmulatorProducer,
)
from L1Trigger.L1CaloTrigger.l1tPhase2GCTEmulatorProducer_cfi import (
    l1tPhase2GCTEmulatorProducer,
)
from L1Trigger.L1CaloTrigger.l1tPhase2L1GCTBToGCTSAdapter_cfi import (
    l1tPhase2L1GCTBToGCTSAdapter,
)
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


def _barrel_adapter(card):
    return l1tPhase2L1GCTBToGCTSAdapter.clone(
        cardIndex=card,
        preIP2=cms.InputTag("l1tPhase2GCTEmulatorProducer", "GCT{}PreIP2".format(card + 1)),
        postIP2=cms.InputTag("l1tPhase2GCTEmulatorProducer", "GCT{}PostIP2".format(card + 1)),
    )


l1tPhase2L1GCTB1ToGCTS = _barrel_adapter(0)
l1tPhase2L1GCTB2ToGCTS = _barrel_adapter(1)
l1tPhase2L1GCTB3ToGCTS = _barrel_adapter(2)


def _side_inputs(side, mhh):
    result = [cms.InputTag(mhh, "LinkOut{}".format(link)) for link in range(3)]
    for card in range(1, 4):
        label = "l1tPhase2L1GCTB{}ToGCTS".format(card)
        result.extend(
            cms.InputTag(label, "LinkOut{}{}".format(side, link))
            for link in range(3)
        )
    return result


l1tPhase2L1GCTSumFullCalo = l1tPhase2L1GCTSumEmulator.clone(
    inputLinks=cms.VInputTag(
        *(
            _side_inputs("Positive", "l1tPhase2L1MHHFromHFAndDMXPositive")
            + _side_inputs("Negative", "l1tPhase2L1MHHFromHFAndDMXNegative")
        )
    )
)

l1tPhase2L1FullCaloMHHGCTSTask = cms.Task(
    l1tPhase2RCTEmulatorProducer,
    l1tPhase2GCTEmulatorProducer,
    l1tPhase2L1GCTB1ToGCTS,
    l1tPhase2L1GCTB2ToGCTS,
    l1tPhase2L1GCTB3ToGCTS,
    l1tPhase2L1HFMHHLinkProducerPositive,
    l1tPhase2L1HFMHHLinkProducerNegative,
    l1tPhase2L1DMXEmulatorPositive,
    l1tPhase2L1DMXEmulatorNegative,
    l1tPhase2L1MHHFromHFAndDMXPositive,
    l1tPhase2L1MHHFromHFAndDMXNegative,
    l1tPhase2L1GCTSumFullCalo,
)
