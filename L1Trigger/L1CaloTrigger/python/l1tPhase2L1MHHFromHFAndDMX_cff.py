import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase2L1DMXEmulator_cfi import (
    l1tPhase2L1DMXEmulator,
)
from L1Trigger.L1CaloTrigger.l1tPhase2L1HFMHHLinkProducer_cfi import (
    l1tPhase2L1HFMHHLinkProducer,
)
from L1Trigger.L1CaloTrigger.l1tPhase2L1MHHEmulator_cfi import (
    l1tPhase2L1MHHEmulator,
)


def _mhh_input_links(hf_label, dmx_label):
    # HF supplies links 0--2.  DMX supplies HGCAL links 8--10 and initialized
    # zero frames for all other positions in the 16-link MHH input contract.
    return cms.VInputTag(
        *[
            cms.InputTag(
                hf_label if link < 3 else dmx_label,
                "LinkIn{}".format(link),
            )
            for link in range(16)
        ]
    )


l1tPhase2L1HFMHHLinkProducerPositive = l1tPhase2L1HFMHHLinkProducer.clone(
    zSide=1
)
l1tPhase2L1HFMHHLinkProducerNegative = l1tPhase2L1HFMHHLinkProducer.clone(
    zSide=-1
)
l1tPhase2L1DMXEmulatorPositive = l1tPhase2L1DMXEmulator.clone(zSide=1)
l1tPhase2L1DMXEmulatorNegative = l1tPhase2L1DMXEmulator.clone(zSide=-1)

l1tPhase2L1MHHFromHFAndDMXPositive = l1tPhase2L1MHHEmulator.clone(
    inputLinks=_mhh_input_links(
        "l1tPhase2L1HFMHHLinkProducerPositive",
        "l1tPhase2L1DMXEmulatorPositive",
    )
)
l1tPhase2L1MHHFromHFAndDMXNegative = l1tPhase2L1MHHEmulator.clone(
    inputLinks=_mhh_input_links(
        "l1tPhase2L1HFMHHLinkProducerNegative",
        "l1tPhase2L1DMXEmulatorNegative",
    )
)

l1tPhase2L1MHHFromHFAndDMXTask = cms.Task(
    l1tPhase2L1HFMHHLinkProducerPositive,
    l1tPhase2L1HFMHHLinkProducerNegative,
    l1tPhase2L1DMXEmulatorPositive,
    l1tPhase2L1DMXEmulatorNegative,
    l1tPhase2L1MHHFromHFAndDMXPositive,
    l1tPhase2L1MHHFromHFAndDMXNegative,
)
