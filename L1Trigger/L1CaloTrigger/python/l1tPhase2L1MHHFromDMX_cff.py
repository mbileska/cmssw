import FWCore.ParameterSet.Config as cms

from L1Trigger.L1CaloTrigger.l1tPhase2L1DMXEmulator_cfi import (
    l1tPhase2L1DMXEmulator,
)
from L1Trigger.L1CaloTrigger.l1tPhase2L1MHHEmulator_cfi import (
    l1tPhase2L1MHHEmulator,
)


def _dmx_input_links(module_label):
    return cms.VInputTag(
        *[
            cms.InputTag(module_label, "LinkIn{}".format(link))
            for link in range(16)
        ]
    )


# The current MHH firmware wrapper implements one eta-side/Card-A data path.
# Clone this pair with zSide=-1 when a second MHH instance is configured.
l1tPhase2L1DMXEmulatorPositive = l1tPhase2L1DMXEmulator.clone(zSide=1)

l1tPhase2L1MHHFromDMX = l1tPhase2L1MHHEmulator.clone(
    inputLinks=_dmx_input_links("l1tPhase2L1DMXEmulatorPositive")
)

l1tPhase2L1MHHFromDMXTask = cms.Task(
    l1tPhase2L1DMXEmulatorPositive,
    l1tPhase2L1MHHFromDMX,
)
