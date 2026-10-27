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


# One DMX/MHH instance is required for each eta side.  The provisional DMX
# producer fills HGCAL links 8--10; HF links 0--2 remain zero until an HF
# source is connected.
l1tPhase2L1DMXEmulatorPositive = l1tPhase2L1DMXEmulator.clone(zSide=1)
l1tPhase2L1DMXEmulatorNegative = l1tPhase2L1DMXEmulator.clone(zSide=-1)

# Keep the original positive-side label for backward compatibility.
l1tPhase2L1MHHFromDMX = l1tPhase2L1MHHEmulator.clone(
    inputLinks=_dmx_input_links("l1tPhase2L1DMXEmulatorPositive")
)
l1tPhase2L1MHHFromDMXNegative = l1tPhase2L1MHHEmulator.clone(
    inputLinks=_dmx_input_links("l1tPhase2L1DMXEmulatorNegative")
)

l1tPhase2L1MHHFromDMXTask = cms.Task(
    l1tPhase2L1DMXEmulatorPositive,
    l1tPhase2L1DMXEmulatorNegative,
    l1tPhase2L1MHHFromDMX,
    l1tPhase2L1MHHFromDMXNegative,
)
