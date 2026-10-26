"""
Description
-----------
CMSSW configuration fragment for the complete GCT IP1 -> IP2 emulator.
The producer consumes the four RCT output products and publishes PreIP1,
PostIP1, PreIP2, and PostIP2 link collections.
"""

import FWCore.ParameterSet.Config as cms

l1tPhase2GCTEmulatorProducer = cms.EDProducer(
    "Phase2L1CaloL1GCTEmulator",
    LinkOut0=cms.InputTag("l1tPhase2RCTEmulatorProducer", "LinkOut0"),
    LinkOut1=cms.InputTag("l1tPhase2RCTEmulatorProducer", "LinkOut1"),
    LinkOut2=cms.InputTag("l1tPhase2RCTEmulatorProducer", "LinkOut2"),
    LinkOut3=cms.InputTag("l1tPhase2RCTEmulatorProducer", "LinkOut3"),
)
