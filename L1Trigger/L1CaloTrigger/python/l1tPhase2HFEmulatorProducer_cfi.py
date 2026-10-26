import FWCore.ParameterSet.Config as cms

l1tPhase2HFEmulatorProducer = cms.EDProducer("Phase2L1CaloL1HFEmulator",
    hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis"),
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(l1tPhase2HFEmulatorProducer,
    hcalDigis = cms.InputTag("mixData", "HcalTriggerPrimitives"),
)