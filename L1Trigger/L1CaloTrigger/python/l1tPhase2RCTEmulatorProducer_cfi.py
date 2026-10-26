import FWCore.ParameterSet.Config as cms

l1tPhase2RCTEmulatorProducer = cms.EDProducer("Phase2L1CaloL1RCTEmulator",
    ecalTPEB = cms.InputTag("simEcalEBTriggerPrimitiveDigis"),
    hcalTP = cms.InputTag("simHcalTriggerPrimitiveDigis"),
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(l1tPhase2RCTEmulatorProducer,
    ecalTPEB = cms.InputTag("DMEcalEBTriggerPrimitiveDigis"),
    hcalTP = "DMHcalTriggerPrimitiveDigis",
)