import FWCore.ParameterSet.Config as cms

l1tPhase2MLonRCTSetupProducer = cms.EDProducer("Phase2L1CaloMLonRCTSetupProducer",
    ecalTPEB = cms.InputTag("simEcalEBTriggerPrimitiveDigis"),
    hcalTP = cms.InputTag("simHcalTriggerPrimitiveDigis"),
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(l1tPhase2MLonRCTSetupProducer,
    ecalTPEB = cms.InputTag("DMEcalEBTriggerPrimitiveDigis"),
    hcalTP = "DMHcalTriggerPrimitiveDigis",
)