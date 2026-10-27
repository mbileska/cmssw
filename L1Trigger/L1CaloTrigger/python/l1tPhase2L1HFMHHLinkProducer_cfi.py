import FWCore.ParameterSet.Config as cms

l1tPhase2L1HFMHHLinkProducer = cms.EDProducer(
    "Phase2L1HFMHHLinkProducer",
    hcalDigis=cms.InputTag("simHcalTriggerPrimitiveDigis"),
    zSide=cms.int32(1),
    etLsb=cms.double(0.5),
    minTowerEt=cms.double(0.0),
    jetSeedThreshold=cms.double(3.0),
    debug=cms.bool(False),
)

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2

premix_stage2.toModify(
    l1tPhase2L1HFMHHLinkProducer,
    hcalDigis=cms.InputTag("mixData", "HcalTriggerPrimitives"),
)
