import FWCore.ParameterSet.Config as cms

l1tPhase2L1DMXEmulator = cms.EDProducer(
    "Phase2L1DMXEmulator",
    hgcalTowers=cms.InputTag("l1tHGCalTowerProducer", "HGCalTowerProcessor"),
    hgcalMulticlusters=cms.InputTag(
        "l1tHGCalBackEndLayer2Producer",
        "HGCalBackendLayer2Processor3DClustering",
    ),
    zSide=cms.int32(1),
    etaMin=cms.double(1.479),
    etaMax=cms.double(3.045),
    etLsb=cms.double(1.0),
    phiOffset=cms.double(0.0),
    minTowerEt=cms.double(1.0),
    minEgEt=cms.double(0.0),
    maxEgHoverE=cms.double(0.3),
    jetSeedThreshold=cms.double(3.0),
    duplicateEGAsIso=cms.bool(False),
    debug=cms.bool(False),
)
