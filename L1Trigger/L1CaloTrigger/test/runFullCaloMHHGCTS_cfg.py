import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.maxEvents = 10
options.outputFile = "fullcalo_mhh_gcts_edm.root"
options.parseArguments()

process = cms.Process("DISTPOC")
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi")
process.HcalTrigTowerGeometryESProducer = cms.ESProducer("HcalTrigTowerGeometryESProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "200X_mcRun4_realistic_v1", "")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.options = cms.untracked.PSet(
    numberOfThreads=cms.untracked.uint32(1),
    numberOfStreams=cms.untracked.uint32(1),
)
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("L1Trigger.L1CaloTrigger.l1tPhase2L1FullCaloMHHGCTS_cff")
process.l1tPhase2RCTEmulatorProducer.ecalTPEB = cms.InputTag("simEcalEBTriggerPrimitiveDigis")
process.l1tPhase2RCTEmulatorProducer.hcalTP = cms.InputTag("simHcalTriggerPrimitiveDigis")
for module in (
    process.l1tPhase2L1HFMHHLinkProducerPositive,
    process.l1tPhase2L1HFMHHLinkProducerNegative,
):
    module.hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis")
for module in (
    process.l1tPhase2L1DMXEmulatorPositive,
    process.l1tPhase2L1DMXEmulatorNegative,
):
    module.hgcalTowers = cms.InputTag("l1tHGCalTowerProducer", "HGCalTowerProcessor")
    module.hgcalMulticlusters = cms.InputTag(
        "l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClustering"
    )

process.run = cms.Path()
process.run.associate(process.l1tPhase2L1FullCaloMHHGCTSTask)
process.output = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string(options.outputFile),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep *_l1tPhase2RCTEmulatorProducer_*_*",
        "keep *_l1tPhase2GCTEmulatorProducer_*_*",
        "keep *_l1tPhase2L1*_*_*",
        "keep *_ak4GenJetsNoNu_*_*",
        "keep *_genMetTrue_*_*",
    ),
)
process.end = cms.EndPath(process.output)
