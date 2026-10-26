import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("L1AlgoTest", eras.Phase2C17I13M9)

# Setup
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load("Configuration.Geometry.GeometryExtendedRun4D110_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.SimL1Emulator_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1))

process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(
        "file:///hdfs/store/user/emettner/Misc/00ea2d85-25bf-4c64-8e71-d3d6726028c4.root",
    ),
    inputCommands=cms.untracked.vstring(
        "keep *",
        "drop l1tTkPrimaryVertexs_*_*_*",
    ),
)


#Global Tag
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "131X_mcRun4_realistic_v6", "")

# HCAL TPG prerequisites
process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")
process.load("CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi")

# HF Emulator
process.load("L1Trigger.L1CaloTrigger.l1tPhase2HFEmulatorProducer_cfi")

process.L1simulation_step = cms.Path(process.l1tPhase2HFEmulatorProducer)
process.HF = cms.Path(process.l1tPhase2HFEmulatorProducer)

#Output (HF link products)
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName=cms.untracked.string("hfemu_out.root"),
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep *_l1tPhase2HFEmulatorProducer_LinkOutIP1PosEtaCh*_*",
        "keep *_l1tPhase2HFEmulatorProducer_LinkOutIP1NegEtaCh*_*",
        "keep *_l1tPhase2HFEmulatorProducer_LinkOutIP2PosEtaCh*_*",
        "keep *_l1tPhase2HFEmulatorProducer_LinkOutIP2NegEtaCh*_*",
    ),
)
process.end = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.HF, process.end)

process.options.numberOfThreads = cms.untracked.uint32(8)
process.options.numberOfStreams = cms.untracked.uint32(0)