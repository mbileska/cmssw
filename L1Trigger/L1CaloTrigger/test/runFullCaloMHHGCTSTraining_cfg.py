"""Training tree for the full-calorimeter DISTPOC EDM output.

Requires mbileska/L1CaloPhase2Analyzer at the pinned revision documented in
DIST_TRIGGER_POC.md to be checked out under CMSSW src/L1Trigger.
"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.maxEvents = -1
options.register("dataset", "UNKNOWN", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Dataset provenance")
options.register("sourceLFN", "UNKNOWN", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Original input LFN")
options.register("ntupleFile", "fullcalo_training.root", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "TFileService output")
options.register("cmsswCommit", "UNKNOWN", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Local CMSSW commit")
options.register("analyzerCommit", "UNKNOWN", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Analyzer commit")
options.parseArguments()

process = cms.Process("DISTTRAIN")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("L1Trigger.L1CaloPhase2Analyzer.mhhGCTSTrainingNtuplizer_cfi")
process.load("L1Trigger.L1CaloTrigger.l1tPhase2L1FullCaloStageNtuplizer_cfi")
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
process.options = cms.untracked.PSet(numberOfThreads=cms.untracked.uint32(1),
                                     numberOfStreams=cms.untracked.uint32(1),
                                     wantSummary=cms.untracked.bool(True))

analyzer = process.mhhGCTSTrainingNtuplizer
# Required source order: MHH followed by B1, B2, B3 for positive eta,
# then the same four sources for negative eta.
gcts_inputs = []
for side, mhh in (("Positive", "l1tPhase2L1MHHFromHFAndDMXPositive"),
                  ("Negative", "l1tPhase2L1MHHFromHFAndDMXNegative")):
    gcts_inputs.extend(cms.InputTag(mhh, "LinkOut" + str(i)) for i in range(3))
    for card in range(1, 4):
        gcts_inputs.extend(
            cms.InputTag("l1tPhase2L1GCTB{}ToGCTS".format(card),
                         "LinkOut{}{}".format(side, i))
            for i in range(3)
        )
analyzer.gctsInputLinks = cms.VInputTag(*gcts_inputs)
analyzer.gctsSideLinks = cms.VInputTag(
    *[cms.InputTag("l1tPhase2L1GCTSumFullCalo", "SumLinkOut" + str(i))
      for i in range(6)]
)
analyzer.gtLinks = cms.VInputTag(
    *[cms.InputTag("l1tPhase2L1GCTSumFullCalo", "LinkOut" + str(i))
      for i in range(6)]
)
analyzer.gctbPresent = cms.bool(True)
analyzer.gctbSourceMask = cms.uint32(7)
analyzer.dataset = cms.string(options.dataset)
analyzer.sourceLFN = cms.string(options.sourceLFN)
analyzer.inputEDM = cms.string(options.inputFiles[0] if options.inputFiles else "UNKNOWN")
analyzer.cmsswRelease = cms.string("CMSSW_20_0_0_patch2")
analyzer.globalTag = cms.string("200X_mcRun4_realistic_v1")
analyzer.cmsswCommit = cms.string(options.cmsswCommit)
analyzer.analyzerCommit = cms.string(options.analyzerCommit)

process.TFileService = cms.Service("TFileService", fileName=cms.string(options.ntupleFile))
process.training = cms.Path(analyzer + process.l1tPhase2L1FullCaloStageNtuplizer)
