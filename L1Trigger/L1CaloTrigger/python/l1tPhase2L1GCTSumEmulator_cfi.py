import FWCore.ParameterSet.Config as cms

phase2L1GCTSumEmulator = cms.EDProducer(
    "Phase2L1GCTSumEmulator",
    debug = cms.bool(False),
    inputLinks = cms.VInputTag(
        cms.InputTag("gctSumTestVectorProducer", "LinkIn0"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn1"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn2"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn3"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn4"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn5"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn6"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn7"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn8"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn9"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn10"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn11"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn12"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn13"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn14"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn15"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn16"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn17"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn18"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn19"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn20"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn21"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn22"),
        cms.InputTag("gctSumTestVectorProducer", "LinkIn23"),
    )
)