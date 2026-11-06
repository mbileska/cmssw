import FWCore.ParameterSet.Config as cms


def _tags(prefix, instances):
    return cms.VInputTag(*(cms.InputTag(prefix, name) for name in instances))


_gct = "l1tPhase2GCTEmulatorProducer"
l1tPhase2L1FullCaloStageNtuplizer = cms.EDAnalyzer(
    "Phase2L1FullCaloStageNtuplizer",
    rctLinks=_tags("l1tPhase2RCTEmulatorProducer", ["LinkOut{}".format(i) for i in range(4)]),
    preIP1Links=_tags(_gct, ["GCT{}SLR{}PreIP1".format(card, slr)
                            for card in range(1, 4) for slr in (3, 1)]),
    postIP1Links=_tags(_gct, ["GCT{}SLR{}PostIP1".format(card, slr)
                             for card in range(1, 4) for slr in (3, 1)]),
    preIP2Links=_tags(_gct, ["GCT{}PreIP2".format(card) for card in range(1, 4)]),
    postIP2Links=_tags(_gct, ["GCT{}PostIP2".format(card) for card in range(1, 4)]),
    adapterFlags=cms.VInputTag(
        *(cms.InputTag("l1tPhase2L1GCTB{}ToGCTS".format(card), "AssumptionFlags")
          for card in range(1, 4))
    ),
)
