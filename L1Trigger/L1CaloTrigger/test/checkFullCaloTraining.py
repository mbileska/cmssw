#!/usr/bin/env python3
"""Check that the training tree contains real barrel inputs and conserved sums."""

import argparse

import ROOT


def signed(value, bits):
    value &= (1 << bits) - 1
    return value - (1 << bits) if value & (1 << (bits - 1)) else value


def source_sums(words):
    first, second = (int(words[0]), int(words[1]))
    return (
        signed(first, 12), signed(first >> 12, 12),
        signed(first >> 24, 12), signed(first >> 36, 12),
        second & 0xFFF, (second >> 12) & 0xFFF, (second >> 24) & 0xFFF,
    )


def side_sums(words):
    first, second = (int(words[0]), int(words[1]))
    return (
        signed(first, 16), signed(first >> 16, 16),
        signed(first >> 32, 16), signed(first >> 48, 16),
        second & 0xFFFF, (second >> 16) & 0xFFFF,
        (second >> 32) & 0xFFFF,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ntuple")
    parser.add_argument("--min-events", type=int, default=1)
    args = parser.parse_args()

    root_file = ROOT.TFile.Open(args.ntuple)
    if not root_file or root_file.IsZombie():
        raise RuntimeError("cannot open training ROOT file")
    tree = root_file.Get("mhhGCTSTrainingNtuplizer/trainingTree")
    if not tree:
        raise RuntimeError("trainingTree is absent")
    stages = root_file.Get("l1tPhase2L1FullCaloStageNtuplizer/stages")
    if not stages or stages.GetEntries() != tree.GetEntries():
        raise RuntimeError("stage tree is absent or event count differs")
    if tree.GetEntries() < args.min_events:
        raise RuntimeError("too few entries: {}".format(tree.GetEntries()))

    active_barrel = False
    for entry in tree:
        if not bool(entry.gctbPresent) or int(entry.gctbSourceMask) != 7:
            raise RuntimeError("event {}: GCTB provenance is not full barrel".format(entry.event))
        if bool(entry.gctsInputDerived):
            raise RuntimeError("event {}: GCTS input was synthesized by the analyzer".format(entry.event))
        inputs = tree.GetLeaf("gcts_input_words")
        sides = tree.GetLeaf("gcts_sum_words")
        for side in range(2):
            summed = [0] * 7
            for source in range(4):
                link = source * 3 + 2
                start = (side * 12 + link) * 9
                words = [inputs.GetValueLong64(start + word) for word in range(9)]
                if any(words[2:]):
                    raise RuntimeError("event {}: reserved source sum words are nonzero".format(entry.event))
                fields = source_sums(words)
                summed = [a + b for a, b in zip(summed, fields)]
                if source > 0 and (words[0] or words[1]):
                    active_barrel = True
            side_start = (side * 3 + 2) * 9
            side_words = [sides.GetValueLong64(side_start + word) for word in range(9)]
            observed = side_sums(side_words)
            expected = tuple(signed(value, 16) if index < 4 else value & 0xFFFF
                             for index, value in enumerate(summed))
            if observed != expected:
                raise RuntimeError(
                    "event {} side {}: GCTS sums {} != four source sums {}".format(
                        entry.event, side, observed, expected
                    )
                )
    if not active_barrel:
        raise RuntimeError("all barrel source sums are zero; full barrel path was not exercised")
    for entry in stages:
        expected_lengths = {
            "rctWords": 4 * 24 * 9,
            "gctbPreIP1Words": 6 * 32 * 9,
            "gctbPostIP1Words": 6 * 21 * 9,
            "gctbPreIP2Words": 3 * 24 * 9,
            "gctbPostIP2Words": 3 * 6 * 9,
            "rctTowerEtGeV": 24 * 102,
            "gctbSTEtGeV": 3 * 2 * 48,
            "gctbEGPtGeV": 3 * 2 * 6,
            "gctbJetPtGeV": 3 * 2 * 6,
            "gctbTauPtGeV": 3 * 2 * 6,
            "gctbAdapterFlags": 3,
        }
        for branch, length in expected_lengths.items():
            if len(getattr(entry, branch)) != length:
                raise RuntimeError("event {}: {} has wrong length".format(entry.event, branch))
        if list(entry.gctbAdapterFlags) != [15, 15, 15]:
            raise RuntimeError("event {}: adapter assumptions changed".format(entry.event))
    print("PASS: {} entries, explicit three-card GCTB inputs, active barrel, and source-to-GCTS sum conservation".format(tree.GetEntries()))


if __name__ == "__main__":
    main()
