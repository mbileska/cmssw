# HF/HGCAL -> MHH -> GCT Sum integration

## Scope

The existing `Phase2L1CaloL1HFEmulator` and its firmware headers are not
modified by this integration.  Its published mixer records cannot be adapted
losslessly after the fact: jets and taus from three phi regions are sorted
together after only a region-local phi value has been stored.  The source
region, and therefore global phi, is absent from the output record.

`Phase2L1HFMHHLinkProducer` is an additive MHH source producer.  It consumes
the same HF trigger primitives, reconstructs the global 24-bin hadronic phi,
and emits three nine-word `std::vector<uint64_t>` products per eta side:

| Product | Contents |
|---|---|
| `LinkIn0` | zero; no HF EG/isolated-EG collection is currently defined |
| `LinkIn1` | six 48-bit jets followed by six 48-bit taus |
| `LinkIn2` | source sums ABI v2 |

HF hadronic eta codes start at 6, immediately after HGCAL codes 0--5.  Object
and sum energies use a 0.5 GeV LSB.  Ex/Ey/total ET and NObj use accepted
physical HF towers.  HTx/HTy/total HT use the accepted reconstructed HF jets.
All Cartesian components are global energy-flow sums, not missing-energy
vectors.

This producer makes the CMSSW integration chain executable without changing
HF-owned code.  It is not a claim of bit-exact agreement with HF IP1: the
authoritative HF-to-MHH output must retain global phi (or its source region),
and its object/sum policy still requires approval from the HF owners.

`l1tPhase2L1MHHFromHFAndDMX_cff.py` builds independent positive- and
negative-eta chains.  HF supplies MHH links 0--2; DMX supplies HGCAL links
8--10 and initialized zero frames for the unused links.

## GCT Sum smoke chain

`l1tPhase2L1MHHGCTSumSmoke_cff.py` connects both MHH outputs to GCT Sum.  It
deliberately supplies zero frames for all three GCTB source triplets per side.
This tests only

```
HF + HGCAL -> MHH(+/-) -> GCT Sum -> six logical GT links
```

It is not a full calorimeter-chain physics configuration.  The current GCTB
emulator still needs separate positive/negative sums-v2 output links,
HTx/HTy/ET/HT/NObj generation, EG-isolation output, and a product-type adapter
before its three regional triplets can replace the smoke-chain zero frames.
TMUX6 is intentionally absent from CMSSW: it serializes and replicates the six
event-level GCT Sum output links but does not change their physics payload.
