# Distributed calorimeter trigger proof of concept

## Status and source basis

This is the `dist-trigger-poc` branch of `https://github.com/mbileska/cmssw.git`.
It starts from `emulator-gcts-mhh-clean` at
`22778e9794fc65ccf0866af773275a2529b75159` and targets
`CMSSW_20_0_0_patch2`, `el9_amd64_gcc13`. The independent emulator source
reviewed on 2026-11-06 is
`cms-caloL1-emulator/cmssw:gct-update-latest-firmware` at
`62ab98fc3dd85dc31a04aa0e48a03818e619b289`. The companion user-owned
analyzer is `mbileska/L1CaloPhase2Analyzer:dist-trigger-poc` at
`d3879c57c5c34d7ef47d15252f6317c2f9297f44`.

The branch contains the current CaloL1 RCT, GCTB IP1/IP2, HF and HGCAL
emulator code plus the user's ABI-v2 MHH/GCTS chain. The latest upstream IP1
and IP2 code and their deterministic tests were imported from the
revision above (only trailing whitespace was removed from
`bitonicSort32IP2_cpp.h`); `Phase2L1CaloL1GCTEmulator.cc` calls the namespaced IP1 entry
point. The existing MHH-compatible DMX implementation is retained because the
upstream DMX still uses the earlier sum layout and a different ET LSB. The
additive HF source producer is retained because the old HF mixer output loses
the source region needed to reconstruct global phi.

This is a research software chain, **not** a claim of firmware-exact,
physics-validated full-calorimeter output. The local development environment
has no CVMFS/CMSSW runtime; this branch has not yet passed a `scram b` or a
full-chain event run. The gates below must pass on an AVX2 Wisconsin Condor
worker before using any output for LNN training.

## Event graph and stored data

```text
ECAL EB TPs + HCAL barrel TPs -> 24 RCT cards -> GCTB IP1 -> GCTB IP2
                                                        | three cards, each eta side
                                                        v
HF TPs -> HF link producer --\                     GCTB-to-GCTS ABI-v2 adapter
                              MHH positive/negative ------------\
HGCAL towers/clusters -> DMX-/                            GCTS -> six GT links
```

`runFullCaloMHHGCTS_cfg.py` writes one replay EDM ROOT file per input shard.
It keeps RCT `LinkOut0..3`, GCTB `PreIP1/PostIP1/PreIP2/PostIP2`, three
adapter triplets per side, HF/DMX/MHH links, GCTS side sums and GT links, and
generator jets/MET. Each frame is 9 x 64-bit words; object frames contain
densely packed 48-bit records. The output module retains the raw products,
which can be replayed without rerunning the emulator.

`runFullCaloMHHGCTSTraining_cfg.py` reads that EDM file and writes a second
ROOT file with two synchronized TTrees:

* `mhhGCTSTrainingNtuplizer/trainingTree`: the user-owned analyzer's HF, DMX,
  MHH, explicit 24 GCTS source links, side outputs, GT outputs, decoded
  quantities in hardware/physical units, truth, and provenance. Its
  `gctbPresent=true`, `gctbSourceMask=7`, `gctsInputDerived=false` distinguish
  this from the older zero-barrel campaign.
* `l1tPhase2L1FullCaloStageNtuplizer/stages`: all RCT and GCTB stage frames,
  decoded RCT tower ET, core GCTB supertower ET, and GCTB IP2 EG/jet/tau pT.
  All decoded pT/ET fields here use a 0.5 GeV LSB. Raw frames remain available
  for redecoding if interface assumptions change.

Flattened stage branch order is:

| Branch | Order | Length per event |
|---|---|---:|
| `rctWords` | link 0..3, card 0..23, word 0..8 | 864 |
| `gctbPreIP1Words` | card 1..3, SLR3 then SLR1, frame 0..31, word 0..8 | 1728 |
| `gctbPostIP1Words` | same SLR order, frame 0..20, word 0..8 | 1134 |
| `gctbPreIP2Words` | card 1..3, frame 0..23, word 0..8 | 648 |
| `gctbPostIP2Words` | card 1..3, frame 0..5, word 0..8 | 162 |
| `rctTowerEtGeV` | card, RCT tower link 1..3, 34 records/link | 2448 |
| `gctbSTEtGeV` | card, positive/negative, eta 0..5, local phi 0..7 | 288 |
| `gctb(EG/Jet/Tau)PtGeV` | card, positive/negative, object 0..5 | 36 each |

## Interface decisions

1. The upstream IP1 update changes PF/EG layouts and sorting, and IP2 now
   runs all six jet and tau iterations. Its source regression tests are built
   by `scram b runtests` along with the pre-existing MHH/GCTS firmware tests.
2. RCT/GCTB routing follows the already checked-in explicit 24-card,
   six-SLR map. GCTB PreIP2 transposes the 12 ST records in each frame because
   IP1 and IP2 use different record orders. Both sides of this permutation
   are persisted.
3. The current GCTB IP2 emits six EG, six jets, six taus per eta side; it
   leaves iso-EG empty, has padded jet/tau coordinates, and supplies only a
   combined legacy sum. The new adapter leaves six iso-EG slots **zero**,
   strips unspecified high bits, subtracts the two-cell IP2 padding from
   jet/tau eta and phi, and emits three source frames per eta side. IP2
   selects fixed EG indices 0,1,8,9,15,22; the adapter restores their
   RCT-card phi offset and reverses negative-eta local phi using the current
   RCT card ordering. It does not silently label an unavailable barrel
   quantity as measured.
4. Barrel Ex/Ey/total ET/NObj are computed from the eight unique ST phi
   columns in GCTB PreIP2 per card/side. HTx/HTy/total HT are computed from
   the selected six IP2 jets. Source sums use signed 12-bit components and
   unsigned 12-bit scalars at 0.5 GeV/LSB; values clamp to those widths.
   NObj currently counts nonzero STs. This is a **POC software policy**, not
   upstream GCTB firmware. The adapter publishes `AssumptionFlags=15` for
   every card: bits 0..3 mean missing iso-EG, derived sums, normalized padded
   coordinates, and assumed phi alignment.
5. RCT card pair zero is at `-pi` in `RCT_IO.h`. The adapter therefore uses
   `phiPhaseRadians=-pi` and 24 ST bins of 15 degrees. The checked-in RCT
   segmentation differs from the interface note's Figure 6 card grid; card
   ownership and especially phi at boundaries require detector-level checks.
   Raw local coordinates are preserved so this can be revisited without
   rerunning the upstream stages. In particular, GCTS currently adds barrel
   card offsets 0/120/240 crystals and 0/8/16 STs, while the checked-in RCT
   card-zero boundary starts at -pi. The resulting non-sum global phi may be
   rotated relative to HGCAL/HF; detector-wide EG/jet/tau phi and barrel to
   endcap stitching are **not yet validated**. The per-source Cartesian
   sums use the RCT physical `-pi` phase so their Ex/Ey are in the physical
   frame; their agreement must be checked separately from packed object phi.
6. GCTS receives MHH + GCTB1 + GCTB2 + GCTB3, independently on positive
   and negative eta. GCTS applies its existing barrel offsets and side
   reductions, then TO GT combines side vectors before MET/MHT magnitude and
   angle. Small MHH/GCTS/GT MET differences can arise from fixed-point
   quantization, final vector magnitude, and clipping. A discrepancy larger
   than the encoded LSB/tolerance requires investigation.

## Required validation gates

* Confirm the chosen `GEN-SIM-DIGI-RAW` file actually contains
  `simEcalEBTriggerPrimitiveDigis`, `simHcalTriggerPrimitiveDigis`,
  `l1tHGCalTowerProducer:HGCalTowerProcessor`, and
  `l1tHGCalBackEndLayer2Producer:HGCalBackendLayer2Processor3DClustering`.
  The previous DY campaign checked only HCAL and HGCAL. The RCT producer now
  gives a clear missing-input exception if ECAL/HCAL products are absent.
* Build on a compatible AVX2 worker with the `CMSSW_20_0_0_patch2` runtime.
  Wisconsin login nodes previously failed in `libActsCore`, and the DAVS
  path previously crashed in Davix. Use Condor and the AAA `root://` PFN.
* First run 10 events, then one complete file, then the campaign. Require
  `cmsRun` status 0, nonempty replay EDM and training ROOT, matching event
  IDs/counts, nonzero barrel activity, source-to-GCTS side sum conservation,
  and plausible detector-coordinate plots before scaling.
* Check phi and eta against generator/reference objects and inspect boundary
  events. The GCTB IP2 EG selection is still fixed-index rather than a final
  pT-ordered/isolation implementation. Current GCTB EG fields lack fully
  validated global coordinates; derived ST sums do not establish final MET
  physics performance.

## Wisconsin run recipe

Set up the CMSSW area and check out this branch under
`$CMSSW_BASE/src/L1Trigger/L1CaloTrigger`. Check out the pinned analyzer under
`$CMSSW_BASE/src/L1Trigger/L1CaloPhase2Analyzer`. On a build-capable host:

```bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$CMSSW_BASE/src"
eval "$(scram runtime -sh)"
scram b -j 8
scram b runtests
```

On an AVX2 Condor worker with a valid transferred CMS proxy, choose an LFN
from DAS and set `PFN="root://cmsxrootd.fnal.gov/${LFN}"`. Preflight with
`edmDumpEventContent "$PFN"` and check the four products above. Then:

```bash
cd "$CMSSW_BASE/src"
cmsRun L1Trigger/L1CaloTrigger/test/runFullCaloMHHGCTS_cfg.py \
  inputFiles="$PFN" maxEvents=10 outputFile=file:fullcalo_edm.root
edmFileUtil -e fullcalo_edm.root
cmsRun L1Trigger/L1CaloTrigger/test/runFullCaloMHHGCTSTraining_cfg.py \
  inputFiles=file:"$PWD/fullcalo_edm.root" maxEvents=-1 \
  ntupleFile=fullcalo_training.root dataset="$DATASET" sourceLFN="$LFN" \
  cmsswCommit="$(git rev-parse HEAD)" \
  analyzerCommit="$(git -C L1Trigger/L1CaloPhase2Analyzer rev-parse HEAD)"
python3 L1Trigger/L1CaloPhase2Analyzer/test/check_mhh_gcts_training_ntuple.py \
  fullcalo_training.root --expected-edm fullcalo_edm.root \
  --expect-full-gctb --require-truth
python3 L1Trigger/L1CaloTrigger/test/checkFullCaloTraining.py \
  fullcalo_training.root --min-events 10
```

`VarParsing` may decorate the EDM filename for finite `maxEvents`; locate
`fullcalo_edm*.root` and pass the actual name to the following commands.
The existing analyzer campaign script validates the older zero-GCTB smoke
mode and must **not** be used unchanged for this full-calorimeter branch.
After the 10-event gate succeeds, `submitFullCaloCampaign.sh` can submit one
Condor shard per DAS file. Start with one file, inspect its two outputs and
validation text, then use `--limit-files 0 --max-events -1` for the full
dataset:

```bash
cd "$CMSSW_BASE/src/L1Trigger/L1CaloTrigger"
./test/submitFullCaloCampaign.sh \
  --dataset '/DYToLL_M-10To50_TuneCP5_14TeV-pythia8/Phase2Fall26DIGI-DIGI_pilot_200X_mcRun4_realistic_v1-v3/GEN-SIM-DIGI-RAW' \
  --limit-files 1 --max-events -1
# After the first full shard passes, repeat with --limit-files 0.
```

The submitter checks the exact reviewed analyzer revision, the current
CMSSW branch name, a valid proxy, and AVX2 Wisconsin worker matching. Its
worker checks input products, the replay EDM products, the existing analyzer
schema/truth checker, and `checkFullCaloTraining.py` before writing a SUCCESS
marker and checksums. A nonzero worker exit holds the Condor job. The
campaign directory is printed by the submitter; monitor from the same login
host used to submit because Wisconsin has multiple schedulers.

Store the final validated ROOT files, logs, exact configs, commit hashes,
LFNs, and SHA256 checksums in backed-up storage; `/scratch` is temporary.
