# Scouting NanoAOD with checked AXOL1TL and CICADA scores

This directory holds what you need to produce HLT-scouting NanoAOD samples with two independent
evaluations of both L1 anomaly scores. Starting from MC `GEN-SIM-RAW`, each job:

1. re-runs the L1 trigger emulator (`L1REPACK:FullMC`) with the 2025 L1 menu, AXOL1TL, and
   `CICADAModel_v2p2p0`,
2. re-runs the HLT (`HLT:GRun` menu of this release),
3. independently rebuilds both model tensors and evaluates both models again,
4. writes scouting NanoAOD (`NANO:@ScoutMonitor`).

Every input event is written out, whether or not it fired a scouting trigger. The main columns:

| Branch | Content |
|---|---|
| `L1AXO_score`, `L1AXO_scoreRecomputed`, `L1AXO_scoreDifference` | uGT score, independent score, and independent minus uGT at BX 0 |
| `nL1AXOInput`, `L1AXOInput_value2` | The 57 exact AXOL1TL inputs; divide `value2` by two before inference |
| `L1CICADA_score`, `L1CICADA_scoreRecomputed`, `L1CICADA_scoreDifference` | Layer-1 score, independent v2p2p0 score, and their difference at BX 0 |
| `nL1CICADAInput`, `L1CICADAInput_et`, `L1CICADAInput_gctEta`, `L1CICADAInput_gctPhi`, `L1CICADAInput_present` | The ordered 252-region CICADA tensor and grid-integrity flag |
| `Background_process`, `Background_sample`, `Background_decay`, `Background_kinematicBin`, `Background_category`, `Background_kinematicValue` | Strict exclusive background label and parent-normalization provenance (background configs only) |
| `L1_AXO_VLoose` … `L1_AXO_VVVTight` | Re-emulated L1 AXO trigger bits |
| `DST_PFScouting_AXO*`, other `DST_*`/`HLT_*` | Re-run HLT decisions |
| `*_pHLT` versions of the trigger bits | HLT decisions from the original MC production, kept for comparison. Paths that exist only in the original menu keep their plain name (for example `DST_PFScouting_AXONominal`); the branch description says `(process: HLT)` for original bits and `(process: reL1reHLT)` for re-run ones. |
| `ScoutingParticle_*` (PF candidates), `ScoutingTrack_*`, `ScoutingPFJet_*`, `ScoutingMuonVtx_*`/`ScoutingMuonNoVtx_*`, `ScoutingElectron_*`, `ScoutingPhoton_*`, `ScoutingPrimaryVertex_*`, `ScoutingMET_*`, `ScoutingRho_*`, … | Scouting objects |
| `genWeight` (and `genEventSumw` in the `Runs` tree) | Generator weights for normalisation |

`nano_branches.txt` is the branch inventory from the earlier AXO-only 2025 test job. Regenerate it
after the first successful closure pilot; it does not yet contain the new independent-score,
tensor, or background-label columns. The original-HLT bits also depend on the menu used to make
the input MC.

---

## Contents of this directory

| Path | What it is |
|---|---|
| `README.md` | This file |
| `makeConfig.sh` | Writes the cmsRun config for one year (Section 3) |
| `crab.py` | CRAB submission wrapper (Section 5) |
| `nano_branches.txt` | Pre-change AXO-only branch inventory; regenerate after the closure pilot |
| `samplelists/` | Input lists for CRAB (Section 6) |
| `samplesearch/` | Scripts used to find the input datasets in DAS, and their outputs (Section 6) |
| `background_categories.json` | Authoritative 37-category partition and optional coarse training groups |
| `background_samples_2025.json` | Parent-dataset IDs and normalization ledger; cross sections must be filled before physics normalization |
| `validateTrialNano.py` | Fails on score disagreement, malformed tensors, duplicate events, invalid labels, or weights |
| `stageBackgrounds.py` | Writes one set of Parquet shards per fine category; training groups remain a reversible mapping |

The code changes this workflow needs live in two other packages of the branch:

| Path (under `$CMSSW_BASE/src`) | What it does |
|---|---|
| `L1Trigger/L1TGlobal/src/GlobalBoard.cc` | Fixes AXO score storage: selects the condition when each menu is initialized, preserves that choice across events, and records only its BX-0 score (including in multi-BX emulation). |
| `L1Trigger/L1TGlobal/plugins/L1TGlobalProducer.cc` | Turns on AXO score saving in `beginRun`, before the trigger conditions are built, and gives the score collection a real BX range. |
| `DPGAnalysis/L1TNanoAOD/plugins/L1AnomalyScoreTableProducer.cc` | Rebuilds the 57/252 inputs, reruns both models, checks closure, and writes both versions |
| `DPGAnalysis/L1TNanoAOD/plugins/BackgroundCategoryTableProducer.cc` | Assigns exactly one fine category and records the configured parent sample ID |
| `DPGAnalysis/L1TNanoAOD/python/l1tnanotables_cff.py` | Defines the anomaly-score and input tables |
| `DPGAnalysis/L1TNanoAOD/python/customiseAXOScoutingNano.py` | The customisation `makeConfig.sh` adds to every config (Section 2) |
| `L1Trigger/L1TCaloLayer1/plugins/L1TCaloSummary.cc` | Zero-initializes the official 252-cell CICADA input for deterministic behavior if a region is absent |

### Exclusive trial design

`background_categories.json` is the source of truth. Its 37 categories are mutually exclusive:

- QCD: six generated-hard-pT intervals: `[15,50)`, `[50,100)`, `[100,200)`, `[200,400)`,
  `[400,800)`, `[800,infinity)`.
- W+jets: direct `e nu`, `mu nu`, or `tau nu`, crossed with LHE-HT intervals `[0,100)`,
  `[100,400)`, `[400,infinity)` (nine categories).
- Z/DY+jets: direct `ee`, `mumu`, `tautau`, or `nunu`, crossed with the same three LHE-HT
  intervals (twelve categories). The listed DY parent is restricted to `MLL>50`; do not mix in
  another mass slice under the same sample ID.
- ttbar: the unordered direct-W pairs `HH`, `HE`, `HM`, `HT`, `EE`, `EM`, `ET`, `MM`, `MT`, `TT`
  (ten categories). A tau is labeled at the W decay, before the tau itself decays.

There is no `other` category. A malformed or unclassifiable event throws and fails the job. The
recommended closure trial target is **5,000 retained events per fine category**. That is 185,000
events once all 37 categories are available, or 170,000 for the current 34-category set before the
pending Z-to-nunu parent is added. Start with a small pilot and scale each parent sample until its least-populated child
reaches 5,000; rare high-HT and `TT` children determine how many input files are needed. Do not
guess a fixed file count before measuring these efficiencies.

The physical staging is also exclusive: `stageBackgrounds.py` writes each event to exactly one
fine-category directory. Training may combine directories through `training_groups` (for example,
all six QCD directories become training label `QCD`) without changing or duplicating the original
fine label.

Normalize by parent sample, not child category:

```
event weight = luminosity * parent cross section * genWeight / parent sum(genWeight)
```

`Background_sample` ties every event to `background_samples_2025.json`. If several official
Z-to-nunu HT slices are needed, each slice gets a separate sample ID, cross section, and sum of
weights. Merge categories only after applying those parent weights.

---

## 0. One-time setup: accounts, grid certificate, storage

You need all of the following before CRAB or Rucio will work. Steps 1–3 can take a few days
because of approvals, so start them early.

1. **CERN computing account and CMS membership.** Your supervisor arranges this.
2. **Grid user certificate.**
   - Go to <https://ca.cern.ch/ca/>, choose *New Grid User Certificate*, set a password and
     download the `.p12` file.
   - Import the `.p12` into your web browser too; you need it for the VO registration below
     and for CMS web services.
   - Install it on lxplus:
     ```bash
     mkdir -p ~/.globus && cd ~/.globus
     # copy myCertificate.p12 here, then:
     openssl pkcs12 -in myCertificate.p12 -clcerts -nokeys -out usercert.pem
     openssl pkcs12 -in myCertificate.p12 -nocerts -out userkey.pem   # choose a PEM pass phrase
     chmod 644 usercert.pem
     chmod 400 userkey.pem
     ```
     If `openssl` complains about an unsupported algorithm, add `-legacy` to both commands.
   - The certificate is valid for one year. Renew it at the same page before it expires.
3. **CMS VO membership.** With the certificate loaded in your browser, register at
   <https://voms-cms-auth.app.cern.ch> and wait for the approval email.
4. **Check that the grid tools recognise you:**
   ```bash
   voms-proxy-init --voms cms --valid 192:00   # asks for the PEM pass phrase
   voms-proxy-info --all                       # should show the /cms VO attribute
   source /cvmfs/cms.cern.ch/common/crab-setup.sh
   crab checkusername                          # should print your CERN username
   ```
5. **Somewhere to write output.** CRAB copies the output to `/store/user/<username>/` at a grid site
   (the "storage site"). Ask your supervisor which site to use. The existing samples are at
   `T2_US_UCSD`, which needs your certificate DN registered by the UCSD site admins. Check that it
   works:
   ```bash
   crab checkwrite --site=T2_US_UCSD     # or the site you were given
   ```
6. **Rucio account** (only for staging data from tape, Section 4). Your Rucio account is your
   CERN username:
   ```bash
   source /cvmfs/cms.cern.ch/rucio/setup-py3.sh
   export RUCIO_ACCOUNT=$USER
   rucio whoami
   rucio list-account-limits $RUCIO_ACCOUNT   # where you have quota to stage data
   ```

You need a valid proxy (`voms-proxy-init ...`) in every new session that uses CRAB, Rucio,
`dasgoclient` or xrootd.

---

## 1. Set up the release

The branch is based on `CMSSW_16_1_0_pre4` (architecture `el8_amd64_gcc13`), so use exactly that
release. Build on an EL8 machine (`lxplus8.cern.ch`). On EL9 lxplus, open an EL8 container first
with `cmssw-el8`. A build made on EL9 for this EL8 release fails on the grid with GLIBC errors.

```bash
ssh <username>@lxplus8.cern.ch          # or: ssh lxplus.cern.ch, then: cmssw-el8
export SCRAM_ARCH=el8_amd64_gcc13
cmsrel CMSSW_16_1_0_pre4
cd CMSSW_16_1_0_pre4/src
cmsenv
git cms-init
git cms-merge-topic mbileska:addaxoscore-review
scram b -j 8
cd DPGAnalysis/L1TNanoAOD/test/axoScoutingNano   # everything below runs from here
```

`git cms-merge-topic` checks out only the packages the branch changes (`DPGAnalysis/L1TNanoAOD`,
`L1Trigger/L1TGlobal`, and `L1Trigger/L1TCaloLayer1`). Don't run `scram b clean` after making local
edits unless you know what you changed.

---

## 2. Which conditions for which year

**The conditions must match the year the MC was simulated for.** The raw data in a `GEN-SIM-RAW`
file encodes that year's detector, so you need the matching era and global tag (GT). Using the 2025
GT on 2023 MC, for example, applies 2025 calibrations and alignment to a 2023 detector, which
shifts the L1 inputs and therefore the AXO score.

The year is in the dataset name: the *campaign* is the start of the second field, and the
processed string tells you the conditions it was made with.

| `makeConfig.sh` year | Input campaign (dataset name contains) | `--era` | `--conditions` (global tag) | Tested |
|---|---|---|---|---|
| `2022preEE` | `Run3Summer22DRPremix` / `mcRun3_2022_realistic_v` | `Run3` | `150X_mcRun3_2022_realistic_v2` | legacy AXO path only; new closure rerun required |
| `2022postEE` | `Run3Summer22EEDRPremix` / `mcRun3_2022_realistic_postEE` | `Run3` | `150X_mcRun3_2022_realistic_postEE_v2` | new closure rerun required |
| `2023preBPix` | `Run3Summer23DRPremix` / `mcRun3_2023_realistic_v` | `Run3_2023` | `150X_mcRun3_2023_realistic_v2` | new closure rerun required |
| `2023postBPix` | `Run3Summer23BPixDRPremix` / `mcRun3_2023_realistic_postBPix` | `Run3_2023` | `150X_mcRun3_2023_realistic_postBPix_v2` | legacy AXO path only; new closure rerun required |
| `2024` | `Run3Winter24Digi` / `mcRun3_2024_realistic` | `Run3_2024` | `150X_mcRun3_2024_realistic_v3` | legacy AXO path only; new closure rerun required |
| `2025` | `Run3Winter25Digi` / `mcRun3_2025_realistic` | `Run3_2025` | `150X_mcRun3_2025_realistic_v12` | legacy AXO path only; new closure rerun required |

The legacy checks established that the era/GT choices and official AXO path run. They did not test
the new independent AXO/CICADA inference or exclusive classifier. Treat every row as unvalidated
for this implementation until Section 3.2 passes on a real input file; never start a large
production from the legacy status alone.

Other campaigns (for example `Run3Winter23Digi`, made with special conditions) are not covered.
Ask before processing them.

The GTs for 2022–2024 are the release's own `autoCond` entries (`auto:phase1_<year>_realistic...`),
and were checked with `conddb list <GT>`. For 2025 we keep `150X_mcRun3_2025_realistic_v12` so new
samples match the ones already produced (the release's `auto:phase1_2025_realistic` is
`151X_mcRun3_2025_realistic_v6`).

The rest is the same for every year:

| Setting | Value | Why |
|---|---|---|
| L1 menu | `--custom_conditions "L1Menu_Collisions2025_v1_2_0_xml,L1TUtmTriggerMenuRcd,,,9999-12-31 23:59:59.000"` | Use the 2025 L1 menu, and with it the 2025 AXO model, in every year. The 2022/2023 menus have no AXO; the 2024 menu has an older model with different seed names. Using one menu gives one score definition across all samples. |
| Steps | `-s L1REPACK:FullMC,HLT:GRun,NANO:@ScoutMonitor` | Re-emulate L1 from the raw detector data, re-run HLT, write scouting NanoAOD |
| Process name | `--process reL1reHLT` | New name, so the original `HLT` results stay readable as `*_pHLT` |
| Modifier `run3_scouting_2025` | added to the `cms.Process(...)` line for years before 2025 | The re-run HLT writes the 2025 scouting format in every year, so the scouting NanoAOD tables must be set up for that format. Without it they look for collections that don't exist and the job stops with `ProductNotFound`. |
| `customiseAXOScoutingNano` | added near the end of the config | Both official and independent anomaly evaluations, exact tensors, optional strict background table, `genWeight`, scouting PF candidate and track tables, `TryToContinue`, and GEM masking |
| Threads | `--nThreads 4` | Must match `--num-cores 4` in the CRAB submission |

`customiseAXOScoutingNano` also handles two things that otherwise stop the job:
- **`TryToContinue = ['FailModule']`.** The HLT menu in this release is newer than the 2025 L1 menu.
  Some HLT paths are seeded by L1 bits that don't exist in that menu, and their seed module throws
  a `FailModule` exception. This setting fails only those paths, not the whole job. Messages like
  `Algorithm L1_XYZ ... cannot be matched to a L1 algo name` in the job log are expected.
- **GEM strip masking off before 2025.** The HLT menu always masks GEM strips, but the 2022/2023
  GTs have no GEM masking conditions. Offline reconstruction only masks from 2025, so the HLT now
  does the same.

---

## 3. Make the config

### 3.1 Config for CRAB

```bash
./makeConfig.sh 2025 --background qcd   --sample-id 0
./makeConfig.sh 2025 --background wjets --sample-id 1
./makeConfig.sh 2025 --background zjets --sample-id 2
./makeConfig.sh 2025 --background ttbar --sample-id 3
```

These write `axoNano_2025_<process>_s<id>_cfg.py`, one explicitly typed config per parent sample.
The sample ID must match `background_samples_2025.json`; a second Z/DY parent must not reuse ID 2.
The script runs `cmsDriver.py` with
the settings above, then makes two edits to the generated file:
1. it adds `run3_scouting_2025` to the `cms.Process(...)` line (years before 2025 only),
2. it adds the `customiseAXOScoutingNano` call.

Both edits are needed. Make configs with this script rather than with your own `cmsDriver.py`
command; if you do use your own, the customisation stops the job with an error that explains what
is missing.

The input file in the config is a placeholder (`file:axoNano_<year>_RAW.root`), which CRAB replaces
with each job's input files. The output is named `out_nano.root`; CRAB uses that name to collect it.

### 3.2 Test locally first

Before submitting, run a few events of one real input file. The file must be on disk (see Section 4
if the dataset is only on tape).

```bash
voms-proxy-init --voms cms --valid 192:00
DATASET=$(head -1 samplelists/backgrounds_2025_wjets.txt)

# Do not assume the first DAS file has a disk replica. Probe a bounded set and
# keep the first file that the global redirector can actually locate.
LFN=
while IFS= read -r candidate; do
  echo "checking ${candidate}"
  if timeout 45 xrdfs root://cms-xrd-global.cern.ch stat "${candidate}" >/dev/null 2>&1; then
    LFN=${candidate}
    break
  fi
done < <(dasgoclient --limit=50 --query="file dataset=${DATASET}")

if [ -z "${LFN}" ]; then
  echo "No readable disk replica found among the first 50 files. Dataset sites:"
  dasgoclient --query="site dataset=${DATASET}" | sort -u
  echo "Stage a file or block as described in Section 4, then retry."
  false
fi

./makeConfig.sh 2025 10 root://cms-xrd-global.cern.ch/${LFN} --background wjets --sample-id 1
cmsRun test_2025_wjets_s1_cfg.py > test_2025_wjets_s1.log 2>&1; echo "exit code: $?"
python3 validateTrialNano.py test_2025_wjets_s1.root --output validation_wjets.json
```

Expected:
- exit code 0
- both score differences within `1e-6`, with all 57 AXO inputs and all 252 CICADA regions present
- every event has one valid `Background_category` and the configured `Background_sample == 1`
- `genWeight` present

The first few minutes go to loading the HLT menu, so 10 events take roughly 5–10 minutes.
An xrootd `3011 No servers are available to read the file` error occurs before event processing and
means that the selected LFN has no reachable disk replica; it is not an AXO/CICADA failure. Check
the individual file with `dasgoclient --query="site file=${LFN}"`, then select another disk-resident
file or stage one as described below.

---

## 4. Input data on tape: stage it with Rucio

Most `GEN-SIM-RAW` lives only on tape. Check where a dataset is:

```bash
dasgoclient --query="site dataset=/X/Y/GEN-SIM-RAW"
```

If it only lists `*_Tape` sites, copy it to disk with a Rucio rule before submitting. CRAB's
automatic tape recall refuses large datasets. Request the copy at a site where you have Rucio quota
(Section 0, step 6), ideally the site where you want the jobs to run.

```bash
source /cvmfs/cms.cern.ch/rucio/setup-py3.sh
export RUCIO_ACCOUNT=$USER
SITE=T2_US_UCSD                          # where you have quota
LIFETIME=2592000                         # 30 days; the disk copy is deleted afterwards

# Small dataset: stage all of it
rucio add-rule cms:/X/Y/GEN-SIM-RAW 1 ${SITE} --lifetime ${LIFETIME}

# Large dataset: stage only some blocks (one block name per line in the file)
dasgoclient --query="block dataset=/X/Y/GEN-SIM-RAW" | head -10 > blocks.txt
while read -r block; do
  rucio add-rule "cms:${block}" 1 ${SITE} --lifetime ${LIFETIME}
done < blocks.txt
```

The block file must have one block per line and end with a newline. A comma-separated file, or one
without a final newline, makes the loop skip blocks silently.

If a rule is refused because it exceeds your quota, re-run the same command with `--ask-approval`.
The request then goes to the site admins, who approve it or not.

Follow the transfer:

```bash
rucio list-rules --account $RUCIO_ACCOUNT        # one line per rule, with its state
rucio rule-info <rule_id>                         # details; state OK = everything is on disk
```

Recalls from tape take from hours to several days. Submit to CRAB only when the rules are `OK`, and
finish processing before the rule lifetime ends.

---

## 5. Submit with CRAB

`crab.py` builds a CRAB config for each line of an input list and submits it. Run it from this
directory, with your proxy and CRAB set up:

```bash
source /cvmfs/cms.cern.ch/common/crab-setup.sh
voms-proxy-init --voms cms --valid 192:00
```

Main options:

| Option | Meaning |
|---|---|
| `-p <cfg>` | The config from Section 3.1 |
| `-i <file>` | Input list: one DAS dataset per line, or a list of files with `--set-input-files` |
| `-t <tag>` | Output dataset tag (for `--set-input-files`, also the CRAB task name) |
| `-o /store/user/<username>/...` | Output directory at the storage site |
| `--site <site>` | Storage site (Section 0, step 5). The default `T3_US_FNALLPC` is probably not yours. |
| `--work-area <dir>` | Local directory for the CRAB task folders. **Use a separate one per year** (see below). |
| `-s FileBased -n N` | File-based splitting with N input files per job (`-n` is honored). Start with `-n 1` for the closure pilot. |
| `--num-cores 4 --max-memory 8000` | Matches `--nThreads 4`. L1+HLT needs this much memory. |
| `--no-publication` | Always use this |
| `--input-blocks "b1,b2,..."` | Process only these blocks (use with data you staged by block) |
| `--max-units N` | Process only N files *in total* |
| `--fnal` | Input files on FNAL EOS (`cmseos.fnal.gov`): run at US sites and read over xrootd |
| `--cern` | Input files on CERN EOS: run at `T2_CH_CERN` |
| `--ignore-locality` | Run at US T2s regardless of where the data is |
| `--dryrun` | Print the CRAB configs without submitting. Do this first. |
| `--status --work-area <dir>` | Status of every task in the work area; resubmits failed jobs |

For DAS datasets, the CRAB task name is built from the dataset's primary name (`-t` is not part
of it). Two datasets with the same primary name, such as the 2023 and 2024 versions of one process,
collide if you submit them into the same `--work-area` (`Working area ... already exists`). That's
why each year gets its own work area below.

A task can have at most 10,000 jobs, so at most 10,000 files. `samplesearch/file_counts.txt` has
the file count of each dataset.

The output ends up at
`<-o>/<PrimaryDataset>/<tag>_<campaign...>/<YYMMDD_HHMMSS>/0000/out_nano_<job>.root`.

### Example workflows

Replace `<username>` and `T2_US_UCSD` with your own. Add `--dryrun` the first time.

**2025 exclusive background closure pilot**
```bash
./makeConfig.sh 2025 --background wjets --sample-id 1
python3 crab.py -p axoNano_2025_wjets_s1_cfg.py -i samplelists/backgrounds_2025_wjets.txt \
  -t AnomalyNano_WJets_trial_v1 -o /store/user/<username>/anomalyNano/2025/WJets/ \
  --site T2_US_UCSD --work-area crab_2025_wjets_trial \
  -s FileBased -n 1 --max-units 10 --num-cores 4 --max-memory 8000 \
  --no-publication --dryrun
```

Repeat with the matching QCD, Z/DY, and ttbar config/list. Remove `--dryrun` only after inspecting
the generated CRAB configs. `crab.py` submits unless `--dryrun` is present.

**2025 flat-pT QCD, staged by block with Rucio (Section 4).** The dataset is about 30M events in
29,900 files (500 files per block), all on tape. Stage and process 10–20 blocks at a time. The 10
blocks in `samplelists/QCD_2025_blocks_v1.txt` are already processed, so take different ones:
```bash
DATASET=$(cat samplelists/QCD_2025_dataset.txt)
dasgoclient --query="block dataset=${DATASET}" | grep -v -F -f samplelists/QCD_2025_blocks_v1.txt \
  | head -10 > samplelists/QCD_2025_blocks_v2.txt
# ... add Rucio rules for these blocks (Section 4) and wait until they are OK ...
./makeConfig.sh 2025 --background qcd --sample-id 0
python3 crab.py -p axoNano_2025_qcd_s0_cfg.py -i samplelists/backgrounds_2025_qcd.txt \
  -t AXOnano_QCD_2025_v2 -o /store/user/<username>/axoNano/QCD_2025/ \
  --site T2_US_UCSD --work-area crab_QCD_2025_v2 \
  --input-blocks "$(paste -sd, samplelists/QCD_2025_blocks_v2.txt)" \
  -s FileBased -n 1 --num-cores 4 --max-memory 8000 --no-publication
```
Commit the new block list to the branch so the next person knows those blocks are done.

**2024 signals**
```bash
./makeConfig.sh 2024
python3 crab.py -p axoNano_2024_cfg.py -i samplelists/datasets_2024.txt \
  -t AXOnano_2024_v1 -o /store/user/<username>/axoNano/2024/ \
  --site T2_US_UCSD --work-area crab_2024 \
  -s FileBased --num-cores 4 --max-memory 8000 --no-publication
```

**2023 post-BPix signals**
```bash
./makeConfig.sh 2023postBPix
python3 crab.py -p axoNano_2023postBPix_cfg.py -i samplelists/datasets_2023postBPix.txt \
  -t AXOnano_2023postBPix_v1 -o /store/user/<username>/axoNano/2023postBPix/ \
  --site T2_US_UCSD --work-area crab_2023postBPix \
  -s FileBased --num-cores 4 --max-memory 8000 --no-publication
```

**2023 pre-BPix and 2022 (pre-/post-EE)**: same pattern with `./makeConfig.sh 2023preBPix`,
`2022preEE` or `2022postEE` and a dataset list you make for that campaign. A list must contain only
datasets of one campaign row of the table in Section 2.

**Private files on FNAL EOS (EMJ signal, 2022)**
```bash
./makeConfig.sh 2022preEE    # assumption: check with the sample producer which GT the DIGI step used
for list in samplelists/emj_private2022/*.txt; do
  name=$(basename ${list} .txt)_private2022
  python3 crab.py -p axoNano_2022preEE_cfg.py -i ${list} --set-input-files \
    -t ${name}_v1 -o /store/user/<username>/axoNano/${name}/ \
    --site T2_US_UCSD --work-area crab_emj \
    -s FileBased --num-cores 4 --max-memory 8000 --no-publication --fnal
done
```
For private samples the year can't be read from a dataset name. Ask whoever produced them which era
and GT they used, and pick the matching row of the table.

### Monitoring

```bash
crab status -d crab_2024/crab_<task name>
crab resubmit -d crab_2024/crab_<task name>              # retry failed jobs
python3 crab.py --status --work-area crab_2024           # all tasks; resubmits failures
python3 crab.py --status --no-resubmit --work-area crab_2024
```

### Validate and stage completed trial output

Run the validator on all outputs from one or more parents. It checks event identity within each
parent sample, exact tensor sizes, a complete CICADA grid, score closure, finite generator weights,
and agreement between the stored category fields and the JSON phase-space boundaries:

```bash
python3 validateTrialNano.py /path/to/out_nano_*.root \
  --expected-categories 0-23,27-36 \
  --min-events-per-populated-category 5000 --output trial_validation.json
```

The expected set above is the current 34-category trial. After adding and validating Z-to-nunu,
replace it with `--expected-categories all`. Add `--require-cross-sections` before any
physics-normalized staging; IDs 24--26 remain empty until the pending Z-to-nunu parent sample is
recorded in `background_samples_2025.json`. With `--min-events-per-populated-category`, the report
also gives each parent's least-populated child, scale factor, and estimated total output-file count.
That estimate is meaningful when the pilot used representative files and CRAB `FileBased -n 1`.

After validation, create physically separate fine-category shards:

```bash
python3 stageBackgrounds.py /path/to/out_nano_*.root \
  --validation-report trial_validation.json --output staged_backgrounds
```

Staging refuses a failed report, a different input list, or a nonempty output directory. It writes
`00_QCD_Pt15to50/`, ..., `36_TTbar_TT/` as applicable plus `manifest.json`. To train
one QCD class, read all category directories listed under `training_groups.QCD`; do not rewrite the
fine `Background_category` field.

---

## 6. Sample lists and search outputs

### `samplelists/`: inputs for `crab.py -i`

| File | Contents | How to use |
|---|---|---|
| `backgrounds_2025_qcd.txt` | Flat-pT QCD parent (sample ID 0) | Use only with the `qcd` config |
| `backgrounds_2025_wjets.txt` | Inclusive W+jets parent (sample ID 1) | Use only with the `wjets` config |
| `backgrounds_2025_zjets.txt` | `MLL>50` DY parent (sample ID 2); Z-to-nunu remains an explicit pending addition | Use only with the `zjets` config |
| `backgrounds_2025_ttbar.txt` | Inclusive ttbar parent (sample ID 3) | Use only with the `ttbar` config |
| `datasets_2025.txt` | 7 `Run3Winter25Digi` datasets: DY→2ℓ, TT, W+jets, VBF H→bb, HHH→4b2τ, gg→HH→2b2τ, VBF H→ττ | `./makeConfig.sh 2025`, then submit as a DAS-dataset list |
| `QCD_2025_dataset.txt` | The flat-pT QCD dataset `QCD_Bin-PT-15to7000_Par-PT-flat2022` (`Run3Winter25Digi`, about 30M events, on tape) | Kept separate because it has to be staged and submitted by block (Section 5) |
| `QCD_2025_blocks_v1.txt` | The 10 QCD blocks (5000 files) already processed, one block per line. Output at `/store/user/mequinna/ntuples/QCD_25/`. | Exclude these when choosing new blocks; add a `_v2`, `_v3`, … file for each new batch |
| `datasets_2024.txt` | 8 `Run3Winter24Digi` signal datasets (H→bb, H→γγ, H→ττ, VBF H→ττ, VBF H→cc, ttH→bb/cc, H→2 long-lived→4b) | `./makeConfig.sh 2024` |
| `datasets_2023postBPix.txt` | 20 `Run3Summer23BPixDRPremix` signal datasets (H→bb, H→γγ (90 GeV), VBF H→inv, long-lived H→4b mass/lifetime points, Higgsino, W→τ→3μ, SUSY bbH→bb mass points, HHH→6b, ttH, VBF H→cc, H→ττ) | `./makeConfig.sh 2023postBPix` |
| `emj_private2022/EMJ_s-channel_mMed-<M>_mDark-20_ctau-<cτ>.txt` | Emerging-jet signal: 12 mass/lifetime points (mMed 100, 200, 500, 750 GeV; cτ 1000, 1500, 2000 mm), privately produced 2022 DIGI files on FNAL EOS, one xrootd path per line | Submit with `--set-input-files --fnal`; check the GT with the producer (Section 5) |

Each `datasets_*.txt` file has one DAS dataset name per line, and all lines in a file are from one
campaign. `crab.py` makes one CRAB task per line.

### `samplesearch/`: how the datasets were found

The list of physics processes we want, signals and SM backgrounds, was searched in DAS: first in
the 2025 campaign, then in 2023 and 2024 for processes that had no 2025 sample.

| File | What it is |
|---|---|
| `samplesearch.sh` | Searches DAS for each process in the `Run3Winter25Digi` campaign. For every dataset found it writes `<process>_files.txt`, a list of all its files as xrootd paths, in the current directory. Those lists are large and not in the branch; regenerate them when needed. |
| `samplesearch_output.txt` | Saved output of `samplesearch.sh`: for each process, `Found: <dataset>` with the file count, or `NOT FOUND`. Most signals were not found in 2025, which is why the next two scripts exist. |
| `samplesearch23.sh`, `samplesearch24.sh` | For the processes missing in 2025, list every `GEN-SIM-RAW` dataset with `2023` or `2024` in its name |
| `output_2023.txt`, `output_2024.txt` | Saved outputs of those two scripts: one `=== <process> ===` header per process, followed by the matching datasets (nothing under a header means no sample). `output_2023.txt` can list several campaigns for one process, for example `Run3Winter23Digi` and `Run3Summer23BPixDRPremix`; the BPix ones went into `datasets_2023postBPix.txt`. |
| `file_counts.txt` | `<files> files \| <events> events \| <dataset>` for every dataset in the lists above, from DAS. Use it to estimate the number of CRAB jobs (one per file) and to stay under 10,000 jobs per task. |

To look for a new process, add a line to the `samples` (or `missing`) array of the right script,
run it with a valid proxy, and add the dataset you choose to the matching `samplelists/datasets_<year>.txt`.
Check its campaign against the table in Section 2.

---

## 7. Check the output

```bash
root -l "root://cms-xrd-global.cern.ch//store/user/<username>/axoNano/.../out_nano_1.root" \
  -e 'Events->Scan("genWeight:L1AXO_score:L1_AXO_Tight", "", "", 10)'
```

- `L1AXO_score`: real values (not all `-999`). Events with `L1_AXO_Tight = 1` should have high scores.
- `genWeight`: present. For the flat-pT QCD sample it spans many orders of magnitude.
- Normalize each `Background_sample` independently using its cross section and summed
  `genEventSumw`; only then merge its fine categories into a training group. The validator checks
  the `Runs` sum against the event-level `genWeight` sum and reports readiness.

Existing samples are under `root://redirector.t2.ucsd.edu:1095//store/user/mequinna/ntuples/`.
Those made before this README existed used the 2025 conditions for every input year, so 2023 and
2024 samples there are not equivalent to new ones made with the table in Section 2.

---

## 8. Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `SUBMITREFUSED` mentioning tape or dataset size | The dataset is on tape. Stage it with Rucio (Section 4) and submit with `--input-blocks`. |
| `Working area ... already exists` | Two datasets with the same primary name in one `--work-area`. Use a separate work area per year. |
| Job log full of `Algorithm L1_... cannot be matched to a L1 algo name` | Expected. These HLT paths are seeded by L1 bits not in the 2025 menu; `TryToContinue` skips them. |
| `RuntimeError: customiseAXOScoutingNano: the process must be built with the run3_scouting_2025 modifier` | The config wasn't made with `makeConfig.sh`. Make it with the script. |
| `ProductNotFound` for `Run3ScoutingMuon` (or other scouting collections) | Same cause: `run3_scouting_2025` is missing from the `cms.Process(...)` line. |
| `No "GEMMaskedStripsRcd" record found` | `customiseAXOScoutingNano` wasn't applied. Make the config with `makeConfig.sh`. |
| `GLIBC_2.xx not found` on the grid | The area was built on EL9. Rebuild on lxplus8 or inside `cmssw-el8`. |
| `L1AXO_score` always `-999` | The L1TGlobal fixes aren't compiled in. Check that `git cms-merge-topic` checked out `L1Trigger/L1TGlobal`, then run `scram b`. |
| Proxy or `checkwrite` errors | Renew the proxy (`voms-proxy-init --voms cms --valid 192:00`); check the certificate expiry date at <https://ca.cern.ch/ca/>. |
