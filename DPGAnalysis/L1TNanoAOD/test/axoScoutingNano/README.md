# Scouting NanoAOD with the AXOL1TL score

This directory holds what you need to produce HLT-scouting NanoAOD samples that include the
AXOL1TL (AXO) L1 anomaly-detection score. Starting from MC `GEN-SIM-RAW`, each job:

1. re-runs the L1 trigger emulator (`L1REPACK:FullMC`) with the 2025 L1 menu and saves the AXO score,
2. re-runs the HLT (`HLT:GRun` menu of this release),
3. writes scouting NanoAOD (`NANO:@ScoutMonitor`).

Every input event is written out, whether or not it fired a scouting trigger. The main columns:

| Branch | Content |
|---|---|
| `L1AXO_score` | Continuous AXOL1TL score at BX 0 (`-999` only if the score could not be read) |
| `L1_AXO_VLoose` … `L1_AXO_VVVTight` | Re-emulated L1 AXO trigger bits |
| `DST_PFScouting_AXO*`, other `DST_*`/`HLT_*` | Re-run HLT decisions |
| `*_pHLT` versions of the trigger bits | HLT decisions from the original MC production, kept for comparison. Paths that exist only in the original menu keep their plain name (for example `DST_PFScouting_AXONominal`); the branch description says `(process: HLT)` for original bits and `(process: reL1reHLT)` for re-run ones. |
| `ScoutingParticle_*` (PF candidates), `ScoutingTrack_*`, `ScoutingPFJet_*`, `ScoutingMuonVtx_*`/`ScoutingMuonNoVtx_*`, `ScoutingElectron_*`, `ScoutingPhoton_*`, `ScoutingPrimaryVertex_*`, `ScoutingMET_*`, `ScoutingRho_*`, … | Scouting objects |
| `genWeight` (and `genEventSumw` in the `Runs` tree) | Generator weights for normalisation |

`nano_branches.txt` lists every branch with its type and description, from a 2025 test job. The
re-run part is the same for every year. The original-HLT bits depend on which menu the input MC
was made with.

---

## Contents of this directory

| Path | What it is |
|---|---|
| `README.md` | This file |
| `makeConfig.sh` | Writes the cmsRun config for one year (Section 3) |
| `crab.py` | CRAB submission wrapper (Section 5) |
| `nano_branches.txt` | All branches of the output NanoAOD |
| `samplelists/` | Input lists for CRAB (Section 6) |
| `samplesearch/` | Scripts used to find the input datasets in DAS, and their outputs (Section 6) |

The code changes this workflow needs live in two other packages of the branch:

| Path (under `$CMSSW_BASE/src`) | What it does |
|---|---|
| `L1Trigger/L1TGlobal/src/GlobalBoard.cc` | Fixes so the AXO score is actually stored. The AXO condition name is now set (an inverted `empty()` check prevented it), and it is no longer cleared every event. |
| `L1Trigger/L1TGlobal/plugins/L1TGlobalProducer.cc` | Turns on AXO score saving in `beginRun`, before the trigger conditions are built, and gives the score collection a real BX range. |
| `DPGAnalysis/L1TNanoAOD/plugins/AXOL1TLScoreTableProducer.cc` | Writes the BX-0 score as the NanoAOD column `L1AXO_score` |
| `DPGAnalysis/L1TNanoAOD/python/l1tnanotables_cff.py` | Defines the `l1EmulAXOTable` module |
| `DPGAnalysis/L1TNanoAOD/python/customiseAXOScoutingNano.py` | The customisation `makeConfig.sh` adds to every config (Section 2) |

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
git cms-merge-topic quinnanm:addaxoscore   # pulls branch addaxoscore from github.com/quinnanm/cmssw
scram b -j 8
cd DPGAnalysis/L1TNanoAOD/test/axoScoutingNano   # everything below runs from here
```

`git cms-merge-topic` checks out only the packages the branch changes (`DPGAnalysis/L1TNanoAOD`
and `L1Trigger/L1TGlobal`). Don't run `scram b clean` after making local edits unless you know what
you changed.

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
| `2022preEE` | `Run3Summer22DRPremix` / `mcRun3_2022_realistic_v` | `Run3` | `150X_mcRun3_2022_realistic_v2` | yes |
| `2022postEE` | `Run3Summer22EEDRPremix` / `mcRun3_2022_realistic_postEE` | `Run3` | `150X_mcRun3_2022_realistic_postEE_v2` | no (same setup as pre-EE, different GT) |
| `2023preBPix` | `Run3Summer23DRPremix` / `mcRun3_2023_realistic_v` | `Run3_2023` | `150X_mcRun3_2023_realistic_v2` | no (same setup as post-BPix, different GT) |
| `2023postBPix` | `Run3Summer23BPixDRPremix` / `mcRun3_2023_realistic_postBPix` | `Run3_2023` | `150X_mcRun3_2023_realistic_postBPix_v2` | yes |
| `2024` | `Run3Winter24Digi` / `mcRun3_2024_realistic` | `Run3_2024` | `150X_mcRun3_2024_realistic_v3` | yes |
| `2025` | `Run3Winter25Digi` / `mcRun3_2025_realistic` | `Run3_2025` | `150X_mcRun3_2025_realistic_v12` | yes, and used for the existing samples |

What "tested" means: the config was made with `makeConfig.sh` in this release, ran on a few events
with exit code 0, and wrote every event with `L1AXO_score` filled and `genWeight` present. The 2025
check used a central `Run3Winter25Digi` file. The existing 2025 production also shows high scores
lining up with `L1_AXO_Tight`. The 2022, 2023 and 2024 checks used TTbar events generated with that
year's GT, because no central `GEN-SIM-RAW` from those years was on disk at the time. In those, the
high-score events fired the AXO L1 bits and `DST_PFScouting_AXO*` paths. Before a large production of a row, especially an untested one, run the local test
in Section 3.2 on a real input file.

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
| `customiseAXOScoutingNano` | added near the end of the config | AXO score on, `L1AXO` table, `genWeight` table, scouting PF candidate and track tables (not part of `@ScoutMonitor`), `TryToContinue`, GEM masking (see below) |
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
./makeConfig.sh 2023postBPix        # or any year from the table in Section 2
```

This writes `axoNano_2023postBPix_cfg.py`, the file you give to CRAB. It runs `cmsDriver.py` with
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
DATASET=$(head -1 samplelists/datasets_2023postBPix.txt)     # any dataset of your year
LFN=$(dasgoclient --query="file dataset=${DATASET}" | head -1)

./makeConfig.sh 2023postBPix 10 root://cms-xrd-global.cern.ch/${LFN}   # writes test_2023postBPix_cfg.py
cmsRun test_2023postBPix_cfg.py > test_2023postBPix.log 2>&1; echo "exit code: $?"
root -l -b -q test_2023postBPix.root \
  -e 'Events->Scan("L1AXO_score:L1_AXO_VLoose:L1_AXO_Tight:DST_PFScouting_AXOVLoose:genWeight")'
```

Expected:
- exit code 0
- `L1AXO_score` real numbers (not `-999`), with high scores tending to have `L1_AXO_Tight = 1`
- `genWeight` present

The first few minutes go to loading the HLT menu, so 10 events take roughly 5–10 minutes.

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
| `-s FileBased` | One input file per job. `-n` is ignored: `crab.py` always uses 1 file per job. |
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

**2025 signals and backgrounds**
```bash
./makeConfig.sh 2025
python3 crab.py -p axoNano_2025_cfg.py -i samplelists/datasets_2025.txt \
  -t AXOnano_2025_v1 -o /store/user/<username>/axoNano/2025/ \
  --site T2_US_UCSD --work-area crab_2025 \
  -s FileBased --num-cores 4 --max-memory 8000 --no-publication
```

**2025 flat-pT QCD, staged by block with Rucio (Section 4).** The dataset is about 30M events in
29,900 files (500 files per block), all on tape. Stage and process 10–20 blocks at a time. The 10
blocks in `samplelists/QCD_2025_blocks_v1.txt` are already processed, so take different ones:
```bash
DATASET=$(cat samplelists/QCD_2025_dataset.txt)
dasgoclient --query="block dataset=${DATASET}" | grep -v -F -f samplelists/QCD_2025_blocks_v1.txt \
  | head -10 > samplelists/QCD_2025_blocks_v2.txt
# ... add Rucio rules for these blocks (Section 4) and wait until they are OK ...
./makeConfig.sh 2025
python3 crab.py -p axoNano_2025_cfg.py -i samplelists/QCD_2025_dataset.txt \
  -t AXOnano_QCD_2025_v2 -o /store/user/<username>/axoNano/QCD_2025/ \
  --site T2_US_UCSD --work-area crab_QCD_2025_v2 \
  --input-blocks "$(paste -sd, samplelists/QCD_2025_blocks_v2.txt)" \
  -s FileBased --num-cores 4 --max-memory 8000 --no-publication
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

---

## 6. Sample lists and search outputs

### `samplelists/`: inputs for `crab.py -i`

| File | Contents | How to use |
|---|---|---|
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
- To normalise, sum `genEventSumw` from the `Runs` tree of all output files.

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
