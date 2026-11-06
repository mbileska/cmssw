#!/usr/bin/env bash
set -euo pipefail

if (($# != 5)); then
  echo "Usage: $0 INDEX LFN DATASET MAX_EVENTS CMSSW_SRC" >&2
  exit 2
fi

printf -v shard '%04d' "$1"
lfn="$2"
dataset="$3"
max_events="$4"
cmssw_src="$5"
package="$cmssw_src/L1Trigger/L1CaloTrigger"
analyzer="$cmssw_src/L1Trigger/L1CaloPhase2Analyzer"
pfn="root://cmsxrootd.fnal.gov/${lfn}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "$cmssw_src"
eval "$(scram runtime -sh)"
cd "${_CONDOR_SCRATCH_DIR:?Run this on a Condor worker}"

grep -qm1 avx2 /proc/cpuinfo
voms-proxy-info -exists -valid 0:30
test -f "$package/test/runFullCaloMHHGCTS_cfg.py"
test -f "$analyzer/python/mhhGCTSFullCaloTrainingNtuplizer_cfi.py"

{
  printf 'source_lfn=%s\nsource_pfn=%s\ndataset=%s\n' "$lfn" "$pfn" "$dataset"
  printf 'host=%s\nstarted_utc=%s\n' "$(hostname -f)" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'cmssw_commit=%s\nanalyzer_commit=%s\n' \
    "$(git -C "$package" rev-parse HEAD)" "$(git -C "$analyzer" rev-parse HEAD)"
  sha256sum "$package/test/runFullCaloMHHGCTS_cfg.py" \
    "$package/test/runFullCaloMHHGCTSTraining_cfg.py"
} > "provenance_${shard}.txt"

edmDumpEventContent "$pfn" > "input_content_${shard}.txt"
for product in simEcalEBTriggerPrimitiveDigis simHcalTriggerPrimitiveDigis \
  l1tHGCalTowerProducer l1tHGCalBackEndLayer2Producer; do
  grep -q "$product" "input_content_${shard}.txt" || {
    echo "Missing input product: $product" >&2
    exit 12
  }
done

requested_edm="fullcalo_${shard}_edm.root"
requested_training="fullcalo_${shard}_training.root"
cmsRun "$package/test/runFullCaloMHHGCTS_cfg.py" \
  inputFiles="$pfn" maxEvents="$max_events" outputFile="file:$requested_edm" \
  > "emulator_${shard}.log" 2>&1

# VarParsing decorates finite-maxEvents output names. Normalize after cmsRun.
actual_edm="$(find . -maxdepth 1 -type f -name "fullcalo_${shard}_edm*.root" -size +100k -print -quit)"
test -n "$actual_edm" || { echo 'Replay EDM is absent or too small' >&2; exit 20; }
if [[ "$actual_edm" != "./$requested_edm" ]]; then
  mv "$actual_edm" "$requested_edm"
fi
edmFileUtil -e "$requested_edm" > "edm_events_${shard}.txt"
edmDumpEventContent "$requested_edm" > "edm_content_${shard}.txt"
for product in l1tPhase2RCTEmulatorProducer l1tPhase2GCTEmulatorProducer \
  l1tPhase2L1GCTB1ToGCTS l1tPhase2L1GCTB2ToGCTS \
  l1tPhase2L1GCTB3ToGCTS l1tPhase2L1GCTSumFullCalo; do
  grep -q "$product" "edm_content_${shard}.txt" || {
    echo "Missing output product: $product" >&2
    exit 21
  }
done

cmsRun "$package/test/runFullCaloMHHGCTSTraining_cfg.py" \
  inputFiles="file:$PWD/$requested_edm" maxEvents=-1 \
  ntupleFile="$requested_training" dataset="$dataset" sourceLFN="$lfn" \
  cmsswCommit="$(git -C "$package" rev-parse HEAD)" \
  analyzerCommit="$(git -C "$analyzer" rev-parse HEAD)" \
  > "training_${shard}.log" 2>&1

python3 "$analyzer/test/check_mhh_gcts_training_ntuple.py" \
  "$requested_training" --expected-edm "$requested_edm" --expect-full-gctb --require-truth \
  > "analyzer_validation_${shard}.txt"
python3 "$package/test/checkFullCaloTraining.py" "$requested_training" \
  --min-events 1 > "fullcalo_validation_${shard}.txt"

sha256sum "$requested_edm" "$requested_training" \
  "emulator_${shard}.log" "training_${shard}.log" \
  "analyzer_validation_${shard}.txt" "fullcalo_validation_${shard}.txt" \
  > "SHA256SUMS_${shard}"
touch "SUCCESS_${shard}"
echo "PASS: full-calorimeter shard $shard"
