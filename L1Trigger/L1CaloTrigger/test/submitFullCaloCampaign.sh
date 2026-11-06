#!/usr/bin/env bash
set -euo pipefail

dataset=''
max_events=10
limit_files=1
tag="dist_poc_$(date +%Y%m%d_%H%M%S)"
while (($#)); do
  case "$1" in
    --dataset) dataset="$2"; shift 2 ;;
    --max-events) max_events="$2"; shift 2 ;;
    --limit-files) limit_files="$2"; shift 2 ;;
    --tag) tag="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

: "${CMSSW_BASE:?Initialize the CMSSW runtime first}"
[[ "$max_events" =~ ^(-1|[1-9][0-9]*)$ ]] || { echo 'Invalid --max-events' >&2; exit 2; }
[[ "$limit_files" =~ ^[0-9]+$ ]] || { echo 'Invalid --limit-files' >&2; exit 2; }
[[ "$tag" =~ ^[A-Za-z0-9._-]+$ ]] || { echo 'Invalid --tag' >&2; exit 2; }
test -n "$dataset" || { echo 'Supply --dataset' >&2; exit 2; }
command -v dasgoclient >/dev/null
command -v condor_submit >/dev/null
voms-proxy-info -exists -valid 1:00 >/dev/null

package="$CMSSW_BASE/src/L1Trigger/L1CaloTrigger"
analyzer="$CMSSW_BASE/src/L1Trigger/L1CaloPhase2Analyzer"
runner="$package/test/runFullCaloShard.sh"
test -x "$runner" || { echo "Missing executable $runner" >&2; exit 3; }
test -f "$analyzer/python/mhhGCTSFullCaloTrainingNtuplizer_cfi.py" || {
  echo 'Pinned mbileska analyzer checkout is missing' >&2
  exit 3
}
test "$(git -C "$package" branch --show-current)" = dist-trigger-poc || {
  echo 'Run from the dist-trigger-poc branch' >&2
  exit 3
}
test "$(git -C "$analyzer" rev-parse HEAD)" = \
  d3879c57c5c34d7ef47d15252f6317c2f9297f44 || {
  echo 'Analyzer is not at the reviewed revision; revalidate before submitting' >&2
  exit 3
}

campaign_dir="/scratch/$USER/$tag"
test ! -e "$campaign_dir" || { echo "Directory exists: $campaign_dir" >&2; exit 4; }
mkdir -p "$campaign_dir/logs"
chmod 700 "$campaign_dir" "$campaign_dir/logs"

dasgoclient --query="file dataset=$dataset" | sort -u > "$campaign_dir/all_lfns.txt"
test -s "$campaign_dir/all_lfns.txt" || { echo 'DAS returned no files' >&2; exit 5; }
if ((limit_files > 0)); then
  head -n "$limit_files" "$campaign_dir/all_lfns.txt" > "$campaign_dir/lfns.txt"
else
  cp "$campaign_dir/all_lfns.txt" "$campaign_dir/lfns.txt"
fi
awk '{print NR-1, $0}' "$campaign_dir/lfns.txt" > "$campaign_dir/jobs.txt"

cat > "$campaign_dir/submit.sub" <<SUBMIT
universe = vanilla
initialdir = $campaign_dir
executable = $runner
arguments = \$(job_index) \$(input_lfn) $dataset $max_events $CMSSW_BASE/src
output = logs/job_\$(job_index).out
error = logs/job_\$(job_index).err
log = logs/condor.log
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
use_x509userproxy = True
request_cpus = 1
request_memory = 4000MB
request_disk = 12000MB
requirements = (TARGET.has_avx2 =?= True) && regexp("hep[.]wisc[.]edu$", TARGET.Machine)
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
queue job_index,input_lfn from jobs.txt
SUBMIT

{
  printf 'dataset=%s\nmax_events_per_file=%s\nselected_files=%s\n' \
    "$dataset" "$max_events" "$(wc -l < "$campaign_dir/lfns.txt")"
  printf 'cmssw_commit=%s\nanalyzer_commit=%s\n' \
    "$(git -C "$package" rev-parse HEAD)" "$(git -C "$analyzer" rev-parse HEAD)"
  printf 'created_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  sha256sum "$package/test/runFullCaloMHHGCTS_cfg.py" \
    "$package/test/runFullCaloMHHGCTSTraining_cfg.py" "$runner"
} > "$campaign_dir/CAMPAIGN_PROVENANCE.txt"

(cd "$campaign_dir" && condor_submit submit.sub | tee submission.txt)
echo "Campaign directory: $campaign_dir"
echo "Monitor on this submit host: condor_q \"$USER\""
