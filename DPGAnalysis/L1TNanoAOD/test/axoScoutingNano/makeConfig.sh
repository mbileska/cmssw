#!/bin/bash
# Make the cmsRun config for MC of one year (see README.md, section 3).
#
#   ./makeConfig.sh <year> [--background <process> --sample-id <id>]
#   ./makeConfig.sh <year> <nEvents> <inputFile> [--background <process> --sample-id <id>]
#
#   <year>: 2022preEE 2022postEE 2023preBPix 2023postBPix 2024 2025
#   <process>: qcd wjets zjets ttbar (omit only for signal/configuration tests)
set -e

usage() {
  echo "usage: $0 <2022preEE|2022postEE|2023preBPix|2023postBPix|2024|2025> [nEvents inputFile] [--background qcd|wjets|zjets|ttbar --sample-id ID]"
}

[ $# -ge 1 ] || { usage; exit 1; }

YEAR=$1
shift
NEVENTS=-1
INPUT=
LOCAL_TEST=0
BACKGROUND=none
SAMPLE_ID=

# Keep the original positional local-test interface, then parse the new background option.
if [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; then
  NEVENTS=$1
  shift
  [ $# -gt 0 ] && [ "${1#--}" = "$1" ] || {
    echo "ERROR: a positional nEvents requires an inputFile"
    usage
    exit 1
  }
  INPUT=$1
  LOCAL_TEST=1
  shift
fi
[[ "${NEVENTS}" =~ ^-?[0-9]+$ ]] || { echo "ERROR: nEvents must be an integer"; exit 1; }
[ "${NEVENTS}" -ge -1 ] || { echo "ERROR: nEvents must be -1 or non-negative"; exit 1; }
while [ $# -gt 0 ]; do
  case "$1" in
    --background)
      [ $# -ge 2 ] || { echo "ERROR: --background needs a value"; exit 1; }
      BACKGROUND=$2
      shift 2
      ;;
    --sample-id)
      [ $# -ge 2 ] || { echo "ERROR: --sample-id needs a value"; exit 1; }
      SAMPLE_ID=$2
      shift 2
      ;;
    *) echo "ERROR: unknown option $1"; exit 1 ;;
  esac
done

case "${BACKGROUND}" in
  none|qcd|wjets|zjets|ttbar) ;;
  *) echo "ERROR: background must be one of: qcd wjets zjets ttbar"; exit 1 ;;
esac
if [ "${BACKGROUND}" != "none" ]; then
  [ -n "${SAMPLE_ID}" ] || { echo "ERROR: --sample-id is required for a background config"; exit 1; }
  [[ "${SAMPLE_ID}" =~ ^[0-9]+$ ]] || { echo "ERROR: --sample-id must be a non-negative integer"; exit 1; }
else
  [ -z "${SAMPLE_ID}" ] || { echo "ERROR: --sample-id requires --background"; exit 1; }
  SAMPLE_ID=0
fi

# era and global tag matching the detector conditions the MC was simulated with
case "${YEAR}" in
  2022preEE)    ERA=Run3;      GT=150X_mcRun3_2022_realistic_v2 ;;
  2022postEE)   ERA=Run3;      GT=150X_mcRun3_2022_realistic_postEE_v2 ;;
  2023preBPix)  ERA=Run3_2023; GT=150X_mcRun3_2023_realistic_v2 ;;
  2023postBPix) ERA=Run3_2023; GT=150X_mcRun3_2023_realistic_postBPix_v2 ;;
  2024)         ERA=Run3_2024; GT=150X_mcRun3_2024_realistic_v3 ;;
  2025)         ERA=Run3_2025; GT=150X_mcRun3_2025_realistic_v12 ;;
  *) usage; exit 1 ;;
esac

SUFFIX=
if [ "${BACKGROUND}" != "none" ]; then
  SUFFIX=_${BACKGROUND}_s${SAMPLE_ID}
fi
if [ "${LOCAL_TEST}" -eq 1 ]; then
  NAME=test_${YEAR}${SUFFIX}
  OUTPUT=file:test_${YEAR}${SUFFIX}.root
  FILEIN=${INPUT}
else
  NAME=axoNano_${YEAR}${SUFFIX}
  OUTPUT=file:out_nano.root   # CRAB collects the output by this name
  FILEIN=file:axoNano_${YEAR}_RAW.root  # CRAB replaces this placeholder
fi

# the 2025 L1 menu (and so the 2025 AXO model) is used for every year
cmsDriver.py "${NAME}" \
  --era "${ERA}" --conditions "${GT}" \
  --custom_conditions "L1Menu_Collisions2025_v1_2_0_xml,L1TUtmTriggerMenuRcd,,,9999-12-31 23:59:59.000" \
  -s L1REPACK:FullMC,HLT:GRun,NANO:@ScoutMonitor \
  --mc --scenario pp --process reL1reHLT \
  --datatier NANOAOD --eventcontent NANOAOD \
  --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run3 \
  --nThreads 4 -n "${NEVENTS}" \
  --filein "${FILEIN}" --fileout "${OUTPUT}" \
  --python_filename "${NAME}_cfg.py" --no_exec

# The re-run HLT writes the 2025 scouting format, so the scouting NanoAOD tables need the run3_scouting_2025
# modifier in every year. cmsDriver's --era does not accept it, so add it to the Process constructor.
if [ "${ERA}" != "Run3_2025" ]; then
  sed -i "s/^process = cms.Process('reL1reHLT',\(.*\))$/from Configuration.Eras.Modifier_run3_scouting_2025_cff import run3_scouting_2025\nprocess = cms.Process('reL1reHLT',\1,run3_scouting_2025)/" "${NAME}_cfg.py"
  grep -q "^process = cms.Process('reL1reHLT',${ERA},run3_scouting_2025)$" "${NAME}_cfg.py" \
    || { echo "ERROR: could not add run3_scouting_2025 to ${NAME}_cfg.py"; exit 1; }
fi

# Add the AXO score, genWeight and related settings after the scouting NanoAOD customisations. This is done here
# rather than with --customise_commands because cmsDriver would run it before the modifier above is added.
CUSTOMISATION="# AXO/CICADA scores, optional exclusive background label, genWeight and safeguards\nfrom DPGAnalysis.L1TNanoAOD.customiseAXOScoutingNano import customiseAXOScoutingNano\nprocess = customiseAXOScoutingNano(process, backgroundProcess='${BACKGROUND}', backgroundSample=${SAMPLE_ID})\n"
if grep -q "^# Add early deletion" "${NAME}_cfg.py"; then
  sed -i "s|^# Add early deletion|${CUSTOMISATION}\n# Add early deletion|" "${NAME}_cfg.py"
else
  printf "\n%b" "${CUSTOMISATION}" >> "${NAME}_cfg.py"
fi
python3 "${NAME}_cfg.py" > /dev/null || { echo "ERROR: ${NAME}_cfg.py does not load"; exit 1; }

echo "Wrote ${NAME}_cfg.py (era ${ERA}, global tag ${GT}, background ${BACKGROUND}, sample ${SAMPLE_ID})"
