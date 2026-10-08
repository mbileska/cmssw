#!/bin/bash
# Make the cmsRun config for MC of one year (see README.md, section 3).
#
#   ./makeConfig.sh <year>                        -> axoNano_<year>_cfg.py for CRAB (all events, CRAB sets the inputs)
#   ./makeConfig.sh <year> <nEvents> <inputFile>  -> test_<year>_cfg.py for a local test, writes test_<year>.root
#
#   <year>: 2022preEE 2022postEE 2023preBPix 2023postBPix 2024 2025
set -e

YEAR=$1
NEVENTS=${2:--1}
INPUT=${3:-}

# era and global tag matching the detector conditions the MC was simulated with
case "${YEAR}" in
  2022preEE)    ERA=Run3;      GT=150X_mcRun3_2022_realistic_v2 ;;
  2022postEE)   ERA=Run3;      GT=150X_mcRun3_2022_realistic_postEE_v2 ;;
  2023preBPix)  ERA=Run3_2023; GT=150X_mcRun3_2023_realistic_v2 ;;
  2023postBPix) ERA=Run3_2023; GT=150X_mcRun3_2023_realistic_postBPix_v2 ;;
  2024)         ERA=Run3_2024; GT=150X_mcRun3_2024_realistic_v3 ;;
  2025)         ERA=Run3_2025; GT=150X_mcRun3_2025_realistic_v12 ;;
  *) echo "usage: $0 <2022preEE|2022postEE|2023preBPix|2023postBPix|2024|2025> [nEvents inputFile]"; exit 1 ;;
esac

if [ -n "${INPUT}" ]; then
  NAME=test_${YEAR}
  OUTPUT=file:test_${YEAR}.root
else
  NAME=axoNano_${YEAR}
  OUTPUT=file:out_nano.root   # CRAB collects the output by this name
fi

# the 2025 L1 menu (and so the 2025 AXO model) is used for every year
cmsDriver.py ${NAME} \
  --era ${ERA} --conditions ${GT} \
  --custom_conditions "L1Menu_Collisions2025_v1_2_0_xml,L1TUtmTriggerMenuRcd,,,9999-12-31 23:59:59.000" \
  -s L1REPACK:FullMC,HLT:GRun,NANO:@ScoutMonitor \
  --mc --scenario pp --process reL1reHLT \
  --datatier NANOAOD --eventcontent NANOAOD \
  --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run3 \
  --nThreads 4 -n ${NEVENTS} \
  --filein "${INPUT}" --fileout ${OUTPUT} \
  --python_filename ${NAME}_cfg.py --no_exec

# The re-run HLT writes the 2025 scouting format, so the scouting NanoAOD tables need the run3_scouting_2025
# modifier in every year. cmsDriver's --era does not accept it, so add it to the Process constructor.
if [ "${ERA}" != "Run3_2025" ]; then
  sed -i "s/^process = cms.Process('reL1reHLT',\(.*\))$/from Configuration.Eras.Modifier_run3_scouting_2025_cff import run3_scouting_2025\nprocess = cms.Process('reL1reHLT',\1,run3_scouting_2025)/" ${NAME}_cfg.py
  grep -q "^process = cms.Process('reL1reHLT',${ERA},run3_scouting_2025)$" ${NAME}_cfg.py \
    || { echo "ERROR: could not add run3_scouting_2025 to ${NAME}_cfg.py"; exit 1; }
fi

# Add the AXO score, genWeight and related settings after the scouting NanoAOD customisations. This is done here
# rather than with --customise_commands because cmsDriver would run it before the modifier above is added.
CUSTOMISATION="# AXO score, genWeight, TryToContinue and GEM masking (DPGAnalysis/L1TNanoAOD/python/customiseAXOScoutingNano.py)\nfrom DPGAnalysis.L1TNanoAOD.customiseAXOScoutingNano import customiseAXOScoutingNano\nprocess = customiseAXOScoutingNano(process)\n"
if grep -q "^# Add early deletion" ${NAME}_cfg.py; then
  sed -i "s|^# Add early deletion|${CUSTOMISATION}\n# Add early deletion|" ${NAME}_cfg.py
else
  printf "\n${CUSTOMISATION}" >> ${NAME}_cfg.py
fi
python3 ${NAME}_cfg.py > /dev/null || { echo "ERROR: ${NAME}_cfg.py does not load"; exit 1; }

echo "Wrote ${NAME}_cfg.py (era ${ERA}, global tag ${GT})"
