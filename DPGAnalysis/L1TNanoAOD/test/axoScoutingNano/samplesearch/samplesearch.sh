# SM backgrounds
declare -A samples
samples["DY"]="/DYto2L-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter25Digi-142X_mcRun3_2025_realistic_v7*/GEN-SIM-RAW"
samples["TT"]="/TT_TuneCP5_13p6TeV_powheg-pythia8/Run3Winter25Digi-142X_mcRun3_2025_realistic_v7*/GEN-SIM-RAW"
samples["WJets"]="/WJetsToLNu_TuneCP5_13p6TeV_madgraphMLM-pythia8/Run3Winter25Digi-142X_mcRun3_2025_realistic_v7*/GEN-SIM-RAW"

# BSM signals
samples["GluGluHToBB"]="/GluGluHToBB_M-125*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["GluGluHToGG_125"]="/GluGluHToGG_M-125*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["GluGluHToGG_90"]="/GluGluHToGG_M-90*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["GluGluHToTauTau"]="/GluGluHToTauTau*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["GluGluHHto2B2W"]="/GluGluHHto2B2WtoLNu2Q*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["HHHTo6B"]="/HHHTo6B*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["HHHto4B2Tau"]="/HHHto4B2Tau*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["HTo2LL4b_1000_450_100000"]="/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["HTo2LL4b_1000_450_10000"]="/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["HTo2LL4b_125_12_900"]="/HTo2LongLivedTo4b_MH-125_MFF-12_CTau-900mm*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["HTo2LL4b_125_25_1500"]="/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-1500mm*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["HTo2LL4b_125_50_3000"]="/HTo2LongLivedTo4b_MH-125_MFF-50_CTau-3000mm*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["SMS_Higgsino"]="/SMS-Higgsino*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["SUSYGGtoBBH_120"]="/SUSYGluGluToBBHToBB_NarrowWidth_M-120*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["SUSYGGtoBBH_1200"]="/SUSYGluGluToBBHToBB_NarrowWidth_M-1200*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["SUSYGGtoBBH_350"]="/SUSYGluGluToBBHToBB_NarrowWidth_M-350*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["SUSYGGtoBBH_600"]="/SUSYGluGluToBBHToBB_NarrowWidth_M-600*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["VBFHToCC"]="/VBFHToCC*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["VBFHToInvisible"]="/VBFHToInvisible*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["VBFHToTauTau"]="/VBFHToTauTau*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["VBFHto2B"]="/VBFHto2B*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["WToTauTo3Mu"]="/WToTauTo3Mu*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["ggXToJpsiJpsi"]="/ggXToJpsiJpsiTo2Mu2E_m7*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["ggXToYY_m18"]="/ggXToYYTo2Mu2E_m18*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["ggXToYY_m26"]="/ggXToYYTo2Mu2E_m26*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["ttHto2B"]="/ttHto2B*/Run3Winter25Digi*/GEN-SIM-RAW"
samples["ttHto2C"]="/ttHto2C*/Run3Winter25Digi*/GEN-SIM-RAW"

for name in "${!samples[@]}"; do
    query="${samples[$name]}"
    echo "Searching for $name: $query"
    # First find the dataset
    dataset=$(dasgoclient --query="dataset=$query" 2>/dev/null | head -1)
    if [ -z "$dataset" ]; then
        echo "  NOT FOUND"
    else
        echo "  Found: $dataset"
        dasgoclient --query="file dataset=$dataset" 2>/dev/null \
            | sed 's|^|root://cms-xrd-global.cern.ch/|' > ${name}_files.txt
        echo "  Files: $(wc -l < ${name}_files.txt)"
    fi
done
