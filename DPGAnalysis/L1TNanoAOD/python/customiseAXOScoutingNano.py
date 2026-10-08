import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Modifier_run3_GEM_2025_cff import run3_GEM_2025
from Configuration.Eras.Modifier_run3_scouting_2025_cff import run3_scouting_2025
from PhysicsTools.NanoAOD.custom_run3scouting_cff import addScoutingParticle, addScoutingTrack


def _missingScouting2025Modifier(process):
    raise RuntimeError(
        "customiseAXOScoutingNano: the process must be built with the run3_scouting_2025 modifier, because the "
        "re-run HLT menu writes the 2025 scouting format. Generate the config with "
        "DPGAnalysis/L1TNanoAOD/test/axoScoutingNano/makeConfig.sh, which adds it for pre-2025 eras.")


def customiseAXOScoutingNano(process):
    """Add the emulated AXOL1TL score and generator weights to a re-L1 + re-HLT scouting NanoAOD job.

    Apply it after the scouting NanoAOD customisations, to a process built with the run3_scouting_2025
    modifier. DPGAnalysis/L1TNanoAOD/test/axoScoutingNano/makeConfig.sh takes care of both.
    """
    # The scouting NanoAOD tables are configured per era for that year's HLT scouting format, but the HLT
    # we re-run always writes the 2025 format; without run3_scouting_2025 the tables look for missing products.
    (~run3_scouting_2025).toModify(process, _missingScouting2025Modifier)

    # scouting PF candidates (ScoutingParticle) and tracks (ScoutingTrack), which @ScoutMonitor does not include
    process = addScoutingTrack(process)
    process = addScoutingParticle(process)

    # emulate the AXOL1TL score in the uGT emulator and store it as a NanoAOD column (L1AXO_score);
    # the parameter is not in the simGtStage2Digis cfi, so it has to be created with an explicit type
    process.simGtStage2Digis.produceAXOL1TLScore = cms.bool(True)
    process.load('DPGAnalysis.L1TNanoAOD.l1tnanotables_cff')
    process.nanoAOD_step += process.l1EmulAXOTable

    # generator weights; the module label ends in "Table", so the NANOAODoutput wildcard already keeps it
    process.load('PhysicsTools.NanoAOD.genWeightsTable_cfi')
    process.nanoAOD_step += process.genWeightsTable

    # HLTL1TSeed throws "FailModule" when an HLT path is seeded by an L1 bit that is not in the L1 menu.
    # The HLT menu in this release is newer than the L1 menu we emulate, so fail only that path, not the job.
    process.options.TryToContinue = ['FailModule']

    # The HLT menu always applies GEM strip masking, but the 2022/2023 MC global tags have no GEM masking
    # conditions. Follow the offline gemRecHits configuration instead: masking only from 2025 on.
    if hasattr(process, 'hltGemRecHits'):
        (~run3_GEM_2025).toModify(process.hltGemRecHits, applyMasking=False)

    return process
