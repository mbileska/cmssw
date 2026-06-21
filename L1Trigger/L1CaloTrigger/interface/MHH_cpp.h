#ifndef L1Trigger_L1CaloTrigger_MHH_cpp
#define L1Trigger_L1CaloTrigger_MHH_cpp

#include "L1Trigger/L1CaloTrigger/interface/MHH_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_MHH_h.h"

namespace p2mhh {

static inline ap_uint<48> pack_gamma_output(const MHHObject& object) {
    #pragma HLS INLINE
    ap_uint<48> out = 0;
    out = ((ap_uint<48>)object.isBarrel << 47) |
          ((ap_uint<48>)object.Phi << 19) |
          ((ap_uint<48>)object.Eta << 12) |
          (ap_uint<48>)object.ET;
    return out;
}

static inline ap_uint<48> pack_had_output(const MHHObject& object) {
    #pragma HLS INLINE
    ap_uint<48> out = 0;
    out = ((ap_uint<48>)object.isBarrel << 47) |
          ((ap_uint<48>)object.PtClusterSeed << 27) |
          ((ap_uint<48>)object.Phi << 18) |
          ((ap_uint<48>)object.Eta << 12) |
          (ap_uint<48>)object.ET;
    return out;
}

static inline ap_int<12> unpack_signed_sum_component(
    ap_uint<48> word,
    ap_uint<2> scenario
) {
    #pragma HLS INLINE
    switch (scenario) {
        case 0: return (ap_int<12>)word.range(11, 0);
        case 1: return (ap_int<12>)word.range(23, 12);
        default: return (ap_int<12>)word.range(35, 24);
    }
}

static inline ap_uint<12> unpack_unsigned_sum_component(ap_uint<48> word) {
    #pragma HLS INLINE
    return word.range(11, 0);
}

static inline ap_uint<48> pack_signed_sum_hypotheses(const ap_int<12> values[3]) {
    #pragma HLS INLINE
    ap_uint<48> out = 0;
    out.range(11, 0) = (ap_uint<12>)values[0];
    out.range(23, 12) = (ap_uint<12>)values[1];
    out.range(35, 24) = (ap_uint<12>)values[2];
    return out;
}

static inline ap_uint<48> pack_unsigned_sum_component(ap_uint<12> value) {
    #pragma HLS INLINE
    ap_uint<48> out = 0;
    out.range(11, 0) = value;
    return out;
}

static inline ap_uint<12> saturating_add(ap_uint<12> a, ap_uint<12> b) {
    #pragma HLS INLINE
    ap_uint<13> s = (ap_uint<13>)a + (ap_uint<13>)b;
    return (s > 0xFFF) ? (ap_uint<12>)0xFFF : (ap_uint<12>)s;
}

static const int MHH_N_ACTIVE_OBJECTS = 12;
static const int MHH_N_SORT_OBJECTS = 32;
static const int MHH_HF_LINK_BASE = 0;
static const int MHH_HGCAL_LINK_BASE = 8;
static const int MHH_CARD_A_OUTPUT_BASE = 0;
static const int MHH_CARD_B_OUTPUT_BASE = 3;

static const ap_uint<10> HGCAL_GAMMA_BOUNDARY_ETA = 84;
static const ap_uint<10> HF_GAMMA_BOUNDARY_ETA = 85;
static const ap_uint<10> HGCAL_HAD_BOUNDARY_ETA = 5;
static const ap_uint<10> HF_HAD_BOUNDARY_ETA = 6;

static const int N_PHI_REGIONS = 8;
static const int GAMMA_PHI_BINS = 360;
static const int GAMMA_REGION_WIDTH = 45;
static const int GAMMA_REGION_OVERLAP = 5;
static const int HAD_PHI_BINS = 24;
static const int HAD_BINS_PER_REGION = 3;
static const int HAD_REGION_OVERLAP = 1;

static inline bool match_dphi_gamma(ap_uint<9> hf_phi, ap_uint<9> hgcal_phi) {
    #pragma HLS INLINE
    ap_int<10> dphi = (ap_int<10>)hf_phi - (ap_int<10>)hgcal_phi;
    if (dphi > 180) dphi -= 360;
    if (dphi < -180) dphi += 360;
    return (dphi >= -1 && dphi <= 1);
}

static inline bool match_dphi_had(ap_uint<9> hf_phi, ap_uint<9> hgcal_phi) {
    #pragma HLS INLINE
    ap_int<10> dphi = (ap_int<10>)hf_phi - (ap_int<10>)hgcal_phi;
    if (dphi > 12) dphi -= 24;
    if (dphi < -12) dphi += 24;
    return (dphi >= -1 && dphi <= 1);
}

static inline void clear_gctvar(MHHObject& object) {
    #pragma HLS INLINE
    object.ET = 0;
    object.Eta = 0;
    object.Phi = 0;
    object.PtClusterSeed = 0;
    object.isBarrel = 0;
}

static inline bool is_active_boundary_flag(
    ap_uint<1> boundary_flag,
    ap_uint<12> et,
    ap_uint<1> cleared
) {
    #pragma HLS INLINE
    return (cleared == 0) && (et != 0) && (boundary_flag == 1);
}

static inline void stitch_state_pair_local(
    ap_uint<12>& hf_et,
    ap_uint<1>& hf_cleared,
    ap_uint<12>& hgcal_et,
    ap_uint<1>& hgcal_cleared
) {
    #pragma HLS INLINE
    ap_uint<12> sum = saturating_add(hf_et, hgcal_et);

    if (hf_et > hgcal_et) {
        hf_et = sum;
        hgcal_et = 0;
        hgcal_cleared = 1;
    } else {
        hgcal_et = sum;
        hf_et = 0;
        hf_cleared = 1;
    }
}

static inline ap_uint<9> normalize_gamma_phi(ap_uint<9> phi) {
    #pragma HLS INLINE
    return (phi >= GAMMA_PHI_BINS) ? (ap_uint<9>)(phi - GAMMA_PHI_BINS) : phi;
}

static inline ap_uint<5> normalize_had_phi(ap_uint<9> phi) {
    #pragma HLS INLINE
    ap_uint<9> p = phi;
    if (p >= 384) p -= 384;
    if (p >= 192) p -= 192;
    if (p >= 96) p -= 96;
    if (p >= 48) p -= 48;
    if (p >= HAD_PHI_BINS) p -= HAD_PHI_BINS;
    return (ap_uint<5>)p;
}

static inline ap_uint<3> gamma_owner_region(ap_uint<9> phi) {
    #pragma HLS INLINE
    ap_uint<9> p = normalize_gamma_phi(phi);
    if (p < (1 * GAMMA_REGION_WIDTH)) return 0;
    if (p < (2 * GAMMA_REGION_WIDTH)) return 1;
    if (p < (3 * GAMMA_REGION_WIDTH)) return 2;
    if (p < (4 * GAMMA_REGION_WIDTH)) return 3;
    if (p < (5 * GAMMA_REGION_WIDTH)) return 4;
    if (p < (6 * GAMMA_REGION_WIDTH)) return 5;
    if (p < (7 * GAMMA_REGION_WIDTH)) return 6;
    return 7;
}

static inline ap_uint<3> had_owner_region(ap_uint<9> phi) {
    #pragma HLS INLINE
    ap_uint<5> p = normalize_had_phi(phi);
    if (p < (1 * HAD_BINS_PER_REGION)) return 0;
    if (p < (2 * HAD_BINS_PER_REGION)) return 1;
    if (p < (3 * HAD_BINS_PER_REGION)) return 2;
    if (p < (4 * HAD_BINS_PER_REGION)) return 3;
    if (p < (5 * HAD_BINS_PER_REGION)) return 4;
    if (p < (6 * HAD_BINS_PER_REGION)) return 5;
    if (p < (7 * HAD_BINS_PER_REGION)) return 6;
    return 7;
}

static inline bool gamma_phi_in_region_overlap(ap_uint<9> phi, ap_uint<3> region) {
    #pragma HLS INLINE
    int start = (int(region) * GAMMA_REGION_WIDTH) - GAMMA_REGION_OVERLAP;
    int end = ((int(region) + 1) * GAMMA_REGION_WIDTH) - 1 + GAMMA_REGION_OVERLAP;
    ap_uint<9> p = normalize_gamma_phi(phi);

    if (start < 0) {
        return (p >= (start + GAMMA_PHI_BINS)) || (p <= end);
    }
    if (end >= GAMMA_PHI_BINS) {
        return (p >= start) || (p <= (end - GAMMA_PHI_BINS));
    }
    return (p >= start) && (p <= end);
}

static inline bool had_phi_in_region_overlap(ap_uint<9> phi, ap_uint<3> region) {
    #pragma HLS INLINE
    int start = (int(region) * HAD_BINS_PER_REGION) - HAD_REGION_OVERLAP;
    int end = ((int(region) + 1) * HAD_BINS_PER_REGION) - 1 + HAD_REGION_OVERLAP;
    ap_uint<5> p = normalize_had_phi(phi);

    if (start < 0) {
        return (p >= (start + HAD_PHI_BINS)) || (p <= end);
    }
    if (end >= HAD_PHI_BINS) {
        return (p >= start) || (p <= (end - HAD_PHI_BINS));
    }
    return (p >= start) && (p <= end);
}

template<bool IS_GAMMA>
static inline void stitch_candidate_pair(
    const MHHObject& hf,
    ap_uint<1> hf_boundary,
    ap_uint<3> hf_owner_region,
    ap_uint<12>& hf_et,
    ap_uint<1>& hf_cleared,
    const MHHObject& hgcal,
    ap_uint<1> hgcal_boundary,
    ap_uint<12>& hgcal_et,
    ap_uint<1>& hgcal_cleared
) {
    #pragma HLS INLINE

    bool in_local_window =
        IS_GAMMA ? gamma_phi_in_region_overlap(hgcal.Phi, hf_owner_region)
                 : had_phi_in_region_overlap(hgcal.Phi, hf_owner_region);

    if (in_local_window &&
        is_active_boundary_flag(hf_boundary, hf_et, hf_cleared) &&
        is_active_boundary_flag(hgcal_boundary, hgcal_et, hgcal_cleared) &&
        (IS_GAMMA ? match_dphi_gamma(hf.Phi, hgcal.Phi)
                  : match_dphi_had(hf.Phi, hgcal.Phi))) {
        stitch_state_pair_local(hf_et, hf_cleared, hgcal_et, hgcal_cleared);
    }
}

template<bool IS_GAMMA>
static void stitch_boundary_collection(
    MHHObject objects[MHH_N_ACTIVE_OBJECTS],
    ap_uint<10> hf_eta,
    ap_uint<10> hgcal_eta
) {
    #pragma HLS INLINE off
    #pragma HLS ARRAY_PARTITION variable=objects complete dim=0
    #pragma HLS aggregate variable=objects compact=bit

    ap_uint<12> hf_et[6];
    ap_uint<12> hgcal_et[6];
    ap_uint<1> hf_cleared[6];
    ap_uint<1> hgcal_cleared[6];
    ap_uint<1> hf_boundary[6];
    ap_uint<1> hgcal_boundary[6];

    #pragma HLS ARRAY_PARTITION variable=hf_et complete dim=0
    #pragma HLS ARRAY_PARTITION variable=hgcal_et complete dim=0
    #pragma HLS ARRAY_PARTITION variable=hf_cleared complete dim=0
    #pragma HLS ARRAY_PARTITION variable=hgcal_cleared complete dim=0
    #pragma HLS ARRAY_PARTITION variable=hf_boundary complete dim=0
    #pragma HLS ARRAY_PARTITION variable=hgcal_boundary complete dim=0

    for (int i = 0; i < 6; ++i) {
        #pragma HLS UNROLL
        hf_et[i] = objects[i].ET;
        hf_cleared[i] = 0;
        hf_boundary[i] = (objects[i].isBarrel == 0) && (objects[i].Eta == hf_eta);

        hgcal_et[i] = objects[6 + i].ET;
        hgcal_cleared[i] = 0;
        hgcal_boundary[i] = (objects[6 + i].isBarrel == 1) && (objects[6 + i].Eta == hgcal_eta);
    }

    for (int hf_idx = 0; hf_idx < 6; ++hf_idx) {
        #pragma HLS UNROLL
        ap_uint<3> owner_region =
            IS_GAMMA ? gamma_owner_region(objects[hf_idx].Phi)
                     : had_owner_region(objects[hf_idx].Phi);
        for (int hgcal_idx = 0; hgcal_idx < 6; ++hgcal_idx) {
            #pragma HLS UNROLL
            stitch_candidate_pair<IS_GAMMA>(
                objects[hf_idx], hf_boundary[hf_idx], owner_region,
                hf_et[hf_idx], hf_cleared[hf_idx],
                objects[6 + hgcal_idx], hgcal_boundary[hgcal_idx],
                hgcal_et[hgcal_idx], hgcal_cleared[hgcal_idx]);
        }
    }

    for (int i = 0; i < 6; ++i) {
        #pragma HLS UNROLL
        if (hf_cleared[i] == 1) {
            clear_gctvar(objects[i]);
        } else {
            objects[i].ET = hf_et[i];
        }

        if (hgcal_cleared[i] == 1) {
            clear_gctvar(objects[6 + i]);
        } else {
            objects[6 + i].ET = hgcal_et[i];
        }
    }
}

static void stitch_gamma_collection(MHHObject objects[MHH_N_ACTIVE_OBJECTS]) {
    #pragma HLS INLINE off
    stitch_boundary_collection<true>(objects, HF_GAMMA_BOUNDARY_ETA, HGCAL_GAMMA_BOUNDARY_ETA);
}

static void stitch_had_collection(MHHObject objects[MHH_N_ACTIVE_OBJECTS]) {
    #pragma HLS INLINE off
    stitch_boundary_collection<false>(objects, HF_HAD_BOUNDARY_ETA, HGCAL_HAD_BOUNDARY_ETA);
}

void updateParams_MHHOutput(
    MHHObject EGs[MHH_N_ACTIVE_OBJECTS],
    MHHObject EGIs[MHH_N_ACTIVE_OBJECTS],
    MHHObject Jets[MHH_N_ACTIVE_OBJECTS],
    MHHObject Taus[MHH_N_ACTIVE_OBJECTS]
) {
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=9
    #pragma HLS ARRAY_PARTITION variable=EGs complete dim=0
    #pragma HLS ARRAY_PARTITION variable=EGIs complete dim=0
    #pragma HLS ARRAY_PARTITION variable=Jets complete dim=0
    #pragma HLS ARRAY_PARTITION variable=Taus complete dim=0
    #pragma HLS aggregate variable=EGs compact=bit
    #pragma HLS aggregate variable=EGIs compact=bit
    #pragma HLS aggregate variable=Jets compact=bit
    #pragma HLS aggregate variable=Taus compact=bit

    stitch_gamma_collection(EGs);
    stitch_gamma_collection(EGIs);
    stitch_had_collection(Jets);
    stitch_had_collection(Taus);
}

void processInputLinks(
    ap_uint<576> link_in[MHH_N_INPUT_LINKS],
    MHHObject EGs[MHH_N_ACTIVE_OBJECTS],
    MHHObject EGIs[MHH_N_ACTIVE_OBJECTS],
    MHHObject Jets[MHH_N_ACTIVE_OBJECTS],
    MHHObject Taus[MHH_N_ACTIVE_OBJECTS],
    MHHSums& Sums
) {
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=9
    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=EGs complete dim=0
    #pragma HLS ARRAY_PARTITION variable=EGIs complete dim=0
    #pragma HLS ARRAY_PARTITION variable=Jets complete dim=0
    #pragma HLS ARRAY_PARTITION variable=Taus complete dim=0
    #pragma HLS aggregate variable=EGs compact=bit
    #pragma HLS aggregate variable=EGIs compact=bit
    #pragma HLS aggregate variable=Jets compact=bit
    #pragma HLS aggregate variable=Taus compact=bit

    for (int source = 0; source < 2; ++source) {
        #pragma HLS UNROLL
        bool hf_source = (source == 0);
        int link_base = hf_source ? MHH_HF_LINK_BASE : MHH_HGCAL_LINK_BASE;
        int object_base = source * 6;

        ap_uint<576> link_A = link_in[link_base + 0];
        for (int j = 0; j < 12; ++j) {
            #pragma HLS UNROLL
            ap_uint<48> raw_48b = link_A.range(48 * j + 47, 48 * j);
            MHHObject tempObj;
            tempObj.getMHHGammas(raw_48b);
            tempObj.isBarrel = hf_source ? 0 : 1;

            if (j < 6) {
                EGs[object_base + j] = tempObj;
            } else {
                EGIs[object_base + (j - 6)] = tempObj;
            }
        }

        ap_uint<576> link_B = link_in[link_base + 1];
        for (int j = 0; j < 12; ++j) {
            #pragma HLS UNROLL
            ap_uint<48> raw_48b = link_B.range(48 * j + 47, 48 * j);
            MHHObject tempObj;
            tempObj.getMHHJetsTaus(raw_48b);
            tempObj.isBarrel = hf_source ? 0 : 1;

            if (j < 6) {
                Jets[object_base + j] = tempObj;
            } else {
                Taus[object_base + (j - 6)] = tempObj;
            }
        }

    }

    ap_uint<576> hf_sums = link_in[MHH_HF_LINK_BASE + 2];
    ap_uint<576> hgcal_sums = link_in[MHH_HGCAL_LINK_BASE + 2];
    ap_uint<48> hf_ex = hf_sums.range(47, 0);
    ap_uint<48> hf_ey = hf_sums.range(95, 48);
    ap_uint<48> hgcal_ex = hgcal_sums.range(47, 0);
    ap_uint<48> hgcal_ey = hgcal_sums.range(95, 48);

    for (int scenario = 0; scenario < 3; ++scenario) {
        #pragma HLS UNROLL
        Sums.Ex[scenario] = unpack_signed_sum_component(hf_ex, scenario) +
                            unpack_signed_sum_component(hgcal_ex, scenario);
        Sums.Ey[scenario] = unpack_signed_sum_component(hf_ey, scenario) +
                            unpack_signed_sum_component(hgcal_ey, scenario);
    }

    Sums.Ht = unpack_unsigned_sum_component(hf_sums.range(143, 96)) +
              unpack_unsigned_sum_component(hgcal_sums.range(143, 96));
    Sums.SumET = unpack_unsigned_sum_component(hf_sums.range(191, 144)) +
                 unpack_unsigned_sum_component(hgcal_sums.range(191, 144));
    Sums.NObj = unpack_unsigned_sum_component(hf_sums.range(239, 192)) +
                unpack_unsigned_sum_component(hgcal_sums.range(239, 192));
}

void sortGCTVars(
    MHHObject EGsVars[MHH_N_SORT_OBJECTS],
    MHHObject EGIsVars[MHH_N_SORT_OBJECTS],
    MHHObject JetsVars[MHH_N_SORT_OBJECTS],
    MHHObject TausVars[MHH_N_SORT_OBJECTS],
    MHHObject sortedEGs[MHH_N_SORT_OBJECTS],
    MHHObject sortedEGIs[MHH_N_SORT_OBJECTS],
    MHHObject sortedJets[MHH_N_SORT_OBJECTS],
    MHHObject sortedTaus[MHH_N_SORT_OBJECTS]
) {
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=9
    #pragma HLS allocation function instances=bitonicSort32 limit=1
    #pragma HLS aggregate variable=EGsVars compact=bit
    #pragma HLS aggregate variable=EGIsVars compact=bit
    #pragma HLS aggregate variable=JetsVars compact=bit
    #pragma HLS aggregate variable=TausVars compact=bit
    #pragma HLS aggregate variable=sortedEGs compact=bit
    #pragma HLS aggregate variable=sortedEGIs compact=bit
    #pragma HLS aggregate variable=sortedJets compact=bit
    #pragma HLS aggregate variable=sortedTaus compact=bit

    bitonicSort32(EGsVars, sortedEGs);
    bitonicSort32(EGIsVars, sortedEGIs);
    bitonicSort32(JetsVars, sortedJets);
    bitonicSort32(TausVars, sortedTaus);
}

void combineObjects(
    MHHObject EGs[MHH_N_ACTIVE_OBJECTS],
    MHHObject EGIs[MHH_N_ACTIVE_OBJECTS],
    MHHObject Jets[MHH_N_ACTIVE_OBJECTS],
    MHHObject Taus[MHH_N_ACTIVE_OBJECTS],
    MHHObject EGsTop6_out[6],
    MHHObject EGIsTop6_out[6],
    MHHObject JetsTop6_out[6],
    MHHObject TausTop6_out[6]
) {
    #pragma HLS INLINE off
    #pragma HLS PIPELINE
    #pragma HLS aggregate variable=EGs compact=bit
    #pragma HLS aggregate variable=EGIs compact=bit
    #pragma HLS aggregate variable=Jets compact=bit
    #pragma HLS aggregate variable=Taus compact=bit
    #pragma HLS aggregate variable=EGsTop6_out compact=bit
    #pragma HLS aggregate variable=EGIsTop6_out compact=bit
    #pragma HLS aggregate variable=JetsTop6_out compact=bit
    #pragma HLS aggregate variable=TausTop6_out compact=bit

    MHHObject EGsVars[MHH_N_SORT_OBJECTS];
    MHHObject EGIsVars[MHH_N_SORT_OBJECTS];
    MHHObject JetsVars[MHH_N_SORT_OBJECTS];
    MHHObject TausVars[MHH_N_SORT_OBJECTS];
    MHHObject sortedEGs[MHH_N_SORT_OBJECTS];
    MHHObject sortedEGIs[MHH_N_SORT_OBJECTS];
    MHHObject sortedJets[MHH_N_SORT_OBJECTS];
    MHHObject sortedTaus[MHH_N_SORT_OBJECTS];

    #pragma HLS aggregate variable=EGsVars compact=bit
    #pragma HLS aggregate variable=EGIsVars compact=bit
    #pragma HLS aggregate variable=JetsVars compact=bit
    #pragma HLS aggregate variable=TausVars compact=bit
    #pragma HLS aggregate variable=sortedEGs compact=bit
    #pragma HLS aggregate variable=sortedEGIs compact=bit
    #pragma HLS aggregate variable=sortedJets compact=bit
    #pragma HLS aggregate variable=sortedTaus compact=bit

    MHHObject dummy;

    for (int i = 0; i < MHH_N_SORT_OBJECTS; ++i) {
        #pragma HLS UNROLL
        if (i < MHH_N_ACTIVE_OBJECTS) {
            EGsVars[i] = EGs[i];
            EGIsVars[i] = EGIs[i];
            JetsVars[i] = Jets[i];
            TausVars[i] = Taus[i];
        } else {
            EGsVars[i] = dummy;
            EGIsVars[i] = dummy;
            JetsVars[i] = dummy;
            TausVars[i] = dummy;
        }
    }

    sortGCTVars(EGsVars, EGIsVars, JetsVars, TausVars,
                sortedEGs, sortedEGIs, sortedJets, sortedTaus);

    for (int i = 0; i < 6; ++i) {
        #pragma HLS UNROLL
        EGsTop6_out[i] = sortedEGs[31 - i];
        EGIsTop6_out[i] = sortedEGIs[31 - i];
        JetsTop6_out[i] = sortedJets[31 - i];
        TausTop6_out[i] = sortedTaus[31 - i];
    }
}

void createOutput(
    MHHObject EGs[MHH_N_ACTIVE_OBJECTS],
    MHHObject EGIs[MHH_N_ACTIVE_OBJECTS],
    MHHObject Jets[MHH_N_ACTIVE_OBJECTS],
    MHHObject Taus[MHH_N_ACTIVE_OBJECTS],
    MHHObject EGsTop6_out[6],
    MHHObject EGIsTop6_out[6],
    MHHObject JetsTop6_out[6],
    MHHObject TausTop6_out[6]
) {
    #pragma HLS INLINE off
    updateParams_MHHOutput(EGs, EGIs, Jets, Taus);
    combineObjects(EGs, EGIs, Jets, Taus, EGsTop6_out, EGIsTop6_out, JetsTop6_out, TausTop6_out);
}

void processOutLinks(
    MHHObject EGsTop6[6],
    MHHObject EGIsTop6[6],
    MHHObject JetsTop6[6],
    MHHObject TausTop6[6],
    const MHHSums& Sums,
    ap_uint<576> link_out[MHH_N_OUTPUT_LINKS]
) {
    #pragma HLS INLINE off
    #pragma HLS PIPELINE II=9

    ap_uint<576> out_link0 = 0;
    ap_uint<576> out_link1 = 0;
    ap_uint<576> out_link2 = 0;

    for (int i = 0; i < 6; ++i) {
        #pragma HLS UNROLL
        out_link0.range(i * 48 + 47, i * 48) = pack_gamma_output(EGsTop6[i]);
        out_link0.range((i + 6) * 48 + 47, (i + 6) * 48) = pack_gamma_output(EGIsTop6[i]);
        out_link1.range(i * 48 + 47, i * 48) = pack_had_output(JetsTop6[i]);
        out_link1.range((i + 6) * 48 + 47, (i + 6) * 48) = pack_had_output(TausTop6[i]);
    }

    // Preserve the three 12-bit Ex/Ey hypotheses required by downstream link C.
    out_link2.range(47, 0) = pack_signed_sum_hypotheses(Sums.Ex);
    out_link2.range(95, 48) = pack_signed_sum_hypotheses(Sums.Ey);
    out_link2.range(143, 96) = pack_unsigned_sum_component(Sums.Ht);
    out_link2.range(191, 144) = pack_unsigned_sum_component(Sums.SumET);
    out_link2.range(239, 192) = pack_unsigned_sum_component(Sums.NObj);

    link_out[MHH_CARD_A_OUTPUT_BASE + 0] = out_link0;
    link_out[MHH_CARD_A_OUTPUT_BASE + 1] = out_link1;
    link_out[MHH_CARD_A_OUTPUT_BASE + 2] = out_link2;
    link_out[MHH_CARD_B_OUTPUT_BASE + 0] = 0;
    link_out[MHH_CARD_B_OUTPUT_BASE + 1] = 0;
    link_out[MHH_CARD_B_OUTPUT_BASE + 2] = 0;
}

void mhh_algo_top(
    ap_uint<576> link_in[MHH_N_INPUT_LINKS],
    ap_uint<576> link_out[MHH_N_OUTPUT_LINKS]
) {
    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
    #pragma HLS PIPELINE II=9
    #pragma HLS INTERFACE ap_ctrl_hs port=return

    MHHObject EGsTop6_out[6], EGs[MHH_N_ACTIVE_OBJECTS], EGIs[MHH_N_ACTIVE_OBJECTS];
    MHHObject Jets[MHH_N_ACTIVE_OBJECTS], Taus[MHH_N_ACTIVE_OBJECTS];
    MHHObject EGIsTop6_out[6], JetsTop6_out[6], TausTop6_out[6];
    MHHSums Sums;

    #pragma HLS aggregate variable=EGsTop6_out compact=bit
    #pragma HLS aggregate variable=EGIsTop6_out compact=bit
    #pragma HLS aggregate variable=JetsTop6_out compact=bit
    #pragma HLS aggregate variable=TausTop6_out compact=bit
    #pragma HLS aggregate variable=EGs compact=bit
    #pragma HLS aggregate variable=EGIs compact=bit
    #pragma HLS aggregate variable=Jets compact=bit
    #pragma HLS aggregate variable=Taus compact=bit
    #pragma HLS ARRAY_PARTITION variable=EGsTop6_out complete
    #pragma HLS ARRAY_PARTITION variable=EGIsTop6_out complete
    #pragma HLS ARRAY_PARTITION variable=JetsTop6_out complete
    #pragma HLS ARRAY_PARTITION variable=TausTop6_out complete

    processInputLinks(link_in, EGs, EGIs, Jets, Taus, Sums);
    createOutput(EGs, EGIs, Jets, Taus, EGsTop6_out, EGIsTop6_out, JetsTop6_out, TausTop6_out);
    processOutLinks(EGsTop6_out, EGIsTop6_out, JetsTop6_out, TausTop6_out, Sums, link_out);
}

}  // namespace p2mhh

#endif

