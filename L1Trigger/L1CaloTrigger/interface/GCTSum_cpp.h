//------------------------------------
// More logic for Phase2L1GCTSumEmulator
// (based heavily on algo_top.cpp in GCT Sum firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_GCTSum_cpp
#define L1Trigger_L1CaloTrigger_GCTSum_cpp

#include <ap_int.h>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>

#include "L1Trigger/L1CaloTrigger/interface/GCTSum_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_GCT_h.h"

namespace p2gctsum {

inline ap_uint<48> pack_gamma_output(const GCTvar& object) {
  ap_uint<48> out = 0;
  out = ((ap_uint<48>)object.isBarrel << 47) | ((ap_uint<48>)object.Phi << 19) | ((ap_uint<48>)object.Eta << 12) |
        (ap_uint<48>)object.ET;
  return out;
}

inline ap_uint<48> pack_had_output(const GCTvar& object) {
  ap_uint<48> out = 0;
  out = ((ap_uint<48>)object.isBarrel << 47) | ((ap_uint<48>)object.PtClusterSeed << 27) |
        ((ap_uint<48>)object.Phi << 18) | ((ap_uint<48>)object.Eta << 12) | (ap_uint<48>)object.ET;
  return out;
}

inline ap_int<12> unpack_signed_sum_component(ap_uint<48> word, ap_uint<2> index) {
  switch (index) {
    case 0:
      return (ap_int<12>)word.range(11, 0);
    case 1:
      return (ap_int<12>)word.range(23, 12);
    default:
      return (ap_int<12>)word.range(35, 24);
  }
}

inline ap_uint<12> unpack_unsigned_sum_component(ap_uint<48> word) {
  return word.range(11, 0);
}

inline ap_uint<48> pack_signed_sum_component(ap_int<16> value) {
  ap_uint<48> out = 0;
  out.range(15, 0) = (ap_uint<16>)value;
  return out;
}

inline ap_uint<48> pack_unsigned_sum_component(ap_uint<16> value) {
  ap_uint<48> out = 0;
  out.range(15, 0) = value;
  return out;
}

inline void processInputLinks(ap_uint<576> link_in[N_INPUT_LINKS],
                              GCTvar EGs[N_GCT_OBJECTS],
                              GCTvar EGIs[N_GCT_OBJECTS],
                              GCTvar Jets[N_GCT_OBJECTS],
                              GCTvar Taus[N_GCT_OBJECTS],
                              GCTsum& Sums) {
  ap_uint<9> gammaOffset[4] = {0, 0, 120, 240};
  ap_uint<9> hadOffset[4] = {0, 0, 8, 16};
  ap_uint<2> sumScenario[4] = {0, 0, 1, 2};

  for (int i = 0; i < N_GCT_CONNECTED; i++) {
    int link_base = 3 * i;
    bool isBarrel = (i > 0);
    ap_uint<9> currentGammaOffset = gammaOffset[i];
    ap_uint<9> currentHadOffset = hadOffset[i];
    ap_uint<2> currentSumScenario = sumScenario[i];

    ap_uint<576> link_A = link_in[link_base + 0];
    for (int j = 0; j < 12; j++) {
      ap_uint<48> raw_48b = link_A.range(48 * j + 47, 48 * j);
      GCTvar tempObj;
      if (isBarrel) {
        tempObj.getGCTvarBarrelGammas(raw_48b, currentGammaOffset);
      } else {
        tempObj.getGCTvarEndcapGammas(raw_48b, 0);
        tempObj.Eta = (ap_uint<10>)tempObj.Eta + 85;
      }

      if (j < 6)
        EGs[i * 6 + j] = tempObj;
      else
        EGIs[i * 6 + (j - 6)] = tempObj;
    }
    ap_uint<576> link_B = link_in[link_base + 1];
    for (int j = 0; j < 12; j++) {
      ap_uint<48> raw_48b = link_B.range(48 * j + 47, 48 * j);
      GCTvar tempObj;
      tempObj.getGCTvarJetsTaus(raw_48b, currentHadOffset, isBarrel);
      if (!isBarrel)
        tempObj.Eta = (ap_uint<10>)tempObj.Eta + 6;

      if (j < 6)
        Jets[i * 6 + j] = tempObj;
      else
        Taus[i * 6 + (j - 6)] = tempObj;
    }

    ap_uint<576> link_C = link_in[link_base + 2];
    ap_uint<48> raw_ex = link_C.range(47, 0);
    ap_uint<48> raw_ey = link_C.range(95, 48);
    ap_uint<48> raw_ht = link_C.range(143, 96);
    ap_uint<48> raw_sumet = link_C.range(191, 144);
    ap_uint<48> raw_nobj = link_C.range(239, 192);

    Sums.Ex += (ap_int<16>)unpack_signed_sum_component(raw_ex, currentSumScenario);
    Sums.Ey += (ap_int<16>)unpack_signed_sum_component(raw_ey, currentSumScenario);
    Sums.Ht += (ap_uint<16>)unpack_unsigned_sum_component(raw_ht);
    Sums.SumET += (ap_uint<16>)unpack_unsigned_sum_component(raw_sumet);
    Sums.NObj += (ap_uint<16>)unpack_unsigned_sum_component(raw_nobj);
  }
}

inline bool is_boundary_barrel_gamma(const GCTvar& o) {
  return (o.ET != 0) && (o.isBarrel == 1) && (o.Eta == BARREL_GAMMA_BOUNDARY_ETA);
}

inline bool is_boundary_endcap_gamma(const GCTvar& o) {
  return (o.ET != 0) && (o.isBarrel == 0) && (o.Eta == ENDCAP_GAMMA_BOUNDARY_ETA);
}

inline bool is_boundary_barrel_had(const GCTvar& o) {
  return (o.ET != 0) && (o.isBarrel == 1) && (o.Eta == BARREL_HAD_BOUNDARY_ETA);
}

inline bool is_boundary_endcap_had(const GCTvar& o) {
  return (o.ET != 0) && (o.isBarrel == 0) && (o.Eta == ENDCAP_HAD_BOUNDARY_ETA);
}

inline void clear_gctvar(GCTvar& o) {
  o.ET = 0;
  o.Eta = 0;
  o.Phi = 0;
  o.PtClusterSeed = 0;
  o.isBarrel = 0;
}

inline void stitch_pair_keep_higher_pt(GCTvar& ec, GCTvar& br) {
  ap_uint<12> ec_et = ec.ET;
  ap_uint<12> br_et = br.ET;
  ap_uint<12> sum = saturatingAdd12(ec_et, br_et);

  if (ec_et > br_et) {
    ec.ET = sum;
    clear_gctvar(br);
  } else {
    br.ET = sum;
    clear_gctvar(ec);
  }
}

inline void stitch_had_pair(GCTvar& ec, GCTvar& br) {
  stitch_pair_keep_higher_pt(ec, br);
}

inline void updateParams_GCTOutput(GCTvar EGs[N_GCT_OBJECTS],
                                   GCTvar EGIs[N_GCT_OBJECTS],
                                   GCTvar Jets[N_GCT_OBJECTS],
                                   GCTvar Taus[N_GCT_OBJECTS]) {
  for (int b = 1; b < 4; ++b) {
    int barrel_base = b * 6;

    for (int ie = 0; ie < 6; ++ie) {
      for (int ib = 0; ib < 6; ++ib) {
        int i_ec = ie;
        int i_br = barrel_base + ib;

        if (is_boundary_endcap_gamma(EGs[i_ec]) && is_boundary_barrel_gamma(EGs[i_br]) &&
            match_dphi_gamma(EGs[i_ec].Phi, EGs[i_br].Phi)) {
          stitch_pair_keep_higher_pt(EGs[i_ec], EGs[i_br]);
        }

        if (is_boundary_endcap_gamma(EGIs[i_ec]) && is_boundary_barrel_gamma(EGIs[i_br]) &&
            match_dphi_gamma(EGIs[i_ec].Phi, EGIs[i_br].Phi)) {
          stitch_pair_keep_higher_pt(EGIs[i_ec], EGIs[i_br]);
        }

        if (is_boundary_endcap_had(Jets[i_ec]) && is_boundary_barrel_had(Jets[i_br]) &&
            match_dphi_had(Jets[i_ec].Phi, Jets[i_br].Phi)) {
          stitch_had_pair(Jets[i_ec], Jets[i_br]);
        }

        if (is_boundary_endcap_had(Taus[i_ec]) && is_boundary_barrel_had(Taus[i_br]) &&
            match_dphi_had(Taus[i_ec].Phi, Taus[i_br].Phi)) {
          stitch_had_pair(Taus[i_ec], Taus[i_br]);
        }
      }
    }
  }
}

inline void sortGCTVars(GCTvar EGsVars[32],
                        GCTvar EGIsVars[32],
                        GCTvar JetsVars[32],
                        GCTvar TausVars[32],
                        GCTvar sortedEGs[32],
                        GCTvar sortedEGIs[32],
                        GCTvar sortedJets[32],
                        GCTvar sortedTaus[32]) {
  bitonicSort32(EGsVars, sortedEGs);
  bitonicSort32(EGIsVars, sortedEGIs);
  bitonicSort32(JetsVars, sortedJets);
  bitonicSort32(TausVars, sortedTaus);
}

inline void combineObjects(GCTvar EGs[24],
                           GCTvar EGIs[24],
                           GCTvar Jets[24],
                           GCTvar Taus[24],
                           GCTvar EGsTop6[6],
                           GCTvar EGIsTop6[6],
                           GCTvar JetsTop6[6],
                           GCTvar TausTop6[6]) {
  GCTvar EGsVars[32];
  GCTvar EGIsVars[32];
  GCTvar JetsVars[32];
  GCTvar TausVars[32];

  GCTvar sortedEGs[32];
  GCTvar sortedEGIs[32];
  GCTvar sortedJets[32];
  GCTvar sortedTaus[32];

  GCTvar dummy;

  for (int i = 0; i < 32; i++) {
    if (i < 24) {
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

  sortGCTVars(EGsVars, EGIsVars, JetsVars, TausVars, sortedEGs, sortedEGIs, sortedJets, sortedTaus);

  for (int i = 0; i < 6; i++) {
    EGsTop6[i] = sortedEGs[31 - i];
    EGIsTop6[i] = sortedEGIs[31 - i];
    JetsTop6[i] = sortedJets[31 - i];
    TausTop6[i] = sortedTaus[31 - i];
  }
}

inline void getHighestInPt(GCTvar EGsTop6[6],
                           GCTvar EGIsTop6[6],
                           GCTvar JetsTop6[6],
                           GCTvar TausTop6[6],
                           GCTvar EGsTop6_reg[6],
                           GCTvar EGIsTop6_reg[6],
                           GCTvar JetsTop6_reg[6],
                           GCTvar TausTop6_reg[6]) {
  for (int i = 0; i < 6; ++i) {
    EGsTop6_reg[i] = EGsTop6[i];
    EGIsTop6_reg[i] = EGIsTop6[i];
    JetsTop6_reg[i] = JetsTop6[i];
    TausTop6_reg[i] = TausTop6[i];
  }
}

inline void processOutLinks(GCTvar EGsTop6[6],
                            GCTvar EGIsTop6[6],
                            GCTvar JetsTop6[6],
                            GCTvar TausTop6[6],
                            const GCTsum& Sums,
                            ap_uint<576> link_out[N_OUTPUT_LINKS]) {
  ap_uint<576> out_link0 = 0;
  ap_uint<576> out_link1 = 0;
  ap_uint<576> out_link2 = 0;

  for (int i = 0; i < 6; i++) {
    out_link0.range(i * 48 + 47, i * 48) = pack_gamma_output(EGsTop6[i]);
  }

  for (int i = 0; i < 6; i++) {
    int slot = i + 6;
    out_link0.range(slot * 48 + 47, slot * 48) = pack_gamma_output(EGIsTop6[i]);
  }

  for (int i = 0; i < 6; i++) {
    out_link1.range(i * 48 + 47, i * 48) = pack_had_output(JetsTop6[i]);
  }

  for (int i = 0; i < 6; i++) {
    int slot = i + 6;
    out_link1.range(slot * 48 + 47, slot * 48) = pack_had_output(TausTop6[i]);
  }

  out_link2.range(47, 0) = pack_signed_sum_component(Sums.Ex);
  out_link2.range(95, 48) = pack_signed_sum_component(Sums.Ey);
  out_link2.range(143, 96) = pack_unsigned_sum_component(Sums.Ht);
  out_link2.range(191, 144) = pack_unsigned_sum_component(Sums.SumET);
  out_link2.range(239, 192) = pack_unsigned_sum_component(Sums.NObj);

  link_out[0] = out_link0;
  link_out[1] = out_link1;
  link_out[2] = out_link2;
}

inline void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]) {
  GCTvar EGs[24], EGIs[24], Jets[24], Taus[24];
  GCTsum Sums;

  GCTvar EGsTop6[6], EGIsTop6[6], JetsTop6[6], TausTop6[6];
  GCTvar EGsTop6_reg[6], EGIsTop6_reg[6], JetsTop6_reg[6], TausTop6_reg[6];

  processInputLinks(link_in, EGs, EGIs, Jets, Taus, Sums);
  updateParams_GCTOutput(EGs, EGIs, Jets, Taus);
  combineObjects(EGs, EGIs, Jets, Taus, EGsTop6, EGIsTop6, JetsTop6, TausTop6);
  getHighestInPt(EGsTop6, EGIsTop6, JetsTop6, TausTop6, EGsTop6_reg, EGIsTop6_reg, JetsTop6_reg, TausTop6_reg);
  processOutLinks(EGsTop6_reg, EGIsTop6_reg, JetsTop6_reg, TausTop6_reg, Sums, link_out);
}

}  // namespace p2gctsum

#endif
