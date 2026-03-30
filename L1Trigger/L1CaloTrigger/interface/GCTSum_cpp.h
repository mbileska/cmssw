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

inline void processInputLinks(ap_uint<576> link_in[N_INPUT_LINKS],
                              GCTvar EGs[N_GCT_OBJECTS],
                              GCTvar EGIs[N_GCT_OBJECTS],
                              GCTvar Jets[N_GCT_OBJECTS],
                              GCTvar Taus[N_GCT_OBJECTS],
                              GCTsum Sums[N_GCT_SUMS]) {
  ap_uint<9> offset[4] = {0, 0, 120, 240};

  for (int i = 0; i < N_GCT_CONNECTED; i++) {
    int link_base = 3 * i;
    bool isBarrel = (i > 0);
    ap_uint<9> current_offset = offset[i];

    ap_uint<576> link_A = link_in[link_base + 0];
    for (int j = 0; j < 12; j++) {
      ap_uint<48> raw_48b = link_A.range(48 * j + 47, 48 * j);
      GCTvar tempObj;
      if (isBarrel)
        tempObj.getGCTvarBarrelGammas(raw_48b, current_offset);
      else
        tempObj.getGCTvarEndcapGammas(raw_48b, current_offset);

      if (j < 6)
        EGs[i * 6 + j] = tempObj;
      else
        EGIs[i * 6 + (j - 6)] = tempObj;
    }
    ap_uint<576> link_B = link_in[link_base + 1];
    for (int j = 0; j < 12; j++) {
      ap_uint<48> raw_48b = link_B.range(48 * j + 47, 48 * j);
      GCTvar tempObj;
      tempObj.getGCTvarJetsTaus(raw_48b, current_offset, isBarrel);

      if (j < 6)
        Jets[i * 6 + j] = tempObj;
      else
        Taus[i * 6 + (j - 6)] = tempObj;
    }

    ap_uint<576> link_C = link_in[link_base + 2];
    ap_uint<48> raw_sum = link_C.range(47, 0);
    Sums[i].getGCTsum(raw_sum);
  }
 }

inline void stitchPair(GCTvar& endcapObj, GCTvar& barrelObj) {
  if (endcapObj.ET == 0 || barrelObj.ET == 0)
    return;
  if (!match_dphi(endcapObj.Phi, barrelObj.Phi))
    return;

  ap_uint<12> stitchedET = saturatingAdd12(endcapObj.ET, barrelObj.ET);

   if (endcapObj.ET > barrelObj.ET) {
     endcapObj.ET = stitchedET;
   } else {
     barrelObj.ET = stitchedET;
   }
}

inline void updateParams_GCTOutput(GCTvar EGs[N_GCT_OBJECTS],
                                   GCTvar EGIs[N_GCT_OBJECTS],
                                   GCTvar Jets[N_GCT_OBJECTS],
                                   GCTvar Taus[N_GCT_OBJECTS]) {
  for (int i = 0; i < 6; i++) {
    stitchPair(EGs[i], EGs[6 + i]);
    stitchPair(EGs[i], EGs[12 + i]);
    stitchPair(EGs[i], EGs[18 + i]);

    stitchPair(EGIs[i], EGIs[6 + i]);
    stitchPair(EGIs[i], EGIs[12 + i]);
    stitchPair(EGIs[i], EGIs[18 + i]);

    stitchPair(Jets[i], Jets[6 + i]);
    stitchPair(Jets[i], Jets[12 + i]);
    stitchPair(Jets[i], Jets[18 + i]);

    stitchPair(Taus[i], Taus[6 + i]);
    stitchPair(Taus[i], Taus[12 + i]);
    stitchPair(Taus[i], Taus[18 + i]);
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
                            GCTsum Sums[4],
                            ap_uint<576> link_out[N_OUTPUT_LINKS]) {
  ap_uint<576> out_link0 = 0;
  ap_uint<576> out_link1 = 0;
  ap_uint<576> out_link2 = 0;

  for (int i = 0; i < 6; i++) {
    out_link0.range(i * 48 + 47, i * 48) = EGsTop6[i].packGamma();
  }

  for (int i = 0; i < 6; i++) {
    int slot = i + 6;
    out_link0.range(slot * 48 + 47, slot * 48) = EGIsTop6[i].packGamma();
  }

  for (int i = 0; i < 6; i++) {
    out_link1.range(i * 48 + 47, i * 48) = JetsTop6[i].packHadron();
  }

  for (int i = 0; i < 6; i++) {
    int slot = i + 6;
    out_link1.range(slot * 48 + 47, slot * 48) = TausTop6[i].packHadron();
  }

  for (int i = 0; i < 4; i++) {
    out_link2.range(i * 48 + 47, i * 48) = Sums[i].pack();
  }

  link_out[0] = out_link0;
  link_out[1] = out_link1;
  link_out[2] = out_link2;
}

inline void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]) {
  GCTvar EGs[24], EGIs[24], Jets[24], Taus[24];
  GCTsum Sums[4];

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
