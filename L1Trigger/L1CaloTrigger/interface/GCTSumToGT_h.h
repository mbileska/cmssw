//------------------------------------
// Header and data formats for GCT SumCard to GT emulator
// (based heavily on algo_top.h in TO_GT_IP firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_GCTSumToGT_h
#define L1Trigger_L1CaloTrigger_GCTSumToGT_h

#include <ap_int.h>
#include <algorithm>
#include <cstdint>
#include <iostream>

#include "L1Trigger/L1CaloTrigger/interface/GCTSum_h.h"

namespace p2gctsumGT {

static constexpr int N_INPUT_LINKS_GT = 6;
static constexpr int N_OUTPUT_LINKS_GT = 6;

using p2gctsum::GCTsum;
using p2gctsum::GCTvar;
typedef ap_uint<10> loop;

class GTvar {
public:
  ap_uint<1> isValid;
  ap_uint<16> ET;
  ap_uint<13> Phi;
  ap_uint<14> Eta;
  ap_uint<4> PtClusterSeed;
  ap_uint<20> Spare;

  GTvar() {
    isValid = 0;
    ET = 0;
    Phi = 0;
    Eta = 0;
    PtClusterSeed = 0;
    Spare = 0;
  }

  void getGTvar(ap_uint<64> i) {
    this->isValid = i.range(0, 0);
    this->ET = i.range(16, 1);
    this->Phi = i.range(29, 17);
    this->Eta = i.range(43, 30);
    this->Spare = i.range(63, 44);
  }

  ap_uint<64> pack() const {
    return (ap_uint<64>)isValid | ((ap_uint<64>)ET << 1) | ((ap_uint<64>)Phi << 17) | ((ap_uint<64>)Eta << 30) |
           ((ap_uint<64>)Spare << 44);
  }

  void convertAndPack(GCTvar& in) {
    this->isValid = 1;
    this->ET = (ap_uint<16>)in.ET << 4;

    ap_uint<13> PhiL = (ap_uint<13>)in.Phi;
    if (in.Phi >= 179 && in.Phi <= 182) {
      PhiL = 178;
    } else if (in.Phi >= 183) {
      PhiL = (360 - in.Phi) | 0x1000;
    }
    this->Phi = (ap_uint<13>(PhiL & 0x1000)) | (ap_uint<13>(PhiL & 0x1FF) * 23);

    this->Eta = (ap_uint<14>(in.Eta & 0x200)) << 4 | (ap_uint<14>(in.Eta & 0x1FF) * 23);

    this->PtClusterSeed = in.PtClusterSeed;
    this->Spare = (ap_uint<20>)in.PtClusterSeed;
  }
};

class GCTtoGT {
public:
  GTvar EGspos[6];
  GTvar EGsneg[6];
  GTvar EGIspos[6];
  GTvar EGIsneg[6];
  GTvar Jetspos[6];
  GTvar Jetsneg[6];
  GTvar Tauspos[6];
  GTvar Tausneg[6];

  ap_uint<64> Sums[4];
  ap_uint<576> linkOutput[6];
  ap_uint<576> link[6];

  GCTtoGT() {
    for (int i = 0; i < 4; i++)
      Sums[i] = 0;
    for (int i = 0; i < 6; i++) {
      linkOutput[i] = 0;
      link[i] = 0;
    }
  }

  void convertObjects(GCTvar _EGspos[6],
                      GCTvar _EGsneg[6],
                      GCTvar _EGIspos[6],
                      GCTvar _EGIsneg[6],
                      GCTvar _Jetspos[6],
                      GCTvar _Jetsneg[6],
                      GCTvar _Tauspos[6],
                      GCTvar _Tausneg[6]) {
    for (int i = 0; i < 6; i++) {
      EGspos[i].convertAndPack(_EGspos[i]);
      EGsneg[i].convertAndPack(_EGsneg[i]);
      EGIspos[i].convertAndPack(_EGIspos[i]);
      EGIsneg[i].convertAndPack(_EGIsneg[i]);
      Jetspos[i].convertAndPack(_Jetspos[i]);
      Jetsneg[i].convertAndPack(_Jetsneg[i]);
      Tauspos[i].convertAndPack(_Tauspos[i]);
      Tausneg[i].convertAndPack(_Tausneg[i]);
    }
  }

  void processSums(GCTsum sumsIn[8]) {
    for (int i = 0; i < 4; i++) {
      ap_uint<12> sumEx = sumsIn[i].Ex + sumsIn[i + 4].Ex;
      ap_uint<12> sumEy = sumsIn[i].Ey + sumsIn[i + 4].Ey;
      ap_uint<12> sumHt = sumsIn[i].Ht + sumsIn[i + 4].Ht;

      this->Sums[i] =
          (ap_uint<64>)1 | ((ap_uint<64>)(sumEx << 4) << 1) | ((ap_uint<64>)(sumEy) << 17) | ((ap_uint<64>)(sumHt << 4) << 30);
    }
  }

  void getcombinedGTfromIP() {
    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 9; j++) {
        ap_uint<10> start = 64 * j;
        ap_uint<10> end = start + 63;
        ap_uint<64> data = 0;

        switch (i) {
          case 0:
            if (j < 6)
              data = EGspos[j].pack();
            else
              data = EGsneg[j - 6].pack();
            break;
          case 1:
            if (j < 3)
              data = EGsneg[j + 3].pack();
            else
              data = EGIspos[j - 3].pack();
            break;
          case 2:
            if (j < 6)
              data = EGIsneg[j].pack();
            else
              data = Jetspos[j - 6].pack();
            break;
          case 3:
            if (j < 3)
              data = Jetspos[j + 3].pack();
            else
              data = Jetsneg[j - 3].pack();
            break;
          case 4:
            if (j < 6)
              data = Tauspos[j].pack();
            else
              data = Tausneg[j - 6].pack();
            break;
          case 5:
            if (j < 3)
              data = Tausneg[j + 3].pack();
            else if (j < 7)
              data = this->Sums[j - 3];
            else
              data = 0;
            break;
        }

        linkOutput[i].range(end, start) = data;
      }
    }
  }

  void putGTtoLink() {
    for (loop i = 0; i < 6; i++) {
      link[i] = linkOutput[i];
    }
  }
};

void algo_top_GT(ap_uint<576> link_in[N_INPUT_LINKS_GT], ap_uint<576> link_out[N_OUTPUT_LINKS_GT]);

}  // namespace p2gctsumGT

#endif
