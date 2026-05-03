//------------------------------------
// Header and data formats for Phase2L1GCTSumEmulator
// (based heavily on algo_top.h in GCT Sum firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_GCTSum_h
#define L1Trigger_L1CaloTrigger_GCTSum_h

#include <ap_int.h>
#include <algorithm>
#include <cstdint>
#include <iostream>

namespace p2gctsum {

static constexpr int N_INPUT_LINKS = 12;
static constexpr int N_OUTPUT_LINKS = 3;
static constexpr int N_GCT_CONNECTED = 4;
static constexpr int N_GCT_OBJECTS_PER_SOURCE = 6;
static constexpr int N_GCT_OBJECTS = 24;
static constexpr int N_GCT_OBJECTS_SORT = 32;
static constexpr int N_GCT_SUMS = 1;

typedef ap_uint<10> loop;
typedef ap_uint<48> xvar;

class GCTvar {
public:
  ap_uint<12> ET;
  ap_uint<10> Eta;
  ap_uint<9> Phi;
  ap_uint<4> PtClusterSeed;
  ap_uint<1> isBarrel;

  GCTvar() { ET = 0; Eta = 0; Phi = 0; PtClusterSeed = 0; isBarrel = 0; }

  GCTvar(const GCTvar& rhs) {
    ET = rhs.ET;
    Eta = rhs.Eta;
    Phi = rhs.Phi;
    PtClusterSeed = rhs.PtClusterSeed;
    isBarrel = rhs.isBarrel;
  }

  GCTvar& operator=(const GCTvar& rhs) {
    this->ET = rhs.ET;
    this->Eta = rhs.Eta;
    this->Phi = rhs.Phi;
    this->PtClusterSeed = rhs.PtClusterSeed;
    this->isBarrel = rhs.isBarrel;
    return *this;
  }

  void unpack(ap_uint<48> i, bool hasSeed) {
    this->ET = i.range(11, 0);
    this->isBarrel = i.range(47, 47);

    if (hasSeed) {
      this->Eta = i.range(17, 12);
      this->Phi = i.range(26, 18);
      this->PtClusterSeed = i.range(30, 27);
    } else {
      this->Eta = i.range(18, 12);
      this->Phi = i.range(27, 19);
      this->PtClusterSeed = 0;
    }
  }

  void getGCTvarBarrelGammas(ap_uint<48> i, ap_uint<9> phiOffset = 0) {
    this->ET = i.range(11, 0);
    this->Eta = i.range(18, 12);
    ap_uint<9> raw_phi = i.range(25, 19);
    this->Phi = (raw_phi + phiOffset) & 0x1FF;
    this->isBarrel = 1;
    this->PtClusterSeed = 0;
  }

  void getGCTvarEndcapGammas(ap_uint<48> i, ap_uint<9> phiOffset = 0) {
    this->ET = i.range(11, 0);
    this->Eta = (ap_uint<10>)i.range(18, 12) + 85;
    ap_uint<9> raw_phi = i.range(27, 19);
    this->Phi = (raw_phi + phiOffset) & 0x1FF;
    this->isBarrel = 0;
    this->PtClusterSeed = 0;
  }

  void getGCTvarJetsTaus(ap_uint<48> i, ap_uint<9> phiOffset, bool barrelFlag) {
    this->ET = i.range(11, 0);
    ap_uint<6> raw_eta = i.range(17, 12);
    ap_uint<10> eta = barrelFlag ? (ap_uint<10>)raw_eta : (ap_uint<10>)(raw_eta + 6);
    this->Eta = eta;
    ap_uint<9> raw_phi = i.range(26, 18);
    this->Phi = (raw_phi + phiOffset) & 0x1FF;
    this->PtClusterSeed = i.range(30, 27);
    this->isBarrel = (barrelFlag ? 1 : 0);
  }

  ap_uint<48> pack() const {
    ap_uint<48> out = 0;

    if (isBarrel) {
      bool isGammaFormat = (Eta > 0x3F);

      if (isGammaFormat) {
        out = ((ap_uint<48>)isBarrel << 47) |
              ((ap_uint<48>)(Phi & 0x7F) << 19) |
              ((ap_uint<48>)Eta << 12) |
              (ap_uint<48>)ET;
      } else {
        out = ((ap_uint<48>)isBarrel << 47) |
              ((ap_uint<48>)PtClusterSeed << 27) |
              ((ap_uint<48>)Phi << 18) |
              ((ap_uint<48>)Eta << 12) |
              (ap_uint<48>)ET;
      }
    } else {
      if (PtClusterSeed == 0 && Eta > 0x3F) {
        out = ((ap_uint<48>)isBarrel << 47) |
              ((ap_uint<48>)Phi << 19) |
              ((ap_uint<48>)Eta << 12) |
              (ap_uint<48>)ET;
      } else {
        out = ((ap_uint<48>)isBarrel << 47) |
              ((ap_uint<48>)PtClusterSeed << 27) |
              ((ap_uint<48>)Phi << 18) |
              ((ap_uint<48>)Eta << 12) |
              (ap_uint<48>)ET;
      }
    }

    return out;
  }
};

class GCTsum {
public:
  ap_int<16> Ex;
  ap_int<16> Ey;
  ap_uint<16> Ht;
  ap_uint<16> SumET;
  ap_uint<16> NObj;

  GCTsum() { Ex = 0; Ey = 0; Ht = 0; SumET = 0; NObj = 0; }

  GCTsum(const GCTsum& rhs) {
    Ex = rhs.Ex;
    Ey = rhs.Ey;
    Ht = rhs.Ht;
    SumET = rhs.SumET;
    NObj = rhs.NObj;
  }

  GCTsum& operator=(const GCTsum& rhs) {
    this->Ex = rhs.Ex;
    this->Ey = rhs.Ey;
    this->Ht = rhs.Ht;
    this->SumET = rhs.SumET;
    this->NObj = rhs.NObj;
    return *this;
  }

  void unpack(ap_uint<576> i) {
    this->Ex = (ap_int<16>)i.range(15, 0);
    this->Ey = (ap_int<16>)i.range(63, 48);
    this->Ht = i.range(111, 96);
    this->SumET = i.range(159, 144);
    this->NObj = i.range(207, 192);
  }
};

static const ap_uint<10> BARREL_GAMMA_BOUNDARY_ETA = 84;
static const ap_uint<10> ENDCAP_GAMMA_BOUNDARY_ETA = 85;
static const ap_uint<10> BARREL_HAD_BOUNDARY_ETA = 5;
static const ap_uint<10> ENDCAP_HAD_BOUNDARY_ETA = 6;

inline bool match_dphi_gamma(ap_uint<9> ephi, ap_uint<9> bphi) {
  ap_int<10> dphi = (ap_int<10>)ephi - (ap_int<10>)bphi;
  if (dphi > 180)
    dphi -= 360;
  if (dphi < -180)
    dphi += 360;
  return (dphi >= -1 && dphi <= 1);
}

inline bool match_dphi_had(ap_uint<9> ephi, ap_uint<9> bphi) {
  ap_int<10> dphi = (ap_int<10>)ephi - (ap_int<10>)bphi;
  if (dphi > 12)
    dphi -= 24;
  if (dphi < -12)
    dphi += 24;
  return (dphi >= -1 && dphi <= 1);
}

void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);

}  // namespace p2gctsum

#endif
