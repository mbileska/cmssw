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
static constexpr int N_GCT_SUMS = 4;

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
    this->Eta = i.range(18, 12);
    ap_uint<9> raw_phi = i.range(27, 19);
    this->Phi = (raw_phi + phiOffset) & 0x1FF;
    this->isBarrel = 0;
    this->PtClusterSeed = 0;
  }

  void getGCTvarJetsTaus(ap_uint<48> i, ap_uint<9> phiOffset, bool barrelFlag) {
    this->ET = i.range(11, 0);
    this->Eta = i.range(17, 12);
    ap_uint<9> raw_phi = i.range(26, 18);
    this->Phi = (raw_phi + phiOffset) & 0x1FF;
    this->PtClusterSeed = i.range(30, 27);
    this->isBarrel = (barrelFlag ? 1 : 0);
  }


  ap_uint<48> packGamma() const {
    ap_uint<48> out = 0;
    out.range(11, 0) = ET;
    out.range(18, 12) = Eta.range(6, 0);
    out.range(27, 19) = Phi;
    out.range(47, 47) = isBarrel;
    return out;
  }

  ap_uint<48> packHadron() const {
    ap_uint<48> out = 0;
    out.range(11, 0) = ET;
    out.range(17, 12) = Eta.range(5, 0);
    out.range(26, 18) = Phi;
    out.range(30, 27) = PtClusterSeed;
    out.range(47, 47) = isBarrel;
    return out;
  }
};

class GCTsum {
public:
  ap_uint<12> Ex;
  ap_uint<12> Ey;
  ap_uint<12> Ht;
  ap_uint<12> Spare;

  GCTsum() { Ex = 0; Ey = 0; Ht = 0; Spare = 0; }

  GCTsum(const GCTsum& rhs) {
    Ex = rhs.Ex;
    Ey = rhs.Ey;
    Ht = rhs.Ht;
    Spare = rhs.Spare;
  }

  GCTsum& operator=(const GCTsum& rhs) {
    this->Ex = rhs.Ex;
    this->Ey = rhs.Ey;
    this->Ht = rhs.Ht;
    this->Spare = rhs.Spare;
    return *this;
  }

  void getGCTsum(ap_uint<48> i) {
    this->Ex = i.range(11, 0);
    this->Ey = i.range(23, 12);
    this->Ht = i.range(35, 24);
    this->Spare = i.range(47, 36);
  }

  ap_uint<48> pack() const {
    ap_uint<48> out = 0;
    out.range(11, 0) = Ex;
    out.range(23, 12) = Ey;
    out.range(35, 24) = Ht;
    out.range(47, 36) = Spare;
    return out;
  }
};

inline ap_uint<12> saturatingAdd12(ap_uint<12> a, ap_uint<12> b) {
  ap_uint<13> sum = (ap_uint<13>)a + (ap_uint<13>)b;
  if (sum > 0xFFF)
    return 0xFFF;
  return sum.range(11, 0);
}

inline ap_uint<9> dphiAbs(ap_uint<9> a, ap_uint<9> b) {
  ap_uint<9> diff1 = (a > b) ? (a - b) : (b - a);
  ap_uint<9> diff2 = 360 - diff1;
  return (diff1 < diff2) ? diff1 : diff2;
}

inline bool match_dphi(ap_uint<9> phiA, ap_uint<9> phiB) { return (dphiAbs(phiA, phiB) <= 1); }

void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);

}  // namespace p2gctsum

#endif
