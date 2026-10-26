/*
 * Description:
 *   CMSSW-safe type and packing definitions for the supplied GCT barrel IP2
 *   firmware.
 *
 *   Source basis:
 *     - attached firmware: algo_topIP2gct.h
 *     - attached firmware: algo_topIP2gct.cpp / gctforgt.cpp
 *
 *   The supplied firmware uses globally-scoped macros/classes whose names
 *   overlap with the existing IP1 firmware.  This adaptation changes only the
 *   C++ integration surface: constants and types are placed in namespace
 *   gctip2 and std::array is used at the CMSSW boundary.  Bit widths, default
 *   values, object packing, input/output counts, and algorithm semantics are
 *   intentionally preserved.
 *
 *   IMPORTANT:
 *   This file does not "complete" behavior that is absent from the attached
 *   firmware.  In particular, the attached source only fills two jet slots,
 *   two tau slots, uses fixed EG source indices, and leaves Ex/Ey at zero.
 */

#ifndef L1Trigger_L1CaloTrigger_algo_topIP2gct_h_h
#define L1Trigger_L1CaloTrigger_algo_topIP2gct_h_h

#include <ap_int.h>

#include <array>
#include <cstddef>

namespace gctip2 {

constexpr std::size_t kInputLinks = 24;
constexpr std::size_t kOutputLinks = 6;
constexpr std::size_t kRegionEta = 6;
constexpr std::size_t kRegionPhi = 8;
constexpr std::size_t kPaddedEta = kRegionEta + 4;
constexpr std::size_t kPaddedPhi = kRegionPhi + 4;
constexpr std::size_t kEcalClustersPerRegion = 72;
constexpr std::size_t kJetsPerRegion = 6;
constexpr std::size_t kTausPerRegion = 6;
constexpr std::size_t kEGsPerRegion = 6;
constexpr std::size_t kObjects48PerLink = 12;

using LinkWord = ap_uint<576>;

struct TowerMask {
  ap_uint<1> energy{0};
  ap_uint<5> eta{31};
  ap_uint<4> phi{15};
};

struct Tau {
  ap_uint<12> energy{0};
  ap_int<6> eta{0};
  ap_int<7> phi{0};
  ap_uint<4> ratio{0};
  // The attached header declares 19 spare bits.  Only bits 31..47 can be
  // serialized into the 48-bit word, exactly as in its shift-based packer.
  ap_uint<19> spare{0};

  ap_uint<48> pack() const {
    ap_uint<48> data = 0;
    data |= static_cast<ap_uint<48>>(energy);
    data |= static_cast<ap_uint<48>>(static_cast<ap_uint<6>>(eta)) << 12;
    data |= static_cast<ap_uint<48>>(static_cast<ap_uint<7>>(phi)) << 18;
    data |= static_cast<ap_uint<48>>(ratio) << 27;
    data |= static_cast<ap_uint<48>>(spare) << 31;
    return data;
  }

  void unpack(ap_uint<48> data) {
    energy = data.range(11, 0);
    eta = data.range(17, 12);
    phi = data.range(24, 18);
    ratio = data.range(30, 27);
    spare = data.range(47, 31);
  }
};

struct Sums {
  ap_uint<12> Ex{0};
  ap_uint<12> Ey{0};
  ap_uint<12> HT{0};
  ap_uint<12> spare{0};

  ap_uint<48> pack() const {
    ap_uint<48> data = 0;
    data |= static_cast<ap_uint<48>>(Ex);
    data |= static_cast<ap_uint<48>>(Ey) << 12;
    data |= static_cast<ap_uint<48>>(HT) << 24;
    data |= static_cast<ap_uint<48>>(spare) << 36;
    return data;
  }

  void unpack(ap_uint<48> data) {
    Ex = data.range(11, 0);
    Ey = data.range(23, 12);
    HT = data.range(35, 24);
    spare = data.range(47, 36);
  }
};

struct Jet {
  ap_uint<12> energy{0};
  ap_int<6> eta{0};
  ap_int<7> phi{0};
  ap_uint<4> ratio{0};
  ap_uint<19> spare{0};

  ap_uint<48> pack() const {
    ap_uint<48> data = 0;
    data |= static_cast<ap_uint<48>>(energy);
    data |= static_cast<ap_uint<48>>(static_cast<ap_uint<6>>(eta)) << 12;
    data |= static_cast<ap_uint<48>>(static_cast<ap_uint<7>>(phi)) << 18;
    data |= static_cast<ap_uint<48>>(ratio) << 27;
    data |= static_cast<ap_uint<48>>(spare) << 31;
    return data;
  }

  void unpack(ap_uint<48> data) {
    energy = data.range(11, 0);
    eta = data.range(17, 12);
    phi = data.range(24, 18);
    ratio = data.range(30, 27);
    spare = data.range(47, 31);
  }
};

struct GctEcalCluster {
  ap_uint<12> energy{0};
  ap_uint<7> eta{0};
  ap_uint<7> phi{0};
  ap_uint<22> spare{0};

  ap_uint<48> pack() const {
    ap_uint<48> data = 0;
    data |= static_cast<ap_uint<48>>(energy);
    data |= static_cast<ap_uint<48>>(eta) << 12;
    data |= static_cast<ap_uint<48>>(phi) << 19;
    data |= static_cast<ap_uint<48>>(spare) << 26;
    return data;
  }

  void unpack(ap_uint<48> data) {
    energy = data.range(11, 0);
    eta = data.range(18, 12);
    phi = data.range(25, 19);
    spare = data.range(47, 26);
  }

  // Literal equivalent of fillgctecalcluster(ap_uint<64>) in the attached
  // firmware: consume only energy/eta/phi and ignore bits above 25.
  void fillFrom64(ap_uint<64> input) {
    energy = input.range(11, 0);
    eta = input.range(18, 12);
    phi = input.range(25, 19);
  }
};

struct STower {
  ap_uint<12> energy{0};
  ap_uint<12> EMenergy{0};
  ap_uint<4> eta{15};
  ap_uint<5> phi{31};
  ap_uint<15> flags{0};

  ap_uint<48> pack() const {
    ap_uint<48> data = 0;
    data |= static_cast<ap_uint<48>>(energy);
    data |= static_cast<ap_uint<48>>(EMenergy) << 12;
    data |= static_cast<ap_uint<48>>(eta) << 24;
    data |= static_cast<ap_uint<48>>(phi) << 28;
    data |= static_cast<ap_uint<48>>(flags) << 33;
    return data;
  }

  void fill(ap_uint<48> input) {
    energy = input.range(11, 0);
    EMenergy = input.range(23, 12);
    eta = input.range(27, 24);
    phi = input.range(32, 28);
    flags = input.range(47, 33);
  }
};

struct Region {
  std::array<GctEcalCluster, kEcalClustersPerRegion> gctecalclusters{};
  std::array<std::array<STower, kPaddedPhi>, kPaddedEta> stowers{};
  std::array<Jet, kJetsPerRegion> jets{};
  std::array<Tau, kTausPerRegion> taus{};
  std::array<GctEcalCluster, kEGsPerRegion> egs{};
  Sums sum{};

  void fillEGLink(LinkWord& output) const {
    output = 0;
    for (std::size_t i = 0; i < kEGsPerRegion; ++i) {
#pragma HLS UNROLL
      output.range(static_cast<int>(48 * i + 47), static_cast<int>(48 * i)) = egs[i].pack();
    }
    // Slots 6..11 remain zero exactly as fillEGlink() in the attached source.
  }

  void fillJetTauLink(LinkWord& output) const {
    output = 0;
    for (std::size_t i = 0; i < kJetsPerRegion; ++i) {
#pragma HLS UNROLL
      output.range(static_cast<int>(48 * i + 47), static_cast<int>(48 * i)) = jets[i].pack();
    }
    for (std::size_t i = 0; i < kTausPerRegion; ++i) {
#pragma HLS UNROLL
      const std::size_t slot = i + kJetsPerRegion;
      output.range(static_cast<int>(48 * slot + 47), static_cast<int>(48 * slot)) = taus[i].pack();
    }
  }

  void fillSumsLink(LinkWord& output) const {
    output = 0;
    output.range(47, 0) = sum.pack();
  }
};

inline void createSums(const std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& positive,
                       const std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& negative,
                       Sums& output);
inline void createJets(const std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& input,
                       std::array<Jet, kJetsPerRegion>& output);
inline void createTaus(std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& towers,
                       std::array<Tau, kTausPerRegion>& output);
inline void createEGs(const std::array<GctEcalCluster, kEcalClustersPerRegion>& input,
                      std::array<GctEcalCluster, kEGsPerRegion>& output);
inline void algoTop(const std::array<LinkWord, kInputLinks>& input,
                    std::array<LinkWord, kOutputLinks>& output);

}  // namespace gctip2

#endif
