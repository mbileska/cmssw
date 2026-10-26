/*
 * Description:
 *   Source-faithful CMSSW adaptation of the attached GCT barrel IP2
 *   gctforgt.cpp helper algorithms.
 *
 *   Only integration mechanics were changed: C arrays are represented with
 *   std::array and symbols live in namespace gctip2.  The algorithm itself is
 *   intentionally NOT upgraded from what is present in the attachment.
 *
 *   Therefore this file deliberately preserves the following attached-firmware
 *   behavior even where the interface PDF describes a more complete target:
 *     - createJets() fills only jets[0] and jets[1];
 *     - createTaus() fills only taus[0] and taus[1];
 *     - createEGs() copies fixed source indices 0,1,8,9,15,22;
 *     - tau ratio is never assigned and therefore remains zero;
 *     - Ex and Ey are explicitly zero;
 *     - HT is the 12-bit wrapped sum of the interior ST energies;
 *     - serialized jet/tau eta and phi retain the padded seed coordinates.
 *
 *   Do not replace these points with PDF-inferred behavior without a newer
 *   authoritative firmware source.
 */

#ifndef L1Trigger_L1CaloTrigger_gctforgtIP2_cpp_h
#define L1Trigger_L1CaloTrigger_gctforgtIP2_cpp_h

#include "L1Trigger/L1CaloTrigger/interface/algo_topIP2gct_h.h"

namespace gctip2 {
namespace detail {

inline STower bestOf2(const STower& first, const STower& second) {
#pragma HLS INLINE
  // Exact strict-'>' tie behavior of the attachment.
  return (first.energy > second.energy) ? first : second;
}

inline void findMaxPhi(const std::array<STower, kPaddedPhi>& towers, STower& maximum) {
#pragma HLS INLINE
  std::array<STower, 6> step1{};
  std::array<STower, 3> step2{};
#pragma HLS ARRAY_PARTITION variable=step1 complete dim=0
#pragma HLS ARRAY_PARTITION variable=step2 complete dim=0

  for (std::size_t k = 0; k < kPaddedPhi; k += 2) {
#pragma HLS UNROLL
    step1[k / 2] = bestOf2(towers[k], towers[k + 1]);
  }
  for (std::size_t k = 0; k < 6; k += 2) {
#pragma HLS UNROLL
    step2[k / 2] = bestOf2(step1[k], step1[k + 1]);
  }
  const STower x1 = bestOf2(step2[0], step2[1]);
  maximum = bestOf2(x1, step2[2]);
}

inline void findMaxEta(const std::array<STower, kRegionEta + 2>& towers, STower& maximum) {
#pragma HLS INLINE
  std::array<STower, 4> step1{};
  std::array<STower, 2> step2{};
#pragma HLS ARRAY_PARTITION variable=step1 complete dim=0
#pragma HLS ARRAY_PARTITION variable=step2 complete dim=0

  for (std::size_t k = 0; k < kRegionEta + 2; k += 2) {
#pragma HLS UNROLL
    step1[k / 2] = bestOf2(towers[k], towers[k + 1]);
  }
  for (std::size_t k = 0; k < 4; k += 2) {
#pragma HLS UNROLL
    step2[k / 2] = bestOf2(step1[k], step1[k + 1]);
  }
  maximum = bestOf2(step2[0], step2[1]);
}

inline void getSeedPosition(const std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& towers,
                            STower& seed) {
#pragma HLS INLINE
  std::array<STower, kRegionEta + 2> etaSlices{};
#pragma HLS ARRAY_PARTITION variable=etaSlices complete dim=0

  // Attached source searches only rows 0..7 (N_GCT_STOWERS_REGION_ETA+2).
  for (std::size_t eta = 0; eta < kRegionEta + 2; ++eta) {
#pragma HLS UNROLL
    std::array<STower, kPaddedPhi> row{};
#pragma HLS ARRAY_PARTITION variable=row complete dim=0
    for (std::size_t phi = 0; phi < kPaddedPhi; ++phi) {
#pragma HLS UNROLL
      row[phi] = towers[eta][phi];
    }
    findMaxPhi(row, etaSlices[eta]);
  }
  findMaxEta(etaSlices, seed);
}

inline void getTau(std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& towers,
                   const STower& seed,
                   Tau& output) {
#pragma HLS INLINE
  std::array<std::array<TowerMask, kPaddedPhi>, kPaddedEta> mask{};
  std::array<std::array<TowerMask, kPaddedPhi>, kPaddedEta> maskDelete{};

  const unsigned int etaSeed = static_cast<unsigned int>(seed.eta);
  const unsigned int phiSeed = static_cast<unsigned int>(seed.phi);

  for (std::size_t eta = 0; eta < kPaddedEta; ++eta) {
#pragma HLS UNROLL
    for (std::size_t phi = 0; phi < kPaddedPhi; ++phi) {
#pragma HLS UNROLL
      if (eta == etaSeed && phi == phiSeed) {
        mask[eta][phi].energy = 1;
        maskDelete[eta][phi].energy = 0;
      } else {
        maskDelete[eta][phi].energy = 1;
      }
    }
  }

  ap_uint<12> value = 0;
  for (std::size_t eta = 0; eta < kRegionEta + 2; ++eta) {
#pragma HLS UNROLL
    for (std::size_t phi = 0; phi < kPaddedPhi; ++phi) {
#pragma HLS UNROLL
      const ap_uint<12> original = towers[eta][phi].energy;
      const ap_uint<12> selected = original * mask[eta][phi].energy;
      value = static_cast<ap_uint<12>>(value + selected);
      towers[eta][phi].energy = original * maskDelete[eta][phi].energy;
    }
  }

  if (seed.eta < 2 || seed.eta > 7 || seed.phi < 2 || seed.phi > 9) {
    value = 0;
  }
  output.energy = value;
  output.eta = seed.eta;
  output.phi = seed.phi;
  // ratio/spare intentionally untouched, matching the attachment.
}

inline void getJet(const std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& towers,
                   const STower& seed,
                   Jet& output) {
#pragma HLS INLINE
  constexpr std::size_t kMaskEta = kPaddedEta + 2;
  constexpr std::size_t kMaskPhi = kPaddedPhi + 2;
  std::array<std::array<TowerMask, kMaskPhi>, kMaskEta> mask{};

  const unsigned int etaSeed = static_cast<unsigned int>(seed.eta);
  const unsigned int phiSeed = static_cast<unsigned int>(seed.phi);

  for (std::size_t eta = 0; eta < kMaskEta; ++eta) {
#pragma HLS UNROLL
    if (eta >= etaSeed && eta <= etaSeed + 2) {
      for (std::size_t phi = 0; phi < kMaskPhi; ++phi) {
#pragma HLS UNROLL
        if (phi >= phiSeed && phi <= phiSeed + 2) {
          mask[eta][phi].energy = 1;
        }
      }
    }
  }

  ap_uint<12> value = 0;
  for (std::size_t eta = 0; eta < kRegionEta + 2; ++eta) {
#pragma HLS UNROLL
    for (std::size_t phi = 0; phi < kPaddedPhi; ++phi) {
#pragma HLS UNROLL
      const ap_uint<12> selected = towers[eta][phi].energy * mask[eta + 1][phi + 1].energy;
      value = static_cast<ap_uint<12>>(value + selected);
    }
  }

  if (seed.eta < 2 || seed.eta > 7 || seed.phi < 2 || seed.phi > 9) {
    value = 0;
  }
  output.energy = value;
  output.eta = seed.eta;
  output.phi = seed.phi;
}

inline void zeroRowSTowers(std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& towers,
                           const STower& seed) {
#pragma HLS INLINE
  constexpr std::size_t kMaskEta = kPaddedEta + 2;
  constexpr std::size_t kMaskPhi = kPaddedPhi + 2;
  std::array<std::array<TowerMask, kMaskPhi>, kMaskEta> mask{};
  std::array<std::array<TowerMask, kMaskPhi>, kMaskEta> inverse{};

  const unsigned int etaSeed = static_cast<unsigned int>(seed.eta);
  const unsigned int phiSeed = static_cast<unsigned int>(seed.phi);

  for (std::size_t eta = 0; eta < kMaskEta; ++eta) {
#pragma HLS UNROLL
    if (eta >= etaSeed && eta <= etaSeed + 2) {
      for (std::size_t phi = 0; phi < kMaskPhi; ++phi) {
#pragma HLS UNROLL
        if (phi >= phiSeed && phi <= phiSeed + 2) {
          mask[eta][phi].energy = 1;
        }
      }
    }
  }

  for (std::size_t eta = 0; eta < kMaskEta; ++eta) {
#pragma HLS UNROLL
    for (std::size_t phi = 0; phi < kMaskPhi; ++phi) {
#pragma HLS UNROLL
      inverse[eta][phi].energy = static_cast<ap_uint<1>>(1 - mask[eta][phi].energy);
    }
  }

  for (std::size_t eta = 0; eta < kRegionEta + 2; ++eta) {
#pragma HLS UNROLL
    for (std::size_t phi = 0; phi < kPaddedPhi; ++phi) {
#pragma HLS UNROLL
      const ap_uint<12> original = towers[eta][phi].energy;
      towers[eta][phi].energy = original * inverse[eta + 1][phi + 1].energy;
    }
  }
}

}  // namespace detail

inline void createJets(const std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& input,
                       std::array<Jet, kJetsPerRegion>& jets) {
#pragma HLS INLINE
  // Region default construction has already zeroed all six slots.  Preserve
  // the attached firmware: only slots 0 and 1 are filled.
  auto towers = input;
  STower seed;

  detail::getSeedPosition(towers, seed);
  detail::getJet(towers, seed, jets[0]);
  detail::zeroRowSTowers(towers, seed);

  detail::getSeedPosition(towers, seed);
  detail::getJet(towers, seed, jets[1]);
}

inline void createTaus(std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& towers,
                       std::array<Tau, kTausPerRegion>& taus) {
#pragma HLS INLINE
  STower seed;

  detail::getSeedPosition(towers, seed);
  detail::getTau(towers, seed, taus[0]);

  detail::getSeedPosition(towers, seed);
  detail::getTau(towers, seed, taus[1]);
}

inline void createSums(const std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& positive,
                       const std::array<std::array<STower, kPaddedPhi>, kPaddedEta>& negative,
                       Sums& output) {
#pragma HLS INLINE
  ap_uint<12> totalEnergy = 0;
  for (std::size_t eta = 0; eta < kRegionEta; ++eta) {
#pragma HLS UNROLL
    for (std::size_t phi = 0; phi < kRegionPhi; ++phi) {
#pragma HLS UNROLL
      totalEnergy = static_cast<ap_uint<12>>(totalEnergy + positive[eta + 2][phi + 2].energy);
      totalEnergy = static_cast<ap_uint<12>>(totalEnergy + negative[eta + 2][phi + 2].energy);
    }
  }
  output.HT = totalEnergy;
  output.Ex = 0;
  output.Ey = 0;
}

inline void createEGs(const std::array<GctEcalCluster, kEcalClustersPerRegion>& input,
                      std::array<GctEcalCluster, kEGsPerRegion>& egs) {
#pragma HLS INLINE
  // Literal fixed-index selection from the attached gctforgt.cpp.
  egs[0] = input[0];
  egs[1] = input[1];
  egs[2] = input[8];
  egs[3] = input[9];
  egs[4] = input[15];
  egs[5] = input[22];
}

}  // namespace gctip2

#endif
