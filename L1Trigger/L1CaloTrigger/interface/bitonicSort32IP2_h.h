#ifndef L1Trigger_L1CaloTrigger_bitonicSort32IP2_h_h
#define L1Trigger_L1CaloTrigger_bitonicSort32IP2_h_h

#include <ap_int.h>

#include <cstddef>

#include "L1Trigger/L1CaloTrigger/interface/algo_topIP2gct_h.h"

namespace gctip2 {

constexpr std::size_t kBitonicSortSize = 32;
using BitonicLoop = ap_uint<6>;

struct BitonicGreaterSmaller {
  GctEcalCluster greater;
  GctEcalCluster smaller;
};

void bitonicSort32(GctEcalCluster input[kBitonicSortSize], GctEcalCluster output[kBitonicSortSize]);

}  // namespace gctip2

#endif
