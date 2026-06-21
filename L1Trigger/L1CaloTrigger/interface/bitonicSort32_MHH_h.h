
#ifndef L1Trigger_L1CaloTrigger_bitonicSort32_MHH_h
#define L1Trigger_L1CaloTrigger_bitonicSort32_MHH_h

#include <iostream>
#include "ap_int.h"
#include "L1Trigger/L1CaloTrigger/interface/MHH_h.h"

namespace p2mhh {

static constexpr int MHH_SORT_SIZE = 32;

using namespace std;

struct MHHSortHandle {
  ap_uint<12> key;  // ET
  ap_uint<5>  idx;  // original index
};

void bitonicSort32(MHHObject in[MHH_SORT_SIZE], MHHObject out[MHH_SORT_SIZE]);

}  // namespace p2mhh

#endif
