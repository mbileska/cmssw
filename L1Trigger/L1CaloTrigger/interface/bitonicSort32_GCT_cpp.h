//------------------------------------
// Bitonic sort implementation for GCT SumCard emulator
// (based heavily on bitonicSort32 in firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_bitonicSort32_GCT_cpp
#define L1Trigger_L1CaloTrigger_bitonicSort32_GCT_cpp

#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_GCT_h.h"

namespace p2gctsum {

inline void compareSwap(GCTvar& a, GCTvar& b, bool dir) {
  bool swap = dir ? gctvar_gt(a, b) : gctvar_gt(b, a);
  if (swap) {
    GCTvar tmp = a;
    a = b;
    b = tmp;
  }
}

inline void bitonicSort32(GCTvar in[32], GCTvar out[32]) {
  GCTvar data[32];
  for (int i = 0; i < 32; ++i) {
    data[i] = in[i];
  }

  for (int k = 2; k <= 32; k <<= 1) {
    for (int j = (k >> 1); j > 0; j >>= 1) {
      for (int i = 0; i < 32; ++i) {
        int ixj = i ^ j;
        if (ixj > i) {
          bool up = ((i & k) == 0);
          compareSwap(data[i], data[ixj], up);
        }
      }
    }
  }

  for (int i = 0; i < 32; ++i) {
    out[i] = data[i];
  }
}

}  // namespace p2gctsum

#endif
