//------------------------------------
// Bitonic sort declarations for GCT SumCard emulator
// (based heavily on bitonicSort32 in firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_bitonicSort32_GCT_h
#define L1Trigger_L1CaloTrigger_bitonicSort32_GCT_h

#include <ap_int.h>
#include "L1Trigger/L1CaloTrigger/interface/GCTSum_h.h"

namespace p2gctsum {

inline bool gctvar_gt(const GCTvar& a, const GCTvar& b) { return (a.ET > b.ET); }

void bitonicSort32(GCTvar in[32], GCTvar out[32]);

}  // namespace p2gctsum

#endif
