#ifndef L1Trigger_L1CaloTrigger_Phase2L1DMXSums_h
#define L1Trigger_L1CaloTrigger_Phase2L1DMXSums_h

#include <ap_int.h>

namespace p2dmx {

  // MHH/GCTB source-sums ABI v2.  The returned 576-bit frame contains two
  // active 64-bit transport words; all reserved bits and words are zero.
  inline ap_uint<576> packSumsV2(
      int ex, int ey, int htx, int hty, unsigned int sumEt, unsigned int ht, unsigned int nObj) {
    ap_uint<576> output = 0;
    output.range(11, 0) = ap_int<12>(ex);
    output.range(23, 12) = ap_int<12>(ey);
    output.range(35, 24) = ap_int<12>(htx);
    output.range(47, 36) = ap_int<12>(hty);
    output.range(75, 64) = ap_uint<12>(sumEt);
    output.range(87, 76) = ap_uint<12>(ht);
    output.range(99, 88) = ap_uint<12>(nObj);
    return output;
  }

}  // namespace p2dmx

#endif
