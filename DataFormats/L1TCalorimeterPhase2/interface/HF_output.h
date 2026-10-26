#ifndef DataFormats_L1TCalorimeterPhase2_HF_output_h
#define DataFormats_L1TCalorimeterPhase2_HF_output_h

#include <ap_int.h>
#include <vector>
#include <array>
#include <cstdint>
#include <cstddef>
#include <map>
#include <string>
#include <algorithm>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace l1tp2 {

  class hfOutputLink {
  public:
    hfOutputLink() { linkDataWords.fill(0ULL); }

    explicit hfOutputLink(ap_uint<576> data) {
      for (std::size_t i = 0; i < linkDataWords.size(); ++i) {
        const int lo = static_cast<int>(64 * i);
        const int hi = lo + 63;
        linkDataWords[i] = static_cast<uint64_t>(data.range(hi, lo));
        }
    }

    ap_uint<576> data() const {
      ap_uint<576> out = 0;
      for (std::size_t i = 0; i < linkDataWords.size(); ++i) {
        const int lo = static_cast<int>(64 * i);
        const int hi = lo + 63;
        out.range(hi, lo) = ap_uint<64>(linkDataWords[i]);
      }
      return out;
    }

  private:
    std::array<uint64_t, 9> linkDataWords;
  };

  // Concrete collection of output objects (with extra tuning information)
  typedef std::vector<l1tp2::hfOutputLink> hfOutputLinkCollection;
}  // namespace l1tp2
#endif