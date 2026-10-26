#ifndef DataFormats_L1TCalorimeterPhase2_RCT_output_h
#define DataFormats_L1TCalorimeterPhase2_RCT_output_h

#include <ap_int.h>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace l1tp2 {

  class rctOutputLink {
  public:
    rctOutputLink() { linkData = (ap_uint<576>)0; }

    rctOutputLink(ap_uint<576> data) { linkData = data; }

    ap_uint<576> data() const { return linkData; }

  private:
    // Information as one value
    ap_uint<576> linkData;
    
  };

  // Concrete collection of output objects (with extra tuning information)
  typedef std::vector<l1tp2::rctOutputLink> rctOutputLinkCollection;
}  // namespace l1tp2
#endif