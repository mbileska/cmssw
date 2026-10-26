#ifndef L1Trigger_L1CaloTrigger_bubl_sorter_h_h
#define L1Trigger_L1CaloTrigger_bubl_sorter_h_h

#include <iostream>
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP21_h.h" // ecalcluster is defined

namespace p2rctIP21 {

static constexpr int N = 32; // Number of clusters

void bubl_sorter(ecalcluster in[N], ecalcluster out[N]);

} // namespace p2rctIP21

#endif
