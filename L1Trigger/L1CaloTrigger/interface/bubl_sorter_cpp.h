#ifndef L1Trigger_L1CaloTrigger_bubl_sorter_cpp_h
#define L1Trigger_L1CaloTrigger_bubl_sorter_cpp_h

#include "L1Trigger/L1CaloTrigger/interface/bubl_sorter_h.h"
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP21_h.h"

namespace p2rctIP21 {

inline void bubl_sorter(ecalcluster in[N], ecalcluster out[N]) {

    ecalcluster temp; // Temporary variable for swapping entire objects

    for (int i = 0; i < N - 1; ++i) {
        for (int j = 0; j < N - i - 1; ++j) {
            if (in[j].energy > in[j + 1].energy) {
                // Swap the entire ecalcluster objects
                temp = in[j];
                in[j] = in[j + 1];
                in[j + 1] = temp;
            }
        }
    }
    for (int i = 0; i < N; ++i) {
        out[i] = in[i];
    }
}

} // namespace p2rctIP21

#endif