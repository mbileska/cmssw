#ifndef _bubl_sorter_H_
#define _bubl_sorter_H_
#include <iostream>
#include "ap_int.h"
#include "algo_topIP1_h.h"


template <typename T, int N>
void bubl_sorter(T (&in)[N], T (&out)[N]) {
    T temp;

    for (loop i = 0; i < N - 1; ++i) {
//    #pragma HLS UNROLL
        for (loop j = 0; j < N - i - 1; ++j) {
            if (in[j].energy < in[j + 1].energy) {
                temp = in[j];
                in[j] = in[j + 1];
                in[j + 1] = temp;
            }
        }
    }

    for (loop i = 0; i < N; ++i) {
        out[i] = in[i];
    }
}

#endif
