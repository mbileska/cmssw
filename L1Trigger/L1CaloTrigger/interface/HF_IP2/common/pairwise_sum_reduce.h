#ifndef PUPPI_PAIRWISE_SUM_REDUCE_H
#define PUPPI_PAIRWISE_SUM_REDUCE_H

#include "hls_reg.h"

template<typename T, unsigned int N>
struct pairwise_sum {
    static T reduce(const T in[N]) {
    
//        static constexpr int halfWidth = N / 2;
        static constexpr int halfWidth = N >>1;
        static constexpr int reducedSize = halfWidth + N % 2;
        T reduced[reducedSize];
        for (int i = 0; i < halfWidth; ++i) {
          reduced[i] = in[2*i] + in[2*i+1];
        }
        if (halfWidth != reducedSize) {
          reduced[reducedSize - 1] = in[N - 1];
        }
        return pairwise_sum<T,reducedSize>::reduce(reduced);
    }
};

template<typename T>
struct pairwise_sum<T,1> {
    static T reduce(const T in[1]) {
        return in[0];
    }
};

template<typename T>
struct pairwise_sum<T,2> {
    static T reduce(const T in[2]) {
        return in[0]+in[1];
    }
};

template<typename T, unsigned int N>
struct pipelined_ternary_sum {
    static T reduce(const T in[N]) {
        static constexpr int N3 = N / 3;
        static constexpr int reducedSize = (N + 2)/ 3;
        T reduced[reducedSize];
        for (int i = 0; i < N3; ++i) {
          reduced[i] = in[3*i] + in[3*i+1] + in[3*i+2];
        }
        if (N - 3*N3 == 2) {
          reduced[reducedSize - 1] = in[N - 1] + in[N - 2];
        } else if (N - 3*N3 == 1) {
          reduced[reducedSize - 1] = in[N - 1];
        }
        return HLS_REG<T>(pipelined_ternary_sum<T,reducedSize>::reduce(reduced));
    }
};

template<typename T>
struct pipelined_ternary_sum<T,1> {
    static T reduce(const T in[1]) {
        return in[0];
    }
};

template<typename T>
struct pipelined_ternary_sum<T,2> {
    static T reduce(const T in[2]) {
        return HLS_REG<T>(in[0]+in[1]);
    }
};

template<typename T>
struct pipelined_ternary_sum<T,3> {
    static T reduce(const T in[3]) {
        return HLS_REG<T>(in[0]+in[1]+in[2]);
    }
};

template<typename T>
struct pipelined_ternary_sum<T,4> {
    static T reduce(const T in[4]) {
        return HLS_REG<T>(in[0]+in[1]+in[2]+in[3]);
    }
};


#endif
