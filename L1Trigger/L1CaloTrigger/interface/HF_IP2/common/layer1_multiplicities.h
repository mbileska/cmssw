#ifndef PUPPI_LAYER1_MULTIPLICITIES_H
#define PUPPI_LAYER1_MULTIPLICITIES_H

#define REG_HF

#define NCALO 8

#define NNEUTRALS 8

#define NALLNEUTRALS NCALO

#define N_SECTORS_PF 24

// Dummy values for barrel compatibility (not used in HF)
#define NMU 1
#define NPV 1
#define NTRACK 1
#define NEMCALO 1
#define NPHOTON NEMCALO
#define NSELCALO 1
#define NTRACK_PER_EMCALO_EGCOMP 1
#define NPUPPIFINALSORTED NNEUTRALS

namespace l1ct {

    template <int N> struct ct_log2_ceil {
      enum { value = ct_log2_ceil<(N >> 1) + (N % 2)>::value + 1 };
    };
    template <> struct ct_log2_ceil<2> {
      enum { value = 1 };
    };
    template <> struct ct_log2_ceil<1> {
      enum { value = 0 };
    };

}

#endif
