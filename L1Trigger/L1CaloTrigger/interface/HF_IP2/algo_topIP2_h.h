#ifndef ALGO_TOPIP2_H
#define ALGO_TOPIP2_H

#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

#include "common/layer1_objs.h"
#include "common/puppi.h"
#include "common/pf.h"
#include "common/layer1_multiplicities.h"

#define N_HF_REGIONS  6
#define NCALO         8
#define LINK_WIDTH    576
#define OBJECT_WIDTH  64

using namespace std;
typedef ap_uint<10> loop;

void algo_topIP2(
    ap_uint<576> link_in[N_HF_REGIONS],
    ap_uint<576> link_out[N_HF_REGIONS]
);

// Direct array interface - for testing and validation
void puppi_top_arrays(
    l1ct::PFRegion   regions[N_HF_REGIONS],
    l1ct::HadCaloObj caloin[N_HF_REGIONS][NCALO],
    l1ct::PuppiObj   puppiOut[N_HF_REGIONS][NNEUTRALS]
);

// Single region processing
void puppi_sector(
    const l1ct::PFRegion   &region,
    const l1ct::HadCaloObj  caloin[NCALO],
    l1ct::PuppiObj          puppiOut[NNEUTRALS]
);

//to create region configuration required for fwdlinpuppi
l1ct::PFRegion createPFRegion(int region);

#endif
