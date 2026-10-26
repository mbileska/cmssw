#ifndef _BITONICSORT16_H_
#define _BITONICSORT16_H_

#include <iostream>
#include "ap_int.h"
#include "algo_topIP1gct_h.h"
//#include "algo_top_parameters.h"
#define Nbclusters 16

using namespace std;

class GreaterSmaller{
public:
    PFcluster greater, smaller;
};

void bitonicSort16AB(PFcluster in[Nbclusters], PFcluster out[Nbclusters]);

#endif
