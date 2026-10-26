#include "bitonicSort16_h.h"

//Main CAE block (compare and exchange)

GreaterSmaller AscendDescend(PFcluster &x, PFcluster &y){

    GreaterSmaller s;

    s.greater = (x.Energy() > y.Energy()) ? x : y;
    s.smaller = (x.Energy() > y.Energy()) ? y : x;

    return s;
}

void FourinSmallFir(PFcluster &x0, PFcluster &x1, PFcluster &x2, PFcluster &x3,
					PFcluster &y0, PFcluster &y1, PFcluster &y2, PFcluster &y3){

    GreaterSmaller res;
        res = AscendDescend(x0, x2);
        y0 = res.smaller; y2 = res.greater;

        res = AscendDescend(x1, x3);
        y1 = res.smaller; y3 = res.greater;
}

void FourinGreatFir(PFcluster &x0, PFcluster &x1, PFcluster &x2, PFcluster &x3,
					PFcluster &y0, PFcluster &y1, PFcluster &y2, PFcluster &y3){

    GreaterSmaller res;
        res = AscendDescend(x0, x2);
        y0 = res.greater; y2 = res.smaller;

        res = AscendDescend(x1, x3);
        y1 = res.greater; y3 = res.smaller;
}

void EightinSmallFir(PFcluster &x0, PFcluster &x1,PFcluster &x2, PFcluster &x3, PFcluster &x4, PFcluster &x5,
                PFcluster &x6, PFcluster &x7, PFcluster &y0, PFcluster &y1, PFcluster &y2, PFcluster &y3, PFcluster &y4, PFcluster &y5, PFcluster &y6, PFcluster &y7){

        GreaterSmaller res;
        res = AscendDescend(x0, x4);
        y0 = res.smaller; y4 = res.greater;

        res = AscendDescend(x1, x5);
        y1 = res.smaller; y5 = res.greater;

        res = AscendDescend(x2, x6);
        y2 = res.smaller; y6 = res.greater;

        res = AscendDescend(x3, x7);
        y3 = res.smaller; y7 = res.greater;

}

void EightinGreatFir(PFcluster &x0, PFcluster &x1, PFcluster &x2, PFcluster &x3, PFcluster &x4, PFcluster &x5,
                PFcluster &x6, PFcluster &x7, PFcluster &y0, PFcluster &y1, PFcluster &y2, PFcluster &y3, PFcluster &y4, PFcluster &y5, PFcluster &y6,PFcluster &y7){

        GreaterSmaller res;
        res = AscendDescend(x0, x4);
        y0 = res.greater; y4 = res.smaller;

        res = AscendDescend(x1, x5);
        y1 = res.greater; y5 = res.smaller;

        res = AscendDescend(x2, x6);
        y2 = res.greater; y6 = res.smaller;

        res = AscendDescend(x3, x7);
        y3 = res.greater; y7= res.smaller;

}

void bitonicSort16AB(PFcluster in[Nbclusters], PFcluster out[Nbclusters]){





GreaterSmaller result;

PFcluster a[Nbclusters], b[Nbclusters], c[Nbclusters], d[Nbclusters], e[Nbclusters], f[Nbclusters], g[Nbclusters], h[Nbclusters], l[Nbclusters];


//Starting of first stage

for(loop i=0; i<Nbclusters/4; i++){

    result = AscendDescend(in[4*i], in[4*i+1]);
    a[4*i] = result.smaller; a[4*i+1] = result.greater; 
}

for(loop i=0; i<Nbclusters/4; i++){

    result = AscendDescend(in[4*i+2], in[4*i+3]);
    a[4*i+2] = result.greater; a[4*i+3] = result.smaller; 
}

//Starting of second stage

for(loop i=0; i<Nbclusters/8; i++){

    FourinSmallFir(a[8*i], a[8*i+1], a[8*i+2], a[8*i+3], b[8*i], b[8*i+1], b[8*i+2], b[8*i+3]);
    FourinGreatFir(a[8*i+4], a[8*i+5], a[8*i+6], a[8*i+7], b[8*i+4], b[8*i+5], b[8*i+6], b[8*i+7]); 
}

//Starting of third stage

for(loop i=0; i<Nbclusters/8; i++){

    result = AscendDescend(b[8*i], b[8*i+1]);
    c[8*i] = result.smaller; c[8*i+1] = result.greater; 
 
    result = AscendDescend(b[8*i+2], b[8*i+3]);
    c[8*i+2] = result.smaller; c[8*i+3] = result.greater;
}

for(loop i=0; i<Nbclusters/8; i++){

    result = AscendDescend(b[8*i+4], b[8*i+5]);
    c[8*i+4] = result.greater; c[8*i+5] = result.smaller; 

    result = AscendDescend(b[8*i+6], b[8*i+7]);
    c[8*i+6] = result.greater; c[8*i+7] = result.smaller;
}

//Starting of fourth stage

EightinSmallFir(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]);
EightinGreatFir(c[8], c[9], c[10], c[11], c[12], c[13], c[14], c[15], d[8], d[9], d[10], d[11], d[12], d[13], d[14], d[15]);

//Starting fifth stage

for(loop i=0; i<Nbclusters/8; i++){

    FourinSmallFir(d[4*i], d[4*i+1], d[4*i+2], d[4*i+3], e[4*i], e[4*i+1], e[4*i+2], e[4*i+3]);
}

for(loop i=0; i<Nbclusters/8; i++){

    FourinGreatFir(d[4*i+8], d[4*i+9], d[4*i+10], d[4*i+11], e[4*i+8], e[4*i+9], e[4*i+10], e[4*i+11]);
}

//Starting sixth stage

for(loop i=0; i<Nbclusters/4; i++)
{

    result = AscendDescend(e[2*i], e[2*i+1]);
    f[2*i] = result.smaller; f[2*i+1] = result.greater;
}

for(loop i=0; i<Nbclusters/4; i++){

    result = AscendDescend(e[2*i+8], e[2*i+9]);
    f[2*i+8] = result.greater; f[2*i+9] = result.smaller;
}

//Starting seventh stage

for(loop i=0; i<Nbclusters/2; i++){

    result = AscendDescend(f[i], f[i+8]);
    g[i] = result.smaller; g[i+8] = result.greater;
}

//starting eighth stage

for(loop i=0; i<Nbclusters/8; i++){

    EightinSmallFir(g[8*i], g[8*i+1], g[8*i+2], g[8*i+3], g[8*i+4], g[8*i+5], g[8*i+6], g[8*i+7],
                h[8*i], h[8*i+1], h[8*i+2], h[8*i+3], h[8*i+4], h[8*i+5], h[8*i+6], h[8*i+7]);
}

//starting ninth stage

for(loop i=0; i<Nbclusters/4; i++){

    FourinSmallFir(h[4*i], h[4*i+1], h[4*i+2], h[4*i+3], l[4*i], l[4*i+1], l[4*i+2], l[4*i+3]);
}

//starting tenth stage

for(loop i=0; i<Nbclusters/2; i++){

    result = AscendDescend(l[2*i], l[2*i+1]);
    out[2*i] = result.smaller; out[2*i+1] = result.greater;
}

//out[0] = in[0] ;

}

