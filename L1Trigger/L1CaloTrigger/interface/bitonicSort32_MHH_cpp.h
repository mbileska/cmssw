#ifndef L1Trigger_L1CaloTrigger_bitonicSort32_MHH_cpp
#define L1Trigger_L1CaloTrigger_bitonicSort32_MHH_cpp

#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_MHH_h.h"

namespace p2mhh {
#include <ap_int.h>



 bool AscendDescendKey( MHHSortHandle &x, MHHSortHandle &y) {
#pragma HLS INLINE off
  return x.key>y.key ;
}
void FourinSmallFirH(MHHSortHandle &x0, MHHSortHandle &x1, MHHSortHandle &x2, MHHSortHandle &x3,
                            MHHSortHandle &y0, MHHSortHandle &y1, MHHSortHandle &y2, MHHSortHandle &y3) {
#pragma HLS INLINE off
  const bool g02 = AscendDescendKey(x0, x2);
  y0 = g02 ? x2 : x0;  y2 = g02 ? x0 : x2;
  const bool g13 = AscendDescendKey(x1, x3);
  y1 = g13 ? x3 : x1;  y3 = g13 ? x1 : x3;
}
void FourinGreatFirH(MHHSortHandle &x0, MHHSortHandle &x1, MHHSortHandle &x2, MHHSortHandle &x3,
                            MHHSortHandle &y0, MHHSortHandle &y1, MHHSortHandle &y2, MHHSortHandle &y3) {
#pragma HLS INLINE off
  const bool g02 = AscendDescendKey(x0, x2);
  y0 = g02 ? x0 : x2;  y2 = g02 ? x2 : x0;
  const bool g13 = AscendDescendKey(x1, x3);
  y1 = g13 ? x1 : x3;  y3 = g13 ? x3 : x1;
}
void EightinSmallFirH(MHHSortHandle &x0, MHHSortHandle &x1, MHHSortHandle &x2, MHHSortHandle &x3,
                             MHHSortHandle &x4, MHHSortHandle &x5, MHHSortHandle &x6, MHHSortHandle &x7,
                             MHHSortHandle &y0, MHHSortHandle &y1, MHHSortHandle &y2, MHHSortHandle &y3,
                             MHHSortHandle &y4, MHHSortHandle &y5, MHHSortHandle &y6, MHHSortHandle &y7) {

  const bool g04 = AscendDescendKey(x0, x4); y0 = g04 ? x4 : x0; y4 = g04 ? x0 : x4;
  const bool g15 = AscendDescendKey(x1, x5); y1 = g15 ? x5 : x1; y5 = g15 ? x1 : x5;
  const bool g26 = AscendDescendKey(x2, x6); y2 = g26 ? x6 : x2; y6 = g26 ? x2 : x6;
  const bool g37 = AscendDescendKey(x3, x7); y3 = g37 ? x7 : x3; y7 = g37 ? x3 : x7;
}
 void EightinGreatFirH(MHHSortHandle &x0, MHHSortHandle &x1, MHHSortHandle &x2, MHHSortHandle &x3,
                             MHHSortHandle &x4, MHHSortHandle &x5, MHHSortHandle &x6, MHHSortHandle &x7,
                             MHHSortHandle &y0, MHHSortHandle &y1, MHHSortHandle &y2, MHHSortHandle &y3,
                             MHHSortHandle &y4, MHHSortHandle &y5, MHHSortHandle &y6, MHHSortHandle &y7) {

  const bool g04 = AscendDescendKey(x0, x4); y0 = g04 ? x0 : x4; y4 = g04 ? x4 : x0;
  const bool g15 = AscendDescendKey(x1, x5); y1 = g15 ? x1 : x5; y5 = g15 ? x5 : x1;
  const bool g26 = AscendDescendKey(x2, x6); y2 = g26 ? x2 : x6; y6 = g26 ? x6 : x2;
  const bool g37 = AscendDescendKey(x3, x7); y3 = g37 ? x3 : x7; y7 = g37 ? x7 : x3;
}
 void SixteenSmallFirH(MHHSortHandle &x0, MHHSortHandle &x1, MHHSortHandle &x2, MHHSortHandle &x3,
                             MHHSortHandle &x4, MHHSortHandle &x5, MHHSortHandle &x6, MHHSortHandle &x7,
                             MHHSortHandle &x8, MHHSortHandle &x9, MHHSortHandle &x10, MHHSortHandle &x11,
                             MHHSortHandle &x12, MHHSortHandle &x13, MHHSortHandle &x14, MHHSortHandle &x15,
                             MHHSortHandle &y0, MHHSortHandle &y1, MHHSortHandle &y2, MHHSortHandle &y3,
                             MHHSortHandle &y4, MHHSortHandle &y5, MHHSortHandle &y6, MHHSortHandle &y7,
                             MHHSortHandle &y8, MHHSortHandle &y9, MHHSortHandle &y10, MHHSortHandle &y11,
                             MHHSortHandle &y12, MHHSortHandle &y13, MHHSortHandle &y14, MHHSortHandle &y15) {

  const bool g0 = AscendDescendKey(x0,x8);   y0 = g0 ? x8  : x0;  y8  = g0 ? x0  : x8;
  const bool g1 = AscendDescendKey(x1,x9);   y1 = g1 ? x9  : x1;  y9  = g1 ? x1  : x9;
  const bool g2 = AscendDescendKey(x2,x10);  y2 = g2 ? x10 : x2;  y10 = g2 ? x2  : x10;
  const bool g3 = AscendDescendKey(x3,x11);  y3 = g3 ? x11 : x3;  y11 = g3 ? x3  : x11;
  const bool g4 = AscendDescendKey(x4,x12);  y4 = g4 ? x12 : x4;  y12 = g4 ? x4  : x12;
  const bool g5 = AscendDescendKey(x5,x13);  y5 = g5 ? x13 : x5;  y13 = g5 ? x5  : x13;
  const bool g6 = AscendDescendKey(x6,x14);  y6 = g6 ? x14 : x6;  y14 = g6 ? x6  : x14;
  const bool g7 = AscendDescendKey(x7,x15);  y7 = g7 ? x15 : x7;  y15 = g7 ? x7  : x15;
}
void SixteenGreatFirH(MHHSortHandle &x0, MHHSortHandle &x1, MHHSortHandle &x2, MHHSortHandle &x3,
                             MHHSortHandle &x4, MHHSortHandle &x5, MHHSortHandle &x6, MHHSortHandle &x7,
                             MHHSortHandle &x8, MHHSortHandle &x9, MHHSortHandle &x10, MHHSortHandle &x11,
                             MHHSortHandle &x12, MHHSortHandle &x13, MHHSortHandle &x14, MHHSortHandle &x15,
                             MHHSortHandle &y0, MHHSortHandle &y1, MHHSortHandle &y2, MHHSortHandle &y3,
                             MHHSortHandle &y4, MHHSortHandle &y5, MHHSortHandle &y6, MHHSortHandle &y7,
                             MHHSortHandle &y8, MHHSortHandle &y9, MHHSortHandle &y10, MHHSortHandle &y11,
                             MHHSortHandle &y12, MHHSortHandle &y13, MHHSortHandle &y14, MHHSortHandle &y15) {

  const bool g0 = AscendDescendKey(x0,x8);   y0 = g0 ? x0  : x8;  y8  = g0 ? x8  : x0;
  const bool g1 = AscendDescendKey(x1,x9);   y1 = g1 ? x1  : x9;  y9  = g1 ? x9  : x1;
  const bool g2 = AscendDescendKey(x2,x10);  y2 = g2 ? x2  : x10; y10 = g2 ? x10 : x2;
  const bool g3 = AscendDescendKey(x3,x11);  y3 = g3 ? x3  : x11; y11 = g3 ? x11 : x3;
  const bool g4 = AscendDescendKey(x4,x12);  y4 = g4 ? x4  : x12; y12 = g4 ? x12 : x4;
  const bool g5 = AscendDescendKey(x5,x13);  y5 = g5 ? x5  : x13; y13 = g5 ? x13 : x5;
  const bool g6 = AscendDescendKey(x6,x14);  y6 = g6 ? x6  : x14; y14 = g6 ? x14 : x6;
  const bool g7 = AscendDescendKey(x7,x15);  y7 = g7 ? x7  : x15; y15 = g7 ? x15 : x7;
}

// ---- top-level: sort handles, then gather payload once ----
void bitonicSort32(MHHObject in[MHH_SORT_SIZE], MHHObject out[MHH_SORT_SIZE]) {
//#pragma HLS PIPELINE
#pragma HLS INLINE off


#pragma HLS ARRAY_PARTITION variable=in
#pragma HLS ARRAY_PARTITION variable=out
#pragma HLS aggregate variable=in compact=bit
#pragma HLS aggregate variable=out compact=bit

  // MHHSortHandle stage buffers
  MHHSortHandle a[MHH_SORT_SIZE], b[MHH_SORT_SIZE], c[MHH_SORT_SIZE], d[MHH_SORT_SIZE], e[MHH_SORT_SIZE], f[MHH_SORT_SIZE], g[MHH_SORT_SIZE], h[MHH_SORT_SIZE], l[MHH_SORT_SIZE], m[MHH_SORT_SIZE], n[MHH_SORT_SIZE], o[MHH_SORT_SIZE], p[MHH_SORT_SIZE], q[MHH_SORT_SIZE], s[MHH_SORT_SIZE], r[MHH_SORT_SIZE];
	#pragma HLS ARRAY_PARTITION variable=a
	#pragma HLS ARRAY_PARTITION variable=b
	#pragma HLS ARRAY_PARTITION variable=c
	#pragma HLS ARRAY_PARTITION variable=d
	#pragma HLS ARRAY_PARTITION variable=e
	#pragma HLS ARRAY_PARTITION variable=f
	#pragma HLS ARRAY_PARTITION variable=g
	#pragma HLS ARRAY_PARTITION variable=h
	#pragma HLS ARRAY_PARTITION variable=l
	#pragma HLS ARRAY_PARTITION variable=m
	#pragma HLS ARRAY_PARTITION variable=n
	#pragma HLS ARRAY_PARTITION variable=o
	#pragma HLS ARRAY_PARTITION variable=p
	#pragma HLS ARRAY_PARTITION variable=q
	#pragma HLS ARRAY_PARTITION variable=s
	#pragma HLS ARRAY_PARTITION variable=r

  // Pack
  for (int i = 0; i < MHH_SORT_SIZE; ++i) {
#pragma HLS UNROLL
    a[i].key = in[i].ET;
    a[i].idx = i;
  }


  // ---------------- Stage 1 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/4; ++i) {
#pragma HLS UNROLL factor=2
    bool g01 = AscendDescendKey(a[4*i], a[4*i+1]);
    b[4*i]   = g01 ? a[4*i+1] : a[4*i];
    b[4*i+1] = g01 ? a[4*i]   : a[4*i+1];
  }
  for (int i = 0; i < MHH_SORT_SIZE/4; ++i) {
#pragma HLS UNROLL factor=2
    bool g23 = AscendDescendKey(a[4*i+2], a[4*i+3]);
    b[4*i+2] = g23 ? a[4*i+2] : a[4*i+3];
    b[4*i+3] = g23 ? a[4*i+3] : a[4*i+2];
  }

  // ---------------- Stage 2 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    FourinSmallFirH(b[8*i+0], b[8*i+1], b[8*i+2], b[8*i+3],
                    c[8*i+0], c[8*i+1], c[8*i+2], c[8*i+3]);
    FourinGreatFirH(b[8*i+4], b[8*i+5], b[8*i+6], b[8*i+7],
                    c[8*i+4], c[8*i+5], c[8*i+6], c[8*i+7]);
  }

  // ---------------- Stage 3 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    bool g01 = AscendDescendKey(c[8*i+0], c[8*i+1]);
    d[8*i+0] = g01 ? c[8*i+1] : c[8*i+0];
    d[8*i+1] = g01 ? c[8*i+0] : c[8*i+1];

     bool g23 = AscendDescendKey(c[8*i+2], c[8*i+3]);
    d[8*i+2] = g23 ? c[8*i+3] : c[8*i+2];
    d[8*i+3] = g23 ? c[8*i+2] : c[8*i+3];
  }
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
     bool g45 = AscendDescendKey(c[8*i+4], c[8*i+5]);
    d[8*i+4] = g45 ? c[8*i+4] : c[8*i+5];
    d[8*i+5] = g45 ? c[8*i+5] : c[8*i+4];

     bool g67 = AscendDescendKey(c[8*i+6], c[8*i+7]);
    d[8*i+6] = g67 ? c[8*i+6] : c[8*i+7];
    d[8*i+7] = g67 ? c[8*i+7] : c[8*i+6];
  }

  // ---------------- Stage 4 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/16; ++i) {
#pragma HLS UNROLL factor=2
    EightinSmallFirH(d[16*i+0], d[16*i+1], d[16*i+2], d[16*i+3],
                     d[16*i+4], d[16*i+5], d[16*i+6], d[16*i+7],
                     e[16*i+0], e[16*i+1], e[16*i+2], e[16*i+3],
                     e[16*i+4], e[16*i+5], e[16*i+6], e[16*i+7]);

    EightinGreatFirH(d[16*i+8], d[16*i+9], d[16*i+10], d[16*i+11],
                     d[16*i+12], d[16*i+13], d[16*i+14], d[16*i+15],
                     e[16*i+8], e[16*i+9], e[16*i+10], e[16*i+11],
                     e[16*i+12], e[16*i+13], e[16*i+14], e[16*i+15]);
  }

  // ---------------- Stage 5 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/16; ++i) {
#pragma HLS UNROLL factor=2
    FourinSmallFirH(e[16*i+0], e[16*i+1], e[16*i+2], e[16*i+3],
                    f[16*i+0], f[16*i+1], f[16*i+2], f[16*i+3]);
    FourinSmallFirH(e[16*i+4], e[16*i+5], e[16*i+6], e[16*i+7],
                    f[16*i+4], f[16*i+5], f[16*i+6], f[16*i+7]);
    FourinGreatFirH(e[16*i+8], e[16*i+9], e[16*i+10], e[16*i+11],
                    f[16*i+8], f[16*i+9], f[16*i+10], f[16*i+11]);
    FourinGreatFirH(e[16*i+12], e[16*i+13], e[16*i+14], e[16*i+15],
                    f[16*i+12], f[16*i+13], f[16*i+14], f[16*i+15]);
  }

  // ---------------- Stage 6 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    bool g0 = AscendDescendKey(f[2*i], f[2*i+1]);
    g[2*i]   = g0 ? f[2*i+1] : f[2*i];
    g[2*i+1] = g0 ? f[2*i]   : f[2*i+1];
  }
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    bool g1 = AscendDescendKey(f[2*i+8], f[2*i+9]);
    g[2*i+8]  = g1 ? f[2*i+8] : f[2*i+9];
    g[2*i+9]  = g1 ? f[2*i+9] : f[2*i+8];
  }
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    bool g2 = AscendDescendKey(f[2*i+16], f[2*i+17]);
    g[2*i+16] = g2 ? f[2*i+17] : f[2*i+16];
    g[2*i+17] = g2 ? f[2*i+16] : f[2*i+17];
  }
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    bool g3 = AscendDescendKey(f[2*i+24], f[2*i+25]);
    g[2*i+24] = g3 ? f[2*i+24] : f[2*i+25];
    g[2*i+25] = g3 ? f[2*i+25] : f[2*i+24];
  }

  // ---------------- Stage 7 ----------------
  SixteenSmallFirH(g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7],
                   g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15],
                   h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7],
                   h[8], h[9], h[10], h[11], h[12], h[13], h[14], h[15]);
  SixteenGreatFirH(g[16], g[17], g[18], g[19], g[20], g[21], g[22], g[23],
                   g[24], g[25], g[26], g[27], g[28], g[29], g[30], g[31],
                   h[16], h[17], h[18], h[19], h[20], h[21], h[22], h[23],
                   h[24], h[25], h[26], h[27], h[28], h[29], h[30], h[31]);

  // ---------------- Stage 8 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/16; ++i) {
#pragma HLS UNROLL factor=2
    EightinSmallFirH(h[8*i+0], h[8*i+1], h[8*i+2], h[8*i+3],
                     h[8*i+4], h[8*i+5], h[8*i+6], h[8*i+7],
                     l[8*i+0], l[8*i+1], l[8*i+2], l[8*i+3],
                     l[8*i+4], l[8*i+5], l[8*i+6], l[8*i+7]);
  }
  for (int i = 0; i < MHH_SORT_SIZE/16; ++i) {
#pragma HLS UNROLL factor=2
    EightinGreatFirH(h[8*i+16], h[8*i+17], h[8*i+18], h[8*i+19],
                     h[8*i+20], h[8*i+21], h[8*i+22], h[8*i+23],
                     l[8*i+16], l[8*i+17], l[8*i+18], l[8*i+19],
                     l[8*i+20], l[8*i+21], l[8*i+22], l[8*i+23]);
  }

  // ---------------- Stage 9 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    FourinSmallFirH(l[4*i+0], l[4*i+1], l[4*i+2], l[4*i+3],
                    m[4*i+0], m[4*i+1], m[4*i+2], m[4*i+3]);
  }
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    FourinGreatFirH(l[4*i+16], l[4*i+17], l[4*i+18], l[4*i+19],
                    m[4*i+16], m[4*i+17], m[4*i+18], m[4*i+19]);
  }

  // ---------------- Stage 10 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/4; ++i) {
#pragma HLS UNROLL factor=2
    bool ga = AscendDescendKey(m[2*i], m[2*i+1]);
    n[2*i]   = ga ? m[2*i+1] : m[2*i];
    n[2*i+1] = ga ? m[2*i]   : m[2*i+1];
  }
  for (int i = 0; i < MHH_SORT_SIZE/4; ++i) {
#pragma HLS UNROLL factor=2
    bool gb = AscendDescendKey(m[2*i+16], m[2*i+17]);
    n[2*i+16] = gb ? m[2*i+16] : m[2*i+17];
    n[2*i+17] = gb ? m[2*i+17] : m[2*i+16];
  }

  // ---------------- Stage 11 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/2; ++i) {
#pragma HLS UNROLL factor=2
     bool gc = AscendDescendKey(n[i], n[i+16]);
    o[i]     = gc ? n[i+16] : n[i];
    o[i+16]  = gc ? n[i]    : n[i+16];
  }

  // ---------------- Stage 12 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/16; ++i) {
#pragma HLS UNROLL factor=2
    SixteenSmallFirH(o[16*i+0],  o[16*i+1],  o[16*i+2],  o[16*i+3],
                     o[16*i+4],  o[16*i+5],  o[16*i+6],  o[16*i+7],
                     o[16*i+8],  o[16*i+9],  o[16*i+10], o[16*i+11],
                     o[16*i+12], o[16*i+13], o[16*i+14], o[16*i+15],
                     p[16*i+0],  p[16*i+1],  p[16*i+2],  p[16*i+3],
                     p[16*i+4],  p[16*i+5],  p[16*i+6],  p[16*i+7],
                     p[16*i+8],  p[16*i+9],  p[16*i+10], p[16*i+11],
                     p[16*i+12], p[16*i+13], p[16*i+14], p[16*i+15]);
  }

  // ---------------- Stage 13 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/8; ++i) {
#pragma HLS UNROLL factor=2
    EightinSmallFirH(p[8*i+0], p[8*i+1], p[8*i+2], p[8*i+3],
                     p[8*i+4], p[8*i+5], p[8*i+6], p[8*i+7],
                     q[8*i+0], q[8*i+1], q[8*i+2], q[8*i+3],
                     q[8*i+4], q[8*i+5], q[8*i+6], q[8*i+7]);
  }

  // ---------------- Stage 14 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/4; ++i) {
#pragma HLS UNROLL factor=2
    FourinSmallFirH(q[4*i+0], q[4*i+1], q[4*i+2], q[4*i+3],
                    s[4*i+0], s[4*i+1], s[4*i+2], s[4*i+3]);
  }

  // ---------------- Stage 15 ----------------
  for (int i = 0; i < MHH_SORT_SIZE/2; ++i) {
    #pragma HLS UNROLL factor=2
    bool gd = AscendDescendKey(s[2*i], s[2*i+1]);
    r[2*i]   = gd ? s[2*i+1] : s[2*i];
    r[2*i+1] = gd ? s[2*i]   : s[2*i+1];
  }

  for (int i = 0; i < MHH_SORT_SIZE; ++i) {
//#pragma HLS PIPELINE II=1
    out[i] = in[r[i].idx];
  }
}

}  // namespace p2mhh

#endif
