//------------------------------------
// Bitonic sort implementation for GCT SumCard emulator
// (based heavily on bitonicSort32 in firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_bitonicSort32_GCT_cpp
#define L1Trigger_L1CaloTrigger_bitonicSort32_GCT_cpp

#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_GCT_h.h"

namespace p2gctsum {

inline bool AscendDescendKey(Handle& x, Handle& y) { return x.key > y.key; }

inline void FourinSmallFirH(Handle& x0, Handle& x1, Handle& x2, Handle& x3, Handle& y0, Handle& y1, Handle& y2, Handle& y3) {
  const bool g02 = AscendDescendKey(x0, x2);
  y0 = g02 ? x2 : x0;
  y2 = g02 ? x0 : x2;
  const bool g13 = AscendDescendKey(x1, x3);
  y1 = g13 ? x3 : x1;
  y3 = g13 ? x1 : x3;
}

inline void FourinGreatFirH(Handle& x0, Handle& x1, Handle& x2, Handle& x3, Handle& y0, Handle& y1, Handle& y2, Handle& y3) {
  const bool g02 = AscendDescendKey(x0, x2);
  y0 = g02 ? x0 : x2;
  y2 = g02 ? x2 : x0;
  const bool g13 = AscendDescendKey(x1, x3);
  y1 = g13 ? x1 : x3;
  y3 = g13 ? x3 : x1;
}

inline void EightinSmallFirH(Handle& x0,
                             Handle& x1,
                             Handle& x2,
                             Handle& x3,
                             Handle& x4,
                             Handle& x5,
                             Handle& x6,
                             Handle& x7,
                             Handle& y0,
                             Handle& y1,
                             Handle& y2,
                             Handle& y3,
                             Handle& y4,
                             Handle& y5,
                             Handle& y6,
                             Handle& y7) {
  const bool g04 = AscendDescendKey(x0, x4);
  y0 = g04 ? x4 : x0;
  y4 = g04 ? x0 : x4;
  const bool g15 = AscendDescendKey(x1, x5);
  y1 = g15 ? x5 : x1;
  y5 = g15 ? x1 : x5;
  const bool g26 = AscendDescendKey(x2, x6);
  y2 = g26 ? x6 : x2;
  y6 = g26 ? x2 : x6;
  const bool g37 = AscendDescendKey(x3, x7);
  y3 = g37 ? x7 : x3;
  y7 = g37 ? x3 : x7;
}

inline void EightinGreatFirH(Handle& x0,
                             Handle& x1,
                             Handle& x2,
                             Handle& x3,
                             Handle& x4,
                             Handle& x5,
                             Handle& x6,
                             Handle& x7,
                             Handle& y0,
                             Handle& y1,
                             Handle& y2,
                             Handle& y3,
                             Handle& y4,
                             Handle& y5,
                             Handle& y6,
                             Handle& y7) {
  const bool g04 = AscendDescendKey(x0, x4);
  y0 = g04 ? x0 : x4;
  y4 = g04 ? x4 : x0;
  const bool g15 = AscendDescendKey(x1, x5);
  y1 = g15 ? x1 : x5;
  y5 = g15 ? x5 : x1;
  const bool g26 = AscendDescendKey(x2, x6);
  y2 = g26 ? x2 : x6;
  y6 = g26 ? x6 : x2;
  const bool g37 = AscendDescendKey(x3, x7);
  y3 = g37 ? x3 : x7;
  y7 = g37 ? x7 : x3;
}

inline void SixteenSmallFirH(Handle& x0,
                             Handle& x1,
                             Handle& x2,
                             Handle& x3,
                             Handle& x4,
                             Handle& x5,
                             Handle& x6,
                             Handle& x7,
                             Handle& x8,
                             Handle& x9,
                             Handle& x10,
                             Handle& x11,
                             Handle& x12,
                             Handle& x13,
                             Handle& x14,
                             Handle& x15,
                             Handle& y0,
                             Handle& y1,
                             Handle& y2,
                             Handle& y3,
                             Handle& y4,
                             Handle& y5,
                             Handle& y6,
                             Handle& y7,
                             Handle& y8,
                             Handle& y9,
                             Handle& y10,
                             Handle& y11,
                             Handle& y12,
                             Handle& y13,
                             Handle& y14,
                             Handle& y15) {
  const bool g0 = AscendDescendKey(x0, x8);
  y0 = g0 ? x8 : x0;
  y8 = g0 ? x0 : x8;
  const bool g1 = AscendDescendKey(x1, x9);
  y1 = g1 ? x9 : x1;
  y9 = g1 ? x1 : x9;
  const bool g2 = AscendDescendKey(x2, x10);
  y2 = g2 ? x10 : x2;
  y10 = g2 ? x2 : x10;
  const bool g3 = AscendDescendKey(x3, x11);
  y3 = g3 ? x11 : x3;
  y11 = g3 ? x3 : x11;
  const bool g4 = AscendDescendKey(x4, x12);
  y4 = g4 ? x12 : x4;
  y12 = g4 ? x4 : x12;
  const bool g5 = AscendDescendKey(x5, x13);
  y5 = g5 ? x13 : x5;
  y13 = g5 ? x5 : x13;
  const bool g6 = AscendDescendKey(x6, x14);
  y6 = g6 ? x14 : x6;
  y14 = g6 ? x6 : x14;
  const bool g7 = AscendDescendKey(x7, x15);
  y7 = g7 ? x15 : x7;
  y15 = g7 ? x7 : x15;
}

inline void SixteenGreatFirH(Handle& x0,
                             Handle& x1,
                             Handle& x2,
                             Handle& x3,
                             Handle& x4,
                             Handle& x5,
                             Handle& x6,
                             Handle& x7,
                             Handle& x8,
                             Handle& x9,
                             Handle& x10,
                             Handle& x11,
                             Handle& x12,
                             Handle& x13,
                             Handle& x14,
                             Handle& x15,
                             Handle& y0,
                             Handle& y1,
                             Handle& y2,
                             Handle& y3,
                             Handle& y4,
                             Handle& y5,
                             Handle& y6,
                             Handle& y7,
                             Handle& y8,
                             Handle& y9,
                             Handle& y10,
                             Handle& y11,
                             Handle& y12,
                             Handle& y13,
                             Handle& y14,
                             Handle& y15) {
  const bool g0 = AscendDescendKey(x0, x8);
  y0 = g0 ? x0 : x8;
  y8 = g0 ? x8 : x0;
  const bool g1 = AscendDescendKey(x1, x9);
  y1 = g1 ? x1 : x9;
  y9 = g1 ? x9 : x1;
  const bool g2 = AscendDescendKey(x2, x10);
  y2 = g2 ? x2 : x10;
  y10 = g2 ? x10 : x2;
  const bool g3 = AscendDescendKey(x3, x11);
  y3 = g3 ? x3 : x11;
  y11 = g3 ? x11 : x3;
  const bool g4 = AscendDescendKey(x4, x12);
  y4 = g4 ? x4 : x12;
  y12 = g4 ? x12 : x4;
  const bool g5 = AscendDescendKey(x5, x13);
  y5 = g5 ? x5 : x13;
  y13 = g5 ? x13 : x5;
  const bool g6 = AscendDescendKey(x6, x14);
  y6 = g6 ? x6 : x14;
  y14 = g6 ? x14 : x6;
  const bool g7 = AscendDescendKey(x7, x15);
  y7 = g7 ? x7 : x15;
  y15 = g7 ? x15 : x7;
}

inline void bitonicSort32(GCTvar in[32], GCTvar out[32]) {
  static constexpr int N = 32;

  // Preserve the firmware's exact comparator network so equal-ET ties keep
  // the same object ordering as the HLS implementation.
  Handle a[N], b[N], c[N], d[N], e[N], f[N], g[N], h[N], l[N], m[N], n[N], o[N], p[N], q[N], s[N], r[N];

  for (int i = 0; i < 32; ++i) {
    a[i].key = in[i].ET;
    a[i].idx = i;
  }

  for (int i = 0; i < N / 4; ++i) {
    bool g01 = AscendDescendKey(a[4 * i], a[4 * i + 1]);
    b[4 * i] = g01 ? a[4 * i + 1] : a[4 * i];
    b[4 * i + 1] = g01 ? a[4 * i] : a[4 * i + 1];
  }
  for (int i = 0; i < N / 4; ++i) {
    bool g23 = AscendDescendKey(a[4 * i + 2], a[4 * i + 3]);
    b[4 * i + 2] = g23 ? a[4 * i + 2] : a[4 * i + 3];
    b[4 * i + 3] = g23 ? a[4 * i + 3] : a[4 * i + 2];
  }

  for (int i = 0; i < N / 8; ++i) {
    FourinSmallFirH(b[8 * i + 0], b[8 * i + 1], b[8 * i + 2], b[8 * i + 3], c[8 * i + 0], c[8 * i + 1], c[8 * i + 2],
                    c[8 * i + 3]);
    FourinGreatFirH(b[8 * i + 4], b[8 * i + 5], b[8 * i + 6], b[8 * i + 7], c[8 * i + 4], c[8 * i + 5], c[8 * i + 6],
                    c[8 * i + 7]);
  }

  for (int i = 0; i < N / 8; ++i) {
    bool g01 = AscendDescendKey(c[8 * i + 0], c[8 * i + 1]);
    d[8 * i + 0] = g01 ? c[8 * i + 1] : c[8 * i + 0];
    d[8 * i + 1] = g01 ? c[8 * i + 0] : c[8 * i + 1];

    bool g23 = AscendDescendKey(c[8 * i + 2], c[8 * i + 3]);
    d[8 * i + 2] = g23 ? c[8 * i + 3] : c[8 * i + 2];
    d[8 * i + 3] = g23 ? c[8 * i + 2] : c[8 * i + 3];
  }
  for (int i = 0; i < N / 8; ++i) {
    bool g45 = AscendDescendKey(c[8 * i + 4], c[8 * i + 5]);
    d[8 * i + 4] = g45 ? c[8 * i + 4] : c[8 * i + 5];
    d[8 * i + 5] = g45 ? c[8 * i + 5] : c[8 * i + 4];

    bool g67 = AscendDescendKey(c[8 * i + 6], c[8 * i + 7]);
    d[8 * i + 6] = g67 ? c[8 * i + 6] : c[8 * i + 7];
    d[8 * i + 7] = g67 ? c[8 * i + 7] : c[8 * i + 6];
  }

  for (int i = 0; i < N / 16; ++i) {
    EightinSmallFirH(d[16 * i + 0],
                     d[16 * i + 1],
                     d[16 * i + 2],
                     d[16 * i + 3],
                     d[16 * i + 4],
                     d[16 * i + 5],
                     d[16 * i + 6],
                     d[16 * i + 7],
                     e[16 * i + 0],
                     e[16 * i + 1],
                     e[16 * i + 2],
                     e[16 * i + 3],
                     e[16 * i + 4],
                     e[16 * i + 5],
                     e[16 * i + 6],
                     e[16 * i + 7]);
    EightinGreatFirH(d[16 * i + 8],
                     d[16 * i + 9],
                     d[16 * i + 10],
                     d[16 * i + 11],
                     d[16 * i + 12],
                     d[16 * i + 13],
                     d[16 * i + 14],
                     d[16 * i + 15],
                     e[16 * i + 8],
                     e[16 * i + 9],
                     e[16 * i + 10],
                     e[16 * i + 11],
                     e[16 * i + 12],
                     e[16 * i + 13],
                     e[16 * i + 14],
                     e[16 * i + 15]);
  }

  for (int i = 0; i < N / 16; ++i) {
    FourinSmallFirH(e[16 * i + 0], e[16 * i + 1], e[16 * i + 2], e[16 * i + 3], f[16 * i + 0], f[16 * i + 1], f[16 * i + 2],
                    f[16 * i + 3]);
    FourinSmallFirH(e[16 * i + 4], e[16 * i + 5], e[16 * i + 6], e[16 * i + 7], f[16 * i + 4], f[16 * i + 5], f[16 * i + 6],
                    f[16 * i + 7]);
    FourinGreatFirH(e[16 * i + 8], e[16 * i + 9], e[16 * i + 10], e[16 * i + 11], f[16 * i + 8], f[16 * i + 9], f[16 * i + 10],
                    f[16 * i + 11]);
    FourinGreatFirH(e[16 * i + 12], e[16 * i + 13], e[16 * i + 14], e[16 * i + 15], f[16 * i + 12], f[16 * i + 13],
                    f[16 * i + 14], f[16 * i + 15]);
  }

  for (int i = 0; i < N / 8; ++i) {
    bool g0 = AscendDescendKey(f[2 * i], f[2 * i + 1]);
    g[2 * i] = g0 ? f[2 * i + 1] : f[2 * i];
    g[2 * i + 1] = g0 ? f[2 * i] : f[2 * i + 1];
  }
  for (int i = 0; i < N / 8; ++i) {
    bool g1 = AscendDescendKey(f[2 * i + 8], f[2 * i + 9]);
    g[2 * i + 8] = g1 ? f[2 * i + 8] : f[2 * i + 9];
    g[2 * i + 9] = g1 ? f[2 * i + 9] : f[2 * i + 8];
  }
  for (int i = 0; i < N / 8; ++i) {
    bool g2 = AscendDescendKey(f[2 * i + 16], f[2 * i + 17]);
    g[2 * i + 16] = g2 ? f[2 * i + 17] : f[2 * i + 16];
    g[2 * i + 17] = g2 ? f[2 * i + 16] : f[2 * i + 17];
  }
  for (int i = 0; i < N / 8; ++i) {
    bool g3 = AscendDescendKey(f[2 * i + 24], f[2 * i + 25]);
    g[2 * i + 24] = g3 ? f[2 * i + 24] : f[2 * i + 25];
    g[2 * i + 25] = g3 ? f[2 * i + 25] : f[2 * i + 24];
  }

  SixteenSmallFirH(g[0],
                   g[1],
                   g[2],
                   g[3],
                   g[4],
                   g[5],
                   g[6],
                   g[7],
                   g[8],
                   g[9],
                   g[10],
                   g[11],
                   g[12],
                   g[13],
                   g[14],
                   g[15],
                   h[0],
                   h[1],
                   h[2],
                   h[3],
                   h[4],
                   h[5],
                   h[6],
                   h[7],
                   h[8],
                   h[9],
                   h[10],
                   h[11],
                   h[12],
                   h[13],
                   h[14],
                   h[15]);
  SixteenGreatFirH(g[16],
                   g[17],
                   g[18],
                   g[19],
                   g[20],
                   g[21],
                   g[22],
                   g[23],
                   g[24],
                   g[25],
                   g[26],
                   g[27],
                   g[28],
                   g[29],
                   g[30],
                   g[31],
                   h[16],
                   h[17],
                   h[18],
                   h[19],
                   h[20],
                   h[21],
                   h[22],
                   h[23],
                   h[24],
                   h[25],
                   h[26],
                   h[27],
                   h[28],
                   h[29],
                   h[30],
                   h[31]);

  for (int i = 0; i < N / 16; ++i) {
    EightinSmallFirH(h[8 * i + 0],
                     h[8 * i + 1],
                     h[8 * i + 2],
                     h[8 * i + 3],
                     h[8 * i + 4],
                     h[8 * i + 5],
                     h[8 * i + 6],
                     h[8 * i + 7],
                     l[8 * i + 0],
                     l[8 * i + 1],
                     l[8 * i + 2],
                     l[8 * i + 3],
                     l[8 * i + 4],
                     l[8 * i + 5],
                     l[8 * i + 6],
                     l[8 * i + 7]);
  }
  for (int i = 0; i < N / 16; ++i) {
    EightinGreatFirH(h[8 * i + 16],
                     h[8 * i + 17],
                     h[8 * i + 18],
                     h[8 * i + 19],
                     h[8 * i + 20],
                     h[8 * i + 21],
                     h[8 * i + 22],
                     h[8 * i + 23],
                     l[8 * i + 16],
                     l[8 * i + 17],
                     l[8 * i + 18],
                     l[8 * i + 19],
                     l[8 * i + 20],
                     l[8 * i + 21],
                     l[8 * i + 22],
                     l[8 * i + 23]);
  }

  for (int i = 0; i < N / 8; ++i) {
    FourinSmallFirH(l[4 * i + 0], l[4 * i + 1], l[4 * i + 2], l[4 * i + 3], m[4 * i + 0], m[4 * i + 1], m[4 * i + 2],
                    m[4 * i + 3]);
  }
  for (int i = 0; i < N / 8; ++i) {
    FourinGreatFirH(l[4 * i + 16], l[4 * i + 17], l[4 * i + 18], l[4 * i + 19], m[4 * i + 16], m[4 * i + 17],
                    m[4 * i + 18], m[4 * i + 19]);
  }

  for (int i = 0; i < N / 4; ++i) {
    bool ga = AscendDescendKey(m[2 * i], m[2 * i + 1]);
    n[2 * i] = ga ? m[2 * i + 1] : m[2 * i];
    n[2 * i + 1] = ga ? m[2 * i] : m[2 * i + 1];
  }
  for (int i = 0; i < N / 4; ++i) {
    bool gb = AscendDescendKey(m[2 * i + 16], m[2 * i + 17]);
    n[2 * i + 16] = gb ? m[2 * i + 16] : m[2 * i + 17];
    n[2 * i + 17] = gb ? m[2 * i + 17] : m[2 * i + 16];
  }

  for (int i = 0; i < N / 2; ++i) {
    bool gc = AscendDescendKey(n[i], n[i + 16]);
    o[i] = gc ? n[i + 16] : n[i];
    o[i + 16] = gc ? n[i] : n[i + 16];
  }

  for (int i = 0; i < N / 16; ++i) {
    SixteenSmallFirH(o[16 * i + 0],
                     o[16 * i + 1],
                     o[16 * i + 2],
                     o[16 * i + 3],
                     o[16 * i + 4],
                     o[16 * i + 5],
                     o[16 * i + 6],
                     o[16 * i + 7],
                     o[16 * i + 8],
                     o[16 * i + 9],
                     o[16 * i + 10],
                     o[16 * i + 11],
                     o[16 * i + 12],
                     o[16 * i + 13],
                     o[16 * i + 14],
                     o[16 * i + 15],
                     p[16 * i + 0],
                     p[16 * i + 1],
                     p[16 * i + 2],
                     p[16 * i + 3],
                     p[16 * i + 4],
                     p[16 * i + 5],
                     p[16 * i + 6],
                     p[16 * i + 7],
                     p[16 * i + 8],
                     p[16 * i + 9],
                     p[16 * i + 10],
                     p[16 * i + 11],
                     p[16 * i + 12],
                     p[16 * i + 13],
                     p[16 * i + 14],
                     p[16 * i + 15]);
  }

  for (int i = 0; i < N / 8; ++i) {
    EightinSmallFirH(p[8 * i + 0],
                     p[8 * i + 1],
                     p[8 * i + 2],
                     p[8 * i + 3],
                     p[8 * i + 4],
                     p[8 * i + 5],
                     p[8 * i + 6],
                     p[8 * i + 7],
                     q[8 * i + 0],
                     q[8 * i + 1],
                     q[8 * i + 2],
                     q[8 * i + 3],
                     q[8 * i + 4],
                     q[8 * i + 5],
                     q[8 * i + 6],
                     q[8 * i + 7]);
  }

  for (int i = 0; i < N / 4; ++i) {
    FourinSmallFirH(q[4 * i + 0], q[4 * i + 1], q[4 * i + 2], q[4 * i + 3], s[4 * i + 0], s[4 * i + 1], s[4 * i + 2],
                    s[4 * i + 3]);
  }

  for (int i = 0; i < N / 2; ++i) {
    bool gd = AscendDescendKey(s[2 * i], s[2 * i + 1]);
    r[2 * i] = gd ? s[2 * i + 1] : s[2 * i];
    r[2 * i + 1] = gd ? s[2 * i] : s[2 * i + 1];
  }

  for (int i = 0; i < N; ++i) {
    out[i] = in[r[i].idx];
  }
}

}  // namespace p2gctsum

#endif
