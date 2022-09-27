#ifndef _PHASE_2_L1_CALO_PFCLUSTER_EMULATOR_H_
#define _PHASE_2_L1_CALO_PFCLUSTER_EMULATOR_H_

#include <cstdlib>

//    eta:  0  1  2  3  4   5  6  7  8   9 10 11 12  13 14 15 16  17 18 19 20
// 0             |                                                     |
// 1             |                                                     |
//               |-----------------------------------------------------|
// 2             |                                                     |
// 3             |                                                     |
// 4             |                                                     |
// 5             |                                                     |
//               | ----------------------------------------------------|
// 6             |                                                     |
// 7             |                                                     |
//
// 8 PFclusters are created in one 21x8 (2+17+2 x 2+4+2)

typedef struct {
float et ;
int eta ;
int phi ;
} GCTpfcluster_t ;

typedef struct {
GCTpfcluster_t GCTpfclusters[8] ;
} GCTPfcluster_t ;

typedef struct {
float et ;
int eta ;
int phi ;
} GCTint_t ;

typedef struct {
   GCTint_t t[8];
} gctEtaStrip_t ;

typedef struct {
   GCTint_t p[19];
} gctEtaStripPeak_t ;

typedef struct {
   gctEtaStrip_t s0;
   gctEtaStrip_t s1;
   gctEtaStrip_t s2;
   gctEtaStrip_t s3;
   gctEtaStrip_t s4;
   gctEtaStrip_t s5;
   gctEtaStrip_t s6;
   gctEtaStrip_t s7;
   gctEtaStrip_t s8;
   gctEtaStrip_t s9;
   gctEtaStrip_t s10;
   gctEtaStrip_t s11;
   gctEtaStrip_t s12;
   gctEtaStrip_t s13;
   gctEtaStrip_t s14;
   gctEtaStrip_t s15;
   gctEtaStrip_t s16;
   gctEtaStrip_t s17;
   gctEtaStrip_t s18;
   gctEtaStrip_t s19;
   gctEtaStrip_t s20;
} region_t;

GCTint_t bestOf2(const GCTint_t& t0, const GCTint_t& t1) {
  GCTint_t x;
  x = (t0.et > t1.et) ? t0 : t1;
  return x;
}

GCTint_t getPeakOfStrip(const gctEtaStrip_t& etaStrip){
  GCTint_t best12 = bestOf2(etaStrip.t[1],etaStrip.t[2]);
  GCTint_t best34 = bestOf2(etaStrip.t[3],etaStrip.t[4]);
  GCTint_t best56 = bestOf2(etaStrip.t[5],etaStrip.t[6]);
  GCTint_t best1234 = bestOf2(best12, best34);
  GCTint_t bestAll = bestOf2(best1234, best56);
  
  return bestAll;
}

GCTint_t getPeakBin(const gctEtaStripPeak_t& etaStripPeak){
  GCTint_t best01 = bestOf2(etaStripPeak.p[0], etaStripPeak.p[1]);
  GCTint_t best23 = bestOf2(etaStripPeak.p[2], etaStripPeak.p[3]);
  GCTint_t best45 = bestOf2(etaStripPeak.p[4], etaStripPeak.p[5]);
  GCTint_t best67 = bestOf2(etaStripPeak.p[6], etaStripPeak.p[7]);
  GCTint_t best89 = bestOf2(etaStripPeak.p[8], etaStripPeak.p[9]);
  GCTint_t best1011 = bestOf2(etaStripPeak.p[10], etaStripPeak.p[11]);
  GCTint_t best1213 = bestOf2(etaStripPeak.p[12], etaStripPeak.p[13]);
  GCTint_t best1415 = bestOf2(etaStripPeak.p[14], etaStripPeak.p[15]);
  GCTint_t best1617 = bestOf2(etaStripPeak.p[16], etaStripPeak.p[17]);
  GCTint_t best0123 = bestOf2(best01, best23);
  GCTint_t best4567 = bestOf2(best45, best67);
  GCTint_t best891011 = bestOf2(best89, best1011);
  GCTint_t best12131415 = bestOf2(best1213, best1415);
  GCTint_t best01234567 = bestOf2(best0123, best4567);
  GCTint_t best01234567891011 = bestOf2(best01234567, best891011);
  GCTint_t best121314151617 = bestOf2(best12131415, best1617);
  GCTint_t best12131415161718 = bestOf2(best121314151617, etaStripPeak.p[18]);
  GCTint_t bestAll = bestOf2(best01234567891011, best12131415161718);

  return bestAll;
}

GCTint_t getPeakPosition(const region_t& region){
  gctEtaStripPeak_t etaPeak;

  etaPeak.p[0] = getPeakOfStrip(region.s1);
  etaPeak.p[1] = getPeakOfStrip(region.s2);
  etaPeak.p[2] = getPeakOfStrip(region.s3);
  etaPeak.p[3] = getPeakOfStrip(region.s4);
  etaPeak.p[4] = getPeakOfStrip(region.s5);
  etaPeak.p[5] = getPeakOfStrip(region.s6);
  etaPeak.p[6] = getPeakOfStrip(region.s7);
  etaPeak.p[7] = getPeakOfStrip(region.s8);
  etaPeak.p[8] = getPeakOfStrip(region.s9);
  etaPeak.p[9] = getPeakOfStrip(region.s10);
  etaPeak.p[10] = getPeakOfStrip(region.s11);
  etaPeak.p[11] = getPeakOfStrip(region.s12);
  etaPeak.p[12] = getPeakOfStrip(region.s13);
  etaPeak.p[13] = getPeakOfStrip(region.s14);
  etaPeak.p[14] = getPeakOfStrip(region.s15);
  etaPeak.p[15] = getPeakOfStrip(region.s16);
  etaPeak.p[16] = getPeakOfStrip(region.s17);
  etaPeak.p[17] = getPeakOfStrip(region.s18);
  etaPeak.p[18] = getPeakOfStrip(region.s19);
  GCTint_t max = getPeakBin(etaPeak);

  return max;
}

region_t initStructure(float temp[21][8]){
  region_t r;

  for(int i=0; i<8; i++){
    r.s0.t[i].et = temp[0][i];
    r.s0.t[i].eta = 0;
    r.s0.t[i].phi = i;
    r.s1.t[i].et = temp[1][i];
    r.s1.t[i].eta = 1;
    r.s1.t[i].phi = i;
    r.s2.t[i].et = temp[2][i];
    r.s2.t[i].eta = 2;
    r.s2.t[i].phi = i;
    r.s3.t[i].et = temp[3][i];
    r.s3.t[i].eta = 3;
    r.s3.t[i].phi = i;
    r.s4.t[i].et = temp[4][i];
    r.s4.t[i].eta = 4;
    r.s4.t[i].phi = i;
    r.s5.t[i].et = temp[5][i];
    r.s5.t[i].eta = 5;
    r.s5.t[i].phi = i;
    r.s6.t[i].et = temp[6][i];
    r.s6.t[i].eta = 6;
    r.s6.t[i].phi = i;
    r.s7.t[i].et = temp[7][i];
    r.s7.t[i].eta = 7;
    r.s7.t[i].phi = i;
    r.s8.t[i].et = temp[8][i];
    r.s8.t[i].eta = 8;
    r.s8.t[i].phi = i;
    r.s9.t[i].et = temp[9][i];
    r.s9.t[i].eta = 9;
    r.s9.t[i].phi = i;
    r.s10.t[i].et = temp[10][i];
    r.s10.t[i].eta = 10;
    r.s10.t[i].phi = i;
    r.s11.t[i].et = temp[11][i];
    r.s11.t[i].eta = 11;
    r.s11.t[i].phi = i;
    r.s12.t[i].et = temp[12][i];
    r.s12.t[i].eta = 12;
    r.s12.t[i].phi = i;
    r.s13.t[i].et = temp[13][i];
    r.s13.t[i].eta = 13;
    r.s13.t[i].phi = i;
    r.s14.t[i].et = temp[14][i];
    r.s14.t[i].eta = 14;
    r.s14.t[i].phi = i;
    r.s15.t[i].et = temp[15][i];
    r.s15.t[i].eta = 15;
    r.s15.t[i].phi = i;
    r.s16.t[i].et = temp[16][i];
    r.s16.t[i].eta = 16;
    r.s16.t[i].phi = i;
    r.s17.t[i].et = temp[17][i];
    r.s17.t[i].eta = 17;
    r.s17.t[i].phi = i;
    r.s18.t[i].et = temp[18][i];
    r.s18.t[i].eta = 18;
    r.s18.t[i].phi = i;
    r.s19.t[i].et = temp[19][i];
    r.s19.t[i].eta = 19;
    r.s19.t[i].phi = i;
    r.s20.t[i].et = temp[20][i];
    r.s20.t[i].eta = 20;
    r.s20.t[i].phi = i;

  }
  return r;
}

float getEt(float temp[21][8], int eta, int phi){
  float et_sumEta[3];

  for(int i=0; i<19; i++){
    for(int j=0; j<6; j++){
      if (i+1 == eta && j+1 == phi){
        for(int k=0; k<3; k++){
          et_sumEta[k] = temp[i+k][j] + temp[i+k][j+1] + temp[i+k][j+2];
        }
      }
    }
  }

  float pfcluster_et = et_sumEta[0] + et_sumEta[1] + et_sumEta[2];

  return pfcluster_et;
}

void RemoveTmp(float temp[21][8], int eta, int phi){
 
  for(int i=0; i<21; i++){
    if(i+1 >= eta && i <= eta+1) {
      for(int j=0; j<8; j++){
        if(j+1 >= phi && j <= phi+1) {
          temp[i][j] = 0;
        }
      }
    }
  }
  return;
}

GCTpfcluster_t recoPfcluster(float temporary[21][8], int etaoffset, int phioffset){
  GCTpfcluster_t pfclusterReturn;

  region_t region;

  region = initStructure(temporary);

  GCTint_t regionMax = getPeakPosition(region);  

  float pfcluster_et = getEt(temporary, regionMax.eta, regionMax.phi);

  RemoveTmp(temporary, regionMax.eta, regionMax.phi);

  if(!(regionMax.eta >= 2 && regionMax.eta <= 18 && regionMax.phi >= 2 && regionMax.phi <= 5)) pfcluster_et = 0;

  pfclusterReturn.et  = pfcluster_et;
  pfclusterReturn.eta = regionMax.eta - 2 + etaoffset;
  pfclusterReturn.phi = regionMax.phi - 2 + phioffset;

  return pfclusterReturn ;
}

GCTPfcluster_t pfcluster(float temporary[21][8], int etaoffset, int phioffset){
  GCTpfcluster_t pfcluster[8];

  for(int i=0; i<8; i++){
    pfcluster[i] = recoPfcluster(temporary, etaoffset, phioffset);
  }

  GCTPfcluster_t GCTPfclusters;

  for(int i=0; i<8; i++){
    GCTPfclusters.GCTpfclusters[i].et = pfcluster[i].et;
    GCTPfclusters.GCTpfclusters[i].eta = pfcluster[i].eta;
    GCTPfclusters.GCTpfclusters[i].phi = pfcluster[i].phi;
  }

  return GCTPfclusters;
}

#endif
