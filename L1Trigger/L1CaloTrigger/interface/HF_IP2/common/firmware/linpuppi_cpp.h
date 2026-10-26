#include "linpuppi.h"
#include <algorithm>
#include <cassert>
#include <cstdio>

#ifndef __SYNTHESIS__
int gdebug_;
void linpuppi_set_debug(bool debug) { gdebug_ = debug; }
#else
void linpuppi_set_debug(bool debug) {}
#endif

using namespace l1ct;
using namespace linpuppi;

inline ap_uint<eta_t::width*2+1> dr2_int_fw(eta_t eta1, phi_t phi1, eta_t eta2, phi_t phi2) {
    ap_int<eta_t::width+1> deta =eta1-eta2;
    ap_int<phi_t::width+1> dphi =phi1-phi2;
    auto deta2 = deta*deta;
    auto dphi2 = dphi*dphi;
//    #pragma HLS RESOURCE variable=deta2 core=Mul_LUT
//    #pragma HLS RESOURCE variable=dphi2 core=Mul_LUT
    return deta2 + dphi2;
}

#define linpuppi_dr2_invert_table_size 512
void _lut_dr2_invert_init(dr2inv_t _table[linpuppi_dr2_invert_table_size]) {
    _table[0] = 1.0f;
    for (int i = 1; i < linpuppi_dr2_invert_table_size; ++i) {
        _table[i] = (1.0f / i);
    }
}

sumTerm_t _lut_dr2_divide(ap_uint<17> num, ap_uint<9> den) {
    static_assert(linpuppi_dr2_invert_table_size == (1 << decltype(den)::width), "Width must be 9 bits");
    static dr2inv_t _table[linpuppi_dr2_invert_table_size];
    _lut_dr2_invert_init(_table);
//    sumTerm_t inv;
//    inv = num* _table[den];
//    return inv;
  return (num * _table[den]);
}


#define fwdlinpuppi_init_x2a_table_size 1024

void fwdlinpuppi_init_x2a_short(x2_t table[fwdlinpuppi_init_x2a_table_size]) {
    const alphaSlope_t alphaSlope = LINPUPPI_alphaSlope * std::log(2);
    const alpha_t alphaZero = LINPUPPI_alphaZero / std::log(2);
    const x2_t logOffset = std::log2(linpuppi::PT2DR2_LSB) - SUM_BITSHIFT;
    for (int i = 0; i < fwdlinpuppi_init_x2a_table_size; ++i) {
        const x2_t C0 = -alphaSlope * alphaZero;
        const x2_t C1 = alphaSlope * logOffset;
        table[i] = C0 + (i > 0 ? x2_t(alphaSlope * alpha_t(std::log2(float(i))) + C1) : x2_t(0));
    }
}

x2_t fwdlinpuppi_calc_x2a_step2(sumTerm_t sum, alphaSlope_t alphaSlope, x2_t alphaCrop, const x2_t table[fwdlinpuppi_init_x2a_table_size]) {
    const int log2lut_bits = 10;
    ap_uint<sumTerm_t::width> rawsum;
    rawsum(sumTerm_t::width-1,0) = sum(sumTerm_t::width-1,0);
    x2_t sumterm = 0;
    auto logarg = rawsum;
    for (int b = decltype(rawsum)::width-1-log2lut_bits; b >= 0; --b) {
        if (rawsum[b+log2lut_bits]) {
            logarg = logarg >> (b + 1);
            sumterm = (b + 1) * alphaSlope;
            break;
        }
    }
    assert(logarg >= 0 && logarg < fwdlinpuppi_init_x2a_table_size);
    x2_t ret = (table[logarg] + sumterm);
    if (ret < -alphaCrop) return -alphaCrop;
    else if (ret > alphaCrop) return alphaCrop;
    else return ret;
}

x2_t fwdlinpuppi_calc_x2a(sumTerm_t sum) {
    static x2_t table[fwdlinpuppi_init_x2a_table_size];
#ifdef __SYNTHESIS__
    fwdlinpuppi_init_x2a_short(table);
#else	// initialize the table only once, otherwise this is really slow
    static bool is_init = false;
    if (!is_init) { fwdlinpuppi_init_x2a_short(table); is_init = true; }
#endif
    // we put a log(2) here since we compute alpha as log2(sum) instead of ln(sum)
    const alphaSlope_t alphaSlope = LINPUPPI_alphaSlope * std::log(2);
    const alpha_t alphaCrop = LINPUPPI_alphaCrop;
    return fwdlinpuppi_calc_x2a_step2(sum, alphaSlope, alphaCrop, table);
}


#define fwdlinpuppi_x2w_table_size 1024

void fwdlinpuppi_init_w(puppiWgt_t table[fwdlinpuppi_x2w_table_size]) {
    const int xavg = fwdlinpuppi_x2w_table_size/2;
    const float X2_LSB = 1.0f/(1<<(x2_t::width-x2_t::iwidth));  // = 1/32 for x2_t=ap_fixed<12,7>
    for (int i = 0; i < fwdlinpuppi_x2w_table_size; ++i) {
        int x2 = i - xavg;
        puppiWgt_t val = 1.0/(1.0 + std::exp(-x2*X2_LSB));
        table[i] = val;
    }
}

void linpuppi_calc_wpt(pt_t pt, x2_t x2, pt_t & puppiPt, puppiWgt_t & wgt) {
    static puppiWgt_t table[fwdlinpuppi_x2w_table_size];
#ifdef __SYNTHESIS__
    fwdlinpuppi_init_w(table);
#else
    static bool is_init = false;
    if (!is_init) { fwdlinpuppi_init_w(table); is_init = true; }
#endif

    const ap_uint<10> xavg = fwdlinpuppi_x2w_table_size>>1;
    ap_int<x2_t::width> x2_raw;
    x2_raw(x2_t::width-1,0) = x2(x2_t::width-1,0);
    ap_uint<10> index = 0;
    if (x2_raw < -xavg) index = 0;
    else if (x2_raw >= xavg) index = fwdlinpuppi_x2w_table_size-1;
    else index = x2_raw + xavg;
    wgt = table[index];
    typedef ap_fixed<pt_t::width, pt_t::iwidth, AP_RND, AP_SAT> pt_rounding_t;
    puppiPt = pt_rounding_t(pt * table[index]);
//    dont pragma
}


void fwdlinpuppiSum(const HadCaloObj caloin[NCALO], sumTerm_t sums[NCALO]) {

    const int DR2MAX = LINPUPPI_DR2MAX;
    const int DR2MIN = LINPUPPI_DR2MIN;
    const int DR2MIN_SHIFT = DR2MIN >> 5;
    const ap_uint<24> PTMAX2_SHIFT = (LINPUPPI_iptMax)*(LINPUPPI_iptMax) >> 5;
    ap_uint<17> pt2_shift[NCALO];
    for (int it = 0; it < NCALO; ++it) {
        auto ipt = Scales::ptToInt(caloin[it].hwPt);
        ap_uint<24> mypt2 = (ipt*ipt)>>5;
//        ap_uint<24> mypt2_ = ipt*ipt;
//        #pragma HLS RESOURCE variable=mypt2_ core=Mul_LUT latency=2
//        ap_uint<24> mypt2 = mypt2_>>5;
        pt2_shift[it] = (mypt2 < PTMAX2_SHIFT ? mypt2 : PTMAX2_SHIFT);
    }

    for (int in = 0; in < NCALO; ++in) {
        sumTerm_t sum = 0;
        for (int it = 0; it < NCALO; ++it) {
            if (it == in) continue;
            auto dr2 = dr2_int_fw(caloin[it].hwEta, caloin[it].hwPhi, caloin[in].hwEta, caloin[in].hwPhi);
            if (dr2 <= DR2MAX) {
                ap_uint<9> dr2short = dr2 >> 5;
                if (dr2short < DR2MIN_SHIFT) dr2short = DR2MIN_SHIFT;
                auto term = _lut_dr2_divide(pt2_shift[it], dr2short);
                sum += term;
            }
        }
        sums[in] = sum;
    }
}



void fwdlinpuppiSum2Pt(const HadCaloObj caloin[NCALO], const sumTerm_t sums[NCALO], pt_t puppiPts[NCALO], puppiWgt_t puppiWgts[NCALO]) {
    const ptSlope_t ptSlopeNe = LINPUPPI_ptSlopeNe;
    const ptSlope_t ptSlopePh = LINPUPPI_ptSlopePh;
    const pt_t ptZeroNe = LINPUPPI_ptZeroNe;
    const pt_t ptZeroPh = LINPUPPI_ptZeroPh;
    const x2_t priorNe = LINPUPPI_priorNe;
    const x2_t priorPh = LINPUPPI_priorPh;
    const pt_t ptCut = LINPUPPI_ptCut;

    x2_t x2a[NCALO], x2ptp[NCALO];

    for (int in = 0; in < NCALO; ++in) {
        x2a[in] = fwdlinpuppi_calc_x2a(sums[in]);
    }

    for (int in = 0; in < NCALO; ++in) {
        if (caloin[in].hwIsEM()) {
            x2ptp[in] = ptSlopePh*(caloin[in].hwPt - ptZeroPh) - priorPh;

        } else {
            x2ptp[in] = ptSlopeNe*(caloin[in].hwPt - ptZeroNe) - priorNe;
        }
    }

    for (int in = 0; in < NCALO; ++in) {
        x2_t x2 = x2a[in] + x2ptp[in];
        linpuppi_calc_wpt(caloin[in].hwPt, x2, puppiPts[in], puppiWgts[in]);
    }
}



void fwdlinpuppiPt(const HadCaloObj caloin[NCALO], pt_t puppiPts[NCALO], puppiWgt_t puppiWgts[NCALO]) {
    sumTerm_t sums[NCALO];
    fwdlinpuppiSum(caloin, sums);
    fwdlinpuppiSum2Pt(caloin, sums, puppiPts, puppiWgts);
}


void fwdlinpuppi(const PFRegion & region, const HadCaloObj caloin[NCALO], PuppiObj pfselne[NNEUTRALS]) {
    pt_t puppiPts[NCALO];
    puppiWgt_t puppiWgts[NCALO];

    fwdlinpuppiPt(caloin, puppiPts, puppiWgts);

    PuppiObj work[NNEUTRALS];

    for (int iout = 0; iout < NNEUTRALS; ++iout) {
        work[iout].clear();
    }

    const pt_t ptCut = LINPUPPI_ptCut;
//    bool valid[NCALO];
//    	#pragma HLS ARRAY_PARTITION variable=valid complete
//        for (int in = 0; in < NCALO; ++in) {
//    		#pragma HLS UNROLL
//            valid[in] = !region.isFiducial(caloin[in]) || puppiPts[in] < ptCut;
//        }
//        for (int in = 0; in < NCALO; ++in) {
//            if (valid[in]) continue;
//            for (int iout = NNEUTRALS-1; iout >= 0; --iout) {
//                if (work[iout].hwPt <= puppiPts[in]) {
//                    if (iout == 0 || work[iout-1].hwPt > puppiPts[in]) {
//                        work[iout].fill(region, caloin[in], puppiPts[in], puppiWgts[in]);
//                    } else {
//                        work[iout] = work[iout-1];
//                    }
//                }
//            }
//        }
        for (int in = 0; in < NCALO; ++in) {
            if (!region.isFiducial(caloin[in]) || puppiPts[in] < ptCut) continue;
            for (int iout = NNEUTRALS-1; iout >= 0; --iout) {
                if (work[iout].hwPt <= puppiPts[in]) {
                    if (iout == 0 || work[iout-1].hwPt > puppiPts[in]) {
                        work[iout].fill(region, caloin[in], puppiPts[in], puppiWgts[in]);
                    } else {
                        work[iout] = work[iout-1];
                    }
                }
            }
        }

    for (int iout = 0; iout < NNEUTRALS; ++iout) {
        pfselne[iout] = work[iout];
    }
}
