#ifndef PUPPI_LINPUPPI_H
#define PUPPI_LINPUPPI_H
#define REG_HF
#define HLS_pipeline_II 1


#include <cmath>
#include "../layer1_objs.h"
#include "../pf.h"
#include "../puppi.h"
#include "../layer1_multiplicities.h"
#include "linpuppi_bits.h"


// neutrals, forward
void fwdlinpuppi(const l1ct::PFRegion & region, const l1ct::HadCaloObj caloin[NCALO], l1ct::PuppiObj pfselne[NNEUTRALS]);
void fwdlinpuppiNoCrop(const l1ct::PFRegion & region, const l1ct::HadCaloObj caloin[NCALO], l1ct::PuppiObj pfallne[NCALO]);

#define LINPUPPI_DATA_SIZE_FWD 64
#define LINPUPPI_NCHANN_FWD_OUTNC (NCALO)
#define LINPUPPI_NCHANN_FWD_OUT (NNEUTRALS)

void packed_fwdlinpuppi(const ap_uint<l1ct::PFRegion::BITWIDTH> input_reg, const ap_uint<l1ct::HadCaloObj::BITWIDTH_SLIM> input_calo[NCALO], ap_uint<LINPUPPI_DATA_SIZE_FWD> output[LINPUPPI_NCHANN_FWD_OUT]);
void packed_fwdlinpuppiNoCrop(const ap_uint<l1ct::PFRegion::BITWIDTH> input_reg, const ap_uint<l1ct::HadCaloObj::BITWIDTH_SLIM> input_calo[NCALO], ap_uint<LINPUPPI_DATA_SIZE_FWD> output[LINPUPPI_NCHANN_FWD_OUTNC]);

void linpuppi_set_debug(bool debug);


//=================================================
// HF Region Parameters
#if defined(REG_HF)

#define LINPUPPI_DR2MAX  4727 // 0.3 cone
#define LINPUPPI_DR2MIN   525 // 0.1 cone
#define LINPUPPI_dzCut     40 // unused

#define LINPUPPI_iptMax    400 // 100.0/LINPUPPI_ptLSB

#define LINPUPPI_ptSlopeNe  0.25
#define LINPUPPI_ptSlopePh  0.25
#define LINPUPPI_ptZeroNe   14.
#define LINPUPPI_ptZeroPh   14.
#define LINPUPPI_alphaSlope 0.6
#define LINPUPPI_alphaZero  9.0
#define LINPUPPI_alphaCrop  4.0
#define LINPUPPI_priorNe    6.0
#define LINPUPPI_priorPh    6.0

#define LINPUPPI_ptCut     10.0

#endif

#endif
