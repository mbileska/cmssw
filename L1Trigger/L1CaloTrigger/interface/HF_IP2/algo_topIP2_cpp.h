#include "algo_topIP2_h.h"
#include "common/firmware/linpuppi.h"
#include "common/layer1_objs.h"
#include "common/puppi.h"
#include "common/pf.h"

// Create PFRegion for each HF region
l1ct::PFRegion createPFRegion(int r) {
    l1ct::PFRegion region;
    region.hwEtaCenter    = l1ct::glbeta_t(6);
    region.hwPhiCenter    = l1ct::glbphi_t(12 * r + 5);
    region.hwEtaHalfWidth = l1ct::eta_t(7);
    region.hwPhiHalfWidth = l1ct::phi_t(6);
    region.hwEtaExtra     = l1ct::eta_t(0);
    region.hwPhiExtra     = l1ct::phi_t(2);
    return region;
}

// Unpack PFCluster into HadCaloObj
inline l1ct::HadCaloObj unpackToHadCalo(const ap_uint<64>& word) {
    l1ct::HadCaloObj h;
    h.hwPt     = word.range(11,  0);
    h.hwEta    = word.range(19, 12);
    h.hwPhi    = word.range(26, 20);
    h.hwHoe    = word.range(30, 27);
    h.hwEmPt   = word.range(42, 31);
    h.hwEmID   = 0;
    h.hwSrrTot = 0;
    h.hwMeanZ  = 0;
    return h;
}

inline ap_uint<64> packPuppiObj(const l1ct::PuppiObj& p) {
    return p.pack();
}

void puppi_sector(
    const l1ct::PFRegion   &region,
    const l1ct::HadCaloObj  caloin[NCALO],
    l1ct::PuppiObj          puppiOut[NNEUTRALS]
) {
    fwdlinpuppi(region, caloin, puppiOut);
}

void algo_topIP2(
    ap_uint<576> link_in[N_HF_REGIONS],
    ap_uint<576> link_out[N_HF_REGIONS]
) {
    for (int r = 0; r < N_HF_REGIONS; r++) {
        l1ct::HadCaloObj hadCalo[NCALO];
        for (int i=0;i<NCALO;i++) {
            hadCalo[i] = unpackToHadCalo(link_in[r].range((i+1)*64-1, i*64));
        }
//        hadCalo[0] = unpackToHadCalo(link_in[r].range( 63,   0));
//        hadCalo[1] = unpackToHadCalo(link_in[r].range(127,  64));
//        hadCalo[2] = unpackToHadCalo(link_in[r].range(191, 128));
//        hadCalo[3] = unpackToHadCalo(link_in[r].range(255, 192));
//        hadCalo[4] = unpackToHadCalo(link_in[r].range(319, 256));
//        hadCalo[5] = unpackToHadCalo(link_in[r].range(383, 320));
//        hadCalo[6] = unpackToHadCalo(link_in[r].range(447, 384));
//        hadCalo[7] = unpackToHadCalo(link_in[r].range(511, 448));

        // Create region config
        l1ct::PFRegion region = createPFRegion(r);

        // Run PUPPI
        l1ct::PuppiObj puppiOut[NNEUTRALS];
        fwdlinpuppi(region, hadCalo, puppiOut);

        // Pack 8 PuppiObj into bits[511:0], spare at bits[575:512] = 0
        link_out[r] = 0;
        link_out[r].range( 63,   0) = packPuppiObj(puppiOut[0]);
        link_out[r].range(127,  64) = packPuppiObj(puppiOut[1]);
        link_out[r].range(191, 128) = packPuppiObj(puppiOut[2]);
        link_out[r].range(255, 192) = packPuppiObj(puppiOut[3]);
        link_out[r].range(319, 256) = packPuppiObj(puppiOut[4]);
        link_out[r].range(383, 320) = packPuppiObj(puppiOut[5]);
        link_out[r].range(447, 384) = packPuppiObj(puppiOut[6]);
        link_out[r].range(511, 448) = packPuppiObj(puppiOut[7]);

#ifndef __SYNTHESIS__
        bool hasInput = false;
        for (int i = 0; i < NCALO; i++) {
            if (hadCalo[i].hwPt > 0) { hasInput = true; break; }
        }
        if (hasInput) {
            std::cout << "=== Region " << r << " ===" << std::endl;
            std::cout << "  etaCenter=" << region.hwEtaCenter
                      << " phiCenter=" << region.hwPhiCenter << std::endl;
            for (int i = 0; i < NCALO; i++) {
                if (hadCalo[i].hwPt > 0)
                    std::cout << "  In[" << i << "]: Pt=" << hadCalo[i].hwPt
                              << " Eta=" << hadCalo[i].hwEta
                              << " Phi=" << hadCalo[i].hwPhi << std::endl;
            }
            for (int i = 0; i < NNEUTRALS; i++) {
                if (puppiOut[i].hwPt > 0)
                    std::cout << "  Out[" << i << "]: Pt=" << puppiOut[i].hwPt
                              << " Eta=" << puppiOut[i].hwEta
                              << " Phi=" << puppiOut[i].hwPhi << std::endl;
            }
        }
#endif
    }
}

// Array interface variant - for testing and validation
void puppi_top_arrays(
    l1ct::PFRegion   regions[N_HF_REGIONS],
    l1ct::HadCaloObj caloin[N_HF_REGIONS][NCALO],
    l1ct::PuppiObj   puppiOut[N_HF_REGIONS][NNEUTRALS]
) {
    for (int r = 0; r < N_HF_REGIONS; r++) {
        fwdlinpuppi(regions[r], caloin[r], puppiOut[r]);
    }
}
