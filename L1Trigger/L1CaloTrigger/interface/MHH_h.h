#ifndef L1Trigger_L1CaloTrigger_MHH_h
#define L1Trigger_L1CaloTrigger_MHH_h

#include <ap_int.h>

namespace p2mhh {

static constexpr int MHH_N_INPUT_LINKS = 16;
static constexpr int MHH_N_OUTPUT_LINKS = 6;

typedef ap_uint<48> xvar;

class MHHObject
{
public:
    ap_uint<12> ET;
    ap_uint<10> Eta;
    ap_uint<9> Phi;
    ap_uint<4> PtClusterSeed;
    ap_uint<1> isBarrel;

    MHHObject() { ET = 0; Eta = 0; Phi = 0; PtClusterSeed = 0; isBarrel = 0; }

    void getMHHGammas(ap_uint<48> i)
    {
        ET = i.range(11, 0);
        Eta = i.range(18, 12);
        Phi = i.range(27, 19);
        PtClusterSeed = 0;
        isBarrel = i[47] ? 1 : 0;
    }

    void getMHHJetsTaus(ap_uint<48> i)
    {
        ET = i.range(11, 0);
        Eta = i.range(17, 12);
        Phi = i.range(26, 18);
        PtClusterSeed = i.range(30, 27);
        isBarrel = i[47] ? 1 : 0;
    }

    MHHObject(const MHHObject& rhs) {
        ET = rhs.ET;
        Eta = rhs.Eta;
        Phi = rhs.Phi;
        PtClusterSeed = rhs.PtClusterSeed;
        isBarrel = rhs.isBarrel;
    }

    MHHObject& operator=(const MHHObject& rhs) {
        ET = rhs.ET;
        Eta = rhs.Eta;
        Phi = rhs.Phi;
        PtClusterSeed = rhs.PtClusterSeed;
        isBarrel = rhs.isBarrel;
        return *this;
    }
};

class MHHSums {
public:
    ap_int<12> Ex[3];
    ap_int<12> Ey[3];
    ap_uint<12> Ht;
    ap_uint<12> SumET;
    ap_uint<12> NObj;

    MHHSums() {
        for (int i = 0; i < 3; ++i) {
            Ex[i] = 0;
            Ey[i] = 0;
        }
        Ht = 0;
        SumET = 0;
        NObj = 0;
    }

    MHHSums(const MHHSums& rhs) {
        for (int i = 0; i < 3; ++i) {
            Ex[i] = rhs.Ex[i];
            Ey[i] = rhs.Ey[i];
        }
        Ht = rhs.Ht;
        SumET = rhs.SumET;
        NObj = rhs.NObj;
    }

    MHHSums& operator=(const MHHSums& rhs) {
        for (int i = 0; i < 3; ++i) {
            Ex[i] = rhs.Ex[i];
            Ey[i] = rhs.Ey[i];
        }
        Ht = rhs.Ht;
        SumET = rhs.SumET;
        NObj = rhs.NObj;
        return *this;
    }
};

void mhh_algo_top(
    ap_uint<576> link_in[MHH_N_INPUT_LINKS],
    ap_uint<576> link_out[MHH_N_OUTPUT_LINKS]
);

}  // namespace p2mhh

#endif

