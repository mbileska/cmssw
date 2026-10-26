//------------------------------------
// IP22 Logic for Phase2L1CaloL1RCTEmulator.cc
// (based heavily on algo_top.h in RCT firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_RCT_IP22_h
#define L1Trigger_L1CaloTrigger_RCT_IP22_h

#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

namespace p2rctIP22 {

static constexpr int N_INPUT_LINKS = 8;
static constexpr int N_OUTPUT_LINKS = 5;

static constexpr int TOWERS_IN_ETA5 = 5;
static constexpr int TOWERS_IN_ETA2 = 2;
static constexpr int TOWERS_IN_PHI = 6;

static constexpr int N_CLUSTERS = 32; 
static constexpr int N_CLUSTERS_OUT = 9; 
static constexpr int N_CLUSTERS_5x6 = 9;
static constexpr int N_CLUSTERS_2x6 = 3;

using namespace std;

typedef ap_uint<8> loop;

class ecalcluster{
    public:
    ap_uint<10> seedEnergy;
    ap_uint<12> energy;
    ap_uint<5> eta;
    ap_uint<5> phi;
    ap_uint<10> et5x5;
    ap_uint<10> et2x5;
    ap_uint<5> timing;
    ap_uint<1> spike;
    ap_uint<1> satur;
    ap_uint<2> brems;
    ap_uint<3> spare;
    ap_uint<64> data;

    ecalcluster(){
        seedEnergy = 0;
        energy = 0;
        eta = 0;
        phi = 0;
        et5x5 = 0;
        et2x5 = 0;
        timing = 0;
        spike = 0;
        satur = 0;
        brems = 0;
        spare = 0;
        data = 0;
    }

    ecalcluster& operator=(const ecalcluster& rhs){
        seedEnergy = rhs.seedEnergy;
        energy = rhs.energy;
        eta = rhs.eta;
        phi = rhs.phi;
        et5x5 = rhs.et5x5;
        et2x5 = rhs.et2x5;
        timing = rhs.timing;
        spike = rhs.spike;
        satur = rhs.satur;
        brems = rhs.brems;
        spare = rhs.spare;
        data = rhs.data;
        return *this;
    }

    ap_uint<64> getecalcluster(void){
      data = (seedEnergy) | 
      (((ap_uint<64>) energy)  << 10) | 
      (((ap_uint<64>) eta)  << 22) | 
      (((ap_uint<64>) phi)  << 27) | 
      (((ap_uint<64>) et5x5)       << 32) |
      (((ap_uint<64>) et2x5)       << 42) |
      (((ap_uint<64>) timing)       << 52) |
      (((ap_uint<64>) spike)       << 57) |
      (((ap_uint<64>) satur)       << 58) |
      (((ap_uint<64>) brems)       << 59) |
      (((ap_uint<64>) spare)       << 61);
    return data ;
    }

    void fillecalcluster(ap_uint<64> i){
    	this->seedEnergy = i.range(9, 0);
    	this->energy = i.range(21, 10);
    	this->eta = i.range(26, 22);
    	this->phi = i.range(31, 27);
    	this->et5x5 = i.range(41, 32);
    	this->et2x5 = i.range(51, 42);
    	this->timing = i.range(56, 52);
    	this->spike = i.range(57, 57);
    	this->satur = i.range(58, 58);
    	this->brems = i.range(60, 59);
    }

 ecalcluster(ap_uint<10> seedEnergy, ap_uint<12> energy, ap_uint<5> eta, ap_uint<5> phi, ap_uint<10> et5x5, ap_uint<10> et2x5, ap_uint<5> timing, ap_uint<1> spike, ap_uint<1> satur, ap_uint<2> brems, ap_uint<3> spare){
        data = (seedEnergy) | 
      (((ap_uint<64>) energy)  << 10) | 
      (((ap_uint<64>) eta)  << 22) | 
      (((ap_uint<64>) phi)  << 27) | 
      (((ap_uint<64>) et5x5)       << 32) |
      (((ap_uint<64>) et2x5)       << 42) |
      (((ap_uint<64>) timing)       << 52) |
      (((ap_uint<64>) spike)       << 57) |
      (((ap_uint<64>) satur)       << 58) |
      (((ap_uint<64>) brems)       << 59) |
      (((ap_uint<64>) spare)       << 61);
    }

//    seedEnergy() {return seedEnergy;}
    ap_uint<12> Energy(void) {return energy;}
    ap_uint<5> Eta(void) {return eta;}
    ap_uint<5> Phi(void) {return phi;}
//    et5x5() {return et5x5;}
//    et2x5() {return et2x5;}
//    timing() {return timing;}
//    spike() {return spike;}
//    satur() {return satur;}
//    brems() {return brems;}
    ap_uint<64> Data(void) {return data;}
	
//    operator uint64_t() {return (ap_uint<64>) data;}


};


class ecaltower{
    public:
    ap_uint<12> energy;
    ap_uint<5> timing;
    ap_uint<1> spike;
    ap_uint<18> data;

    ecaltower(){
        energy = 0;
        timing = 0;
        spike = 0;
        data = 0;
    }

    ecaltower& operator=(const ecaltower& rhs){
        energy = rhs.energy;
        timing = rhs.timing;
        spike = rhs.spike;
        data = rhs.data;
        return *this;
    }

    ap_uint<18> getecaltower(void){
      data = (((ap_uint<18>) energy)  | 
      (((ap_uint<18>) timing) << 12) |
      (((ap_uint<18>) spike)  << 17)) ;
    return data ;
    }

    void fillecaltower(ap_uint<18> i){
    	this->energy = i.range(11, 0);
    	this->timing = i.range(16, 12);
    	this->spike = i.range(17, 17);
    }

};

//void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);
void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);


} // namespace p2rctIP22

#endif
