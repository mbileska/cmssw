//------------------------------------
// IP1 5x6 Logic for Phase2L1CaloL1RCTEmulator.cc
// (based heavily on algo_top.h in RCT firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_RCT_IP1_5x6_h
#define L1Trigger_L1CaloTrigger_RCT_IP1_5x6_h

#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

namespace p2rctIP1_5x6 {

static constexpr int N_INPUT_LINKS = 30;
static constexpr int N_OUTPUT_LINKS = 2;

static constexpr int CRYSTALS_IN_PHI = 30;
static constexpr int CRYSTALS_IN_ETA = 25;
static constexpr int CRYSTALS_IN_ETA2 = 7;
static constexpr int CRYSTALS_IN_ETA23 = 11;

static constexpr int TOWERS_IN_ETA = 5;
static constexpr int TOWERS_IN_PHI = 6;
static constexpr int N_CLUSTERS = 3;

static constexpr int SEED_THRESHOLD = 2;
//static constexpr int SEED_THRESHOLD = 1023; // dont do any clustering

using namespace std;

typedef ap_uint<6> loop;

class ecalcrystalmask{
    public:
//    ap_uint<10> energy;
    ap_uint<1> energy;

    ecalcrystalmask(){
        energy = 0;
    }

    ecalcrystalmask(const ecalcrystalmask& rhs){
//    cout << " copy " << endl ;
    energy=rhs.energy;
    }

    ecalcrystalmask& operator=(const ecalcrystalmask& rhs){
//    cout << " assign " << endl ;
    this->energy=rhs.energy;
        return *this;
    }
};

class ecalcrystal{
    public:
    ap_uint<10> energy;
    ap_uint<5> timing;
    ap_uint<1> spike;
    ap_uint<5> eta;
    ap_uint<5> phi;

    ecalcrystal(){
        energy = 0;
        timing = 0;
        spike = 0;
        eta = 0;
        phi = 0;
    }

    inline ap_uint<16> getecalcrystal(void){
    	ap_uint<16> data;
        data  = 
	((ap_uint<16>)energy & 0x3FF) |
	(((ap_uint<16>)timing)<<10 & 0x1F) |
	(((ap_uint<16>)spike)<<15 & 0x1) ;
	return data ;
    }

    ecalcrystal(ap_uint<16> i){
    	this->energy = i.range(9, 0);
    	this->timing = i.range(14, 10);
    	this->spike = i.range(15, 15);
    }

    inline void setecalcrystal(ap_uint<16> i, ap_uint<5> k, ap_uint<5> j){
    	this->energy = i.range(9, 0);
    	this->timing = i.range(14, 10);
    	this->spike = i.range(15, 15);
    	this->eta = k;
    	this->phi = j;
    }

    ecalcrystal(const ecalcrystal& rhs){
//    cout << " copy " << endl ;
    energy=rhs.energy;
    timing=rhs.timing;
    spike=rhs.spike;
    eta=rhs.eta;
    phi=rhs.phi;
    }

    ecalcrystal& operator=(const ecalcrystal& rhs){
//    cout << " assign " << endl ;
    this->energy=rhs.energy;
    this->timing=rhs.timing;
    this->spike=rhs.spike;
    this->eta=rhs.eta;
    this->phi=rhs.phi;
        return *this;
    }

    inline ap_uint<10> Energy(void) {return energy;}
    inline ap_uint<5> Eta(void) {return eta;}
    inline ap_uint<5> Phi(void) {return phi;}
};


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

    inline ap_uint<64> getecalcluster(void){
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

    inline void fillecalcluster(ap_uint<64> i){
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
    inline ap_uint<12> Energy(void) {return energy;}
    inline ap_uint<5> Eta(void) {return eta;}
    inline ap_uint<5> Phi(void) {return phi;}
//    et5x5() {return et5x5;}
//    et2x5() {return et2x5;}
//    timing() {return timing;}
//    spike() {return spike;}
//    satur() {return satur;}
//    brems() {return brems;}
    inline ap_uint<64> Data(void) {return data;}
	
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

    inline ap_uint<18> getecaltower(void){
      data = (((ap_uint<18>) energy)  | 
      (((ap_uint<18>) timing) << 12) |
      (((ap_uint<18>) spike)  << 17)) ;
    return data ;
    }

};

void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);

} // namespace p2rctIP1_5x6

#endif
