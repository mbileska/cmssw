//------------------------------------
// IP21 Logic for Phase2L1CaloL1RCTEmulator.cc
// (based heavily on algo_top.h in RCT firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_RCT_IP21_h
#define L1Trigger_L1CaloTrigger_RCT_IP21_h

#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

namespace p2rctIP21 {

static constexpr int N_INPUT_LINKS = 8;
static constexpr int N_OUTPUT_LINKS = 8;

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
        this->spare = i.range(63, 61);
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

class rctecalcluster{
    public:
    ap_uint<10> seedEnergy;
    ap_uint<12> energy;
    ap_uint<7> eta;
    ap_uint<5> phi;
    ap_uint<10> et5x5;
    ap_uint<7> wps;
    ap_uint<5> timing;
    ap_uint<1> spike;
    ap_uint<1> satur;
    ap_uint<2> brems;
    ap_uint<3> spare;
    ap_uint<64> data;

    rctecalcluster(){
        seedEnergy = 0;
        energy = 0;
        eta = 0;
        phi = 0;
        et5x5 = 0;
        wps = 0;
        timing = 0;
        spike = 0;
        satur = 0;
        brems = 0;
        spare = 0;
        data = 0;
    }

    rctecalcluster& operator=(const rctecalcluster& rhs){
        seedEnergy = rhs.seedEnergy;
        energy = rhs.energy;
        eta = rhs.eta;
        phi = rhs.phi;
        et5x5 = rhs.et5x5;
        wps = rhs.wps;
        timing = rhs.timing;
        spike = rhs.spike;
        satur = rhs.satur;
        brems = rhs.brems;
        spare = rhs.spare;
        data = rhs.data;
        return *this;
    }

    ap_uint<64> getrctecalcluster(void){
      data = (seedEnergy) | 
      (((ap_uint<64>) energy)  << 10) | 
      (((ap_uint<64>) eta)  << 22) | 
      (((ap_uint<64>) phi)  << 29) | 
      (((ap_uint<64>) et5x5)       << 34) |
      (((ap_uint<64>) wps)       << 44) |
      (((ap_uint<64>) timing)       << 51) |
      (((ap_uint<64>) spike)       << 56) |
      (((ap_uint<64>) satur)       << 57) |
      (((ap_uint<64>) brems)       << 58) |
      (((ap_uint<64>) spare)       << 60);
    return data ;
    }

    void fillrctecalcluster(ap_uint<64> i){
    	this->seedEnergy = i.range(9, 0);
    	this->energy = i.range(21, 10);
    	this->eta = i.range(28, 22);
    	this->phi = i.range(33, 29);
    	this->et5x5 = i.range(43, 34);
    	this->wps = i.range(50, 44);
    	this->timing = i.range(55, 51);
    	this->spike = i.range(56, 56);
    	this->satur = i.range(57, 57);
    	this->brems = i.range(59, 58);
        this->spare = i.range(62,60);
    }
   
    void fillrctecalcluster2(ecalcluster i, ap_uint<7> eta,ap_uint<7> SS_fun[25]){
   
	ap_uint<7> WPs=1 ;
	//ap_uint<7> SS_fun[25] = {126, 125, 125, 124, 123, 123, 122, 122, 121, 121, 121, 121, 121, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120} ;

	for(loop j=0; j<25; j++){
	bool one = ( i.energy >= ((loop)j<<2)) ? true : false ;
	bool two = ( i.energy < ((loop)j<<2)+4) ? true : false ;
	ap_uint<18> left = ((ap_uint<18>)i.et2x5<<7) ;
	ap_uint<18> right = 0  ;

		for(loop k=0; k<7; k++){
		ap_uint<18> tmp = 0  ;
        	if ((SS_fun[j]&(1<<k)) != 0) tmp = ((ap_uint<18>)i.et5x5<<k);
                right += tmp ;}

	bool three = ( left < right ) ? true : false ;
	if ( one && two && three)  WPs=0;
	}
	//WPs=1;

    	this->seedEnergy = i.seedEnergy;
    	this->energy = i.energy;
    	this->eta = eta;
    	this->phi = i.phi;
    	this->et5x5 = i.et5x5;
    	this->wps = WPs ;
    	this->timing = i.timing;
    	this->spike = i.spike;
    	this->satur = i.satur;
    	this->brems = i.brems;
        this->spare = i.spare;
    }


 rctecalcluster(ap_uint<10> seedEnergy, ap_uint<12> energy, ap_uint<7> eta, ap_uint<5> phi, ap_uint<10> et5x5, ap_uint<7> wps, ap_uint<5> timing, ap_uint<1> spike, ap_uint<1> satur, ap_uint<2> brems, ap_uint<4> spare){
        data = (seedEnergy) | 
      (((ap_uint<64>) energy)  << 10) | 
      (((ap_uint<64>) eta)  << 22) | 
      (((ap_uint<64>) phi)  << 29) | 
      (((ap_uint<64>) et5x5)       << 34) |
      (((ap_uint<64>) wps)       << 44) |
      (((ap_uint<64>) timing)       << 51) |
      (((ap_uint<64>) spike)       << 56) |
      (((ap_uint<64>) satur)       << 57) |
      (((ap_uint<64>) brems)       << 58) |
      (((ap_uint<64>) spare)       << 60);
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



//void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);
void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS], ap_uint<7> ss_fun0, ap_uint<7> ss_fun1, ap_uint<7> ss_fun2, ap_uint<7> ss_fun3,
        ap_uint<7> ss_fun4, ap_uint<7> ss_fun5, ap_uint<7> ss_fun6, ap_uint<7> ss_fun7,
        ap_uint<7> ss_fun8, ap_uint<7> ss_fun9, ap_uint<7> ss_fun10, ap_uint<7> ss_fun11,
        ap_uint<7> ss_fun12, ap_uint<7> ss_fun13, ap_uint<7> ss_fun14, ap_uint<7> ss_fun15,
        ap_uint<7> ss_fun16, ap_uint<7> ss_fun17, ap_uint<7> ss_fun18, ap_uint<7> ss_fun19,
        ap_uint<7> ss_fun20, ap_uint<7> ss_fun21, ap_uint<7> ss_fun22, ap_uint<7> ss_fun23,
        ap_uint<7> ss_fun24);


} // namespace p2rctIP21

#endif
