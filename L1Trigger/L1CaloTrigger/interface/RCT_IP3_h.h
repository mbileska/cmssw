//------------------------------------
// IP3 Logic for Phase2L1CaloL1RCTEmulator.cc
// (based heavily on algo_top.h in RCT firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_RCT_IP3_h
#define L1Trigger_L1CaloTrigger_RCT_IP3_h

#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

namespace p2rctIP3 {

static constexpr int N_INPUT_LINKS_ECAL = 5;
static constexpr int N_INPUT_LINKS_HCAL = 4;
static constexpr int N_INPUT_LINKS = 9;
static constexpr int N_OUTPUT_LINKS = 4;

static constexpr int N_TOWERS_IN_ETA5 = 5;
static constexpr int N_TOWERS_IN_ETA2 = 2;
static constexpr int N_TOWERS_IN_ETA = 17;
static constexpr int N_TOWERS_IN_PHI = 6;

static constexpr int N_CLUSTERS = 9;

static constexpr int N_TOWERS_HCAL_IN_PHI = 4;
static constexpr int N_TOWERS_HCAL_IN_ETA = 16;

static constexpr int N_TOWERS_IN_HCAL_REGION = 32;

using namespace std;

typedef ap_uint<8> loop;

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

 rctecalcluster(ap_uint<10> seedEnergy, ap_uint<12> energy, ap_uint<7> eta, ap_uint<5> phi, ap_uint<10> et5x5, ap_uint<7> wps, ap_uint<5> timing, ap_uint<1> spike, ap_uint<1> satur, ap_uint<2> brems, ap_uint<3> spare){
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
    ap_uint<7> Eta(void) {return eta;}
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


class hcaltower{
    public:
    ap_uint<10> energy;
    ap_uint<6> fb;
    ap_uint<16> data;

    hcaltower(){
        energy = 0;
        fb = 0;
        data = 0;
    }

    hcaltower& operator=(const hcaltower& rhs){
        energy = rhs.energy;
        fb = rhs.fb;
        data = rhs.data;
        return *this;
    }

    ap_uint<16> gethcaltower(void){
      data = ((ap_uint<16>) energy)  | 
      (((ap_uint<16>) fb) << 10) ;
    return data ;
    }

    void fillhcaltower(ap_uint<16> i){
    	this->energy = i.range(9, 0);
    	this->fb = i.range(15, 10);
    }

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

class rcttower{
    public:
    ap_uint<10> energy;
    ap_uint<4> hoe;
    ap_uint<2> flags;
    ap_uint<16> data;

    rcttower(){
        energy = 0;
        hoe = 0;
        flags = 0;
        data = 0;
    }

    rcttower& operator=(const rcttower& rhs){
        energy = rhs.energy;
        hoe = rhs.hoe;
        flags = rhs.flags;
        data = rhs.data;
        return *this;
    }

    ap_uint<16> getrcttower(void){
      data = (((ap_uint<16>) energy)  | 
      (((ap_uint<16>) hoe) << 10) |
      (((ap_uint<16>) flags)  << 14)) ;
    return data ;
    }

    void fillrcttower(ap_uint<16> i){
    	this->energy = i.range(9, 0);
    	this->hoe = i.range(13, 10);
    	this->flags = i.range(15, 14);
    }

    void uploadECAL(ecaltower ecaltower){
        this->energy = ((ap_uint<10>)ecaltower.energy);
        if(ecaltower.energy >= 1023) this->energy = 1023 ;
        if(ecaltower.timing >= 5) this->flags = 1 ;
    }

    ap_uint<4> gethoe(ap_uint<10> ECAL, ap_uint<10> HCAL){

	ap_uint<4> hoeOut ;
	ap_uint<1> hoeLSB = 0 ;
	ap_uint<3> hoe = 0 ;
	ap_uint<10> A;
	ap_uint<10> B;

	A = (ECAL > HCAL)?ECAL:HCAL;
	B = (ECAL > HCAL)?HCAL:ECAL;
 
	//if( ECAL == 0 || HCAL == 0 || HCAL >= ECAL) hoeLSB = 0 ;
    if(HCAL >= ECAL) hoeLSB = 0 ;
    else hoeLSB = 1 ;
//	 if( A > B ){
        for(loop i=0; i<7; i++){
	if (A > (B<<(7-i))) { hoe = 7-i ; break ;}
	else  ;
	}
//	 }

	hoeOut = hoeLSB | (hoe << 1) ;
	return hoeOut ;
   }

    void uploadHCAL(hcaltower &hcaltower){
	this->hoe = gethoe(this->energy, hcaltower.energy) ;
	ap_uint<11> tmp = hcaltower.energy + this->energy;
    this->energy = (tmp > 1023) ? (ap_uint<10>) 0x3FF : (ap_uint<10>) tmp;
	ap_uint<2> tmp3 = this->flags ;
	if(hcaltower.fb != 0) this->flags = tmp3 + 2 ;
    }
};

void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);


} // namespace p2rctIP3

#endif
