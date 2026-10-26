#ifndef _ALGO_TOPIP1_H_
#define _ALGO_TOPIP1_H_

#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>


#define N_INPUT_LINKS  18 
#define N_OUTPUT_LINKS_CL1  6
#define N_OUTPUT_LINKS_MIX  3

#define TOWERS_ETA 12
#define TOWERS_PHI 72

#define N_HF_REGIONS 24 

#define N_HF_TOWERS_REGION_ETA 6
#define N_HF_TOWERS_REGION_PHI 6
#define N_HF_PFCLUSTERS_REGION 4 

#define N_SECTORS_PF 6

#define STOWERS_ETA_R 4
#define STOWERS_PHI_R 8
#define STOWERS_ETA 4
#define STOWERS_PHI 24
#define N_STOWERS 96
#define N_HF_STOWERS_REGION 4

#define TEN 10
#define FIVE 5 

#define N_JETS 6
#define N_TAUS 6


using namespace std;

typedef ap_uint<10> loop;

class hftower{
    public:
    ap_uint<8> energy;
    ap_uint<2> fb;
    ap_uint<5> eta;
    ap_uint<7> phi;

    hftower(){
        energy = 0;
        fb = 0;
	eta = 0;
	phi = 0;
    }

    ap_uint<10> gettower(void){
    	ap_uint<10> data;
        data  = 
	((ap_uint<10>)energy & 0xFF) |
	((ap_uint<10>)fb<<8 & 0x300) ;
	return data ;
    }

    void fillhftower(ap_uint<10> i){
    	this->energy = i.range(7, 0);
    	this->fb = i.range(9, 8);
    }

    hftower(ap_uint<10> i){
    	this->energy = i.range(7, 0);
    	this->fb = i.range(9, 8);
    }

    void hftower12(ap_uint<10> i){
	ap_uint<8> en = i.range(7, 0);
    	this->energy = en >> 1 - en >> 2;
    	this->fb = i.range(9, 8);
    }

    void hftower13(ap_uint<10> i){
	ap_uint<8> en = i.range(7, 0);
    	this->energy = en - en >> 1 - en >> 2;
    	this->fb = i.range(9, 8);
    }

    void hftower14(ap_uint<10> i){
	ap_uint<8> en = i.range(7, 0);
    	this->energy = en >> 2;
    	this->fb = i.range(9, 8);
    }

    hftower(const hftower& rhs){
    energy=rhs.energy;
    fb=rhs.fb;
    eta=rhs.eta;
    phi=rhs.phi;
    }

    hftower& operator=(const hftower& rhs){
    this->energy=rhs.energy;
    this->fb=rhs.fb;
    this->eta=rhs.eta;
    this->phi=rhs.phi;
    return *this;
    }

    ap_uint<8> Energy(void) {return energy;}
    ap_uint<2> Fb(void) {return fb;}
};


class wedge{
    public:
    hftower hftowers[48] ;

    void fillwedge(ap_uint<220> i) {

	for(loop j=0; j<10; j++){
	ap_uint<10> value = (i>>(10*j))&0x3FF;
	ap_uint<8> energy1  = (value&0xFF)>>1;
	ap_uint<8> energytmp  = value&0xFF;
	ap_uint<8> energy2  = energytmp - energy1;
	ap_uint<2> fb1  = (value&0x300)>>8 ;
        this->hftowers[j].energy = energy1 ;
        this->hftowers[j].fb = fb1 ;
	this->hftowers[12+j].energy = energy2 ;
	this->hftowers[12+j].fb = fb1 ;
	value = (i>>(10*j+110))&0x3FF;
	energy1  = (value&0xFF)>>1;
	energytmp  = value&0xFF;
	energy2  = energytmp - energy1;
	fb1  = (value&300)>>8 ;
        this->hftowers[24+j].energy = energy1 ;
        this->hftowers[24+j].fb = fb1 ;
	this->hftowers[36+j].energy = energy2 ;
	this->hftowers[36+j].fb = fb1 ;
	}

	ap_uint<10> value = (i>>100)&0x3FF;
	ap_uint<8> energy1  = (value&0xFF)>>2;
	ap_uint<8> energytmp  = value&0xFF;
	ap_uint<8> energy2  = (energytmp - 2*energy1)>>1;
	ap_uint<8> energy3  = (energytmp - energy2-2*energy1);
	ap_uint<2> fb1  = (value&300)>>8 ;
        this->hftowers[10].energy = energy1 ;
        this->hftowers[10].fb = fb1 ;
        this->hftowers[22].energy = energy2 ;
        this->hftowers[22].fb = fb1 ;
        this->hftowers[34].energy = energy3 ;
        this->hftowers[34].fb = fb1 ;
	this->hftowers[46].energy = energy1 ;
	this->hftowers[46].fb = fb1 ;

	value = (i>>210)&0x3FF;
	energy1  = (value&0xFF)>>2;
	energytmp  = value&0xFF;
	energy2  = (energytmp - 2*energy1)>>1;
	energy3  = (energytmp - energy2-2*energy1);
	fb1  = (value&300)>>8 ;
        this->hftowers[11].energy = energy1 ;
        this->hftowers[11].fb = fb1 ;
        this->hftowers[23].energy = energy2 ;
        this->hftowers[23].fb = fb1 ;
        this->hftowers[35].energy = energy3 ;
        this->hftowers[35].fb = fb1 ;
	this->hftowers[47].energy = energy1 ;
	this->hftowers[47].fb = fb1 ;

    }
};

class towermask{
    public:
//    ap_uint<10> energy;
    ap_uint<1> energy;
    ap_uint<5> eta;
    ap_uint<4> phi;

    towermask(){
        energy = 0;
        eta = 31;
	phi = 15;
    }

    towermask(const towermask& rhs){
//    cout << " copy " << endl ;
    energy=rhs.energy;
    eta=rhs.eta;
    phi=rhs.phi;
    }

    towermask& operator=(const towermask& rhs){
//    cout << " assign " << endl ;
    this->energy=rhs.energy;
    this->eta=rhs.eta;
    this->phi=rhs.phi;
        return *this;
    }
};


class tau{
    public:
    ap_uint<12> energy;
    ap_int<8> eta;
    ap_int<7> phi;
    ap_uint<4> flags;
    ap_uint<12> ECAL;
    ap_uint<12> seedet ;
    ap_uint<9> spare;
    ap_uint<64> data;

    tau(){
    energy=0;
    eta=0;
    phi=0;
    flags=0;
    ECAL=0;
    seedet=0;
    spare=0;
    data=0;
    }

    tau& operator=(const tau& rhs){
        energy = rhs.energy;
        eta = rhs.eta;
        phi = rhs.phi;
        flags = rhs.flags;
        ECAL = rhs.ECAL;
        seedet = rhs.seedet;
        spare = rhs.spare;
        data = rhs.data;
        return *this;
    }

    ap_uint<64> getTau(void){
      data = (energy) | 
      (((ap_uint<64>) eta)  << 12) | 
      (((ap_uint<64>) phi)  << 20) | 
      (((ap_uint<64>) flags)       << 27) |
      (((ap_uint<64>) ECAL)       << 31) |
      (((ap_uint<64>) seedet)       << 43) |
      (((ap_uint<64>) spare)       << 55);
    return data ;
    }

    void fillTau(ap_uint<64> i){
    	this->energy = i.range(11, 0);
    	this->eta = i.range(19, 12);
    	this->phi = i.range(26, 20);
    	this->flags = i.range(30, 27);
    	this->ECAL = i.range(42, 31);
    	this->seedet = i.range(54, 43);
    	this->spare = i.range(63, 55);
    }

    ap_uint<12> Energy(void) {return energy;}
    ap_uint<8> Eta(void) {return eta;}
    ap_uint<7> Phi(void) {return phi;}
    ap_uint<64> Data(void) {return data;}

};


class sums{
    public:
    ap_uint<12> energy;
    ap_int<8> eta;
    ap_int<7> phi;
    ap_uint<4> flags;
    ap_uint<12> ECAL;
    ap_uint<12> seedet ;
    ap_uint<9> spare;
    ap_uint<64> data;

    sums(){
    energy=0;
    eta=0;
    phi=0;
    flags=0;
    ECAL=0;
    seedet=0;
    spare=0;
    data=0;
    }

    sums& operator=(const sums& rhs){
        energy = rhs.energy;
        eta = rhs.eta;
        phi = rhs.phi;
        flags = rhs.flags;
        ECAL = rhs.ECAL;
        seedet = rhs.seedet;
        spare = rhs.spare;
        data = rhs.data;
        return *this;
    }

    ap_uint<64> getSum(void){
      data = (energy) | 
      (((ap_uint<64>) eta)  << 12) | 
      (((ap_uint<64>) phi)  << 20) | 
      (((ap_uint<64>) flags)       << 27) |
      (((ap_uint<64>) ECAL)       << 31) |
      (((ap_uint<64>) seedet)       << 43) |
      (((ap_uint<64>) spare)       << 55);
    return data ;
    }

    void fillsums(ap_uint<64> i){
    	this->energy = i.range(11, 0);
    	this->eta = i.range(19, 12);
    	this->phi = i.range(26, 20);
    	this->flags = i.range(30, 27);
    	this->ECAL = i.range(42, 31);
    	this->seedet = i.range(54, 43);
    	this->spare = i.range(63, 55);
    }

    ap_uint<12> Energy(void) {return energy;}
    ap_uint<8> Eta(void) {return eta;}
    ap_uint<7> Phi(void) {return phi;}
    ap_uint<64> Data(void) {return data;}

};

class jet{
    public:
    ap_uint<12> energy;
    ap_int<8> eta;
    ap_int<7> phi;
    ap_uint<4> flags;
    ap_uint<12> ECAL;
    ap_uint<12> seedet ;
    ap_uint<9> spare;
    ap_uint<64> data;

    jet(){
    energy=0;
    eta=0;
    phi=0;
    flags=0;
    ECAL=0;
    seedet=0;
    spare=0;
    data=0;
    }

    jet& operator=(const jet& rhs){
        energy = rhs.energy;
        eta = rhs.eta;
        phi = rhs.phi;
        flags = rhs.flags;
        ECAL = rhs.ECAL;
        seedet = rhs.seedet;
        spare = rhs.spare;
        data = rhs.data;
        return *this;
    }

    ap_uint<64> getJet(void){
      data = (energy) | 
      (((ap_uint<64>) eta)  << 12) | 
      (((ap_uint<64>) phi)  << 20) | 
      (((ap_uint<64>) flags)       << 27) |
      (((ap_uint<64>) ECAL)       << 31) |
      (((ap_uint<64>) seedet)       << 43) |
      (((ap_uint<64>) spare)       << 55);
    return data ;
    }

    void fillJet(ap_uint<64> i){
    	this->energy = i.range(11, 0);
    	this->eta = i.range(19, 12);
    	this->phi = i.range(26, 20);
    	this->flags = i.range(30, 27);
    	this->ECAL = i.range(42, 31);
    	this->seedet = i.range(54, 43);
    	this->spare = i.range(63, 55);
    }

    ap_uint<12> Energy(void) {return energy;}
    ap_uint<8> Eta(void) {return eta;}
    ap_uint<7> Phi(void) {return phi;}
    ap_uint<64> Data(void) {return data;}

};

class stower{
    public:
    ap_uint<12> energy;
    ap_uint<12> SeedEnergy;
    ap_uint<4> eta ;
    ap_uint<5> phi ;
    ap_uint<15> flags;
    ap_uint<48> data;

    stower(){
        energy = 0;
        SeedEnergy = 0;
     	eta = 15 ;
   	phi = 31 ;
        flags = 0;
        data = 0;
    }

    stower(const stower& rhs){
        energy = rhs.energy;
        SeedEnergy = rhs.SeedEnergy;
        flags = rhs.flags;
        data = rhs.data;
	eta = rhs.eta ;
	phi = rhs.phi ;
    }

    stower& operator=(const stower& rhs){
        energy = rhs.energy;
        SeedEnergy = rhs.SeedEnergy;
        flags = rhs.flags;
        data = rhs.data;
	eta = rhs.eta ;
	phi = rhs.phi ;
        return *this;
    }

    ap_uint<48> getstower(void){
      data = (((ap_uint<48>) energy)  | 
      (((ap_uint<48>) SeedEnergy) << 12) |
      (((ap_uint<48>) eta) << 24) |
      (((ap_uint<48>) phi) << 28) |
      (((ap_uint<48>) flags)  << 33)) ;
    return data ;
    }

    void fillstower(ap_uint<48> i){
    	this->energy = i.range(11, 0);
    	this->SeedEnergy = i.range(23, 12);
    	this->eta = i.range(27, 24);
    	this->phi = i.range(32, 28);
    	this->flags = i.range(47, 33);
    }

};

class PFcluster{
    public:
    ap_uint<12> energy;
    ap_int<8> eta;
    ap_uint<7> phi;
    ap_uint<4> hoe;
    ap_uint<12> ECAL;
    ap_uint<12> HCAL;
    ap_uint<9> spare;
    ap_uint<64> data;

    PFcluster(){
    energy=0;
    eta=0;
    phi=0;
    hoe=0;
    ECAL=0;
    HCAL=0;
    spare=0;
    data=0;
    }

    PFcluster& operator=(const PFcluster& rhs){
        energy = rhs.energy;
        eta = rhs.eta;
        phi = rhs.phi;
        hoe = rhs.hoe;
        ECAL = rhs.ECAL;
        HCAL = rhs.HCAL;
        spare = rhs.spare;
        data = rhs.data;
        return *this;
    }

    ap_uint<64> getPFcluster(void){
      data = (energy) | 
      (((ap_uint<64>) eta)  << 12) | 
      (((ap_uint<64>) phi)  << 20) | 
      (((ap_uint<64>) hoe)       << 27) |
      (((ap_uint<64>) ECAL)       << 31) |
      (((ap_uint<64>) HCAL)       << 43) |
      (((ap_uint<64>) spare)       << 55);
    return data ;
    }

    void fillPFcluster(ap_uint<64> i){
    	this->energy = i.range(11, 0);
    	this->eta = i.range(19, 12);
    	this->phi = i.range(26, 20);
    	this->hoe = i.range(30, 27);
    	this->ECAL = i.range(42, 31);
    	this->HCAL = i.range(54, 43);
    	this->spare = i.range(63, 55);
    }

    ap_uint<12> Energy(void) {return energy;}
    ap_uint<8> Eta(void) {return eta;}
    ap_uint<7> Phi(void) {return phi;}
    ap_uint<64> Data(void) {return data;}
	

};

class hfregion{
public:
hftower hftowers[N_HF_TOWERS_REGION_ETA+4][N_HF_TOWERS_REGION_PHI+4];
//stower stowers[2][2];
//PFcluster pfclusters[N_HF_PFCLUSTERS_REGION] ;
};

class stregion{
public:
stower stowers[2][2];
//PFcluster pfclusters[N_HF_PFCLUSTERS_REGION] ;
};

class pfregion{
public:
//stower stowers[2][2];
PFcluster pfclusters[N_HF_PFCLUSTERS_REGION] ;
};


void getPFClusters(hftower towers[N_HF_TOWERS_REGION_ETA+4][N_HF_TOWERS_REGION_PHI+4], PFcluster pfclusters[N_HF_PFCLUSTERS_REGION]) ;
void createSums(stregion region[N_HF_REGIONS] , sums& sum) ;

void createJets(stower stowers[STOWERS_ETA_R][STOWERS_PHI_R+4], jet jets[5]) ;
void createTaus(stower stowers[STOWERS_ETA_R][STOWERS_PHI_R+4], tau taus[5]) ;

void algo_topIP1(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS_CL1+N_OUTPUT_LINKS_MIX]);

#endif

