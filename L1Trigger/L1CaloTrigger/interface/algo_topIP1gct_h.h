#ifndef _ALGO_TOPIP2_H_
#define _ALGO_TOPIP2_H_



#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

#define N_RCT_REGIONS 8 // 4+4 RCT cards are connected to 1 SLR

#define N_INPUT_LINKS  32 // 4 links per RCTcard 4*8 
#define N_OUTPUT_LINKS  21 // [linkA(EG)+linkB(PFClusters)]x2[eta+;eta-]x6[TMI]

#define TOWERS_IN_ETA5 5
#define TOWERS_IN_ETA2 2
#define TOWERS_IN_PHI 6

#define N_CLUSTERS 32 
#define N_CLUSTERS_OUT 9 
#define N_CLUSTERS_5x6 9
#define N_CLUSTERS_2x6 3

#define N_RCT_TOWERS_REGION13 34 // 2 links 2x17 towers 
#define N_RCT_TOWERS_REGION 102 // 2 links 2x17 towers 
#define N_RCT_CLUSTERS_REGION 9 // 
#define N_GCT_TOWERS_REGION_ETA 17
#define N_GCT_TOWERS_REGION_PHI 6
#define N_GCT_PFCLUSTERS_REGION 12

#define N_STOWERS_REGION 12

#define TEN 10
#define FIVE 5 

using namespace std;

typedef ap_uint<8> loop;

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

class PFcluster{
    public:
    ap_uint<12> energy;
    ap_int<8> eta;
    ap_int<7> phi;
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

    void getHoE(){

	ap_uint<1> hoeLSB = 0 ;
	ap_uint<3> he = 0 ;
	ap_uint<12> A;
	ap_uint<12> B;

	A = (ECAL > HCAL)?ECAL:HCAL;
	B = (ECAL > HCAL)?HCAL:ECAL;
 
	if( HCAL >= ECAL) hoeLSB = 0 ;
	else hoeLSB = 1 ;
	 /*if( A > B ){
 	  if(A > 2*B) he = 0b001 ;
  	  if(A > 4*B) he = 0b010 ;
  	  if(A > 8*B) he = 0b011 ;
  	  if(A > 16*B) he = 0b100 ;
  	  if(A > 32*B) he = 0b101 ;
  	  if(A > 64*B) he = 0b110 ;
 	  if(A > 128*B) he = 0b111 ;
 	 }*/
	/*if (A > B) {
	  if (A > 128 * B) he = 0b111;
	  else if (A > 64 * B) he = 0b110;
	  else if (A > 32 * B) he = 0b101;
	  else if (A > 16 * B) he = 0b100;
	  else if (A > 8 * B) he = 0b011;
	  else if (A > 4 * B) he = 0b010;
	  else if (A > 2 * B) he = 0b001;
	  else he = 0b000;
	}*/

	if (A > (B << 7)) he = 0b111;
	else if (A > (B << 6)) he = 0b110;
	else if (A > (B << 5)) he = 0b101;
	else if (A > (B << 4)) he = 0b100;
	else if (A > (B << 3)) he = 0b011;
	else if (A > (B << 2)) he = 0b010;
	else if (A > (B << 1)) he = 0b001;
	else he = 0b000;

	hoe = hoeLSB | (he << 1) ;
   }

    ap_uint<12> Energy(void) {return energy;}
    ap_uint<8> Eta(void) {return eta;}
    ap_uint<7> Phi(void) {return phi;}
    ap_uint<64> Data(void) {return data;}
	

};


class EGcluster{
    public:
    ap_uint<12> energy;
    ap_int<7> eta;
    ap_int<7> phi;
    ap_uint<4> hoe;
    ap_uint<2> hoeWP;
    ap_uint<3> iso;
    ap_uint<2> isoWP;
    ap_uint<6> FB;
    ap_uint<5> timing;
    ap_uint<2> shapeWP;
    ap_uint<2> brems;
    ap_uint<12> spare;
    ap_uint<64> data;

    EGcluster(){
    energy=0;
    eta=0;
    phi=0;
    hoe=0;
    hoeWP=0;
    iso=0;
    isoWP=0;
    FB=0;
    timing=0;
    shapeWP=0;
    brems=0;
    spare=0;
    data=0;
    }

    EGcluster& operator=(const EGcluster& rhs){
        energy = rhs.energy;
        eta = rhs.eta;
        phi = rhs.phi;
        hoe = rhs.hoe;
        hoeWP = rhs.hoeWP;
        iso = rhs.iso;
        isoWP = rhs.isoWP;
        FB = rhs.FB;
        timing = rhs.timing;
        shapeWP = rhs.shapeWP;
        brems = rhs.brems;
        spare = rhs.spare;
        data = rhs.data;
        return *this;
    }

    ap_uint<64> getEGcluster(void){
      data = (energy) | 
      (((ap_uint<64>) eta)  << 12) | 
      (((ap_uint<64>) phi)  << 19) | 
      (((ap_uint<64>) hoe)       << 26) |
      (((ap_uint<64>) hoeWP)       << 30) |
      (((ap_uint<64>) iso)       << 32) |
      (((ap_uint<64>) isoWP)       << 35) |
      (((ap_uint<64>) FB)       << 37) |
      (((ap_uint<64>) timing)       << 43) |
      (((ap_uint<64>) shapeWP)       << 48) |
      (((ap_uint<64>) brems)       << 50) |
      (((ap_uint<64>) spare)       << 52);
    return data ;
    }

    void fillEGcluster(ap_uint<64> i){
    	this->energy = i.range(11, 0);
    	this->eta = i.range(18, 12);
    	this->phi = i.range(25, 19);
    	this->hoe = i.range(29, 26);
    	this->hoeWP = i.range(31, 30);
    	this->iso= i.range(34, 32);
    	this->isoWP= i.range(36, 35);
    	this->FB = i.range(42, 37);
    	this->timing = i.range(47, 43);
    	this->shapeWP = i.range(49, 48);
    	this->brems = i.range(51, 50);
    	this->spare = i.range(63, 52);
    }

    ap_uint<12> Energy(void) {return energy;}
    ap_uint<8> Eta(void) {return eta;}
    ap_uint<7> Phi(void) {return phi;}
    ap_uint<64> Data(void) {return data;}
	

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
    ap_uint<4> spare;
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


class SuperTower{
    public:
    ap_uint<12> energy;
    ap_uint<12> EMenergy;
    ap_uint<4> eta ;
    ap_uint<5> phi ;
    ap_uint<15> flags;
    ap_uint<48> data;

    SuperTower(){
        energy = 0;
        EMenergy = 0;
     	eta = 15 ;
   	phi = 31 ;
        flags = 0;
        data = 0;
    }

    SuperTower(const SuperTower& rhs){
        energy = rhs.energy;
        EMenergy = rhs.EMenergy;
        flags = rhs.flags;
        data = rhs.data;
	eta = rhs.eta ;
	phi = rhs.phi ;
    }

    SuperTower& operator=(const SuperTower& rhs){
        energy = rhs.energy;
        EMenergy = rhs.EMenergy;
        flags = rhs.flags;
        data = rhs.data;
	eta = rhs.eta ;
	phi = rhs.phi ;
        return *this;
    }

    ap_uint<48> getST(void){
      data = (((ap_uint<48>) energy)  | 
      (((ap_uint<48>) EMenergy) << 12) |
      (((ap_uint<48>) eta) << 24) |
      (((ap_uint<48>) phi) << 28) |
      (((ap_uint<48>) flags)  << 33)) ;
    return data ;
    }

    void fillST(ap_uint<48> i){
    	this->energy = i.range(11, 0);
    	this->EMenergy = i.range(23, 12);
    	this->eta = i.range(27, 24);
    	this->phi = i.range(32, 28);
    	this->flags = i.range(47, 33);
    }


};

class rcttower{
    public:
    ap_uint<10> energy;
    ap_uint<4> hoe;
    ap_uint<2> flags;
    ap_uint<16> data;
    ap_uint<5> eta ;
    ap_uint<4> phi ;

    rcttower(){
        energy = 0;
        hoe = 0;
        flags = 0;
        data = 0;
     	eta = 31 ;
   	phi = 15 ;
    }

    rcttower(const rcttower& rhs){
        energy = rhs.energy;
        hoe = rhs.hoe;
        flags = rhs.flags;
        data = rhs.data;
	eta = rhs.eta ;
	phi = rhs.phi ;
    }

    rcttower& operator=(const rcttower& rhs){
        energy = rhs.energy;
        hoe = rhs.hoe;
        flags = rhs.flags;
        data = rhs.data;
	eta = rhs.eta ;
	phi = rhs.phi ;
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

    void getEH(ap_uint<10> &ECAL, ap_uint<10> &HCAL){

	ap_uint<10> A = energy>>(hoe>>1) ;
	ap_uint<10> B = energy - A;

	ECAL = ((hoe&0x1) == 1)?B:A;
	HCAL = ((hoe&0x1) == 1)?A:B;
 
   }

};

class fulltower{
    public:
    ap_uint<12> energy;
    ap_uint<4> hoe;
    ap_uint<2> flags;
    ap_uint<12> ECAL;
    ap_uint<12> HCAL;
    ap_uint<12> EGclusterEnergy;
    ap_uint<5> eta ;
    ap_uint<4> phi ;
    ap_uint<1> spare ;
    ap_uint<64> data;

    fulltower(){
        energy = 0;
        hoe = 0;
        flags = 0;
    	ECAL=0;
    	HCAL=0;
    	EGclusterEnergy=0;
	eta = 0;
        phi = 0;
        spare = 0;
        data = 0;

    }

    fulltower& operator=(const fulltower& rhs){
        energy = rhs.energy;
        hoe = rhs.hoe;
        flags = rhs.flags;
        ECAL = rhs.ECAL;
        HCAL = rhs.HCAL;
        EGclusterEnergy = rhs.EGclusterEnergy;
        eta = rhs.eta;
        phi = rhs.phi;
        spare = rhs.spare;
        data = rhs.data;
        return *this;
    }

    fulltower(const fulltower& rhs){
        energy = rhs.energy;
        hoe = rhs.hoe;
        flags = rhs.flags;
        ECAL = rhs.ECAL;
        HCAL = rhs.HCAL;
        EGclusterEnergy = rhs.EGclusterEnergy;
        eta = rhs.eta;
        phi = rhs.phi;
        spare = rhs.spare;
        data = rhs.data;
    }


    ap_uint<64> getfulltower(void){
      data = ((ap_uint<64>) energy)  | 
      (((ap_uint<64>) hoe) << 12) |
      (((ap_uint<64>) flags)  << 16) |
      (((ap_uint<64>) ECAL)  << 18) |
      (((ap_uint<64>) HCAL)  << 30) |
      (((ap_uint<64>) EGclusterEnergy)  << 42) |
      (((ap_uint<64>) eta) << 54) |
      (((ap_uint<64>) phi) << 59) |
      (((ap_uint<64>) spare)  << 63) ;
    return data ;
    }

    void fillfulltower(ap_uint<64> i){
    	this->energy = i.range(11, 0);
    	this->hoe = i.range(15, 12);
    	this->flags = i.range(17, 16);
    	this->ECAL = i.range(29, 18);
    	this->HCAL = i.range(41,30);
    	this->EGclusterEnergy = i.range(53, 42);
    	this->eta = i.range(58, 54);
    	this->phi = i.range(62, 59);
    	this->spare = i.range(63, 63);
    }

    void copyfromrcttower(rcttower rcttower){
	ap_uint<10> A = rcttower.energy>>(rcttower.hoe>>1) ;
	ap_uint<10> B = rcttower.energy - A;

	ECAL = ((hoe&0x1) == 1)?B:A;
	HCAL = ((hoe&0x1) == 1)?A:B;

	energy = rcttower.energy ;
	hoe = 0 ;
	flags = rcttower.flags ;
   }

    void copyfromrctecalcluster(rctecalcluster rctecalcluster, towermask mask){
    	ap_uint<12> A = EGclusterEnergy + rctecalcluster.energy*mask.energy ;
    	ap_uint<12> B = energy + rctecalcluster.energy*mask.energy ;

	energy = B ;
	EGclusterEnergy = A ;
	hoe = 0 ;
   }

    void getHoE(){

	ap_uint<1> hoeLSB = 0 ;
	ap_uint<3> he = 0 ;
	ap_uint<12> A;
	ap_uint<12> B;

	A = (ECAL > HCAL)?ECAL:HCAL;
	B = (ECAL > HCAL)?HCAL:ECAL;
 
	if( HCAL >= ECAL) hoeLSB = 0 ;
	else hoeLSB = 1 ;
	 if( A > B ){
 	  if(A > 2*B) he = 0b001 ;
  	  if(A > 4*B) he = 0b010 ;
  	  if(A > 8*B) he = 0b011 ;
  	  if(A > 16*B) he = 0b100 ;
  	  if(A > 32*B) he = 0b101 ;
  	  if(A > 64*B) he = 0b110 ;
 	  if(A > 128*B) he = 0b111 ;
 	 }
	hoe = hoeLSB | (he << 1) ;
   }
};



class rctregion{

public:
rctecalcluster rctecalclusters[N_RCT_CLUSTERS_REGION] ;
//rcttower rcttowers[N_RCT_TOWERS_REGION] ;
rcttower rcttowers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4] ;

SuperTower supertowers[N_STOWERS_REGION] ;

EGcluster egclusters[N_RCT_CLUSTERS_REGION] ;

PFcluster pfclusters[N_GCT_PFCLUSTERS_REGION] ;

towermask Grid[N_RCT_CLUSTERS_REGION] ;

 void calculateST(){

	for(loop i=0; i<3; i++){

		for(loop k=0; k<3; k++){

			for(loop m=0; m<N_STOWERS_REGION/2; m++){

	    supertowers[2*m].energy += rcttowers[i+2+3*m][k+2].energy ;
	    supertowers[2*m+1].energy += rcttowers[i+2+3*m][k+5].energy ;
	}}} ;
}


 void completeEG(){

	for(loop i=0; i<N_RCT_CLUSTERS_REGION; i++){

	    egclusters[i].energy = rctecalclusters[i].energy ;
	    egclusters[i].eta = rctecalclusters[i].eta ;
	    egclusters[i].phi = rctecalclusters[i].phi ;
	    egclusters[i].brems = rctecalclusters[i].brems ;
	}} ;


 void fillSTlink(ap_uint<576> &link_out){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<N_STOWERS_REGION; i++)
          {    

		start=i*48 ; end=start+47; 
		link_out.range(end, start) = supertowers[i].getST() ; 
	  }
	}


 void fillEGlink(ap_uint<576> &link_out){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<N_RCT_CLUSTERS_REGION; i++)
          {    

		start=i*64 ; end=start+63; 
		link_out.range(end, start) = egclusters[i].getEGcluster() ; 
	  }
	}

 void fillEGlink1n8(ap_uint<576> &link_out1){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<8; i++)
          {    

		start=(i+1)*64 ; end=start+63; 
		link_out1.range(end, start) = egclusters[i].getEGcluster() ; 
	  }
	}

 void fillEGlink8n1(ap_uint<576> &link_out1){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<8; i++)
          {    

		start=i*64 ; end=start+63; 
		link_out1.range(end, start) = egclusters[i].getEGcluster() ; 
	  }
	}

 void fillEGlink17(ap_uint<576> &link_out1, ap_uint<576> &link_out2){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<1; i++)
          {    

		start=(i+8)*64 ; end=start+63; 
		link_out1.range(end, start) = egclusters[i].getEGcluster() ; 
	  }

        for(loop i=1; i<8; i++)
          {    

		start=(i-1)*64 ; end=start+63; 
		link_out2.range(end, start) = egclusters[i].getEGcluster() ; 
	  }
	}

 void fillEGlink26(ap_uint<576> &link_out1, ap_uint<576> &link_out2){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<2; i++)
          {    

		start=(i+7)*64 ; end=start+63; 
		link_out1.range(end, start) = egclusters[i].getEGcluster() ; 
	  }

        for(loop i=2; i<8; i++)
          {    

		start=(i-2)*64 ; end=start+63; 
		link_out2.range(end, start) = egclusters[i].getEGcluster() ; 
	  }
	}

 void fillPFlink93(ap_uint<576> &link_out1, ap_uint<576> &link_out2){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<9; i++)
          {    

		start=i*64 ; end=start+63; 
		link_out1.range(end, start) = pfclusters[i].getPFcluster() ; 
	  }

        for(loop i=9; i<12; i++)
          {    

		start=(i-9)*64 ; end=start+63; 
		link_out2.range(end, start) = pfclusters[i].getPFcluster() ; 
	  }
}

 void fillPFlink66(ap_uint<576> &link_out1, ap_uint<576> &link_out2){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<6; i++)
          {    

		start=(3+i)*64 ; end=start+63; 
		link_out1.range(end, start) = pfclusters[i].getPFcluster() ; 
	  }

        for(loop i=6; i<12; i++)
          {    

		start=(i-6)*64 ; end=start+63; 
		link_out2.range(end, start) = pfclusters[i].getPFcluster() ; 
	  }
}

 void fillPFlink39(ap_uint<576> &link_out1, ap_uint<576> &link_out2){
 
        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop i=0; i<3; i++)
          {    

		start=(i+6)*64 ; end=start+63; 
		link_out1.range(end, start) = pfclusters[i].getPFcluster() ; 
	  }

        for(loop i=0; i<9; i++)
          {    

		start=i*64 ; end=start+63; 
		link_out2.range(end, start) = pfclusters[i].getPFcluster() ; 
	  }
}


} ;

void createPFClusters(rcttower regiongcttowers[22][10], PFcluster pfclusters[N_GCT_PFCLUSTERS_REGION]) ;
void algo_topIP1gct(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);

#endif

