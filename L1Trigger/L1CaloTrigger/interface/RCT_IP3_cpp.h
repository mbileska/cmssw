//------------------------------------
// More IP3 Logic for Phase2L1CaloL1RCTEmulator.cc
// (based heavily on algo_top.cpp in RCT firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_RCT_IP3_cpp
#define L1Trigger_L1CaloTrigger_RCT_IP3_cpp

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
// #include <cstdint>
#include <cstdlib>  // for rand
#include <iostream>
#include <fstream>
#include <memory>

#include <bitset>

// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

// ECAL TPs
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"

// HCAL TPs
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"

// Output tower collection
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloTower.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedClusterCorrelator.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedTowerCorrelator.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/DigitizedClusterGT.h"

#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "L1Trigger/L1CaloTrigger/interface/ParametricCalibration.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1CaloTrigger/interface/Phase2L1CaloEGammaUtils.h"
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1RCT.h"
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1GCT.h"

// RCT CMSSW - firmware interface
#include "L1Trigger/L1CaloTrigger/interface/RCT_IO.h"

// RCT IP3 header files
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP3_h.h"

///////////////////////////////////////////////////////////////////////////

// Define algo_top for RCT IP3

namespace p2rctIP3 {

inline void processOutLinks(rctecalcluster RCTECALClusters[N_CLUSTERS],
                     rcttower RCTTowersPHI[N_TOWERS_IN_PHI][N_TOWERS_IN_ETA],
                     ap_uint<576> link_out[N_OUTPUT_LINKS]){

ap_uint<10> start ;
ap_uint<10> end ;
ap_uint<10> start1 ;
ap_uint<10> end1 ;

for(loop i=0 ; i<N_CLUSTERS; i++){
	start=i*64 ; end=start+63;
	link_out[0].range(end, start) = RCTECALClusters[i].getrctecalcluster() ;
}

for(loop i=0; i<N_TOWERS_IN_ETA; i++){
	start=i*16 ; end=start+15;
	start1=start+272 ; end1=start1+15;
	link_out[1].range(end, start) = RCTTowersPHI[0][i].getrcttower() ;
	link_out[1].range(end1, start1) = RCTTowersPHI[1][i].getrcttower() ;
	link_out[2].range(end, start) = RCTTowersPHI[2][i].getrcttower() ;
	link_out[2].range(end1, start1) = RCTTowersPHI[3][i].getrcttower() ;
	link_out[3].range(end, start) = RCTTowersPHI[4][i].getrcttower() ;
	link_out[3].range(end1, start1) = RCTTowersPHI[5][i].getrcttower() ;
	 //cout << hex <<link_out[3]<<endl;
	 //cout<<"link out 3 ";
}

}

inline void getECAL(ap_uint<576> link_in[N_INPUT_LINKS],ecaltower ECALTowersSLR3[N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI],ecaltower ECALTowersSLR2[N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI],ecaltower ECALTowersSLR1[N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI],ecaltower ECALTowersSLR0[N_TOWERS_IN_ETA2*N_TOWERS_IN_PHI]){

for(loop i=0; i<N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI; i++){
    ap_uint<10> start   = i*18;
    ap_uint<10> end = start + 17;
    ECALTowersSLR3[i].fillecaltower(link_in[1].range(end, start));
    //cout << hex << link_in[1]<<endl;
    //cout <<"tower link 1"<<endl;


  }
for(loop i=0; i<N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI; i++){
    ap_uint<10> start   = i*18;
    ap_uint<10> end = start + 17;
    ECALTowersSLR2[i].fillecaltower(link_in[2].range(end, start));
   // cout << hex << link_in[2]<<endl;
    //cout <<"tower link 2"<<endl;

  }
for(loop i=0; i<N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI; i++){
    ap_uint<10> start   = i*18;
    ap_uint<10> end = start + 17;
    ECALTowersSLR1[i].fillecaltower(link_in[3].range(end, start));
    //cout << hex << link_in[3]<<endl;
    //cout <<"tower link 3"<<endl;

  }
for(loop i=0; i<N_TOWERS_IN_ETA2*N_TOWERS_IN_PHI; i++){
    ap_uint<10> start   = i*18;
    ap_uint<10> end = start + 17;
    ECALTowersSLR0[i].fillecaltower(link_in[4].range(end, start));
   // cout << hex << link_in[4];
   // cout <<"tower link 4"<<endl;


  }
//cout<<ECALTowersSLR0[0].energy<<endl;
  // cout<<"zero zero en";
}

inline void processInputLinks(ap_uint<576> link_in[N_INPUT_LINKS], rctecalcluster RCTECALClusters[N_CLUSTERS], rcttower RCTTowersPHI[N_TOWERS_IN_PHI][N_TOWERS_IN_ETA], bool HCALsecondhalfstarts){

ecaltower ECALTowersSLR3[N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI];
ecaltower ECALTowersSLR2[N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI];
ecaltower ECALTowersSLR1[N_TOWERS_IN_ETA5*N_TOWERS_IN_PHI];
ecaltower ECALTowersSLR0[N_TOWERS_IN_ETA2*N_TOWERS_IN_PHI];


hcaltower HCALTowersR0[N_TOWERS_IN_HCAL_REGION];
hcaltower HCALTowersR1[N_TOWERS_IN_HCAL_REGION];
hcaltower HCALTowersR2[N_TOWERS_IN_HCAL_REGION];
hcaltower HCALTowersR3[N_TOWERS_IN_HCAL_REGION];

for(loop i=0; i<N_CLUSTERS; i++){
    ap_uint<10> start   = i*64;
    ap_uint<10> end = start + 63;
    RCTECALClusters[i].fillrctecalcluster(link_in[0].range(end, start));

  }



getECAL(link_in,ECALTowersSLR3, ECALTowersSLR2, ECALTowersSLR1, ECALTowersSLR0);

//cout<<"ECALTowersSLR0[0].energy :";
//cout<<ECALTowersSLR0[0].energy<<endl;


//cout<<"eta::";
//cout<<RCTECALClusters[0].eta<<endl;

for(loop i=0; i<N_TOWERS_IN_HCAL_REGION; i++){
    ap_uint<10> start   = i*16;
    ap_uint<10> end = start + 15;
    HCALTowersR0[i].fillhcaltower(link_in[5].range(end, start));
  }
for(loop i=0; i<N_TOWERS_IN_HCAL_REGION; i++){
    ap_uint<10> start   = i*16;
    ap_uint<10> end = start + 15;
    HCALTowersR1[i].fillhcaltower(link_in[6].range(end, start));
  }
for(loop i=0; i<N_TOWERS_IN_HCAL_REGION; i++){
    ap_uint<10> start   = i*16;
    ap_uint<10> end = start + 15;
    HCALTowersR2[i].fillhcaltower(link_in[7].range(end, start));
  }
for(loop i=0; i<N_TOWERS_IN_HCAL_REGION; i++){
    ap_uint<10> start   = i*16;
    ap_uint<10> end = start + 15;
    HCALTowersR3[i].fillhcaltower(link_in[8].range(end, start));
  }

for(loop j=0; j<N_TOWERS_IN_PHI; j++){
for(loop i=0; i<N_TOWERS_IN_ETA5; i++){
RCTTowersPHI[j][i+0].uploadECAL(ECALTowersSLR3[i*6+j]) ;
}}
for(loop j=0; j<N_TOWERS_IN_PHI; j++){
for(loop i=0; i<N_TOWERS_IN_ETA5; i++){
RCTTowersPHI[j][i+5].uploadECAL(ECALTowersSLR2[i*6+j]) ;
}}
for(loop j=0; j<N_TOWERS_IN_PHI; j++){
for(loop i=0; i<N_TOWERS_IN_ETA5; i++){
RCTTowersPHI[j][i+10].uploadECAL(ECALTowersSLR1[i*6+j]) ;
}}

for(loop j=0; j<N_TOWERS_IN_PHI; j++){
for(loop i=0; i<N_TOWERS_IN_ETA2; i++){
RCTTowersPHI[j][i+15].uploadECAL(ECALTowersSLR0[i*6+j]) ;
}}

if(!HCALsecondhalfstarts){
for(loop j=0; j<N_TOWERS_HCAL_IN_PHI; j++){
for(loop i=0; i<N_TOWERS_HCAL_IN_ETA/2; i++){
//    #pragma HLS UNROLL
RCTTowersPHI[j][i].uploadHCAL(HCALTowersR0[i*4+j]) ;
}}

for(loop j=0; j<N_TOWERS_HCAL_IN_PHI; j++){
for(loop i=0; i<N_TOWERS_HCAL_IN_ETA/2; i++){
//    #pragma HLS UNROLL
RCTTowersPHI[j][i+8].uploadHCAL(HCALTowersR1[i*4+j]) ;
}}

for(loop j=0; j<N_TOWERS_HCAL_IN_PHI/2; j++){
for(loop i=0; i<N_TOWERS_HCAL_IN_ETA/2; i++){
//    #pragma HLS UNROLL
RCTTowersPHI[j+4][i].uploadHCAL(HCALTowersR2[i*4+j]) ;
}}
for(loop j=0; j<N_TOWERS_HCAL_IN_PHI/2; j++){
for(loop i=0; i<N_TOWERS_HCAL_IN_ETA/2; i++){
//    #pragma HLS UNROLL
RCTTowersPHI[j+4][i+8].uploadHCAL(HCALTowersR3[i*4+j]) ;
}}

}
else{
for(loop j=0; j<N_TOWERS_HCAL_IN_PHI/2; j++){
for(loop i=0; i<N_TOWERS_HCAL_IN_ETA/2; i++){
//    #pragma HLS UNROLL
RCTTowersPHI[j][i].uploadHCAL(HCALTowersR2[i*4+j+2]);
}}
for(loop j=0; j<N_TOWERS_HCAL_IN_PHI/2; j++){
for(loop i=0; i<N_TOWERS_HCAL_IN_ETA/2; i++){
//    #pragma HLS UNROLL
RCTTowersPHI[j][i+8].uploadHCAL(HCALTowersR3[i*4+j+2]);
}}
for(loop j=0; j<N_TOWERS_HCAL_IN_PHI; j++){
for(loop i=0; i<N_TOWERS_HCAL_IN_ETA/2; i++){
//    #pragma HLS UNROLL
RCTTowersPHI[j+2][i].uploadHCAL(HCALTowersR0[i*4+j]) ;
}}
for(loop j=0; j<N_TOWERS_HCAL_IN_PHI; j++){
for(loop i=0; i<N_TOWERS_HCAL_IN_ETA/2; i++){
//    #pragma HLS UNROLL
RCTTowersPHI[j+2][i+8].uploadHCAL(HCALTowersR1[i*4+j]) ;
}}
}
}

inline void algo_top(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]){
//#pragma HLS latency min=145

rctecalcluster RCTECALClusters[N_CLUSTERS] ;

rcttower RCTTowersPHI[N_TOWERS_IN_PHI][N_TOWERS_IN_ETA];

//reading clusters and towers
// order of HCAL tower

bool HCALsecondhalfstarts = 0 ;

processInputLinks(link_in, RCTECALClusters, RCTTowersPHI, HCALsecondhalfstarts) ;
//out<< RCTTowersPHI[0][0].energy<<endl;
//cout<<"RCTTowersPHI[0][0]";
//cout<< RCTTowersPHI[5][0].energy<<endl;
//cout<<"RCTTowersPHI[5][0]";
/*---------------------------------link 0------------------------------------*/
        
link_out[0] = 0;
link_out[1] = 0;
link_out[2] = 0;
link_out[3] = 0;

processOutLinks(RCTECALClusters, RCTTowersPHI, link_out);

}

// Copy of algo_top, but with HCALsecondhalfstarts = 1 (to be called on certain RCT cards, depending on orientation of HCAL links)
inline void algo_top_HCALsecondhalfstarts(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]){
//#pragma HLS latency min=145

rctecalcluster RCTECALClusters[N_CLUSTERS] ;

rcttower RCTTowersPHI[N_TOWERS_IN_PHI][N_TOWERS_IN_ETA];

//reading clusters and towers
// order of HCAL tower

bool HCALsecondhalfstarts = 1 ;

processInputLinks(link_in, RCTECALClusters, RCTTowersPHI, HCALsecondhalfstarts) ;
//out<< RCTTowersPHI[0][0].energy<<endl;
//cout<<"RCTTowersPHI[0][0]";
//cout<< RCTTowersPHI[5][0].energy<<endl;
//cout<<"RCTTowersPHI[5][0]";
/*---------------------------------link 0------------------------------------*/
        
link_out[0] = 0;
link_out[1] = 0;
link_out[2] = 0;
link_out[3] = 0;

processOutLinks(RCTECALClusters, RCTTowersPHI, link_out);

}

} // namespace p2rctIP3

#endif
