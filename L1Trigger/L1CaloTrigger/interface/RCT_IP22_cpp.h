//------------------------------------
// More IP22 Logic for Phase2L1CaloL1RCTEmulator.cc
// (based heavily on IP22.cpp in RCT firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_RCT_IP22_cpp
#define L1Trigger_L1CaloTrigger_RCT_IP22_cpp

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

// RCT IP22 header files
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP22_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bubl_sorter_h.h"

//////////////////////////////////////////////////////////////////////////

// Define algo_top for RCT IP22

namespace p2rctIP22 {

inline void processOutLinks(
   ecaltower ECALTowersSLR0[TOWERS_IN_ETA2 * TOWERS_IN_PHI],
   ecaltower ECALTowersSLR1[TOWERS_IN_ETA5 * TOWERS_IN_PHI],
   ecaltower ECALTowersSLR2[TOWERS_IN_ETA5 * TOWERS_IN_PHI],
   ecaltower ECALTowersSLR3[TOWERS_IN_ETA5 * TOWERS_IN_PHI],
   ap_uint < 576 > link_out[N_OUTPUT_LINKS]) {

   ap_uint < 10 > start;
   ap_uint < 10 > end;
//packing towers for IP3
   for (loop i = 0; i < TOWERS_IN_ETA5 * TOWERS_IN_PHI; i++) {
      start = i * 18;
      end = start + 17;
      link_out[1].range(end, start) = ECALTowersSLR3[i].getecaltower();
      link_out[2].range(end, start) = ECALTowersSLR2[i].getecaltower();
      link_out[3].range(end, start) = ECALTowersSLR1[i].getecaltower();
   }

   for (loop i = 0; i < TOWERS_IN_ETA2 * TOWERS_IN_PHI; i++) {
      start = i * 18;
      end = start + 17;
      link_out[4].range(end, start) = ECALTowersSLR0[i].getecaltower();
   }

}

inline void processInputLinks(ap_uint < 576 > link_in[N_INPUT_LINKS], ecaltower ECALTowersSLR3[TOWERS_IN_ETA5 * TOWERS_IN_PHI], ecaltower ECALTowersSLR2[TOWERS_IN_ETA5 * TOWERS_IN_PHI], ecaltower ECALTowersSLR1[TOWERS_IN_ETA5 * TOWERS_IN_PHI], ecaltower ECALTowersSLR0[TOWERS_IN_ETA2 * TOWERS_IN_PHI], ecalcluster ECALClustersRejected[27]) {

   for (loop i = 0; i < TOWERS_IN_ETA5 * TOWERS_IN_PHI; i++) {
      ap_uint < 10 > start = i * 18;
      ap_uint < 10 > end = start + 17;
      ECALTowersSLR3[i].fillecaltower(link_in[4].range(end, start));
   }
   for (loop i = 0; i < TOWERS_IN_ETA5 * TOWERS_IN_PHI; i++) {
      ap_uint < 10 > start = i * 18;
      ap_uint < 10 > end = start + 17;
      ECALTowersSLR2[i].fillecaltower(link_in[5].range(end, start));
   }
   for (loop i = 0; i < TOWERS_IN_ETA5 * TOWERS_IN_PHI; i++) {
      ap_uint < 10 > start = i * 18;
      ap_uint < 10 > end = start + 17;
      ECALTowersSLR1[i].fillecaltower(link_in[6].range(end, start));
   }

   for (loop i = 0; i < TOWERS_IN_ETA2 * TOWERS_IN_PHI; i++) {
      ap_uint < 10 > start = i * 18;
      ap_uint < 10 > end = start + 17;
      ECALTowersSLR0[i].fillecaltower(link_in[7].range(end, start));
   }

//unpack rejected clusters received from IP21
   for (loop i = 0; i < 9; i++) {
      ap_uint < 10 > start = i * 64;
      ap_uint < 10 > end = start + 63;
      ECALClustersRejected[i].fillecalcluster(link_in[1].range(end, start));
      ECALClustersRejected[i+9].fillecalcluster(link_in[2].range(end, start));
      ECALClustersRejected[i + 18].fillecalcluster(link_in[3].range(end, start));
   }

}
inline void putbackintower2x6(ecalcluster ECALCluster, ecaltower ECALTowers[TOWERS_IN_ETA2 * TOWERS_IN_PHI]) {
   loop id = (ECALCluster.eta) * 6 + (ECALCluster.phi);

   for (loop i = 0; i < TOWERS_IN_ETA2 * TOWERS_IN_PHI; i++) {
      ap_uint < 12 > sum = ECALTowers[i].energy + ECALCluster.energy;
      if (id == i) {
         ECALTowers[i].energy = sum;
      }
   }
}

inline void putbackintower(ecalcluster ECALCluster, ecaltower ECALTowers[TOWERS_IN_ETA5 * TOWERS_IN_PHI]) {

   loop id = (ECALCluster.eta) * 6 + (ECALCluster.phi);

   for (loop i = 0; i < TOWERS_IN_ETA5 * TOWERS_IN_PHI; i++) {
      ap_uint < 12 > sum = ECALTowers[i].energy + ECALCluster.energy;
      if (id == i) {
         ECALTowers[i].energy = sum;
      }
   }

}
inline void dummy_pass_through(ap_uint<576> in, ap_uint<576> &out) {
out = in;
}
inline void algo_top(ap_uint < 576 > link_in[N_INPUT_LINKS], ap_uint < 576 > link_out[N_OUTPUT_LINKS]) {

   ecalcluster ECALClustersRejected[27];
   ecaltower ECALTowersSLR3[TOWERS_IN_ETA5 * TOWERS_IN_PHI];
   ecaltower ECALTowersSLR2[TOWERS_IN_ETA5 * TOWERS_IN_PHI];
   ecaltower ECALTowersSLR1[TOWERS_IN_ETA5 * TOWERS_IN_PHI];
   ecaltower ECALTowersSLR0[TOWERS_IN_ETA2 * TOWERS_IN_PHI];

   //reading clusters and towers

   processInputLinks(link_in, ECALTowersSLR3, ECALTowersSLR2, ECALTowersSLR1, ECALTowersSLR0, ECALClustersRejected);

   ecaltower ECALTowersSLR3C[TOWERS_IN_ETA5 * TOWERS_IN_PHI];
   ecaltower ECALTowersSLR2C[TOWERS_IN_ETA5 * TOWERS_IN_PHI];
   ecaltower ECALTowersSLR1C[TOWERS_IN_ETA5 * TOWERS_IN_PHI];
   ecaltower ECALTowersSLR0C[TOWERS_IN_ETA2 * TOWERS_IN_PHI];
   for (loop i = 0; i < 21; i++) {

      if (ECALClustersRejected[i].spare == 0) putbackintower2x6(ECALClustersRejected[i], ECALTowersSLR0C);
      if (ECALClustersRejected[i].spare == 1) putbackintower(ECALClustersRejected[i], ECALTowersSLR1C);
      if (ECALClustersRejected[i].spare == 2) putbackintower(ECALClustersRejected[i], ECALTowersSLR2C);
      if (ECALClustersRejected[i].spare == 3) putbackintower(ECALClustersRejected[i], ECALTowersSLR3C);
   }

   for (loop i = 0; i < TOWERS_IN_ETA5 * TOWERS_IN_PHI; i++) {
      ap_uint < 12 > energy3 = ECALTowersSLR3[i].energy + ECALTowersSLR3C[i].energy;
      ECALTowersSLR3[i].energy = energy3;
      ap_uint < 12 > energy2 = ECALTowersSLR2[i].energy + ECALTowersSLR2C[i].energy;
      ECALTowersSLR2[i].energy = energy2;
      ap_uint < 12 > energy1 = ECALTowersSLR1[i].energy + ECALTowersSLR1C[i].energy;
      ECALTowersSLR1[i].energy = energy1;
   }

   for (loop i = 0; i < TOWERS_IN_ETA2 * TOWERS_IN_PHI; i++) {
      ap_uint < 12 > energy = ECALTowersSLR0[i].energy + ECALTowersSLR0C[i].energy;
      ECALTowersSLR0[i].energy = energy;
   }

   /*---------------------------------link 0------------------------------------*/

   link_out[0] = 0;
   link_out[1] = 0;
   link_out[2] = 0;
   link_out[3] = 0;
   link_out[4] = 0;

   processOutLinks(ECALTowersSLR0, ECALTowersSLR1, ECALTowersSLR2, ECALTowersSLR3, link_out);
   dummy_pass_through(link_in[0], link_out[0]);

}

} // namespace p2rctIP22

#endif
