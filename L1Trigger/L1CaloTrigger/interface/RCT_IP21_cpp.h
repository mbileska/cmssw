//------------------------------------
// More IP21 Logic for Phase2L1CaloL1RCTEmulator.cc
// (based heavily on algo_top.cpp in RCT firmware repo)
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_RCT_IP21_cpp
#define L1Trigger_L1CaloTrigger_RCT_IP21_cpp

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

// RCT IP21 header files and data formats
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP21_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bubl_sorter_h.h"

//////////////////////////////////////////////////////////////////////////

// Define algo_top for RCT IP21

namespace p2rctIP21 {

inline void processOutLinks(ecalcluster ECALClustersSorted[N_CLUSTERS],
   ap_uint < 576 > link_out[N_OUTPUT_LINKS], ap_uint < 7 > SS_fun[25]) {

   ap_uint < 10 > start;
   ap_uint < 10 > end;

   rctecalcluster RCTECALClusters[N_CLUSTERS_OUT];

//threshhold calibration for top 9 clusters
   for (loop i = 0; i < N_CLUSTERS_OUT; i++) {
      ap_uint < 7 > slr = ECALClustersSorted[i + 23].spare;
      //ap_uint < 7 > etatmp = 3 - ECALClustersSorted[i + 23].eta;
      ap_uint < 7 > etatmp = ECALClustersSorted[i + 23].eta;

      ap_uint < 7 >  rctEta = etatmp + ((ap_uint < 7 > ) slr << 4) + ((ap_uint < 7 > ) slr << 3) + slr;
      //rctEta = etatmp + (slr == 0 ? 0 : (slr == 1 ? 9 : (slr == 2 ? 24 : 59)));



      RCTECALClusters[i].fillrctecalcluster2(ECALClustersSorted[i + 23], rctEta, SS_fun);
   }


//packing top 9 clusters to send to IP3
   for (loop i = 0; i < N_CLUSTERS_OUT; i++) {
      start = i * 64;
      end = start + 63;
      //cout<<RCTECALClusters[i].energy<<endl; //<<-------------------------------------------------------------
      link_out[0].range(end, start) = RCTECALClusters[i].getrctecalcluster();
   }

loop m=0;
//packing rejected clusters to send to IP22
   for (loop i = 2; i < 11; i++) {
      start = m * 64;
      end = start + 63;
      link_out[1].range(end, start) = ECALClustersSorted[i].getecalcluster();
      link_out[2].range(end, start) = ECALClustersSorted[i+9].getecalcluster();
      link_out[3].range(end, start) = ECALClustersSorted[i+18].getecalcluster();
      m++;
   }




}

inline void processInputLinks(ap_uint < 576 > link_in[N_INPUT_LINKS], ecalcluster ECALClustersSLR3[N_CLUSTERS_5x6], ecalcluster ECALClustersSLR2[N_CLUSTERS_5x6], ecalcluster ECALClustersSLR1[N_CLUSTERS_5x6], ecalcluster ECALClustersSLR0[N_CLUSTERS_2x6]) {

   for (loop i = 0; i < N_CLUSTERS_5x6; i++) {
      ap_uint < 10 > start = i * 64;
      ap_uint < 10 > end = start + 63;
      ECALClustersSLR3[i].fillecalcluster(link_in[0].range(end, start));
   }
   for (loop i = 0; i < N_CLUSTERS_5x6; i++) {
      ap_uint < 10 > start = i * 64;
      ap_uint < 10 > end = start + 63;
      ECALClustersSLR2[i].fillecalcluster(link_in[1].range(end, start));
   }
   for (loop i = 0; i < N_CLUSTERS_5x6; i++) {
      ap_uint < 10 > start = i * 64;
      ap_uint < 10 > end = start + 63;
      ECALClustersSLR1[i].fillecalcluster(link_in[2].range(end, start));
   }

   for (loop i = 0; i < N_CLUSTERS_2x6; i++) {
      ap_uint < 10 > start = i * 64;
      ap_uint < 10 > end = start + 63;
      ECALClustersSLR0[i].fillecalcluster(link_in[3].range(end, start));
   }

}

inline void getSums(ecalcluster low, ecalcluster high, ap_uint < 12 > & EnergySum, ap_uint < 10 > & EnergySum5x5, ap_uint < 10 > & EnergySum2x5) {

   ap_uint < 12 > EnergyH = high.energy;
   ap_uint < 10 > EnergyH5x5 = high.et5x5;
   ap_uint < 10 > EnergyH2x5 = high.et2x5;

   ap_uint < 12 > EnergyL = low.energy;
   ap_uint < 10 > EnergyL5x5 = low.et5x5;
   ap_uint < 10 > EnergyL2x5 = low.et2x5;

   EnergySum = EnergyH + EnergyL;
   EnergySum5x5 = EnergyH5x5 + EnergyL5x5;
   EnergySum2x5 = EnergyH2x5 + EnergyL2x5;

}

inline void stitchClusters2x6(ecalcluster ECALClustersLowEta[N_CLUSTERS_5x6], ecalcluster ECALClustersHighEta[N_CLUSTERS_2x6]) {

   /* low Eta we use last 3 clusters, high Eta first 3, we match those clusters that
   are at eta=24 Low with eta=0 High if the dPhi<=2 */

   for (loop j = 6; j < 9; j++) {
      /* low Eta */
      for (loop i = 0; i < 3; i++) {
         /* high Eta */

         bool etaH = (ECALClustersHighEta[i].eta == 0) ? true : false;
         bool etaL = (ECALClustersLowEta[j].eta == 24) ? true : false;

         ap_uint < 5 > PhiH = ECALClustersHighEta[i].phi;
         ap_uint < 5 > PhiL = ECALClustersLowEta[j].phi;

         ap_uint < 5 > dPhi = (PhiL > PhiH) ? (PhiL - PhiH) : (PhiH - PhiL);
         bool dphi = (dPhi <= 2) ? true : false;

         bool stitch = etaL && etaH && dphi;

         ap_uint < 12 > EnergyH = ECALClustersHighEta[i].energy;
         ap_uint < 12 > EnergyL = ECALClustersLowEta[j].energy;

         ap_uint < 12 > EnergySum;
         ap_uint < 10 > EnergySum5x5;
         ap_uint < 10 > EnergySum2x5;

         getSums(ECALClustersLowEta[j], ECALClustersHighEta[i], EnergySum, EnergySum5x5, EnergySum2x5);

         if (stitch) {
            if (EnergyH > EnergyL) {
               ECALClustersHighEta[i].energy = EnergySum;
               ECALClustersHighEta[i].et5x5 = EnergySum5x5;
               ECALClustersHighEta[i].et2x5 = EnergySum2x5;
               ECALClustersLowEta[j].energy = 0;
            } else {
               ECALClustersHighEta[i].energy = 0;
               ECALClustersLowEta[j].energy = EnergySum;
               ECALClustersLowEta[j].et5x5 = EnergySum5x5;
               ECALClustersLowEta[j].et2x5 = EnergySum2x5;
            }
         }

      }
   }

}

inline void stitchClusters(ecalcluster ECALClustersLowEta[N_CLUSTERS_5x6], ecalcluster ECALClustersHighEta[N_CLUSTERS_5x6]) {

   /* low Eta we use last 3 clusters, high Eta first 3, we match those clusters that
   are at eta=24 Low with eta=0 High if the dPhi<=2 */

   for (loop j = 6; j < 9; j++) {
      /* low Eta */
      for (loop i = 0; i < 3; i++) {
         /* high Eta */

         bool etaH = (ECALClustersHighEta[i].eta == 0) ? true : false;
         bool etaL = (ECALClustersLowEta[j].eta == 24) ? true : false;

         ap_uint < 5 > PhiH = ECALClustersHighEta[i].phi;
         ap_uint < 5 > PhiL = ECALClustersLowEta[j].phi;

         ap_uint < 5 > dPhi = (PhiL > PhiH) ? (PhiL - PhiH) : (PhiH - PhiL);
         bool dphi = (dPhi <= 2) ? true : false;

         bool stitch = etaL && etaH && dphi;

         ap_uint < 12 > EnergyH = ECALClustersHighEta[i].energy;
         ap_uint < 12 > EnergyL = ECALClustersLowEta[j].energy;

         ap_uint < 12 > EnergySum;
         ap_uint < 10 > EnergySum5x5;
         ap_uint < 10 > EnergySum2x5;

         getSums(ECALClustersLowEta[j], ECALClustersHighEta[i], EnergySum, EnergySum5x5, EnergySum2x5);

         if (stitch) {
            if (EnergyH > EnergyL) {
               ECALClustersHighEta[i].energy = EnergySum;
               ECALClustersHighEta[i].et5x5 = EnergySum5x5;
               ECALClustersHighEta[i].et2x5 = EnergySum2x5;
               ECALClustersLowEta[j].energy = 0;
            } else {
               ECALClustersHighEta[i].energy = 0;
               ECALClustersLowEta[j].energy = EnergySum;
               ECALClustersLowEta[j].et5x5 = EnergySum5x5;
               ECALClustersLowEta[j].et2x5 = EnergySum2x5;
            }
         }

      }
   }

}


inline void algo_top(ap_uint < 576 > link_in[N_INPUT_LINKS], ap_uint < 576 > link_out[N_OUTPUT_LINKS],
   ap_uint < 7 > ss_fun0, ap_uint < 7 > ss_fun1, ap_uint < 7 > ss_fun2, ap_uint < 7 > ss_fun3,
   ap_uint < 7 > ss_fun4, ap_uint < 7 > ss_fun5, ap_uint < 7 > ss_fun6, ap_uint < 7 > ss_fun7,
   ap_uint < 7 > ss_fun8, ap_uint < 7 > ss_fun9, ap_uint < 7 > ss_fun10, ap_uint < 7 > ss_fun11,
   ap_uint < 7 > ss_fun12, ap_uint < 7 > ss_fun13, ap_uint < 7 > ss_fun14, ap_uint < 7 > ss_fun15,
   ap_uint < 7 > ss_fun16, ap_uint < 7 > ss_fun17, ap_uint < 7 > ss_fun18, ap_uint < 7 > ss_fun19,
   ap_uint < 7 > ss_fun20, ap_uint < 7 > ss_fun21, ap_uint < 7 > ss_fun22, ap_uint < 7 > ss_fun23,
   ap_uint < 7 > ss_fun24) {

   ecalcluster ECALClusters[N_CLUSTERS];
   ecalcluster ECALClustersSorted[N_CLUSTERS];
   ecalcluster ECALClustersSLR3[N_CLUSTERS_5x6];
   ecalcluster ECALClustersSLR2[N_CLUSTERS_5x6];
   ecalcluster ECALClustersSLR1[N_CLUSTERS_5x6];
   ecalcluster ECALClustersSLR0[N_CLUSTERS_2x6];
   ap_uint < 7 > SS_fun[25];

   // Pack the individual ss_fun variables into the array
   SS_fun[0] = ss_fun0;
   SS_fun[1] = ss_fun1;
   SS_fun[2] = ss_fun2;
   SS_fun[3] = ss_fun3;
   SS_fun[4] = ss_fun4;
   SS_fun[5] = ss_fun5;
   SS_fun[6] = ss_fun6;
   SS_fun[7] = ss_fun7;
   SS_fun[8] = ss_fun8;
   SS_fun[9] = ss_fun9;
   SS_fun[10] = ss_fun10;
   SS_fun[11] = ss_fun11;
   SS_fun[12] = ss_fun12;
   SS_fun[13] = ss_fun13;
   SS_fun[14] = ss_fun14;
   SS_fun[15] = ss_fun15;
   SS_fun[16] = ss_fun16;
   SS_fun[17] = ss_fun17;
   SS_fun[18] = ss_fun18;
   SS_fun[19] = ss_fun19;
   SS_fun[20] = ss_fun20;
   SS_fun[21] = ss_fun21;
   SS_fun[22] = ss_fun22;
   SS_fun[23] = ss_fun23;
   SS_fun[24] = ss_fun24;

   //reading clusters

   processInputLinks(link_in, ECALClustersSLR3, ECALClustersSLR2, ECALClustersSLR1, ECALClustersSLR0);

   stitchClusters(ECALClustersSLR3, ECALClustersSLR2);
   stitchClusters2x6(ECALClustersSLR1, ECALClustersSLR0);
   stitchClusters(ECALClustersSLR2, ECALClustersSLR1);

   for (loop i = 0; i < N_CLUSTERS_5x6; i++) {
      ECALClusters[i] = ECALClustersSLR3[i];
      ECALClusters[i].spare = 3;
   }
   for (loop i = 0; i < N_CLUSTERS_5x6; i++) {
      ECALClusters[N_CLUSTERS_5x6 + i] = ECALClustersSLR2[i];
      ECALClusters[N_CLUSTERS_5x6 + i].spare = 2;
   }
   for (loop i = 0; i < N_CLUSTERS_5x6; i++) {
      ECALClusters[2 * N_CLUSTERS_5x6 + i] = ECALClustersSLR1[i];
      ECALClusters[2 * N_CLUSTERS_5x6 + i].spare = 1;
   }
   for (loop i = 0; i < N_CLUSTERS_2x6; i++) {
      ECALClusters[3 * N_CLUSTERS_5x6 + i] = ECALClustersSLR0[i];
      ECALClusters[3 * N_CLUSTERS_5x6 + i].spare = 0;
   }

   //bitonicSort32(ECALClusters, ECALClustersSorted);

   bubl_sorter(ECALClusters, ECALClustersSorted);

   /*---------------------------------link 0------------------------------------*/

   link_out[0] = 0;
   link_out[1] = 0;
   link_out[2] = 0;
   link_out[3] = 0;
   link_out[4] = 0;
   link_out[5] = 0;
   link_out[6] = 0;
   link_out[7] = 0;

   processOutLinks(ECALClustersSorted, link_out, SS_fun);

   link_out[4] = link_in[4];
   link_out[5] = link_in[5];
   link_out[6] = link_in[6];
   link_out[7] = link_in[7];

   // Set last bit of links 4-7 equal to the first bit of link 0
   ap_uint<1> first_bit = link_out[0].range(0, 0);
   link_out[4].range(575, 575) = first_bit;
   link_out[5].range(575, 575) = first_bit;
   link_out[6].range(575, 575) = first_bit;
   link_out[7].range(575, 575) = first_bit;


}

} // namespace p2rctIP21

#endif
