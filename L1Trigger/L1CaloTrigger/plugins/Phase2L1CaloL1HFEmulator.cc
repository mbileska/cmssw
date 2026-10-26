/* 
 * Description: Phase 2 HF Emulator
 */

// system include files
//#include <ap_int.h>
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


// HF CMSSW - firmware interface
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/HF_output.h"


//HF IP1 header files and data formats
#include "L1Trigger/L1CaloTrigger/interface/HF_IP1/algo_topIP1_h.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP1/algo_topIP1_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP1/hfcommon_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP1/bubl_sorter_h.h"

//HF IP2 header files and data formats
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/algo_topIP2_h.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/algo_topIP2_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/bit_encoding.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/datatypes.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/hls_reg.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/layer1_multiplicities.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/layer1_objs.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/pairwise_sum_reduce.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/pf.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/puppi.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/firmware/linpuppi.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/firmware/linpuppi_bits.h"
#include "L1Trigger/L1CaloTrigger/interface/HF_IP2/common/firmware/linpuppi_cpp.h"
////////////////////////////////////////////////////////////////////////////////

// Declare the Phase2L1CaloL1HFEmulator class and its methods 




class Phase2L1CaloL1HFEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1CaloL1HFEmulator(const edm::ParameterSet&);
  ~Phase2L1CaloL1HFEmulator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hfToken_;
};


////////////////////////////////////////////////////////////////////////////////


//Phase2L1CaloL1HFEmulator intializer, destructor, and produce methods

Phase2L1CaloL1HFEmulator::Phase2L1CaloL1HFEmulator(const edm::ParameterSet& iConfig)
    : hfToken_(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("hcalDigis"))) {
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh0");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh1");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh2");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh3");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh4");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh5");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh6");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh7");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1PosEtaCh8");

  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh0");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh1");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh2");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh3");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh4");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh5");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh6");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh7");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP1NegEtaCh8");

  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2PosEtaCh0");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2PosEtaCh1");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2PosEtaCh2");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2PosEtaCh3");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2PosEtaCh4");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2PosEtaCh5");

  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2NegEtaCh0");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2NegEtaCh1");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2NegEtaCh2");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2NegEtaCh3");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2NegEtaCh4");
  produces<l1tp2::hfOutputLinkCollection>("LinkOutIP2NegEtaCh5");
}


void Phase2L1CaloL1HFEmulator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Initialize the 18 input links for the HF emulator
  ap_uint<576> link_in_pos[N_INPUT_LINKS] = {};
  ap_uint<576> link_in_neg[N_INPUT_LINKS] = {};

  edm::Handle<HcalTrigPrimDigiCollection> hfHandle;
  iEvent.getByToken(hfToken_, hfHandle);

  ap_uint<10> tp_pos[TOWERS_ETA][TOWERS_PHI] = {};
  ap_uint<10> tp_neg[TOWERS_ETA][TOWERS_PHI] = {};
  for (const auto& hit : *hfHandle.product()) {
    int ieta = hit.id().ieta();
    int iphi = hit.id().iphi();
    if (abs(ieta) < 30 || abs(ieta) > 41) continue; // only HF
    if (hit.id().version() != 1) continue; // only version 1
    ap_uint<10> val = ap_uint<10>(hit.SOI_compressedEt()) | (ap_uint<10>(hit.SOI_fineGrain(0)) << 8) |
                       (ap_uint<10>(hit.SOI_fineGrain(1)) << 9); // get the 10-bit energy value
    if (ieta > 0)
      tp_pos[ieta - 30][iphi - 1] = val;
    else
      tp_neg[-ieta - 30][iphi - 1] = val;
  }

  // link inde = 3 * sector + chunk, with per-link bitfields [A(0...109), B(110...219)]
  for (int sector = 0; sector < N_INPUT_LINKS / 3; ++sector) {
    for (int i = 0; i < TOWERS_ETA - 2; ++i) {
      const ap_uint<10> startA = i * 10;
      const ap_uint<10> endA = startA + 9;
      const ap_uint<10> startB = startA + 110;
      const ap_uint<10> endB = startB + 9;

      for (int chunk = 0; chunk < 3; ++chunk) {
        const int linkIdx = 3 * sector + chunk;
        const int phiBase = 12 * sector + 4 * chunk;
        link_in_pos[linkIdx].range(endA, startA) = tp_pos[i][phiBase]; // A side
        link_in_pos[linkIdx].range(endB, startB) = tp_pos[i][phiBase + 2]; // B side
        link_in_neg[linkIdx].range(endA, startA) = tp_neg[i][phiBase];
        link_in_neg[linkIdx].range(endB, startB) = tp_neg[i][phiBase + 2];
      }
    }

    for (int chunk = 0; chunk < 3; ++chunk) {
      const int linkIdx = 3 * sector + chunk;
      const int phiBase = 12 * sector + 4 * chunk;
      link_in_pos[linkIdx].range(109, 100) = tp_pos[10][phiBase];  // row 10, A side
      link_in_pos[linkIdx].range(219, 210) = tp_pos[11][phiBase];  // row 11, B side
      link_in_neg[linkIdx].range(109, 100) = tp_neg[10][phiBase];  // row 10, A side
      link_in_neg[linkIdx].range(219, 210) = tp_neg[11][phiBase];  // row 11, B side
    }
  }

  
  /////////////// Run the HF Emulator firmware ///////////////////////////

  std::cout << "Starting the HF Emulator..." << std::endl;

  // Run algo_top (firmware code)
  ap_uint<576> link_out_ip1_pos[N_OUTPUT_LINKS_CL1 + N_OUTPUT_LINKS_MIX];
  ap_uint<576> link_out_ip1_neg[N_OUTPUT_LINKS_CL1 + N_OUTPUT_LINKS_MIX];
  algo_topIP1(link_in_pos, link_out_ip1_pos);
  algo_topIP1(link_in_neg, link_out_ip1_neg);

  ap_uint<576> link_in_ip2_pos[N_HF_REGIONS];
  ap_uint<576> link_out_ip2_pos[N_HF_REGIONS];
  ap_uint<576> link_in_ip2_neg[N_HF_REGIONS];
  ap_uint<576> link_out_ip2_neg[N_HF_REGIONS];

  // Copy the output of IP1 to the input of IP2
  for (int i = 0; i < N_HF_REGIONS; i++) {
    link_in_ip2_pos[i] = link_out_ip1_pos[i];
    link_in_ip2_neg[i] = link_out_ip1_neg[i];
  }
  // Run the IP2 emulator
  algo_topIP2(link_in_ip2_pos, link_out_ip2_pos);
  algo_topIP2(link_in_ip2_neg, link_out_ip2_neg);


  // put IP1 into output collections
  for (int i = 0; i < N_OUTPUT_LINKS_CL1 + N_OUTPUT_LINKS_MIX; ++i) {
    auto col_pos = std::make_unique<l1tp2::hfOutputLinkCollection>();
    auto col_neg = std::make_unique<l1tp2::hfOutputLinkCollection>();
    col_pos->push_back(l1tp2::hfOutputLink(link_out_ip1_pos[i]));
    col_neg->push_back(l1tp2::hfOutputLink(link_out_ip1_neg[i]));
    iEvent.put(std::move(col_pos), "LinkOutIP1PosEtaCh" + std::to_string(i));
    iEvent.put(std::move(col_neg), "LinkOutIP1NegEtaCh" + std::to_string(i));
  }

  // Put IP2 into output collections
  for (int i = 0; i < N_HF_REGIONS; ++i) {
    auto col_pos = std::make_unique<l1tp2::hfOutputLinkCollection>();
    auto col_neg = std::make_unique<l1tp2::hfOutputLinkCollection>();
    col_pos->push_back(l1tp2::hfOutputLink(link_out_ip2_pos[i]));
    col_neg->push_back(l1tp2::hfOutputLink(link_out_ip2_neg[i]));
    iEvent.put(std::move(col_pos), "LinkOutIP2PosEtaCh" + std::to_string(i));
    iEvent.put(std::move(col_neg), "LinkOutIP2NegEtaCh" + std::to_string(i));
  }
}


///////////////////////////////////////////////////////////////////////////////////


void Phase2L1CaloL1HFEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("hcalDigis", edm::InputTag("simHcalTriggerPrimitiveDigis"));
  descriptions.addWithDefaultLabel(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1CaloL1HFEmulator);