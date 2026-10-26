/* 
 * Description: Generates info for the CaloL1 ML on RCT project
 *              Based on Phase2 CaloL1 emulator RCT code
 */

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

// RCT IP1 header files and data formats
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP1_5x6_h.h"
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP1_5x6_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP1_2x6_h.h"
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP1_2x6_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/RCT_IP3_h.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/RCT_output.h"

///////////////////////////////////////////////////////////////////////////////

// Declare the Phase2L1CaloMLonRCTSetupProducer class and its methods

class Phase2L1CaloMLonRCTSetupProducer : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1CaloMLonRCTSetupProducer(const edm::ParameterSet&);
  ~Phase2L1CaloMLonRCTSetupProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPEBToken_;
  edm::EDGetTokenT<edm::SortedCollection<HcalTriggerPrimitiveDigi>> hcalTPToken_;
  edm::ESGetToken<CaloTPGTranscoder, CaloTPGRecord> decoderTag_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryTag_;
  const CaloSubdetectorGeometry* ebGeometry;
  const CaloSubdetectorGeometry* hbGeometry;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> hbTopologyTag_;
  const HcalTopology* hcTopology_;
};

//////////////////////////////////////////////////////////////////

// Phase2L1CaloMLonRCTSetupProducer initializer, destructor, and produce methods

Phase2L1CaloMLonRCTSetupProducer::Phase2L1CaloMLonRCTSetupProducer(const edm::ParameterSet& iConfig)
    : ecalTPEBToken_(consumes<EcalEBTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalTPEB"))),
      hcalTPToken_(consumes<edm::SortedCollection<HcalTriggerPrimitiveDigi>>(iConfig.getParameter<edm::InputTag>("hcalTP"))),
      decoderTag_(esConsumes<CaloTPGTranscoder, CaloTPGRecord>(edm::ESInputTag("", ""))),
      caloGeometryTag_(esConsumes<CaloGeometry, CaloGeometryRecord>(edm::ESInputTag("", ""))),
      hbTopologyTag_(esConsumes<HcalTopology, HcalRecNumberingRecord>(edm::ESInputTag("", ""))) {
  
  produces<l1tp2::rctOutputLinkCollection>("EGammaSLR3");
  produces<l1tp2::rctOutputLinkCollection>("EGammaSLR2");
  produces<l1tp2::rctOutputLinkCollection>("EGammaSLR1");
  produces<l1tp2::rctOutputLinkCollection>("EGammaSLR0");
  produces<l1tp2::rctOutputLinkCollection>("ECALUnclusteredSLR3");
  produces<l1tp2::rctOutputLinkCollection>("ECALUnclusteredSLR2");
  produces<l1tp2::rctOutputLinkCollection>("ECALUnclusteredSLR1");
  produces<l1tp2::rctOutputLinkCollection>("ECALUnclusteredSLR0");
  produces<l1tp2::rctOutputLinkCollection>("HCAL8");
  produces<l1tp2::rctOutputLinkCollection>("HCAL7");
  produces<l1tp2::rctOutputLinkCollection>("HCAL6");
  produces<l1tp2::rctOutputLinkCollection>("HCAL5");
}

void Phase2L1CaloMLonRCTSetupProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Output collections
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outEGammaSLR3(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outEGammaSLR2(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outEGammaSLR1(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outEGammaSLR0(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outECALUnclusteredSLR3(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outECALUnclusteredSLR2(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outECALUnclusteredSLR1(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outECALUnclusteredSLR0(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outHCAL8(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outHCAL7(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outHCAL6(make_unique<l1tp2::rctOutputLinkCollection>());
  std::unique_ptr<l1tp2::rctOutputLinkCollection> link_outHCAL5(make_unique<l1tp2::rctOutputLinkCollection>());

  std::cout << "Starting the ML on RCT setup producer..." << std::endl;

  // Detector geometry
  const auto& caloGeometry = iSetup.getData(caloGeometryTag_);
  ebGeometry = caloGeometry.getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  hbGeometry = caloGeometry.getSubdetectorGeometry(DetId::Hcal, HcalBarrel);

  const auto& hbTopology = iSetup.getData(hbTopologyTag_);
  hcTopology_ = &hbTopology;
  HcalTrigTowerGeometry theTrigTowerGeometry(hcTopology_);

  const auto& decoder = iSetup.getData(decoderTag_);

  //***************************************************//
  // Get the ECAL hits
  //***************************************************//
  edm::Handle<EcalEBTrigPrimDigiCollection> pcalohits;
  iEvent.getByToken(ecalTPEBToken_, pcalohits);

  // Initialize array of energies to all 0s
  float crystalEnergies[p2rctIO::N_CARDS][p2rctIO::N_TOWERS_ETA*p2rctIO::CRYSTALS_IN_TOWER_ETA][p2rctIO::N_TOWERS_PHI*p2rctIO::CRYSTALS_IN_TOWER_PHI] = {{{ 0 }}};

  for (const auto& hit : *pcalohits.product()) {
    if (hit.encodedEt() > 0)  // hit.encodedEt() returns an int corresponding to 2x the crystal Et
    {
      // Et is 10 bit, by keeping the ADC saturation Et at 120 GeV it means that you have to multiply by 0.125 (input LSB)
      float et = hit.encodedEt() * 0.125;
      if (et < 0.5) {
        continue;  // Reject hits with < 500 MeV ET
      }

      // Get cell coordinates and info
      auto cell = ebGeometry->getGeometry(hit.id());

      p2rctIO::SimpleCaloHit ehit;
      ehit.setPosition(GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z()));
      ehit.setEt(et);

      // Find the card that this hit is in
      for(int cc = 0; cc < p2rctIO::N_CARDS; cc++) {
        if (ehit.isInCard(cc)) {

          // Get the crystal iEta and iPhi, relative to the bottom left corner of the card
          int local_iEta = ehit.crystalLocaliEta(cc);
          int local_iPhi = ehit.crystalLocaliPhi(cc);

          crystalEnergies[cc][local_iEta][local_iPhi] = crystalEnergies[cc][local_iEta][local_iPhi] + ehit.et();

        }
      }
      // break;
    }
  }

  //***************************************************//
  // Get the HCAL hits
  //***************************************************//
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi>> hbhecoll;
  iEvent.getByToken(hcalTPToken_, hbhecoll);

  // Initialize array of energies to all 0s
  float HCALtowerEnergies[p2rctIO::N_REGIONS_ETA][p2rctIO::N_REGIONS_PHI][p2rctIO::TOWERS_IN_REGION_ETA][p2rctIO::TOWERS_IN_REGION_PHI] = {{{{ 0 }}}};
  // Initialize array of feature bits to all 0s as ints, convert back to bitset later
  int HCALtowerFeatureBits[p2rctIO::N_REGIONS_ETA][p2rctIO::N_REGIONS_PHI][p2rctIO::TOWERS_IN_REGION_ETA][p2rctIO::TOWERS_IN_REGION_PHI] = {{{{ 0 }}}};

  for (const auto& hit : *hbhecoll.product()) {
    float et = decoder.hcaletValue(hit.id(), hit.t0());
    // same thing as SOI_compressedEt() in HcalTriggerPrimitiveDigi.h///
    if (et <= 0)
      continue;

    if (!(hcTopology_->validHT(hit.id()))) {
      LogError("Phase2L1CaloL1RCTEmulator")
          << " -- Hcal hit DetID not present in HCAL Geom: " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloL1RCTEmulator");
      continue;
    }
    const std::vector<HcalDetId>& hcId = theTrigTowerGeometry.detIds(hit.id());
    if (hcId.empty()) {
      LogError("Phase2L1CaloL1RCTEmulator") << "Cannot find any HCalDetId corresponding to " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloL1RCTEmulator");
      continue;
    }
    if (hcId[0].subdetId() > 1) {
      continue;
    }
    GlobalVector hcal_tp_position = GlobalVector(0., 0., 0.);
    for (const auto& hcId_i : hcId) {
      if (hcId_i.subdetId() > 1) {
        continue;
      }
      // get the first HCAL TP/ cell
      auto cell = hbGeometry->getGeometry(hcId_i);
      if (cell == nullptr) {
        continue;
      }
      GlobalVector tmpVector = GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());
      hcal_tp_position = tmpVector;

      break;
    }
    p2rctIO::SimpleCaloHit hhit;
    hhit.setPosition(hcal_tp_position);
    hhit.setEt(et);

    // Calculate feature bits as int to store
    int fbits = 0;
    for(int iBit=0; iBit < p2rctIO::N_BITS_HCAL_FBITS; iBit++) {
      fbits = fbits + (int)hit.SOI_fineGrain(iBit) * std::pow(2, iBit);
    }
    hhit.setFeatureBits(fbits);

    int BCP_iEta = hhit.BCPcardiEta();
    int BCP_iPhi = hhit.BCPcardiPhi();
    int iEtaInBCP = hhit.iEtaInBCPcard();
    int iPhiInBCP = hhit.iPhiInBCPcard();

    if(HCALtowerEnergies[BCP_iEta][BCP_iPhi][iEtaInBCP][iPhiInBCP] > 0){
      std::cout << "Adding energy to an HCAL tower that already has non-zero energy!" << std::endl;
    }
    HCALtowerEnergies[BCP_iEta][BCP_iPhi][iEtaInBCP][iPhiInBCP] = HCALtowerEnergies[BCP_iEta][BCP_iPhi][iEtaInBCP][iPhiInBCP] + hhit.et();
    HCALtowerFeatureBits[BCP_iEta][BCP_iPhi][iEtaInBCP][iPhiInBCP] = hhit.featurebits();

  }

  // Iterate through cells of crystalEnergies and fill a vector of p2rctIO::RCTcardECAL objects

  std::vector<p2rctIO::RCTcardECAL> cardsECAL(p2rctIO::N_CARDS);
  for(int cc=0; cc < p2rctIO::N_CARDS; cc++) {
    for(int iEtaCrystalCard=0; iEtaCrystalCard < p2rctIO::N_TOWERS_ETA*p2rctIO::CRYSTALS_IN_TOWER_ETA; iEtaCrystalCard++) {
      for(int iPhiCrystalCard=0; iPhiCrystalCard < p2rctIO::N_TOWERS_PHI*p2rctIO::CRYSTALS_IN_TOWER_PHI; iPhiCrystalCard++) {
        float thisCrystalEnergy = crystalEnergies[cc][iEtaCrystalCard][iPhiCrystalCard];
        cardsECAL[cc].addHit(thisCrystalEnergy, 0.0, 0, iEtaCrystalCard, iPhiCrystalCard);
      }
    }
  }

  // Iterate through towers of HCALtowerEnergies and fill an array of HCAL links via a BCPcardsHCAL object

  p2rctIO::BCPcardsHCAL cardsHCAL;
  for(int BCPiEta=0; BCPiEta < p2rctIO::N_REGIONS_ETA; BCPiEta++) {
    for(int BCPiPhi=0; BCPiPhi < p2rctIO::N_REGIONS_PHI; BCPiPhi++) {
      for(int iEtaInBCP=0; iEtaInBCP < p2rctIO::TOWERS_IN_REGION_ETA; iEtaInBCP++) {
        for(int iPhiInBCP=0; iPhiInBCP < p2rctIO::TOWERS_IN_REGION_PHI; iPhiInBCP++) {
          float thisTowerEnergy = HCALtowerEnergies[BCPiEta][BCPiPhi][iEtaInBCP][iPhiInBCP];
          int thisFeatureBits = HCALtowerFeatureBits[BCPiEta][BCPiPhi][iEtaInBCP][iPhiInBCP];
          cardsHCAL.addHit(thisTowerEnergy, thisFeatureBits, BCPiEta, BCPiPhi, iEtaInBCP, iPhiInBCP);
        }
      }
    }
  }

  // Loop through cards
  for(int cc=0; cc < p2rctIO::N_CARDS; cc++) {
    //////////////////////////// IP1 ////////////////////////////

    // Separate out the hits into the 5x6 and 2x6 areas
    std::vector<ap_uint<576>> link_in_SLR3_vec = cardsECAL[cc].getAx6(5, 0);
    std::vector<ap_uint<576>> link_in_SLR2_vec = cardsECAL[cc].getAx6(5, 5);
    std::vector<ap_uint<576>> link_in_SLR1_vec = cardsECAL[cc].getAx6(5, 10);
    std::vector<ap_uint<576>> link_in_SLR0_vec = cardsECAL[cc].getAx6(2, 15);

    ap_uint<576>* link_in_SLR3 = &link_in_SLR3_vec[0];
    ap_uint<576>* link_in_SLR2 = &link_in_SLR2_vec[0];
    ap_uint<576>* link_in_SLR1 = &link_in_SLR1_vec[0];
    ap_uint<576>* link_in_SLR0 = &link_in_SLR0_vec[0];

    // Initialize output links for IP1
    ap_uint<576> link_outIP1_SLR3[p2rctIP1_5x6::N_OUTPUT_LINKS];
    ap_uint<576> link_outIP1_SLR2[p2rctIP1_5x6::N_OUTPUT_LINKS];
    ap_uint<576> link_outIP1_SLR1[p2rctIP1_5x6::N_OUTPUT_LINKS];
    ap_uint<576> link_outIP1_SLR0[p2rctIP1_2x6::N_OUTPUT_LINKS];

    // Apply IP1 algo_top
    p2rctIP1_5x6::algo_top(link_in_SLR3, link_outIP1_SLR3);
    p2rctIP1_5x6::algo_top(link_in_SLR2, link_outIP1_SLR2);
    p2rctIP1_5x6::algo_top(link_in_SLR1, link_outIP1_SLR1);
    p2rctIP1_2x6::algo_top(link_in_SLR0, link_outIP1_SLR0);

    // Move cluster outputs into output vectors
    link_outEGammaSLR3->push_back(link_outIP1_SLR3[0]);
    link_outEGammaSLR2->push_back(link_outIP1_SLR2[0]);
    link_outEGammaSLR1->push_back(link_outIP1_SLR1[0]);
    link_outEGammaSLR0->push_back(link_outIP1_SLR0[0]);

    // Move towers of unclustered ECAL energy into output vectors
    link_outECALUnclusteredSLR3->push_back(link_outIP1_SLR3[1]);
    link_outECALUnclusteredSLR2->push_back(link_outIP1_SLR2[1]);
    link_outECALUnclusteredSLR1->push_back(link_outIP1_SLR1[1]);
    link_outECALUnclusteredSLR0->push_back(link_outIP1_SLR0[1]);

    // Declare link_inIP3, which contains HCAL links
    ap_uint<576> link_inIP3[p2rctIP3::N_INPUT_LINKS];

    // HCAL links
    bool secondhalfstarts = (((cc + 3) % 4) > 1); //True for cards 0,3,4,7,etc.
    int BCPiPhiLow = (3*(cc / 2)) / 2;
    int BCPiPhiHigh = BCPiPhiLow + 1;
    if((cc % 2) == 1) {
      // Positive eta
      int BCPiEtaLow = 2;
      int BCPiEtaHigh = 3;

      if (!secondhalfstarts) {
        link_inIP3[5] = cardsHCAL.getLink(BCPiEtaLow,BCPiPhiLow).Data();
        link_inIP3[6] = cardsHCAL.getLink(BCPiEtaHigh,BCPiPhiLow).Data();
        link_inIP3[7] = cardsHCAL.getLink(BCPiEtaLow,BCPiPhiHigh).Data();
        link_inIP3[8] = cardsHCAL.getLink(BCPiEtaHigh,BCPiPhiHigh).Data();
      }
      else {
        link_inIP3[5] = cardsHCAL.getLink(BCPiEtaLow,BCPiPhiHigh).Data();
        link_inIP3[6] = cardsHCAL.getLink(BCPiEtaHigh,BCPiPhiHigh).Data();
        link_inIP3[7] = cardsHCAL.getLink(BCPiEtaLow,BCPiPhiLow).Data();
        link_inIP3[8] = cardsHCAL.getLink(BCPiEtaHigh,BCPiPhiLow).Data();
      }
    }
    else {
      // Negative eta
      int BCPiEtaLow = 0;
      int BCPiEtaHigh = 1;

      if (!secondhalfstarts) {
        link_inIP3[5] = cardsHCAL.getLink(BCPiEtaHigh,BCPiPhiHigh).Data();
        link_inIP3[6] = cardsHCAL.getLink(BCPiEtaLow,BCPiPhiHigh).Data();
        link_inIP3[7] = cardsHCAL.getLink(BCPiEtaHigh,BCPiPhiLow).Data();
        link_inIP3[8] = cardsHCAL.getLink(BCPiEtaLow,BCPiPhiLow).Data();
      }
      else {
        link_inIP3[5] = cardsHCAL.getLink(BCPiEtaHigh,BCPiPhiLow).Data();
        link_inIP3[6] = cardsHCAL.getLink(BCPiEtaLow,BCPiPhiLow).Data();
        link_inIP3[7] = cardsHCAL.getLink(BCPiEtaHigh,BCPiPhiHigh).Data();
        link_inIP3[8] = cardsHCAL.getLink(BCPiEtaLow,BCPiPhiHigh).Data();
      }
    }

    // Move outputs into output vectors
    link_outHCAL8->push_back(link_inIP3[8]);
    link_outHCAL7->push_back(link_inIP3[7]);
    link_outHCAL6->push_back(link_inIP3[6]);
    link_outHCAL5->push_back(link_inIP3[5]);

  }

  iEvent.put(std::move(link_outEGammaSLR3), "EGammaSLR3");
  iEvent.put(std::move(link_outEGammaSLR2), "EGammaSLR2");
  iEvent.put(std::move(link_outEGammaSLR1), "EGammaSLR1");
  iEvent.put(std::move(link_outEGammaSLR0), "EGammaSLR0");
  iEvent.put(std::move(link_outECALUnclusteredSLR3), "ECALUnclusteredSLR3");
  iEvent.put(std::move(link_outECALUnclusteredSLR2), "ECALUnclusteredSLR2");
  iEvent.put(std::move(link_outECALUnclusteredSLR1), "ECALUnclusteredSLR1");
  iEvent.put(std::move(link_outECALUnclusteredSLR0), "ECALUnclusteredSLR0");
  iEvent.put(std::move(link_outHCAL8), "HCAL8");
  iEvent.put(std::move(link_outHCAL7), "HCAL7");
  iEvent.put(std::move(link_outHCAL6), "HCAL6");
  iEvent.put(std::move(link_outHCAL5), "HCAL5");
}

//////////////////////////////////////////////////////////////////////////

void Phase2L1CaloMLonRCTSetupProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("ecalTPEB", edm::InputTag("simEcalEBTriggerPrimitiveDigis"));
  desc.add<edm::InputTag>("hcalTP", edm::InputTag("simHcalTriggerPrimitiveDigis"));
  {
    edm::ParameterSetDescription psd0;
    psd0.add<std::vector<double>>("etaBins",
                                  {
                                      0.087,
                                      0.174,
                                      0.261,
                                      0.348,
                                      0.435,
                                      0.522,
                                      0.609,
                                      0.696,
                                      0.783,
                                      0.87,
                                      0.957,
                                      1.044,
                                      1.131,
                                      1.218,
                                      1.305,
                                      1.392,
                                      1.479,
                                  });
    psd0.add<std::vector<double>>("ptBins",
                                  {
                                      12,
                                      20,
                                      30,
                                      40,
                                      55,
                                      90,
                                      1000000.0,
                                  });
    psd0.add<std::vector<double>>("scale",
                                  {
                                      1.298,  1.287,
                                      1.309,  1.298,
                                      1.309,  1.309,
                                      1.309,  1.298,
                                      1.309,  1.298,
                                      1.309,  1.309,
                                      1.309,  1.32,
                                      1.309,  1.32,
                                      1.309,  1.1742,
                                      1.1639, 1.1639,
                                      1.1639, 1.1639,
                                      1.1639, 1.1639,
                                      1.1742, 1.1742,
                                      1.1639, 1.1639,
                                      1.1742, 1.1639,
                                      1.1639, 1.1742,
                                      1.1742, 1.1536000000000002,
                                      1.11,   1.11,
                                      1.11,   1.11,
                                      1.11,   1.11,
                                      1.11,   1.11,
                                      1.11,   1.11,
                                      1.11,   1.11,
                                      1.11,   1.11,
                                      1.11,   1.11,
                                      1.1,    1.09,
                                      1.09,   1.09,
                                      1.09,   1.09,
                                      1.09,   1.09,
                                      1.09,   1.09,
                                      1.09,   1.09,
                                      1.09,   1.09,
                                      1.09,   1.09,
                                      1.09,   1.09,
                                      1.07,   1.07,
                                      1.07,   1.07,
                                      1.07,   1.07,
                                      1.07,   1.08,
                                      1.07,   1.07,
                                      1.08,   1.08,
                                      1.07,   1.08,
                                      1.08,   1.08,
                                      1.08,   1.06,
                                      1.06,   1.06,
                                      1.06,   1.05,
                                      1.05,   1.06,
                                      1.06,   1.06,
                                      1.06,   1.06,
                                      1.06,   1.06,
                                      1.06,   1.06,
                                      1.06,   1.06,
                                      1.04,   1.04,
                                      1.04,   1.04,
                                      1.05,   1.04,
                                      1.05,   1.05,
                                      1.05,   1.05,
                                      1.05,   1.05,
                                      1.05,   1.05,
                                      1.05,   1.05,
                                      1.05,
                                  });
    desc.add<edm::ParameterSetDescription>("calib", psd0);
  }
  descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1CaloMLonRCTSetupProducer);
