/* 
 * Description: Phase 2 RCT Layer 1 emulator: create ECAL crystal collections
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
// #include <cstdint>
#include <cstdlib> // for rand
#include <iostream>
#include <fstream>
#include <memory>

// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

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
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "L1Trigger/L1CaloTrigger/interface/ParametricCalibration.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Phase2L1CaloEGammaEmulator.h"
#include "Phase2L1RCT.h"
#include "Phase2L1GCT.h"
#include "Phase2L1GCT_algo.h"


// Declare the Phase2L1CaloEGammaEmulator class and its methods

class Phase2L1CaloEGammaEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1CaloEGammaEmulator(const edm::ParameterSet&);
  ~Phase2L1CaloEGammaEmulator() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  // bool passes_ss(float pt, float ss);
  // bool passes_photon(float pt, float pss);
  // bool passes_iso(float pt, float iso);
  // bool passes_looseTkss(float pt, float ss);
  // bool passes_looseTkiso(float pt, float iso);

  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPEBToken_;
  edm::EDGetTokenT<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcalTPToken_;
  edm::ESGetToken<CaloTPGTranscoder, CaloTPGRecord> decoderTag_;

  l1tp2::ParametricCalibration calib_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryTag_;
  const CaloSubdetectorGeometry* ebGeometry;
  const CaloSubdetectorGeometry* hbGeometry;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> hbTopologyTag_;
  const HcalTopology* hcTopology_;
  
};

//////////////////////////////////////////////////////////////////////////

// Phase2L1CaloEGammaEmulator initializer, destructor, and produce methods

Phase2L1CaloEGammaEmulator::Phase2L1CaloEGammaEmulator(const edm::ParameterSet & iConfig)
    : ecalTPEBToken_(consumes<EcalEBTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalTPEB"))),
      hcalTPToken_(
          consumes<edm::SortedCollection<HcalTriggerPrimitiveDigi> >(iConfig.getParameter<edm::InputTag>("hcalTP"))),
      decoderTag_(esConsumes<CaloTPGTranscoder, CaloTPGRecord>(edm::ESInputTag("", ""))),
      calib_(iConfig.getParameter<edm::ParameterSet>("calib")),
      caloGeometryTag_(esConsumes<CaloGeometry, CaloGeometryRecord>(edm::ESInputTag("", ""))),
      hbTopologyTag_(esConsumes<HcalTopology, HcalRecNumberingRecord>(edm::ESInputTag("", ""))) {
  produces<l1tp2::CaloCrystalClusterCollection>( "RCT" );
  produces<l1tp2::CaloCrystalClusterCollection>( "GCT" );
  produces<l1tp2::CaloTowerCollection>( "RCT" );
  produces<l1tp2::CaloTowerCollection>( "GCT" );
  produces<l1tp2::CaloTowerCollection>( "GCTFullTowers" );
}

Phase2L1CaloEGammaEmulator::~Phase2L1CaloEGammaEmulator() {}



void Phase2L1CaloEGammaEmulator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Detector geometry
  const auto& caloGeometry = iSetup.getData(caloGeometryTag_);
  ebGeometry = caloGeometry.getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  hbGeometry = caloGeometry.getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  const auto& hbTopology = iSetup.getData(hbTopologyTag_);
  hcTopology_ = &hbTopology;
  HcalTrigTowerGeometry theTrigTowerGeometry(hcTopology_);

  const auto& decoder = iSetup.getData(decoderTag_);

  //***************************************************// 
  // Declare RCT output collections
  //***************************************************// 

  auto L1EGXtalClusters = std::make_unique<l1tp2::CaloCrystalClusterCollection>();
  auto L1CaloTowers = std::make_unique<l1tp2::CaloTowerCollection>();

  //***************************************************//
  // Get the ECAL hits
  //***************************************************// 
  edm::Handle<EcalEBTrigPrimDigiCollection> pcalohits;   
  iEvent.getByToken(ecalTPEBToken_, pcalohits);

  std::vector<SimpleCaloHit> ecalhits;

  for (const auto& hit : *pcalohits.product()) {
    if (hit.encodedEt() > 0)  // hit.encodedEt() returns an int corresponding to 2x the crystal Et
    {
      // Et is 10 bit, by keeping the ADC saturation Et at 120 GeV it means that you have to divide by 8
      float et = hit.encodedEt() / 8.;
      if (et < cut_500_MeV) {
        continue;  // Reject hits with < 500 MeV ET 
      }

      // Get cell coordinates and info
      auto cell = ebGeometry->getGeometry(hit.id());
      
      std::cout << "Found ECAL cell/hit with et " << et << "and coordinates " << cell->getPosition().x() << "," 
          << cell->getPosition().y() << "," 
          << cell->getPosition().z() << " and ET (GeV) " 
          << et << std::endl;

      SimpleCaloHit ehit;
      ehit.setId(hit.id());
      ehit.setPosition(GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z()));
      ehit.setEnergy(et);
      ehit.setEt_uint((ap_uint<10>) hit.encodedEt());  // also save the uint Et
      ehit.setPt();
      ecalhits.push_back(ehit);	
    }
  }

  // // Debugging only: inject hits in problematic towers
  // std::vector<GlobalVector> myGlobalVectors = {
  //   GlobalVector(93.8306, 88.8616, -1.5158), // (-0.0436499, 0.741765)
  //     GlobalVector( 89.2289, 93.825, 38.6819), // not problematic: at (eta, phi) (0.30555, 0.829031)
  //   GlobalVector(84.1384, 98.077, -3.91203), // problematic: at (eta, phi) (-0.0436499, 0.829031)
  //   GlobalVector(80.3284, 101.263, -11.1033), // (-0.0436499, 0.916298)
  //   GlobalVector(89.2799, 93.7701, -53.0322)
  //                                              };
  // for (auto globVec : myGlobalVectors) {
  //   SimpleCaloHit ehit;
  //   ehit.setPosition(globVec);
  //   ehit.setEnergy(0.5);
  //   ehit.setEt_uint(4);
  //   ehit.setPt();
  //   ecalhits.push_back(ehit);
  // }
  //***************************************************// 
  // Get the HCAL hits
  //***************************************************// 
  std::vector<SimpleCaloHit> hcalhits;
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;
  iEvent.getByToken(hcalTPToken_, hbhecoll);
  
  for (const auto& hit : *hbhecoll.product()) {
    float et = decoder.hcaletValue(hit.id(), hit.t0());
    ap_uint<10> encodedEt = hit.t0().compressedEt(); 
    // same thing as SOI_compressedEt() in HcalTriggerPrimitiveDigi.h///
    if (et <= 0)
      continue;
    
    if (!(hcTopology_->validHT(hit.id()))) {
      LogError("Phase2L1CaloEGammaEmulator")  << " -- Hcal hit DetID not present in HCAL Geom: " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloEGammaEmulator");
      continue;
    }
    const std::vector<HcalDetId>& hcId = theTrigTowerGeometry.detIds(hit.id());
    if (hcId.empty()) {
      LogError("Phase2L1CaloEGammaEmulator") << "Cannot find any HCalDetId corresponding to " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloEGammaEmulator");
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
      
      // std::cout << "Found HCAL cell/TP with coordinates " << cell->getPosition().x() << ","
      //  		<< cell->getPosition().y() << ","
      //  		<< cell->getPosition().z() << " and ET (GeV) " << et
      // 		<< ", encoded Et " << encodedEt << std::endl;
      break;
    }
    SimpleCaloHit hhit;
    hhit.setId(hit.id());
    hhit.setIdHcal(hit.id());
    hhit.setPosition(hcal_tp_position);
    hhit.setEnergy(et);
    hhit.setPt();
    hhit.setEt_uint(encodedEt);
    hcalhits.push_back(hhit);
  }

  //***************************************************// 
  // Initialize necessary arrays for tower and clusters
  //***************************************************// 

  // L1 Outputs definition: Arrays that use firmware convention for indexing
  tower_t towerHCALCard[n_towers_cardEta][n_towers_cardPhi][n_towers_halfPhi]; // 17x4x36 array (not to be confused with the 12x1 array of ap_uints, towerEtHCAL                                               
  tower_t towerECALCard[n_towers_cardEta][n_towers_cardPhi][n_towers_halfPhi];
  // There is one vector of clusters per card (up to 12 clusters before stitching across ECAL regions)
  std::vector<Cluster> cluster_list[n_towers_halfPhi];
  // After merging/stitching the clusters, we only take the 8 highest pt per card                      
  std::vector<Cluster> cluster_list_merged[n_towers_halfPhi];

  //***************************************************// 
  // Fill RCT ECAL regions with ECAL hits
  //***************************************************// 
  for (int cc = 0; cc < n_towers_halfPhi; ++cc) {  // Loop over 36 L1 cards
   
    card rctCard;
    rctCard.setIdx(cc);
        
    for (const auto& hit : ecalhits) {

      // Check if the hit is in cards 0-35
      if ((getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc)) &&
      	  (getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc)) &&
      	  (getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc)) &&
      	  (getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc))) {


        // Get the crystal eta and phi, relative to the bottom left corner of the card 
        // (0 up to 17*5, 0 up to 4*5) 
        int local_iEta = getCrystal_local_iEta(hit.position().eta(), cc);
        int local_iPhi = getCrystal_local_iPhi(hit.position().phi(), cc);

        // Region number (0-5) depends only on the crystal iEta in the card
        int regionNumber = getRegionNumber(local_iEta);
	
        // Tower eta and phi index inside the card (17x4)
        int inCard_tower_iEta = int(local_iEta / CRYSTALS_IN_TOWER_ETA); 
        int inCard_tower_iPhi = int(local_iPhi / CRYSTALS_IN_TOWER_PHI);
	
        // Tower eta and phi index inside the region (3x4)
        int inRegion_tower_iEta = inCard_tower_iEta % TOWER_IN_ETA;
        int inRegion_tower_iPhi = inCard_tower_iPhi % TOWER_IN_PHI;
	
        // Crystal eta and phi index inside the 3x4 region (15x20)
        int inRegion_crystal_iEta = local_iEta % (TOWER_IN_ETA * CRYSTALS_IN_TOWER_ETA);
        int inRegion_crystal_iPhi = local_iPhi;

        // Crystal eta and phi index inside the tower (5x5)
        int inLink_crystal_iEta = (inRegion_crystal_iEta % CRYSTALS_IN_TOWER_ETA);
        int inLink_crystal_iPhi = (inRegion_crystal_iPhi % CRYSTALS_IN_TOWER_PHI);
          
        // Add the crystal energy to the rctCard
        region3x4& myRegion = rctCard.getRegion3x4(regionNumber);
        linkECAL& myLink = myRegion.getLinkECAL(inRegion_tower_iEta, inRegion_tower_iPhi);
	      myLink.addCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi, hit.et_uint());

        std::cout << "[Building ECAL regions:] Added (including calibration here for comparison) " << hit.pt() * calib_(0, std::abs(hit.position().eta())) << " GeV "
            << "at hit.position().eta() and .phi() (" << hit.position().eta() << ", "<< hit.position().phi() << ")"
            << " and GlobalVector position " << hit.position().x() << ", "<< hit.position().y() << ", " << hit.position().z() 
            << " to tower at ("
            << getTowerEta_fromAbsID(getAbsID_iEta_fromFirmwareCardTowerLink(cc, inCard_tower_iEta, inCard_tower_iPhi)) << ", "
            << getTowerPhi_fromAbsID(getAbsID_iPhi_fromFirmwareCardTowerLink(cc, inCard_tower_iEta, inCard_tower_iPhi)) << ")" << ", for an uncalibrated energy of "
            << myLink.getCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi)/8.0 << " GeV"
            << std::endl;
	
      } 
    } 

    //***************************************************// 
    // Build RCT towers from HCAL hits
    //***************************************************// 
    for (const auto& hit : hcalhits) {
      if (getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc) &&
          getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc) &&
          getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc) &&
          getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc) && hit.pt() > 0) {

      	// Get crystal eta and phi, relative to the bottom left corner of the card 
      	// (0 up to 17*5, 0 up to 4*5) 
      	int local_iEta = getCrystal_local_iEta(hit.position().eta(), cc);
      	int local_iPhi = getCrystal_local_iPhi(hit.position().phi(), cc);

        // Region (0-5) the hit falls into 
        int regionNumber = getRegionNumber(local_iEta);
        
        // Tower eta and phi index inside the card (17x4)
        int inCard_tower_iEta = int(local_iEta / CRYSTALS_IN_TOWER_ETA); 
        int inCard_tower_iPhi = int(local_iPhi / CRYSTALS_IN_TOWER_PHI);
        
        // Tower eta and phi index inside the region (3x4)
        int inRegion_tower_iEta = inCard_tower_iEta % TOWER_IN_ETA;
        int inRegion_tower_iPhi = inCard_tower_iPhi % TOWER_IN_PHI;

        // Access the right HCAL region and tower and increment the ET
        towers3x4& myTowers3x4 = rctCard.getTowers3x4(regionNumber);
        towerHCAL& myTower = myTowers3x4.getTowerHCAL(inRegion_tower_iEta, inRegion_tower_iPhi);
        myTower.addEt(hit.et_uint());
      }
    }

    // // Debugging
    // for (int idxRegion = 0; idxRegion < N_REGIONS_PER_CARD; idxRegion++) {
    //   region3x4& myRegion = rctCard.getRegion3x4(idxRegion);
    //   // In each 3x4 region, loop through the links (one link = one tower)
    //   for (int iLinkEta = 0; iLinkEta < TOWER_IN_ETA; iLinkEta++) {
    //     for (int iLinkPhi = 0; iLinkPhi < TOWER_IN_PHI; iLinkPhi++) {
    //       linkECAL& myLink = myRegion.getLinkECAL(iLinkEta, iLinkPhi);
    //       int towerAbsiEta = (getCrystal_iEta_fromCardRegionInfo(cc, idxRegion, iLinkEta, 0)) / 5; // crystal idx doesn't matter
    //       int towerAbsiPhi = (getCrystal_iPhi_fromCardRegionInfo(cc, idxRegion, iLinkPhi, 0)) / 5; // crystal idx doesn't matter

    //       float myTowerEnergy = 0.0;
    //       for (int i = 0; i < 5; i++) {
    //         for (int j = 0; j < 5; j++) {
    //           myTowerEnergy += myLink.getCrystalE(i, j)/8.0;
    //         }
    //       }
    //       // std::cout << "After towerECALCard initialized: tower at ("
    //       //           << getTowerEta_fromAbsID(towerAbsiEta) << ", "
    //       //           << getTowerPhi_fromAbsID(towerAbsiPhi) << "), contains energy: " << myTowerEnergy << std::endl;
    //     }
    //   }
    // }
    //***************************************************// 
    // Make clusters in each ECAL region independently
    //***************************************************// 
   for (int idxRegion = 0; idxRegion < N_REGIONS_PER_CARD; idxRegion++) {
      
      // ECAL crystals array
      crystal temporary[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI]; 
      // HCAL towers in 3x4 region
      ap_uint<12> towerEtHCAL[TOWER_IN_ETA * TOWER_IN_PHI];   
      
      region3x4& myRegion = rctCard.getRegion3x4(idxRegion);
      towers3x4& myTowers = rctCard.getTowers3x4(idxRegion);

    	//	std::cout << std::endl << "[----] DOING CARD " << cc << " AND REGION IDX " << myRegion.getIdx() << std::endl;

      // In each 3x4 region, loop through the links (one link = one tower)
      for (int iLinkEta = 0; iLinkEta < TOWER_IN_ETA; iLinkEta++) {
        for (int iLinkPhi = 0; iLinkPhi < TOWER_IN_PHI; iLinkPhi++) {
          
          // Get the ECAL link (one link per tower) 
          linkECAL& myLink = myRegion.getLinkECAL(iLinkEta, iLinkPhi);

          // // Debugging only:
          // int towerAbsiEta = (getCrystal_iEta_fromCardRegionInfo(cc, idxRegion, iLinkEta, 0)) / 5; // crystal idx doesn't matter
          // int towerAbsiPhi = (getCrystal_iPhi_fromCardRegionInfo(cc, idxRegion, iLinkPhi, 0)) / 5; // crystal idx doesn't matter
          // float myTowerEnergy = 0.0;
          // for (int i = 0; i < 5; i++) {
          //   for (int j = 0; j < 5; j++) {
          //     myTowerEnergy += myLink.getCrystalE(i, j)/8.0;
          //   }
          // }
          // std::cout << "(Before packaging into temp array): in card " << cc << " region " << idxRegion << ", tower at ("
          //           << getTowerEta_fromAbsID(towerAbsiEta) << ", "
          //           << getTowerPhi_fromAbsID(towerAbsiPhi) << "), contains energy: " << myTowerEnergy << std::endl;
          // end of debugging block
          
          // We have an array of 3x4 links/towers, each link/tower is 5x5 in crystals. We need to convert this to a 15x20 of crystals
          int ref_iEta = (iLinkEta * CRYSTALS_IN_TOWER_ETA); 
          int ref_iPhi = (iLinkPhi * CRYSTALS_IN_TOWER_PHI); 

          // In the link, get the crystals (5x5 in each link) 
          for (int iEta = 0; iEta < CRYSTALS_IN_TOWER_ETA; iEta++) {           	    
            for (int iPhi = 0; iPhi < CRYSTALS_IN_TOWER_PHI; iPhi++) {    
        
              // Et as unsigned int
              ap_uint<10> uEnergy = myLink.getCrystalE(iEta, iPhi);

              // Fill the 'temporary' array with a crystal object 
              temporary[ref_iEta + iEta][ref_iPhi + iPhi] = crystal(uEnergy);
	          }
	        } // end of loop over crystals

          // HCAL tower ET
          towerHCAL& myTower = myTowers.getTowerHCAL(iLinkEta, iLinkPhi);
          towerEtHCAL[(iLinkEta * TOWER_IN_PHI) + iLinkPhi] = myTower.getEt();
	      }
	    }

      // Debugging only:
      std::cout << "........ Before clustering, card " << cc << " and region " << idxRegion;
      float remainingTowerEnergies = 0.0;
      for (int i = 0; i < 15; i++) {
        for (int j = 0; j < 20; j++) {
          remainingTowerEnergies += temporary[i][j].energy;
        }
      }
      std::cout << ", has unclustered energy (GeV) " << remainingTowerEnergies/8.0 << std::endl;
	
      // Iteratively find four clusters and remove them from 'temporary' as we go, and fill cluster_list
      for (int c = 0; c < N_CLUSTERS_PER_REGION; c++) {
        std::cout << "... Iteratively finding clusters: loop number " << c << " in card " << cc << " and region " << idxRegion << std::endl;
        Cluster newCluster = getClusterFromRegion3x4(temporary); // remove cluster from 'temporary' 
        newCluster.setRegionIdx(idxRegion);                      // add the region number 
        if (newCluster.clusterEnergy() > 0) {
          cluster_list[cc].push_back(newCluster);                // do not push back 0-energy clusters
        }
        // Debugging: iteratively print out temp. array here too
        std::cout << "........ End of loop number " << c << " in card " << cc << " and region " << idxRegion << std::endl;
        float remainingE = 0.0;
        for (int i = 0; i < 15; i++) {
          for (int j = 0; j < 20; j++) {
            remainingE += temporary[i][j].energy;
          }
        }
        std::cout << ", region has unclustered energy (GeV) " << remainingE/8.0 << std::endl;
        // end of debugging block
      }
	
      // Create towers using remaining ECAL energy, and the HCAL towers were already calculated in towersEtHCAL[12] 
      ap_uint<12> towerEtECAL[12];
      getECALTowersEt(temporary, towerEtECAL);

      // // Debugging only
      // float twelveListRemainingE = 0.0;

      // Fill towerHCALCard and towerECALCard arrays 
      for (int i = 0; i < 12; i++) {
        // Get the tower's indices in a (17x4) card
        int iEta = (idxRegion * TOWER_IN_ETA) + (i / TOWER_IN_PHI);   
        int iPhi = (i % TOWER_IN_PHI);

        // If the region number is 5 (i.e. only 2x4 in towers, skip the third row) N_REGIONS_PER_CARD = 6. i.e. we do not want to consider
        // i = 8, 9, 10, 11
        if ( (idxRegion == (N_REGIONS_PER_CARD - 1)) && (i > 7) ) {
          continue;
        }
        towerHCALCard[iEta][iPhi][cc] = tower_t(towerEtHCAL[i], 0); // tower_t initializer takes an ap-uint<12> for the energy
        towerECALCard[iEta][iPhi][cc] = tower_t(towerEtECAL[i], 0);

        // std::cout << "[Inside filling:] card " << cc << " region " << idxRegion << ", i = " << i << ", iEta = " << iEta << ", iPhi = " << iPhi << ", tower has towerEtECAL[i] = " 
        //           << towerEtECAL[i] << " mapping to towerECALCard[iEta][iPhi][cc] = " << towerECALCard[iEta][iPhi][cc].et() << std::endl;

        // // Debugging only
        // twelveListRemainingE += (towerEtECAL[i]/8.0);

      }
      // std::cout << "... After mysteriously repackaging into an array, the total energy of card " << cc << " region " << idxRegion << " is " << twelveListRemainingE << std::endl;

    }


    //-------------------------------------------// 
    // Stitching across ECAL regions             //
    //-------------------------------------------// 
    
    // TESTING ONLY: use dummy clusters                                                                                          
    // for (int i = 0; i < 12; i++) {                                                                                            
    //   Cluster dummy = Cluster((ap_uint<12>) i+5, i%12, i%4, 0, 0, 0);                                                   
    //   test_cluster.push_back(dummy);                                                                                         
    // }                    

    // TESTING ONLY: Put two fake clusters on either side of a boundary: e.g. tower eta 15, tower phi 1, crystal eta 0, crystal phi 3
    // and the other directly below it in eta, but either 0 or 1 away in crystal phi (first try with same phi)
    // for (int i = 0; i < 5; i++) {
    //   int dummyTowerPhi = (i % 4);
    //   int dummyCrystalPhi = (i % 5);
    //   Cluster dummy1 = Cluster((ap_uint<12>) 25, 0, dummyTowerPhi, 0, dummyCrystalPhi, 0);
    //   dummy1.regionIdx = (i+1);
    //   Cluster dummy2 = Cluster((ap_uint<12>) 20, 2, dummyTowerPhi, 4, dummyCrystalPhi + 1, 1); 
    //   dummy2.regionIdx = i;
    //   test_cluster.push_back(dummy1);
    //   test_cluster.push_back(dummy2);
    // }
    // START HERE: Replace cluster_list[cc] with test_cluster to do a test  

   // std::cout << "Card " << cc << ": BEFORE stitching: " << std::endl;
    
    // for ( Cluster c: cluster_list[cc] ) {
    //   c.printInfo("Before stitching");
    // }
    int nRegionBoundariesEta = (N_REGIONS_PER_CARD - 1); // 6 regions -> 5 boundaries to check
    // Upper and lower boundaries respectively, to check for stitching along
    int towerEtaBoundaries[nRegionBoundariesEta][2] = 
      {	{15, 14},
        {12, 11},
        {9, 8},
        {6, 5},
        {3, 2}      };

    for (int iBound = 0; iBound < nRegionBoundariesEta; iBound++) {
      stitchClusterOverRegionBoundary(cluster_list[cc], towerEtaBoundaries[iBound][0], towerEtaBoundaries[iBound][1], cc);
    }
    
    // std::cout << "Card " << cc << ": AFTER stitching:" << std::endl;
    // for ( Cluster c : cluster_list[cc] ) {
    //   c.printInfo("After stitching");
    // }

    //--------------------------------------------------------------------------------//  
    // Sort the clusters, take the 8 with highest pT, and return extras to tower 
    //--------------------------------------------------------------------------------//  
    if (cluster_list[cc].empty()) {
      // << "No clusters found in card " << cc << std::endl;
    }
    else {
      std::sort(cluster_list[cc].begin(), cluster_list[cc].end(), compareClusterET);
      
      // If there are more than eight clusters, return the unused energy to the towers
      for (unsigned int kk = n_clusters_4link; kk < cluster_list[cc].size(); ++kk) {
        Cluster cExtra = cluster_list[cc][kk];
        if (cExtra.clusterEnergy() > 0) {

          // Increment tower ET 
          // Get tower (eta, phi) (up to (17, 4)) in the RCT card
          int whichTowerEtaInCard = ((cExtra.region() * TOWER_IN_ETA) + cExtra.towerEta());
          int whichTowerPhiInCard = cExtra.towerPhi(); 
          ap_uint<12> oldTowerEt = towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc].et();
          ap_uint<12> newTowerEt = (oldTowerEt + cExtra.clusterEnergy());
          ap_uint<4>  hoe        = towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc].hoe();
          towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc] = tower_t(newTowerEt, hoe);

          std::cout << "[Adding clusters back to towers:] Added extra cluster # " << kk << ": energy " << cExtra.clusterEnergy()/8.0 << " GeV, at ("
                    << getEta_fromCrystaliEta( getCrystal_iEta_fromCardRegionInfo(cc, cExtra.region(), cExtra.towerEta(), cExtra.clusterEta()) ) << ", "
                    << getPhi_fromCrystaliPhi( getCrystal_iEta_fromCardRegionInfo(cc, cExtra.region(), cExtra.towerPhi(), cExtra.clusterPhi()) )
                    << " to tower with real (iEta, iPhi) "
                    << getTowerEta_fromAbsID(getAbsID_iEta_fromFirmwareCardTowerLink(cc, whichTowerEtaInCard, whichTowerPhiInCard)) << ", "
                    << getTowerPhi_fromAbsID(getAbsID_iPhi_fromFirmwareCardTowerLink(cc, whichTowerEtaInCard, whichTowerPhiInCard)) << ")"
                    << " with original energy " << oldTowerEt/8.0 << " and new energy " << newTowerEt/8.0
                    << std::endl;
          
        }
      }
      
      // Save up to eight clusters: loop over cluster_list
      for (unsigned int kk = 0; kk < cluster_list[cc].size(); ++kk) {
        if (kk >= n_clusters_4link) continue;
        if (cluster_list[cc][kk].clusterEnergy() > 0) {
          cluster_list_merged[cc].push_back(cluster_list[cc][kk]);
        }
      }
      // std::cout << "Sorted, up to 8 clusters: ";
      // for ( Cluster c : cluster_list_merged[cc] ) {
      //   c.printInfo("cluster_list_merged before calibration");
      // }
      // std::cout << std::endl;
      
    }

    //-------------------------------------------//                          
    // Calibrate clusters                    
    //-------------------------------------------//   
    for (auto &c : cluster_list_merged[cc]) {
      float realEta = getEta_fromCrystaliEta(getCrystal_iEta_fromCardRegionInfo(cc, 
                                                                                c.region(),
                                                                                c.towerEta(),
                                                                                c.clusterEta()));
      
      // std::cout << "Calib with pT " << c.getPt() << " and eta " << realEta << ": Old pT: "
		  //           << c.getPt() << " (uint:" << c.clusterEnergy() << ")" << std::endl;

      c.calib = calib_(c.getPt(), std::abs(realEta));
      c.applyCalibration(c.calib);

      // std::cout << ">>> New pT after c.applyCalibration: " << c.getPt() << " (uint:" << c.clusterEnergy() << ")" << std::endl;
      // c.printInfo("cluster_list_merged after calibration");
      
    }

    //-------------------------------------------//                          
    // Cluster shower shape flags
    //-------------------------------------------// 
    // Shower shape flags
    for (auto &c : cluster_list_merged[cc]) { 
      c.is_ss   = passes_ss(c.getPt(), (c.getEt2x5() / c.getEt5x5()) );
      c.is_looseTkss = passes_looseTkss(c.getPt(), (c.getEt2x5() / c.getEt5x5() ));
    }

    //-------------------------------------------//                          
    // Calibrate towers                   
    //-------------------------------------------// 
    for (int ii = 0; ii < n_towers_cardEta; ++ii) { // 17 towers per card in eta
      for (int jj = 0; jj < n_towers_cardPhi; ++jj) { // 4 towers per card in phi
        float tRealEta = getTowerEta_fromAbsID(getAbsID_iEta_fromFirmwareCardTowerLink(cc, ii, jj));  // real eta of center of tower
        double tCalib = calib_(0, tRealEta);                                    // calibration factor
        bool verbose = false;
        towerECALCard[ii][jj][cc].applyCalibration(tCalib, verbose);
      }
    }

    //-------------------------------------------//                          
    // Calculate tower HoE                    
    //-------------------------------------------// 
    // std::cout << "Computing HoE for towers: " << std::endl;
    for (int ii = 0; ii < n_towers_cardEta; ++ii) { // 17 towers per card in eta
      for (int jj = 0; jj < n_towers_cardPhi; ++jj) { // 4 towers per card in phi   
        ap_uint<12> ecalEt = towerECALCard[ii][jj][cc].et();
        ap_uint<12> hcalEt = towerHCALCard[ii][jj][cc].et();
        towerECALCard[ii][jj][cc].getHoverE(ecalEt, hcalEt); 
      }
    }

    //-----------------------------------------------------------//
    // Produce output collections for event display and analyzer
    //-----------------------------------------------------------//
    for (auto & c : cluster_list_merged[cc]) { 
      // c.printInfo("Start of loop to output collections");                                                    

      float realEta = getEta_fromCrystaliEta( getCrystal_iEta_fromCardRegionInfo(cc, c.region(), c.towerEta(), c.clusterEta()) );
      float realPhi = getPhi_fromCrystaliPhi( getCrystal_iPhi_fromCardRegionInfo(cc, c.region(), c.towerPhi(), c.clusterPhi()) );

      reco::Candidate::PolarLorentzVector p4calibrated(c.getPt(), realEta, realPhi, 0.);

      // Constructor definition at: https://github.com/cms-l1t-offline/cmssw/blob/l1t-phase2-v3.3.11/DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h#L34
      l1tp2::CaloCrystalCluster cluster(p4calibrated,
                                        c.getPt(), // use float
                                        0,  // float h over e
                                        0,  // float iso
                                        0,  // DetId seedCrystal 
                                        0,  // puCorrPt
                                        c.getBrems(), // 0, 1, or 2 (as computed in firmware)
                                        0,            // et2x2 (not calculated)
                                        c.getEt2x5(), // et2x5 (as computed in firmware, save float)
                                        0,            // et3x5 (not calculated)
                                        c.getEt5x5(), // et5x5 (as computed in firmware, save float)
                                        c.getIsSS(), // standalone WP
                                        c.getIsSS(), // electronWP98
                                        false, // photonWP80
                                        c.getIsSS(), // electronWP90
                                        c.getIsLooseTkss(), // looseL1TkMatchWP
                                        c.getIsSS()  // stage2effMatch
                                        );

      std::map<std::string, float> params;
      params["standaloneWP_showerShape"] = c.getIsSS();
      params["trkMatchWP_showerShape"]   = c.getIsLooseTkss();
      cluster.setExperimentalParams(params);
      
      L1EGXtalClusters->push_back(cluster);
    }  
    // Output tower collections
    for (int ii = 0; ii < n_towers_cardEta; ++ii) { // 17 towers per card in eta
      for (int jj = 0; jj < n_towers_cardPhi; ++jj) { // 4 towers per card in phi

        l1tp2::CaloTower l1CaloTower;
        // Divide by 8.0 to get ET as float (GeV)
        l1CaloTower.setEcalTowerEt(towerECALCard[ii][jj][cc].et()/8.0);
        // HCAL TPGs encoded ET: multiply by the LSB (0.5) to convert to GeV
        float hcalLSB = 0.5;
        l1CaloTower.setHcalTowerEt(towerHCALCard[ii][jj][cc].et() * hcalLSB);
        int absToweriEta = getAbsID_iEta_fromFirmwareCardTowerLink(cc, ii, jj);
        int absToweriPhi = getAbsID_iPhi_fromFirmwareCardTowerLink(cc, ii, jj);
        l1CaloTower.setTowerIEta(absToweriEta);
        l1CaloTower.setTowerIPhi(absToweriPhi);
        l1CaloTower.setTowerEta(getTowerEta_fromAbsID(absToweriEta));
        l1CaloTower.setTowerPhi(getTowerPhi_fromAbsID(absToweriPhi));
        
        L1CaloTowers->push_back(l1CaloTower);
      }
    }


  } // end of loop over cards

  iEvent.put(std::move(L1EGXtalClusters), "RCT");
  iEvent.put(std::move(L1CaloTowers),     "RCT");
  std::cout << "Finished producing RCT emulator for this event" << std::endl;
  
  //*******************************************************************
  // Do GCT geometry for inputs
  //*******************************************************************

  // Loop over GCT cards (three of them)
  GCTcard_t gctCards[N_GCTCARDS];
  GCTtoCorr_t gctToCorr[N_GCTCARDS];

  for (unsigned int gcc = 0; gcc < N_GCTCARDS; gcc++) {
    // Each GCT card encompasses 16 RCT cards, listed in 
    // GCTcardtoRCTcardnumber[3][16]. 
    std::cout << "[Phase2L1CaloEGammaEmulator.cc] Packaging GCT: Starting GCT Card " << gcc << "..." << std::endl;
    // i goes from 0 to <16
    for (int i = 0; i < (N_RCTCARDS_PHI * 2); i++) {

      unsigned int rcc = GCTcardtoRCTcardnumber[gcc][i];
      std::cout << "[Phase2L1CaloEGammaEmulator.cc] Packaging GCT ... Starting RCT Card: " << rcc << "..." << std::endl;

      // Positive eta? Fist row is in positive eta
      bool isPositiveEta = (i < N_RCTCARDS_PHI);
      
      // Sum tower ECAL and HCAL energies
      // 17 towers per link
      for (int iTower = 0; iTower < N_GCTTOWERS_FIBER; iTower++) {
        // 4 links per card
        for (int iLink = 0; iLink < n_links_card; iLink++) {
          tower_t t0_ecal = towerECALCard[iTower][iLink][rcc];
          tower_t t0_hcal = towerHCALCard[iTower][iLink][rcc];
          RCTtower_t t;
          t.et = t0_ecal.et() + convertHcalETtoEcalET(t0_hcal.et());
          t.hoe = t0_ecal.hoe();
          // Not needed for GCT firmware but will be written into GCT CMSSW outputs : 12 bits each 
          t.ecalEt = t0_ecal.et();
          t.hcalEt = t0_hcal.et();
          // std::cout << "tower et and hoe: " << t.et << ", " << t.hoe 
          //           << " from ECAL (uint) " << t.ecalEt << " and HCAL (using HCAL uint) " << t.hcalEt << std::endl;

          if (isPositiveEta) {  
            gctCards[gcc].RCTcardEtaPos[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower] = t;
          } else {
            gctCards[gcc].RCTcardEtaNeg[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower] = t;
          }

          // std::cout << "Stashing tower with ECAL energy " << t0_ecal.et()/8.0 << " GeV with real (eta, phi): ("
          //           << getTowerEta_fromAbsID(getAbsID_iEta_fromFirmwareCardTowerLink(rcc, iTower, iLink)) << ", "
          //           << getTowerPhi_fromAbsID(getAbsID_iPhi_fromFirmwareCardTowerLink(rcc, iTower, iLink)) << ")"
          //           << std::endl;

	      }
      } 
      
      // Distribute at most 8 RCT clusters across four links: convert to GCT coordinates
      for (size_t iCluster = 0; (iCluster < cluster_list_merged[rcc].size()) && (iCluster < (N_RCTGCT_FIBERS * N_RCTCLUSTERS_FIBER)); iCluster++) {
        Cluster c0 = cluster_list_merged[rcc][iCluster];
        RCTcluster_t c;
        c.et     = c0.clusterEnergy();
        
        // tower Eta: c0.towerEta() refers to the tower iEta INSIDE the region, so we need to convert to tower iEta inside the card
        c.towEta = (c0.region() * TOWER_IN_ETA) + c0.towerEta();     
        c.towPhi = c0.towerPhi();
        c.crEta  = c0.clusterEta();
        c.crPhi  = c0.clusterPhi();
        c.et5x5  = c0.uint_et5x5(); 
        c.et2x5  = c0.uint_et2x5(); 
        c.is_ss        = c0.getIsSS();      
        c.is_looseTkss = c0.getIsLooseTkss(); 
	
        unsigned int iIdxInGCT = i % N_RCTCARDS_PHI;  
        unsigned int iLinkC = iCluster % N_RCTGCT_FIBERS;
        unsigned int iPosC  = iCluster / N_RCTGCT_FIBERS;
        // std::cout << c.et << ", "
        //           << "accessing link " << iCluster % N_RCTGCT_FIBERS << " "
        //           << "and position "   << iCluster / N_RCTGCT_FIBERS << " "
        //           << "isPositiveEta "  << isPositiveEta << " " 
        //           << "with iIdxInGCT " << iIdxInGCT << " " 
        //           << "with uint et5x5, et2x5 " << c.et5x5 << ", " << c.et2x5 << " " 
        //           << "with shower shape flags ss/looseTkss " << c.is_ss << ", " << c.is_looseTkss << std::endl;
	
        if (isPositiveEta) {
          gctCards[gcc].RCTcardEtaPos[iIdxInGCT].RCTtoGCTfiber[iLinkC].RCTclusters[iPosC] = c;
        } else {
          gctCards[gcc].RCTcardEtaNeg[iIdxInGCT].RCTtoGCTfiber[iLinkC].RCTclusters[iPosC] = c;
        }
      }
      // If there were fewer than eight clusters, make sure the remaining fiber clusters are zero'd out.
      for (size_t iZeroCluster = cluster_list_merged[rcc].size(); iZeroCluster < (N_RCTGCT_FIBERS * N_RCTCLUSTERS_FIBER); iZeroCluster++) {
	      unsigned int iIdxInGCT = i % N_RCTCARDS_PHI;
        unsigned int iLinkC = iZeroCluster % N_RCTGCT_FIBERS;
        unsigned int iPosC  = iZeroCluster / N_RCTGCT_FIBERS;

	      RCTcluster_t cZero;
        cZero.et     = 0;
        cZero.towEta = 0;
        cZero.towPhi = 0;
        cZero.crEta  = 0;
        cZero.crPhi  = 0;
        cZero.et5x5  = 0;
        cZero.et2x5  = 0;
        cZero.is_ss  = false;
        cZero.is_looseTkss = false;
	      if (isPositiveEta) {
          gctCards[gcc].RCTcardEtaPos[iIdxInGCT].RCTtoGCTfiber[iLinkC].RCTclusters[iPosC] = cZero;
        } 
        else {
          gctCards[gcc].RCTcardEtaNeg[iIdxInGCT].RCTtoGCTfiber[iLinkC].RCTclusters[iPosC] = cZero;
        }
      }
    }
  } // end of loop over initializing GCT cards
	
  // Sanity check: go back into the GCT card and see if we get values
  for (unsigned int gcc = 0; gcc < N_GCTCARDS; gcc++) {
    //    std::cout << "GCT: Starting Card " << gcc << "..." << std::endl;
    
    for (int i = 0; i < (N_RCTCARDS_PHI); i++) {
      //       std::cout << "Inside GCT: card (out of 8 in the positive side) " << i << std::endl;
      // for (int iLink = 0; iLink < n_links_card; iLink++) {
      //   for (int iTower = 0; iTower < N_GCTTOWERS_FIBER; iTower++) {
      //     	std::cout << gctCards[gcc].RCTcardEtaPos[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower].et << " = (ECAL) "
      //               << gctCards[gcc].RCTcardEtaPos[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower].ecalEt << " + (HCAL) "
      //               << gctCards[gcc].RCTcardEtaPos[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower].hcalEt << std::endl;
	    //     }
      //   }
      }
    // std::cout << std::endl;
    // std::cout << "Clusters: ";
    for (int iLink = 0; iLink < n_links_card; iLink++) {
      for (int iCluster = 0; iCluster < N_RCTCLUSTERS_FIBER; iCluster++) {
        // std::cout << gctCards[gcc].RCTcardEtaPos[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTclusters[iCluster].et << ", ";
      }
    }
    
    std::cout << std::endl;
    for (int i = N_RCTCARDS_PHI; i < (N_RCTCARDS_PHI * 2); i++) {
      // std::cout<< "Inside GCT: card (out of 8 in the negative side) " << i % N_RCTCARDS_PHI << std::endl;
	    // for (int iLink = 0; iLink < n_links_card; iLink++) {
	    //   for (int iTower = 0; iTower < N_GCTTOWERS_FIBER; iTower++) {
	    //     std::cout << gctCards[gcc].RCTcardEtaNeg[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower].et << " = (ECAL) "
      //               << gctCards[gcc].RCTcardEtaNeg[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower].ecalEt << " + (HCAL in HCAL convention) "
      //               << gctCards[gcc].RCTcardEtaNeg[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower].hcalEt << std::endl;
	    //   }
	    // }
      //	std::cout << std::endl;
      //	std::cout << "Clusters: ";
      for (int iLink = 0; iLink < n_links_card; iLink++) {
        for (int iCluster = 0; iCluster < N_RCTCLUSTERS_FIBER; iCluster++) {
          //  std::cout << gctCards[gcc].RCTcardEtaNeg[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTclusters[iCluster].et << ", ";
        }
      }
    }
    
    std::cout << std::endl;
  }	   
  std::cout << std::endl;

  //----------------------------------------------------
  // Output collections for the GCT emulator
  //----------------------------------------------------   
  auto L1GCTClusters   = std::make_unique<l1tp2::CaloCrystalClusterCollection>();
  auto L1GCTTowers     = std::make_unique<l1tp2::CaloTowerCollection>();
  auto L1GCTFullTowers = std::make_unique<l1tp2::CaloTowerCollection>();

  // Apply the GCT firmware code to each GCT
  for (unsigned int gcc = 0; gcc < N_GCTCARDS; gcc++) {
    // getClustersTowers
    // getClustersCombined
    // getFullTowers
    algo_top(gctCards[gcc], gctToCorr[gcc], gcc, L1GCTClusters, L1GCTTowers, L1GCTFullTowers, calib_);
  }

  // Check that L1GCTClusters has stuff in it
  std::cout << "Checking that the GCT output collections has stuff in them..." << std::endl;
  for (auto const& c: *L1GCTClusters ) {

    std::cout << "GCT cluster energy " << c.pt() << std::endl;
  }
  // for (auto const& c: *L1GCTTowers ) {
    
  //   std::cout << "GCT tower energy " << c.ecalTowerEt() + c.hcalTowerEt()
	//       << " at real eta, phi: (" << c.towerEta() << ", " << c.towerPhi() << ")" << std::endl;
    
  // }

  // std::cout << "Checking that the GCT Full Towers have entries..." << std::endl;
  // for (auto const& c: *L1GCTFullTowers ) {
    
  //   std::cout << "GCT FULL tower energy " << c.ecalTowerEt() 
	//             << " at real eta, phi: (" << c.towerEta() << ", " << c.towerPhi() << ")" << std::endl;
    
  // }
  
  
  iEvent.put(std::move(L1GCTClusters), "GCT");
  iEvent.put(std::move(L1GCTTowers),   "GCT");
  iEvent.put(std::move(L1GCTFullTowers), "GCTFullTowers");
} 



////////////////////////////////////////////////////////////////////////// 

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1CaloEGammaEmulator);

