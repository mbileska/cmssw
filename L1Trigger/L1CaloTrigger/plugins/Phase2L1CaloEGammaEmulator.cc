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
  edm::ESHandle<CaloTPGTranscoder> decoder_;

  l1tp2::ParametricCalibration calib_;

  edm::ESHandle<CaloGeometry> caloGeometry_;
  const CaloSubdetectorGeometry* ebGeometry;
  const CaloSubdetectorGeometry* hbGeometry;
  edm::ESHandle<HcalTopology> hbTopology;
  const HcalTopology* hcTopology_;
  
};

//////////////////////////////////////////////////////////////////////////

// Phase2L1CaloEGammaEmulator initializer, destructor, and produce methods

Phase2L1CaloEGammaEmulator::Phase2L1CaloEGammaEmulator(const edm::ParameterSet & iConfig)
  : ecalTPEBToken_(consumes<EcalEBTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("ecalTPEB"))),
    hcalTPToken_(
		 consumes<edm::SortedCollection<HcalTriggerPrimitiveDigi> >(iConfig.getParameter<edm::InputTag>("hcalTP"))),
  calib_(iConfig.getParameter<edm::ParameterSet>("calib")) {
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
  iSetup.get<CaloGeometryRecord>().get(caloGeometry_);
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
  hbGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Hcal, HcalBarrel);
  iSetup.get<HcalRecNumberingRecord>().get(hbTopology);
  hcTopology_ = hbTopology.product();
  HcalTrigTowerGeometry theTrigTowerGeometry(hcTopology_);
  iSetup.get<CaloTPGRecord>().get(decoder_);

  //***************************************************// 
  // Declare RCT output collections
  //***************************************************// 

  auto L1EGXtalClusters = std::make_unique<l1tp2::CaloCrystalClusterCollection>();
  auto L1CaloTowers = std::make_unique<l1tp2::CaloTowerCollection>();

  //***************************************************//
  // Get the ECAL hits
  //***************************************************// 
  edm::Handle<EcalEBTrigPrimDigiCollection> pcalohits;   // from l.354
  iEvent.getByToken(ecalTPEBToken_, pcalohits);

  std::vector<SimpleCaloHit> ecalhits;

  for (const auto& hit : *pcalohits.product()) {
    if (hit.encodedEt() > 0)  // hit.encodedEt() returns an int corresponding to 2x the crystal Et
      {
	// Et is 10 bit, by keeping the ADC saturation Et at 120 GeV it means that you have to divide by 8
	float et = hit.encodedEt() / 8.;
	if (et < cut_500_MeV)
	  continue;  // keep the 500 MeV ET Cut

	// std::cout << "ECAL hit et: " << et << std::endl;

	// Get cell coordinates and info
	auto cell = ebGeometry->getGeometry(hit.id());
	
	// std::cout << "Found ECAL cell/hit with coordinates " << cell->getPosition().x() << "," 
	// 	  << cell->getPosition().y() << "," 
	// 	  << cell->getPosition().z() << " and ET (GeV) " 
	// 	  << et << std::endl;

	SimpleCaloHit ehit;
	ehit.setId(hit.id());
	ehit.setPosition(GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z()));
	ehit.setEnergy(et);
	ehit.setEt_uint((ap_uint<10>) hit.encodedEt());  // also save the uint Et
	ehit.setPt();
	ecalhits.push_back(ehit);

	
      }
  }

  //***************************************************// 
  // Get the HCAL hits
  //***************************************************// 
  std::vector<SimpleCaloHit> hcalhits;
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hbhecoll;
  iEvent.getByToken(hcalTPToken_, hbhecoll);
  
  for (const auto& hit : *hbhecoll.product()) {
    float et = decoder_->hcaletValue(hit.id(), hit.t0());
    ap_uint<10> encodedEt = hit.t0().compressedEt(); 
    // same thing as SOI_compressedEt() in HcalTriggerPrimitiveDigi.h///
    if (et <= 0)
      continue;
    
    if (!(hcTopology_->validHT(hit.id()))) {
      LogError("Phase2L1CaloEGammaEmulator")
  	<< " -- Hcal hit DetID not present in HCAL Geom: " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloEGammaEmulator");
      continue;
    }
    const std::vector<HcalDetId>& hcId = theTrigTowerGeometry.detIds(hit.id());
    if (hcId.empty()) {
      LogError("Phase2L1CaloEGammaEmulator")
  	<< "Cannot find any HCalDetId corresponding to " << hit.id() << std::endl;
      throw cms::Exception("Phase2L1CaloEGammaEmulator");
      continue;
    }
    if (hcId[0].subdetId() > 1)
      continue;
    GlobalVector hcal_tp_position = GlobalVector(0., 0., 0.);
    for (const auto& hcId_i : hcId) {
      if (hcId_i.subdetId() > 1)
        continue;
      // get the first HCAL TP/ cell
      auto cell = hbGeometry->getGeometry(hcId_i);
      if (cell == nullptr)
  	continue;
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
  // Declare old-CMSSW outputs
  //***************************************************// 
  //
  // L1 Outputs definition (using previous emulator by Cecile) - NOT firmware convention
  // Firmware convention -> CMSSW indexing conversion is done near the end.
  // 36 L1 cards, each with 4x17 towers. All using CMSSW indexing convention, NOT firmware convention
  //
  ap_uint<12> ECAL_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  ap_uint<12> HCAL_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  int iEta_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];  
  int iPhi_tower_L1Card[n_links_card][n_towers_per_link][n_towers_halfPhi];
  // 36 L1 cards with 4 links, each with 2 clusters (up to 8 per card)
  ap_uint<12> energy_cluster_L1Card[n_links_card][n_clusters_link][n_towers_halfPhi]; 
  int crystalID_cluster_L1Card[n_links_card][n_clusters_link][n_towers_halfPhi];  // range: [0, 5*5) 
  int towerID_cluster_L1Card[n_links_card][n_clusters_link][n_towers_halfPhi];    // range: [0, 17*4)

  //
  // L1 Outputs definition: Arrays that use firmware convention for indexing
  //
  tower_t towerHCALCard[n_towers_cardEta][n_towers_cardPhi][n_towers_halfPhi]; // 17x4x36 array (not to be confused with the 12x1 array of ap_uints, towerEtHCAL                                               
  tower_t towerECALCard[n_towers_cardEta][n_towers_cardPhi][n_towers_halfPhi];
  // There is one vector of clusters per card (up to 12 clusters before stitching across ECAL regions)
  std::vector<Cluster> cluster_list[n_towers_halfPhi];
  // After merging/stitching the clusters, we only take the 8 highest pt per card                      
  std::vector<Cluster> cluster_list_merged[n_towers_halfPhi];


  // Zero out the L1 outputs
  for (int ii = 0; ii < n_links_card; ++ii) {
    for (int jj = 0; jj < n_towers_per_link; ++jj) {
      for (int ll = 0; ll < n_towers_halfPhi; ++ll) {
        ECAL_tower_L1Card[ii][jj][ll] = 0;
        HCAL_tower_L1Card[ii][jj][ll] = 0;
        iPhi_tower_L1Card[ii][jj][ll] = -999;
        iEta_tower_L1Card[ii][jj][ll] = -999;
      }
    }
  }

  for (int ii = 0; ii < n_links_card; ++ii) {
    for (int jj = 0; jj < n_clusters_link; ++jj) {
      for (int ll = 0; ll < n_towers_halfPhi; ++ll) {
	energy_cluster_L1Card[ii][jj][ll]    = 0;
	towerID_cluster_L1Card[ii][jj][ll]   = 0;
        crystalID_cluster_L1Card[ii][jj][ll] = 0;
	
      }
    }
  }

  //***************************************************// 
  // Fill RCT emulator with ECAL hits
  //***************************************************// 
  int theCard = 28; // Debugging only
  
  for (int cc = 0; cc < n_towers_halfPhi; ++cc) {  // Loop over 36 L1 cards
    // Debugging purposes only: TEMP: only do one card
    // if (cc != theCard) continue;
    
    // Initialize variables
    card rctCard;
    rctCard.setIdx(cc);
    
    // Debugging only: inject four hits at specific locations
    // int nHits = 4;
    // int arr_local_iEta[nHits] = {0, 1, 8, 9};
    // int arr_local_iPhi[nHits] = {0, 0, 10, 10};
    
    //    for (int j = 0; j < nHits; j++) {
    for (const auto& hit : ecalhits) {

      // Debugging only: do not use this block
      //  if (true) {

      // Check if the hit is in cards 0-35
      if ((getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc)) &&
      	  (getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc)) &&
      	  (getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc)) &&
      	  (getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc))) {
	
	// Get the crystal eta and phi, relative to the bottom left corner of the card 
	// (0 up to 17*5, 0 up to 4*5) 
	int local_iEta = getCrystal_local_iEta(hit.position().eta(), cc);
        int local_iPhi = getCrystal_local_iPhi(hit.position().phi(), cc);
      
	// int local_iEta = arr_local_iEta[j];
	// int local_iPhi = arr_local_iPhi[j];
	
	// Once we have the iEta and iPhi of the crystal relative to the bottom left corner,
	// everything else is determined.
	
	// Figure out what region (0-5) the hit falls into 
	// The region number (0-5) depends only on the local crystal iEta
	int regionNumber = getRegionNumber(local_iEta);
	
	// Get the tower eta and phi index inside the card (17x4)
	int inCard_tower_iEta = int(local_iEta / CRYSTALS_IN_TOWER_ETA); 
	int inCard_tower_iPhi = int(local_iPhi / CRYSTALS_IN_TOWER_PHI);
	
	// Get the tower eta and phi index inside the region (3x4)
	int inRegion_tower_iEta = inCard_tower_iEta % TOWER_IN_ETA;
	int inRegion_tower_iPhi = inCard_tower_iPhi % TOWER_IN_PHI;
	
	// Within the region, figure out which crystal and link the hit falls into                                         
        // Get the crystal eta and phi index inside the 3x4 region (15x20)
	int inRegion_crystal_iEta = local_iEta % (TOWER_IN_ETA * CRYSTALS_IN_TOWER_ETA);
	int inRegion_crystal_iPhi = local_iPhi;

	// Get the crystal eta and phi index inside the link (5x5). E.g. crystal (7, 9) is crystal number (2, 4) in its link
	int inLink_crystal_iEta = (inRegion_crystal_iEta % CRYSTALS_IN_TOWER_ETA);
	int inLink_crystal_iPhi = (inRegion_crystal_iPhi % CRYSTALS_IN_TOWER_PHI);

	
	// std::cout << "local_iEta, local_iPhi, regionNumber: "
	// 	  << local_iEta << ", " 
	// 	  << local_iPhi << ", "
	// 	    << regionNumber << std::endl;
	// std::cout << "inRegion_crystal_iEta, inRegion_crystal_iPhi (expecting 15x20), regionNumber (expecting 0-5): " << std::endl;
	// << inRegion_crystal_iEta << ", "
	// << inRegion_crystal_iPhi << ", "
	// << regionNumber << std::endl;
	  
	// Access the right region -> link -> crystal and increment the energy
	if (regionNumber < N_REGIONS_PER_CARD) {
	  region3x4& myRegion = rctCard.getRegion3x4(regionNumber);

	  //	  std::cout << "inRegion_tower_iEta, inRegion_tower_iPhi (expecting 3x4): " << inRegion_tower_iEta << ", "
	  //		    << inRegion_tower_iPhi << std::endl;
	  
	  // Get the right link
	  linkECAL& myLink = myRegion.getLinkECAL(inRegion_tower_iEta, inRegion_tower_iPhi);
	  
	  // Add the energy to the right crystal 
	  //	  std::cout << "inLink_crystal_iEta, inLink_crystal_iPhi (expecting 5x5): " 
	  //		    << inLink_crystal_iEta << ", " << inLink_crystal_iPhi << std::endl;

	  myLink.addCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi, hit.et_uint());

	  // TEMP: add a dummy energy
	  // myLink.addCrystalE(inLink_crystal_iEta, inLink_crystal_iPhi, 10);
	  
	}

	if (cc == 28 && regionNumber == 3) {
	  std::cout << "Card: " << cc 
		    << ", hit (Eta, phi, et): "
		    << hit.position().eta() << ", " << hit.position().phi() << ", " << hit.et_uint() << ". "
		    << "(iEta, iPhi): " 
		    << getCrystal_iEta(hit.position().eta()) << ", "
		    << getCrystal_iPhi(hit.position().phi()) << ". "
		    << "(Crystal eta, phi): "
		    << getEta_fromCrystaliEta(getCrystal_iEta(hit.position().eta())) << ", "
		    << getPhi_fromCrystaliPhi(getCrystal_iPhi(hit.position().phi())) << ". "
		    // << "(Go back to iEta, iPhi in card: ) "
		    // << getCrystal_iEta_fromCardRegionInfo(cc, regionNumber, inRegion_tower_iEta, inLink_crystal_iEta) << ", "
		    // << getCrystal_iPhi_fromCardRegionInfo(cc, regionNumber, inRegion_tower_iPhi, inLink_crystal_iPhi) << "."
		  // << "local_iEta/iPhi: " << local_iEta << ", " << local_iPhi << ", " 
		  // << "inCard_tower_iEta/iPhi: " << inCard_tower_iEta     << ", " << inCard_tower_iPhi << ", "
		  // << "region/inRegion_tower_iEta/iPhi: " << regionNumber << ", " << inRegion_tower_iEta << ", "
		  // << inRegion_tower_iPhi << ", "
		  // << "inLink crystal iEta/iPhi: " << inLink_crystal_iEta << ", " << inLink_crystal_iPhi 
		    << std::endl;

	}
      }      
  } // end of loop over ECAL hits

    // Also initialize tower (iEta, iPhi) coordinates (code lifted from 
    // https://github.com/cms-l1t-offline/cmssw/blob/25a1610b718c4cf94c33afb6e23767b5d3a677d7/L1Trigger/L1CaloTrigger/plugins/L1EGammaCrystalsEmulatorProducer.cc#L717-L727), changing the const variable names and getTowers_absEtaID function name
    // n.b. iEta_tower_L1Card and iPhi_tower_L1Card (the CMSSW emulator Layer 1 outputs) use a different convention
    // so later on when we compute a 17x4 tower array, it is rotated for the negative eta cards.
    static constexpr float tower_width = 0.0873;
    for (int jj = 0; jj < n_links_card ; ++jj) {
      for (int ii = 0; ii < n_towers_per_link; ++ii) {
	float phi = getPhiMin_card_emulator(cc) * tower_width / CRYSTALS_IN_TOWER_PHI - M_PI + (jj + 0.5) * tower_width;
	float eta = getEtaMin_card_emulator(cc) * tower_width / CRYSTALS_IN_TOWER_ETA - n_towers_cardEta * tower_width +
	  (ii + 0.5) * tower_width;
	iEta_tower_L1Card[jj][ii][cc] = getTower_absEtaID(eta);
	iPhi_tower_L1Card[jj][ii][cc] = getTower_absPhiID(phi);
	// std::cout << "Real (eta, phi): " << eta << ", " << phi << "; " 
	// 	  << "iEta_tower_L1Card[" << jj << "][" << ii 
	// 	  << "][" << cc << "] = " << getTower_absEtaID(eta) << ", "
	// 	  << "iPhi_tower_L1Card[" << jj << "][" << ii
	// 	  << "][" << cc << "] = " << getTower_absPhiID(phi) << std::endl;
      }
    }
    

    //***************************************************// 
    // Build RCT towers from HCAL
    //***************************************************// 

    // Same idea as the ECAL RCT geometry, except we only care about the ET in towers 

    // Loop over hcal hits to get the HCAL towers.

    // TEMP: don't do HCAL hits
    //    for (int k = 0; k < 1; k++) {
    for (const auto& hit : hcalhits) {

      // TEMP: do not use this conditional
      //      if (true) {
      if (getCrystal_iPhi(hit.position().phi()) <= getCard_iPhiMax(cc) &&
          getCrystal_iPhi(hit.position().phi()) >= getCard_iPhiMin(cc) &&
          getCrystal_iEta(hit.position().eta()) <= getCard_iEtaMax(cc) &&
          getCrystal_iEta(hit.position().eta()) >= getCard_iEtaMin(cc) && hit.pt() > 0) {

      	// HCAL: Get the crystal eta and phi, relative to the bottom left corner of the card 
      	// (0 up to 17*5, 0 up to 4*5) 
      	int local_iEta = getCrystal_local_iEta(hit.position().eta(), cc);
      	int local_iPhi = getCrystal_local_iPhi(hit.position().phi(), cc);

	// int local_iEta = 0;
	// int local_iPhi = 0;
	
	// HCAL: Once we have the iEta and iPhi of the crystal relative to the bottom left corner,
	// everything else is determined.
	
	// HCAL: Figure out what region (0-5) the hit falls into 
	// The region number (0-5) depends only on the local crystal iEta
	int regionNumber = getRegionNumber(local_iEta);
	
	// HCAL: Get the tower eta and phi index inside the card (17x4)
	int inCard_tower_iEta = int(local_iEta / CRYSTALS_IN_TOWER_ETA); 
	int inCard_tower_iPhi = int(local_iPhi / CRYSTALS_IN_TOWER_PHI);
	
	// HCAL: Get the tower eta and phi index inside the region (3x4)
	int inRegion_tower_iEta = inCard_tower_iEta % TOWER_IN_ETA;
	int inRegion_tower_iPhi = inCard_tower_iPhi % TOWER_IN_PHI;

	// std::cout << "HCAL hit eta/phi : "  << hit.position().eta() << ", " << hit.position().phi() << ", region: " << regionNumber << ", "
	// 	  << "inRegion_tower_iEta and iPhi : " << inRegion_tower_iEta << ", " << inRegion_tower_iPhi 
	// 	  << ", energy " << hit.et_uint() 
	// 	  << std::endl;

	// Access the right HCAL region -> tower and increment the ET
	if (regionNumber < N_REGIONS_PER_CARD) {
	  towers3x4& myTowers3x4 = rctCard.getTowers3x4(regionNumber);
	  towerHCAL& myTower = myTowers3x4.getTowerHCAL(inRegion_tower_iEta, inRegion_tower_iPhi);

	  // TEMP: add a dummy energy
	  // myTower.addEt(1);
	  myTower.addEt(hit.et_uint());

	}
      }
    } // end of loop over hcal hits

    //***************************************************// 
    // In each ECAL region, perform clustering
    //***************************************************// 

    // Dummy for testing
    std::vector<Cluster> test_cluster;

    for (int idxRegion = 0; idxRegion < N_REGIONS_PER_CARD; idxRegion++) {
      
      crystal temporary[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI];       // ECAL crystal array (will be changed)
      ap_uint<12> towerEtHCAL[TOWER_IN_ETA * TOWER_IN_PHI];    // HCAL tower ET in the 3x4 region

      if (cc > -1) {
	region3x4& myRegion = rctCard.getRegion3x4(idxRegion);
	towers3x4& myTowers = rctCard.getTowers3x4(idxRegion);

	//	std::cout << std::endl << "[----] DOING CARD " << cc << " AND REGION IDX " << myRegion.getIdx() << std::endl;

	// In each 3x4 region, loop through the links (one link per tower)
	for (int iLinkEta = 0; iLinkEta < TOWER_IN_ETA; iLinkEta++) {
	  for (int iLinkPhi = 0; iLinkPhi < TOWER_IN_PHI; iLinkPhi++) {
	    
	    // Get the ECAL link (one link per tower) 
	    linkECAL& myLink = myRegion.getLinkECAL(iLinkEta, iLinkPhi);
	    
	    // Each link has a different bottom left corner in the (iEta, iPhi) of the ECAL region.
	    // This will help us fill the ECAL region 'temporary' array, going link by link 
	    int ref_iEta = (iLinkEta * CRYSTALS_IN_TOWER_ETA); 
	    int ref_iPhi = (iLinkPhi * CRYSTALS_IN_TOWER_PHI); 

	    // In the link, get the crystals (5x5 in each link) 
	    for (int iEta = 0; iEta < CRYSTALS_IN_TOWER_ETA; iEta++) {           	    
	      for (int iPhi = 0; iPhi < CRYSTALS_IN_TOWER_PHI; iPhi++) {    
		
		// Et as unsigned int
		ap_uint<10> uEnergy = myLink.getCrystalE(iEta, iPhi);

		// if (uEnergy > 0) {
		//   std::cout << "energy>0: Accessing temporary array: " << (ref_iEta + iEta) << ", " << (ref_iPhi + iPhi) 
		// 	    << ", writing energy (uint:) " << uEnergy << std::endl;
		// }

		// Fill the 'temporary' array with a crystal object 
		temporary[ref_iEta + iEta][ref_iPhi + iPhi] = crystal(uEnergy);
	      }
	    } // end of loop over crystals

	    ///// Get HCAL tower ET
	    towerHCAL& myTower = myTowers.getTowerHCAL(iLinkEta, iLinkPhi);
	    //  std::cout << "towerEtHCAL: " << myTower.getEt() << ", ";
	    towerEtHCAL[(iLinkEta * TOWER_IN_PHI) + iLinkPhi] = myTower.getEt();
	  }
	}
	
	// Iteratively find four clusters and remove them from 'temporary' as we go, and fill cluster_list
	for (int c = 0; c < N_CLUSTERS_PER_REGION; c++) {
	  Cluster newCluster = getClusterFromRegion3x4(temporary); // remove cluster from 'temporary' 
	  newCluster.setRegionIdx(idxRegion);                      // add the region number 
	  if (newCluster.clusterEnergy() > 0) {
	    cluster_list[cc].push_back(newCluster);                // do not push back 0-energy clusters
	  }
	}
	
	// Create towers using remaining ECAL energy, and the HCAL towers were already calculated in towersEtHCAL[12] 
	ap_uint<12> towerEtECAL[12];
	getECALTowersEt(temporary, towerEtECAL);
	
	// Fill in towerHCALCard and towerECALCard arrays (for whole detector)
	for (int i = 0; i < 12; i++) {
	  // Get the tower's indices in a (17x4) card
	  int iEta = (idxRegion * TOWER_IN_ETA) + (i / TOWER_IN_PHI);   
	  int iPhi = (i % TOWER_IN_PHI);
	  //	  std::cout << "(" << iEta << "," << iPhi << ")";
	  // towerHCALCard[iEta][iPhi][cc] = tower_t(towerEtHCAL[i], 0, 0);
	  // towerECALCard[iEta][iPhi][cc] = tower_t(towerEtECAL[i], 0, 0);
	  towerHCALCard[iEta][iPhi][cc] = tower_t(towerEtHCAL[i], 0);
	  towerECALCard[iEta][iPhi][cc] = tower_t(towerEtECAL[i], 0);
	}

      } // end of "if"

    } // end of loop over regions

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

       std::cout << "Card " << cc << ": BEFORE stitching: " << std::endl;
    
    for (unsigned int kk = 0; kk < cluster_list[cc].size(); ++kk) {
      Cluster c = cluster_list[cc][kk];
      std::cout << cluster_list[cc][kk].clusterEnergy()
    		<< " (" << cluster_list[cc][kk].towerEtaInCard()
    		<< ", " << cluster_list[cc][kk].towerPhi() << ") "
		<< " (" << cluster_list[cc][kk].clusterEta() 
		<< ", " << cluster_list[cc][kk].clusterPhi() << ") ";
    }
    int nRegionBoundariesEta = (N_REGIONS_PER_CARD - 1); // 6 regions -> 5 boundaries to check
    int towerEtaBoundaries[nRegionBoundariesEta][2] = 
      {	{15, 14},
	{12, 11},
	{9, 8},
	{6, 5},
	{3, 2}      };

    for (int iBound = 0; iBound < nRegionBoundariesEta; iBound++) {
      // std::cout << "   Will check for stitching along " << towerEtaBoundaries[iBound][0]
      //           << " and " << towerEtaBoundaries[iBound][1] << std::endl;
      // First number in the tuple is the upper boundary, second number is the lower boundary
      stitchClusterOverRegionBoundary(cluster_list[cc], towerEtaBoundaries[iBound][0], towerEtaBoundaries[iBound][1]);
    }
    
    std::cout << "Card " << cc << ": AFTER stitching:" << std::endl;
    for (unsigned int kk = 0; kk < cluster_list[cc].size(); ++kk) {
      Cluster c = cluster_list[cc][kk];
      std::cout << cluster_list[cc][kk].clusterEnergy()
		<< " (" << cluster_list[cc][kk].towerEtaInCard()
		<< ", " << cluster_list[cc][kk].towerPhi() << ")" 
		<< " (" << cluster_list[cc][kk].clusterEta()
		<< ", " << cluster_list[cc][kk].clusterPhi() <<") ";
    }

    //-------------------------------------------//  
    // Sort the clusters:                        //
    // Take only the 8 highest ET clusters,      //
    // add leftover clusters back to ECAL towers //
    //-------------------------------------------//
    std::cout << "Card " << cc << ": unsorted ET: "; 
    
    std::cout << "Up to twelve clusters going in: printing ET and tower (17x4):";
    if (cluster_list[cc].empty()) {
      std::cout << "The vector is empty; skip this step" << std::endl;
    }
    else {
      //      std::cout << "There will be " << cluster_list[cc].size() << " elements to loop over" << std::endl;
      // for (unsigned int kk = 0; kk < cluster_list[cc].size(); ++kk) {
      // 	Cluster c = cluster_list[cc][kk];
      // 	std::cout << cluster_list[cc][kk].clusterEnergy() 
      // 		  << " (" << cluster_list[cc][kk].towerEtaInCard() 
      // 		  << ", " << cluster_list[cc][kk].towerPhi() << ") "; 
      
      std::cout << std::endl;
      std::sort(cluster_list[cc].begin(), cluster_list[cc].end(), compareClusterET);
      
      // If there are more than eight clusters, return the unused energy to the towers
      for (unsigned int kk = n_clusters_4link; kk < cluster_list[cc].size(); ++kk) {
	Cluster cExtra = cluster_list[cc][kk];
	if (cExtra.clusterEnergy() > 0) {
	  //	  std::cout << "Extra cluster # " << kk << ": energy (GeV) " << cExtra.clusterEnergy()
	  //		    << std::endl;
	  // Increment tower ET in towerECALCard[17][4] which uses RCT HLS geometry
	  // Get tower (eta, phi) (up to (17, 4)) in the RCT card
	  int whichTowerEtaInCard = ((cExtra.region() * TOWER_IN_ETA) + cExtra.towerEta());
	  int whichTowerPhiInCard = cExtra.towerPhi(); 
	  ap_uint<12> oldTowerEt = towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc].et();
	  ap_uint<12> newTowerEt = (oldTowerEt + cExtra.clusterEnergy());
	  ap_uint<4>  hoe        = towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc].hoe();
	  // ap_uint<4>  satur      = towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc].satur();
	  // towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc] = tower_t(newTowerEt, hoe, satur);
	  towerECALCard[whichTowerEtaInCard][whichTowerPhiInCard][cc] = tower_t(newTowerEt, hoe);
	  
	  // std::cout << "... Adding to card eta (0-17) " << whichTowerEtaInCard
	  // 	    << ", card phi (0-4) " << whichTowerPhiInCard
	  // 	    << " old energy : " << oldTowerEt
	  // 	    << " new energy : " << newTowerEt
	  // 	    << std::endl;
	}
      }
      
      // Build the sorted cluster list IFF there were clusters 
      std::cout << "Building sorted cluster list: ";
      if (cluster_list[cc].empty()) {
	//	std::cout << "No clusters: do not build sorted cluster list" << std::endl;
      }
      else {
	// Save up to eight clusters: loop over cluster_list
	for (unsigned int kk = 0; kk < cluster_list[cc].size(); ++kk) {
	  if (kk >= n_clusters_4link) continue;
	  if (cluster_list[cc][kk].clusterEnergy() > 0) 
	    cluster_list_merged[cc].push_back(cluster_list[cc][kk]);
	}
	std::cout << "Sorted, up to 8 clusters: ";
	for (unsigned int kk = 0; kk < cluster_list_merged[cc].size(); ++kk) {
	  std::cout << cluster_list_merged[cc][kk].clusterEnergy() 
		    << " (" << cluster_list_merged[cc][kk].towerEtaInCard() 
		    << ", " << cluster_list_merged[cc][kk].towerPhi() << ") ";
	}
	std::cout << std::endl;
      }
    }
    // END HERE for str-replace testing

    //-------------------------------------------//                          
    // Calibration:                                                                                    
    //   - Applied to the cluster_list_merged                                 
    //     (calib_(cluster_list[cc][jj].cpt, std::abs(cluster_list[cc][jj].craweta_)))     
    //   - Also applied to the ECAL hits that were not clustered, that go into ECAL_tower_L1Card
    //     (calib_(0, std::abs(hit.position().eta()))                           
    //-------------------------------------------//   
    
    // Calibrate clusters
    for (auto & c : cluster_list_merged[cc]) {
      float realEta = getEta_fromCrystaliEta( getCrystal_iEta_fromCardRegionInfo(cc, 
										 c.region(),
										 c.towerEta(),
										 c.clusterEta()));
      
      c.calib = calib_(c.getPt(), std::abs(realEta));
      std::cout << "Calib with pT " << c.getPt() << " and eta " << realEta << ": " 
                << c.getCalib() << std::endl;
      std::cout << ">>> Old pT: "
		<< c.getPt() << " (uint:" <<c.clusterEnergy() << ")" << std::endl;
      c.applyCalibration(c.calib);

      std::cout << ">>> New pT after c.applyCalibration: " 
		<< c.getPt() << " (uint:" << c.clusterEnergy() << ")" << std::endl;
      
    }

    // Calibrate ECAL tower ET
    // Access towerECALCard[17][4] (array of tower_t) which uses RCT HLS geometry
    for (int ii = 0; ii < n_towers_cardPhi; ++ii) { // 4 towers per card in phi
      for (int jj = 0; jj < n_towers_cardEta; ++jj) { // 17 towers per card in eta
	float tRealEta = getTowerEta_fromAbsID(iEta_tower_L1Card[ii][jj][cc]);  // real eta of center of tower
	double tCalib = calib_(0, tRealEta);                                    // calibration factor
	towerECALCard[jj][ii][cc].applyCalibration(tCalib);
	
      }
    }

    // Compute HoE for towers
    std::cout << "Computing HoE for towers: " << std::endl;
    for (int ii = 0; ii < n_towers_cardPhi; ++ii) { // 4 towers per card in phi                                                                            
      for (int jj = 0; jj < n_towers_cardEta; ++jj) { // 17 towers per card in eta   
	// std::cout << (int) towerECALCard[jj][ii][cc].hoe();
	ap_uint<12> ecalEt = towerECALCard[jj][ii][cc].et();
	ap_uint<12> hcalEt = towerHCALCard[jj][ii][cc].et();
	towerECALCard[jj][ii][cc].getHoverE(ecalEt, hcalEt); 
	// DUMMY TEST
	// towerECALCard[jj][ii][cc].getHoverE(rand() % 100 + 1, rand() % 100 + 1);
	// std::cout << " -> " << (int) towerECALCard[jj][ii][cc].hoe() << ", ";
      }
    }
    std::cout << std::endl;
					      
    //-------------------------------------------------------------------------------//
    // Write the L1 outputs in the style of previous CMSSW
    //-------------------------------------------------------------------------------//

    std::cout << "Sanity check: Card " << cc << ": SORTED ET: (if >0 clusters exist, print stuff)";
    
    // (May be deprecated:) If the cluster list isn't empty, write out the old CMSSW (by Cecile)
    // -style arrays energy_cluster_L1Card, crystalID_cluster_L1Card, and towerID_clusterL1Card
    // which treats (0, 0) as always the top left corner of the RCT card if you look at the
    // usual RCT diagram. (aka negative eta cards are NOT rotated)
    if (!cluster_list_merged[cc].empty()) {
      for (unsigned int jj = 0; jj < cluster_list_merged[cc].size(); ++jj) {
	
	Cluster c = cluster_list_merged[cc][jj];
	// std::cout << c.clusterEnergy() << " (" << c.clusterEta() << ", " << c.clusterPhi() << ") "
	// 	  << ", setting crystalID_cluster_L1Card[" << jj % n_links_card << "]["
	// 	  << jj / n_links_card << "][" << cc << "] = " << getCrystalIDInTower_emulator(c.clusterEta(), c.clusterPhi(), cc) 
	// 	  << ", setting towerID_cluster_L1Card[" << jj % n_links_card << "]["
	// 	  << jj / n_links_card << "][" << cc << "] = " << getTowerID_emulator(c.towerEta(), c.towerPhi(), cc, c.region())
	// 	  << std::endl;
	
	energy_cluster_L1Card[jj % n_links_card][jj / n_links_card][cc] =
	  c.clusterEnergy();
	crystalID_cluster_L1Card[jj % n_links_card][jj / n_links_card][cc] = 
	  getCrystalIDInTower_emulator(c.clusterEta(), c.clusterPhi(), cc);
	towerID_cluster_L1Card[jj % n_links_card][jj / n_links_card][cc] =
	  getTowerID_emulator(c.towerEta(), c.towerPhi(), cc, c.region());
      }
    }
    
    std::cout << std::endl;

    // (May be deprecated:) Fill out ECAL_tower_L1Card and HCAL_tower_L1Card, which is in the style of the old CMSSW (by Cecile). The firmware 17x4 array treats the "bottom left" corner of the card as (0, 0)
    // (rotating the negative eta cards so that the endcap region is pointing up)
    // while the old CMSSW (by Cecile) treats the 4x17 array as always starting in the top left corner if we look
    // at the usual RCT diagram.
    //    std::cout << "Re-package towerECALCard into old CMSSW emulator geometry" << std::endl;
    for (int i = 0; i < n_towers_cardEta; i++) {
      for (int j = 0; j < n_towers_cardPhi; j++ ) {
        
    	// n.b. L1 output is 4*17*36, hence [j][i] instead of [i][j] on the L.H.S.    
    	if ((cc % 2) == 1) { // if cc is odd (positive eta)
    	  ECAL_tower_L1Card[j][i][cc] = towerECALCard[i][j][cc].et();
    	  HCAL_tower_L1Card[j][i][cc] = towerHCALCard[i][j][cc].et();
    	}
    	else {  // if cc is even (negative eta), we need to rotate the coordinates
    	  //	  std::cout << "writing to [" << (3-j) << "][" << (16-i) << "] ";
    	  ECAL_tower_L1Card[3-j][16-i][cc] = towerECALCard[i][j][cc].et();
    	  HCAL_tower_L1Card[3-j][16-i][cc] = towerHCALCard[i][j][cc].et();
    	}
      }
    }
    std::cout << std::endl;


    //-----------------------------------------------------------//
    // Produce output collections for event display and analyzer
    //-----------------------------------------------------------//
    for (auto & c : cluster_list_merged[cc]) {                                                      
      std::cout << c.clusterEnergy() 
		<< " c.region: " << c.region() << "," 
		<< " c.towerEta and Phi:  (" << c.towerEta() << ", " << c.towerPhi() << "), "
		<< " c.clusterEta and Phi: (" << c.clusterEta() << ", " << c.clusterPhi() << ") ";    

      float realEta = getEta_fromCrystaliEta( getCrystal_iEta_fromCardRegionInfo(cc, c.region(), c.towerEta(), c.clusterEta()) );
      float realPhi = getPhi_fromCrystaliPhi( getCrystal_iPhi_fromCardRegionInfo(cc, c.region(), c.towerPhi(), c.clusterPhi()) );
      
      std::cout << "Real eta, phi " << " (" << realEta << ", " << realPhi << ") ";
      c.calib = calib_(c.clusterEnergy()/8.0, std::abs(realEta));
      std::cout << "Calib with pT " << c.getPt() << " and eta " << realEta << ": "
		<< c.getCalib() << std::endl;
      reco::Candidate::PolarLorentzVector p4calibrated(c.getPt(), // use float
						       realEta,
						       realPhi,
						       0.);
      
      // Constructor definition at: https://github.com/cms-l1t-offline/cmssw/blob/l1t-phase2-v3.3.11/DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h#L34
      // showerShape_cluster_L1Card: depends on passes_ss(cpt, 2x5/5x5)
      // equivalent to cshowershape_. equivalent to showerShape_cluster_L2Card
      // is_ss
      bool is_ss = passes_ss(c.getPt(),
			     (c.getEt2x5() / c.getEt5x5()) );
      c.is_ss        = is_ss;
      
      // showerShapeLooseTk_cluster_L1Card: depends on passes_looseTkss(cpt, 2x5/5x5)
      // equivalent to cphotonshowershape_. equivalent to photonShowerShape_cluster_L2Card
      // is_looseTkss
      bool is_looseTkss = passes_looseTkss(c.getPt(), 
					   (c.getEt2x5() / c.getEt5x5() ));  // calibration cancels
      //      std::cout << "is_ss: " << is_ss << ", is_looseTkss: " << is_looseTkss << std::endl;
      c.is_looseTkss = is_looseTkss;

      // std::cout << "is_ss: " << is_ss << ", is_looseTkss: " << is_looseTkss << std::endl;
      // std::cout << "c.getIsSS(): " << c.getIsSS() << ", c.getIsLooseTkss(): " << c.getIsLooseTkss() << std::endl;
      
      // Using some dummy/stand-in variables for the flags, see params[] for flags instead
      // isolation is not computed at this stage 
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

      // RCT flags
      std::map<std::string, float> params;
      params["standaloneWP_showerShape"] = is_ss;
      params["trkMatchWP_showerShape"]   = is_looseTkss;
      cluster.setExperimentalParams(params);
      
      L1EGXtalClusters->push_back(cluster);
    }  // end of loop over clusters                                                                                                           
    // Produce output tower collections
    for (int ii = 0; ii < n_towers_cardPhi; ++ii) { // 4 towers per card in phi
      for (int jj = 0; jj < n_towers_cardEta; ++jj) { // 17 towers per card in eta

	l1tp2::CaloTower l1CaloTower;
	// Divide by 8.0 to get ET as float (GeV)
	l1CaloTower.setEcalTowerEt(ECAL_tower_L1Card[ii][jj][cc]/8.0);
	// HCAL TPGs encoded ET: multiply by the LSB (0.5) to convert to GeV
	float hcalLSB = 0.5;
	l1CaloTower.setHcalTowerEt(HCAL_tower_L1Card[ii][jj][cc] * hcalLSB);
	l1CaloTower.setTowerIEta(getToweriEta_fromAbsID(iEta_tower_L1Card[ii][jj][cc]));
        l1CaloTower.setTowerIPhi(getToweriPhi_fromAbsID(iPhi_tower_L1Card[ii][jj][cc]));
        l1CaloTower.setTowerEta(getTowerEta_fromAbsID(iEta_tower_L1Card[ii][jj][cc]));
        l1CaloTower.setTowerPhi(getTowerPhi_fromAbsID(iPhi_tower_L1Card[ii][jj][cc]));
	
	L1CaloTowers->push_back(l1CaloTower);
      }
    }


  } // end of loop over cards

  iEvent.put(std::move(L1EGXtalClusters), "RCT");
  iEvent.put(std::move(L1CaloTowers),     "RCT");
  std::cout << "Finished producing RCT emulator for this event" << std::endl;
  
  //*******************************************************************
  //*************** Do GCT geometry for inputs ************************
  //*******************************************************************

  // Loop over GCT cards (three of them)
  GCTcard_t gctCards[N_GCTCARDS];
  GCTtoCorr_t gctToCorr[N_GCTCARDS];

  for (unsigned int gcc = 0; gcc < N_GCTCARDS; gcc++) {
    // Each GCT card encompasses 16 RCT cards, listed in 
    // GCTcardtoRCTcardnumber[3][16]. 
    //    std::cout << "GCT: Starting Card " << gcc << "..." << std::endl;
    // i goes from 0 to <16
    for (int i = 0; i < (N_RCTCARDS_PHI * 2); i++) {

      unsigned int rcc = GCTcardtoRCTcardnumber[gcc][i];
      //      std::cout << "... Starting RCT Card: " << rcc << "..." << std::endl;

      // Positive eta? Fist row is in positive eta
      bool isPositiveEta = (i < N_RCTCARDS_PHI);
      
      for (int iLink = 0; iLink < n_links_card; iLink++) {

	// Get the towers (to-do: H/E has been calculated at this point already in the towerECALCard array)
	// Add the tower ECAL and tower HCAL energies!
	// std::cout << "Adding ECAL and HCAL towers.. ";
	for (int iTower = 0; iTower < N_GCTTOWERS_FIBER; iTower++) {
	  tower_t t0_ecal = towerECALCard[iTower][iLink][rcc];
	  tower_t t0_hcal = towerHCALCard[iTower][iLink][rcc];
	  RCTtower_t t;
	  t.et = t0_ecal.et() + convertHcalETtoEcalET(t0_hcal.et());
	  t.hoe = t0_ecal.hoe();
	  // Not needed for GCT firmware but will be written into GCT CMSSW outputs : 12 bits each 
	  t.ecalEt = t0_ecal.et();
	  t.hcalEt = t0_hcal.et(); 
	  // std::cout << "tower et and hoe: " << t.et << ", " << t.hoe;

	  if (isPositiveEta) {  
	    gctCards[gcc].RCTcardEtaPos[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower] = t;
	  } else {
	    gctCards[gcc].RCTcardEtaNeg[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower] = t;
	  }
	} // end of loop over 17 towers in one link
      } // end of loop over 4 links in one card
      
      // Get the RCT clusters, and distribute them across four links, ensuring that all values are
      // correctly converted to their GCT equivalents.
      // There are at most eight clusters per card
      for (size_t iCluster = 0; 
	   (iCluster < cluster_list_merged[rcc].size()) && 
	     (iCluster < (N_RCTGCT_FIBERS * N_RCTCLUSTERS_FIBER));
	   iCluster++) {
	Cluster c0 = cluster_list_merged[rcc][iCluster];
	RCTcluster_t c;
	c.et     = c0.clusterEnergy();
	
	c.towEta = (c0.region() * TOWER_IN_ETA) + c0.towerEta();     
	// towerEta is unusual because the class 'Cluster' stores the region number (region())
	// and towerEta refers to the tower iEta INSIDE the region, but RCTTcluster_t stores no
	// information about the region number and assumes tower iEta goes from 0 to 17.
	c.towPhi = c0.towerPhi();
	c.crEta  = c0.clusterEta();
	c.crPhi  = c0.clusterPhi();
	c.et5x5  = c0.uint_et5x5();  // newly added
	c.et2x5  = c0.uint_et2x5();  // newly added
	c.is_ss        = c0.getIsSS();        // newly added
	c.is_looseTkss = c0.getIsLooseTkss(); // newly added
	
	unsigned int iIdxInGCT = i % N_RCTCARDS_PHI;  
	unsigned int iLinkC = iCluster % N_RCTGCT_FIBERS;
	unsigned int iPosC  = iCluster / N_RCTGCT_FIBERS;
	std::cout << c.et << ", "
		  << "accessing link " << iCluster % N_RCTGCT_FIBERS << " "
		  << "and position "   << iCluster / N_RCTGCT_FIBERS << " "
		  << "isPositiveEta "  << isPositiveEta << " " 
		  << "with iIdxInGCT " << iIdxInGCT << " " 
		  << "with uint et5x5, et2x5 " << c.et5x5 << ", " << c.et2x5 << " " 
		  << "with shower shape flags ss/looseTkss " << c.is_ss << ", " << c.is_looseTkss << std::endl;
	
	if (isPositiveEta) {
	  gctCards[gcc].RCTcardEtaPos[iIdxInGCT].RCTtoGCTfiber[iLinkC].RCTclusters[iPosC] = c;
	} else {
	  gctCards[gcc].RCTcardEtaNeg[iIdxInGCT].RCTtoGCTfiber[iLinkC].RCTclusters[iPosC] = c;
	}
      }
      // If there were fewer than eight clusters, make sure the remaining fiber clusters are zero'd out.
      for (size_t iZeroCluster = cluster_list_merged[rcc].size();
	   iZeroCluster < (N_RCTGCT_FIBERS * N_RCTCLUSTERS_FIBER);
	   iZeroCluster++) {
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
        } else {
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

      for (int iLink = 0; iLink < n_links_card; iLink++) {
	for (int iTower = 0; iTower < N_GCTTOWERS_FIBER; iTower++) {
	  //  std::cout << gctCards[gcc].RCTcardEtaPos[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower].et << ", ";
	}
      }
      // std::cout << std::endl;
      //      std::cout << "Clusters: ";
      for (int iLink = 0; iLink < n_links_card; iLink++) {
	for (int iCluster = 0; iCluster < N_RCTCLUSTERS_FIBER; iCluster++) {
	  // std::cout << gctCards[gcc].RCTcardEtaPos[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTclusters[iCluster].et << ", ";
	}
      }
    }
    std::cout << std::endl;
    for (int i = N_RCTCARDS_PHI; i < (N_RCTCARDS_PHI * 2); i++) {
      //      std::cout<< "Inside GCT: card (out of 8 in the negative side) " << i % N_RCTCARDS_PHI << std::endl;
	for (int iLink = 0; iLink < n_links_card; iLink++) {
	  for (int iTower = 0; iTower < N_GCTTOWERS_FIBER; iTower++) {
	    //    std::cout << gctCards[gcc].RCTcardEtaNeg[i % N_RCTCARDS_PHI].RCTtoGCTfiber[iLink].RCTtowers[iTower].et << ", ";
	  }
	}
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
    algo_top(gctCards[gcc], gctToCorr[gcc], gcc, L1GCTClusters, L1GCTTowers, L1GCTFullTowers);
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

  std::cout << "Checking that the GCT Full Towers have entries..." << std::endl;
  for (auto const& c: *L1GCTFullTowers ) {
    
    std::cout << "GCT FULL tower energy " << c.ecalTowerEt() 
	            << " at real eta, phi: (" << c.towerEta() << ", " << c.towerPhi() << ")" << std::endl;
    
  }
  
  
  iEvent.put(std::move(L1GCTClusters), "GCT");
  iEvent.put(std::move(L1GCTTowers),   "GCT");
  iEvent.put(std::move(L1GCTFullTowers), "GCTFullTowers");
} 



////////////////////////////////////////////////////////////////////////// 

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1CaloEGammaEmulator);

