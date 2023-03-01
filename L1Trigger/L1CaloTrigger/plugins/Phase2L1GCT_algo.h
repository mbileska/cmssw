#include "Phase2L1RCT.h"
#include "Phase2L1GCT.h"

#ifndef _PHASE_2_L1_GCT_ALGO_H_
#define _PHASE_2_L1_GCT_ALGO_H_

/* 
 * Get GCT cluster c's iEta (global iEta convention).
 * Use with getEta_fromCrystaliEta from Phase2L1RCT.h to convert from GCT cluster iEta to real eta.
 */
int getCluster_global_cluster_iEta(unsigned int nGCTCard, GCTcluster_t c) {

  // First get the "iEta/iPhi" in the GCT card. i.e. in the diagram where the barrel
  // is split up into three GCT cards, (iEta, iPhi) = (0, 0) is the top left corner
  // of the GCT card. 
  int iEta_in_gctCard;

  if (c.towEtaNeg) {
    // Negative eta: c.towEta and c.crEta count outwards from the real eta = 0 center line, so to convert to the barrel diagram global iEta
    // (global iEta = 0 from LHS of page), do (17*5 - 1) minus the GCT value.
    // e.g. If in GCT, a negative card's cluster had iEta = 84, this would be global iEta = 0.
    iEta_in_gctCard = ((N_GCTTOWERS_FIBER * CRYSTALS_IN_TOWER_ETA - 1) - ((c.towEta * CRYSTALS_IN_TOWER_ETA) + c.crEta));
  }
  else {
    // c.towEta and c.crEta count outwards from the real eta = 0 center line, so for positive
    // eta we need to add the 17*5 offset so that positive eta 0+epsilon starts at 17*5.
    // e.g. If in GCT, a positive card's cluster had iEta = 0, this would be global iEta = 85.
    // e.g. If in GCT, a positive card's cluster had iEta = 84, this would be global iEta = 169.
    iEta_in_gctCard = ((N_GCTTOWERS_FIBER * CRYSTALS_IN_TOWER_ETA) + ((c.towEta * CRYSTALS_IN_TOWER_ETA) + c.crEta));
  }

  // Last, convert to the global iEta/iPhi in the barrel region. For eta there is nothing to
  // change.
  int iEta_in_barrel = iEta_in_gctCard;
  
  return iEta_in_barrel;
}

/* 
 * Get GCT cluster c's iPhi (global convention).
 * Use with getPhi_fromCrystaliPhi from Phase2L1RCT.h to convert from GCT cluster to real phi.
 * If returnGlobalGCTiPhi is true (Default value) then return the iPhi in the entire GCT barrel. Otherwise
 * just return the iPhi in the current GCT card.
 */
int getCluster_global_cluster_iPhi(unsigned int nGCTCard, GCTcluster_t c, bool returnGlobalGCTiPhi = true) {

  assert(nGCTCard <= 2); 

  int iPhi_in_gctCard = ((c.towPhi * CRYSTALS_IN_TOWER_PHI) + c.crPhi);

  // If we should return the global GCT iPhi, get the iPhi offset due to the number of the GCT card
  int iPhi_card_offset = 0;
  if (returnGlobalGCTiPhi) {
    if      (nGCTCard == 0) iPhi_card_offset = GCTCARD_0_TOWER_IPHI_OFFSET * CRYSTALS_IN_TOWER_PHI; 
    else if (nGCTCard == 1) iPhi_card_offset = GCTCARD_1_TOWER_IPHI_OFFSET * CRYSTALS_IN_TOWER_PHI;   
    else if (nGCTCard == 2) iPhi_card_offset = GCTCARD_2_TOWER_IPHI_OFFSET * CRYSTALS_IN_TOWER_PHI;
  }

  // Detector wraps around in phi: modulo number of crystals in phi (n_towers_Phi = 72)
  int iPhi_in_barrel = (iPhi_card_offset + iPhi_in_gctCard) % (n_towers_Phi * CRYSTALS_IN_TOWER_PHI);

  return iPhi_in_barrel;
}


/*
 * Each crystal falls in a tower: get this tower iEta in the GCT card (same as global) given the cluster's info. 
 */
 int getCluster_global_tower_iEta(unsigned int nGCTCard, GCTcluster_t c) {
  int crystaliEta_in_GCT_card = getCluster_global_cluster_iEta(nGCTCard, c);
  return (int) (crystaliEta_in_GCT_card / 5);
 }

 /*
  * Each crystal falls in a tower: get this tower iPhi in the GCT card (same as global) given the cluster's info. 
  */ 
int getCluster_global_tower_iPhi(unsigned int nGCTCard, GCTcluster_t c, bool returnGlobalGCTiPhi = true) {
  int crystaliPhi_in_GCT_card = getCluster_global_cluster_iPhi(nGCTCard, c, returnGlobalGCTiPhi);
  return (int) (crystaliPhi_in_GCT_card / 5);
}


/* 
 * Correlator fiber convention -> Global GCT convention
 * Get tower's global (iEta) from the GCTCorrFiber index [0, 64) and the tower's postion in the fiber [0, 17).
 * Recall that GCTCorrFiber is [0, 32) for negative eta and [32, 64) for positive eta. The tower's position in the fiber [0, 17)
 * always counts outwards from real eta = 0.
 * Use in conjunction with (float) getTowerEta_fromAbsID(int id) from Phase2L1RCT.h to get a tower's real eta.
 */ 
int getTower_global_toweriEta(unsigned int nGCTCard, unsigned int gctCorrFiberIdx, unsigned int posInFiber) {

  (void) nGCTCard; // not needed

  bool isTowerInPositiveEta = (gctCorrFiberIdx < N_GCTPOSITIVE_FIBERS); // N_GCTPOSITIVE_FIBERS = 32
  
  int global_toweriEta; 
  if (isTowerInPositiveEta) {  
    // e.g. For positive eta, posInFiber = 0 is at real eta = 0, so global tower iEta is 0 + 17 = 17
    global_toweriEta = (N_GCTTOWERS_FIBER + posInFiber); // N_GCTTOWERS_FIBER = 17
  } 
  else {
    // e.g. For negative eta, posInFiber = 0 is at real eta = 0, and global tower iEta is 17 - 1 - 0 = 16
    // posInFiber = 16 is at real eta = -1.4841, and global tower iEta is 17 - 1 - 16 = 0.
    global_toweriEta = (N_GCTTOWERS_FIBER - 1 - posInFiber); 
  }
  return global_toweriEta;
}

/* 
 * Correlator fiber convention -> Global GCT convention
 * Get tower's global (iPhi) from the GCT card number (0, 1, 2), and the GCTCorrFiber index [0, 64).
 * GCTCorrFiber is [0, 32) for negative eta and [32, 64) for positive eta. In the phi direction, fiber index #0 has the same phi 
 * as fiber index #32, so only the (fiber index modulo 32) matters for the phi direction. 
 * The tower's position in the fiber doesn't matter; in each fiber the phi is the same. 
 * Use in conjunction with (float) getTowerPhi_fromAbsID(int id) from Phase2L1RCT.h to get a tower's real phi.
 */
int getTower_global_toweriPhi(unsigned int nGCTCard, unsigned int gctCorrFiberIdx, unsigned int posInFiber) {

  (void) posInFiber; // not needed

  int global_tower_iPhi;
  
  unsigned int effectiveFiberIdx = (gctCorrFiberIdx % N_GCTPOSITIVE_FIBERS);  // N_GCTPOSITIVE_FIBERS = 32

  assert(nGCTCard <= 2);  // Make sure the card number is valid
  int toweriPhi_card_offset = 0;   
  if      (nGCTCard == 0) toweriPhi_card_offset = GCTCARD_0_TOWER_IPHI_OFFSET;
  else if (nGCTCard == 1) toweriPhi_card_offset = GCTCARD_1_TOWER_IPHI_OFFSET;
  else if (nGCTCard == 2) toweriPhi_card_offset = GCTCARD_2_TOWER_IPHI_OFFSET;

  global_tower_iPhi = (toweriPhi_card_offset + effectiveFiberIdx) % (n_towers_Phi);  //  as explained above, effectiveFiberIdx is [0, 32). n_towers_Phi = 72

  return global_tower_iPhi;
}

/*
 * GCT card convention -> Global GCT convention
 * Get tower's global iEta from the tower's iEta inside the GCT card (0-34).
 */
 int getTower_global_toweriEta_fromGCTcard(unsigned int nGCTCard, unsigned int iEtaInGCTCard) {
  (void) nGCTCard;
  int global_iEta = iEtaInGCTCard;
  return global_iEta;
 }

 /*
  * GCT card convention -> Global GCT convention
  * Get tower's global iPhi from the GCT card number (0, 1, or 2) and the tower's iPhi inside the GCT card (0-32).
  */
int getTower_global_toweriPhi_fromGCTcard(unsigned int nGCTCard, unsigned int iPhiInGCTCard) {
  assert(nGCTCard <= 2);  // Make sure the card number is valid
  int toweriPhi_card_offset = 0;   
  if      (nGCTCard == 0) toweriPhi_card_offset = GCTCARD_0_TOWER_IPHI_OFFSET;
  else if (nGCTCard == 1) toweriPhi_card_offset = GCTCARD_1_TOWER_IPHI_OFFSET;
  else if (nGCTCard == 2) toweriPhi_card_offset = GCTCARD_2_TOWER_IPHI_OFFSET;

  int global_iPhi = (toweriPhi_card_offset + iPhiInGCTCard) % (n_towers_Phi);  //   n_towers_Phi = 72

  return global_iPhi;
}


GCTcard_t getClustersCombined(const GCTcard_t& GCTcard, unsigned int nGCTCard){

  GCTcard_t GCTcombinedClusters ;

  // Initialize the output
  for(int i = 0; i < N_RCTCARDS_PHI; i++){
    for(int j = 0; j < N_RCTGCT_FIBERS; j++){
      for(int k = 0; k < N_RCTCLUSTERS_FIBER; k++){
        GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k] = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k];
	      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k] = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k];
      }
    }
  }
  // we will store new et in the GCTcombinedClusters, 0'ing lower clusters after stitching, dont need to care about other variables they stay the 
  // same as input for now at least
  // we combine even phi boundaries positive eta, when combined the lowest et is set to 0
  
  for(int i=0; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
        ap_uint<15> eta1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
        ap_uint<15> phi1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;

        for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
          for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
            ap_uint<15> eta2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
            ap_uint<15> phi2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;

            ap_uint<15> dPhi ; dPhi=((5 - phi1) + phi2) ;    
            ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;

            // Proximity stitching: even phi boundaries, positive eta
            if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3 && phi1 == 4){
              if( GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 && phi2 == 0) {
                if(dEta < 2) {
                  ap_uint<12> one = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two){ 
                    // std::cout<< "merging.." << std::endl;
                    GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;		    
                    GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                  }
                  else {
                    //  std::cout<< "merging.." << std::endl;
                    GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                    GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                  }
                }
              }
            }
            // Brems correction: even phi boundaries, positive eta. Note no requirement on crystal phi in below line
            if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3) {
              if( GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 ) {
                if( (dPhi <= 5) && (dEta < 2) ) {  
                  ap_uint<12> one = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two){ 
                    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (two > (0.10 * one)) {
                      std::cout<< "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;		    
                      GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                    }
                    else {
                      std::cout << "energy of 'two' is insufficient for merge" << std::endl;
                    }
                  }
                  else {
                    // Test: only stitch if energy is >10% of the other
                    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (one > (0.10 * two)) {
                      std::cout<< "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                      GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                    }
                    else {
                      std::cout << "energy of 'one' is insufficiennt for merge" << std::endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  // now we combine odd phi boundaries positive eta
    
  for(int i=1; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
        ap_uint<15> eta1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
        ap_uint<15> phi1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;

        for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){ 
          for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
            ap_uint<15> eta2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
            ap_uint<15> phi2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;

            ap_uint<15> dEta = (eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
            ap_uint<15> dPhi = ((5 - phi1) + phi2) ; 

            // Proximity stitching: odd phi boundaries, positive eta
            if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3 && phi1 == 4) {
              if(GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 && phi2 == 0) {
                if(dEta <2) {
                  ap_uint<12> one = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two){
                    // std::cout<< "merging.." << std::endl;
                    GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
                    GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                  }
                  else {
                    // std::cout<< "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                      GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                  }
                }
              }
            }

            // Brems correction: odd phi boundaries, positive eta. Clusters only need to be opposite in towers
            if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3) {
              if( GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 ) {
                if( (dPhi <= 5) && (dEta < 2) ) {         
                  ap_uint<12> one = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two){
                    // Test: only stitch if energy of one is >10% of the other                                                           
                    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (two > (0.10 * one)) {
                      std::cout<< "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
                      GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                    }
                    else {
                      std::cout << "energy of 'two' is insufficient for merge" << std::endl;
                    }
                  }
                  else {
                    // Test: only stitch if energy of one is >10% of the other          
                    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (one > (0.10 * two)) { 
                      std::cout<< "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                      GCTcombinedClusters.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                    }
                    else {
                      std::cout << "energy of 'one' is insufficient for merge" << std::endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  // Debugging printouts
  std::cout << ">>> GCT card " << nGCTCard << ": print all negative eta clusters:" << std::endl;
  for (int i=0; i < N_RCTCARDS_PHI; i++) {
    std::cout << "  ... RCT Card " << i << std::endl;
    for (int j=0; j < N_RCTGCT_FIBERS; j++) {
        for (int k=0; k < N_RCTCLUSTERS_FIBER; k++) {
          std::cout << "Fiber " << j << ", cluster num. " << k << ": energy " << GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et/8.0
                    << " with tower phi " << GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi 
                    << " and crystal phi" << GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi << std::endl;
      }
    }
  }
    
  // repeat above steps for negative eta
  // even phi boundaries
    
  for(int i=0; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
        ap_uint<15> eta1 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
        ap_uint<15> phi1 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;

        for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
          for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){

            ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
            ap_uint<15> phi2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;

            ap_uint<15> dPhi = ((5 - phi2) + phi1) ;  // reversed for negative eta
            ap_uint<15> dEta = (eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;

            // Proximity stitching: clusters need to be exactly across each other in crystals
            if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0 && phi1 == 0 ) {
              if(GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 && phi2 == 4 ) {
                if(dEta < 2) {
                  ap_uint<12> one = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two){ 
                      //std::cout<< "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
                      GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                  }                
                  else {
                      //std::cout << "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                      GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                  }
                }
              }
            }

            // Brems correction: even phi boundaries, negative eta: clusters only need to be across each other in towers
            if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0) {
              if(GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3) {
                if( (dPhi <= 5) && (dEta < 2) ) {   
                  ap_uint<12> one = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two){ 
                    // Test: only stitch if energy of one is >10% of the other                                                         
                    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (two > (0.10 * one)) {
                      std::cout<< "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
                      GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                    }
                    else { 
                      std::cout << "Energy of 'two' was insufficient for merge" << std::endl; 
                      }
                  }
                  else {
                    // Test: 
                    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (one > (0.10 * two)) {
                      std::cout << "merging.." << std::endl;
                      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                      GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                    }
                    else { 
                      std::cout <<"Energy of 'one' was insufficient for merge" <<std::endl; 
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
    
  // now we combine odd boundaries 
    
  for(int i=1; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
        ap_uint<15> eta1 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
        ap_uint<15> phi1 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;

        for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
          for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){

            ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
            ap_uint<15> phi2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;

            ap_uint<15> dPhi = ((5 - phi2) + phi1) ; // reversed compared to positive eta
            ap_uint<15> dEta = (eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;

            // Proximity stitching: the clusters need to be directly across the RCT card boundary in crystals
            if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0 && phi1 == 0) { 
              if(GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 && phi2 == 4 ) {
                if(dEta < 2) {        
                  ap_uint<12> one = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two) { 
                    // std::cout << "merging..." << std::endl;
                    GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
                    GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                  }
                  else {
                    // std::cout << "merging..." << std::endl;
                    GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                    GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                  }
                }
              }
            }

            // Brems stitching: odd boundaries, negative eta. The clusters should be in the towers saddling the boundary 
            if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0) {
              if( GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 ) {
                if( (dPhi <= 5) && (dEta < 2) ) {        
                  ap_uint<12> one = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two){ 
                    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (two > (0.10 * one)) {
                      std::cout << "merging..." << std::endl;
                      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
                      GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                    }
                    else {
                      std::cout << "energy insufficient to merge" << std::endl;
                    }
                  }
                  else {
                    // Test: only stitch if energy of one is >10% of the other                                                           
                    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (one > (0.10 * two)) {
                      std::cout << "merging..." << std::endl;
                      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                      GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                    }
                    else {
                      std::cout <<"energy insufficient to merge" << std::endl;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  
   // we need to store what we did before we start phi stitching

  GCTcard_t GCTout ;
  for(int i = 0; i < N_RCTCARDS_PHI; i++){
    for(int j = 0; j < N_RCTGCT_FIBERS; j++){
      for(int k = 0; k < N_RCTCLUSTERS_FIBER; k++){
        GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k] = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k];
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k] = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k];
      }
    }
  }
          
  // now we combine eta boundaries, just positive and negative eta 
     
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
        ap_uint<15> phi1 = (i*4+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi)*5+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;
        ap_uint<15> eta1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
        if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta == 0 && eta1 == 0 ) {
          for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
            for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
              ap_uint<15> phi2 = (i*4+(3-GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi))*5+(4-GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi) ;
              ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
              if( GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].towEta == 0 && eta2 == 0 ) {
                ap_uint<15> dPhi ; dPhi=(phi1 > phi2)?(phi1-phi2):(phi2-phi1) ;
                if( dPhi < 2 ) {
                  ap_uint<12> one = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
                  ap_uint<12> two = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
                  if (one > two){ 
                    GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
                    GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
                  }
                  else {
                    GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
                    GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return GCTout ;
}

/*
 * Populate a GCTinternal_t struct (consisting of 64 fibers, each fiber has clusters and towers) by converting RCT clusters and towers to GCT notation.
 */ 

GCTinternal_t getClustersTowers(const GCTcard_t& GCTcard, unsigned int nGCTCard){

  GCTcard_t GCTcombinedClusters ;
  GCTinternal_t GCTout ;

  // here we will stitch the clusters in phi and eta
  GCTcombinedClusters = getClustersCombined(GCTcard, nGCTCard) ;

  // create internal structure of GCT card
  // we start from RCT card 0 - it is overlap with other GCT card and fill structure that we will use to send data to Correlator
  // we only need to care about clusters et in combinrdClusters, since the rest remains unchanged wrt input, the cluster that we set to 0
  // remain in the data at the same place , it will just get 0 et now
  // we need to code Positive and Negative Eta differently !  For negative Eta link 0 for each RCT
  // region becomes 3 in GCT output, the RCT card is rotated around 0:0 point of the card 
  // First 16 fibers - positive Eta , second 16 - negative. Eta coded 0...16 and towEtaNeg = 0 or 1 for clusters ; 
  // Phi is coded 0...15 , in case if whole card 0...33 and subdevision 1/5 in crPhi and crEta 0...4 for 
  // position in tower
  //
  // towers are put in link starting from eta=0, the link number defines Eta negative or positive and Phi position of tower.
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
        // positive eta: cluster eta/phi info
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].et    = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].et2x5 = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 ;
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].et5x5 = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 ;
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].is_ss = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss ;
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].is_looseTkss = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss;

        // positive eta: cluster eta/phi info
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towEtaNeg  = 0 ;
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towEta  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta  ;
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].towPhi  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi+i*4  ;
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].crEta  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta  ;
        GCTout.GCTCorrfiber[i*4+j].GCTclusters[k].crPhi  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi  ;

        // negative eta: cluster energy info
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].et    = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].et2x5 = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].et5x5 = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].is_ss = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].is_looseTkss = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss;
        
        // negative eta: cluster eta/phi info
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].towEtaNeg  = 1 ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].towEta  = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towEta  ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].towPhi  = (3-GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi)+i*4  ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].crEta  = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crEta  ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTclusters[k].crPhi  = (4-GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi)  ;

      }
      for(int k=0; k<N_RCTTOWERS_FIBER; k++){
        GCTout.GCTCorrfiber[i*4+j].GCTtowers[k].et     = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTtowers[k].et  ;
        GCTout.GCTCorrfiber[i*4+j].GCTtowers[k].hoe    = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTtowers[k].hoe ;
        GCTout.GCTCorrfiber[i*4+j].GCTtowers[k].ecalEt = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTtowers[k].ecalEt ;
        GCTout.GCTCorrfiber[i*4+j].GCTtowers[k].hcalEt = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTtowers[k].hcalEt ; 
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTtowers[k].et     = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTtowers[k].et  ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTtowers[k].hoe    = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTtowers[k].hoe ; 
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTtowers[k].ecalEt = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTtowers[k].ecalEt ;
        GCTout.GCTCorrfiber[i*4+(3-j)+N_GCTPOSITIVE_FIBERS].GCTtowers[k].hcalEt = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTtowers[k].hcalEt ; 
      }
    }
  }
  return GCTout ;
}

/*
 * Return full towers with the tower energy (i.e. unclustered energy) and cluster energy added together.
 */
GCTintTowers_t  getFullTowers(const GCTinternal_t& GCTinternal) {
  GCTintTowers_t GCTintTowers;
  
  // Positive eta
  for(int i=0; i<N_GCTPOSITIVE_FIBERS; i=i+4){
    for(int i1=0; i1<4; i1++){
      for(int k=0; k<N_GCTTOWERS_FIBER; k++){
	      ap_uint<15> phi = i+i1 ; 
	      ap_uint<15> eta = N_GCTETA/2 + k ; 
	      GCTintTowers.GCTtower[eta][phi].et  = GCTinternal.GCTCorrfiber[phi].GCTtowers[k].et ;
        GCTintTowers.GCTtower[eta][phi].hoe = GCTinternal.GCTCorrfiber[phi].GCTtowers[k].hoe ;
	      for(int ic1=0; ic1<4; ic1++){
	        for(int jc=0; jc<N_GCTCLUSTERS_FIBER; jc++){
	          ap_uint<15> eta1 = N_GCTETA/2 + GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].towEta ; 
	          ap_uint<15> phi1 = GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].towPhi ; 
	          if( eta == eta1 && phi == phi1) {
              GCTintTowers.GCTtower[eta][phi].et = (GCTintTowers.GCTtower[eta][phi].et + GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].et) ;
            }
          }
        }
      }
    }
  }

  // Negative eta
  for(int i=N_GCTPOSITIVE_FIBERS; i<N_GCTINTERNAL_FIBERS; i=i+4){
    for(int i1=0; i1<4; i1++){
      for(int k=0; k<N_GCTTOWERS_FIBER; k++){
	      ap_uint<15> eta = N_GCTETA/2 - k - 1 ; 
	      ap_uint<15> phi = i+i1-N_GCTPOSITIVE_FIBERS ; 
	      GCTintTowers.GCTtower[eta][phi].et  = GCTinternal.GCTCorrfiber[i+i1].GCTtowers[k].et ;
        GCTintTowers.GCTtower[eta][phi].hoe = GCTinternal.GCTCorrfiber[i+i1].GCTtowers[k].hoe ;
	      for(int ic1=0; ic1<4; ic1++){
	        for(int jc=0; jc<N_GCTCLUSTERS_FIBER; jc++){
	          ap_uint<15> eta1 = N_GCTETA/2 - 1 - GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].towEta ;
	          ap_uint<15> phi1 = GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].towPhi ;
	          if( eta == eta1 && phi == phi1) {
              GCTintTowers.GCTtower[eta][phi].et = (GCTintTowers.GCTtower[eta][phi].et + GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].et);
            }
          }
        }
      }
    }
  }
  
  return GCTintTowers ;
}

/*
 * Compute isolation (sum of unclustered energy in 7x7 window IN TOWERS) for a single cluster in a GCTinternal_t struct,
 * at iFiber and iCluster. Needs the full GCTinternal_t to access tower energies for isolation sum.
 */ 
void computeIso(GCTinternal_t& GCTinternal, unsigned int iFiber, unsigned int iCluster, int nGCTCard, l1tp2::ParametricCalibration calib_) {


  // We will only save clusters with > 0 GeV, so only need to do this for clusters with >0 energy 
  if (GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster].et == 0) {
    GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster].iso = 0;
    return;
  }

  ap_uint<12> uint_isolation = 0;

  bool getGlobal_iPhi = false;   // for the phi function: do not add the GCT card off-set, so we remain in the
  // gct local card iEta/iPhi
  int toweriEta_in_GCT_card = getCluster_global_tower_iEta(nGCTCard, GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster]);
  int toweriPhi_in_GCT_card = getCluster_global_tower_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster], getGlobal_iPhi);

  // If cluster is in the overlap region, do not compute isolation 
  bool inOverlapWithAnotherGCTCard = ( ((toweriPhi_in_GCT_card >= 0) && (toweriPhi_in_GCT_card < 4)) || ((toweriPhi_in_GCT_card >= 28) && (toweriPhi_in_GCT_card < 32)) );
  if (inOverlapWithAnotherGCTCard) {
    GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster].iso = 0;
    return;
  }

  // Size 5x5 in towers: include the overlap-region-between-GCT-cards-if-applicable. In eta direction, the min and max towers (inclusive!) are:
  int isoWindow_toweriEta_in_GCT_card_min = std::max(0, toweriEta_in_GCT_card - 2);
  int isoWindow_toweriEta_in_GCT_card_max = std::min(toweriEta_in_GCT_card + 2, N_GCTETA - 1);  // N_GCTETA = 34
  // e.g. if our window is centered at tower_iEta = 5, we want to sum towers_iEta 3, 4, (5), 6, 7, inclusive 
  // e.g. if our window is near the boundary, tower_iEta = 32, we want to sum towers_iEta 30, 31, (32), 33
  // inclusive (but there are only N_GCTETA = 34 towers, so we stop at tower_iEta = 33)
  
  // in phi direction, the min and max towers (inclusive!) are:
  int isoWindow_toweriPhi_in_GCT_card_min = std::max(0, toweriPhi_in_GCT_card - 2);
  int isoWindow_toweriPhi_in_GCT_card_max = std::min(toweriPhi_in_GCT_card + 2, N_GCTPHI - 1);  
  
  // Keep track of the number of towers we summed over
  int nTowersSummed = 0;

  std::cout << "[Computing isolation...]" << std::endl;
  std::cout << " - Current cluster: " << std::endl;
  std::cout << "    pT " <<  GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster].et/8.0 << " GeV" << " at ("
            << getEta_fromCrystaliEta(getCluster_global_cluster_iEta(nGCTCard, GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster]))
            << ", "
            << getPhi_fromCrystaliPhi(getCluster_global_cluster_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster]))
            << ") " 
            << " (cluster global iEta/iPhi: " << getCluster_global_cluster_iEta(nGCTCard, GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster]) << ", " 
            << getCluster_global_cluster_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster])  << ")"
            << std::endl;

  std::cout << " - Isolation components: " << std::endl;

  // First add any nearby clusters to the isolation
  for (unsigned int candFiber = 0; candFiber < N_GCTINTERNAL_FIBERS; candFiber++) {
    for (unsigned int candCluster = 0; candCluster < N_GCTCLUSTERS_FIBER; candCluster++ ) {

      // Do not double-count the cluster we are calculating the isolation for
      if ( !((candFiber == iFiber) && (candCluster == iCluster))) {

        // Only consider clusters with et > 0 for isolation sum
        if (GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster].et > 0) {

          // Get the candidate cluster's tower iEta and iPhi in GCT card
          bool getGlobal_iPhi = false;   // for the phi function: do not add the GCT card off-set (so we remain in the gct local card iEta/iPhi)

          int candidate_toweriEta = getCluster_global_tower_iEta(nGCTCard, GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster]);
          int candidate_toweriPhi = getCluster_global_tower_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster], getGlobal_iPhi);

          // If the tower that the candidate cluster is in, is within a 5x5 window, add the candidate cluster energy's to the isolation as a proxy for the ECAL energy
          if ( ( (candidate_toweriEta >= isoWindow_toweriEta_in_GCT_card_min) && (candidate_toweriEta <= isoWindow_toweriEta_in_GCT_card_max) ) &&
                ( (candidate_toweriPhi >= isoWindow_toweriPhi_in_GCT_card_min) && (candidate_toweriPhi <= isoWindow_toweriPhi_in_GCT_card_max) ) ) {
            
            uint_isolation += GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster].et;

            bool getGlobalPhi = true; // just for print-out statements
            std::cout << "    Adding cluster pT " << GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster].et/8.0 << " GeV " << " at (" 
                      << getEta_fromCrystaliEta(getCluster_global_cluster_iEta(nGCTCard, GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster]))
                      << ", "
                      << getPhi_fromCrystaliPhi(getCluster_global_cluster_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster], getGlobalPhi))
                      << ") "
                      << " (cluster global iEta/iPhi: " << getCluster_global_cluster_iEta(nGCTCard, GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster]) << ", " 
                      << getCluster_global_cluster_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[candFiber].GCTclusters[candCluster], getGlobalPhi)  << ")"
                      << std::endl;
          }
        }
      }
    }
  }



  //  From "tower index in GCT card", get which fiber it is in (out of 64 fibers), and which tower it is inside the fiber (out of 17 towers)
  for (int iEta = isoWindow_toweriEta_in_GCT_card_min; iEta <= isoWindow_toweriEta_in_GCT_card_max; iEta++) {
    for (int iPhi = isoWindow_toweriPhi_in_GCT_card_min; iPhi <= isoWindow_toweriPhi_in_GCT_card_max; iPhi++) {
      
      nTowersSummed += 1;
      
      int indexInto64Fibers;
      int indexInto17TowersInFiber; 
      
      bool isTowerInPositiveEta = (iEta >= N_GCTTOWERS_FIBER); 
      if (isTowerInPositiveEta) { 
        // phi index is simple (e.g. if real phi = +80 degrees, iPhi in GCT = 31)
        indexInto64Fibers = iPhi; 
        // if real eta = 1.47, iEta in GCT card = 33. If real eta = 0.0, iEta in GCT = 17, so iEta in fiber = 17%17 = 0.
        indexInto17TowersInFiber = (iEta % 17); 
      }
      else { 
        // add offset (e.g. if real phi = +80 degrees, iPhi in GCT = 31, and my index into GCT fibers 31 + 32 = 63)
        indexInto64Fibers = (iPhi + N_GCTPOSITIVE_FIBERS); 
        // e.g.  if real eta = 0, iEta in GCT card = 16, i.e. our index into the GCT fiber is 16-16 = 0
        indexInto17TowersInFiber = (16 - iEta); 
      }
      
      // std::cout << "... indexInto64Fibers: " << indexInto64Fibers << std::endl;
      // std::cout << "... indexInto17TowersInFiber: " << indexInto17TowersInFiber << std::endl;
      
      ap_uint<12> towerEt = GCTinternal.GCTCorrfiber[indexInto64Fibers].GCTtowers[indexInto17TowersInFiber].et;
   //   ap_uint<12> hcalEtInEcalConvention =  convertHcalETtoEcalET(GCTinternal.GCTCorrfiber[indexInto64Fibers].GCTtowers[indexInto17TowersInFiber].hcalEt);
      ap_uint<12> ecalEt = GCTinternal.GCTCorrfiber[indexInto64Fibers].GCTtowers[indexInto17TowersInFiber].ecalEt;

      std::cout << "    Adding tower pT " << ecalEt/8.0 << " GeV" << " at (" 
                << getTowerEta_fromAbsID(getTower_global_toweriEta(nGCTCard, indexInto64Fibers, indexInto17TowersInFiber)) << ", "
                << getTowerPhi_fromAbsID(getTower_global_toweriPhi(nGCTCard, indexInto64Fibers, indexInto17TowersInFiber)) << ") "
                << "(corresponds to calib_(0, " << std::abs(getTowerEta_fromAbsID(getTower_global_toweriEta(nGCTCard, indexInto64Fibers, indexInto17TowersInFiber))) << "): "
                << calib_(0, std::abs(getTowerEta_fromAbsID(getTower_global_toweriEta(nGCTCard, indexInto64Fibers, indexInto17TowersInFiber))))
                << ")"
                << std::endl;
	    uint_isolation += ecalEt;
    }
  }
  
  // Scale the isolation sum up if we summed over fewer than (5x5) = 25 towers
  float scaleFactor = ((float) (N_GCTTOWERS_CLUSTER_ISO_ONESIDE * N_GCTTOWERS_CLUSTER_ISO_ONESIDE) / (float) nTowersSummed);

   std::cout << std::endl;
  std::cout << "    Total: " << uint_isolation/8.0 << " GeV" << std::endl;

  std::cout << " - Isolation: " << std::endl;
  std::cout << "    Summed over " << nTowersSummed << " towers" << std::endl;
  std::cout << "    Scaling by (25.0/" << nTowersSummed << ")" << std::endl;
  std::cout << "    Total, scaled: " << (uint_isolation/8.0) * scaleFactor << std::endl;
  std::cout << "    Relative isolation: " << (uint_isolation/8.0) * scaleFactor / (GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster].et/8.0)  << std::endl;
  std::cout << std::endl;


  uint_isolation = (ap_uint<12>) (((float) uint_isolation) * scaleFactor);
  
  // Set the iso in the cluster
  GCTinternal.GCTCorrfiber[iFiber].GCTclusters[iCluster].iso = uint_isolation;

}

/*
 * Compute relative isolation for a GCT cluster and set its flags in-place, and return 1.
 */ 
int computeRelIsoAndFlags(GCTcluster_t &cluster) {

  float relative_iso = 0;
  if (cluster.et > 0) {
	  relative_iso = ( ((float) cluster.iso / 8 ) / ( (float) cluster.et / 8 ) );
  }
  cluster.relIso  = relative_iso;
  cluster.is_iso        = passes_iso(cluster.et/8.0, cluster.relIso);
  cluster.is_looseTkiso = passes_looseTkiso(cluster.et/8.0, cluster.relIso);
  
  return 1;
}

/*
 * Compute cluster isolation for an entire GCTinternal_t struct given the GCT card index as well. Returns 1.
 */
int computeClusterIsolationsForGCTCard(GCTinternal_t &gctInternal, int nGCTCard, l1tp2::ParametricCalibration calib_) {
  for (unsigned int iFiber = 0; iFiber < N_GCTINTERNAL_FIBERS; iFiber++) {
    for (unsigned int iCluster = 0; iCluster < N_GCTCLUSTERS_FIBER; iCluster++ ) {

      computeIso(gctInternal, iFiber, iCluster, nGCTCard, calib_);
      computeRelIsoAndFlags(gctInternal.GCTCorrfiber[iFiber].GCTclusters[iCluster]);
     
    }
  }
  return 1;
}

/* 
 * algo_top: First two arguments are the same as in the original firmware.
 * nGCTCard is 0, 1, or 2 (needed for getting the cluster real eta/phis for CMSSW collections).
 * gctClusters is the CMSSW-style output collection of clusters.
 * gctTowers is the CMSSW-style output collection of towers.
 */

void algo_top(const GCTcard_t& GCTcard, GCTtoCorr_t& GCTtoCorr,
              unsigned int nGCTCard,
              std::unique_ptr<l1tp2::CaloCrystalClusterCollection> const& gctClusters,
              std::unique_ptr<l1tp2::CaloTowerCollection> const& gctTowers,
              std::unique_ptr<l1tp2::CaloTowerCollection> const& gctFullTowers,
              l1tp2::ParametricCalibration calib_) {
  
  //-------------------------//
  // Initialize the GCT area 
  //-------------------------//
  GCTinternal_t GCTinternal = getClustersTowers(GCTcard, nGCTCard);

  //------------------------------------------------//
  // Combine towers and clusters to get full towers
  //------------------------------------------------//
  GCTintTowers_t GCTintTowers = getFullTowers(GCTinternal);

  //---------------------------//
  // Compute cluster isolation
  //--------------------------//
  computeClusterIsolationsForGCTCard(GCTinternal, nGCTCard, calib_); 

  //-----------------------------------------------------------------------------------------------------------------------//
  // Output to correlator, positive eta. Skip overlap region, i.e. fibers i = 0, 1, 2, 3, and i = 28, 29, 30, 31.
  //-----------------------------------------------------------------------------------------------------------------------//
  for(int i = 4; i < (N_GCTPOSITIVE_FIBERS-N_RCTGCT_FIBERS); i++){
    for(int k = 0; k < N_GCTCLUSTERS_FIBER; k++){

      // Tower iEta is from 0-16 where iEta = 0 is real eta = 0 (indexes increase towards larger abs(eta).
      // Tower iPhi is from 0-3 where iPhi = 0 is the leftmost (in a 'sideways' diagram of the GCT card like on the TWiki). 

      GCTcluster_t posCluster = GCTinternal.GCTCorrfiber[i].GCTclusters[k];

      // Write the GCT internal clusters to the output to correlator: all fields are the same with the exception of towPhi, which
      // needs to be subtracted by 4 becauuse the output to correlator does NOT include the overlap region.
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k] = posCluster;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towPhi  =  posCluster.towPhi-4 ;

      // Get the real eta, phi using two helper functions
      int crystaliEta_in_barrel = getCluster_global_cluster_iEta(nGCTCard, posCluster);
      int crystaliPhi_in_barrel = getCluster_global_cluster_iPhi(nGCTCard, posCluster);
      float realEta = getEta_fromCrystaliEta(crystaliEta_in_barrel);
      float realPhi = getPhi_fromCrystaliPhi(crystaliPhi_in_barrel);

      reco::Candidate::PolarLorentzVector p4cluster(posCluster.et/8.0,
                                                    realEta,
                                                    realPhi,
                                                    0.);
      l1tp2::CaloCrystalCluster cluster(p4cluster, 
                                        posCluster.et/8.0,   // convert to float
                                        posCluster.iso/8.0,  // TECHNICALLY INCORRECT: for testing only, supposed to be float h over e                              
                                        posCluster.relIso,		   // for consistency with the old emulator, in this field save (iso energy sum)/(cluster energy)
                                        0,  // DetId seedCrystal                              
                                        0,  // puCorrPt                                           
                                        0,  // 0, 1, or 2 (as computed in firmware)                
                                        0,  // et2x2 (not calculated)                             
                                        posCluster.et2x5/8.0,  // et2x5 (as computed in firmware, save float)           
                                        0,  // et3x5 (not calculated)                             
                                        posCluster.et5x5/8.0,   // et5x5 (as computed in firmware, save float)  
                                        posCluster.is_ss,  // standalone WP: not computed
                                        posCluster.is_ss, // electronWP98: not computed 
                                        false, // is_photon in Cecile's emulator, photonWP80: not computed
                                        posCluster.is_ss, // electronWP90: not computed
                                        posCluster.is_looseTkss, // looseL1TkMatchWP
                                        posCluster.is_ss  // stage2effMatch: not computed
                                        );

        // Flags
        std::map<std::string, float> params;
        params["standaloneWP_showerShape"] = posCluster.is_ss;
        params["standaloneWP_isolation"]   = posCluster.is_iso;
        params["trkMatchWP_showerShape"]   = posCluster.is_looseTkss;
        params["trkMatchWP_isolation"]     = posCluster.is_looseTkiso;
        cluster.setExperimentalParams(params);

        if (cluster.pt() > 0.0) {
          gctClusters->push_back(cluster);
          // std::cout << "--- cluster pT, global iEta, iPhi and real eta, phi: "
          //           << posCluster.et/8.0  << ", "
          //           << "(" << crystaliEta_in_barrel << ", " << crystaliPhi_in_barrel << "), "
          //           << "(" << realEta << "," << realPhi << ")"
          //           << " with relative isolation " << posCluster.relIso
          //           << std::endl;
          // std::cout << "    with the GCTinternal values: " << std::endl;
          // printGCTClusterInfo(posCluster, "positive cluster writeout");
	      }
    }
    // Positive eta towers : push back to CMSSW collection
    for(int k=0; k<N_GCTTOWERS_FIBER; k++){
      // std::cout<< "Accessing positive eta: GCTCorrfiber " << i-4
      // 	       << " , GCTtowers " << k
      // 	       << " , energy " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].et
	    //    << " , hoe " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].hoe 
	    //    << " , ecalEt " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].ecalEt
      //          << " , hcalEt " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].hcalEt
	    //    << std::endl;
      GCTtoCorr.GCTCorrfiber[i-4].GCTtowers[k] = GCTinternal.GCTCorrfiber[i].GCTtowers[k];

      l1tp2::CaloTower l1CaloTower;
      l1CaloTower.setEcalTowerEt(GCTinternal.GCTCorrfiber[i].GCTtowers[k].ecalEt/8.0); // float: ECAL divide by 8.0
      float hcalLSB = 0.5;
      l1CaloTower.setHcalTowerEt(GCTinternal.GCTCorrfiber[i].GCTtowers[k].hcalEt * hcalLSB);  // float: HCAL multiply by LSB 
      int global_toweriEta = getTower_global_toweriEta( nGCTCard, i, k );
      int global_toweriPhi = getTower_global_toweriPhi( nGCTCard, i, k );
      l1CaloTower.setTowerIEta( global_toweriEta );
      l1CaloTower.setTowerIPhi( global_toweriPhi );
      l1CaloTower.setTowerEta( getTowerEta_fromAbsID( global_toweriEta ) );
      l1CaloTower.setTowerPhi( getTowerPhi_fromAbsID( global_toweriPhi ) );
      
      gctTowers->push_back(l1CaloTower);
    }
  }
  //-----------------------------------------------------------------------------------------------------------------------//
  // Output to correlator: In negative eta, the overlap region to skip is fibers 32, 33, 34, 35, and 61, 62, 63, 64.
  //-----------------------------------------------------------------------------------------------------------------------//
  for (int i = (N_GCTPOSITIVE_FIBERS+N_RCTGCT_FIBERS); i < (N_GCTINTERNAL_FIBERS-N_RCTGCT_FIBERS); i++) {
    for (int k = 0; k < N_GCTCLUSTERS_FIBER; k++) {

      GCTcluster_t negCluster = GCTinternal.GCTCorrfiber[i].GCTclusters[k];

      // Write the GCT internal clusters to the output to correlator: all fields are the same with the exception of towPhi, which
      // needs to be subtracted by 4 becauuse the output to correlator does NOT include the overlap region.
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k] = negCluster;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].towPhi  =  negCluster.towPhi-4 ;

      // Get the real eta, phi using two helper functions
      int globaliEta = getCluster_global_cluster_iEta(nGCTCard, negCluster);
      int globaliPhi = getCluster_global_cluster_iPhi(nGCTCard, negCluster);
      float realEta = getEta_fromCrystaliEta(globaliEta);
      float realPhi = getPhi_fromCrystaliPhi(globaliPhi);
      
      reco::Candidate::PolarLorentzVector p4cluster(negCluster.et/8.0,
                                                    realEta,
                                                    realPhi,
                                                    0.);
      l1tp2::CaloCrystalCluster cluster(p4cluster, 
                                        negCluster.et/8.0, // conver to float
                                        negCluster.iso/8.0,  // TECHNICALLY INCORRECT, for studies only, supposed to be .hoe (float h over e)                              
                                        negCluster.relIso, // follow old emulator's convention: save relative iso in this field
                                        0,  // DetId seedCrystal                              
                                        0,  // puCorrPt                                           
                                        0,  // 0, 1, or 2 (as computed in firmware)                
                                        0,  // et2x2 (not calculated)                             
                                        negCluster.et2x5/8.0,  // et2x5 (as computed in firmware, save float)           
                                        0,  // et3x5 (not calculated)                             
                                        negCluster.et5x5/8.0,   // et5x5 (as computed in firmware, save float)  
                                        negCluster.is_ss,  // standalone WP
                                        negCluster.is_ss, // electronWP98: not computed 
                                        false, // photonWP80: not computed
                                        negCluster.is_ss, // electronWP90: not computed
                                        negCluster.is_looseTkss, // looseL1TkMatchWP
                                        negCluster.is_ss  // stage2effMatch: not computed
                                        );

      // Experimental parameters
      std::map<std::string, float> params;
      params["standaloneWP_showerShape"] = negCluster.is_ss;
      params["standaloneWP_isolation"]   = negCluster.is_iso;
      params["trkMatchWP_showerShape"]   = negCluster.is_looseTkss;
      params["trkMatchWP_isolation"]     = negCluster.is_looseTkiso;
      cluster.setExperimentalParams(params);

      // Push back negative eta clusters to the output CMSSW collection
      if (cluster.pt() > 0.0) {
	        gctClusters->push_back(cluster);
	        // std::cout << "--- cluster pT, global iEta, iPhi and real eta, phi: "
          //           << negCluster.et/8.0  << ", "
          //           << "(" << globaliEta << ", "  << globaliPhi << "), "
          //           << "(" << realEta << "," << realPhi << ")"
          //           << " with relative isolation " << negCluster.relIso
          //           << std::endl;
          // std::cout << "... with the GCTinternal values: " << std::endl;
          // printGCTClusterInfo(negCluster, "negative cluster writeout");
	      }
      
    }
    // Negative eta towers : push back to CMSSW as well 
    for(int k=0; k<N_GCTTOWERS_FIBER; k++){
      // std::cout<< "Accessing negative eta: GCTCorrfiber " << i
      // 	       << " , GCTtowers " << k
      // 	       << " , energy " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].et 
	    //    << " , hoe " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].hoe 
	    //    << " , ecalEt " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].ecalEt 
	    //    << " , hcalEt " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].hcalEt 
	    //    << std::endl;
      GCTtoCorr.GCTCorrfiber[i-12].GCTtowers[k] = GCTinternal.GCTCorrfiber[i].GCTtowers[k];

      l1tp2::CaloTower l1CaloTower;
      l1CaloTower.setEcalTowerEt(GCTinternal.GCTCorrfiber[i].GCTtowers[k].ecalEt/8.0); // float: ECAL divide by 8.0
      float hcalLSB = 0.5;
      l1CaloTower.setHcalTowerEt(GCTinternal.GCTCorrfiber[i].GCTtowers[k].hcalEt * hcalLSB);  // float: HCAL multiply by LSB 
      int global_toweriEta = getTower_global_toweriEta( nGCTCard, i, k );
      int global_toweriPhi = getTower_global_toweriPhi( nGCTCard, i, k );
      l1CaloTower.setTowerIEta( global_toweriEta );
      l1CaloTower.setTowerIPhi( global_toweriPhi );
      l1CaloTower.setTowerEta( getTowerEta_fromAbsID( global_toweriEta ) );
      l1CaloTower.setTowerPhi( getTowerPhi_fromAbsID( global_toweriPhi ) );
      
      gctTowers->push_back(l1CaloTower);
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------//
  // CMSSW outputs for GCT Full Towers (clusters + towers) output for PFClusters.
  //-----------------------------------------------------------------------------------------------------------------------//
  for (unsigned int iEta = 0; iEta < N_GCTETA; iEta++) {
    for (unsigned int iPhi = 0; iPhi < N_GCTPHI; iPhi++) {

      l1tp2::CaloTower fullTower;
      // Store total Et (HCAL+ECAL) in the ECAL ET member.
      fullTower.setEcalTowerEt(GCTintTowers.GCTtower[iEta][iPhi].et/8.0); // convert to float by dividing by 8.0
      // Leave HCAL Tower Et = 0.
      int global_toweriEta = getTower_global_toweriEta_fromGCTcard(nGCTCard, iEta);
      int global_toweriPhi = getTower_global_toweriPhi_fromGCTcard(nGCTCard, iPhi);
      fullTower.setTowerIEta( global_toweriEta );
      fullTower.setTowerIPhi( global_toweriPhi );
      fullTower.setTowerEta( getTowerEta_fromAbsID( global_toweriEta ) );
      fullTower.setTowerPhi( getTowerPhi_fromAbsID( global_toweriPhi ) );
      
      gctFullTowers->push_back(fullTower);
    }
  
  }

}

#endif
