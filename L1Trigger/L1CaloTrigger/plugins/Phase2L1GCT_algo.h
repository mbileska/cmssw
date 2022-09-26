#include "Phase2L1RCT.h"
#include "Phase2L1GCT.h"

#ifndef _PHASE_2_L1_GCT_ALGO_H_
#define _PHASE_2_L1_GCT_ALGO_H_

// Get GCT cluster c's iEta (global iEta convention).
// Use with getEta_fromCrystaliEta from Phase2L1RCT.h to convert from GCT cluster
// to real eta.
int getCluster_global_iEta(unsigned int nGCTCard, GCTcluster_t c) {

  // First get the "iEta/iPhi" in the GCT card. i.e. in the diagram where the barrel
  // is split up into three GCT cards, (iEta, iPhi) = (0, 0) is the top left corner
  // of the GCT card. 
  int iEta_in_gctCard;

  // towEtaNeg = 1 (true) for negative eta
  if (c.towEtaNeg) {
    // c.towEta and c.crEta count outwards from the real eta = 0 center line, so for negative
    // eta, to convert to the barrel diagram global iEta (global iEta = 0 from LHS of page),
    // do (17*5 - 1) minus the GCT value.
    // e.g. If in GCT, a negative card's cluster had iEta = 84, this would be global iEta = 0.
    iEta_in_gctCard = ((17 * 5 - 1) - ((c.towEta * 5) + c.crEta));
  }
  else {
    // c.towEta and c.crEta count outwards from the real eta = 0 center line, so for positive
    // eta we need to add the 17*5 offset so that positive eta 0+epsilon starts at 17*5.
    // e.g. If in GCT, a positive card's cluster had iEta = 0, this would be global iEta = 85.
    // e.g. If in GCT, a positive card's cluster had iEta = 84, this would be global iEta = 169.
    iEta_in_gctCard = ((17 * 5) + ((c.towEta * 5) + c.crEta));
  }

  // Last, convert to the global iEta/iPhi in the barrel region. For eta there is nothing to
  // change.
  int iEta_in_barrel = iEta_in_gctCard;
  
  return iEta_in_barrel;
}

// Get GCT cluster c's iPhi (global convention).
// Use with getPhi_fromCrystaliPhi from Phase2L1RCT.h to convert from GCT cluster to real phi.
// If returnGlobalGCTiPhi is true (Default value) then return the iPhi in the entire GCT barrel. Otherwise
// just return the iPhi in the current GCT card.
int getCluster_global_iPhi(unsigned int nGCTCard, GCTcluster_t c, bool returnGlobalGCTiPhi = true) {
  
  // First get the "iEta/iPhi" in the GCT card. i.e. in the diagram where the barrel
  // is split up into three GCT cards, (iEta, iPhi) = (0, 0) is the top left corner
  // of the GCT card.  

  // Luckily, in the GCT algo convention and the global convention, iPhi always increases from the
  // top of the page to the bottom of page in the barrel diagram.
  int iPhi_in_gctCard = ((c.towPhi * 5) + c.crPhi);
  /* std::cout << ">>>>>>> getCluster_global_iPhi: invoked with inputs nGCTCard " << nGCTCard << ", cluster c.towPhi: " << c.towPhi << ", c.crPhi: "  */
  /* 	    << c.crPhi << std::endl; */
  /* std::cout << ">>>>>>> (contd.) iPhi_in_gctCard: " << iPhi_in_gctCard << std::endl; */

  // Last, convert to the global iEta/iPhi in the barrel region. For phi, we need to add the offset
  // of the GCT card in the phi direction, and modulo with the total number of crystals in the barrel
  // in the phi direction, since it wraps around.
  assert(nGCTCard <= 2);  // Make sure the card number is valid
  int iPhi_card_offset = 0;

  // (default behavior) If we should return the global GCT iPhi, get the iPhi offset due to the number of the GCT card
  if (returnGlobalGCTiPhi) {
    if      (nGCTCard == 0) iPhi_card_offset = GCTCARD_0_TOWER_IPHI_OFFSET * 5;  // tower #20, and five crystals per tower
    else if (nGCTCard == 1) iPhi_card_offset = GCTCARD_1_TOWER_IPHI_OFFSET * 5;   
    else if (nGCTCard == 2) iPhi_card_offset = GCTCARD_2_TOWER_IPHI_OFFSET * 5;
  }
  // Else, treat it as no offset due to GCT card number
  // std::cout << ">>>>>>> (contd.) returnGlobalGCTiPhi? " << returnGlobalGCTiPhi << ", giving iPhi_card_offset " << iPhi_card_offset << std::endl;
    

  int iPhi_in_barrel = (iPhi_card_offset + iPhi_in_gctCard) % (n_towers_Phi * 5); // detector wraps around in phi: n_towers_Phi = 72
  // std::cout << ">>>>>>> (contd.) total with offset (even if offset is 0): " << iPhi_card_offset + iPhi_in_gctCard << ", returning " << iPhi_in_barrel  << std::endl;

  return iPhi_in_barrel;
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

  int global_toweriEta; 
  
  bool isTowerInPositiveEta = (gctCorrFiberIdx < N_GCTPOSITIVE_FIBERS); // N_GCTPOSITIVE_FIBERS = 32
  
  if (isTowerInPositiveEta) {  
    // e.g. For positive eta, posInFiber = 0 is at real eta = 0, so global tower iEta is 0 + 17 = 17
    global_toweriEta = (N_GCTTOWERS_FIBER + posInFiber); // N_GCTTOWERS_FIBER = 17
  } else {
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

  int global_toweriPhi;
  
  unsigned int effectiveFiberIdx = (gctCorrFiberIdx % N_GCTPOSITIVE_FIBERS);  // N_GCTPOSITIVE_FIBERS = 32

  assert(nGCTCard <= 2);  // Make sure the card number is valid
  int toweriPhi_card_offset = 0;   
  if      (nGCTCard == 0) toweriPhi_card_offset = GCTCARD_0_TOWER_IPHI_OFFSET;
  else if (nGCTCard == 1) toweriPhi_card_offset = GCTCARD_1_TOWER_IPHI_OFFSET;
  else if (nGCTCard == 2) toweriPhi_card_offset = GCTCARD_2_TOWER_IPHI_OFFSET;

  global_toweriPhi = (toweriPhi_card_offset + effectiveFiberIdx) % (n_towers_Phi);  //  as explained above, effectiveFiberIdx is [0, 32). n_towers_Phi = 72

  return global_toweriPhi;
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


GCTcard_t getClustersCombined(const GCTcard_t& GCTcard){

  GCTcard_t GCTcombinedClusters ;

  // first we initialize the output == input 
  //
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	
	// et2x5, et5x5, is_ss, and is_looseTkss remain the same whether or not the cluster is stitched across GCT card boundaries
	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5 ; 

	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5 ;

	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss ;

	GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss ;
	GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss ;
      }}}

  // we will store new et in the GCTcombinedClusters, 0'ing lower clusters after stiching, dont need to care about other variabls they stay the 
  // same as input for now at least
  // we combine even phi boudaries positive eta, when combined the lowest et is set to 0
  
  for(int i=0; i<N_RCTCARDS_PHI-1; i=i+2){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	ap_uint<15> eta1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towEta*5+GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crEta ;
	ap_uint<15> phi1 = GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].crPhi ;
	// if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3 && phi1 == 4){
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3) {
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      // if( GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 && phi2 == 0) {
	      if( GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 ) {
		ap_uint<15> dPhi ; dPhi=((5 - phi1) + phi2) ;    
	        ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		// if( dEta < 2 ) {
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
		}}
	    }}
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
	// if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3 && phi1 == 4) {
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 3){ 
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      // if(GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0 && phi2 == 0) {
	      if(GCTcard.RCTcardEtaPos[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 0) {  
		ap_uint<15> dPhi ; dPhi=((5 - phi1) + phi2) ; 
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
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
		}}
	    }}
	}
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
	// if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0 && phi1 == 0 ) {
	if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0) {
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      // if(GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 && phi2 == 4 ) {
	      if(GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3) {
		ap_uint<15> dPhi ; dPhi=((5 - phi2) + phi1) ;  // reversed for negative eta
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
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
		    else { std::cout << "Energy of 'two' was insufficient for merge" << std::endl; }
		  }
		  else {
		    // Test: 
		    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (one > (0.10 * two)) {
		      std::cout << "merging.." << std::endl;
		      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = 0 ;
		      GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = one + two ; 
		    }
		    else { std::cout <<"Energy of 'one' was insufficient for merge" <<std::endl; 
		    }
		  }
		}}
	    }}
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
	
	// if(GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0 && phi1 == 0 ) {
	if(GCTcard.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].towPhi == 0){ 
	  for(int j1=0; j1<N_RCTGCT_FIBERS; j1++){
	    for(int k1=0; k1<N_RCTCLUSTERS_FIBER; k1++){
	      ap_uint<15> eta2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towEta*5+GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crEta ;
	      ap_uint<15> phi2 = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].crPhi ;
	      // if( GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 && phi2 == 4 ) {
	      if( GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].towPhi == 3 ) {
		ap_uint<15> dPhi ; dPhi=((5 - phi2) + phi1) ; // reversed compared to positive eta
		ap_uint<15> dEta ; dEta=(eta1 > eta2)?(eta1-eta2):(eta2-eta1) ;
		if( (dPhi <= 5) && (dEta < 2) ) {        
		  ap_uint<12> one = GCTcard.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et ;
		  ap_uint<12> two = GCTcard.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et ;
		  if (one > two){ 
		    // Test: only stitch if energy of one is >10% of the other                                                         
		    std::cout << "Comparing 'one' and 'two': " << one << ", " << two << std::endl;
                    if (two > (0.10 * one)) {
		      std::cout << "merging..." << std::endl;
		      GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et = one + two ;
		      GCTcombinedClusters.RCTcardEtaNeg[i+1].RCTtoGCTfiber[j1].RCTclusters[k1].et = 0 ; 
		    }
		    else { std::cout << "energy insufficient to merge" << std::endl;
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
		    else { std::cout <<"energy insufficient to merge" << std::endl;
		    }
		  }}
	      }}
	  }
	}
      }
    }
  }
  //
 
  // we need to store what we did before we start phi stiching
  //
  GCTcard_t GCTout ;
  for(int i=0; i<N_RCTCARDS_PHI; i++){
    for(int j=0; j<N_RCTGCT_FIBERS; j++){
      for(int k=0; k<N_RCTCLUSTERS_FIBER; k++){
	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et  ;
	
	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5  ;
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et2x5  ;

	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5  ;
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].et5x5  ;

	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss  ;
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_ss  ;
	
	GCTout.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss  = GCTcombinedClusters.RCTcardEtaPos[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss  ;
        GCTout.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss  = GCTcombinedClusters.RCTcardEtaNeg[i].RCTtoGCTfiber[j].RCTclusters[k].is_looseTkss  ;
      }}}
          
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
		}}
	    }}
	}
      }
    }
  }
  return GCTout ;
}


GCTinternal_t getClustersTowers(const GCTcard_t& GCTcard){

  GCTcard_t GCTcombinedClusters ;
  GCTinternal_t GCTout ;

  // here we will stitch the clusters in phi and eta
  GCTcombinedClusters = getClustersCombined(GCTcard) ;

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
  //
  //
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
    }}
  return GCTout ;
}

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
	          if( eta == eta1 && phi == phi1) GCTintTowers.GCTtower[eta][phi].et  = GCTintTowers.GCTtower[eta][phi].et + GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].et ;
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
	          if( eta == eta1 && phi == phi1) GCTintTowers.GCTtower[eta][phi].et  = GCTintTowers.GCTtower[eta][phi].et + GCTinternal.GCTCorrfiber[i+ic1].GCTclusters[jc].et ;
          }
        }
      }
    }
  }
  
  return GCTintTowers ;
}

//
// Compute isolation (sum of unclustered energy in 7x7 window IN TOWERS) for one cluster in place, where i is the fiber index 
// and k is the cluster-in-fiber index for the cluster. 
// nGCTCard is GCT card 0/1/2 (though it is not needed).
// 
void compute_isolation_for_one_cluster(GCTinternal_t& GCTinternal, int i, int k, int nGCTCard) {

  // We will only save clusters with > 0 GeV, so only need to do this for clusters with >0 energy 
  if (GCTinternal.GCTCorrfiber[i].GCTclusters[k].et == 0) {
    // iso for these clusters is zero
    // std::cout << "Cluster with zero energy: zero iso" << std::endl;
    GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso = 0;
    return;
  }
  
  std::cout << ">>> Calculating isolation..." << std::endl;
  ap_uint<12> uint_isolation = 0;

  /* std::cout << ">>> In iso calculation: GCTinternal.GCTCorrfiber[i].GCTclusters[k]: "  */
  /* 	    << "c.towEta: "  */
  /* 	    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta << ", "  */
  /* 	    << "c.crEta: "  */
  /* 	    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta << ", "  */
  /* 	    << "c.towPhi: "  */
  /* 	    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi << ", " */
  /* 	    << "c.crPhi: "  */
  /* 	    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi << std::endl; */

  bool getGlobal_iPhi = false;   // for the phi function: do not add the GCT card off-set, so we remain in the
  // gct local card iEta/iPhi
  int crystaliEta_in_GCT_card = getCluster_global_iEta(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
  int crystaliPhi_in_GCT_card = getCluster_global_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k], getGlobal_iPhi );
      
  int toweriEta_in_GCT_card = (int) (crystaliEta_in_GCT_card / 5);
  int toweriPhi_in_GCT_card = (int) (crystaliPhi_in_GCT_card / 5);
      
  std::cout << ">>> cluster's crystal ieta/iphi in GCT card: "
	    << crystaliEta_in_GCT_card << ", "
	    << crystaliPhi_in_GCT_card << ", "
	    << ">>> cluster's tower ieta/iphi in GCT card: " 
	    << toweriEta_in_GCT_card << ","
	    << toweriPhi_in_GCT_card << std::endl;

  // Is the cluster in a RCT card which overlaps with other GCT cards?
  bool inOverlapWithAnotherGCTCard = ( ((toweriPhi_in_GCT_card >= 0) && (toweriPhi_in_GCT_card < 4)) || ((toweriPhi_in_GCT_card >= 28) && (toweriPhi_in_GCT_card < 32)) );
      
  // If cluster is in the overlap region, do not compute isolation 
  if (inOverlapWithAnotherGCTCard) {
    GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso = 0;
    return;
  }

  // Size 7x7 in towers: include the overlap-region-between-GCT-cards-if-applicable
  // in eta direction, the min and max towers (inclusive!) are:
  int isoWindow_toweriEta_in_GCT_card_min = std::max(0, toweriEta_in_GCT_card - 3);
  int isoWindow_toweriEta_in_GCT_card_max = std::min(toweriEta_in_GCT_card + 3, N_GCTETA - 1);  // N_GCTETA = 34
  // e.g. if our window is centered at tower_iEta = 5, we want to sum towers_iEta 2, 3, 4, (5), 6, 7, 8, inclusive 
  // e.g. if our window is near the boundary, tower_iEta = 31, we want to sum towers_iEta 28, 29, 30, (31), 32, 33
  // inclusive (but there are only N_GCTETA = 34 towers, so we stop at tower_iEta = 33)
  
  // in phi direction, the min and max towers (inclusive!) are:
  int isoWindow_toweriPhi_in_GCT_card_min = std::max(0, toweriPhi_in_GCT_card - 3);
  int isoWindow_toweriPhi_in_GCT_card_max = std::min(toweriPhi_in_GCT_card + 3, N_GCTPHI - 1);  
  
  
  /* std::cout << ">>> window min/max eta: "  */
  /*   << isoWindow_toweriEta_in_GCT_card_min << ", " << isoWindow_toweriEta_in_GCT_card_max */
  /*   << std::endl; */
  /* std::cout << ">>> window min/max phi: "  */
  /*   << isoWindow_toweriPhi_in_GCT_card_min << ", " << isoWindow_toweriPhi_in_GCT_card_max  */
  /*   << std::endl; */
  
  // Sanity check: print all towers
  for (int iFiber = 0; iFiber < 64; iFiber++) {
    std::cout << "(fiber " << iFiber << "): ";
    for (int iTower = 0; iTower < 17; iTower++) {
      std::cout << GCTinternal.GCTCorrfiber[iFiber].GCTtowers[iTower].et << ", " ;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  // Translate this "tower index in GCT card" into something sensible for accessing the tower Ets

  // Keep track of the number of towers we summed over
  int nTowersSummed = 0;
  
  // For each tower we need....
  for (int iEta = isoWindow_toweriEta_in_GCT_card_min; 
       iEta <= isoWindow_toweriEta_in_GCT_card_max;
       iEta++) {
    
    for (int iPhi = isoWindow_toweriPhi_in_GCT_card_min;
	 iPhi <= isoWindow_toweriPhi_in_GCT_card_max;
	 iPhi++) {
      
      nTowersSummed++;
      
      // std::cout << "(iEta, iPhi): " << iEta << ", " << iPhi << std::endl;
      
      // Declare indices for accessing the tower Et
      int myIndexIntoGCT_64Fibers;
      int myIndexIntoGCT_Fiber_17Towers; 
      
      bool isTowerInPositiveEta = (iEta >= N_GCTTOWERS_FIBER); 
      if (isTowerInPositiveEta) { myIndexIntoGCT_64Fibers = iPhi; } // positive eta: phi index is simple
      // pos eta: e.g. if real phi = +80 degrees, iPhi in GCT = 31
      else                      { myIndexIntoGCT_64Fibers = (iPhi + N_GCTPOSITIVE_FIBERS); } // neg eta: add offset
      // neg eta: e.g. if real phi = +80 degrees, iPhi in GCT = 31, and my index into GCT fibers 31 + 32 = 63
      
      if (isTowerInPositiveEta) { myIndexIntoGCT_Fiber_17Towers = (iEta % 17); } // pos eta: if real eta = 1.47, iEta in GCT card = 33. If real eta = 0.0, iEta in GCT = 17, so iEta in fiber = 17%17 = 0.
      else                      { myIndexIntoGCT_Fiber_17Towers = (16 - iEta); }  // neg eta: if real eta = 0, iEta in GCT card = 16, i.e. our index into the GCT fiber is 16-16 = 0
      
      // std::cout << "... myIndexIntoGCT_64Fibers: " << myIndexIntoGCT_64Fibers << std::endl;
      // std::cout << "... myIndexIntoGCT_Fiber_17Towers: " << myIndexIntoGCT_Fiber_17Towers << std::endl;
      
      // Increment uint_isolation 
      ap_uint<12> myTowerEt = GCTinternal.GCTCorrfiber[myIndexIntoGCT_64Fibers].GCTtowers[myIndexIntoGCT_Fiber_17Towers].et;
      if (myTowerEt != 0) { 
	std::cout << "... myTowerEt (as float, non-zero): " << myTowerEt/8.0 << ", "
		  << "for 'tower-index-in-GCT-card' " << iEta << "," << iPhi << ", " 
		  << "and myindexes " << myIndexIntoGCT_64Fibers << "(fiber) and " 
		  << myIndexIntoGCT_Fiber_17Towers << "(tower-in-fiber)" 
		  << std::endl; 
	uint_isolation += myTowerEt;
      }
      
    }
  }
  
  // If we summed fewer than (7x7) = 49 towers because the cluster was at the edge of the permissible region,
  // scale up the isolation sum.
  float scaleFactor = ((float) (N_GCTTOWERS_CLUSTER_ISO_ONESIDE * N_GCTTOWERS_CLUSTER_ISO_ONESIDE) / (float) nTowersSummed);
  std::cout << "--> Summed over " << nTowersSummed << " towers: scaling iso " << uint_isolation 
	    << " by " << scaleFactor << " to get " << (uint_isolation * scaleFactor)
	    << std::endl;
  uint_isolation = (ap_uint<12>) (((float) uint_isolation) * scaleFactor);
  
  // Set the iso in the cluster
  GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso = uint_isolation;
  std::cout << "end of isolation calculation: (in GeV): " << uint_isolation / 8.0 
	    << ". Saved (uint) as: " <<  GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso
	    << std::endl;

}

// algo_top: First two arguments are the same as in the original firmware.
// nGCTCard is 0, 1, or 2 (needed for getting the cluster real eta/phis for CMSSW collections).
// gctClusters is the CMSSW-style collection of clusters, to be used in the analyzer.
// gctTowers is the CMSSW-style collection of towers, to be used in the analyzer.

void algo_top(const GCTcard_t& GCTcard, GCTtoCorr_t& GCTtoCorr,
	      unsigned int nGCTCard,
	      std::unique_ptr<l1tp2::CaloCrystalClusterCollection> const& gctClusters,
	      std::unique_ptr<l1tp2::CaloTowerCollection> const& gctTowers,
        std::unique_ptr<l1tp2::CaloTowerCollection> const& gctFullTowers) {
  
  GCTinternal_t GCTinternal ;

  // Fill the GCT area with proper et and eta/phi; the eta is 0....16 and towEtaNeg = 0 or 1
  // Phi 0....15 since this is half GCT card and 0...3 is overlap currently on one side, crEta/Phi 0...4
  GCTinternal = getClustersTowers(GCTcard) ;
  // Combine towers and clusters to get full towers
  GCTintTowers_t GCTintTowers;
  GCTintTowers = getFullTowers(GCTinternal) ;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Output to correlator, positive eta. Skip overlap region, i.e. fibers i = 0, 1, 2, 3, and i = 28, 29, 30, 31.
  ////////////////////////////////////////////////////////////////////////////////////////////////

  for(int i=4; i<(N_GCTPOSITIVE_FIBERS-N_RCTGCT_FIBERS); i++){
    for(int k=0; k<N_GCTCLUSTERS_FIBER; k++){

      // Comments from Stephanie:
      // Use GCTinternal for GCT clusters. 
      // Indexing: etas of towers go from 0-16 where iEta = 0 is real eta = 0 (indexes increase towards
      // larger abs(eta). Phis of towers go from 0-3 where iPhi = 0 is the leftmost (if you look
      // at a 'sideways' diagram of the GCT card like on the TWiki). 
      // 

      // Compute isolation 
      compute_isolation_for_one_cluster(GCTinternal, i, k, nGCTCard);
      // Compute isolation flags (inputs are: et (as float), and then relative isolation (iso/et))
      float relative_iso;
      if (GCTinternal.GCTCorrfiber[i].GCTclusters[k].et > 0) {
	        relative_iso = ( ((float) GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso / 8 ) / ( (float) GCTinternal.GCTCorrfiber[i].GCTclusters[k].et / 8 ) );
	        std::cout << "relative iso is " << (float) GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso/8 
		                << " divided by "     << (float) GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8
		                 << " giving " << relative_iso << std::endl;
      }
      else {
	// std::cout << "Cluster energy is 0: relative_iso set to 0" << std::endl;
	      relative_iso = 0;
      }
      bool is_iso        = passes_iso(        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0, relative_iso);
      bool is_looseTkiso = passes_looseTkiso( GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0, relative_iso);
      GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_iso        = is_iso;
      GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkiso = is_looseTkiso;

      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].et  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et   ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towEtaNeg  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg  ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].towPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi-4 ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].crEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].crPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].iso    = GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso ;   // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].et2x5  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5 ; // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].et5x5  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5 ; // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].is_ss  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss ; // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].is_looseTkss = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss ; // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].is_iso = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_iso ; // new
      GCTtoCorr.GCTCorrfiber[i-4].GCTclusters[k].is_looseTkiso = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkiso; // new
      // Get the real eta, phi using two helper functions
      int crystaliEta_in_barrel = getCluster_global_iEta(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
      int crystaliPhi_in_barrel = getCluster_global_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
      float realEta = getEta_fromCrystaliEta(crystaliEta_in_barrel);
      float realPhi = getPhi_fromCrystaliPhi(crystaliPhi_in_barrel);

      reco::Candidate::PolarLorentzVector p4cluster(GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0,
						    realEta,
						    realPhi,
						    0.);
      l1tp2::CaloCrystalCluster cluster(p4cluster, 
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0,   // convert to float
					0,  // float h over e                              
					relative_iso,		   // for consistency with the old emulator, in this field save (iso energy sum)/(cluster energy)
                                        0,  // DetId seedCrystal                              
                                        0,  // puCorrPt                                           
                                        0,  // 0, 1, or 2 (as computed in firmware)                
                                        0,  // et2x2 (not calculated)                             
                                        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5/8.0,  // et2x5 (as computed in firmware, save float)           
                                        0,  // et3x5 (not calculated)                             
                                        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5/8.0,   // et5x5 (as computed in firmware, save float)  
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss,  // standalone WP: not computed
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss, // electronWP98: not computed 
					false, // is_photon in Cecile's emulator, photonWP80: not computed
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss, // electronWP90: not computed
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss, // looseL1TkMatchWP
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss  // stage2effMatch: not computed
                                        );

      // Flags
      std::map<std::string, float> params;
      params["standaloneWP_showerShape"] = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss;
      params["standaloneWP_isolation"]   = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_iso;
      params["trkMatchWP_showerShape"]   = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss;
      params["trkMatchWP_isolation"]     = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkiso;
      cluster.setExperimentalParams(params);

      if (cluster.pt() > 0.0) {
	      gctClusters->push_back(cluster);
	      std::cout << "--- cluster pT, global iEta, iPhi and real eta, phi: "
		              << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0 
		              << ", ("
		              << crystaliEta_in_barrel
		              << ", "
	                << crystaliPhi_in_barrel
		              << "), ("
		              << realEta 
		              << "," 
		              << realPhi
		              << ")"
		              << " with relative isolation " << relative_iso
	                << std::endl;
	      std::cout << "    with the GCTinternal values: " 
		  << "towEtaNeg: " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg << ", "
		  << "towEta: "    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta    << ", "
		  << "towPhi: "    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi    << ", "
		  << "crEta: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta   << ", "
		  << "crPhi: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi   << ", "
	          << "iso: "     << GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso/8.0 << ", " 
		  << "et2x5: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5/8.0 << ", "
		  << "et5x5: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5/8.0 << ", "
		  << "is_ss: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss << ", "
		  << "is_looseTkss" << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss << ", "
		  << "is_iso: " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_iso << ", "
		  << "is_looseTkiso: " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkiso << std::endl;
	
	
      }

    }
    // Positive eta towers : push back to CMSSW collection
    for(int k=0; k<N_GCTTOWERS_FIBER; k++){
      std::cout<< "Accessing positive eta: GCTCorrfiber " << i-4
      	       << " , GCTtowers " << k
      	       << " , energy " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].et
	       << " , hoe " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].hoe 
	       << " , ecalEt " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].ecalEt
               << " , hcalEt " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].hcalEt
	       << std::endl;
      GCTtoCorr.GCTCorrfiber[i-4].GCTtowers[k].et     = GCTinternal.GCTCorrfiber[i].GCTtowers[k].et ;
      GCTtoCorr.GCTCorrfiber[i-4].GCTtowers[k].hoe    = GCTinternal.GCTCorrfiber[i].GCTtowers[k].hoe ; 
      GCTtoCorr.GCTCorrfiber[i-4].GCTtowers[k].ecalEt = GCTinternal.GCTCorrfiber[i].GCTtowers[k].ecalEt ; 
      GCTtoCorr.GCTCorrfiber[i-4].GCTtowers[k].hcalEt = GCTinternal.GCTCorrfiber[i].GCTtowers[k].hcalEt ; 

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
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Output to correlator: In negative eta, the overlap region to skip is fibers 32, 33, 34, 35, and 61, 62, 63, 64.
  //////////////////////////////////////////////////////////////////////////////////////////////// 
  for(int i=(N_GCTPOSITIVE_FIBERS+N_RCTGCT_FIBERS); i<(N_GCTINTERNAL_FIBERS-N_RCTGCT_FIBERS); i++){
    for(int k=0; k<N_GCTCLUSTERS_FIBER; k++){
      /* std::cout << "Accessing negative eta: GCTCorrfiber " << i-12 */
      /*           << " , GCTclusters " << k  */
      /*           << " , energy " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et << std::endl; */
      // Isolation
      compute_isolation_for_one_cluster(GCTinternal, i, k, nGCTCard);

      // Isolation flags (inputs are: et (as float), and then relative isolation (iso/et))                                                               
      float relative_iso;
      if (GCTinternal.GCTCorrfiber[i].GCTclusters[k].et > 0) {
        relative_iso = ( ((float) GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso / 8 ) / ( (float) GCTinternal.GCTCorrfiber[i].GCTclusters[k].et / 8 ) );
	std::cout << "relative iso is "<< (float) GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso/8
                  << " divided by "     << (float) GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8
                  << " giving " << relative_iso << std::endl;
      }
      else {
	      // std::cout << "Cluster energy is 0: relative_iso set to 0" << std::endl;
        relative_iso = 0;
      }
      // Build negative eta clusters
      bool is_iso        = passes_iso(        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0, relative_iso);
      bool is_looseTkiso = passes_looseTkiso( GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0, relative_iso);
      GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_iso        = is_iso;
      GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkiso = is_looseTkiso;

      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].et  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et   ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].towEtaNeg  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg  ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].towEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].towPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi-4 ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].crEta  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].crPhi  =  GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].iso    = GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso;    // new
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].et2x5  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5 ; // new                            
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].et5x5  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5 ; // new                         
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].is_ss  = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss ; // new                
      GCTtoCorr.GCTCorrfiber[i-12].GCTclusters[k].is_looseTkss = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss ; // new    

      // Get the real eta, phi using two helper functions
      int globaliEta = getCluster_global_iEta(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
      int globaliPhi = getCluster_global_iPhi(nGCTCard, GCTinternal.GCTCorrfiber[i].GCTclusters[k]);
      float realEta = getEta_fromCrystaliEta(globaliEta);
      float realPhi = getPhi_fromCrystaliPhi(globaliPhi);
      
      reco::Candidate::PolarLorentzVector p4cluster(GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0,
						    realEta,
						    realPhi,
						    0.);
      l1tp2::CaloCrystalCluster cluster(p4cluster, 
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0, // conver to float
					0,  // float h over e                              
					relative_iso, // follow old emulator's convention: save relative iso in this field
                                        0,  // DetId seedCrystal                              
                                        0,  // puCorrPt                                           
                                        0,  // 0, 1, or 2 (as computed in firmware)                
                                        0,  // et2x2 (not calculated)                             
                                        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5/8.0,  // et2x5 (as computed in firmware, save float)           
                                        0,  // et3x5 (not calculated)                             
                                        GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5/8.0,   // et5x5 (as computed in firmware, save float)  
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss,  // standalone WP
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss, // electronWP98: not computed 
					false, // photonWP80: not computed
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss, // electronWP90: not computed
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss, // looseL1TkMatchWP
					GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss  // stage2effMatch: not computed
                                        );

      // Experimental parameters
      std::map<std::string, float> params;
      params["standaloneWP_showerShape"] = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss;
      params["standaloneWP_isolation"]   = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_iso;
      params["trkMatchWP_showerShape"]   = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss;
      params["trkMatchWP_isolation"]     = GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkiso;
      cluster.setExperimentalParams(params);

      // Push back negative eta clusters to the output CMSSW collection
      if (cluster.pt() > 0.0) {
	        gctClusters->push_back(cluster);
	        std::cout << "--- cluster pT, global iEta, iPhi and real eta, phi: "
		  << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et/8.0 
		  << ", ("
		  << globaliEta 
		  << ", "
		  << globaliPhi
		  << "), ("
		  << realEta 
		  << "," 
		  << realPhi
		  << ")"
		  << " with relative isolation " 
		  << relative_iso
		  << std::endl;
	std::cout << "    with the GCTinternal values: "
                  << "towEtaNeg: " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEtaNeg << ", "
                  << "towEta: "    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towEta    << ", "
                  << "towPhi: "    << GCTinternal.GCTCorrfiber[i].GCTclusters[k].towPhi    << ", "
                  << "crEta: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crEta   << ", "
                  << "crPhi: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].crPhi   << ", "
		  << "iso: "     << GCTinternal.GCTCorrfiber[i].GCTclusters[k].iso   << ", "
		  << "et2x5: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et2x5 << ", "
                  << "et5x5: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].et5x5 << ", "
                  << "is_ss: "   << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_ss << ", "
                  << "is_looseTkss: " << GCTinternal.GCTCorrfiber[i].GCTclusters[k].is_looseTkss << std::endl;
	
      }
      
    }
    // Negative eta towers : push back to CMSSW as well 
    for(int k=0; k<N_GCTTOWERS_FIBER; k++){
      std::cout<< "Accessing negative eta: GCTCorrfiber " << i
      	       << " , GCTtowers " << k
      	       << " , energy " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].et 
	       << " , hoe " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].hoe 
	       << " , ecalEt " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].ecalEt 
	       << " , hcalEt " << GCTinternal.GCTCorrfiber[i].GCTtowers[k].hcalEt 
	       << std::endl;
      GCTtoCorr.GCTCorrfiber[i-12].GCTtowers[k].et     = GCTinternal.GCTCorrfiber[i].GCTtowers[k].et ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTtowers[k].hoe    = GCTinternal.GCTCorrfiber[i].GCTtowers[k].hoe ;
      GCTtoCorr.GCTCorrfiber[i-12].GCTtowers[k].ecalEt = GCTinternal.GCTCorrfiber[i].GCTtowers[k].ecalEt ; 
      GCTtoCorr.GCTCorrfiber[i-12].GCTtowers[k].hcalEt = GCTinternal.GCTCorrfiber[i].GCTtowers[k].hcalEt ; 

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

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // Create GCT Full Towers (clusters + towers) output to CMSSW.
  ////////////////////////////////////////////////////////////////////////////////////////////////
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
