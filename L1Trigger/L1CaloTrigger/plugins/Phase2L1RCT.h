
#ifndef _PHASE_2_L1_RCT_H_
#define _PHASE_2_L1_RCT_H_

static constexpr int n_towers_Eta = 34;
static constexpr int n_towers_Phi = 72;
static constexpr int n_towers_halfPhi = 36;
static constexpr int n_towers_cardEta = 17;   // new: equivalent to n_towers_per_link
static constexpr int n_towers_cardPhi = 4;    
static constexpr int n_crystals_cardEta = (n_towers_Eta * n_towers_cardEta);
static constexpr int n_crystals_cardPhi = (n_towers_Phi * n_towers_cardPhi);

// outputs
static constexpr int n_links_card = 4;        // 4 links per card
static constexpr int n_clusters_link = 2;     // 2 clusters sent out in each link
static constexpr int n_clusters_4link = 8;    // 8 clusters sent out in 4 links
static constexpr int n_towers_per_link = 17;  

static constexpr int CRYSTALS_IN_TOWER_ETA = 5;
static constexpr int CRYSTALS_IN_TOWER_PHI = 5;

static constexpr int TOWER_IN_ETA = 3;      // number of towers in eta, in one 3x4 region (barrel)
static constexpr int TOWER_IN_PHI = 4;      // number of towers in phi, in one 3x4 region (barrel)

//static constexpr int TOWER_IN_ETA_OVERLAP = 2; // number of towers in eta, in one 2x4 region (overlap)
//static constexpr int TOWER_IN_PHI_OVERLAP = 4; // number of towers in phi, in one 2x4 region (overlap)

static constexpr int CRYSTAL_IN_ETA = 15;   // number of crystals in eta, in one 3x4 region (barrel)
static constexpr int CRYSTAL_IN_PHI = 20;   // number of crystals in phi, in one 3x4 region (barrel)

static constexpr float ECAL_eta_range = 1.4841;
static constexpr float half_crystal_size = 0.00873;

static constexpr float slideIsoPtThreshold = 80;
static constexpr float a0_80 = 0.85, a1_80 = 0.0080, a0 = 0.21;                        // passes_iso
static constexpr float b0 = 0.38, b1 = 1.9, b2 = 0.05;                                 // passes_looseTkiso
static constexpr float c0_ss = 0.94, c1_ss = 0.052, c2_ss = 0.044;                     // passes_ss
static constexpr float d0 = 0.96, d1 = 0.0003;                                         // passes_photon
static constexpr float e0_looseTkss = 0.944, e1_looseTkss = 0.65, e2_looseTkss = 0.4;  // passes_looseTkss
static constexpr float cut_500_MeV = 0.5;

static constexpr int N_CLUSTERS_PER_REGION = 4;       // number of clusters per ECAL region
static constexpr int N_REGIONS_PER_CARD = 6;          // number of ECAL regions per card

// absolute IDs range from 0-33
// 0-16 are iEta -17 to -1
// 17 to 33 are iEta 1 to 17
static constexpr int toweriEta_fromAbsoluteID_shift = 16;

// absolute IDs range from 0-71.
// To align with detector tower IDs (1 - n_towers_Phi)
// shift all indices by 37 and loop over after 72
static constexpr int toweriPhi_fromAbsoluteID_shift_lowerHalf = 37;
static constexpr int toweriPhi_fromAbsoluteID_shift_upperHalf = 35;


// Convert HCAL ET to ECAL ET convention 
ap_uint<12> convertHcalETtoEcalET(ap_uint<12> HCAL) {
  float hcalLSB = 0.5;
  float hcalEtAsFloat = HCAL * hcalLSB;
  return (ap_uint<12> (hcalEtAsFloat * 8));
}

////////////////////////////////////////////////////////////////////////// 
// RCT Cards: indexing helper functions
////////////////////////////////////////////////////////////////////////// 

// Assert that the card index is within bounds. (Valid cc: 0 to 35, since there are 36 RCT cards)
bool isValidCard(int cc) {
  return ((cc > -1) && (cc < 36));
}

//////////////////////////////////////////////////////////////////////////
// RCT Cards: need to know their min/max crystal boundaries.
//////////////////////////////////////////////////////////////////////////
// For a card (ranging from 0 to 35, since there are 36 cards), return the iEta of the crystal with max iEta.
// This represents the card boundaries in eta (identical to getEtaMax_card in the original emulator)
int getCard_iEtaMax(int cc) {
  assert(isValidCard(cc));
  
  int etamax = 0;
  if (cc % 2 == 0)   // Even card: negative eta
    etamax = (n_towers_cardEta * CRYSTALS_IN_TOWER_ETA - 1);  // First eta half. 5 crystals in eta in 1 tower.
  else               // Odd card: positive eta
    etamax = (n_towers_Eta * CRYSTALS_IN_TOWER_ETA - 1);
  return etamax;
}

// Same as above but for minimum iEta.
int getCard_iEtaMin(int cc) {
  int etamin = 0;
  if (cc % 2 == 0)  // Even card: negative eta 
    etamin = (0);  
  else                // Odd card: positive eta
    etamin = (n_towers_cardEta * CRYSTALS_IN_TOWER_ETA);
  return etamin;
}

// Same as above but for maximum iPhi. 
int getCard_iPhiMax(int cc) {
  int phimax = ((cc / 2) + 1) * 4 * CRYSTALS_IN_TOWER_PHI - 1;
  return phimax;
}

// Same as above but for minimum iPhi.
int getCard_iPhiMin(int cc) {
  int phimin = (cc / 2) * 4 * CRYSTALS_IN_TOWER_PHI;
  return phimin;
}

// Given the RCT card number (0-35), get the crystal iEta of the "bottom left" corner
int getCard_refCrystal_iEta(int cc) {

  if ((cc % 2) == 1) {  // if cc is odd (positive eta)
    return (17 * CRYSTALS_IN_TOWER_ETA);
  }
  else {   // if cc is even (negative eta) the bottom left corner is further in eta, hence +4
    return ((16 * CRYSTALS_IN_TOWER_ETA) + 4);
  }
}

// Given the RCT card number (0-35), get the global crystal iPhi of the "bottom left" corner (0- 71*5)
int getCard_refCrystal_iPhi(int cc) {

  if ((cc % 2) == 1) { 
    // if cc is odd: positive eta   
    return int(cc / 2) * TOWER_IN_PHI * CRYSTALS_IN_TOWER_PHI; 
  } 
  else {  
    // if cc is even, the bottom left corner is further in phi, hence the +4 and -1
    return (((int(cc / 2) * TOWER_IN_PHI) + 4) * CRYSTALS_IN_TOWER_PHI) - 1;
  }
}

////////////////////////////////////////////////////////////////////////// 
// Towers: Go from real (eta, phi) to tower absolute ID
////////////////////////////////////////////////////////////////////////// 

/* 
 * For a real eta, get the tower absolute Eta index (possible values are 0-33, since there
 * are 34 towers in eta. (Adapted from getTower_absoluteEtaID)
 */
int getTower_absEtaID(float eta) {
  float size_cell = 2 * ECAL_eta_range / n_towers_Eta;
  int etaID = int((eta + ECAL_eta_range) / size_cell);
  return etaID;
}

/* 
 * Same as above, but for phi.
 * Possible values range from 0-71 (Adapted from getTower_absolutePhiID)
 */
int getTower_absPhiID(float phi) {
  float size_cell = 2 * M_PI / n_towers_Phi;
  int phiID = int((phi + M_PI) / size_cell);
  return phiID;
}

////////////////////////////////////////////////////////////////////////// 
// Towers: Go from firmware specifics (RCT card, tower number in link, and link number (all firmware convention))
//         to tower absolute ID.
////////////////////////////////////////////////////////////////////////// 

/* 
 * Get the global tower iEta (0-31) from the firmware card, tower number (0-16), and link (0-3). Respects the fact that
 * in the firmware, negative eta cards are "rotated" (link 0, tower 0) starts in the "top right" corner if we
 * look at a diagram of the barrel region.
 */
int getAbsID_iEta_fromFirmwareCardTowerLink(int nCard, int nTower, int nLink) {
  // iEta only depends on the tower position in the link
  (void) nCard;
  (void) nLink; 
  if ((nCard % 2) == 1) { // if cc is odd (positive eta), e.g. nTower = 0 will correspond to absolute iEta ID 17.
    return n_towers_per_link + nTower;
  }
  else { // if cc is even (negative eta): e.g. nTower = 0 will correspond to absolute iEta ID 16.
    return (16-nTower);
  }
}

/*
 * Get the global tower iPhi (0-71) from the firmware card, tower number (0-16), and link (0-3).
 */
int getAbsID_iPhi_fromFirmwareCardTowerLink(int nCard, int nTower, int nLink) {
  // iPhi only depends on the card and link number
  (void) nTower;
  if ((nCard % 2) == 1) { // if cc is odd (positive eta),
    // e.g. cc=3, link #2, global iPhi = int(3/2) * 4 + 2 = 6
    return (int(nCard/2) * TOWER_IN_PHI) + nLink;
  }
  else { // if cc is even (negative eta)
    // e.g. cc=4, link #2, global iPhi = int(4/2) * 4 + (4 - 2 - 1)
    //                                 = 2*4 + 1
    //                                 = 9
    // minus one is because TOWER_IN_PHI is 4
    return (int(nCard/2) * TOWER_IN_PHI) + (TOWER_IN_PHI - nLink - 1);
  }
}

////////////////////////////////////////////////////////////////////////// 
// Towers: Go from absolute ID, back to real eta and phi.
////////////////////////////////////////////////////////////////////////// 

/*
 * From the tower absolute ID in eta (0-33), get the real eta of the tower center
 * Same as getTowerEta_fromAbsoluteID in previous CMSSW emulator
 */
float getTowerEta_fromAbsID(int id) {
  float size_cell = 2 * ECAL_eta_range / n_towers_Eta;
  float eta = (id * size_cell) - ECAL_eta_range + 0.5 * size_cell;
  return eta;
}

/*
 * From the tower absolute ID in phi (0-71), get the real phi of the tower center
 * Same as getTowerPhi_fromAbsoluteID in previous CMSSW emulator
 */
float getTowerPhi_fromAbsID(int id) {
  float size_cell = 2 * M_PI / n_towers_Phi;
  float phi = (id * size_cell) - M_PI + 0.5 * size_cell;
  return phi;
}

//////////////////////////////////////////////////////////////////////////
// Crystals: go from real eta and phi to crystal iEta/iPhi in detector.
//////////////////////////////////////////////////////////////////////////
/* 
 * Get crystal's global iEta from real eta. (identical to getCrystal_etaID in L1EGammaCrystalsEmulatorProducer.cc)
 * This "global" iEta ranges from 0 to (33*5) since there are 34 towers in eta in the full detector, 
 * each with five crystals in eta.
 */
int getCrystal_iEta(float eta) {
  float size_cell = 2 * ECAL_eta_range / (CRYSTALS_IN_TOWER_ETA * n_towers_Eta);
  int iEta = int((eta + ECAL_eta_range) / size_cell);
  return iEta;
}

/*
 * Get crystal's global iPhi from real phi. (identical to getCrystal_phiID in L1EGammaCrystalsEmulatorProducer.cc)
 * This "global" iPhi ranges from 0 to (71*5) since there are 72 towers in phi in the full detector,
 * each with five crystals in eta.
 */
int getCrystal_iPhi(float phi) {
  float size_cell = 2 * M_PI / (CRYSTALS_IN_TOWER_PHI * n_towers_Phi);
  int iPhi = int((phi + M_PI) / size_cell);
  return iPhi;
}

////////////////////////////////////////////////////////////////////////// 
// Crystals: going from real eta/phi to firmware conventions, given the RCT card.
////////////////////////////////////////////////////////////////////////// 

/*
 * For a crystal with real eta, and falling in card cc, get its local iEta 
 * relative to the bottom left corner of the card (possible local iEta ranges from 0 to 17 * 5,
 * since in one card, there are 17 towers in eta, each with 5 crystals in eta.
 */
int getCrystal_local_iEta(float hitEta, int cc) {
  return abs(getCard_refCrystal_iEta(cc) - getCrystal_iEta(hitEta));
}

/*
 * Same as above, but for iPhi (possible local iPhi ranges from 0 to (3*5), since in one card,
 * there are 4 towers in phi, each with 5 crystals in phi.
 */
int getCrystal_local_iPhi(float hitPhi, int cc) {
  return abs(getCard_refCrystal_iPhi(cc) - getCrystal_iPhi(hitPhi));
}

/* 
 * Get the RCT card region that a crystal is in, given the "local" iEta of the crystal 
 * 0 is region closest to eta = 0. Regions 0, 1, 2, 3, 4 are in the barrel, Region 5 is in overlap
 */
int getRegionNumber(int local_iEta) {
  int no = int(local_iEta / (TOWER_IN_ETA * CRYSTALS_IN_TOWER_ETA));
  assert(no < 6); 
  return no;
}

//////////////////////////////////////////////////////////////////////////
// Clusters: going from firmware conventions to global iEta/iPhi.
//////////////////////////////////////////////////////////////////////////
/* 
 * From a crystal in card cc, region regionIdx, tower number (0-2 for a 3x4 region) and crystal number (in the tower)
 * return the crystal iEta (0-33*5) in the entire card. Use this with getEta_fromCrystal_iEta to go from a cluster
 * to the real eta of its position.
 */
int getCrystal_iEta_fromCardRegionInfo(int cc, int regionIdx, int towerEta_in_region, int crystalEta_in_tower) {
  
  int crystalEta_in_card = ((regionIdx * TOWER_IN_ETA * CRYSTALS_IN_TOWER_ETA) + (towerEta_in_region * CRYSTALS_IN_TOWER_ETA)
			    + crystalEta_in_tower);
  if ((cc % 2) == 1) {
    // if card is odd (positive eta)
    return (getCard_refCrystal_iEta(cc) + crystalEta_in_card);
  }
  else {
    // if card is even (negative eta)
    return (getCard_refCrystal_iEta(cc) - crystalEta_in_card);
  }
}

/* 
 * Same thing as above but in phi.
 */
int getCrystal_iPhi_fromCardRegionInfo(int cc, int regionIdx, int towerPhi_in_region, int crystalPhi_in_tower) {
  
  int crystalPhi_in_card = (towerPhi_in_region * CRYSTALS_IN_TOWER_PHI) + crystalPhi_in_tower;
  if ((cc % 2) == 1) {
    // if card is odd (positive eta)   
    return (getCard_refCrystal_iPhi(cc) + crystalPhi_in_card);
  }
  else {
    // if card is even (negative eta)
    return (getCard_refCrystal_iPhi(cc) - crystalPhi_in_card);
  }
}

//////////////////////////////////////////////////////////////////////////
// Crystals: go from crystal global iEta/iPhi in detector back to real eta, phi.
//////////////////////////////////////////////////////////////////////////

/* 
 * Get crystal's real eta from global iEta (0-33*5).
 */
float getEta_fromCrystaliEta(int iEta) {
  float size_cell = 2 * ECAL_eta_range / (CRYSTALS_IN_TOWER_ETA * n_towers_Eta);
  return iEta * size_cell - ECAL_eta_range + half_crystal_size;
}

/* 
 * Get crystal's real eta from global iPhi (0-71*5).
 */
float getPhi_fromCrystaliPhi(int iPhi) {
  float size_cell = 2 * M_PI / (CRYSTALS_IN_TOWER_PHI * n_towers_Phi);
  return iPhi * size_cell - M_PI + half_crystal_size;
}


////////////////////////////////////////////////////////////////////////// 
// Cluster flags: taken from previous emulator
//  https://github.com/skkwan/cmssw/blob/devel-Phase2RCTCluster/L1Trigger/L1CaloTrigger/plugins/L1EGammaCrystalsEmulatorProducer.cc#L1240-L1280
////////////////////////////////////////////////////////////////////////// 

bool passes_iso(float pt, float iso) {
  bool is_iso = true;
  if (pt < slideIsoPtThreshold) {
    if (!((a0_80 - a1_80 * pt) > iso))
      is_iso = false;
  } else {
    if (iso > a0)
      is_iso = false;
  }
  if (pt > 130)
    is_iso = true;
  return is_iso;
}

bool passes_looseTkiso(float pt, float iso) {
  bool is_iso = (b0 + b1 * std::exp(-b2 * pt) > iso);
  if (pt > 130)
    is_iso = true;
  return is_iso;
}

bool passes_ss(float pt, float ss) {
  bool is_ss = ((c0_ss + c1_ss * std::exp(-c2_ss * pt)) <= ss);
  if (pt > 130)
    is_ss = true;
  return is_ss;
}

bool passes_looseTkss(float pt, float ss) {
  bool is_ss = ((e0_looseTkss - e1_looseTkss * std::exp(-e2_looseTkss * pt)) <= ss);
  if (pt > 130)
    is_ss = true;
  return is_ss;
}



//////////////////////////////////////////////////////////////////////////  
// Other emulator helper functions
//////////////////////////////////////////////////////////////////////////

/*
 * Declare the SimpleCaloHit class (taken from previous emulator)
 */

class SimpleCaloHit {
 private:
  float pt_ = 0;
  float energy_ = 0.;
  ap_uint<10> et_uint_;
  GlobalVector position_;  // As opposed to GlobalPoint, so we can add them (for weighted average)
  HcalDetId id_hcal_;
  EBDetId id_;
  
 public:
  // tool functions
  inline void setPt() { pt_ = (position_.mag2() > 0) ? energy_ * sin(position_.theta()) : 0; };
  inline void setEnergy(float et) { energy_ = et / sin(position_.theta()); };
  inline void setEt_uint(ap_uint<10> et_uint) { et_uint_ = et_uint; }
  inline void setPosition(const GlobalVector& pos) { position_ = pos; };
  inline void setIdHcal(const HcalDetId& idhcal) { id_hcal_ = idhcal; };
  inline void setId(const EBDetId& id) { id_ = id; };

  inline float pt() const { return pt_; };
  inline float energy() const { return energy_; };
  inline ap_uint<10> et_uint() const { return et_uint_; };
  inline const GlobalVector& position() const { return position_; };
  inline const EBDetId& id() const { return id_; };
};

/*******************************************************************/

/*
 * linkECAL class: represents one ECAL link (one tower: 5x5 crystals)
 */

class linkECAL {
 private:
  ap_uint<10> crystalE[CRYSTALS_IN_TOWER_ETA][CRYSTALS_IN_TOWER_PHI];

 public:
  // constructor                                                                                               
  linkECAL() { crystalE[CRYSTALS_IN_TOWER_ETA][CRYSTALS_IN_TOWER_PHI] = {};  }

  // Set members
  inline void zeroOut() {  // zero out the crystalE array
    for (int i = 0; i < CRYSTALS_IN_TOWER_ETA; i++) {
      for (int j = 0; j < CRYSTALS_IN_TOWER_PHI; j++) { 
	crystalE[i][j] = 0;
      }}};
  inline void setCrystalE(int iEta, int iPhi, ap_uint<10> energy) {
    assert(iEta < CRYSTALS_IN_TOWER_ETA); assert(iPhi < CRYSTALS_IN_TOWER_PHI); crystalE[iEta][iPhi] = energy; };
  inline void addCrystalE(int iEta, int iPhi, ap_uint<10> energy) { 
    assert(iEta < CRYSTALS_IN_TOWER_ETA); assert(iPhi < CRYSTALS_IN_TOWER_PHI); crystalE[iEta][iPhi] += energy; };

  // Access members
  inline ap_uint<10> getCrystalE(int iEta, int iPhi) { assert(iEta < 5); assert(iPhi < 5); return crystalE[iEta][iPhi]; };

};

/*******************************************************************/

/*
 * region3x4 class: represents one 3x4 ECAL region. The region stores no
 *                  information about which card it is located in.
 *                  idx: 0-4. Region 0 is the one closest to eta = 0, counting outwards in eta  
 */

class region3x4 {
 private:
  int idx_ = -1; 
  linkECAL linksECAL[TOWER_IN_ETA][TOWER_IN_PHI];  // 3x4 in towers

 public:
  // constructor                                                                            
  region3x4() { idx_ = -1; }

  // copy constructor
  region3x4(const region3x4& other) {
    idx_ = other.idx_;
    for (int i = 0; i < TOWER_IN_ETA; i++) {
      for (int j = 0; j < TOWER_IN_PHI; j++ ) {
        linksECAL[i][j] = other.linksECAL[i][j];
      }
    }
  }
  
  // overload operator= to use copy constructor                                                            
  region3x4 operator=(const region3x4& other) {
    region3x4 newRegion(other);
    return newRegion;
  };

  // set members
  inline void zeroOut() { 
    for (int i = 0; i < TOWER_IN_ETA; i++) { 
      for (int j = 0; j < TOWER_IN_PHI; j++) { 
	linksECAL[i][j].zeroOut();
      }}
  };
  inline void setIdx(int idx) { idx_ = idx; };

  // get members
  inline float getIdx() const { return idx_; };
  inline linkECAL& getLinkECAL (int iEta, int iPhi) { return linksECAL[iEta][iPhi]; };
};

/*******************************************************************/

/*
 * towerHCAL class: represents one HCAL tower
 */

class towerHCAL {
 private:
  ap_uint<10> et;
  ap_uint<6> fb;

 public:
  // constructor
  towerHCAL() { et = 0;  fb = 0; };

  // copy constructor
  towerHCAL(const towerHCAL &other) {
    et = other.et; fb = other.fb;   };

  // set members
  inline void zeroOut() { et = 0; fb = 0; };
  inline void addEt(ap_uint<10> newEt) { et += newEt; };

  // get members
  inline ap_uint<10> getEt() { return et; };

};

/*******************************************************************/

/*
 * towers3x4 class: represents 3x4 array of HCAL towers. idx = 0, 1, ... 4 are the barrel gion
 */

class towers3x4 {
 private:
  int idx_ = -1; 
  towerHCAL towersHCAL[TOWER_IN_ETA][TOWER_IN_PHI];  // 3x4 in towers

 public:
  // constructor                                                                            
  towers3x4() { idx_ = -1; };

  // copy constructor
  towers3x4(const towers3x4& other) {
    idx_ = other.idx_;
    for (int i = 0; i < TOWER_IN_ETA; i++) {
      for (int j = 0; j < TOWER_IN_PHI; j++ ) {
	towersHCAL[i][j] = other.towersHCAL[i][j]; }}; };

  // set members
  inline void zeroOut() { 
    for (int i = 0; i < TOWER_IN_ETA; i++) { 
      for (int j = 0; j < TOWER_IN_PHI; j++) { 
	towersHCAL[i][j].zeroOut();      }}   };
  inline void setIdx(int idx) { idx_ = idx; };

  // get members
  inline float getIdx() const { return idx_; };
  inline towerHCAL& getTowerHCAL ( int iEta, int iPhi ){ return towersHCAL[iEta][iPhi]; }; 
  
};

/*******************************************************************/

/* 
 * card class: represents one RCT card. Each card has five 3x4 regions and one 2x4 region,
 *             which is represented by a 3x4 region with its third row zero'd out.
 *             idx 0-35: odd values of cardIdx span eta = 0 to eta = 1.41 
 *                       even values of cardIdx span eta = -1.41 to eta = 0
 *             The realEta and realPhi arrays store the (eta, phi) of the center of the towers.
 */

class card {
 private:
  int idx_ = -1 ; 
  region3x4 card3x4Regions[N_REGIONS_PER_CARD];
  towers3x4 card3x4Towers[N_REGIONS_PER_CARD];
  
 public:
  // constructor
  card() {  
    idx_ = -1; 
    for (int i = 0; i < N_REGIONS_PER_CARD; i++) {
      card3x4Regions[i].setIdx(i);  
      card3x4Regions[i].zeroOut();
      card3x4Towers[i].setIdx(i);
      card3x4Towers[i].zeroOut();
    }};
  
  // copy constructor
  card(const card& other) {
    idx_ = other.idx_;
    for (int i = 0; i < N_REGIONS_PER_CARD; i++) {
      card3x4Regions[i] = other.card3x4Regions[i]; 
      card3x4Towers[i]  = other.card3x4Towers[i];  }};

  // overload operator= to use copy constructor
  card operator=(const card& other) {
    card newCard(other);
    return newCard;  };
  
  // set members
  inline void setIdx(int idx) { idx_ = idx; };
  inline void zeroOut() { 
    for (int i = 0; i < N_REGIONS_PER_CARD; i++) { card3x4Regions[i].zeroOut(); card3x4Towers[i].zeroOut();    }; };

  // get members
  inline float getIdx() const { return idx_; };
  inline region3x4& getRegion3x4(int idx) { assert(idx < N_REGIONS_PER_CARD); return card3x4Regions[idx]; }
  inline towers3x4& getTowers3x4(int idx) { assert(idx < N_REGIONS_PER_CARD); return card3x4Towers[idx]; }

};

/*******************************************************************/

/*
 *  crystal class: 
 */

class crystal{
 public:
  ap_uint<10> energy;   // formerly ap_uint<10>
  //  uint8_t  timing; // formerly ap_uint<4>

  crystal(){
    energy = 0;
    //    timing = 0;
  }

  crystal(ap_uint<10> energy){  // To-do: add timing information
    this->energy = energy;
    //    this->timing = 0; 
  }

  crystal& operator=(const crystal& rhs){
    energy = rhs.energy;
    //    timing = rhs.timing;
    return *this;
  }
};

/*
 * crystalMax class: adapted from algo_top.h
 */ 
class crystalMax{
 public:
  ap_uint<10> energy;
  uint8_t phiMax;
  uint8_t etaMax;

  crystalMax(){
    energy = 0;
    phiMax = 0;
    etaMax = 0;
  }

  crystalMax& operator=(const crystalMax& rhs){
    energy = rhs.energy;
    phiMax = rhs.phiMax;
    etaMax = rhs.etaMax;
    return *this;
  }
};
/*******************************************************************/

// Copied from algo_top.h, with the exception of changing out the ap uints to ap_uint<10>

typedef struct {
  ap_uint<10> energy;
  ap_uint<5> phi;
  ap_uint<5> eta;
} ecaltp_t ;

typedef struct {
  ecaltp_t cr0;
  ecaltp_t cr1 ;
  ecaltp_t cr2 ;
  ecaltp_t cr3 ;
  ecaltp_t cr4 ;
  ecaltp_t cr5 ;
  ecaltp_t cr6 ;
  ecaltp_t cr7 ;
  ecaltp_t cr8 ;
  ecaltp_t cr9 ;
  ecaltp_t cr10 ;
  ecaltp_t cr11 ;
  ecaltp_t cr12 ;
  ecaltp_t cr13 ;
  ecaltp_t cr14 ;
  ecaltp_t cr15 ;
  ecaltp_t cr16 ;
  ecaltp_t cr17 ;
  ecaltp_t cr18 ;
  ecaltp_t cr19 ;
} etaStrip_t ;

typedef struct {
  etaStrip_t etaStrip0 ;
  etaStrip_t etaStrip1 ;
  etaStrip_t etaStrip2 ;
  etaStrip_t etaStrip3 ;
  etaStrip_t etaStrip4 ;
  etaStrip_t etaStrip5 ;
  etaStrip_t etaStrip6 ;
  etaStrip_t etaStrip7 ;
  etaStrip_t etaStrip8 ;
  etaStrip_t etaStrip9 ;
  etaStrip_t etaStrip10 ;
  etaStrip_t etaStrip11 ;
  etaStrip_t etaStrip12 ;
  etaStrip_t etaStrip13 ;
  etaStrip_t etaStrip14 ;
} ecalRegion_t ;


typedef struct {
  ecaltp_t pk0;
  ecaltp_t pk1;
  ecaltp_t pk2;
  ecaltp_t pk3;
  ecaltp_t pk4;
  ecaltp_t pk5;
  ecaltp_t pk6;
  ecaltp_t pk7;
  ecaltp_t pk8;
  ecaltp_t pk9;
  ecaltp_t pk10;
  ecaltp_t pk11;
  ecaltp_t pk12;
  ecaltp_t pk13;
  ecaltp_t pk14;
} etaStripPeak_t ;

class tower_t {
 public:
  ap_uint<16> data;

  tower_t(){
    data = 0;
  }
  tower_t& operator=(const tower_t& rhs){
    data = rhs.data;
    return *this;
  }

  tower_t(ap_uint<12> et, ap_uint<4> hoe) {
    data = (et) | 
      (((ap_uint<16>) hoe)  << 12) ; 
    //      (((ap_uint<16>) satur)  << 15);
  }

  ap_uint<12> et() {return (data & 0xFFF);}
  // ap_uint<3> hoe() {return ((data >> 12) & 0x7);}
  ap_uint<4> hoe() {return ((data >> 12) & 0xF); }
  // ap_uint<1> satur() {return ((data >> 15) & 0x2);}

  float getEt() { return (float) et() / 8.0; }
  
  operator uint16_t() {return (uint16_t) data;}
  
  // Only for ECAL towers! Apply calibration and modify the et() value.
  void applyCalibration(float factor, bool verbose = false) {  
    if (verbose) {
      std::cout << "    applyCalibration: old data as hex " << std::hex << data << std::endl;
      std::cout << std::dec;
    }
    // Get the new pT as a float
    float newEt = getEt() * factor;

    if (verbose) {
      std::cout << "                      uncalibrated Et (unsigned int): " << et() << ", float: " << getEt() 
                << ", calibration factor: " << factor << ", giving new Et (float): " << newEt << ". ";
    }
    // Convert the new pT to an unsigned int (16 bits so we can take the logical OR with the bit mask later)
    ap_uint<16> newEt_uint = (ap_uint<16>) (int) (newEt * 8.0);
    // Make sure the first four bits are zero
    newEt_uint = (newEt_uint & 0x0FFF); 

    if (verbose) {
      std::cout << "      Multiplying by 8.0 to convert to unsigned int: " << (int) newEt_uint << std::endl; 
    }
    // Modify 'data'
    ap_uint<16> bitMask = 0xF000;  // last twelve digits are zero
    data = (data & bitMask);                 // zero out the last twelve digits
    data = (data | newEt_uint);              // write in the new ET
    
    // std::cout << "      New Et as bits: " << et() << ", old data as hex " << std::hex << data << std::endl;
    // std::cout << std::dec;
  }

  // For towers: Calculate H/E ratio given the ECAL and HCAL energies and modify the hoe() value.
  void getHoverE(ap_uint<12> ECAL, ap_uint<12> HCAL_inHcalConvention){
    // Convert HCAL ET to ECAL ET convention 
    ap_uint<12> HCAL = convertHcalETtoEcalET(HCAL_inHcalConvention);

    // std::cout << " getHoverE with ECAL (int) " << (int) ECAL << ", HCAL (int) " << (int) HCAL_inHcalConvention
    //           << " (which, in ECAL ET convention, is): " << (int) HCAL << std::endl;
    // std::cout << " old data as int: " << (int) data << std::endl;      
    ap_uint<4> hoeOut ;
    ap_uint<1> hoeLSB = 0 ;
    ap_uint<4> hoe = 0 ;
    ap_uint<12> A;
    ap_uint<12> B;

    A = (ECAL > HCAL)?ECAL:HCAL;
    B = (ECAL > HCAL)?HCAL:ECAL;
 
    if( ECAL == 0 || HCAL == 0 || HCAL >= ECAL) hoeLSB = 0 ;
    else hoeLSB = 1 ;
    if( A > B ){
      if(A > 2*B) hoe = 0x1 ;
      if(A > 4*B) hoe = 0x2 ;
      if(A > 8*B) hoe = 0x3 ;
      if(A > 16*B) hoe = 0x4 ;
      if(A > 32*B) hoe = 0x5 ;
      if(A > 64*B) hoe = 0x6 ;
      if(A > 128*B) hoe = 0x7 ;
    }
    hoeOut = hoeLSB | (hoe << 1) ;
    
    // std::cout << "  hoeOut as int: " << (int) hoeOut << std::endl;
    ap_uint<16> hoeOutLong = ( ( ((ap_uint<16>) hoeOut) << 12 ) | 0x0000 ) ;  // e.g. 0b ____ 0000 0000 0000 where ___ are the hoe digits
    // std::cout << "   hoeOutLong as int: " << (int) hoeOutLong << std::endl;

    // Take the logical OR to preserve the saturation and tower ET bits
    ap_uint<16> bitMask = 0x0FFF;  // 0000 1111 1111 1111 : zero-out the HoE bits
    data = (data & bitMask);       // zero-out the HoE bits
    data = (data | hoeOutLong);    // write in the new HoE bits
    // std::cout << "    new data as int: " << (int) data << std::endl;
  }
  
};


//--------------------------------------------------------// 

ecalRegion_t initStructure(crystal temporary[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI]){

  ap_uint<5> Eta = 0x0 ;
  ap_uint<5> Phi = 0x0 ;

  ecalRegion_t out;

  out.etaStrip0.cr0.energy=temporary[Eta+0][Phi+0].energy; out.etaStrip0.cr0.eta=0; out.etaStrip0.cr0.phi=0;
  out.etaStrip0.cr1.energy=temporary[Eta+0][Phi+1].energy; out.etaStrip0.cr1.eta=0; out.etaStrip0.cr1.phi=1;
  out.etaStrip0.cr2.energy=temporary[Eta+0][Phi+2].energy; out.etaStrip0.cr2.eta=0; out.etaStrip0.cr2.phi=2;
  out.etaStrip0.cr3.energy=temporary[Eta+0][Phi+3].energy; out.etaStrip0.cr3.eta=0; out.etaStrip0.cr3.phi=3;
  out.etaStrip0.cr4.energy=temporary[Eta+0][Phi+4].energy; out.etaStrip0.cr4.eta=0; out.etaStrip0.cr4.phi=4;
  out.etaStrip0.cr5.energy=temporary[Eta+0][Phi+5].energy; out.etaStrip0.cr5.eta=0; out.etaStrip0.cr5.phi=5;
  out.etaStrip0.cr6.energy=temporary[Eta+0][Phi+6].energy; out.etaStrip0.cr6.eta=0; out.etaStrip0.cr6.phi=6;
  out.etaStrip0.cr7.energy=temporary[Eta+0][Phi+7].energy; out.etaStrip0.cr7.eta=0; out.etaStrip0.cr7.phi=7;
  out.etaStrip0.cr8.energy=temporary[Eta+0][Phi+8].energy; out.etaStrip0.cr8.eta=0; out.etaStrip0.cr8.phi=8;
  out.etaStrip0.cr9.energy=temporary[Eta+0][Phi+9].energy; out.etaStrip0.cr9.eta=0; out.etaStrip0.cr9.phi=9;
  out.etaStrip0.cr10.energy=temporary[Eta+0][Phi+10].energy; out.etaStrip0.cr10.eta=0; out.etaStrip0.cr10.phi=10;
  out.etaStrip0.cr11.energy=temporary[Eta+0][Phi+11].energy; out.etaStrip0.cr11.eta=0; out.etaStrip0.cr11.phi=11;
  out.etaStrip0.cr12.energy=temporary[Eta+0][Phi+12].energy; out.etaStrip0.cr12.eta=0; out.etaStrip0.cr12.phi=12;
  out.etaStrip0.cr13.energy=temporary[Eta+0][Phi+13].energy; out.etaStrip0.cr13.eta=0; out.etaStrip0.cr13.phi=13;
  out.etaStrip0.cr14.energy=temporary[Eta+0][Phi+14].energy; out.etaStrip0.cr14.eta=0; out.etaStrip0.cr14.phi=14;
  out.etaStrip0.cr15.energy=temporary[Eta+0][Phi+15].energy; out.etaStrip0.cr15.eta=0; out.etaStrip0.cr15.phi=15;
  out.etaStrip0.cr16.energy=temporary[Eta+0][Phi+16].energy; out.etaStrip0.cr16.eta=0; out.etaStrip0.cr16.phi=16;
  out.etaStrip0.cr17.energy=temporary[Eta+0][Phi+17].energy; out.etaStrip0.cr17.eta=0; out.etaStrip0.cr17.phi=17;
  out.etaStrip0.cr18.energy=temporary[Eta+0][Phi+18].energy; out.etaStrip0.cr18.eta=0; out.etaStrip0.cr18.phi=18;
  out.etaStrip0.cr19.energy=temporary[Eta+0][Phi+19].energy; out.etaStrip0.cr19.eta=0; out.etaStrip0.cr19.phi=19;

  out.etaStrip1.cr0.energy=temporary[Eta+1][Phi+0].energy; out.etaStrip1.cr0.eta=1; out.etaStrip1.cr0.phi=0;
  out.etaStrip1.cr1.energy=temporary[Eta+1][Phi+1].energy; out.etaStrip1.cr1.eta=1; out.etaStrip1.cr1.phi=1;
  out.etaStrip1.cr2.energy=temporary[Eta+1][Phi+2].energy; out.etaStrip1.cr2.eta=1; out.etaStrip1.cr2.phi=2;
  out.etaStrip1.cr3.energy=temporary[Eta+1][Phi+3].energy; out.etaStrip1.cr3.eta=1; out.etaStrip1.cr3.phi=3;
  out.etaStrip1.cr4.energy=temporary[Eta+1][Phi+4].energy; out.etaStrip1.cr4.eta=1; out.etaStrip1.cr4.phi=4;
  out.etaStrip1.cr5.energy=temporary[Eta+1][Phi+5].energy; out.etaStrip1.cr5.eta=1; out.etaStrip1.cr5.phi=5;
  out.etaStrip1.cr6.energy=temporary[Eta+1][Phi+6].energy; out.etaStrip1.cr6.eta=1; out.etaStrip1.cr6.phi=6;
  out.etaStrip1.cr7.energy=temporary[Eta+1][Phi+7].energy; out.etaStrip1.cr7.eta=1; out.etaStrip1.cr7.phi=7;
  out.etaStrip1.cr8.energy=temporary[Eta+1][Phi+8].energy; out.etaStrip1.cr8.eta=1; out.etaStrip1.cr8.phi=8;
  out.etaStrip1.cr9.energy=temporary[Eta+1][Phi+9].energy; out.etaStrip1.cr9.eta=1; out.etaStrip1.cr9.phi=9;
  out.etaStrip1.cr10.energy=temporary[Eta+1][Phi+10].energy; out.etaStrip1.cr10.eta=1; out.etaStrip1.cr10.phi=10;
  out.etaStrip1.cr11.energy=temporary[Eta+1][Phi+11].energy; out.etaStrip1.cr11.eta=1; out.etaStrip1.cr11.phi=11;
  out.etaStrip1.cr12.energy=temporary[Eta+1][Phi+12].energy; out.etaStrip1.cr12.eta=1; out.etaStrip1.cr12.phi=12;
  out.etaStrip1.cr13.energy=temporary[Eta+1][Phi+13].energy; out.etaStrip1.cr13.eta=1; out.etaStrip1.cr13.phi=13;
  out.etaStrip1.cr14.energy=temporary[Eta+1][Phi+14].energy; out.etaStrip1.cr14.eta=1; out.etaStrip1.cr14.phi=14;
  out.etaStrip1.cr15.energy=temporary[Eta+1][Phi+15].energy; out.etaStrip1.cr15.eta=1; out.etaStrip1.cr15.phi=15;
  out.etaStrip1.cr16.energy=temporary[Eta+1][Phi+16].energy; out.etaStrip1.cr16.eta=1; out.etaStrip1.cr16.phi=16;
  out.etaStrip1.cr17.energy=temporary[Eta+1][Phi+17].energy; out.etaStrip1.cr17.eta=1; out.etaStrip1.cr17.phi=17;
  out.etaStrip1.cr18.energy=temporary[Eta+1][Phi+18].energy; out.etaStrip1.cr18.eta=1; out.etaStrip1.cr18.phi=18;
  out.etaStrip1.cr19.energy=temporary[Eta+1][Phi+19].energy; out.etaStrip1.cr19.eta=1; out.etaStrip1.cr19.phi=19;

  out.etaStrip2.cr0.energy=temporary[Eta+2][Phi+0].energy; out.etaStrip2.cr0.eta=2; out.etaStrip2.cr0.phi=0;
  out.etaStrip2.cr1.energy=temporary[Eta+2][Phi+1].energy; out.etaStrip2.cr1.eta=2; out.etaStrip2.cr1.phi=1;
  out.etaStrip2.cr2.energy=temporary[Eta+2][Phi+2].energy; out.etaStrip2.cr2.eta=2; out.etaStrip2.cr2.phi=2;
  out.etaStrip2.cr3.energy=temporary[Eta+2][Phi+3].energy; out.etaStrip2.cr3.eta=2; out.etaStrip2.cr3.phi=3;
  out.etaStrip2.cr4.energy=temporary[Eta+2][Phi+4].energy; out.etaStrip2.cr4.eta=2; out.etaStrip2.cr4.phi=4;
  out.etaStrip2.cr5.energy=temporary[Eta+2][Phi+5].energy; out.etaStrip2.cr5.eta=2; out.etaStrip2.cr5.phi=5;
  out.etaStrip2.cr6.energy=temporary[Eta+2][Phi+6].energy; out.etaStrip2.cr6.eta=2; out.etaStrip2.cr6.phi=6;
  out.etaStrip2.cr7.energy=temporary[Eta+2][Phi+7].energy; out.etaStrip2.cr7.eta=2; out.etaStrip2.cr7.phi=7;
  out.etaStrip2.cr8.energy=temporary[Eta+2][Phi+8].energy; out.etaStrip2.cr8.eta=2; out.etaStrip2.cr8.phi=8;
  out.etaStrip2.cr9.energy=temporary[Eta+2][Phi+9].energy; out.etaStrip2.cr9.eta=2; out.etaStrip2.cr9.phi=9;
  out.etaStrip2.cr10.energy=temporary[Eta+2][Phi+10].energy; out.etaStrip2.cr10.eta=2; out.etaStrip2.cr10.phi=10;
  out.etaStrip2.cr11.energy=temporary[Eta+2][Phi+11].energy; out.etaStrip2.cr11.eta=2; out.etaStrip2.cr11.phi=11;
  out.etaStrip2.cr12.energy=temporary[Eta+2][Phi+12].energy; out.etaStrip2.cr12.eta=2; out.etaStrip2.cr12.phi=12;
  out.etaStrip2.cr13.energy=temporary[Eta+2][Phi+13].energy; out.etaStrip2.cr13.eta=2; out.etaStrip2.cr13.phi=13;
  out.etaStrip2.cr14.energy=temporary[Eta+2][Phi+14].energy; out.etaStrip2.cr14.eta=2; out.etaStrip2.cr14.phi=14;
  out.etaStrip2.cr15.energy=temporary[Eta+2][Phi+15].energy; out.etaStrip2.cr15.eta=2; out.etaStrip2.cr15.phi=15;
  out.etaStrip2.cr16.energy=temporary[Eta+2][Phi+16].energy; out.etaStrip2.cr16.eta=2; out.etaStrip2.cr16.phi=16;
  out.etaStrip2.cr17.energy=temporary[Eta+2][Phi+17].energy; out.etaStrip2.cr17.eta=2; out.etaStrip2.cr17.phi=17;
  out.etaStrip2.cr18.energy=temporary[Eta+2][Phi+18].energy; out.etaStrip2.cr18.eta=2; out.etaStrip2.cr18.phi=18;
  out.etaStrip2.cr19.energy=temporary[Eta+2][Phi+19].energy; out.etaStrip2.cr19.eta=2; out.etaStrip2.cr19.phi=19;

  out.etaStrip3.cr0.energy=temporary[Eta+3][Phi+0].energy; out.etaStrip3.cr0.eta=3; out.etaStrip3.cr0.phi=0;
  out.etaStrip3.cr1.energy=temporary[Eta+3][Phi+1].energy; out.etaStrip3.cr1.eta=3; out.etaStrip3.cr1.phi=1;
  out.etaStrip3.cr2.energy=temporary[Eta+3][Phi+2].energy; out.etaStrip3.cr2.eta=3; out.etaStrip3.cr2.phi=2;
  out.etaStrip3.cr3.energy=temporary[Eta+3][Phi+3].energy; out.etaStrip3.cr3.eta=3; out.etaStrip3.cr3.phi=3;
  out.etaStrip3.cr4.energy=temporary[Eta+3][Phi+4].energy; out.etaStrip3.cr4.eta=3; out.etaStrip3.cr4.phi=4;
  out.etaStrip3.cr5.energy=temporary[Eta+3][Phi+5].energy; out.etaStrip3.cr5.eta=3; out.etaStrip3.cr5.phi=5;
  out.etaStrip3.cr6.energy=temporary[Eta+3][Phi+6].energy; out.etaStrip3.cr6.eta=3; out.etaStrip3.cr6.phi=6;
  out.etaStrip3.cr7.energy=temporary[Eta+3][Phi+7].energy; out.etaStrip3.cr7.eta=3; out.etaStrip3.cr7.phi=7;
  out.etaStrip3.cr8.energy=temporary[Eta+3][Phi+8].energy; out.etaStrip3.cr8.eta=3; out.etaStrip3.cr8.phi=8;
  out.etaStrip3.cr9.energy=temporary[Eta+3][Phi+9].energy; out.etaStrip3.cr9.eta=3; out.etaStrip3.cr9.phi=9;
  out.etaStrip3.cr10.energy=temporary[Eta+3][Phi+10].energy; out.etaStrip3.cr10.eta=3; out.etaStrip3.cr10.phi=10;
  out.etaStrip3.cr11.energy=temporary[Eta+3][Phi+11].energy; out.etaStrip3.cr11.eta=3; out.etaStrip3.cr11.phi=11;
  out.etaStrip3.cr12.energy=temporary[Eta+3][Phi+12].energy; out.etaStrip3.cr12.eta=3; out.etaStrip3.cr12.phi=12;
  out.etaStrip3.cr13.energy=temporary[Eta+3][Phi+13].energy; out.etaStrip3.cr13.eta=3; out.etaStrip3.cr13.phi=13;
  out.etaStrip3.cr14.energy=temporary[Eta+3][Phi+14].energy; out.etaStrip3.cr14.eta=3; out.etaStrip3.cr14.phi=14;
  out.etaStrip3.cr15.energy=temporary[Eta+3][Phi+15].energy; out.etaStrip3.cr15.eta=3; out.etaStrip3.cr15.phi=15;
  out.etaStrip3.cr16.energy=temporary[Eta+3][Phi+16].energy; out.etaStrip3.cr16.eta=3; out.etaStrip3.cr16.phi=16;
  out.etaStrip3.cr17.energy=temporary[Eta+3][Phi+17].energy; out.etaStrip3.cr17.eta=3; out.etaStrip3.cr17.phi=17;
  out.etaStrip3.cr18.energy=temporary[Eta+3][Phi+18].energy; out.etaStrip3.cr18.eta=3; out.etaStrip3.cr18.phi=18;
  out.etaStrip3.cr19.energy=temporary[Eta+3][Phi+19].energy; out.etaStrip3.cr19.eta=3; out.etaStrip3.cr19.phi=19;

  out.etaStrip4.cr0.energy=temporary[Eta+4][Phi+0].energy; out.etaStrip4.cr0.eta=4; out.etaStrip4.cr0.phi=0;
  out.etaStrip4.cr1.energy=temporary[Eta+4][Phi+1].energy; out.etaStrip4.cr1.eta=4; out.etaStrip4.cr1.phi=1;
  out.etaStrip4.cr2.energy=temporary[Eta+4][Phi+2].energy; out.etaStrip4.cr2.eta=4; out.etaStrip4.cr2.phi=2;
  out.etaStrip4.cr3.energy=temporary[Eta+4][Phi+3].energy; out.etaStrip4.cr3.eta=4; out.etaStrip4.cr3.phi=3;
  out.etaStrip4.cr4.energy=temporary[Eta+4][Phi+4].energy; out.etaStrip4.cr4.eta=4; out.etaStrip4.cr4.phi=4;
  out.etaStrip4.cr5.energy=temporary[Eta+4][Phi+5].energy; out.etaStrip4.cr5.eta=4; out.etaStrip4.cr5.phi=5;
  out.etaStrip4.cr6.energy=temporary[Eta+4][Phi+6].energy; out.etaStrip4.cr6.eta=4; out.etaStrip4.cr6.phi=6;
  out.etaStrip4.cr7.energy=temporary[Eta+4][Phi+7].energy; out.etaStrip4.cr7.eta=4; out.etaStrip4.cr7.phi=7;
  out.etaStrip4.cr8.energy=temporary[Eta+4][Phi+8].energy; out.etaStrip4.cr8.eta=4; out.etaStrip4.cr8.phi=8;
  out.etaStrip4.cr9.energy=temporary[Eta+4][Phi+9].energy; out.etaStrip4.cr9.eta=4; out.etaStrip4.cr9.phi=9;
  out.etaStrip4.cr10.energy=temporary[Eta+4][Phi+10].energy; out.etaStrip4.cr10.eta=4; out.etaStrip4.cr10.phi=10;
  out.etaStrip4.cr11.energy=temporary[Eta+4][Phi+11].energy; out.etaStrip4.cr11.eta=4; out.etaStrip4.cr11.phi=11;
  out.etaStrip4.cr12.energy=temporary[Eta+4][Phi+12].energy; out.etaStrip4.cr12.eta=4; out.etaStrip4.cr12.phi=12;
  out.etaStrip4.cr13.energy=temporary[Eta+4][Phi+13].energy; out.etaStrip4.cr13.eta=4; out.etaStrip4.cr13.phi=13;
  out.etaStrip4.cr14.energy=temporary[Eta+4][Phi+14].energy; out.etaStrip4.cr14.eta=4; out.etaStrip4.cr14.phi=14;
  out.etaStrip4.cr15.energy=temporary[Eta+4][Phi+15].energy; out.etaStrip4.cr15.eta=4; out.etaStrip4.cr15.phi=15;
  out.etaStrip4.cr16.energy=temporary[Eta+4][Phi+16].energy; out.etaStrip4.cr16.eta=4; out.etaStrip4.cr16.phi=16;
  out.etaStrip4.cr17.energy=temporary[Eta+4][Phi+17].energy; out.etaStrip4.cr17.eta=4; out.etaStrip4.cr17.phi=17;
  out.etaStrip4.cr18.energy=temporary[Eta+4][Phi+18].energy; out.etaStrip4.cr18.eta=4; out.etaStrip4.cr18.phi=18;
  out.etaStrip4.cr19.energy=temporary[Eta+4][Phi+19].energy; out.etaStrip4.cr19.eta=4; out.etaStrip4.cr19.phi=19;

  out.etaStrip5.cr0.energy=temporary[Eta+5][Phi+0].energy; out.etaStrip5.cr0.eta=5; out.etaStrip5.cr0.phi=0;
  out.etaStrip5.cr1.energy=temporary[Eta+5][Phi+1].energy; out.etaStrip5.cr1.eta=5; out.etaStrip5.cr1.phi=1;
  out.etaStrip5.cr2.energy=temporary[Eta+5][Phi+2].energy; out.etaStrip5.cr2.eta=5; out.etaStrip5.cr2.phi=2;
  out.etaStrip5.cr3.energy=temporary[Eta+5][Phi+3].energy; out.etaStrip5.cr3.eta=5; out.etaStrip5.cr3.phi=3;
  out.etaStrip5.cr4.energy=temporary[Eta+5][Phi+4].energy; out.etaStrip5.cr4.eta=5; out.etaStrip5.cr4.phi=4;
  out.etaStrip5.cr5.energy=temporary[Eta+5][Phi+5].energy; out.etaStrip5.cr5.eta=5; out.etaStrip5.cr5.phi=5;
  out.etaStrip5.cr6.energy=temporary[Eta+5][Phi+6].energy; out.etaStrip5.cr6.eta=5; out.etaStrip5.cr6.phi=6;
  out.etaStrip5.cr7.energy=temporary[Eta+5][Phi+7].energy; out.etaStrip5.cr7.eta=5; out.etaStrip5.cr7.phi=7;
  out.etaStrip5.cr8.energy=temporary[Eta+5][Phi+8].energy; out.etaStrip5.cr8.eta=5; out.etaStrip5.cr8.phi=8;
  out.etaStrip5.cr9.energy=temporary[Eta+5][Phi+9].energy; out.etaStrip5.cr9.eta=5; out.etaStrip5.cr9.phi=9;
  out.etaStrip5.cr10.energy=temporary[Eta+5][Phi+10].energy; out.etaStrip5.cr10.eta=5; out.etaStrip5.cr10.phi=10;
  out.etaStrip5.cr11.energy=temporary[Eta+5][Phi+11].energy; out.etaStrip5.cr11.eta=5; out.etaStrip5.cr11.phi=11;
  out.etaStrip5.cr12.energy=temporary[Eta+5][Phi+12].energy; out.etaStrip5.cr12.eta=5; out.etaStrip5.cr12.phi=12;
  out.etaStrip5.cr13.energy=temporary[Eta+5][Phi+13].energy; out.etaStrip5.cr13.eta=5; out.etaStrip5.cr13.phi=13;
  out.etaStrip5.cr14.energy=temporary[Eta+5][Phi+14].energy; out.etaStrip5.cr14.eta=5; out.etaStrip5.cr14.phi=14;
  out.etaStrip5.cr15.energy=temporary[Eta+5][Phi+15].energy; out.etaStrip5.cr15.eta=5; out.etaStrip5.cr15.phi=15;
  out.etaStrip5.cr16.energy=temporary[Eta+5][Phi+16].energy; out.etaStrip5.cr16.eta=5; out.etaStrip5.cr16.phi=16;
  out.etaStrip5.cr17.energy=temporary[Eta+5][Phi+17].energy; out.etaStrip5.cr17.eta=5; out.etaStrip5.cr17.phi=17;
  out.etaStrip5.cr18.energy=temporary[Eta+5][Phi+18].energy; out.etaStrip5.cr18.eta=5; out.etaStrip5.cr18.phi=18;
  out.etaStrip5.cr19.energy=temporary[Eta+5][Phi+19].energy; out.etaStrip5.cr19.eta=5; out.etaStrip5.cr19.phi=19;

  out.etaStrip6.cr0.energy=temporary[Eta+6][Phi+0].energy; out.etaStrip6.cr0.eta=6; out.etaStrip6.cr0.phi=0;
  out.etaStrip6.cr1.energy=temporary[Eta+6][Phi+1].energy; out.etaStrip6.cr1.eta=6; out.etaStrip6.cr1.phi=1;
  out.etaStrip6.cr2.energy=temporary[Eta+6][Phi+2].energy; out.etaStrip6.cr2.eta=6; out.etaStrip6.cr2.phi=2;
  out.etaStrip6.cr2.energy=temporary[Eta+6][Phi+2].energy; out.etaStrip6.cr2.eta=6; out.etaStrip6.cr2.phi=2;
  out.etaStrip6.cr3.energy=temporary[Eta+6][Phi+3].energy; out.etaStrip6.cr3.eta=6; out.etaStrip6.cr3.phi=3;
  out.etaStrip6.cr4.energy=temporary[Eta+6][Phi+4].energy; out.etaStrip6.cr4.eta=6; out.etaStrip6.cr4.phi=4;
  out.etaStrip6.cr5.energy=temporary[Eta+6][Phi+5].energy; out.etaStrip6.cr5.eta=6; out.etaStrip6.cr5.phi=5;
  out.etaStrip6.cr6.energy=temporary[Eta+6][Phi+6].energy; out.etaStrip6.cr6.eta=6; out.etaStrip6.cr6.phi=6;
  out.etaStrip6.cr7.energy=temporary[Eta+6][Phi+7].energy; out.etaStrip6.cr7.eta=6; out.etaStrip6.cr7.phi=7;
  out.etaStrip6.cr8.energy=temporary[Eta+6][Phi+8].energy; out.etaStrip6.cr8.eta=6; out.etaStrip6.cr8.phi=8;
  out.etaStrip6.cr9.energy=temporary[Eta+6][Phi+9].energy; out.etaStrip6.cr9.eta=6; out.etaStrip6.cr9.phi=9;
  out.etaStrip6.cr10.energy=temporary[Eta+6][Phi+10].energy; out.etaStrip6.cr10.eta=6; out.etaStrip6.cr10.phi=10;
  out.etaStrip6.cr11.energy=temporary[Eta+6][Phi+11].energy; out.etaStrip6.cr11.eta=6; out.etaStrip6.cr11.phi=11;
  out.etaStrip6.cr12.energy=temporary[Eta+6][Phi+12].energy; out.etaStrip6.cr12.eta=6; out.etaStrip6.cr12.phi=12;
  out.etaStrip6.cr13.energy=temporary[Eta+6][Phi+13].energy; out.etaStrip6.cr13.eta=6; out.etaStrip6.cr13.phi=13;
  out.etaStrip6.cr14.energy=temporary[Eta+6][Phi+14].energy; out.etaStrip6.cr14.eta=6; out.etaStrip6.cr14.phi=14;
  out.etaStrip6.cr15.energy=temporary[Eta+6][Phi+15].energy; out.etaStrip6.cr15.eta=6; out.etaStrip6.cr15.phi=15;
  out.etaStrip6.cr16.energy=temporary[Eta+6][Phi+16].energy; out.etaStrip6.cr16.eta=6; out.etaStrip6.cr16.phi=16;
  out.etaStrip6.cr17.energy=temporary[Eta+6][Phi+17].energy; out.etaStrip6.cr17.eta=6; out.etaStrip6.cr17.phi=17;
  out.etaStrip6.cr18.energy=temporary[Eta+6][Phi+18].energy; out.etaStrip6.cr18.eta=6; out.etaStrip6.cr18.phi=18;
  out.etaStrip6.cr19.energy=temporary[Eta+6][Phi+19].energy; out.etaStrip6.cr19.eta=6; out.etaStrip6.cr19.phi=19;

  out.etaStrip7.cr0.energy=temporary[Eta+7][Phi+0].energy; out.etaStrip7.cr0.eta=7; out.etaStrip7.cr0.phi=0;
  out.etaStrip7.cr1.energy=temporary[Eta+7][Phi+1].energy; out.etaStrip7.cr1.eta=7; out.etaStrip7.cr1.phi=1;
  out.etaStrip7.cr2.energy=temporary[Eta+7][Phi+2].energy; out.etaStrip7.cr2.eta=7; out.etaStrip7.cr2.phi=2;
  out.etaStrip7.cr3.energy=temporary[Eta+7][Phi+3].energy; out.etaStrip7.cr3.eta=7; out.etaStrip7.cr3.phi=3;
  out.etaStrip7.cr4.energy=temporary[Eta+7][Phi+4].energy; out.etaStrip7.cr4.eta=7; out.etaStrip7.cr4.phi=4;
  out.etaStrip7.cr5.energy=temporary[Eta+7][Phi+5].energy; out.etaStrip7.cr5.eta=7; out.etaStrip7.cr5.phi=5;
  out.etaStrip7.cr6.energy=temporary[Eta+7][Phi+6].energy; out.etaStrip7.cr6.eta=7; out.etaStrip7.cr6.phi=6;
  out.etaStrip7.cr7.energy=temporary[Eta+7][Phi+7].energy; out.etaStrip7.cr7.eta=7; out.etaStrip7.cr7.phi=7;
  out.etaStrip7.cr8.energy=temporary[Eta+7][Phi+8].energy; out.etaStrip7.cr8.eta=7; out.etaStrip7.cr8.phi=8;
  out.etaStrip7.cr9.energy=temporary[Eta+7][Phi+9].energy; out.etaStrip7.cr9.eta=7; out.etaStrip7.cr9.phi=9;
  out.etaStrip7.cr10.energy=temporary[Eta+7][Phi+10].energy; out.etaStrip7.cr10.eta=7; out.etaStrip7.cr10.phi=10;
  out.etaStrip7.cr11.energy=temporary[Eta+7][Phi+11].energy; out.etaStrip7.cr11.eta=7; out.etaStrip7.cr11.phi=11;
  out.etaStrip7.cr12.energy=temporary[Eta+7][Phi+12].energy; out.etaStrip7.cr12.eta=7; out.etaStrip7.cr12.phi=12;
  out.etaStrip7.cr13.energy=temporary[Eta+7][Phi+13].energy; out.etaStrip7.cr13.eta=7; out.etaStrip7.cr13.phi=13;
  out.etaStrip7.cr14.energy=temporary[Eta+7][Phi+14].energy; out.etaStrip7.cr14.eta=7; out.etaStrip7.cr14.phi=14;
  out.etaStrip7.cr15.energy=temporary[Eta+7][Phi+15].energy; out.etaStrip7.cr15.eta=7; out.etaStrip7.cr15.phi=15;
  out.etaStrip7.cr16.energy=temporary[Eta+7][Phi+16].energy; out.etaStrip7.cr16.eta=7; out.etaStrip7.cr16.phi=16;
  out.etaStrip7.cr17.energy=temporary[Eta+7][Phi+17].energy; out.etaStrip7.cr17.eta=7; out.etaStrip7.cr17.phi=17;
  out.etaStrip7.cr18.energy=temporary[Eta+7][Phi+18].energy; out.etaStrip7.cr18.eta=7; out.etaStrip7.cr18.phi=18;
  out.etaStrip7.cr19.energy=temporary[Eta+7][Phi+19].energy; out.etaStrip7.cr19.eta=7; out.etaStrip7.cr19.phi=19;

  out.etaStrip8.cr0.energy=temporary[Eta+8][Phi+0].energy; out.etaStrip8.cr0.eta=8; out.etaStrip8.cr0.phi=0;
  out.etaStrip8.cr1.energy=temporary[Eta+8][Phi+1].energy; out.etaStrip8.cr1.eta=8; out.etaStrip8.cr1.phi=1;
  out.etaStrip8.cr2.energy=temporary[Eta+8][Phi+2].energy; out.etaStrip8.cr2.eta=8; out.etaStrip8.cr2.phi=2;
  out.etaStrip8.cr3.energy=temporary[Eta+8][Phi+3].energy; out.etaStrip8.cr3.eta=8; out.etaStrip8.cr3.phi=3;
  out.etaStrip8.cr4.energy=temporary[Eta+8][Phi+4].energy; out.etaStrip8.cr4.eta=8; out.etaStrip8.cr4.phi=4;
  out.etaStrip8.cr5.energy=temporary[Eta+8][Phi+5].energy; out.etaStrip8.cr5.eta=8; out.etaStrip8.cr5.phi=5;
  out.etaStrip8.cr6.energy=temporary[Eta+8][Phi+6].energy; out.etaStrip8.cr6.eta=8; out.etaStrip8.cr6.phi=6;
  out.etaStrip8.cr7.energy=temporary[Eta+8][Phi+7].energy; out.etaStrip8.cr7.eta=8; out.etaStrip8.cr7.phi=7;
  out.etaStrip8.cr8.energy=temporary[Eta+8][Phi+8].energy; out.etaStrip8.cr8.eta=8; out.etaStrip8.cr8.phi=8;
  out.etaStrip8.cr9.energy=temporary[Eta+8][Phi+9].energy; out.etaStrip8.cr9.eta=8; out.etaStrip8.cr9.phi=9;
  out.etaStrip8.cr10.energy=temporary[Eta+8][Phi+10].energy; out.etaStrip8.cr10.eta=8; out.etaStrip8.cr10.phi=10;
  out.etaStrip8.cr11.energy=temporary[Eta+8][Phi+11].energy; out.etaStrip8.cr11.eta=8; out.etaStrip8.cr11.phi=11;
  out.etaStrip8.cr12.energy=temporary[Eta+8][Phi+12].energy; out.etaStrip8.cr12.eta=8; out.etaStrip8.cr12.phi=12;
  out.etaStrip8.cr13.energy=temporary[Eta+8][Phi+13].energy; out.etaStrip8.cr13.eta=8; out.etaStrip8.cr13.phi=13;
  out.etaStrip8.cr14.energy=temporary[Eta+8][Phi+14].energy; out.etaStrip8.cr14.eta=8; out.etaStrip8.cr14.phi=14;
  out.etaStrip8.cr15.energy=temporary[Eta+8][Phi+15].energy; out.etaStrip8.cr15.eta=8; out.etaStrip8.cr15.phi=15;
  out.etaStrip8.cr16.energy=temporary[Eta+8][Phi+16].energy; out.etaStrip8.cr16.eta=8; out.etaStrip8.cr16.phi=16;
  out.etaStrip8.cr17.energy=temporary[Eta+8][Phi+17].energy; out.etaStrip8.cr17.eta=8; out.etaStrip8.cr17.phi=17;
  out.etaStrip8.cr18.energy=temporary[Eta+8][Phi+18].energy; out.etaStrip8.cr18.eta=8; out.etaStrip8.cr18.phi=18;
  out.etaStrip8.cr19.energy=temporary[Eta+8][Phi+19].energy; out.etaStrip8.cr19.eta=8; out.etaStrip8.cr19.phi=19;

  out.etaStrip9.cr0.energy=temporary[Eta+9][Phi+0].energy; out.etaStrip9.cr0.eta=9; out.etaStrip9.cr0.phi=0;
  out.etaStrip9.cr1.energy=temporary[Eta+9][Phi+1].energy; out.etaStrip9.cr1.eta=9; out.etaStrip9.cr1.phi=1;
  out.etaStrip9.cr2.energy=temporary[Eta+9][Phi+2].energy; out.etaStrip9.cr2.eta=9; out.etaStrip9.cr2.phi=2;
  out.etaStrip9.cr3.energy=temporary[Eta+9][Phi+3].energy; out.etaStrip9.cr3.eta=9; out.etaStrip9.cr3.phi=3;
  out.etaStrip9.cr4.energy=temporary[Eta+9][Phi+4].energy; out.etaStrip9.cr4.eta=9; out.etaStrip9.cr4.phi=4;
  out.etaStrip9.cr5.energy=temporary[Eta+9][Phi+5].energy; out.etaStrip9.cr5.eta=9; out.etaStrip9.cr5.phi=5;
  out.etaStrip9.cr6.energy=temporary[Eta+9][Phi+6].energy; out.etaStrip9.cr6.eta=9; out.etaStrip9.cr6.phi=6;
  out.etaStrip9.cr7.energy=temporary[Eta+9][Phi+7].energy; out.etaStrip9.cr7.eta=9; out.etaStrip9.cr7.phi=7;
  out.etaStrip9.cr8.energy=temporary[Eta+9][Phi+8].energy; out.etaStrip9.cr8.eta=9; out.etaStrip9.cr8.phi=8;
  out.etaStrip9.cr9.energy=temporary[Eta+9][Phi+9].energy; out.etaStrip9.cr9.eta=9; out.etaStrip9.cr9.phi=9;
  out.etaStrip9.cr10.energy=temporary[Eta+9][Phi+10].energy; out.etaStrip9.cr10.eta=9; out.etaStrip9.cr10.phi=10;
  out.etaStrip9.cr11.energy=temporary[Eta+9][Phi+11].energy; out.etaStrip9.cr11.eta=9; out.etaStrip9.cr11.phi=11;
  out.etaStrip9.cr12.energy=temporary[Eta+9][Phi+12].energy; out.etaStrip9.cr12.eta=9; out.etaStrip9.cr12.phi=12;
  out.etaStrip9.cr13.energy=temporary[Eta+9][Phi+13].energy; out.etaStrip9.cr13.eta=9; out.etaStrip9.cr13.phi=13;
  out.etaStrip9.cr14.energy=temporary[Eta+9][Phi+14].energy; out.etaStrip9.cr14.eta=9; out.etaStrip9.cr14.phi=14;
  out.etaStrip9.cr15.energy=temporary[Eta+9][Phi+15].energy; out.etaStrip9.cr15.eta=9; out.etaStrip9.cr15.phi=15;
  out.etaStrip9.cr16.energy=temporary[Eta+9][Phi+16].energy; out.etaStrip9.cr16.eta=9; out.etaStrip9.cr16.phi=16;
  out.etaStrip9.cr17.energy=temporary[Eta+9][Phi+17].energy; out.etaStrip9.cr17.eta=9; out.etaStrip9.cr17.phi=17;
  out.etaStrip9.cr18.energy=temporary[Eta+9][Phi+18].energy; out.etaStrip9.cr18.eta=9; out.etaStrip9.cr18.phi=18;
  out.etaStrip9.cr19.energy=temporary[Eta+9][Phi+19].energy; out.etaStrip9.cr19.eta=9; out.etaStrip9.cr19.phi=19;

  out.etaStrip10.cr0.energy=temporary[Eta+10][Phi+0].energy; out.etaStrip10.cr0.eta=10; out.etaStrip10.cr0.phi=0;
  out.etaStrip10.cr1.energy=temporary[Eta+10][Phi+1].energy; out.etaStrip10.cr1.eta=10; out.etaStrip10.cr1.phi=1;
  out.etaStrip10.cr2.energy=temporary[Eta+10][Phi+2].energy; out.etaStrip10.cr2.eta=10; out.etaStrip10.cr2.phi=2;
  out.etaStrip10.cr3.energy=temporary[Eta+10][Phi+3].energy; out.etaStrip10.cr3.eta=10; out.etaStrip10.cr3.phi=3;
  out.etaStrip10.cr4.energy=temporary[Eta+10][Phi+4].energy; out.etaStrip10.cr4.eta=10; out.etaStrip10.cr4.phi=4;
  out.etaStrip10.cr5.energy=temporary[Eta+10][Phi+5].energy; out.etaStrip10.cr5.eta=10; out.etaStrip10.cr5.phi=5;
  out.etaStrip10.cr6.energy=temporary[Eta+10][Phi+6].energy; out.etaStrip10.cr6.eta=10; out.etaStrip10.cr6.phi=6;
  out.etaStrip10.cr7.energy=temporary[Eta+10][Phi+7].energy; out.etaStrip10.cr7.eta=10; out.etaStrip10.cr7.phi=7;
  out.etaStrip10.cr8.energy=temporary[Eta+10][Phi+8].energy; out.etaStrip10.cr8.eta=10; out.etaStrip10.cr8.phi=8;
  out.etaStrip10.cr9.energy=temporary[Eta+10][Phi+9].energy; out.etaStrip10.cr9.eta=10; out.etaStrip10.cr9.phi=9;
  out.etaStrip10.cr10.energy=temporary[Eta+10][Phi+10].energy; out.etaStrip10.cr10.eta=10; out.etaStrip10.cr10.phi=10;
  out.etaStrip10.cr11.energy=temporary[Eta+10][Phi+11].energy; out.etaStrip10.cr11.eta=10; out.etaStrip10.cr11.phi=11;
  out.etaStrip10.cr12.energy=temporary[Eta+10][Phi+12].energy; out.etaStrip10.cr12.eta=10; out.etaStrip10.cr12.phi=12;
  out.etaStrip10.cr13.energy=temporary[Eta+10][Phi+13].energy; out.etaStrip10.cr13.eta=10; out.etaStrip10.cr13.phi=13;
  out.etaStrip10.cr14.energy=temporary[Eta+10][Phi+14].energy; out.etaStrip10.cr14.eta=10; out.etaStrip10.cr14.phi=14;
  out.etaStrip10.cr15.energy=temporary[Eta+10][Phi+15].energy; out.etaStrip10.cr15.eta=10; out.etaStrip10.cr15.phi=15;
  out.etaStrip10.cr16.energy=temporary[Eta+10][Phi+16].energy; out.etaStrip10.cr16.eta=10; out.etaStrip10.cr16.phi=16;
  out.etaStrip10.cr17.energy=temporary[Eta+10][Phi+17].energy; out.etaStrip10.cr17.eta=10; out.etaStrip10.cr17.phi=17;
  out.etaStrip10.cr18.energy=temporary[Eta+10][Phi+18].energy; out.etaStrip10.cr18.eta=10; out.etaStrip10.cr18.phi=18;
  out.etaStrip10.cr19.energy=temporary[Eta+10][Phi+19].energy; out.etaStrip10.cr19.eta=10; out.etaStrip10.cr19.phi=19;

  out.etaStrip11.cr0.energy=temporary[Eta+11][Phi+0].energy; out.etaStrip11.cr0.eta=11; out.etaStrip11.cr0.phi=0;
  out.etaStrip11.cr1.energy=temporary[Eta+11][Phi+1].energy; out.etaStrip11.cr1.eta=11; out.etaStrip11.cr1.phi=1;
  out.etaStrip11.cr2.energy=temporary[Eta+11][Phi+2].energy; out.etaStrip11.cr2.eta=11; out.etaStrip11.cr2.phi=2;
  out.etaStrip11.cr3.energy=temporary[Eta+11][Phi+3].energy; out.etaStrip11.cr3.eta=11; out.etaStrip11.cr3.phi=3;
  out.etaStrip11.cr4.energy=temporary[Eta+11][Phi+4].energy; out.etaStrip11.cr4.eta=11; out.etaStrip11.cr4.phi=4;
  out.etaStrip11.cr5.energy=temporary[Eta+11][Phi+5].energy; out.etaStrip11.cr5.eta=11; out.etaStrip11.cr5.phi=5;
  out.etaStrip11.cr6.energy=temporary[Eta+11][Phi+6].energy; out.etaStrip11.cr6.eta=11; out.etaStrip11.cr6.phi=6;
  out.etaStrip11.cr7.energy=temporary[Eta+11][Phi+7].energy; out.etaStrip11.cr7.eta=11; out.etaStrip11.cr7.phi=7;
  out.etaStrip11.cr8.energy=temporary[Eta+11][Phi+8].energy; out.etaStrip11.cr8.eta=11; out.etaStrip11.cr8.phi=8;
  out.etaStrip11.cr9.energy=temporary[Eta+11][Phi+9].energy; out.etaStrip11.cr9.eta=11; out.etaStrip11.cr9.phi=9;
  out.etaStrip11.cr10.energy=temporary[Eta+11][Phi+10].energy; out.etaStrip11.cr10.eta=11; out.etaStrip11.cr10.phi=10;
  out.etaStrip11.cr11.energy=temporary[Eta+11][Phi+11].energy; out.etaStrip11.cr11.eta=11; out.etaStrip11.cr11.phi=11;
  out.etaStrip11.cr12.energy=temporary[Eta+11][Phi+12].energy; out.etaStrip11.cr12.eta=11; out.etaStrip11.cr12.phi=12;
  out.etaStrip11.cr13.energy=temporary[Eta+11][Phi+13].energy; out.etaStrip11.cr13.eta=11; out.etaStrip11.cr13.phi=13;
  out.etaStrip11.cr14.energy=temporary[Eta+11][Phi+14].energy; out.etaStrip11.cr14.eta=11; out.etaStrip11.cr14.phi=14;
  out.etaStrip11.cr15.energy=temporary[Eta+11][Phi+15].energy; out.etaStrip11.cr15.eta=11; out.etaStrip11.cr15.phi=15;
  out.etaStrip11.cr16.energy=temporary[Eta+11][Phi+16].energy; out.etaStrip11.cr16.eta=11; out.etaStrip11.cr16.phi=16;
  out.etaStrip11.cr17.energy=temporary[Eta+11][Phi+17].energy; out.etaStrip11.cr17.eta=11; out.etaStrip11.cr17.phi=17;
  out.etaStrip11.cr18.energy=temporary[Eta+11][Phi+18].energy; out.etaStrip11.cr18.eta=11; out.etaStrip11.cr18.phi=18;
  out.etaStrip11.cr19.energy=temporary[Eta+11][Phi+19].energy; out.etaStrip11.cr19.eta=11; out.etaStrip11.cr19.phi=19;

  out.etaStrip12.cr0.energy=temporary[Eta+12][Phi+0].energy; out.etaStrip12.cr0.eta=12; out.etaStrip12.cr0.phi=0;
  out.etaStrip12.cr1.energy=temporary[Eta+12][Phi+1].energy; out.etaStrip12.cr1.eta=12; out.etaStrip12.cr1.phi=1;
  out.etaStrip12.cr2.energy=temporary[Eta+12][Phi+2].energy; out.etaStrip12.cr2.eta=12; out.etaStrip12.cr2.phi=2;
  out.etaStrip12.cr3.energy=temporary[Eta+12][Phi+3].energy; out.etaStrip12.cr3.eta=12; out.etaStrip12.cr3.phi=3;
  out.etaStrip12.cr4.energy=temporary[Eta+12][Phi+4].energy; out.etaStrip12.cr4.eta=12; out.etaStrip12.cr4.phi=4;
  out.etaStrip12.cr5.energy=temporary[Eta+12][Phi+5].energy; out.etaStrip12.cr5.eta=12; out.etaStrip12.cr5.phi=5;
  out.etaStrip12.cr6.energy=temporary[Eta+12][Phi+6].energy; out.etaStrip12.cr6.eta=12; out.etaStrip12.cr6.phi=6;
  out.etaStrip12.cr7.energy=temporary[Eta+12][Phi+7].energy; out.etaStrip12.cr7.eta=12; out.etaStrip12.cr7.phi=7;
  out.etaStrip12.cr8.energy=temporary[Eta+12][Phi+8].energy; out.etaStrip12.cr8.eta=12; out.etaStrip12.cr8.phi=8;
  out.etaStrip12.cr9.energy=temporary[Eta+12][Phi+9].energy; out.etaStrip12.cr9.eta=12; out.etaStrip12.cr9.phi=9;
  out.etaStrip12.cr10.energy=temporary[Eta+12][Phi+10].energy; out.etaStrip12.cr10.eta=12; out.etaStrip12.cr10.phi=10;
  out.etaStrip12.cr11.energy=temporary[Eta+12][Phi+11].energy; out.etaStrip12.cr11.eta=12; out.etaStrip12.cr11.phi=11;
  out.etaStrip12.cr12.energy=temporary[Eta+12][Phi+12].energy; out.etaStrip12.cr12.eta=12; out.etaStrip12.cr12.phi=12;
  out.etaStrip12.cr13.energy=temporary[Eta+12][Phi+13].energy; out.etaStrip12.cr13.eta=12; out.etaStrip12.cr13.phi=13;
  out.etaStrip12.cr14.energy=temporary[Eta+12][Phi+14].energy; out.etaStrip12.cr14.eta=12; out.etaStrip12.cr14.phi=14;
  out.etaStrip12.cr15.energy=temporary[Eta+12][Phi+15].energy; out.etaStrip12.cr15.eta=12; out.etaStrip12.cr15.phi=15;
  out.etaStrip12.cr16.energy=temporary[Eta+12][Phi+16].energy; out.etaStrip12.cr16.eta=12; out.etaStrip12.cr16.phi=16;
  out.etaStrip12.cr17.energy=temporary[Eta+12][Phi+17].energy; out.etaStrip12.cr17.eta=12; out.etaStrip12.cr17.phi=17;
  out.etaStrip12.cr18.energy=temporary[Eta+12][Phi+18].energy; out.etaStrip12.cr18.eta=12; out.etaStrip12.cr18.phi=18;
  out.etaStrip12.cr19.energy=temporary[Eta+12][Phi+19].energy; out.etaStrip12.cr19.eta=12; out.etaStrip12.cr19.phi=19;

  out.etaStrip13.cr0.energy=temporary[Eta+13][Phi+0].energy; out.etaStrip13.cr0.eta=13; out.etaStrip13.cr0.phi=0;
  out.etaStrip13.cr1.energy=temporary[Eta+13][Phi+1].energy; out.etaStrip13.cr1.eta=13; out.etaStrip13.cr1.phi=1;
  out.etaStrip13.cr2.energy=temporary[Eta+13][Phi+2].energy; out.etaStrip13.cr2.eta=13; out.etaStrip13.cr2.phi=2;
  out.etaStrip13.cr3.energy=temporary[Eta+13][Phi+3].energy; out.etaStrip13.cr3.eta=13; out.etaStrip13.cr3.phi=3;
  out.etaStrip13.cr4.energy=temporary[Eta+13][Phi+4].energy; out.etaStrip13.cr4.eta=13; out.etaStrip13.cr4.phi=4;
  out.etaStrip13.cr5.energy=temporary[Eta+13][Phi+5].energy; out.etaStrip13.cr5.eta=13; out.etaStrip13.cr5.phi=5;
  out.etaStrip13.cr6.energy=temporary[Eta+13][Phi+6].energy; out.etaStrip13.cr6.eta=13; out.etaStrip13.cr6.phi=6;
  out.etaStrip13.cr7.energy=temporary[Eta+13][Phi+7].energy; out.etaStrip13.cr7.eta=13; out.etaStrip13.cr7.phi=7;
  out.etaStrip13.cr8.energy=temporary[Eta+13][Phi+8].energy; out.etaStrip13.cr8.eta=13; out.etaStrip13.cr8.phi=8;
  out.etaStrip13.cr9.energy=temporary[Eta+13][Phi+9].energy; out.etaStrip13.cr9.eta=13; out.etaStrip13.cr9.phi=9;
  out.etaStrip13.cr10.energy=temporary[Eta+13][Phi+10].energy; out.etaStrip13.cr10.eta=13; out.etaStrip13.cr10.phi=10;
  out.etaStrip13.cr11.energy=temporary[Eta+13][Phi+11].energy; out.etaStrip13.cr11.eta=13; out.etaStrip13.cr11.phi=11;
  out.etaStrip13.cr12.energy=temporary[Eta+13][Phi+12].energy; out.etaStrip13.cr12.eta=13; out.etaStrip13.cr12.phi=12;
  out.etaStrip13.cr13.energy=temporary[Eta+13][Phi+13].energy; out.etaStrip13.cr13.eta=13; out.etaStrip13.cr13.phi=13;
  out.etaStrip13.cr14.energy=temporary[Eta+13][Phi+14].energy; out.etaStrip13.cr14.eta=13; out.etaStrip13.cr14.phi=14;
  out.etaStrip13.cr15.energy=temporary[Eta+13][Phi+15].energy; out.etaStrip13.cr15.eta=13; out.etaStrip13.cr15.phi=15;
  out.etaStrip13.cr16.energy=temporary[Eta+13][Phi+16].energy; out.etaStrip13.cr16.eta=13; out.etaStrip13.cr16.phi=16;
  out.etaStrip13.cr17.energy=temporary[Eta+13][Phi+17].energy; out.etaStrip13.cr17.eta=13; out.etaStrip13.cr17.phi=17;
  out.etaStrip13.cr18.energy=temporary[Eta+13][Phi+18].energy; out.etaStrip13.cr18.eta=13; out.etaStrip13.cr18.phi=18;
  out.etaStrip13.cr19.energy=temporary[Eta+13][Phi+19].energy; out.etaStrip13.cr19.eta=13; out.etaStrip13.cr19.phi=19;
  
  out.etaStrip14.cr0.energy=temporary[Eta+14][Phi+0].energy; out.etaStrip14.cr0.eta=14; out.etaStrip14.cr0.phi=0;
  out.etaStrip14.cr1.energy=temporary[Eta+14][Phi+1].energy; out.etaStrip14.cr1.eta=14; out.etaStrip14.cr1.phi=1;
  out.etaStrip14.cr2.energy=temporary[Eta+14][Phi+2].energy; out.etaStrip14.cr2.eta=14; out.etaStrip14.cr2.phi=2;
  out.etaStrip14.cr3.energy=temporary[Eta+14][Phi+3].energy; out.etaStrip14.cr3.eta=14; out.etaStrip14.cr3.phi=3;
  out.etaStrip14.cr4.energy=temporary[Eta+14][Phi+4].energy; out.etaStrip14.cr4.eta=14; out.etaStrip14.cr4.phi=4;
  out.etaStrip14.cr5.energy=temporary[Eta+14][Phi+5].energy; out.etaStrip14.cr5.eta=14; out.etaStrip14.cr5.phi=5;
  out.etaStrip14.cr6.energy=temporary[Eta+14][Phi+6].energy; out.etaStrip14.cr6.eta=14; out.etaStrip14.cr6.phi=6;
  out.etaStrip14.cr7.energy=temporary[Eta+14][Phi+7].energy; out.etaStrip14.cr7.eta=14; out.etaStrip14.cr7.phi=7;
  out.etaStrip14.cr8.energy=temporary[Eta+14][Phi+8].energy; out.etaStrip14.cr8.eta=14; out.etaStrip14.cr8.phi=8;
  out.etaStrip14.cr9.energy=temporary[Eta+14][Phi+9].energy; out.etaStrip14.cr9.eta=14; out.etaStrip14.cr9.phi=9;
  out.etaStrip14.cr10.energy=temporary[Eta+14][Phi+10].energy; out.etaStrip14.cr10.eta=14; out.etaStrip14.cr10.phi=10;
  out.etaStrip14.cr11.energy=temporary[Eta+14][Phi+11].energy; out.etaStrip14.cr11.eta=14; out.etaStrip14.cr11.phi=11;
  out.etaStrip14.cr12.energy=temporary[Eta+14][Phi+12].energy; out.etaStrip14.cr12.eta=14; out.etaStrip14.cr12.phi=12;
  out.etaStrip14.cr13.energy=temporary[Eta+14][Phi+13].energy; out.etaStrip14.cr13.eta=14; out.etaStrip14.cr13.phi=13;
  out.etaStrip14.cr14.energy=temporary[Eta+14][Phi+14].energy; out.etaStrip14.cr14.eta=14; out.etaStrip14.cr14.phi=14;
  out.etaStrip14.cr15.energy=temporary[Eta+14][Phi+15].energy; out.etaStrip14.cr15.eta=14; out.etaStrip14.cr15.phi=15;
  out.etaStrip14.cr16.energy=temporary[Eta+14][Phi+16].energy; out.etaStrip14.cr16.eta=14; out.etaStrip14.cr16.phi=16;
  out.etaStrip14.cr17.energy=temporary[Eta+14][Phi+17].energy; out.etaStrip14.cr17.eta=14; out.etaStrip14.cr17.phi=17;
  out.etaStrip14.cr18.energy=temporary[Eta+14][Phi+18].energy; out.etaStrip14.cr18.eta=14; out.etaStrip14.cr18.phi=18;
  out.etaStrip14.cr19.energy=temporary[Eta+14][Phi+19].energy; out.etaStrip14.cr19.eta=14; out.etaStrip14.cr19.phi=19;
  
  return out;
}

//--------------------------------------------------------// 

// Compare two ecaltp_t and return the one with the larger pT.
ecaltp_t bestOf2(const ecaltp_t& ecaltp0, const ecaltp_t& ecaltp1) {
  ecaltp_t x;
  x = (ecaltp0.energy > ecaltp1.energy)?ecaltp0:ecaltp1;

  return x;
}

//--------------------------------------------------------// 

// For a given etaStrip_t, find the ecaltp_t (out of 20 of them) with the largest pT, using pairwise comparison
ecaltp_t getPeakBin20N(const etaStrip_t& etaStrip){

  ecaltp_t best01 = bestOf2(etaStrip.cr0,etaStrip.cr1) ;
  ecaltp_t best23 = bestOf2(etaStrip.cr2,etaStrip.cr3) ;
  ecaltp_t best45 = bestOf2(etaStrip.cr4,etaStrip.cr5) ;
  ecaltp_t best67 = bestOf2(etaStrip.cr6,etaStrip.cr7) ;
  ecaltp_t best89 = bestOf2(etaStrip.cr8,etaStrip.cr9) ;
  ecaltp_t best1011 = bestOf2(etaStrip.cr10,etaStrip.cr11) ;
  ecaltp_t best1213 = bestOf2(etaStrip.cr12,etaStrip.cr13) ;
  ecaltp_t best1415 = bestOf2(etaStrip.cr14,etaStrip.cr15) ;
  ecaltp_t best1617 = bestOf2(etaStrip.cr16,etaStrip.cr17) ;
  ecaltp_t best1819 = bestOf2(etaStrip.cr18,etaStrip.cr19) ;

  ecaltp_t best0123 = bestOf2(best01,best23) ;
  ecaltp_t best4567 = bestOf2(best45,best67) ;
  ecaltp_t best891011 = bestOf2(best89,best1011) ;
  ecaltp_t best12131415 = bestOf2(best1213,best1415) ;
  ecaltp_t best16171819 = bestOf2(best1617,best1819) ;

  ecaltp_t best01234567 = bestOf2(best0123,best4567) ;
  ecaltp_t best89101112131415 = bestOf2(best891011,best12131415) ;

  ecaltp_t best0to15 = bestOf2(best01234567,best89101112131415) ;
  ecaltp_t bestOf20 = bestOf2(best0to15,best16171819) ;

  return bestOf20 ;
}

//--------------------------------------------------------//

// For a given etaStripPeak_t (representing the 15 crystals, one per row in eta, not necessarily with the same phi),
// return the crystal with the highest pT).
// This is used in getClusterPosition to get the cluster seed.

crystalMax getPeakBin15N(const etaStripPeak_t& etaStrip){
  crystalMax x;

  ecaltp_t best01 = bestOf2(etaStrip.pk0,etaStrip.pk1) ;
  ecaltp_t best23 = bestOf2(etaStrip.pk2,etaStrip.pk3) ;
  ecaltp_t best45 = bestOf2(etaStrip.pk4,etaStrip.pk5) ;
  ecaltp_t best67 = bestOf2(etaStrip.pk6,etaStrip.pk7) ;
  ecaltp_t best89 = bestOf2(etaStrip.pk8,etaStrip.pk9) ;
  ecaltp_t best1011 = bestOf2(etaStrip.pk10,etaStrip.pk11) ;
  ecaltp_t best1213 = bestOf2(etaStrip.pk12,etaStrip.pk13) ;

  ecaltp_t best0123 = bestOf2(best01,best23) ;
  ecaltp_t best4567 = bestOf2(best45,best67) ;
  ecaltp_t best891011 = bestOf2(best89,best1011) ;
  ecaltp_t best121314 = bestOf2(best1213,etaStrip.pk14) ;

  ecaltp_t best01234567 = bestOf2(best0123,best4567);
  ecaltp_t best891011121314 = bestOf2(best891011,best121314) ;

  ecaltp_t bestOf15 = bestOf2(best01234567,best891011121314) ;

  x.energy = bestOf15.energy ;
  x.etaMax = bestOf15.eta ;
  x.phiMax = bestOf15.phi ;

  /* std::cout << "[---] getPeakBin15N: bestOf15.energy, bestOf15.etaMax, bestOf15.phiMax: "  */
  /*       << bestOf15.energy << ", "  */
  /*       << bestOf15.eta    << ", " */
  /* 	    << bestOf15.phi    << std::endl; */
  // << x.energy << ", "                                                                                            
  // << x.etaMax    << ", "                                                                                            
  // << x.phiMax    << std::endl;     

  return x ;
}

//--------------------------------------------------------// 

// Take a 3x4 ECAL region (i.e. 15x20 in crystals, add crystal energies in squares of 5x5, giving
// 3x4 = 12 ECAL tower sums.) Store these 12 values in towerEt.

void getECALTowersEt(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<12> towerEt[12]){

  ap_uint<10> temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI] ;
  ap_uint<12> towerEtN[3][4][5] ;
  for (int i = 0; i < CRYSTAL_IN_ETA; i++){
    for(int k = 0; k < CRYSTAL_IN_PHI; k++){
      temp[i][k] = tempX[i][k].energy ;
    }}
    
  for(int i=0; i<CRYSTAL_IN_ETA; i=i+5){
    for(int k=0; k<CRYSTAL_IN_PHI; k=k+5){
      towerEtN[i/5][k/5][0] = temp[i][k] + temp[i][k+1] + temp[i][k+2] + temp[i][k+3] + temp[i][k+4] ;
      towerEtN[i/5][k/5][1] = temp[i+1][k] + temp[i+1][k+1] + temp[i+1][k+2] + temp[i+1][k+3] + temp[i+1][k+4] ;
      towerEtN[i/5][k/5][2] = temp[i+2][k] + temp[i+2][k+1] + temp[i+2][k+2] + temp[i+2][k+3] + temp[i+2][k+4] ;
      towerEtN[i/5][k/5][3] = temp[i+3][k] + temp[i+3][k+1] + temp[i+3][k+2] + temp[i+3][k+3] + temp[i+3][k+4] ;
      towerEtN[i/5][k/5][4] = temp[i+4][k] + temp[i+4][k+1] + temp[i+4][k+2] + temp[i+4][k+3] + temp[i+4][k+4] ;
    }}

  towerEt[0]= towerEtN[0][0][0] + towerEtN[0][0][1] + towerEtN[0][0][2] + towerEtN[0][0][3] + towerEtN[0][0][4] ;
  towerEt[1]= towerEtN[0][1][0] + towerEtN[0][1][1] + towerEtN[0][1][2] + towerEtN[0][1][3] + towerEtN[0][1][4] ;
  towerEt[2]= towerEtN[0][2][0] + towerEtN[0][2][1] + towerEtN[0][2][2] + towerEtN[0][2][3] + towerEtN[0][2][4] ;
  towerEt[3]= towerEtN[0][3][0] + towerEtN[0][3][1] + towerEtN[0][3][2] + towerEtN[0][3][3] + towerEtN[0][3][4] ;
  towerEt[4]= towerEtN[1][0][0] + towerEtN[1][0][1] + towerEtN[1][0][2] + towerEtN[1][0][3] + towerEtN[1][0][4] ;
  towerEt[5]= towerEtN[1][1][0] + towerEtN[1][1][1] + towerEtN[1][1][2] + towerEtN[1][1][3] + towerEtN[1][1][4] ;
  towerEt[6]= towerEtN[1][2][0] + towerEtN[1][2][1] + towerEtN[1][2][2] + towerEtN[1][2][3] + towerEtN[1][2][4] ;
  towerEt[7]= towerEtN[1][3][0] + towerEtN[1][3][1] + towerEtN[1][3][2] + towerEtN[1][3][3] + towerEtN[1][3][4] ;
  towerEt[8]= towerEtN[2][0][0] + towerEtN[2][0][1] + towerEtN[2][0][2] + towerEtN[2][0][3] + towerEtN[2][0][4] ;
  towerEt[9]= towerEtN[2][1][0] + towerEtN[2][1][1] + towerEtN[2][1][2] + towerEtN[2][1][3] + towerEtN[2][1][4] ;
  towerEt[10]= towerEtN[2][2][0] + towerEtN[2][2][1] + towerEtN[2][2][2] + towerEtN[2][2][3] + towerEtN[2][2][4] ;
  towerEt[11]= towerEtN[2][3][0] + towerEtN[2][3][1] + towerEtN[2][3][2] + towerEtN[2][3][3] + towerEtN[2][3][4] ;

  std::cout << "getECALTowersEt (after subtracting clusters), total energy in this 3x4 region: ";
  ap_uint<12> totalEt;
  for (int i = 0; i < 12; i++) {
    totalEt += towerEt[i];
  }
  std::cout << totalEt/8.0 << " GeV" << std::endl;
    
}

//--------------------------------------------------------//  

class clusterInfo{
 public:
  ap_uint<10> seedEnergy;
  ap_uint<15> energy;
  ap_uint<15> et5x5;
  ap_uint<15> et2x5;
  ap_uint<5> phiMax;
  ap_uint<5> etaMax;
  ap_uint<2> brems;

  clusterInfo(){
    seedEnergy = 0;
    energy = 0;
    et5x5 = 0;
    et2x5 = 0;
    phiMax = 0;
    etaMax = 0;
    brems = 0;
  }

  clusterInfo& operator=(const clusterInfo& rhs){
    seedEnergy = rhs.seedEnergy;
    energy = rhs.energy;
    et5x5 = rhs.et5x5;
    et2x5 = rhs.et2x5;
    phiMax = rhs.phiMax;
    etaMax = rhs.etaMax;
    brems = rhs.brems;
    return *this;
  }
};

//--------------------------------------------------------//

class Cluster{
 public:
  ap_uint<28> data;
  int regionIdx;    // region is 0 through 4 in barrel, 5 for overlap in barrel+endcap
  float calib;     //  ECAL energy calibration factor
  ap_uint<2> brems;       // brems flag
  ap_uint<15> et5x5;      // et5x5 sum in towers around cluster
  ap_uint<15> et2x5;      // et2x5 sum in towers around cluster
  bool is_ss;         // shower shape flag
  bool is_looseTkss;  // loose Tk shower shape flag

  Cluster(){
    data = 0;
    regionIdx = -1;
    calib = 1.0;
    brems = 0;
    et5x5 = 0;
    et2x5 = 0;
    is_ss = false;
    is_looseTkss = false;
  }

  Cluster(ap_uint<12> clusterEnergy, ap_uint<5> towerEta, ap_uint<2> towerPhi, ap_uint<3> clusterEta, ap_uint<3> clusterPhi, ap_uint<3> satur, 
          ap_uint<15> clusterEt5x5 = 0, ap_uint<15> clusterEt2x5 = 0, ap_uint<2> clusterBrems = 0, 
          bool cluster_is_ss = false, bool cluster_is_looseTkss = false, int clusterRegionIdx = 0) {
        data = (clusterEnergy) | 
            (((ap_uint<32>) towerEta)  << 12) |
            (((ap_uint<32>) towerPhi)  << 17) |
            (((ap_uint<32>) clusterEta)  << 19) |
            (((ap_uint<32>) clusterPhi) << 22) |
            (((ap_uint<32>) satur)       << 25);
    regionIdx = clusterRegionIdx,
    et5x5 = clusterEt5x5;
    et2x5 = clusterEt2x5;
    brems = clusterBrems;
    is_ss = cluster_is_ss;
    is_looseTkss = cluster_is_looseTkss;
  }

  Cluster& operator=(const Cluster& rhs){
    data      = rhs.data;
    regionIdx = rhs.regionIdx;
    calib     = rhs.calib;
    brems     = rhs.brems;
    et5x5     = rhs.et5x5;
    et2x5     = rhs.et2x5;
    is_ss     = rhs.is_ss;
    is_looseTkss = rhs.is_looseTkss;
    return *this;
  }

  void setRegionIdx(int regIdx) { regionIdx = regIdx; }  // Newly added

  ap_uint<12> clusterEnergy() {return (data & 0xFFF);}
  ap_uint<5> towerEta() {return ((data >> 12) & 0x1F);}  // goes from 0 to 3 (need region for full info)
  ap_uint<2> towerPhi() {return ((data >> 17) & 0x3);}
  ap_uint<3> clusterEta() {return ((data >> 19) & 0x7);}
  ap_uint<3> clusterPhi() {return ((data >> 22) & 0x7);}
  ap_uint<3> satur() {return ((data >> 25) & 0x7);}
  ap_uint<15> uint_et2x5(){ return et2x5; }
  ap_uint<15> uint_et5x5(){ return et5x5; }
  
  operator uint32_t() {return (ap_uint<28>) data;}
  int region() { return regionIdx; }                       // Newly added
  int getBrems() { return (int) brems; }              
  float getCalib() { return (float) calib; }                    
  float getPt()    { return ((float) clusterEnergy() / 8.0); }   // Return pT as a float
  float getEt5x5() { return ((float) et5x5 / 8.0); }             // Return ET5x5 as a float
  float getEt2x5() { return ((float) et2x5 / 8.0); }             // Return ET2x5 as a float
  
  int towerEtaInCard() { return ((int) (region() * TOWER_IN_ETA) + towerEta() ); } 

  bool getIsSS() { return is_ss; }
  bool getIsLooseTkss() { return is_looseTkss; }

  /*
   * Apply calibration (float) to the pT in-place.
   */
  void applyCalibration(float factor) {  
    
    float newPt = getPt() * factor;
    // Convert the new pT to an unsigned int, 28 bits to take the logical OR with the mask
    ap_uint<28> newPt_uint = (ap_uint<28>) (int) (newPt * 8.0);
    // Make sure that the new pT only takes up the last twelve bits
    newPt_uint = (newPt_uint & 0x0000FFF);
    
    // std::cout << "   newPt_uint: " << (int) newPt_uint << std::endl;
    // Modify 'data'
    ap_uint<28> bitMask = 0xFFFF000;  // last twelve digits are zero
    data = (data & bitMask);                 // zero out the last twelve digits
    data = (data | newPt_uint);              // write in the new ET
    // std::cout << "   new: " << (int) data << std::endl;
  }

  /*
   * Print a cluster's information.
   */  
  void printInfo(std::string description) {
    std::cout << "[Cluster::printInfo:] [" << description << "]:" 
              << "pT (float) " << getPt() << ", "
              << "region: " << region() << "," 
              << "towerEta and Phi:  " << towerEta() << ", " << towerPhi() << "), "
              << "clusterEta and Phi: (" << clusterEta() << ", " << clusterPhi() << ") "
              << "et2x5 and et5x5 (float): " << getEt2x5() << ", " << getEt5x5() << ")" 
              << std::endl;
  }

  // Get real eta
  float realEta(int cc) {
    return getEta_fromCrystaliEta(getCrystal_iEta_fromCardRegionInfo(cc, region(), towerEta(), clusterEta()) );
  }

  // Get real phi
  float realPhi(int cc) {
    return getPhi_fromCrystaliPhi(getCrystal_iPhi_fromCardRegionInfo(cc, region(), towerPhi(), clusterPhi()) );
  }




};

//--------------------------------------------------------// 

/*
 * Compare the ET of two clusters (pass this to std::sort to get clusters sorted in decreasing ET).
 */

bool compareClusterET(Cluster& lhs, Cluster& rhs) {
  return ( lhs.clusterEnergy() > rhs.clusterEnergy() );
}

//--------------------------------------------------------//

clusterInfo getClusterPosition(const ecalRegion_t& ecalRegion){
  etaStripPeak_t etaStripPeak;
  clusterInfo cluster;

  etaStripPeak.pk0 = getPeakBin20N(ecalRegion.etaStrip0);
  etaStripPeak.pk1 = getPeakBin20N(ecalRegion.etaStrip1);
  etaStripPeak.pk2 = getPeakBin20N(ecalRegion.etaStrip2);
  etaStripPeak.pk3 = getPeakBin20N(ecalRegion.etaStrip3);
  etaStripPeak.pk4 = getPeakBin20N(ecalRegion.etaStrip4);
  etaStripPeak.pk5 = getPeakBin20N(ecalRegion.etaStrip5);
  etaStripPeak.pk6 = getPeakBin20N(ecalRegion.etaStrip6);
  etaStripPeak.pk7 = getPeakBin20N(ecalRegion.etaStrip7);
  etaStripPeak.pk8 = getPeakBin20N(ecalRegion.etaStrip8);
  etaStripPeak.pk9 = getPeakBin20N(ecalRegion.etaStrip9);
  etaStripPeak.pk10 = getPeakBin20N(ecalRegion.etaStrip10);
  etaStripPeak.pk11 = getPeakBin20N(ecalRegion.etaStrip11);
  etaStripPeak.pk12 = getPeakBin20N(ecalRegion.etaStrip12);
  etaStripPeak.pk13 = getPeakBin20N(ecalRegion.etaStrip13);
  etaStripPeak.pk14 = getPeakBin20N(ecalRegion.etaStrip14);

  crystalMax peakIn15;
  peakIn15 = getPeakBin15N(etaStripPeak);
  
  cluster.seedEnergy = peakIn15.energy;
  cluster.energy = 0;
  cluster.etaMax = peakIn15.etaMax;
  cluster.phiMax = peakIn15.phiMax;
  cluster.brems = 0;
  cluster.et5x5 = 0;
  cluster.et2x5 = 0;

  return cluster;
}

//--------------------------------------------------------//

/*
 * Return initialized cluster with specified Et, eta, phi, with all other fields (saturation, Et2x5, Et5x5, brems, flags initialized to 0/ false).
 */ 
Cluster packCluster(ap_uint<15>& clusterEt, ap_uint<5>& etaMax_t, ap_uint<5>& phiMax_t){

  std::cout << "[-->] packCluster: clusterEt (GeV) " << clusterEt/8.0 << ", eta and phi: " << etaMax_t << ", " << phiMax_t << std::endl;
  ap_uint<12> peggedEt;
  Cluster pack;

  ap_uint<5> towerEta = (etaMax_t)/5;
  ap_uint<2> towerPhi = (phiMax_t)/5;
  ap_uint<3> clusterEta = etaMax_t - 5*towerEta;
  ap_uint<3> clusterPhi = phiMax_t - 5*towerPhi;

  peggedEt = (clusterEt > 0xFFF)? (ap_uint<12>)0xFFF : (ap_uint<12>) clusterEt;

  pack = Cluster(peggedEt, towerEta, towerPhi, clusterEta, clusterPhi, 0);

  return pack;
}


//--------------------------------------------------------// 

// Given the cluster seed_eta, seed_phi, and brems, remove the cluster energy
// from the given crystal array temp. Functionally identical to "RemoveTmp".

void removeClusterFromCrystal(crystal temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi, ap_uint<2> brems) {

  // Zero out the crystal energies in a 3 (eta) by 5 (phi) window (the clusters are 3x5 in crystals)
  for (int i = 0; i < CRYSTAL_IN_ETA; i++){
    for (int k = 0; k < CRYSTAL_IN_PHI; k++){
      if ((i >= (seed_eta-1)) && (i <= (seed_eta+1)) && (k >= (seed_phi-2)) && (k <= (seed_phi+2))) {
        temp[i][k].energy = 0;
      }
    }
  }

  // If brems flag is 1, *also* zero the energies in the 3x5 window to the "left" of the cluster 
  // N.B. in the positive eta cards, "left" in the region = towards negative phi, 
  // but for negative eta cards, everything is flipped, so "left" in the region" = towards positive phi
  if (brems == 1) {
    for (int i = 0; i < CRYSTAL_IN_ETA; i++){
      for (int k = 0; k < CRYSTAL_IN_PHI; k++){
	      if ((i >= (seed_eta-1)) && (i <= (seed_eta+1)) && (k >= (seed_phi-2-5)) && (k <= (seed_phi+2-5))) {
          temp[i][k].energy = 0;
        } 
      }
    }
  }  
  // If brems flag is 2, *also* zero the energies in the 3x5 window to the "right" of the cluster
  // N.B. in the positive eta cards, "right" in the region = towards POSITIVE phi,
  // but for negative eta cards, everything is flipped, so "right" in the region = towards NEGATIVE phi
  else if (brems == 2) {
    for (int i = 0; i < CRYSTAL_IN_ETA; i++){
      for (int k = 0; k < CRYSTAL_IN_PHI; k++){
	      if ((i >= (seed_eta-1)) && (i <= (seed_eta+1)) && (k >= (seed_phi-2+5)) && (k <= (seed_phi+2+5))) {
          temp[i][k].energy = 0;
        }
      }
    }
  }
  
}

//--------------------------------------------------------//

// Given a 15x20 crystal tempX, and a seed with seed_eta and seed_phi, return a clusterInfo containing
// the cluster energy for a positive bremmstrahulung shift 

clusterInfo getBremsValuesPos(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi ){

  ap_uint<12> temp[CRYSTAL_IN_ETA+2][CRYSTAL_IN_PHI+4];
  ap_uint<12> phi0eta[3], phi1eta[3], phi2eta[3], phi3eta[3], phi4eta[3];
  ap_uint<12> eta_slice[3];
  clusterInfo cluster_tmp;
  
  // Set all entries in a new ((15+2)x(20+4)) array to be zero.
  for (int i = 0; i < (CRYSTAL_IN_ETA + 2); i++) {
    for (int j = 0; j < (CRYSTAL_IN_PHI + 4); j++) {
      temp[i][j] = 0;
    }
  }
  
  // Read the energies of the input crystal tempX into the slightly larger array temp, with an offset so temp is tempX
  // except shifted +1 in eta, and -3 in phi. 
  for (int i = 0; i < (CRYSTAL_IN_ETA); i++) {
    for (int j = 0; j < (CRYSTAL_IN_PHI-3); j++){
      // std::cout << "getBremsValuesPos: i, j " << i << ", " << j << ", [i+1][j-3]: " << i+1 << ", " << j-3 << std::endl;
      temp[i+1][j] = tempX[i][j+3].energy ;
      if (tempX[i][j].energy > 0) {
        // std::cout << "For seed_eta " << seed_eta << " seed phi " << seed_phi << ", " 
        //     << "getBremsValuesPos: tempX[" << i << "][" << j+3 << "].energy: " << tempX[i][j+3].energy
        //     << ", temp[" << i+1 << "][" << j << "]" << temp[i+1][j] << std::endl;
      }
    }
  }
  
  ap_uint<6> seed_eta1, seed_phi1;
  seed_eta1 = seed_eta; //to start from corner
  seed_phi1 = seed_phi; //to start from corner

  // now we are in the left bottom corner 
  // Loop over the shifted array, and at the original location of the seed (seed_eta1/seed_phi1), 
  // read a 3 (eta) x 5 (phi) rectangle of crystals where the original location of the seed is in the bottom left corner
  for (int j = 0; j < CRYSTAL_IN_ETA; j++) {
    if (j == seed_eta1) {
      for (int k = 0; k < CRYSTAL_IN_PHI; k++) {
        if (k == seed_phi1) {
          // Same eta as the seed, read next five crystals in phi
          // std::cout << "j, k = " << j << ", " << k << std::endl;
          
          phi0eta[0] = temp[j][k];
          phi1eta[0] = temp[j][k+1];
          phi2eta[0] = temp[j][k+2];
          phi3eta[0] = temp[j][k+3];
          phi4eta[0] = temp[j][k+4];

          // for (int kk = k; kk < (k+5); kk++) {
          //   std::cout << "--> phi{0-4}eta[0]:" << temp[j][kk] << std::endl;
          // }

          // +1 eta from the seed, read next five crystals in phi
          phi0eta[1] = temp[j+1][k];
          phi1eta[1] = temp[j+1][k+1];
          phi2eta[1] = temp[j+1][k+2];
          phi3eta[1] = temp[j+1][k+3];
          phi4eta[1] = temp[j+1][k+4];

          // for (int kk =k; kk <(k+5); kk++) {
          //   std::cout << "--> phi{0-4}eta[1]:"<< temp[j+1][kk] << std::endl;
          // }
                              
          // +2 eta from the seed, read next five crystals in phi
          phi0eta[2] = temp[j+2][k];
          phi1eta[2] = temp[j+2][k+1];
          phi2eta[2] = temp[j+2][k+2];
          phi3eta[2] = temp[j+2][k+3];
          phi4eta[2] = temp[j+2][k+4];

          // for (int kk = k; kk < (k+5); kk++) {
          //   std::cout << "--> phi{0-4}eta[2]:"<< temp[j+2][kk] << std::endl;
          // }
	        continue;
        }
      }
    }
  }

  // Add up the energies in this 3x5 of crystals, initialize a cluster_tmp, and return it
  for (int i = 0; i < 3; i++) { 
    eta_slice[i] = phi0eta[i] + phi1eta[i] + phi2eta[i] + phi3eta[i] + phi4eta[i];
    // std::cout << ">>> i = " << i << ": "
	  //           << phi0eta[i] << " " << phi1eta[i] << " " << phi2eta[i] << " " << phi3eta[i] << " " << phi4eta[i] << std::endl;
  }
  cluster_tmp.energy = (eta_slice[0] + eta_slice[1] + eta_slice[2]);

  // std::cout << "getBremsValuesPos: energy, seed eta/phi = " << cluster_tmp.energy << ", " << seed_eta << ", " << seed_phi << std::endl;
  return cluster_tmp;

}

//--------------------------------------------------------//  

// Given a 15x20 crystal tempX, and a seed with seed_eta and seed_phi, return a clusterInfo containing                                     
// the cluster energy for a *negative* bremmstrahlung shift 

clusterInfo getBremsValuesNeg(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi){

  ap_uint<12> temp[CRYSTAL_IN_ETA+2][CRYSTAL_IN_PHI+4];
  ap_uint<12> phi0eta[3], phi1eta[3], phi2eta[3], phi3eta[3], phi4eta[3];

  ap_uint<12> eta_slice[3];

  clusterInfo cluster_tmp;

  // Initialize all entries in a new ((15+2)x(20+4)) array to be zero.                                                                  
  for (int i = 0; i < (CRYSTAL_IN_ETA + 2); i++) {
    for (int j = 0; j < (CRYSTAL_IN_PHI + 4); j++) {
      temp[i][j] = 0;
    }
  }

  // Read the energies of the input crystal tempX into the slightly larger array temp, with an offset so temp is tempX
  // except shifted in +1 in eta and +7 in phi
  for (int i = 0; i < (CRYSTAL_IN_ETA); i++) {
    for (int j = 0; j < (CRYSTAL_IN_PHI - 1); j++) {
      temp[i+1][j+7] = tempX[i][j].energy; // fixed from firmware
      if (tempX[i][j].energy > 0) {
	      // std::cout << "For seed_eta " << seed_eta << " seed phi " << seed_phi << ", "
		    //           << "getBremsValuesNeg: tempX[" << i << "][" << j << "].energy: " << tempX[i][j].energy << std::endl;
      }
    }
  }

  ap_uint<6> seed_eta1, seed_phi1;
  seed_eta1 = seed_eta; //to start from corner
  seed_phi1 = seed_phi; //to start from corner

  // Loop over the shifted array, and at the original location of the seed (seed_eta1/seed_phi1), 
  // read a 3 (eta) x 5 (phi) rectangle of crystals where the original location of the seed is in the bottom left corner
  for (int j = 0; j < CRYSTAL_IN_ETA; j++) {
    if (j == seed_eta1) {
      for (int k = 0; k < CRYSTAL_IN_PHI; k++) {
        if (k == seed_phi1) {
          // Same eta as the seed, read next five crystals in phi
          phi0eta[0] = temp[j][k];
          phi1eta[0] = temp[j][k+1];
          phi2eta[0] = temp[j][k+2];
          phi3eta[0] = temp[j][k+3];
          phi4eta[0] = temp[j][k+4];
              
          // +1 eta from the seed, read next five crystals in phi
          phi0eta[1] = temp[j+1][k];
          phi1eta[1] = temp[j+1][k+1];
          phi2eta[1] = temp[j+1][k+2];
          phi3eta[1] = temp[j+1][k+3];
          phi4eta[1] = temp[j+1][k+4];
                              
          // +2 eta from the seed, read next five crystals in phi
          phi0eta[2] = temp[j+2][k];
          phi1eta[2] = temp[j+2][k+1];
          phi2eta[2] = temp[j+2][k+2];
          phi3eta[2] = temp[j+2][k+3];
          phi4eta[2] = temp[j+2][k+4];              
          continue;
        }
      }
    }
  }
  
  // Add up the energies in this 3x5 of crystals, initialize a cluster_tmp, and return it
  for (int i = 0; i < 3; i++) { 
    eta_slice[i] = phi0eta[i] + phi1eta[i] + phi2eta[i] + phi3eta[i] + phi4eta[i];
  }
  cluster_tmp.energy = (eta_slice[0] + eta_slice[1] + eta_slice[2]);
  
  // std::cout << "getBremsValuesNeg: energy, seed eta/phi = " << cluster_tmp.energy << ", " << seed_eta << ", " << seed_phi << std::endl;
  
  return cluster_tmp;
  
}


//--------------------------------------------------------//                                                                           

// Given a 15x20 crystal tempX, and a seed with seed_eta and seed_phi, return a clusterInfo containing                              
// the cluster energy (central value)

clusterInfo getClusterValues(crystal tempX[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI], ap_uint<5> seed_eta,  ap_uint<5> seed_phi){

  ap_uint<12> temp[CRYSTAL_IN_ETA+4][CRYSTAL_IN_PHI+4];
  ap_uint<12> phi0eta[5], phi1eta[5], phi2eta[5], phi3eta[5], phi4eta[5];
  ap_uint<12> eta_slice[5];
  ap_uint<12> et2x5_1Tot, et2x5_2Tot, etSum2x5;
  ap_uint<12> et5x5Tot;

  clusterInfo cluster_tmp;
  // Initialize empty (15+4)x(20+4) array
  for (int i = 0; i < (CRYSTAL_IN_ETA + 4); i++){
    for (int k = 0; k< (CRYSTAL_IN_PHI + 4); k++){
      temp[i][k] = 0 ;
    }
  }

  // Copy input array energies into temp array with +2 eta and +2 phi offset.
  // std::cout << ">>> [INFO: Phase2L1RCT.h: getClusterValues]: entries in the input 15x20 array are: ";
  for (int i = 0; i < (CRYSTAL_IN_ETA); i++){
    for (int k = 0; k < (CRYSTAL_IN_PHI); k++){
      temp[i+2][k+2] = tempX[i][k].energy ;
      // std::cout << tempX[i][k].energy << ", ";
    }
  }
 // std::cout << std::endl;

  ap_uint<6> seed_eta1, seed_phi1;
  seed_eta1 = seed_eta; //to start from corner
  seed_phi1 = seed_phi; //to start from corner
  
  // now we are in the left bottom corner 
  // Loop over the shifted array, and at the original location of the seed (seed_eta1/seed_phi1), 
  // read a 5 (eta) x 5 (phi) rectangle of crystals where the original location of the seed is in the bottom left corner
  for (int j = 0; j < CRYSTAL_IN_ETA; j++) {
    if (j == seed_eta1) {
      for (int k = 0; k < CRYSTAL_IN_PHI; k++) {
	      if (k == seed_phi1) {
	      // Same eta as the seed, read next five crystals in phi
	      phi0eta[0] = temp[j][k];
	      phi1eta[0] = temp[j][k+1];
	      phi2eta[0] = temp[j][k+2];
	      phi3eta[0] = temp[j][k+3];
	      phi4eta[0] = temp[j][k+4];
	      
	      // +1 eta from the seed, read next five crystals in phi
	      phi0eta[1] = temp[j+1][k];
	      phi1eta[1] = temp[j+1][k+1];
	      phi2eta[1] = temp[j+1][k+2];
	      phi3eta[1] = temp[j+1][k+3];
	      phi4eta[1] = temp[j+1][k+4];

	      // +2 eta from the seed, read next five crystals in phi
	      phi0eta[2] = temp[j+2][k];
	      phi1eta[2] = temp[j+2][k+1];
	      phi2eta[2] = temp[j+2][k+2];
	      phi3eta[2] = temp[j+2][k+3];
	      phi4eta[2] = temp[j+2][k+4];

	      // +3 eta from the seed, read next five crystals in phi 
	      phi0eta[3] = temp[j+3][k];
	      phi1eta[3] = temp[j+3][k+1];
	      phi2eta[3] = temp[j+3][k+2];
	      phi3eta[3] = temp[j+3][k+3];
	      phi4eta[3] = temp[j+3][k+4];

	      // +4 eta from the seed, read next five crystals in phi
	      phi0eta[4] = temp[j+4][k];
	      phi1eta[4] = temp[j+4][k+1];
	      phi2eta[4] = temp[j+4][k+2];
	      phi3eta[4] = temp[j+4][k+3];
	      phi4eta[4] = temp[j+4][k+4];
	      
	      continue;
        }
      }
    }
  }

  // Add the first three eta strips into the cluster energy
  for (int i = 0; i < 5; i++) {
    eta_slice[i] = phi0eta[i] + phi1eta[i] + phi2eta[i] + phi3eta[i] + phi4eta[i];
  }
  cluster_tmp.energy = (eta_slice[1] + eta_slice[2] + eta_slice[3]);

    // Print the energy in the slices 
  //std::cout << ">>> [INFO: Phase2L1RCT.h: getClusterValues]: The energy in the five eta slices is: ";
  for (unsigned int i = 0; i < 5; i++) {
   // std::cout << "Slice " << i << ": " << eta_slice[i] << ", ";
  }
  // std::cout << std::endl;

  // Get the energy totals in the 5x5 and also in two 2x5 
  et5x5Tot   = (eta_slice[0] + eta_slice[1] + eta_slice[2] + eta_slice[3] + eta_slice[4]);
  et2x5_1Tot = (eta_slice[1] + eta_slice[2]);
  et2x5_2Tot = (eta_slice[2] + eta_slice[3]);

  // std::cout << ">>> [INFO: Phase2L1RCT.h: getClusterValues]: et2x5_1Tot (slice 1+2) " << et2x5_1Tot << " and et2x5_2Tot (slice 2+3): " << et2x5_2Tot << std::endl;

  if (et2x5_1Tot >= et2x5_2Tot) etSum2x5 = et2x5_1Tot;
  else                          etSum2x5 = et2x5_2Tot;

  cluster_tmp.et5x5 = et5x5Tot;
  cluster_tmp.et2x5 = etSum2x5;

  // std::cout << ">>> [INFO: Phase2L1RCT.h: getClusterValues:] output energy, et5x5Tot (eta slices 0+1+2+3+4), etSum2x5 (eta slices 1+2 or 2+3 whichever is greater) = " << cluster_tmp.energy << ", " 
	//           << cluster_tmp.et5x5 << ", " << cluster_tmp.et2x5 << std::endl;

  return cluster_tmp;
}

//--------------------------------------------------------// 

// Functionally identical to "getRegion3x4" from algo_top.cc (Renamed to avoid confusion with the card class method)
// In 15x20 crystal array temp, return the next cluster, and remove the cluster's energy
// from the crystal array.

Cluster getClusterFromRegion3x4(crystal temp[CRYSTAL_IN_ETA][CRYSTAL_IN_PHI]){

  Cluster returnCluster;
  clusterInfo cluster_tmp;
  clusterInfo cluster_tmpCenter;
  clusterInfo cluster_tmpBneg;
  clusterInfo cluster_tmpBpos;

  ecalRegion_t ecalRegion;
  ecalRegion = initStructure(temp);

  cluster_tmp = getClusterPosition(ecalRegion);

  float seedEnergyFloat = cluster_tmp.seedEnergy/8.0;

  // Do not make cluster if seed is less than 1.0 GeV
  if (seedEnergyFloat < 1.0) {  
    cluster_tmp.energy = 0;
    cluster_tmp.phiMax = 0;
    cluster_tmp.etaMax = 0;
    return packCluster(cluster_tmp.energy, cluster_tmp.phiMax, cluster_tmp.etaMax);
  }
  
  ap_uint<5> seed_phi = cluster_tmp.phiMax;
  ap_uint<5> seed_eta = cluster_tmp.etaMax;
  
  cluster_tmpCenter = getClusterValues(temp, seed_eta, seed_phi);
  cluster_tmpBneg   = getBremsValuesNeg(temp, seed_eta, seed_phi); 
  cluster_tmpBpos   = getBremsValuesPos(temp, seed_eta, seed_phi);

  cluster_tmp.energy = cluster_tmpCenter.energy;
  cluster_tmp.brems = 0;

  // Create a cluster 
  if ((cluster_tmpBneg.energy > cluster_tmpCenter.energy/8) && (cluster_tmpBneg.energy > cluster_tmpBpos.energy)) {   
    cluster_tmp.energy = (cluster_tmpCenter.energy + cluster_tmpBneg.energy);
    // std::cout << "getClusterFromRegion3x4: Brems in negative phi direction: set cluster ET to " << cluster_tmp.energy << std::endl;
    cluster_tmp.brems = 1; }
  else if(cluster_tmpBpos.energy > cluster_tmpCenter.energy/8) {
    cluster_tmp.energy = (cluster_tmpCenter.energy + cluster_tmpBpos.energy);
    // std::cout << "getClusterFromRegion3x4: Brems in positive phi direction: set cluster ET to " << cluster_tmp.energy << std::endl;
    cluster_tmp.brems = 2; }
  
  returnCluster = packCluster(cluster_tmp.energy, cluster_tmp.etaMax, cluster_tmp.phiMax);
  removeClusterFromCrystal(temp, seed_eta, seed_phi, cluster_tmp.brems);
  
  // Add clusterInfo members to the output cluster members
  returnCluster.brems = cluster_tmp.brems; 
  returnCluster.et5x5 = cluster_tmpCenter.et5x5;  // get et5x5 from the center value
  returnCluster.et2x5 = cluster_tmpCenter.et2x5;  // get et2x5 from the center value
 //  std::cout << "[INFO: Phase2L1RCT.h: getClusterFromRegion3x4():] et2x5 and et5x5 are " << returnCluster.et2x5 << ", " << returnCluster.et5x5 << std::endl;
  return returnCluster;
  
}

//--------------------------------------------------------// 

// Stitch clusters in cluster_list across the boundary specified by 
// towerEtaUpper and towerEtaLower (using RCT card notation). Modifies the input vector
// (passed by reference). If two clusters are combined, modify the higher-energy cluster and
// zero out the energy of the smaller-energy cluster.
// cc is the RCT card number (for print-out statements only).
// Returns 1.

int stitchClusterOverRegionBoundary(std::vector<Cluster>& cluster_list, int towerEtaUpper, int towerEtaLower, int cc) {

  (void) cc; // for printout statements

  int crystalEtaUpper = 0;
  int crystalEtaLower = 4;

  //std::cout << "Comparing for towerEta boundaries " << towerEtaUpper << " and " << towerEtaLower << std::endl;
  for (int i = 0; i < (int) cluster_list.size(); i++) {
    for (int j = 0; j < (int) cluster_list.size(); j++) {
      // Do not double-count
      if (i == j) continue; 

      Cluster c1 = cluster_list[i];
      Cluster c2 = cluster_list[j];

      // Use the .towerEtaInCard() method to get the tower eta in the entire RCT card
      if ( (c1.clusterEnergy() > 0) && (c1.towerEtaInCard() == towerEtaUpper) && (c1.clusterEta() == crystalEtaUpper) ) {
        if ( (c2.clusterEnergy() > 0) && (c2.towerEtaInCard() == towerEtaLower) && (c2.clusterEta() == crystalEtaLower ) ) {
            // std::cout << "  Found two clusters along eta border (not checked if near in phi): "
            //           << " upper has energy "     << c1.clusterEnergy()
            //           << " and lower has energy " << c2.clusterEnergy() << std::endl;
	          ap_uint<5> phi1 = c1.towerPhi()*5 + c1.clusterPhi();
            ap_uint<5> phi2 = c2.towerPhi()*5 + c2.clusterPhi();
            ap_uint<5> dPhi; 
	          dPhi = (phi1 > phi2)?(phi1-phi2):(phi2-phi1);
	  
            if( dPhi < 2 ) {
              ap_uint<15> totalEnergy = c1.clusterEnergy() + c2.clusterEnergy();
              ap_uint<15> totalEt2x5 = c1.uint_et2x5() + c2.uint_et2x5();
              ap_uint<15> totalEt5x5 = c1.uint_et5x5() + c2.uint_et5x5();
              // std::cout << "  Total energy is " << totalEnergy/8.0 << std::endl;

                // Initialize a cluster with the larger cluster's position and total energy
              if (c1.clusterEnergy() > c2.clusterEnergy()){
                c1 = Cluster(totalEnergy, c1.towerEta(), c1.towerPhi(), c1.clusterEta(), c1.clusterPhi(), c1.satur(), totalEt5x5, totalEt2x5, c1.getBrems(), c1.getIsSS(), c1.getIsLooseTkss(), c1.region());
                c2 = Cluster(0, c2.towerEta(), c2.towerPhi(), c2.clusterEta(), c2.clusterPhi(), 0, 0, 0, 0, false, false, c2.region());
                cluster_list[i] = c1;
                cluster_list[j] = c2;
                // std::cout << "Stitched:  cluster_list[i] and [j] have pT, eta, phi : " << cluster_list[i].getPt() << "(" << cluster_list[i].realEta(cc) << "," << cluster_list[i].realPhi(cc) << ")" 
                //     << ", " << cluster_list[j].getPt() << "(" << cluster_list[j].realEta(cc) << "," << cluster_list[j].realPhi(cc) << ")" <<  std::endl;
              }
              else {
                // Analogous to above portion
                c1 = Cluster(0, c1.towerEta(), c1.towerPhi(), c1.clusterEta(), c1.clusterPhi(), 0, 0, 0, 0, false, false, c1.region());
                c2 = Cluster(totalEnergy, c2.towerEta(), c2.towerPhi(), c2.clusterEta(), c2.clusterPhi(), c2.satur(), totalEt5x5, totalEt2x5, c2.getBrems(), c2.getIsSS(), c2.getIsLooseTkss(), c2.region());
                cluster_list[i] = c1;
                cluster_list[j] = c2;
                // std::cout << "Stitched:  cluster_list[i] and [j] have pT, eta, phi : " << cluster_list[i].getPt() << "(" << cluster_list[i].realEta(cc) << "," << cluster_list[i].realPhi(cc) << ")" 
                //     << ", " << cluster_list[j].getPt() << "(" << cluster_list[j].realEta(cc) << "," << cluster_list[j].realPhi(cc) << ")" <<  std::endl;
              }
	          }	  
          }
        }  
	  }	  
  }
      
  return 1;
}


//--------------------------------------------------------//  

#endif
