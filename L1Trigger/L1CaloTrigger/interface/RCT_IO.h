//------------------------------------
// IO Logic for Phase2L1CaloL1RCTEmulator.cc.
// Interface between CMSSW-provided trigger primitives and the input links
// that the firmware code takes as inputs.
//------------------------------------
#ifndef L1Trigger_L1CaloTrigger_RCT_IO
#define L1Trigger_L1CaloTrigger_RCT_IO

#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

#include <bitset>

namespace p2rctIO {

static constexpr int N_BITS_CRYSTAL = 16;
static constexpr int N_BITS_ENERGY = 10;
static constexpr int N_BITS_TIMING = 5;
static constexpr int N_BITS_SPIKE = 1;
static constexpr int N_BITS_TOWER_HCAL = 16;
static constexpr int N_BITS_HCAL_FBITS = 6;
static constexpr float LSB_ENERGY = 0.5;

static constexpr int N_TOWERS_PHI = 6;
static constexpr int N_TOWERS_ETA = 17;

static constexpr int N_CARDS = 24;

static constexpr int CRYSTALS_IN_TOWER_ETA = 5;
static constexpr int CRYSTALS_IN_TOWER_PHI = 5;

static constexpr int TOWERS_IN_REGION_ETA = 8;
static constexpr int TOWERS_IN_REGION_PHI = 4;
static constexpr int REGIONS_IN_CARD_ETA = 2;
static constexpr int REGIONS_IN_CARD_PHI = 2;
static constexpr int N_REGIONS_PHI = 18;
static constexpr int N_REGIONS_ETA = 4;

static constexpr float ECAL_ETA_RANGE = 1.4841;
static constexpr float HCAL_ETA_RANGE = 1.3968;

//////////////////////////////////////////////////////////////////////////
  // RCT: indexing helper functions
  //////////////////////////////////////////////////////////////////////////

  // Assert that the card index is within bounds. (Valid cc: 0 to N_CARDS, since there are N_CARDS RCT cards)
  inline bool isValidCard(int cc) { return ((cc > -1) && (cc < N_CARDS)); }

  // RCT Cards: need to know their min/max crystal boundaries.

  // For a card (ranging from 0 to N_CARDS, since there are N_CARDS cards), return the iEta of the crystal with max iEta.
  // This represents the card boundaries in eta (identical to getEtaMax_card in the original emulator)
  inline int getCard_iEtaMax(int cc) {
    assert(isValidCard(cc));

    int etamax = 0;
    if (cc % 2 == 0)                                            // Even card: negative eta
      etamax = (N_TOWERS_ETA * CRYSTALS_IN_TOWER_ETA - 1);  // First eta half. 5 crystals in eta in 1 tower.
    else                                                        // Odd card: positive eta
      etamax = (2 * N_TOWERS_ETA * CRYSTALS_IN_TOWER_ETA - 1);
    return etamax;
  }

  // Same as above but for minimum iEta.
  inline int getCard_iEtaMin(int cc) {
    int etamin = 0;
    if (cc % 2 == 0)  // Even card: negative eta
      etamin = (0);
    else  // Odd card: positive eta
      etamin = (N_TOWERS_ETA * CRYSTALS_IN_TOWER_ETA);
    return etamin;
  }

  // Same as above but for maximum iPhi.
  inline int getCard_iPhiMax(int cc) {
    int phimax = ((cc / 2) + 1) * N_TOWERS_PHI * CRYSTALS_IN_TOWER_PHI - 1;
    return phimax;
  }

  // Same as above but for minimum iPhi.
  inline int getCard_iPhiMin(int cc) {
    int phimin = (cc / 2) * N_TOWERS_PHI * CRYSTALS_IN_TOWER_PHI;
    return phimin;
  }

  // Given the RCT card number (0 through N_CARDS-1), get the global crystal iEta of the "bottom left" corner
  inline int getCard_refCrystal_iEta(int cc) {
    if ((cc % 2) == 1) {  // if cc is odd (positive eta)
      return (N_TOWERS_ETA * CRYSTALS_IN_TOWER_ETA);
    } else {  // if cc is even (negative eta) the bottom left corner is further in eta
      return (N_TOWERS_ETA * CRYSTALS_IN_TOWER_ETA - 1);
    }
  }

  // Given the RCT card number (0 through N_CARDS-1), get the global crystal iPhi of the "bottom left" corner (0- 71*5)
  inline int getCard_refCrystal_iPhi(int cc) {
    if ((cc % 2) == 1) {
      // if cc is odd: positive eta
      return int(cc / 2) * N_TOWERS_PHI * CRYSTALS_IN_TOWER_PHI;
    } else {
      // if cc is even, the bottom left corner is further in phi, hence the +1 and -1
      return ((int(cc / 2)+1) * N_TOWERS_PHI * CRYSTALS_IN_TOWER_PHI) - 1;
    }
  }

  // Given the RCT card number (0 through N_CARDS-1), get the global HCAL tower iEta of the "bottom left" corner
  inline int getCard_refHCALTower_iEta(int cc) {
    if ((cc % 2) == 1) {  // if cc is odd (positive eta)
      return (TOWERS_IN_REGION_ETA * REGIONS_IN_CARD_ETA); // towers in BCP region * BCP regions in RCT card = #towers in RCT card
    } else {  // if cc is even (negative eta) the bottom left corner is further in eta
      return (TOWERS_IN_REGION_ETA * REGIONS_IN_CARD_ETA - 1);
    }
  }

  // Given the RCT card number (0 through N_CARDS-1), get the global HCAL tower iPhi of the "bottom left" corner
  inline int getCard_refHCALTower_iPhi(int cc) {
    if ((cc % 2) == 1) {
      // if cc is odd: positive eta
      return int(cc / 2) * N_TOWERS_PHI;
    } else {
      // if cc is even, the bottom left corner is further in phi, hence the +1 and -1
      return ((int(cc / 2)+1) * N_TOWERS_PHI) - 1;
    }
  }

class linkECAL {
    private:
    ap_uint<576> data;

    public:
    // constructor
    linkECAL() {
        data = (ap_uint<576>)0;
    }

    inline ap_uint<576> Data(void) {return data;}

    inline void setCrystal(float energy, float timing, int spike, int iEta, int iPhi) {
        // Clear the 16 bits for this crystal
        clearCrystal(iEta, iPhi);

        // Build the 16 bits for this crystal
        ap_uint<16> crystalData;
        crystalData = 
        ((ap_uint<16>)(energy/LSB_ENERGY) & 0x3FF) |
        (((ap_uint<16>)(timing)<<N_BITS_ENERGY) & 0x1F) |
        (((ap_uint<16>)(spike)<<(N_BITS_ENERGY+N_BITS_TIMING)) & 0x1) ;

        // Put this crystal's 16 bits into the 576 bit data
        int startId = iEta*CRYSTALS_IN_TOWER_PHI+iPhi;
        int start = startId * N_BITS_CRYSTAL;
        ap_uint<576> shiftedCrystalData = (ap_uint<576>)crystalData<<start;

        this->data = this->data | shiftedCrystalData;
    }

    inline void clearCrystal(int iEta, int iPhi) {
        // Find first bit corresponding to crystal at iEta, iPhi
        int startId = iEta*CRYSTALS_IN_TOWER_PHI+iPhi;
        int start = startId * N_BITS_CRYSTAL;

        ap_uint<576> spaceHolder = (ap_uint<576>)(0xFFFF)<<start; //16 1s, shifted to start at start

        this->data = this->data & ~spaceHolder; //Turn the 16 bits for this crystal to 0
    }

    inline void clearAllCrystals(void) {
        for(int i=0; i<CRYSTALS_IN_TOWER_ETA; i++) {
            for(int j=0; j<CRYSTALS_IN_TOWER_PHI; j++) {
                clearCrystal(i,j);
            }
        }
    }
};

class linkHCAL {
    private:
    ap_uint<576> data;

    public:
    // constructor
    linkHCAL() {
        data = (ap_uint<576>)0;
    }

    inline ap_uint<576> Data(void) {return data;}

    inline void setTower(float energy, int features, int iEta, int iPhi) {
        // Clear the 16 bits for this tower
        clearTower(iEta, iPhi);

        // Build the 16 bits for this tower
        ap_uint<16> towerData;
        towerData = 
        ((ap_uint<16>)(energy/LSB_ENERGY) & 0x3FF) |
        (((ap_uint<16>)(features) & 0x3F)<<N_BITS_ENERGY) ;

        // Put this tower's 16 bits into the 576 bit data
        int startId = iEta*TOWERS_IN_REGION_PHI+iPhi;
        int start = startId * N_BITS_TOWER_HCAL;
        ap_uint<576> shiftedTowerData = (ap_uint<576>)towerData<<start;

        this->data = this->data | shiftedTowerData;
    }

    inline void clearTower(int iEta, int iPhi) {
        // Find first bit corresponding to crystal at iEta, iPhi
        int startId = iEta*TOWERS_IN_REGION_PHI+iPhi;
        int start = startId * N_BITS_TOWER_HCAL;

        ap_uint<576> spaceHolder = (ap_uint<576>)(0xFFFF)<<start; //16 1s, shifted to start at start

        this->data = this->data & ~spaceHolder; //Turn the 16 bits for this crystal to 0
    }

    inline void clearAllTowers(void) {
        for(int i=0; i<TOWERS_IN_REGION_ETA; i++) {
            for(int j=0; j<TOWERS_IN_REGION_PHI; j++) {
                clearTower(i,j);
            }
        }
    }
};

class RCTcardECAL {
    private:
    linkECAL links[N_TOWERS_ETA][N_TOWERS_PHI];

    public:
    //constructor
    RCTcardECAL() {
        for(int iEta=0; iEta<N_TOWERS_ETA; iEta++) {
            for(int iPhi=0; iPhi<N_TOWERS_PHI; iPhi++) {
                links[iEta][iPhi] = linkECAL();
            }
        }
    }

    inline void addHit(float energy, float timing, int spike, int iEtaCrystalCard, int iPhiCrystalCard) {
        // Find which link(=tower) this hit is in
        int iEtaTowerCard = iEtaCrystalCard/CRYSTALS_IN_TOWER_ETA;
        int iPhiTowerCard = iPhiCrystalCard/CRYSTALS_IN_TOWER_PHI;

        // Find which crystal within that tower it is
        int iEtaCrystalTower = iEtaCrystalCard%CRYSTALS_IN_TOWER_ETA;
        int iPhiCrystalTower = iPhiCrystalCard%CRYSTALS_IN_TOWER_PHI;

        // Update that crystal within that link
        links[iEtaTowerCard][iPhiTowerCard].setCrystal(energy, timing, spike, iEtaCrystalTower, iPhiCrystalTower);
        // if (energy > 0.0) {
          // std::cout << "Link output: " << (bitset<576>)links[iEtaTowerCard][iPhiTowerCard].Data() << std::endl;
        // }
    }

    inline linkECAL getLink(int iEtaTowerCard, int iPhiTowerCard) {
      // Return the link at iEtaTowerCard, iPhiTowerCard. int accessible via .Data() method
      return links[iEtaTowerCard][iPhiTowerCard];
    }

    inline std::vector<ap_uint<576>> getAx6(int A, int iEtaOffset) {
      // Return A*6 words as the input links for an Ax6 region, with an offset of rows in iEta
      std::vector<ap_uint<576>> out;

      for(int iEta=0; iEta<A; iEta++) {
        for(int iPhi=0; iPhi<6; iPhi++) {
          out.push_back(getLink(iEta+iEtaOffset, iPhi).Data());
        }
      }
      return out;
    }
};

class BCPcardsHCAL {
    private:
    linkHCAL links[N_REGIONS_ETA][N_REGIONS_PHI];

    public:
    //constructor
    BCPcardsHCAL() {
        for(int iEta=0; iEta<N_REGIONS_ETA; iEta++) {
            for(int iPhi=0; iPhi<N_REGIONS_PHI; iPhi++) {
                links[iEta][iPhi] = linkHCAL();
            }
        }
    }

    inline void addHit(float energy, int features, int BCPiEta, int BCPiPhi, int iEtaTowerCard, int iPhiTowerCard) {
        // Update that tower within that link
        links[BCPiEta][BCPiPhi].setTower(energy, features, iEtaTowerCard, iPhiTowerCard);
    }

    inline linkHCAL getLink(int BCPiEta, int BCPiPhi) {
      // Return the link at BCPiEta, BCPiPhi. int accessible via .Data() method
      return links[BCPiEta][BCPiPhi];
    }
};

/* 
* Represents one input HCAL or ECAL hit.
*/
class SimpleCaloHit {
    private:
    float et_ = 0.;
    GlobalVector position_;  // As opposed to GlobalPoint, so we can add them (for weighted average)
    int featurebits_;

    public:
    // tool functions
    inline void setEt(float et) { et_ = et; };
    inline void setPosition(const GlobalVector& pos) { position_ = pos; };
    inline void setFeatureBits(int bits) { featurebits_ = bits; };

    inline float et() const { return et_; };
    inline const GlobalVector& position() const { return position_; };
    inline int featurebits() const { return featurebits_; };

    /* 
       * Get crystal's iEta from real eta. (identical to getCrystal_etaID in L1EGammaCrystalsEmulatorProducer.cc)
       * This "global" iEta ranges from 0 to (34*5)-1 since there are 34 towers in eta in the full detector, 
       * each with five crystals in eta.
       */
    int crystaliEta(void) const {
      float size_cell = ECAL_ETA_RANGE / (CRYSTALS_IN_TOWER_ETA * N_TOWERS_ETA);
      int iEta = int((position().eta() + ECAL_ETA_RANGE) / size_cell);
      return iEta;
    }

    /* 
       * Get crystal's iPhi from real phi. (identical to getCrystal_phiID in L1EGammaCrystalsEmulatorProducer.cc)
       * This "global" iPhi ranges from 0 to (71*5) since there are 72 towers in phi in the full detector, each with five crystals in phi.
       */
    int crystaliPhi(void) const {
      float phi = position().phi();
      float size_cell = 2 * M_PI / (CRYSTALS_IN_TOWER_PHI * N_TOWERS_PHI * N_CARDS / 2);
      int iPhi = int((phi + M_PI) / size_cell);
      return iPhi;
    }

    /* 
       * Get HCAL tower's iEta from real eta.
       * This "global" iEta ranges from 0 to 31 since there are 32 towers in eta in the full detector
       */
    int HCALtoweriEta(void) const {
      float size_cell = HCAL_ETA_RANGE / (TOWERS_IN_REGION_ETA * REGIONS_IN_CARD_ETA);
      int iEta = int((position().eta() + HCAL_ETA_RANGE) / size_cell);
      return iEta;
    }

    /* 
       * Get tower's iPhi from real phi.
       * This "global" iPhi ranges from 0 to 6*12-1 since there are 12 cards spanning phi, each with 6 towers in phi.
       */
    int HCALtoweriPhi(void) const {
      float phi = position().phi();
      float size_cell = 2 * M_PI / (N_TOWERS_PHI * N_CARDS / 2);
      int iPhi = int((phi + M_PI) / size_cell);
      return iPhi;
    }

    /*
       * Check if it falls within the boundary of a card.
       */
    bool isInCard(int cc) const {
      // Written using ECAL crystal coordinates, but works equally for HCAL
      return (crystaliPhi() <= getCard_iPhiMax(cc) && crystaliPhi() >= getCard_iPhiMin(cc) &&
              crystaliEta() <= getCard_iEtaMax(cc) && crystaliEta() >= getCard_iEtaMin(cc));
    };

    /*
      * For a crystal with real eta, and falling in card cc, get its local iEta 
      * relative to the bottom left corner of the card (possible local iEta ranges from 0 to N_TOWERS_ETA * CRYSTALS_IN_TOWER_ETA,
      * since in one card, there are N_TOWERS_ETA towers in eta, each with CRYSTALS_IN_TOWER_ETA crystals in eta.
      */
    int crystalLocaliEta(int cc) const { return abs(getCard_refCrystal_iEta(cc) - crystaliEta()); }

    /*
      * Same as above, but for iPhi (possible local iPhi ranges from 0 to (N_TOWERS_PHI*CRYSTALS_IN_TOWER_PHI), since in one card,
      * there are N_TOWERS_PHI towers in phi, each with CRYSTALS_IN_TOWER_PHI crystals in phi.
      */
    int crystalLocaliPhi(int cc) const { return abs(getCard_refCrystal_iPhi(cc) - crystaliPhi()); }

    /*
      * For a tower with real eta, and falling in card cc, get its local iEta 
      * relative to the bottom left corner of the card (possible local iEta ranges from 0 to
      * TOWERS_IN_REGION_ETA * REGIONS_IN_CARD_ETA.
      */
    int HCALtowerLocaliEta(int cc) const { return abs(getCard_refHCALTower_iEta(cc) - HCALtoweriEta()); }

    // Same as above, but for iPhi
    int HCALtowerLocaliPhi(int cc) const { return abs(getCard_refHCALTower_iPhi(cc) - HCALtoweriPhi()); }

    // Find the iEta for the HCAL BCP card this is in (0 through 3)
    int BCPcardiEta(void) const {
      int hit_iEta = HCALtoweriEta();
      int out = hit_iEta / TOWERS_IN_REGION_ETA;
      return out;
    }

    // Same as above, but for iPhi (0 through 17)
    int BCPcardiPhi(void) const {
      int hit_iPhi = HCALtoweriPhi();
      int out = hit_iPhi / TOWERS_IN_REGION_PHI;
      return out;
    }

    // Find the iEta of this hit within its HCAL BCP card (0 through 7)
    int iEtaInBCPcard(void) const {
      int hit_iEta = HCALtoweriEta();
      int BCP_iEta = BCPcardiEta();
      int out;
      if (BCP_iEta >= N_REGIONS_ETA/2) {
        out = hit_iEta - BCP_iEta*TOWERS_IN_REGION_ETA; //positive eta
      }
      else {
        out = (BCP_iEta+1)*TOWERS_IN_REGION_ETA - 1 - hit_iEta;
      }
      return out;
    }

    // Same as above, but for iPhi (0 through 3)
    int iPhiInBCPcard(void) const {
      int hit_iPhi = HCALtoweriPhi();
      int BCP_iPhi = BCPcardiPhi();
      int BCP_iEta = BCPcardiEta();
      int out;
      if (BCP_iEta >= N_REGIONS_ETA/2) {
        out = hit_iPhi - BCP_iPhi*TOWERS_IN_REGION_PHI; //positive eta
      }
      else {
        out = (BCP_iPhi+1)*TOWERS_IN_REGION_PHI - 1 - hit_iPhi;
      }
      return out;
    }
};

inline void printIP3OutputClusters(ap_uint<576> Data) {
  ap_uint<10> seed_pt;
  ap_uint<12> pt;
  ap_uint<7> eta;
  ap_uint<5> phi;
  ap_uint<10> et5x5;
  ap_uint<7> wps;
  ap_uint<5> timing;
  ap_uint<1> spike;
  ap_uint<1> satur;
  ap_uint<2> brems;
  ap_uint<4> spare;

  for(int i=0; i<9; i++){
    int start = i*64;
    std::cout << "Cluster " << i << std::endl;

    seed_pt = Data.range(start+9,start);
    std::cout << "seed_pt: " << seed_pt*0.5 << std::endl;
    pt = Data.range(start+21,start+10);
    std::cout << "pt: " << pt*0.5 << std::endl;
    eta = Data.range(start+28,start+22);
    std::cout << "eta: " << eta << std::endl;
    phi = Data.range(start+33,start+29);
    std::cout << "phi: " << phi << std::endl;
    et5x5 = Data.range(start+43,start+34);
    std::cout << "et5x5: " << et5x5*2.0 << std::endl;
    wps = Data.range(start+50,start+44);
    std::cout << "wps: " << wps << std::endl;
    timing = Data.range(start+55,start+51);
    std::cout << "timing: " << timing << std::endl;
    spike = Data.range(start+56,start+56);
    std::cout << "spike: " << spike << std::endl;
    satur = Data.range(start+57,start+57);
    std::cout << "satur: " << satur << std::endl;
    brems = Data.range(start+59,start+58);
    std::cout << "brems: " << brems << std::endl;
    spare = Data.range(start+63,start+60);
    std::cout << "spare: " << spare << std::endl;
  }
}

inline void printIP3OutputTowers(ap_uint<576> Data) {
  ap_uint<10> et;
  ap_uint<4> hoe;
  ap_uint<2> fb;

  std::cout << "Lower iPhi in this link:" << std::endl;
  for(int i=0; i<17; i++) {
    int start = i*16;
    if(Data.range(start+9,start) > 0){
    std::cout << "Tower iEta " << i << std::endl;

    et = Data.range(start+9,start);
    std::cout << "et: " << et*0.5 << std::endl;
    hoe = Data.range(start+13,start+10);
    std::cout << "hoe: " << hoe << std::endl;
    fb = Data.range(start+15,start+14);
    std::cout << "fb: " << fb << std::endl;
    }
  }

  std::cout << "Higher iPhi in this link:" << std::endl;
  for(int i=0; i<17; i++) {
    int start = i*16 + 272;
    if(Data.range(start+9,start) > 0){
    std::cout << "Tower iEta " << i << std::endl;

    et = Data.range(start+9,start);
    std::cout << "et: " << et*0.5 << std::endl;
    hoe = Data.range(start+13,start+10);
    std::cout << "hoe: " << hoe << std::endl;
    fb = Data.range(start+15,start+14);
    std::cout << "fb: " << fb << std::endl;
    }
  }
}

} // namespace p2rctIO

#endif
