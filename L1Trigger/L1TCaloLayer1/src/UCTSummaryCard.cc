#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

#include "L1Trigger/L1TCaloLayer1/interface/UCTSummaryCard.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTObject.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTLayer1.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTCrate.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTCard.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTGeometry.hh"
#include "L1Trigger/L1TCaloLayer1/interface/UCTLogging.hh"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

using namespace l1tcalo;

UCTSummaryCard::UCTSummaryCard(const std::vector<std::vector<std::vector<uint32_t> > >* l,
                               uint32_t jetSeedIn,
                               uint32_t tauSeedIn,
                               double tauIsolationFactorIn,
                               uint32_t eGammaSeedIn,
                               double eGammaIsolationFactorIn)
    : pumLUT(l),
      jetSeed(jetSeedIn),
      tauSeed(tauSeedIn),
      tauIsolationFactor(tauIsolationFactorIn),
      eGammaSeed(tauSeedIn),
      eGammaIsolationFactor(tauIsolationFactorIn) {
  UCTGeometry g;
  sinPhi[0] = 0;
  cosPhi[0] = 1;
  for (int iPhi = 1; iPhi <= 72; iPhi++) {
    sinPhi[iPhi] = sin(g.getUCTTowerPhi(iPhi));
    cosPhi[iPhi] = cos(g.getUCTTowerPhi(iPhi));
  }
}

UCTSummaryCard::~UCTSummaryCard() {
  for (uint32_t i = 0; i < regions.size(); i++) {
    if (regions[i] != nullptr)
      delete regions[i];
  }
}

bool UCTSummaryCard::process() {
  clearEvent();
  UCTGeometry g;
  uint32_t etValue = 0;
  uint32_t htValue = 0;
  int sumEx = 0;
  int sumEy = 0;
  int sumHx = 0;
  int sumHy = 0;
  int etaMin = -NRegionsInCard;
  int etaMax = NRegionsInCard;
  // Determine pumLevel using only the central regions
  pumLevel = 0;
  for (int iEta = etaMin; iEta <= etaMax; iEta++) {  // Ranging between -7 to 7
    if (iEta == 0)
      continue;                                                 // Note that eta == 0 is illegal
    for (uint32_t iPhi = 0; iPhi < MaxUCTRegionsPhi; iPhi++) {  // Note that phi ranges from 0 to 17
      const UCTRegion* uctRegion = getRegion(iEta, iPhi);
      uint32_t et = uctRegion->et();
      if (et > 0)
        pumLevel++;
    }
  }
  uint32_t totalRegionCount = 2 * etaMax * MaxUCTRegionsPhi;
  uint32_t pumBinSize = totalRegionCount / pumLUT->size();
  pumBin = floor(pumLevel / pumBinSize);
  if (pumBin >= pumLUT->size())
    pumBin = pumLUT->size() - 1;  // Max index
  // We walk the eta-phi plane looping over all regions.
  // to make global objects like TotalET, HT, MET, MHT
  // subtracting pileup along the way
  // For compact objects we use processRegion(regionIndex)
  uint32_t pileup = 0;
  uint32_t et = 0;
  for (uint32_t side = 0; side < 2; side++) {
    bool negativeSide = true;
    if (side == 1)
      negativeSide = false;
    for (uint32_t crate = 0; crate < g.getNCrates(); crate++) {
      for (uint32_t card = 0; card < g.getNCards(); card++) {
        for (uint32_t region = 0; region < NRegionsInCard; region++) {
          int iEta = g.getUCTRegionEtaIndex(negativeSide, region);
          uint32_t iPhi = g.getUCTRegionPhiIndex(crate, card);
          UCTRegionIndex regionIndex(iEta, iPhi);
          processRegion(regionIndex);
          //const UCTRegion* uctRegion = getRegion(iEta, iPhi);
          const UCTRegion* uctRegion = new UCTRegion(crate, card, negativeSide, region, 1); 

          //regions.push_back(new UCTRegion(crate, card, true, rgn, fwVersion));
          //regions.push_back(new UCTRegion(crate, card, false, rgn, fwVersion));
          if (uctRegion == nullptr) {
            continue;
          }
          et = uctRegion->et();
          uint32_t pileupInRegion = (*pumLUT)[pumBin][side][region];
          pileup += pileupInRegion;
          if (pileupInRegion < et)
            et -= pileupInRegion;
          else
            et = 0;
          UCTTowerIndex hitTower = g.getUCTTowerIndexFromL1CaloRegion(regionIndex, uctRegion->rawData());
          int hitCaloPhi = hitTower.second;
          sumEx += ((int)(((double)et) * cosPhi[hitCaloPhi]));
          sumEy += ((int)(((double)et) * sinPhi[hitCaloPhi]));
          etValue += et;
          if (et > 10) {
            sumHx += ((int)(((double)et) * cosPhi[hitCaloPhi]));
            sumHy += ((int)(((double)et) * sinPhi[hitCaloPhi]));
            htValue += et;
          }
          delete uctRegion;
        }
      }
    }
  }
  uint32_t metSquare = sumEx * sumEx + sumEy * sumEy;
  uint32_t metValue = sqrt((double)metSquare);
  double metPhi = (atan2(sumEy, sumEx) * 180. / 3.1415927) + 180.;
  int metIPhi = (int)(72. * (metPhi / 360.));
  uint32_t mhtSquare = sumHx * sumHx + sumHy * sumHy;
  uint32_t mhtValue = sqrt((double)mhtSquare);
  double mhtPhi = (atan2(sumHy, sumHx) * 180. / 3.1415927) + 180.;
  int mhtIPhi = (int)(72. * (mhtPhi / 360.));

  ET = new UCTObject(UCTObject::ET, etValue, 0, metIPhi, pileup, 0, 0);
  HT = new UCTObject(UCTObject::HT, htValue, 0, mhtIPhi, pileup, 0, 0);
  MET = new UCTObject(UCTObject::MET, metValue, 0, metIPhi, pileup, 0, 0);
  MHT = new UCTObject(UCTObject::MHT, mhtValue, 0, mhtIPhi, pileup, 0, 0);

  // Then sort the candidates for output usage
  emObjs.sort();
  isoEMObjs.sort();
  tauObjs.sort();
  isoTauObjs.sort();
  centralJetObjs.sort();
  forwardJetObjs.sort();
  boostedJetObjs.sort();
  // Cool we never fail :)
  return true;
}

bool UCTSummaryCard::processRegion(UCTRegionIndex center) {
  // We process the region looking at nearest neighbor data
  // We should never need beyond nearest-neighbor for most
  // objects - eGamma, tau or jet

  std::vector<uint32_t> boostedJetTowers;
  boostedJetTowers.clear();
  for (size_t iPhi = 0; iPhi < 12; iPhi++) {
    for (size_t iEta = 0; iEta < 12; iEta++) {
      boostedJetTowers.push_back(0);
    }
  }

  std::vector<uint32_t> boostedJetRegionET, boostedJetRegionTauVeto, boostedJetRegionEGammaVeto;
  boostedJetRegionET.clear();
  boostedJetRegionTauVeto.clear();
  boostedJetRegionEGammaVeto.clear();
  for (size_t i = 0; i < 9; i++) {
    boostedJetRegionET.push_back(0);
    boostedJetRegionTauVeto.push_back(0);
    boostedJetRegionEGammaVeto.push_back(0);
  }

  UCTGeometryExtended g;

  const UCTRegion* cRegion = getRegion(center.first, center.second);
  if (cRegion == nullptr) {
    return false;
  }

  // Central jets do not use HF - note border exclusion
  bool centralRegion = false;
  if (cRegion->getRegion() < (NRegionsInCard - 1))
    centralRegion = true;

  uint32_t centralET = cRegion->et();
  uint32_t centralPU = std::min(centralET, (*pumLUT)[pumBin][0][cRegion->getRegion()]);
  bitset<4> cActiveTowerEta = 0;
  bitset<4> cActiveTowerPhi = 0;
  bool cActiveTowers[4][4];
  //bitset<4> cActiveTowerEta = cRegion->activeTowerEta();
  //bitset<4> cActiveTowerPhi = cRegion->activeTowerPhi();
  //std::cout<<"cActiveTowerEta: "<<cActiveTowerEta<<"\t"<<"cActiveTowerPhi: "<<cActiveTowerPhi<<std::endl;
  if (!cRegion->isNegativeEta())
    centralPU = std::min(centralET, (*pumLUT)[pumBin][1][cRegion->getRegion()]);
  centralET -= centralPU;
  UCTTowerIndex centralHitTower = g.getUCTTowerIndexFromL1CaloRegion(center, cRegion->rawData());
  int nTauLike = 0;
  int nEGammaLike = 0;
  bool centralIsTauLike = cRegion->isTauLike();
  if (centralIsTauLike)
    nTauLike++;
  bool centralIsEGammaLike = cRegion->isEGammaLike();
  if (centralIsEGammaLike)
    nEGammaLike++;
  int hitCaloEta = centralHitTower.first;
  int hitCaloPhi = centralHitTower.second;
  boostedJetRegionET[4] = centralET;
  boostedJetRegionTauVeto[4] = centralIsTauLike;
  boostedJetRegionEGammaVeto[4] = centralIsEGammaLike;
  //get the towers from central
  int etaOffset = 4;
  int phiOffset = 4;

  //vector<UCTTower*> ctowers;
  //ctowers.clear();
  //for (uint32_t iEta = 0; iEta < 4; iEta++) {
  //  for (uint32_t iPhi = 0; iPhi < 4; iPhi++) {
  //    ctowers.push_back(new UCTTower(cRegion->getCrate(), cRegion->getCard(), cRegion->isNegativeEta(), cRegion->getRegion(), iEta, iPhi, 3));
  //  }
  //}

  //vector<UCTTower*> ctowers = cRegion->regionTowers();
  //for(size_t iPhi = 0; iPhi < 4; iPhi++){
  //  for(size_t iEta = 0; iEta < 4; iEta++) {
  //    if(iEta*4+iPhi < ctowers.size()) {
  //      if(center.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = ctowers[iEta*4+iPhi]->et();
  //      else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = ctowers[iEta*4+iPhi]->et();

  //      uint32_t towerET = ctowers[iEta*4+iPhi]->et();
  //      std::cout<<"towerET: "<<towerET<<std::endl;
  //      if(towerET > centralET*0.015875) {
  //        cActiveTowers[iEta][iPhi] = true;
  //      }
  //      else
  //        cActiveTowers[iEta][iPhi] = false;
  //    }
      //std::cout<<cActiveTowers[iEta][iPhi]<<std::endl;
  //  }
  //}
  //for(uint32_t iEta = 0; iEta < 4; iEta++) {
  //  bool activeStrip = false;
  //  for(uint32_t iPhi = 0; iPhi < 4; iPhi++) {
  //    if(cActiveTowers[iEta][iPhi]) activeStrip = true;
  //  }
  //  if(activeStrip) cActiveTowerEta |= (0x1 << iEta);
  //}
  //for(uint32_t iPhi = 0; iPhi < 4; iPhi++) {
  //  bool activeStrip = false;
  //  for(uint32_t iEta = 0; iEta < 4; iEta++) {
  //    if(cActiveTowers[iEta][iPhi]) activeStrip = true;
  //  }
  //  if(activeStrip) cActiveTowerPhi |= (0x1 << iPhi);
  //}
  //std::cout<<"cActiveTowerEta: "<<cActiveTowerEta<<"\t"<<"cActiveTowerPhi: "<<cActiveTowerPhi<<std::endl;

  UCTRegionIndex northIndex = g.getUCTRegionNorth(center);
  const UCTRegion* northRegion = getRegion(northIndex.first, northIndex.second);
  uint32_t northET = 0;
  uint32_t northPU = 0;
  UCTTowerIndex northHitTower;
  bool northIsTauLike = false;
  bool northIsEGammaLike = false;
  bitset<4> nActiveTowerEta = 0;
  bitset<4> nActiveTowerPhi = 0;
  if (northRegion != nullptr) {
    //nActiveTowerEta = northRegion->activeTowerEta();
    //nActiveTowerPhi = northRegion->activeTowerPhi();
    northET = northRegion->et();
    northPU = std::min(northET, (*pumLUT)[pumBin][0][northRegion->getRegion()]);
    if (!northRegion->isNegativeEta())
      northPU = std::min(northET, (*pumLUT)[pumBin][1][northRegion->getRegion()]);
    northET -= northPU;
    northHitTower = g.getUCTTowerIndexFromL1CaloRegion(northIndex, northRegion->rawData());
    northIsTauLike = northRegion->isTauLike();
    if (northIsTauLike)
      nTauLike++;
    northIsEGammaLike = northRegion->isEGammaLike();
    if (northIsEGammaLike)
      nEGammaLike++;
    boostedJetRegionET[3] = northET;
    boostedJetRegionTauVeto[3] = northIsTauLike;
    boostedJetRegionEGammaVeto[3] = northIsEGammaLike;
    //get the towers from North
    int etaOffset = 4;
    int phiOffset = 8;
    for(uint32_t iPhi = 0; iPhi < 4; iPhi++){
      for(uint32_t iEta = 0; iEta < 4; iEta++) {
        if(iEta*4+iPhi < northRegion->regionTowers().size()) {
           if(northIndex.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = northRegion->regionTowers()[iEta*4+iPhi]->et();
           else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = northRegion->regionTowers()[iEta*4+iPhi]->et();
        }
      }
    }
  }

  UCTRegionIndex southIndex = g.getUCTRegionSouth(center);
  const UCTRegion* southRegion = getRegion(southIndex.first, southIndex.second);
  uint32_t southET = 0;
  uint32_t southPU = 0;
  UCTTowerIndex southHitTower;
  bool southIsTauLike = false;
  bool southIsEGammaLike = false;
  bitset<4> sActiveTowerEta = 0;
  bitset<4> sActiveTowerPhi = 0;
  if (southRegion != nullptr) {
    //sActiveTowerEta = southRegion->activeTowerEta();
    //sActiveTowerPhi = southRegion->activeTowerPhi();
    southET = southRegion->et();
    southPU = std::min(southET, (*pumLUT)[pumBin][0][southRegion->getRegion()]);
    if (!southRegion->isNegativeEta())
      southPU = std::min(southET, (*pumLUT)[pumBin][1][southRegion->getRegion()]);
    southET -= southPU;
    southHitTower = g.getUCTTowerIndexFromL1CaloRegion(southIndex, southRegion->rawData());
    southIsTauLike = southRegion->isTauLike();
    if (southIsTauLike)
      nTauLike++;
    southIsEGammaLike = southRegion->isEGammaLike();
    if (southIsEGammaLike)
      nEGammaLike++;
    boostedJetRegionET[5] = southET;
    boostedJetRegionTauVeto[5] = southIsTauLike;
    boostedJetRegionEGammaVeto[5] = southIsEGammaLike;
    //get the towers from south
    int etaOffset = 4;
    int phiOffset = 0;
    for(uint32_t iPhi = 0; iPhi < 4; iPhi++){
      for(uint32_t iEta = 0; iEta < 4; iEta++) {
        if(iEta*4+iPhi < southRegion->regionTowers().size()) {
          if(southIndex.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = southRegion->regionTowers()[iEta*4+iPhi]->et();
          else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = southRegion->regionTowers()[iEta*4+iPhi]->et();
        }
      }
    }
  }

  UCTRegionIndex westIndex = g.getUCTRegionWest(center);
  const UCTRegion* westRegion = getRegion(westIndex.first, westIndex.second);
  uint32_t westET = 0;
  uint32_t westPU = 0;
  UCTTowerIndex westHitTower;
  bool westIsTauLike = false;
  bool westIsEGammaLike = false;
  bitset<4> wActiveTowerEta = 0;
  bitset<4> wActiveTowerPhi = 0;
  if (westRegion != nullptr) {
    //wActiveTowerEta = westRegion->activeTowerEta();
    //wActiveTowerPhi = westRegion->activeTowerPhi();
    westET = westRegion->et();
    westPU = std::min(westET, (*pumLUT)[pumBin][0][westRegion->getRegion()]);
    if (!westRegion->isNegativeEta())
      westPU = std::min(westET, (*pumLUT)[pumBin][1][westRegion->getRegion()]);
    westET -= westPU;
    westHitTower = g.getUCTTowerIndexFromL1CaloRegion(westIndex, westRegion->rawData());
    westIsTauLike = westRegion->isTauLike();
    if (westIsTauLike)
      nTauLike++;
    westIsEGammaLike = westRegion->isEGammaLike();
    if (westIsEGammaLike)
      nEGammaLike++;
    boostedJetRegionET[7] = westET;
    boostedJetRegionTauVeto[7] = westIsEGammaLike;
    boostedJetRegionEGammaVeto[7] = westIsEGammaLike;
    //get the towers from west
    int etaOffset = 0;
    int phiOffset = 4;
    for(uint32_t iPhi = 0; iPhi < 4; iPhi++){
      for(uint32_t iEta = 0; iEta < 4; iEta++) {
        if(iEta*4+iPhi < westRegion->regionTowers().size()) {
          if(westIndex.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = westRegion->regionTowers()[iEta*4+iPhi]->et();
          else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = westRegion->regionTowers()[iEta*4+iPhi]->et();
        }
      }
    }
  }

  UCTRegionIndex eastIndex = g.getUCTRegionEast(center);
  const UCTRegion* eastRegion = getRegion(eastIndex.first, eastIndex.second);
  uint32_t eastET = 0;
  uint32_t eastPU = 0;
  UCTTowerIndex eastHitTower;
  bool eastIsTauLike = false;
  bool eastIsEGammaLike = false;
  bitset<4> eActiveTowerEta = 0;
  bitset<4> eActiveTowerPhi = 0;
  if (eastRegion != nullptr) {
    //eActiveTowerEta = eastRegion->activeTowerEta();
    //eActiveTowerPhi = eastRegion->activeTowerPhi();
    eastET = eastRegion->et();
    eastPU = std::min(eastET, (*pumLUT)[pumBin][0][eastRegion->getRegion()]);
    if (!eastRegion->isNegativeEta())
      eastPU = std::min(eastET, (*pumLUT)[pumBin][1][eastRegion->getRegion()]);
    eastET -= eastPU;
    eastHitTower = g.getUCTTowerIndexFromL1CaloRegion(eastIndex, eastRegion->rawData());
    eastIsTauLike = eastRegion->isTauLike();
    if (eastIsTauLike)
      nTauLike++;
    eastIsEGammaLike = eastRegion->isEGammaLike();
    if (eastIsEGammaLike)
      nEGammaLike++;
    boostedJetRegionET[1] = eastET;
    boostedJetRegionTauVeto[1] = eastIsTauLike;
    boostedJetRegionEGammaVeto[1] = eastIsEGammaLike;
    //get the towers from east
    int etaOffset = 8;
    int phiOffset = 4;
    for(uint32_t iPhi = 0; iPhi < 4; iPhi++){
      for(uint32_t iEta = 0; iEta < 4; iEta++) {
        if(iEta*4+iPhi < eastRegion->regionTowers().size()) {
          if(eastIndex.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = eastRegion->regionTowers()[iEta*4+iPhi]->et();
          else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = eastRegion->regionTowers()[iEta*4+iPhi]->et();
        }
      }
    }
  }

  UCTRegionIndex nwIndex = g.getUCTRegionNW(center);
  const UCTRegion* nwRegion = getRegion(nwIndex.first, nwIndex.second);
  uint32_t nwET = 0;
  uint32_t nwPU = 0;
  UCTTowerIndex nwHitTower;
  bool nwIsTauLike = false;
  bool nwIsEGammaLike = false;
  bitset<4> nwActiveTowerEta = 0;
  bitset<4> nwActiveTowerPhi = 0;
  if (nwRegion != nullptr) {
    //nwActiveTowerEta = nwRegion->activeTowerEta();
    //nwActiveTowerPhi = nwRegion->activeTowerPhi();
    nwET = nwRegion->et();
    nwPU = std::min(nwET, (*pumLUT)[pumBin][0][nwRegion->getRegion()]);
    if (!nwRegion->isNegativeEta())
      nwPU = std::min(nwET, (*pumLUT)[pumBin][1][nwRegion->getRegion()]);
    nwET -= nwPU;
    nwHitTower = g.getUCTTowerIndexFromL1CaloRegion(nwIndex, nwRegion->rawData());
    nwIsTauLike = nwRegion->isTauLike();
    if (nwIsTauLike)
      nTauLike++;
    nwIsEGammaLike = nwRegion->isEGammaLike();
    if (nwIsEGammaLike)
      nEGammaLike++;
    boostedJetRegionET[6] = nwET;
    boostedJetRegionTauVeto[6] = nwIsTauLike;
    boostedJetRegionEGammaVeto[6] = nwIsEGammaLike;
    //get the towers from north west
    int etaOffset = 0;
    int phiOffset = 8;
    for(uint32_t iPhi = 0; iPhi < 4; iPhi++){
      for(uint32_t iEta = 0; iEta < 4; iEta++) {
        if(iEta*4+iPhi < nwRegion->regionTowers().size()) {
          if(nwIndex.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = nwRegion->regionTowers()[iEta*4+iPhi]->et();
          else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = nwRegion->regionTowers()[iEta*4+iPhi]->et();
        }
      }
    }
  }

  UCTRegionIndex neIndex = g.getUCTRegionNE(center);
  const UCTRegion* neRegion = getRegion(neIndex.first, neIndex.second);
  uint32_t neET = 0;
  uint32_t nePU = 0;
  UCTTowerIndex neHitTower;
  bool neIsTauLike = false;
  bool neIsEGammaLike = false;
  bitset<4> neActiveTowerEta = 0;
  bitset<4> neActiveTowerPhi = 0;
  if (neRegion != nullptr) {
    //neActiveTowerEta = neRegion->activeTowerEta();
    //neActiveTowerPhi = neRegion->activeTowerPhi();
    neET = neRegion->et();
    nePU = std::min(neET, (*pumLUT)[pumBin][0][neRegion->getRegion()]);
    if (!neRegion->isNegativeEta())
      nePU = std::min(neET, (*pumLUT)[pumBin][1][neRegion->getRegion()]);
    neET -= nePU;
    neHitTower = g.getUCTTowerIndexFromL1CaloRegion(neIndex, neRegion->rawData());
    neIsTauLike = neRegion->isTauLike();
    if (neIsTauLike)
      nTauLike++;
    neIsEGammaLike = neRegion->isEGammaLike();
    if (neIsEGammaLike)
      nEGammaLike++;
    boostedJetRegionET[0] = neET;
    boostedJetRegionTauVeto[0] = neIsTauLike;
    boostedJetRegionEGammaVeto[0] = neIsEGammaLike;
    //get the towers from north east
    int etaOffset = 8;
    int phiOffset = 8;
    for(uint32_t iPhi = 0; iPhi < 4; iPhi++){
      for(uint32_t iEta = 0; iEta < 4; iEta++) {
        if(iEta*4+iPhi < neRegion->regionTowers().size()) {
           if(neIndex.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = neRegion->regionTowers()[iEta*4+iPhi]->et();
           else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = neRegion->regionTowers()[iEta*4+iPhi]->et();
        }
      }
    }
  }

  UCTRegionIndex swIndex = g.getUCTRegionSW(center);
  const UCTRegion* swRegion = getRegion(swIndex.first, swIndex.second);
  uint32_t swET = 0;
  uint32_t swPU = 0;
  UCTTowerIndex swHitTower;
  bool swIsTauLike = false;
  bool swIsEGammaLike = false;
  bitset<4> swActiveTowerEta = 0;
  bitset<4> swActiveTowerPhi = 0;
  if (swRegion != nullptr) {
    //swActiveTowerEta = swRegion->activeTowerEta();
    //swActiveTowerPhi = swRegion->activeTowerPhi();
    swET = swRegion->et();
    swPU = std::min(swET, (*pumLUT)[pumBin][0][swRegion->getRegion()]);
    if (!swRegion->isNegativeEta())
      swPU = std::min(swET, (*pumLUT)[pumBin][1][swRegion->getRegion()]);
    swET -= swPU;
    swHitTower = g.getUCTTowerIndexFromL1CaloRegion(swIndex, swRegion->rawData());
    swIsTauLike = swRegion->isTauLike();
    if (swIsTauLike)
      nTauLike++;
    swIsEGammaLike = swRegion->isEGammaLike();
    if (swIsEGammaLike)
      nEGammaLike++;
    boostedJetRegionET[8] = swET;
    boostedJetRegionTauVeto[8] = swIsTauLike;
    boostedJetRegionEGammaVeto[8] = swIsEGammaLike;
    //get the towers from south west
    int etaOffset = 0;
    int phiOffset = 0;
    for(uint32_t iPhi = 0; iPhi < 4; iPhi++){
      for(uint32_t iEta = 0; iEta < 4; iEta++) {
        if(iEta*4+iPhi < swRegion->regionTowers().size()) {
          if(swIndex.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = swRegion->regionTowers()[iEta*4+iPhi]->et();
          else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = swRegion->regionTowers()[iEta*4+iPhi]->et();
        }
      }
    }
  }

  UCTRegionIndex seIndex = g.getUCTRegionSE(center);
  const UCTRegion* seRegion = getRegion(seIndex.first, seIndex.second);
  uint32_t seET = 0;
  uint32_t sePU = 0;
  UCTTowerIndex seHitTower;
  bool seIsTauLike = false;
  bool seIsEGammaLike = false;
  bitset<4> seActiveTowerEta = 0;
  bitset<4> seActiveTowerPhi = 0;
  if (seRegion != nullptr) {
    //seActiveTowerEta = seRegion->activeTowerEta();
    //seActiveTowerPhi = seRegion->activeTowerPhi();
    seET = seRegion->et();
    sePU = std::min(seET, (*pumLUT)[pumBin][0][seRegion->getRegion()]);
    if (!seRegion->isNegativeEta())
      sePU = std::min(seET, (*pumLUT)[pumBin][1][seRegion->getRegion()]);
    seET -= sePU;
    seHitTower = g.getUCTTowerIndexFromL1CaloRegion(seIndex, seRegion->rawData());
    seIsTauLike = seRegion->isTauLike();
    if (seIsTauLike)
      nTauLike++;
    seIsEGammaLike = seRegion->isEGammaLike();
    if (seIsEGammaLike)
      nEGammaLike++;
    boostedJetRegionET[2] = seET;
    boostedJetRegionTauVeto[2] = seIsTauLike;
    boostedJetRegionEGammaVeto[2] = seIsEGammaLike;
    //get the towers from south east
    int etaOffset = 8;
    int phiOffset = 0;
    for(uint32_t iPhi = 0; iPhi < 4; iPhi++){
      for(uint32_t iEta = 0; iEta < 4; iEta++) {
        if(iEta*4+iPhi < seRegion->regionTowers().size()) {
          if(seIndex.first < 0) boostedJetTowers[(etaOffset+3-iEta)*12+(phiOffset+iPhi)] = seRegion->regionTowers()[iEta*4+iPhi]->et();
          else boostedJetTowers[(etaOffset+iEta)*12+(phiOffset+iPhi)] = seRegion->regionTowers()[iEta*4+iPhi]->et();
        }
      }
    }
  }

  uint32_t et3x3 = centralET + northET + nwET + westET + swET + southET + seET + eastET + neET;
  if (et3x3 > 0x3FF)
    et3x3 = 0x3FF;  // Peg at 10-bits

  uint32_t pu3x3 = centralPU + northPU + nwPU + westPU + swPU + southPU + sePU + eastPU + nePU;

  // Jet - a 3x3 object with center greater than a seed and all neighbors
  uint32_t jetET = et3x3;

  if (centralET >= northET && centralET >= nwET && centralET >= westET && centralET >= swET && centralET > southET &&
      centralET > seET && centralET > eastET && centralET > neET && centralET > jetSeed) {
    if (centralRegion)
      centralJetObjs.push_back(new UCTObject(UCTObject::jet, jetET, hitCaloEta, hitCaloPhi, pu3x3, 0, et3x3));
    else
      forwardJetObjs.push_back(new UCTObject(UCTObject::jet, jetET, hitCaloEta, hitCaloPhi, pu3x3, 0, et3x3));
  }

  auto boostedJet = new UCTObject(UCTObject::jet, jetET, hitCaloEta, hitCaloPhi, pu3x3, 0, et3x3);
  boostedJet->setNTaus(nTauLike);
  //boostedJet->setNEGammas(nEGammaLike);
  boostedJet->setBoostedJetRegionET(boostedJetRegionET);
  boostedJet->setBoostedJetRegionTauVeto(boostedJetRegionTauVeto);
  boostedJet->setBoostedJetRegionEGammaVeto(boostedJetRegionEGammaVeto);

  bitset<4> wEta = nwActiveTowerEta | wActiveTowerEta | swActiveTowerEta ;
  bitset<4> cEta =  nActiveTowerEta | cActiveTowerEta |  sActiveTowerEta ;
  bitset<4> eEta = neActiveTowerEta | eActiveTowerEta | seActiveTowerEta ;
  bitset<4> nPhi = nwActiveTowerPhi | nActiveTowerPhi | neActiveTowerPhi ;
  bitset<4> cPhi =  wActiveTowerPhi | cActiveTowerPhi |  eActiveTowerPhi ;
  bitset<4> sPhi = swActiveTowerPhi | sActiveTowerPhi | seActiveTowerPhi ;
  bitset<12> eta((string)(wEta.to_string() + cEta.to_string() + eEta.to_string()));
  bitset<12> phi((string)(nPhi.to_string() + cPhi.to_string() + sPhi.to_string()));
  boostedJet->setActiveTowerEta(eta);
  boostedJet->setActiveTowerPhi(phi);
  boostedJet->setBoostedJetTowers(boostedJetTowers);

  //std::cout<<"patterns: "<<boostedJet->activeTowerEta()<<"\t"<<boostedJet->activeTowerPhi()<<std::endl;
  
  //std::cout<<boostedJet->boostedJetRegionET()[0]<<"\t"<<boostedJet->boostedJetRegionET()[1]<<"\t"<<boostedJet->boostedJetRegionET()[2]<<"\t"<<boostedJet->boostedJetRegionET()[3]<<std::endl;

  boostedJetObjs.push_back(boostedJet);

  // tau Object - a single region or a 2-region sum, where the neighbor with lower ET is located using matching hit calo towers

  if (centralRegion && centralIsTauLike && centralET > tauSeed) {
    uint32_t tauET = centralET;
    uint32_t tauPU = centralPU;
    int neighborMatchCount = 0;
    //check to see if we are on the edge of the calorimeter
    if (!g.isEdgeTower(centralHitTower)) {
      //Check to see if the neighbor regions are tau like and if central ET is greater
      //if region is tau like and a neighbor AND with less energy, set it to 0.
      if (g.areNeighbors(centralHitTower, northHitTower) && northIsTauLike && centralET >= northET) {
        tauET += northET;
        tauPU += northPU;
        neighborMatchCount++;
      } else if (g.areNeighbors(centralHitTower, northHitTower) && northIsTauLike && centralET < northET) {
        tauET = 0;
      }
      if (g.areNeighbors(centralHitTower, southHitTower) && southIsTauLike && centralET > southET) {
        tauET += southET;
        tauPU += southPU;
        neighborMatchCount++;
      } else if (g.areNeighbors(centralHitTower, southHitTower) && southIsTauLike && centralET <= southET) {
        tauET = 0;
      }
      if (g.areNeighbors(centralHitTower, westHitTower) && westIsTauLike && centralET >= westET) {
        tauET += westET;
        tauPU += westPU;
        neighborMatchCount++;
      } else if (g.areNeighbors(centralHitTower, westHitTower) && westIsTauLike && centralET < westET) {
        tauET = 0;
      }
      if (g.areNeighbors(centralHitTower, eastHitTower) && eastIsTauLike && centralET > eastET) {
        tauET += eastET;
        tauPU += eastPU;
        neighborMatchCount++;
      } else if (g.areNeighbors(centralHitTower, eastHitTower) && eastIsTauLike && centralET <= eastET) {
        tauET = 0;
      }
      if (neighborMatchCount > 2) {
        tauET = 0;
      }
    }
    if (tauET != 0) {
      tauObjs.push_back(new UCTObject(UCTObject::tau, tauET, hitCaloEta, hitCaloPhi, tauPU, 0xDEADBEEF, et3x3));
      // Subtract footprint
      uint32_t isolation = 0;
      if (et3x3 > tauET)
        isolation = et3x3 - tauET;
      if (isolation < ((uint32_t)(tauIsolationFactor * (double)tauET))) {
        isoTauObjs.push_back(new UCTObject(UCTObject::isoTau, tauET, hitCaloEta, hitCaloPhi, pu3x3, isolation, et3x3));
      }
    }
  }

  // eGamma Object - This is a sad story
  // eGamma should be in 2-3 contiguous towers, but we have no bandwidth to get a second cluster from cards
  // so we use essentially the same clustering as for tau, but demand that energy is almost all in the ECAL
  // pileup subtraction is critical to not overshoot - further this should work better for isolated eGamma
  // a single region or a 2-region sum, where the neighbor with lower ET is located using matching hit calo towers

  if (centralRegion && centralIsEGammaLike && centralET > eGammaSeed) {
    uint32_t eGammaET = centralET;
    uint32_t eGammaPU = centralPU;
    int neighborMatchCount = 0;
    if (!g.isEdgeTower(centralHitTower)) {
      if (g.areNeighbors(centralHitTower, northHitTower) && northIsEGammaLike && centralET >= northET) {
        eGammaET += northET;
        eGammaPU += northPU;
        neighborMatchCount++;
      } else if (g.areNeighbors(centralHitTower, northHitTower) && northIsEGammaLike && centralET < northET) {
        eGammaET = 0;
      }
      if (g.areNeighbors(centralHitTower, southHitTower) && southIsEGammaLike && centralET > southET) {
        eGammaET += southET;
        eGammaPU += southPU;
        neighborMatchCount++;
      } else if (g.areNeighbors(centralHitTower, southHitTower) && southIsEGammaLike && centralET <= southET) {
        eGammaET = 0;
      }
      if (g.areNeighbors(centralHitTower, westHitTower) && westIsEGammaLike && centralET >= westET) {
        eGammaET += westET;
        eGammaPU += westPU;
        neighborMatchCount++;
      } else if (g.areNeighbors(centralHitTower, westHitTower) && westIsEGammaLike && centralET < westET) {
        eGammaET = 0;
      }
      if (g.areNeighbors(centralHitTower, eastHitTower) && eastIsEGammaLike && centralET > eastET) {
        eGammaET += eastET;
        eGammaPU += eastPU;
        neighborMatchCount++;
      } else if (g.areNeighbors(centralHitTower, eastHitTower) && eastIsEGammaLike && centralET <= eastET) {
        eGammaET = 0;
      }
      if (neighborMatchCount > 2) {
        eGammaET = 0;
      }
    }
    if (eGammaET != 0) {
      emObjs.push_back(new UCTObject(UCTObject::eGamma, eGammaET, hitCaloEta, hitCaloPhi, eGammaPU, 0xDEADBEEF, et3x3));
      uint32_t isolation = 0;
      if (et3x3 > eGammaET)
        isolation = et3x3 - eGammaET;
      if (isolation < ((uint32_t)(eGammaIsolationFactor * (double)eGammaET))) {
        isoEMObjs.push_back(
            new UCTObject(UCTObject::isoEGamma, eGammaET, hitCaloEta, hitCaloPhi, pu3x3, isolation, et3x3));
      }
    }
  }

  return true;
}

bool UCTSummaryCard::clearEvent() {
  emObjs.clear();
  isoEMObjs.clear();
  tauObjs.clear();
  isoTauObjs.clear();
  centralJetObjs.clear();
  forwardJetObjs.clear();
  boostedJetObjs.clear();
  return true;
}

bool UCTSummaryCard::clearRegions() {
  regions.clear();
  return true;
}

const UCTRegion* UCTSummaryCard::getRegion(int regionEtaIndex, uint32_t regionPhiIndex) const {
  if (regionEtaIndex == 0 || (uint32_t)std::abs(regionEtaIndex) > NRegionsInCard ||
      regionPhiIndex >= MaxUCTRegionsPhi) {
    return nullptr;
  }
  UCTGeometry g;
  UCTRegionIndex r = UCTRegionIndex(regionEtaIndex, regionPhiIndex);
  UCTTowerIndex t = g.getUCTTowerIndex(r);
  uint32_t absCaloEta = std::abs(t.first);
  uint32_t absCaloPhi = std::abs(t.second);
  bool negativeEta = false;
  if (t.first < 0)
    negativeEta = true;
  uint32_t nCrate = g.getCrate(absCaloEta, absCaloPhi);
  uint32_t nCard = g.getCard(absCaloEta, absCaloPhi);
  uint32_t nRegion = g.getRegion(absCaloEta, absCaloPhi);
  uint32_t i = ((NCardsInCrate * NRegionsInCard * nCrate) + (NRegionsInCard * nCard) + nRegion) *
               2;  // Correct index for region collection vector of size 3*6*7*2
  if (!negativeEta)
    i++;
  if (i > regions.size()) {
    edm::LogError("L1TCaloSummary") << "UCTSummaryCard: Incorrect region requested -- bailing" << std::endl;
    exit(1);
  }
  return regions[i];
}

bool UCTSummaryCard::setRegionData(std::vector<UCTRegion*> inputRegions) {
  for (long unsigned int i = 0; i < inputRegions.size(); i++) {
    regions.push_back(inputRegions[i]);
  }
  return true;
}

void UCTSummaryCard::print() {
  if (cardSummary > 0)
    edm::LogInfo("L1TCaloSummary") << "UCTSummaryCard: result = " << cardSummary << std::endl;
}
