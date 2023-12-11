// -*- C++ -*-
//
// Package:    L1Trigger/L1CaloTrigger
// Class:      Phase2L1CaloJetEmulator
//
/**\class Phase2L1CaloJetEmulator Phase2L1CaloJetEmulator.cc L1Trigger/L1CaloTrigger/plugins/Phase2L1CaloJetEmulator.cc

 Description: Producing GCT calo jets using GCT barrel, HGCal and HF towers, based on firmware logic.

 Implementation:
     Depends on producers for CaloTowerCollection, HGCalTowerBxCollection and HcalTrigPrimDigiCollection.
*/
//
// Original Author:  Pallabi Das
//         Created:  Tue, 11 Apr 2023 11:27:33 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1TCalorimeterPhase2/interface/CaloCrystalCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloTower.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/CaloPFCluster.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/Phase2L1CaloJet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1THGCal/interface/HGCalTower.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#include <ap_int.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1CaloJetEmulator.h"
#include "TF1.h"

//
// class declaration
//

class Phase2L1CaloJetEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1CaloJetEmulator(const edm::ParameterSet&);
  ~Phase2L1CaloJetEmulator() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  float get_jet_pt_calibration(float &jet_pt, float &jet_eta) const;
  float get_tau_pt_calibration(float &tau_pt, float &tau_eta) const;

  // ----------member data ---------------------------
  edm::EDGetTokenT<l1tp2::CaloTowerCollection> caloTowerToken_;
  edm::EDGetTokenT<l1t::HGCalTowerBxCollection> hgcalTowerToken_;
  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hfToken_;
  edm::ESGetToken<CaloTPGTranscoder, CaloTPGRecord> decoderTag_;
  std::vector<edm::ParameterSet> nHits_to_nvtx_params;
  std::vector<edm::ParameterSet> nvtx_to_PU_sub_params;
  std::map<std::string, TF1> nHits_to_nvtx_funcs;
  std::map<std::string, TF1> hgcalEM_nvtx_to_PU_sub_funcs;
  std::map<std::string, TF1> hgcalHad_nvtx_to_PU_sub_funcs;
  std::map<std::string, TF1> hf_nvtx_to_PU_sub_funcs;
  std::map<std::string, std::map<std::string, TF1> > all_nvtx_to_PU_sub_funcs;

  // For fetching jet pt calibrations
  std::vector<double> jetPtBins;
  std::vector<double> absEtaBinsBarrel;
  std::vector<double> jetCalibrationsBarrel;
  std::vector<double> absEtaBinsHGCal;
  std::vector<double> jetCalibrationsHGCal;
  std::vector<double> absEtaBinsHF;
  std::vector<double> jetCalibrationsHF;

  // For fetching tau pt calibrations
  std::vector<double> tauPtBins;
  std::vector<double> tauAbsEtaBinsBarrel;
  std::vector<double> tauCalibrationsBarrel;
  std::vector<double> tauAbsEtaBinsHGCal;
  std::vector<double> tauCalibrationsHGCal;

  // For storing jet calibrations 
  std::vector<std::vector<double>> calibrationsBarrel;
  std::vector<std::vector<double>> calibrationsHGCal;
  std::vector<std::vector<double>> calibrationsHF;

  // For storing tau calibrations
  std::vector<std::vector<double>> tauPtCalibrationsBarrel;
  std::vector<std::vector<double>> tauPtCalibrationsHGCal;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Phase2L1CaloJetEmulator::Phase2L1CaloJetEmulator(const edm::ParameterSet& iConfig)
    : caloTowerToken_(consumes<l1tp2::CaloTowerCollection>(iConfig.getParameter<edm::InputTag>("gctFullTowers"))),
      hgcalTowerToken_(consumes<l1t::HGCalTowerBxCollection>(iConfig.getParameter<edm::InputTag>("hgcalTowers"))),
      hfToken_(consumes<HcalTrigPrimDigiCollection>(iConfig.getParameter<edm::InputTag>("hcalDigis"))),
      decoderTag_(esConsumes<CaloTPGTranscoder, CaloTPGRecord>(edm::ESInputTag("", ""))),
      nHits_to_nvtx_params(iConfig.getParameter<std::vector<edm::ParameterSet> >("nHits_to_nvtx_params")),
      nvtx_to_PU_sub_params(iConfig.getParameter<std::vector<edm::ParameterSet> >("nvtx_to_PU_sub_params")), 
      jetPtBins(iConfig.getParameter<std::vector<double>>("jetPtBins")),
      absEtaBinsBarrel(iConfig.getParameter<std::vector<double>>("absEtaBinsBarrel")),
      jetCalibrationsBarrel(iConfig.getParameter<std::vector<double>>("jetCalibrationsBarrel")),
      absEtaBinsHGCal(iConfig.getParameter<std::vector<double>>("absEtaBinsHGCal")),
      jetCalibrationsHGCal(iConfig.getParameter<std::vector<double>>("jetCalibrationsHGCal")),
      absEtaBinsHF(iConfig.getParameter<std::vector<double>>("absEtaBinsHF")),
      jetCalibrationsHF(iConfig.getParameter<std::vector<double>>("jetCalibrationsHF")),
      tauPtBins(iConfig.getParameter<std::vector<double>>("tauPtBins")),
      tauAbsEtaBinsBarrel(iConfig.getParameter<std::vector<double>>("tauAbsEtaBinsBarrel")),
      tauCalibrationsBarrel(iConfig.getParameter<std::vector<double>>("tauCalibrationsBarrel")),
      tauAbsEtaBinsHGCal(iConfig.getParameter<std::vector<double>>("tauAbsEtaBinsHGCal")),
      tauCalibrationsHGCal(iConfig.getParameter<std::vector<double>>("tauCalibrationsHGCal")) {
  for (uint i = 0; i < nHits_to_nvtx_params.size(); i++) {
    edm::ParameterSet* pset = &nHits_to_nvtx_params.at(i);
    std::string calo = pset->getParameter<std::string>("fit");
    nHits_to_nvtx_funcs[calo.c_str()] = TF1(calo.c_str(), "[0] + [1] * x");
    nHits_to_nvtx_funcs[calo.c_str()].SetParameter(0, pset->getParameter<std::vector<double> >("params").at(0));
    nHits_to_nvtx_funcs[calo.c_str()].SetParameter(1, pset->getParameter<std::vector<double> >("params").at(1));
  }
  all_nvtx_to_PU_sub_funcs["hgcalEM"] = hgcalEM_nvtx_to_PU_sub_funcs;
  all_nvtx_to_PU_sub_funcs["hgcalHad"] = hgcalHad_nvtx_to_PU_sub_funcs;
  all_nvtx_to_PU_sub_funcs["hf"] = hf_nvtx_to_PU_sub_funcs;
  for (uint i = 0; i < nvtx_to_PU_sub_params.size(); i++) {
    edm::ParameterSet* pset = &nvtx_to_PU_sub_params.at(i);
    std::string calo = pset->getParameter<std::string>("calo");
    std::string iEta = pset->getParameter<std::string>("iEta");
    double p1 = pset->getParameter<std::vector<double> >("params").at(0);
    double p2 = pset->getParameter<std::vector<double> >("params").at(1);

    all_nvtx_to_PU_sub_funcs[calo.c_str()][iEta.c_str()] = TF1(calo.c_str(), "[0] + [1] * x");
    all_nvtx_to_PU_sub_funcs[calo.c_str()][iEta.c_str()].SetParameter(0, p1);
    all_nvtx_to_PU_sub_funcs[calo.c_str()][iEta.c_str()].SetParameter(1, p2);
  }

  // Fill the jet pt calibration 2D vector
  // Dimension 1 is AbsEta bin
  // Dimension 2 is jet pT bin which is filled with the actual callibration value
  // size()-1 b/c the inputs have lower and upper bounds
  // Do Barrel, then HGCal, then HF
  int index = 0;
  for (unsigned int abs_eta = 0; abs_eta < absEtaBinsBarrel.size() - 1; abs_eta++) {
    std::vector<double> pt_bin_calibs;
    for (unsigned int pt = 0; pt < jetPtBins.size() - 1; pt++) {
      pt_bin_calibs.push_back(jetCalibrationsBarrel.at(index));
      index++;
    }
    calibrationsBarrel.push_back(pt_bin_calibs);
  }

  index = 0;
  for (unsigned int abs_eta = 0; abs_eta < absEtaBinsHGCal.size() - 1; abs_eta++) {
    std::vector<double> pt_bin_calibs;
    for (unsigned int pt = 0; pt < jetPtBins.size() - 1; pt++) {
      pt_bin_calibs.push_back(jetCalibrationsHGCal.at(index));
      index++;
    }
    calibrationsHGCal.push_back(pt_bin_calibs);
  }

  index = 0;
  for (unsigned int abs_eta = 0; abs_eta < absEtaBinsHF.size() - 1; abs_eta++) {
    std::vector<double> pt_bin_calibs;
    for (unsigned int pt = 0; pt < jetPtBins.size() - 1; pt++) {
      pt_bin_calibs.push_back(jetCalibrationsHF.at(index));
      index++;
    }
    calibrationsHF.push_back(pt_bin_calibs);
  }

  // Fill the tau pt calibration 2D vector
  // Dimension 1 is AbsEta bin
  // Dimension 2 is tau pT bin which is filled with the actual callibration value
  // size()-1 b/c the inputs have lower and upper bounds (except L1EG b/c that is a cound)
  // Do Barrel, then HGCal
  //
  // Note to future developers: be very concious of the order in which the calibrations are printed
  // out in tool which makse the cfg files.  You need to match that exactly when loading them and
  // using the calibrations below.
  index = 0;
  for (unsigned int abs_eta = 0; abs_eta < tauAbsEtaBinsBarrel.size() - 1; abs_eta++) {
    std::vector<double> pt_bin_calibs;
    for (unsigned int pt = 0; pt < tauPtBins.size() - 1; pt++) {
      pt_bin_calibs.push_back(tauCalibrationsBarrel.at(index));
      index++;
    }
    tauPtCalibrationsBarrel.push_back(pt_bin_calibs);
  }

  index = 0;
  for (unsigned int abs_eta = 0; abs_eta < tauAbsEtaBinsHGCal.size() - 1; abs_eta++) {
    std::vector<double> pt_bin_calibs;
    for (unsigned int pt = 0; pt < tauPtBins.size() - 1; pt++) {
      pt_bin_calibs.push_back(tauCalibrationsHGCal.at(index));
      index++;
    }
    tauPtCalibrationsHGCal.push_back(pt_bin_calibs);
  }

  produces<l1tp2::Phase2L1CaloJetCollection>("GCTJet");
}

Phase2L1CaloJetEmulator::~Phase2L1CaloJetEmulator() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void Phase2L1CaloJetEmulator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  std::unique_ptr<l1tp2::Phase2L1CaloJetCollection> jetCands(make_unique<l1tp2::Phase2L1CaloJetCollection>());

  edm::Handle<std::vector<l1tp2::CaloTower>> caloTowerCollection;
  if (!iEvent.getByToken(caloTowerToken_, caloTowerCollection))
    edm::LogError("Phase2L1CaloJetEmulator") << "Failed to get towers from caloTowerCollection!";

  iEvent.getByToken(caloTowerToken_, caloTowerCollection);
  float GCTintTowers[nBarrelEta][nBarrelPhi];
  float realEta[nBarrelEta][nBarrelPhi];
  float realPhi[nBarrelEta][nBarrelPhi];
  for (const l1tp2::CaloTower& i : *caloTowerCollection) {

    int ieta = i.towerIEta();
    int iphi = i.towerIPhi();
    if (i.ecalTowerEt() > 1.) GCTintTowers[ieta][iphi] = i.ecalTowerEt(); // suppress <= 1 GeV towers
    else GCTintTowers[ieta][iphi] = 0;
    realEta[ieta][iphi] = i.towerEta();
    realPhi[ieta][iphi] = i.towerPhi();
  }

  //Assign ETs to each eta-half of the barrel region (17x72 --> 18x72 to be able to make 3x3 super towers)
  //Then we create 6x24 super towers in each eta-half
  //Find 6 jets in each eta-half

  float temporary[nBarrelEta/2][nBarrelPhi];

  edm::Handle<l1t::HGCalTowerBxCollection> hgcalTowerCollection;
  if (!iEvent.getByToken(hgcalTowerToken_, hgcalTowerCollection))
    edm::LogError("Phase2L1CaloJetEmulator") << "Failed to get towers from hgcalTowerCollection!";
  l1t::HGCalTowerBxCollection hgcalTowerColl;
  iEvent.getByToken(hgcalTowerToken_, hgcalTowerCollection);
  hgcalTowerColl = (*hgcalTowerCollection.product());

  edm::Handle<HcalTrigPrimDigiCollection> hfHandle;
  if (!iEvent.getByToken(hfToken_, hfHandle))
    edm::LogError("Phase2L1CaloJetEmulator") << "Failed to get HcalTrigPrimDigi for HF!";
  iEvent.getByToken(hfToken_, hfHandle);

  int i_hgcalEM_hits_leq_threshold = 0;
  int i_hgcalHad_hits_leq_threshold = 0;
  int i_hf_hits_leq_threshold = 0;
  for (auto it = hgcalTowerColl.begin(0); it != hgcalTowerColl.end(0); it++) {
    if (it->etEm() <= 1.75 && it->etEm() >= 1.25) {
      i_hgcalEM_hits_leq_threshold++;
    }
    if (it->etHad() <= 1.25 && it->etHad() >= 0.75) {
      i_hgcalHad_hits_leq_threshold++;
    }
  }
  const auto& decoder = iSetup.getData(decoderTag_);
  for (const auto& hit : *hfHandle.product()) {
    double et = decoder.hcaletValue(hit.id(), hit.t0());
    if (abs(hit.id().ieta()) < l1t::CaloTools::kHFBegin) continue;
    if (abs(hit.id().ieta()) > l1t::CaloTools::kHFEnd) continue;
    if (et <= 15.0 && et >= 10.0) i_hf_hits_leq_threshold++;
  }

  double hgcalEM_nvtx = nHits_to_nvtx_funcs["hgcalEM"].Eval(i_hgcalEM_hits_leq_threshold);
  if (hgcalEM_nvtx < 0) hgcalEM_nvtx = 0;
  double hgcalHad_nvtx = nHits_to_nvtx_funcs["hgcalHad"].Eval(i_hgcalHad_hits_leq_threshold);
  if (hgcalHad_nvtx < 0) hgcalHad_nvtx = 0;
  double hf_nvtx = nHits_to_nvtx_funcs["hf"].Eval(i_hf_hits_leq_threshold);
  if (hf_nvtx < 0) hf_nvtx = 0;
  double EstimatedNvtx = (hgcalEM_nvtx + hgcalHad_nvtx + hf_nvtx) / 3.;
    
  // HGCal info 
  float hgcalTowers[nHgcalEta][nHgcalPhi];
  float hgcalEta[nHgcalEta][nHgcalPhi];
  float hgcalPhi[nHgcalEta][nHgcalPhi];

  for (int iphi = 0; iphi < nHgcalPhi; iphi++) {
    for (int ieta = 0; ieta < nHgcalEta; ieta++) {
      hgcalTowers[ieta][iphi] = 0;
      if (ieta < nHgcalEta/2) hgcalEta[ieta][iphi] = -3.045 + ieta*0.087 + 0.0435;
      else hgcalEta[ieta][iphi] = 1.479 + (ieta-nHgcalEta/2)*0.087 + 0.0435;
      hgcalPhi[ieta][iphi] = - M_PI + (iphi*M_PI/36) + (M_PI/72);
    }
  }

  for (auto it = hgcalTowerColl.begin(0); it != hgcalTowerColl.end(0); it++) {
    float eta = it->eta();
    int ieta;
    if (eta < 0) ieta = 19 - it->id().iEta();
    else ieta = 20 + it->id().iEta();
    if (eta > 1.479) ieta = ieta - 4;
    int iphi = it->id().iPhi();
    float hgcal_etEm = it->etEm(); 
    float hgcal_etHad = it->etHad();
    if (abs(eta) <= 1.8) {
      hgcal_etEm = it->etEm() - all_nvtx_to_PU_sub_funcs["hgcalEM"]["er1p4to1p8"].Eval(EstimatedNvtx);
      hgcal_etHad = it->etHad() - all_nvtx_to_PU_sub_funcs["hgcalHad"]["er1p4to1p8"].Eval(EstimatedNvtx);
    }
    if (abs(eta) <= 2.1 && abs(eta) > 1.8) {
      hgcal_etEm = it->etEm() - all_nvtx_to_PU_sub_funcs["hgcalEM"]["er1p8to2p1"].Eval(EstimatedNvtx);
      hgcal_etHad = it->etHad() - all_nvtx_to_PU_sub_funcs["hgcalHad"]["er1p8to2p1"].Eval(EstimatedNvtx);
    }
    if (abs(eta) <= 2.4 && abs(eta) > 2.1) {
      hgcal_etEm = it->etEm() - all_nvtx_to_PU_sub_funcs["hgcalEM"]["er2p1to2p4"].Eval(EstimatedNvtx);
      hgcal_etHad = it->etHad() - all_nvtx_to_PU_sub_funcs["hgcalHad"]["er2p1to2p4"].Eval(EstimatedNvtx);
    }
    if (abs(eta) <= 2.7 && abs(eta) > 2.4) {
      hgcal_etEm = it->etEm() - all_nvtx_to_PU_sub_funcs["hgcalEM"]["er2p4to2p7"].Eval(EstimatedNvtx);
      hgcal_etHad = it->etHad() - all_nvtx_to_PU_sub_funcs["hgcalHad"]["er2p4to2p7"].Eval(EstimatedNvtx);
    }
    if (abs(eta) <= 3.1 && abs(eta) > 2.7) {
      hgcal_etEm = it->etEm() - all_nvtx_to_PU_sub_funcs["hgcalEM"]["er2p7to3p1"].Eval(EstimatedNvtx);
      hgcal_etHad = it->etHad() - all_nvtx_to_PU_sub_funcs["hgcalHad"]["er2p7to3p1"].Eval(EstimatedNvtx);
    }
    if (hgcal_etEm < 0) hgcal_etEm = 0;
    if (hgcal_etHad < 0) hgcal_etHad = 0;
    if (hgcal_etEm + hgcal_etHad > 1. && abs(eta) > 1.479) hgcalTowers[ieta][iphi] = hgcal_etEm + hgcal_etHad; // suppress <= 1 GeV towers
    //std::cout<<(it->etEm()+it->etHad())<<"\t"<<(hgcal_etEm + hgcal_etHad)<<std::endl;
  }

  //Assign ETs to each eta-half of the endcap region (18x72)
  //Then we create 6x24 super towers in each eta-half
  //Find 6 jets in each eta-half

  float temporary_hgcal[nHgcalEta/2][nHgcalPhi];

  // HF info
  float hfTowers[nHfEta][nHfPhi];
  float hfEta[nHfEta][nHfPhi];
  float hfPhi[nHfEta][nHfPhi];

  for (int iphi = 0; iphi < nHfPhi; iphi++) {
    for (int ieta = 0; ieta < nHfEta; ieta++) {
      hfTowers[ieta][iphi] = 0;
      int temp = ieta;
      if (ieta < 12) temp = ieta - 41;
      else temp = ieta - 12 + 30;
      hfEta[ieta][iphi] = l1t::CaloTools::towerEta(temp);
      hfPhi[ieta][iphi] = - M_PI + (iphi*M_PI/36) + (M_PI/72);
    }
  }

  //const auto& decoder = iSetup.getData(decoderTag_);
  for (const auto& hit : *hfHandle.product()) {
    double et = decoder.hcaletValue(hit.id(), hit.t0());
    int ieta = 0;
    if (abs(hit.id().ieta()) < l1t::CaloTools::kHFBegin) continue;
    if (abs(hit.id().ieta()) > l1t::CaloTools::kHFEnd) continue;
    if (hit.id().ieta() <= -30) {
      ieta = hit.id().ieta() + 41;
    }
    else if (hit.id().ieta() >= 30) {
      ieta = 12 + (hit.id().ieta() - 30);
    }
    int iphi = 0;
    if (hit.id().iphi() <= 36) iphi = hit.id().iphi() + 35;
    else if (hit.id().iphi() > 36) iphi = hit.id().iphi() - 37;
    if (abs(hit.id().ieta()) <= 33 && abs(hit.id().ieta()) >= 29) et = et - all_nvtx_to_PU_sub_funcs["hf"]["er29to33"].Eval(EstimatedNvtx);
    if (abs(hit.id().ieta()) <= 37 && abs(hit.id().ieta()) >= 34) et = et - all_nvtx_to_PU_sub_funcs["hf"]["er34to37"].Eval(EstimatedNvtx);
    if (abs(hit.id().ieta()) <= 41 && abs(hit.id().ieta()) >= 38) et = et - all_nvtx_to_PU_sub_funcs["hf"]["er38to41"].Eval(EstimatedNvtx);
    if (et < 0) et = 0;
    if (et > 1.) hfTowers[ieta][iphi] = et; // suppress <= 1 GeV towers
  }

  //Assign ETs to each eta-half of the endcap region (12x72)
  //Then we create 4x24 super towers in each eta-half
  //Find 6 jets in each eta-half

  float temporary_hf[nHfEta/2][nHfPhi];

  //Begin creating jets

  vector<l1tp2::Phase2L1CaloJet> halfBarrelJets, halfHgcalJets, halfHfJets;
  halfBarrelJets.clear(); halfHgcalJets.clear(); halfHfJets.clear();
  vector<l1tp2::Phase2L1CaloJet> allJets;
  allJets.clear();

  for (int k = 0; k < 2; k++) {
    halfBarrelJets.clear();
    halfHgcalJets.clear();
    halfHfJets.clear();
    jetInfo jet[30] ;

    // BARREL
    for (int iphi = 0; iphi < nBarrelPhi; iphi++) {
      for (int ieta = 0; ieta < nBarrelEta/2; ieta++) {
        if (k == 0) temporary[ieta][iphi] = GCTintTowers[ieta][iphi];
        else temporary[ieta][iphi] = GCTintTowers[nBarrelEta/2 + ieta][iphi];
      }
    }

    GCTsupertower_t tempST[nSTEta][nSTPhi];
    makeST(temporary, tempST);

    for (int i = 0; i < 10; i++) {
      jet[i] = getRegion(tempST);
      l1tp2::Phase2L1CaloJet tempJet;
      int gctjeteta = jet[i].etaCenter;
      int gctjetphi = jet[i].phiCenter;
      tempJet.setJetIEta(gctjeteta+k*nBarrelEta/2);
      tempJet.setJetIPhi(gctjetphi);
      float jeteta = realEta[gctjeteta+k*nBarrelEta/2][gctjetphi];
      float jetphi = realPhi[gctjeteta+k*nBarrelEta/2][gctjetphi];
      tempJet.setJetEta(jeteta);
      tempJet.setJetPhi(jetphi);
      tempJet.setJetEt(get_jet_pt_calibration(jet[i].energy, jeteta));
      tempJet.setTauEt(get_tau_pt_calibration(jet[i].tauEt, jeteta));
      tempJet.setTowerEt(jet[i].energyMax);
      int gcttowereta = jet[i].etaMax;
      int gcttowerphi = jet[i].phiMax;
      tempJet.setTowerIEta(gcttowereta+k*nBarrelEta/2);
      tempJet.setTowerIPhi(gcttowerphi);
      float towereta = realEta[gcttowereta+k*nBarrelEta/2][gcttowerphi];
      float towerphi = realPhi[gcttowereta+k*nBarrelEta/2][gcttowerphi];
      tempJet.setTowerEta(towereta);
      tempJet.setTowerPhi(towerphi);
      reco::Candidate::PolarLorentzVector tempJetp4;
      tempJetp4.SetPt(tempJet.jetEt());
      tempJetp4.SetEta(tempJet.jetEta());
      tempJetp4.SetPhi(tempJet.jetPhi());
      tempJetp4.SetM(0.);
      tempJet.setP4(tempJetp4);

      if (jet[i].energy > 0.) halfBarrelJets.push_back(tempJet);
    }

    // ENDCAP
    for (int iphi = 0; iphi < nHgcalPhi; iphi++) {
      for (int ieta = 0; ieta < nHgcalEta/2; ieta++) {
        if (k == 0) temporary_hgcal[ieta][iphi] = hgcalTowers[ieta][iphi];
        else temporary_hgcal[ieta][iphi] = hgcalTowers[nHgcalEta/2 + ieta][iphi];
      }
    }

    GCTsupertower_t tempST_hgcal[nSTEta][nSTPhi];
    makeST_hgcal(temporary_hgcal, tempST_hgcal);
    for (int i = 10; i < 20; i++) {
      jet[i] = getRegion(tempST_hgcal);
      l1tp2::Phase2L1CaloJet tempJet;
      int hgcaljeteta = jet[i].etaCenter;
      int hgcaljetphi = jet[i].phiCenter;
      tempJet.setJetIEta(hgcaljeteta+k*nHgcalEta/2);
      tempJet.setJetIPhi(hgcaljetphi);
      float jeteta = hgcalEta[hgcaljeteta+k*nHgcalEta/2][hgcaljetphi];
      float jetphi = hgcalPhi[hgcaljeteta+k*nHgcalEta/2][hgcaljetphi];
      tempJet.setJetEta(jeteta);
      tempJet.setJetPhi(jetphi);
      tempJet.setJetEt(get_jet_pt_calibration(jet[i].energy, jeteta));
      tempJet.setTauEt(get_tau_pt_calibration(jet[i].tauEt, jeteta));
      tempJet.setTowerEt(jet[i].energyMax);
      int hgcaltowereta = jet[i].etaMax;
      int hgcaltowerphi = jet[i].phiMax;
      tempJet.setTowerIEta(hgcaltowereta+k*nHgcalEta/2);
      tempJet.setTowerIPhi(hgcaltowerphi);
      float towereta = hgcalEta[hgcaltowereta+k*nHgcalEta/2][hgcaltowerphi];
      float towerphi = hgcalPhi[hgcaltowereta+k*nHgcalEta/2][hgcaltowerphi];
      tempJet.setTowerEta(towereta);
      tempJet.setTowerPhi(towerphi);
      reco::Candidate::PolarLorentzVector tempJetp4;
      tempJetp4.SetPt(tempJet.jetEt());
      tempJetp4.SetEta(tempJet.jetEta());
      tempJetp4.SetPhi(tempJet.jetPhi());
      tempJetp4.SetM(0.);
      tempJet.setP4(tempJetp4);

      if (jet[i].energy > 0.) halfHgcalJets.push_back(tempJet);
    }

    // HF
    for (int iphi = 0; iphi < nHfPhi; iphi++) {
      for (int ieta = 0; ieta < nHfEta/2; ieta++) {
        if (k == 0) temporary_hf[ieta][iphi] = hfTowers[ieta][iphi];
        else temporary_hf[ieta][iphi] = hfTowers[nHfEta/2 + ieta][iphi];
      }
    }

    GCTsupertower_t tempST_hf[nSTEta][nSTPhi];
    makeST_hf(temporary_hf, tempST_hf);
    for (int i = 20; i < 30; i++) {
      jet[i] = getRegion(tempST_hf);
      l1tp2::Phase2L1CaloJet tempJet;
      int hfjeteta = jet[i].etaCenter;
      int hfjetphi = jet[i].phiCenter;
      tempJet.setJetIEta(hfjeteta+k*nHfEta/2);
      tempJet.setJetIPhi(hfjetphi);
      float jeteta = hfEta[hfjeteta+k*nHfEta/2][hfjetphi];
      float jetphi = hfPhi[hfjeteta+k*nHfEta/2][hfjetphi];
      tempJet.setJetEta(jeteta);
      tempJet.setJetPhi(jetphi);
      tempJet.setJetEt(get_jet_pt_calibration(jet[i].energy, jeteta));
      tempJet.setTauEt(get_tau_pt_calibration(jet[i].tauEt, jeteta));
      tempJet.setTowerEt(jet[i].energyMax);
      int hftowereta = jet[i].etaMax;
      int hftowerphi = jet[i].phiMax;
      tempJet.setTowerIEta(hftowereta+k*nHfEta/2);
      tempJet.setTowerIPhi(hftowerphi);
      float towereta = hfEta[hftowereta+k*nHfEta/2][hftowerphi];
      float towerphi = hfPhi[hftowereta+k*nHfEta/2][hftowerphi];
      tempJet.setTowerEta(towereta);
      tempJet.setTowerPhi(towerphi);
      reco::Candidate::PolarLorentzVector tempJetp4;
      tempJetp4.SetPt(tempJet.jetEt());
      tempJetp4.SetEta(tempJet.jetEta());
      tempJetp4.SetPhi(tempJet.jetPhi());
      tempJetp4.SetM(0.);
      tempJet.setP4(tempJetp4);

      if (jet[i].energy > 0.) halfHfJets.push_back(tempJet);
    }

    // Stitching:
    // if the jet eta is at the boundary: for HB it should be within 0-1 in -ve eta, 32-33 in +ve eta; for HE it should be within 0-1/16-17 in -ve eta, 34-35/18-19 in +ve eta; for HF it should be within 10-11 in -ve eta, 12-13 in +ve eta
    // then get the phi of that jet and check if there is a neighbouring jet with the same phi, then merge to the jet that has higher ET
    // in both eta/phi allow a maximum of one tower between jet centers for stitching

    for (size_t i = 0; i < halfHgcalJets.size(); i++) {
      if (halfHgcalJets.at(i).jetIEta() > 15 && halfHgcalJets.at(i).jetIEta() < 20) {
        float hgcal_ieta = k*nBarrelEta + halfHgcalJets.at(i).jetIEta();
        for (size_t j = 0; j < halfBarrelJets.size(); j++) {
          float barrel_ieta = nHgcalEta/2 + halfBarrelJets.at(j).jetIEta();
          if (abs(barrel_ieta - hgcal_ieta) < 3 && abs(halfBarrelJets.at(j).jetIPhi() - halfHgcalJets.at(i).jetIPhi()) < 3) {
            float totalet = halfBarrelJets.at(j).jetEt() + halfHgcalJets.at(i).jetEt();
	    float totalTauEt = halfBarrelJets.at(j).tauEt() + halfHgcalJets.at(i).tauEt();
            if (halfBarrelJets.at(j).jetEt() > halfHgcalJets.at(i).jetEt()) {
              halfHgcalJets.at(i).setJetEt(0.);
	      halfHgcalJets.at(i).setTauEt(0.);
              halfBarrelJets.at(j).setJetEt(totalet);
	      halfBarrelJets.at(j).setTauEt(totalTauEt);
              reco::Candidate::PolarLorentzVector tempJetp4;
              tempJetp4.SetPt(totalet);
              tempJetp4.SetEta(halfBarrelJets.at(j).jetEta());
              tempJetp4.SetPhi(halfBarrelJets.at(j).jetPhi());
              tempJetp4.SetM(0.);
              halfBarrelJets.at(j).setP4(tempJetp4);
            }
            else {
              halfHgcalJets.at(i).setJetEt(totalet);
	      halfHgcalJets.at(i).setTauEt(totalTauEt);
              halfBarrelJets.at(j).setJetEt(0.);
	      halfBarrelJets.at(j).setTauEt(0.);
              reco::Candidate::PolarLorentzVector tempJetp4;
              tempJetp4.SetPt(totalet);
              tempJetp4.SetEta(halfHgcalJets.at(i).jetEta());
              tempJetp4.SetPhi(halfHgcalJets.at(i).jetPhi());
              tempJetp4.SetM(0.);
              halfHgcalJets.at(i).setP4(tempJetp4);
            }
          }
        }
      }
      else if (halfHgcalJets.at(i).jetIEta() < 2 || halfHgcalJets.at(i).jetIEta() > 33) {
        float hgcal_ieta = k*nBarrelEta + nHfEta/2 + halfHgcalJets.at(i).jetIEta();
        for (size_t j = 0; j < halfHfJets.size(); j++) {
          float hf_ieta = k*nBarrelEta + k*nHgcalEta + halfHfJets.at(j).jetIEta();
          if (abs(hgcal_ieta - hf_ieta) < 3 && abs(halfHfJets.at(j).jetIPhi() - halfHgcalJets.at(i).jetIPhi()) < 3) {
            float totalet = halfHfJets.at(j).jetEt() + halfHgcalJets.at(i).jetEt();
	    float totalTauEt = halfHfJets.at(j).tauEt() + halfHgcalJets.at(i).tauEt();
            if (halfHfJets.at(j).jetEt() > halfHgcalJets.at(i).jetEt()) {
              halfHgcalJets.at(i).setJetEt(0.);
	      halfHgcalJets.at(i).setTauEt(0.);
              halfHfJets.at(j).setJetEt(totalet);
	      halfHfJets.at(j).setTauEt(totalTauEt);
              reco::Candidate::PolarLorentzVector tempJetp4;
              tempJetp4.SetPt(totalet);
              tempJetp4.SetEta(halfHfJets.at(j).jetEta());
              tempJetp4.SetPhi(halfHfJets.at(j).jetPhi());
              tempJetp4.SetM(0.);
              halfHfJets.at(j).setP4(tempJetp4);
            }
            else {
              halfHgcalJets.at(i).setJetEt(totalet);
	      halfHgcalJets.at(i).setTauEt(totalTauEt);
              halfHfJets.at(j).setJetEt(0.);
	      halfHfJets.at(j).setTauEt(0.);
              reco::Candidate::PolarLorentzVector tempJetp4;
              tempJetp4.SetPt(totalet);
              tempJetp4.SetEta(halfHgcalJets.at(i).jetEta());
              tempJetp4.SetPhi(halfHgcalJets.at(i).jetPhi());
              tempJetp4.SetM(0.);
              halfHgcalJets.at(i).setP4(tempJetp4);
            }
          }
        }
      }
    }

    std::sort(halfBarrelJets.begin(), halfBarrelJets.end(), compareByEt);
    for (size_t i = 0; i < halfBarrelJets.size(); i++) {
      if (halfBarrelJets.at(i).jetEt() > 0. && i < 6) allJets.push_back(halfBarrelJets.at(i));
    }

    std::sort(halfHgcalJets.begin(), halfHgcalJets.end(), compareByEt);
    for (size_t i = 0; i < halfHgcalJets.size(); i++) {
      if (halfHgcalJets.at(i).jetEt() > 0. && i < 6) allJets.push_back(halfHgcalJets.at(i));
    }

    std::sort(halfHfJets.begin(), halfHfJets.end(), compareByEt);
    for (size_t i = 0; i < halfHfJets.size(); i++) {
      if (halfHfJets.at(i).jetEt() > 0. && i < 6) allJets.push_back(halfHfJets.at(i));
    }

  }

  std::sort(allJets.begin(), allJets.end(), compareByEt);
  for (size_t i = 0; i < allJets.size(); i++) {
    jetCands->push_back(allJets.at(i));
  }  

  iEvent.put(std::move(jetCands), "GCTJet");
}

// ------------ method called when starting to processes a run  ------------
/*
void
Phase2L1CaloJetEmulator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
Phase2L1CaloJetEmulator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
Phase2L1CaloJetEmulator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
Phase2L1CaloJetEmulator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phase2L1CaloJetEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// Apply calibrations to HCAL energy based on Jet Eta, Jet pT
float Phase2L1CaloJetEmulator::get_jet_pt_calibration(float &jet_pt, float &jet_eta) const {
  float abs_eta = std::abs(jet_eta);
  float tmp_jet_pt = jet_pt;
  if (tmp_jet_pt > 499)
    tmp_jet_pt = 499;

  // Different indices sizes in different calo regions.
  // Barrel...
  size_t eta_index = 0;
  size_t pt_index = 0;
  float calib = 1.0;
  if (abs_eta <= 1.5) {
    // Start loop checking 2nd value
    for (unsigned int i = 1; i < absEtaBinsBarrel.size(); i++) {
      if (abs_eta <= absEtaBinsBarrel.at(i))
        break;
      eta_index++;
    }

    // Start loop checking 2nd value
    for (unsigned int i = 1; i < jetPtBins.size(); i++) {
      if (tmp_jet_pt <= jetPtBins.at(i))
        break;
      pt_index++;
    }
    calib = calibrationsBarrel[eta_index][pt_index];
  }                         // end Barrel
  else if (abs_eta <= 3.0)  // HGCal
    {
      // Start loop checking 2nd value
      for (unsigned int i = 1; i < absEtaBinsHGCal.size(); i++) {
	if (abs_eta <= absEtaBinsHGCal.at(i))
	  break;
	eta_index++;
      }

      // Start loop checking 2nd value
      for (unsigned int i = 1; i < jetPtBins.size(); i++) {
	if (tmp_jet_pt <= jetPtBins.at(i))
	  break;
	pt_index++;
      }
      calib = calibrationsHGCal[eta_index][pt_index];
    }     // end HGCal
  else  // HF
    {
      // Start loop checking 2nd value
      for (unsigned int i = 1; i < absEtaBinsHF.size(); i++) {
	if (abs_eta <= absEtaBinsHF.at(i))
	  break;
	eta_index++;
      }

      // Start loop checking 2nd value
      for (unsigned int i = 1; i < jetPtBins.size(); i++) {
	if (tmp_jet_pt <= jetPtBins.at(i))
	  break;
	pt_index++;
      }
      calib = calibrationsHF[eta_index][pt_index];
    }  // end HF

  return jet_pt*calib;
}

// Apply calibrations to tau pT based on L1EG info, EM Fraction, Tau Eta, Tau pT
float Phase2L1CaloJetEmulator::get_tau_pt_calibration(float &tau_pt, float &tau_eta) const {
  float abs_eta = std::abs(tau_eta);
  float tmp_tau_pt = tau_pt;
  if (tmp_tau_pt > 199)
    tmp_tau_pt = 199;

  // Different indices sizes in different calo regions.
  // Barrel...
  size_t eta_index = 0;
  size_t pt_index = 0;
  float calib = 1.0;
  // HERE
  if (abs_eta <= 1.5) {
    // Start loop checking 2nd value
    for (unsigned int i = 1; i < tauAbsEtaBinsBarrel.size(); i++) {
      if (abs_eta <= tauAbsEtaBinsBarrel.at(i))
        break;
      eta_index++;
    }

    // Start loop checking 2nd value
    for (unsigned int i = 1; i < tauPtBins.size(); i++) {
      if (tmp_tau_pt <= tauPtBins.at(i))
        break;
      pt_index++;
    }
    calib = tauPtCalibrationsBarrel[eta_index][pt_index];
  }                         // end Barrel
  else if (abs_eta <= 3.0)  // HGCal
    {
      // Start loop checking 2nd value
      for (unsigned int i = 1; i < tauAbsEtaBinsHGCal.size(); i++) {
	if (abs_eta <= tauAbsEtaBinsHGCal.at(i))
	  break;
	eta_index++;
      }

      // Start loop checking 2nd value
      for (unsigned int i = 1; i < tauPtBins.size(); i++) {
	if (tmp_tau_pt <= tauPtBins.at(i))
	  break;
	pt_index++;
      }
      calib = tauPtCalibrationsHGCal[eta_index][pt_index];
    }  // end HGCal
  else
    return calib;

  return tau_pt*calib;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1CaloJetEmulator);
