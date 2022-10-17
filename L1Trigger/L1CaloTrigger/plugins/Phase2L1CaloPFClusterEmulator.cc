// -*- C++ -*-
//
// Package:    L1Trigger/L1CaloTrigger
// Class:      Phase2L1CaloPFClusterEmulator
//
/**\class Phase2L1CaloPFClusterEmulator Phase2L1CaloPFClusterEmulator.cc L1Trigger/L1CaloTrigger/plugins/Phase2L1CaloPFClusterEmulator.cc

 Description: Creates 3x3 PF clusters from GCTintTowers to be sent to correlator. Follows firmware logic, creates 8 clusters per (2+17+2)x(2+4+2).

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Pallabi Das
//         Created:  Wed, 21 Sep 2022 14:54:20 GMT
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
#include "DataFormats/L1Trigger/interface/EGamma.h"

#include <ap_int.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include "Phase2L1CaloPFClusterEmulator.h"

//#include "Phase2L1RCT.h"
//#include "Phase2L1GCT.h"
//#include "Phase2L1GCT_algo.h"


//
// class declaration
//

class Phase2L1CaloPFClusterEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1CaloPFClusterEmulator(const edm::ParameterSet&);
  ~Phase2L1CaloPFClusterEmulator();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  //vector<l1tp2::CaloTower>              "Phase2L1CaloEGammaEmulatorProducer"   "GCTFullTowers"   "L1AlgoTest"  
  edm::EDGetTokenT<l1tp2::CaloTowerCollection> caloTowerToken_;
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
Phase2L1CaloPFClusterEmulator::Phase2L1CaloPFClusterEmulator(const edm::ParameterSet& iConfig) :
  caloTowerToken_(consumes<l1tp2::CaloTowerCollection>( edm::InputTag("Phase2L1CaloEGammaEmulatorProducer","GCTFullTowers","L1AlgoTest")))
{
  produces<l1tp2::CaloPFClusterCollection>( "GCTPFCluster" ) ;
  //register your products
/* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
*/
  //now do what ever other initialization is needed
}

Phase2L1CaloPFClusterEmulator::~Phase2L1CaloPFClusterEmulator() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void Phase2L1CaloPFClusterEmulator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  std::unique_ptr<l1tp2::CaloPFClusterCollection> pfclusterCands(new l1tp2::CaloPFClusterCollection);

  edm::Handle<std::vector<l1tp2::CaloTower>> caloTowerCollection;
  if(!iEvent.getByToken(caloTowerToken_, caloTowerCollection)) edm::LogError("Phase2L1CaloPFClusterEmulator") << "Failed to get towers from caloTowerCollection!" ;

  iEvent.getByToken(caloTowerToken_, caloTowerCollection);
  float GCTintTowers[34][72];
  float realEta[34][72];
  float realPhi[34][72];
  for (const l1tp2::CaloTower &i : *caloTowerCollection) {
    //std::cout << "Phase2L1CaloPFClusterEmulator: GCT FULL tower energy " << i.ecalTowerEt()
    //          << " at eta, phi: (" << i.towerIEta() << ", " << i.towerIPhi() << ")" << std::endl; 

    int ieta = i.towerIEta(); int iphi = i.towerIPhi();
    GCTintTowers[ieta][iphi] = i.ecalTowerEt();
    realEta[ieta][iphi] = i.towerEta();
    realPhi[ieta][iphi] = i.towerPhi();

    //l1tp2::CaloPFCluster l1CaloPFCluster;
    //l1CaloPFCluster.setClusterEt(i.ecalTowerEt());
    //pfclusterCands->push_back(l1CaloPFCluster);
  }

  float regions[36][21][8];

  for(int ieta = 0; ieta < 21; ieta++){
    for(int iphi = 0; iphi < 8; iphi++){
      for(int k = 0; k < 36; k++){
        regions[k][ieta][iphi] = 0.;
      }
    }
  }

  //Assign ETs to 36 21x8 regions

  for(int ieta = 0; ieta < 21; ieta++){
    for(int iphi = 0; iphi < 8; iphi++){
      if(ieta > 1){
        if(iphi > 1) regions[0][ieta][iphi] = GCTintTowers[ieta-2][iphi-2];
        regions[2][ieta][iphi] = GCTintTowers[ieta-2][iphi+2];
        regions[4][ieta][iphi] = GCTintTowers[ieta-2][iphi+6];
        regions[6][ieta][iphi] = GCTintTowers[ieta-2][iphi+10];
        regions[8][ieta][iphi] = GCTintTowers[ieta-2][iphi+14];
        regions[10][ieta][iphi] = GCTintTowers[ieta-2][iphi+18];
        regions[12][ieta][iphi] = GCTintTowers[ieta-2][iphi+22];
        regions[14][ieta][iphi] = GCTintTowers[ieta-2][iphi+26];
        regions[16][ieta][iphi] = GCTintTowers[ieta-2][iphi+30];
        regions[18][ieta][iphi] = GCTintTowers[ieta-2][iphi+34];
        regions[20][ieta][iphi] = GCTintTowers[ieta-2][iphi+38];
        regions[22][ieta][iphi] = GCTintTowers[ieta-2][iphi+42];
        regions[24][ieta][iphi] = GCTintTowers[ieta-2][iphi+46];
        regions[26][ieta][iphi] = GCTintTowers[ieta-2][iphi+50];
        regions[28][ieta][iphi] = GCTintTowers[ieta-2][iphi+54];
        regions[30][ieta][iphi] = GCTintTowers[ieta-2][iphi+58];
        regions[32][ieta][iphi] = GCTintTowers[ieta-2][iphi+62];
        if(iphi < 6) regions[34][ieta][iphi] = GCTintTowers[ieta-2][iphi+66];
      }
    if(ieta < 19){
        if(iphi > 1) regions[1][ieta][iphi] = GCTintTowers[ieta+15][iphi-2];
        regions[3][ieta][iphi] = GCTintTowers[ieta+15][iphi+2];
        regions[5][ieta][iphi] = GCTintTowers[ieta+15][iphi+6];
        regions[7][ieta][iphi] = GCTintTowers[ieta+15][iphi+10];
        regions[9][ieta][iphi] = GCTintTowers[ieta+15][iphi+14];
        regions[11][ieta][iphi] = GCTintTowers[ieta+15][iphi+18];
        regions[13][ieta][iphi] = GCTintTowers[ieta+15][iphi+22];
        regions[15][ieta][iphi] = GCTintTowers[ieta+15][iphi+26];
        regions[17][ieta][iphi] = GCTintTowers[ieta+15][iphi+30];
        regions[19][ieta][iphi] = GCTintTowers[ieta+15][iphi+34];
        regions[21][ieta][iphi] = GCTintTowers[ieta+15][iphi+38];
        regions[23][ieta][iphi] = GCTintTowers[ieta+15][iphi+42];
        regions[25][ieta][iphi] = GCTintTowers[ieta+15][iphi+46];
        regions[27][ieta][iphi] = GCTintTowers[ieta+15][iphi+50];
        regions[29][ieta][iphi] = GCTintTowers[ieta+15][iphi+54];
        regions[31][ieta][iphi] = GCTintTowers[ieta+15][iphi+58];
        regions[33][ieta][iphi] = GCTintTowers[ieta+15][iphi+62];
        if(iphi < 6) regions[35][ieta][iphi] = GCTintTowers[ieta+15][iphi+66];
      }
    }
  }

  float temporary[21][8];
  int etaoffset = 0;
  int phioffset = 0;

  //Use same code from firmware for finding clusters
  for(int k = 0;  k < 36; k++){
    for(int ieta = 0; ieta < 21; ieta++){
      for(int iphi = 0; iphi < 8; iphi++){
        temporary[ieta][iphi] = regions[k][ieta][iphi];
      }
    }
    if(k % 2 ==  0) etaoffset = 0;
    else etaoffset = 17;
    if(k > 1 && k % 2 ==  0) phioffset = phioffset+4;

    GCTPfcluster_t tempPfclusters;
    tempPfclusters = pfcluster(temporary, etaoffset, phioffset);  

    for(int i=0; i<8; i++){
      int gcteta = tempPfclusters.GCTpfclusters[i].eta;
      int gctphi = tempPfclusters.GCTpfclusters[i].phi;
      float towereta = realEta[gcteta][gctphi];
      float towerphi = realPhi[gcteta][gctphi];
      //std::cout << "Phase2L1CaloPFClusterEmulator: GCT Pfclusters energy " << tempPfclusters.GCTpfclusters[i].et
      //          << " at ieta, iphi: (" << tempPfclusters.GCTpfclusters[i].eta << ", " << tempPfclusters.GCTpfclusters[i].phi << ")"
      //          << " at eta, phi: (" << towereta << ", " << towerphi << ")" << std::endl;
      l1tp2::CaloPFCluster l1CaloPFCluster;
      l1CaloPFCluster.setClusterEt(tempPfclusters.GCTpfclusters[i].et);
      l1CaloPFCluster.setClusterIEta(gcteta);
      l1CaloPFCluster.setClusterIPhi(gctphi);
      l1CaloPFCluster.setClusterEta(towereta);
      l1CaloPFCluster.setClusterPhi(towerphi);
      pfclusterCands->push_back(l1CaloPFCluster);
      //allPFClusters.push_back(tempPfclusters.GCTpfclusters[i]);
    }
  }

  //cout<<"PFClusters size: "<<pfclusterCands->size()<<endl;
  iEvent.put(std::move(pfclusterCands), "GCTPFCluster");
    
/*

  //float et_c, et_n, et_s, et_e, et_w, et_ne, et_nw, et_se, et_sw;
  vector<float> allPFClusters;
  allPFClusters.clear();
  for(int i=0; i<8; i++){
    cout<<"Pfclusters1: "<<Pfclusters1.GCTpfclusters[i].et<<endl;
     allPFClusters.push_back(Pfclusters1.GCTpfclusters[i].et);
  }

  float pfclusteret_3x3[34][72];
  for(int ieta = 0; ieta < 34; ieta++){
    for(int iphi = 0; iphi < 72; iphi++){
      et_c = GCTintTowers[ieta][iphi];

      if(iphi == 0) et_n = GCTintTowers[ieta][71];
      else if(iphi > 0) et_n = GCTintTowers[ieta][iphi-1];

      if(iphi == 71) et_s = GCTintTowers[ieta][0];
      else if(iphi < 71) et_s = GCTintTowers[ieta][iphi+1];

      if(ieta == 33) et_e = 0;
      else if(ieta < 33) et_e = GCTintTowers[ieta+1][iphi];

      if(ieta == 0) et_w = 0;
      else if(ieta > 0) et_w = GCTintTowers[ieta-1][iphi];

      if(ieta == 33) et_ne = 0;
      else if(iphi == 0 && ieta < 33) et_ne = GCTintTowers[ieta+1][71];
      else if(iphi > 0 && ieta < 33) et_ne = GCTintTowers[ieta+1][iphi-1];

      if(ieta == 0) et_nw = 0;
      else if(iphi == 0 && ieta > 0) et_nw = GCTintTowers[ieta-1][71];
      else if(iphi > 0 && ieta > 0) et_nw = GCTintTowers[ieta-1][iphi-1];

      if(ieta == 33) et_se = 0;
      else if(iphi == 71 && ieta < 33) et_se = GCTintTowers[ieta+1][0];
      else if(iphi < 71 && ieta < 33) et_ne = GCTintTowers[ieta+1][iphi+1];

      if(ieta == 0) et_sw = 0;
      else if(iphi == 71 && ieta > 0) et_sw = GCTintTowers[ieta-1][0];
      else if(iphi < 71 && ieta > 0) et_sw = GCTintTowers[ieta-1][iphi+1];
     
      float et_3x3 = et_c + et_n + et_s + et_e + et_w + et_ne + et_nw + et_se + et_sw;
      if(et_c >= et_n && et_c >= et_s && et_c >= et_e && et_c >= et_w && et_c >= et_ne && et_c >= et_nw && et_c >= et_se && et_c >= sw) pfclusteret_3x3[ieta][iphi] = et_3x3;
      else pfclusteret_3x3[ieta][iphi] = 0;

      allPFClusters.push_back(et_3x3);
    }
  }
  if(allPFClusters.size() > 1){ std::sort(allPFClusters.begin(),allPFClusters.end(),std::greater<float>()); }
  for(long unsigned int i = 0; i < allPFClusters.size(); i++){
    std::cout << "Phase2L1CaloPFClusterEmulator: GCT PF cluster energy " << allPFClusters.at(i) << std::endl;
    l1tp2::CaloPFCluster l1CaloPFCluster;
    l1CaloPFCluster.setClusterEt(allPFClusters.at(i));
    pfclusterCands->push_back(l1CaloPFCluster);
  //  if(i < 48) ...
  }
*/

/* This is an event example
  //Read 'ExampleData' from the Event
  ExampleData const& in = iEvent.get(inToken_);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  iEvent.put(std::make_unique<ExampleData2>(in));
*/

/* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
*/
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void Phase2L1CaloPFClusterEmulator::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void Phase2L1CaloPFClusterEmulator::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
Phase2L1CaloPFClusterEmulator::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
Phase2L1CaloPFClusterEmulator::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
Phase2L1CaloPFClusterEmulator::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
Phase2L1CaloPFClusterEmulator::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phase2L1CaloPFClusterEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1CaloPFClusterEmulator);
