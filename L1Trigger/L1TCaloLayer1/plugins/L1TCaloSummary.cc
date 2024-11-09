// -*- C++ -*-
//
// Package:    L1Trigger/L1TCaloSummary
// Class:      L1TCaloSummary
//
/**\class L1TCaloSummary L1TCaloSummary.cc L1Trigger/L1TCaloSummary/plugins/L1TCaloSummary.cc

   Description: The package L1Trigger/L1TCaloSummary is prepared for monitoring the CMS Layer-1 Calorimeter Trigger.

   Implementation:
   It prepares region objects and puts them in the event
*/
//
// Original Author:  Sridhara Dasu
//         Created:  Sat, 14 Nov 2015 14:18:27 GMT
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

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "L1Trigger/L1TCaloLayer1/src/UCTLayer1.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCrate.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTCard.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTTower.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"

#include "L1Trigger/L1TCaloLayer1/src/UCTObject.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTSummaryCard.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometryExtended.hh"

#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/CICADA.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "L1Trigger/L1TCaloLayer1/src/UCTLogging.hh"
#include <bitset>

#include <string>
#include <sstream>

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"


//Anomaly detection includes
#include "ap_fixed.h"
#include "hls4ml/emulator.h"

using namespace l1tcalo;
using namespace l1extra;
using namespace std;

//
// class declaration
//

template <class INPUT, class OUTPUT>
class L1TCaloSummary : public edm::stream::EDProducer<> {
public:
  explicit L1TCaloSummary(const edm::ParameterSet&);
  ~L1TCaloSummary() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  //void beginJob() override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  //void endJob() override;

  void print();

  // ----------member data ---------------------------

  uint32_t nPumBins;

  std::vector<std::vector<std::vector<uint32_t>>> pumLUT;

  double caloScaleFactor;

  uint32_t jetSeed;
  uint32_t tauSeed;
  float tauIsolationFactor;
  uint32_t eGammaSeed;
  double eGammaIsolationFactor;
  double boostedJetPtFactor;

  bool verbose;
  int fwVersion;

  edm::EDGetTokenT<L1CaloRegionCollection> regionToken;

  UCTLayer1* layer1;

  hls4mlEmulator::ModelLoader loader;
  std::shared_ptr<hls4mlEmulator::Model> model;

  bool overwriteWithTestPatterns;
  std::vector<edm::ParameterSet> testPatterns;

  tensorflow::Options options;
  tensorflow::MetaGraphDef* metaGraph;
  tensorflow::Session* session; 
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
template <class INPUT, class OUTPUT>
L1TCaloSummary<INPUT, OUTPUT>::L1TCaloSummary(const edm::ParameterSet& iConfig)
    : nPumBins(iConfig.getParameter<unsigned int>("nPumBins")),
      pumLUT(nPumBins, std::vector<std::vector<uint32_t>>(2, std::vector<uint32_t>(13))),
      caloScaleFactor(iConfig.getParameter<double>("caloScaleFactor")),
      jetSeed(iConfig.getParameter<unsigned int>("jetSeed")),
      tauSeed(iConfig.getParameter<unsigned int>("tauSeed")),
      tauIsolationFactor(iConfig.getParameter<double>("tauIsolationFactor")),
      eGammaSeed(iConfig.getParameter<unsigned int>("eGammaSeed")),
      eGammaIsolationFactor(iConfig.getParameter<double>("eGammaIsolationFactor")),
      boostedJetPtFactor(iConfig.getParameter<double>("boostedJetPtFactor")),
      verbose(iConfig.getParameter<bool>("verbose")),
      fwVersion(iConfig.getParameter<int>("firmwareVersion")),
      regionToken(consumes<L1CaloRegionCollection>(edm::InputTag("simCaloStage2Layer1Digis"))),
      loader(hls4mlEmulator::ModelLoader(iConfig.getParameter<string>("CICADAModelVersion"))),
      overwriteWithTestPatterns(iConfig.getParameter<bool>("useTestPatterns")),
      testPatterns(iConfig.getParameter<std::vector<edm::ParameterSet>>("testPatterns")) {
  std::vector<double> pumLUTData;
  char pumLUTString[10];
  for (uint32_t pumBin = 0; pumBin < nPumBins; pumBin++) {
    for (uint32_t side = 0; side < 2; side++) {
      if (side == 0)
        sprintf(pumLUTString, "pumLUT%2.2dp", pumBin);
      else
        sprintf(pumLUTString, "pumLUT%2.2dn", pumBin);
      pumLUTData = iConfig.getParameter<std::vector<double>>(pumLUTString);
      for (uint32_t iEta = 0; iEta < std::max((uint32_t)pumLUTData.size(), MaxUCTRegionsEta); iEta++) {
        pumLUT[pumBin][side][iEta] = (uint32_t)round(pumLUTData[iEta] / caloScaleFactor);
      }
      if (pumLUTData.size() != (MaxUCTRegionsEta))
        edm::LogError("L1TCaloSummary") << "PUM LUT Data size integrity check failed; Expected size = "
                                        << MaxUCTRegionsEta << "; Provided size = " << pumLUTData.size()
                                        << "; Will use what is provided :(" << std::endl;
    }
  }
  produces<L1JetParticleCollection>("Boosted");

  //anomaly trigger loading
  model = loader.load_model();
  produces<l1t::CICADABxCollection>("CICADAScore");

  std::string fullPathToModel("/afs/hep.wisc.edu/user/mbileska/100try/CMSSW_14_1_0_pre5/src/L1Trigger/L1TCaloLayer1/data/model/");
  metaGraph = tensorflow::loadMetaGraphDef(fullPathToModel);
  session = tensorflow::createSession(metaGraph, fullPathToModel, options);
}


// Destructor
template <class INPUT, class OUTPUT>
L1TCaloSummary<INPUT, OUTPUT>::~L1TCaloSummary()
{
  delete metaGraph;
  metaGraph = nullptr;
  tensorflow::closeSession(session);
  session=nullptr;
}

//
// member functions
//
std::vector<float> eta_low_bound = {-2.5 , -2.172 , -1.392 , -1.044, -0.696, -0.348, 0.001, 0.348, 0.696, 1.044, 1.392 , 1.74, 2.172, 2.5};
std::array<float, 18> central_phi = {{0.000, 0.349, 0.698, 1.047, 1.396, 1.744, 2.093, 2.442, 2.791, -3.14159, -2.791, -2.442, -2.093, -1.744, -1.396, -1.047, -0.698, -0.349}};

float Rct_to_eta(float ieta) {
    if (ieta < 0 || ieta >= eta_low_bound.size() - 1) {
        return -999; 
    }
    float eta_lower = eta_low_bound[ieta];
    float eta_upper = eta_low_bound[ieta + 1];
    return (eta_lower + eta_upper) / 2.0;
}

float Rct_to_phi(float iphi) {
    if (iphi < 0 || iphi >= central_phi.size()) {
        return -999;  
    }
    return central_phi[iphi];  
}
// ------------ method called to produce the data  ------------
template <class INPUT, class OUTPUT>
void L1TCaloSummary<INPUT, OUTPUT>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::unique_ptr<L1JetParticleCollection> bJetCands(new L1JetParticleCollection);

  std::unique_ptr<l1t::CICADABxCollection> CICADAScore = std::make_unique<l1t::CICADABxCollection>();
  CICADAScore->setBXRange(-2, 2);

  UCTGeometry g;

  // Here we read region data from the region collection created by L1TCaloLayer1 instead of
  // independently creating regions from TPGs for processing by the summary card. This results
  // in a single region vector of size 252 whereas from independent creation we had 3*6 vectors
  // of size 7*2. Indices are mapped in UCTSummaryCard accordingly.
  UCTSummaryCard summaryCard =
      UCTSummaryCard(&pumLUT, jetSeed, tauSeed, tauIsolationFactor, eGammaSeed, eGammaIsolationFactor);
  std::vector<UCTRegion*> inputRegions;
  inputRegions.clear();
  edm::Handle<std::vector<L1CaloRegion>> regionCollection;
  if (!iEvent.getByToken(regionToken, regionCollection))
    edm::LogError("L1TCaloSummary") << "UCT: Failed to get regions from region collection!";
  iEvent.getByToken(regionToken, regionCollection);
  //Model input
  //This is done as a flat vector input, but future versions may involve 2D input
  //This will have to be handled later
  INPUT modelInput[252];
  tensorflow::Tensor input(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, 18, 14, 1}));
  auto inputTensorMapped = input.tensor<float, 4>();

  for (const L1CaloRegion& i : *regionCollection) {
    UCTRegionIndex r = g.getUCTRegionIndexFromL1CaloRegion(i.gctEta(), i.gctPhi());
    UCTTowerIndex t = g.getUCTTowerIndexFromL1CaloRegion(r, i.raw());
    uint32_t absCaloEta = std::abs(t.first);
    uint32_t absCaloPhi = std::abs(t.second);
    bool negativeEta = false;
    if (t.first < 0)
      negativeEta = true;
    uint32_t crate = g.getCrate(t.first, t.second);
    uint32_t card = g.getCard(t.first, t.second);
    uint32_t region = g.getRegion(absCaloEta, absCaloPhi);
    UCTRegion* test = new UCTRegion(crate, card, negativeEta, region, fwVersion);
    test->setRegionSummary(i.raw());
    inputRegions.push_back(test);
    //This *should* fill the tensor in the proper order to be fed to the anomaly model
    //We take 4 off of the GCT eta/iEta.
    //iEta taken from this ranges from 4-17, (I assume reserving lower and higher for forward regions)
    //So our first index, index 0, is technically iEta=4, and so-on.
    //CICADA reads this as a flat vector


    //IMPORTANT CHANGE
    modelInput[14 * i.gctPhi() + (i.gctEta() - 4)] = i.et();
    inputTensorMapped(0, i.gctPhi(), i.gctEta() - 4, 0) = i.et();



  }
  // Check if we're using test patterns. If so, we overwrite the inputs with a test pattern
  if (overwriteWithTestPatterns) {
    unsigned int evt = iEvent.id().event();
    unsigned int totalTestPatterns = testPatterns.size();
    unsigned int patternElement = evt % totalTestPatterns;
    const edm::ParameterSet& element = testPatterns.at(patternElement);
    std::stringstream inputStream;
    std::string PhiRowString;

    edm::LogWarning("L1TCaloSummary") << "Overwriting existing CICADA input with test pattern!\n";

    for (unsigned short int iPhi = 1; iPhi <= 18; ++iPhi) {
      PhiRowString = "";
      std::stringstream PhiRowStringStream;
      PhiRowStringStream << "iPhi_" << iPhi;
      PhiRowString = PhiRowStringStream.str();
      std::vector<unsigned int> phiRow = element.getParameter<std::vector<unsigned int>>(PhiRowString);
      for (unsigned short int iEta = 1; iEta <= 14; ++iEta) {
        modelInput[14 * (iPhi - 1) + (iEta - 1)] = phiRow.at(iEta - 1);
        inputStream << phiRow.at(iEta - 1) << " ";
      }
      inputStream << "\n";
    }
    edm::LogInfo("L1TCaloSummary") << "Input Stream:\n" << inputStream.str();
  }

  //Extract model output
  OUTPUT modelResult[1] = {
      OUTPUT("0.0", 10)};  //the 10 here refers to the fact that we read in "0.0" as a decimal number
  model->prepare_input(modelInput);
  model->predict();
  model->read_result(modelResult);

  CICADAScore->push_back(0, modelResult[0].to_float());

  if (overwriteWithTestPatterns)
    edm::LogInfo("L1TCaloSummary") << "Test Pattern Output: " << CICADAScore->at(0, 0);

  summaryCard.setRegionData(inputRegions);

  if (!summaryCard.process()) {
    edm::LogError("L1TCaloSummary") << "UCT: Failed to process summary card" << std::endl;
    exit(1);
  }

  double pt = 0;
  double eta = -999.;
  double phi = -999.;
  double mass = 0.0;

  auto input_tensor_name = metaGraph->signature_def().at("serving_default").inputs().at("INPUT").name();
  auto phi_output_name = metaGraph->signature_def().at("serving_default").outputs().at("phi_output").name();
  auto eta_output_name = metaGraph->signature_def().at("serving_default").outputs().at("eta_output").name();
  auto reconstructed_output_name = metaGraph->signature_def().at("serving_default").outputs().at("reconstructed_output").name();

// edm::LogError("L1TCaloSummary") << "Input Tensor Name: " << input_tensor_name;
// edm::LogError("L1TCaloSummary") << "Phi Output Name: " << phi_output_name;
// edm::LogError("L1TCaloSummary") << "Eta Output Name: " << eta_output_name;


  std::vector<tensorflow::Tensor> outputs;

  // tensorflow::run(session, {{"serving_default_INPUT", input}}, {"StatefulPartitionedCall:0"}, &outputs);
  tensorflow::run(session, 
                {{"serving_default_INPUT", input}},  
                {"StatefulPartitionedCall:0","StatefulPartitionedCall:1"},  
                &outputs);
// if (outputs.size() == 2) {
//      edm::LogError("L1TCaloSummary") << "Error: Expected 2 outputs, but got " << outputs.size() << std::endl;
//     edm::LogError("L1TCaloSummary") << "output0: " << outputs[0];
//     edm::LogError("L1TCaloSummary") << "outputMappedEta: " << outputs[1];
//     return;
// }
  auto outputMappedEta = outputs[0].flat<float>();  // Now use flat<float> correctly for a 1D tensor
  auto outputMappedPhi = outputs[1].flat<float>();  // Now use flat<float> correctly for a 1D tensor

//   edm::LogError("L1TCaloSummary") << "outputMappedEta: " << outputMappedEta;
// edm::LogError("L1TCaloSummary") << "outputMappedEta(0): " << outputMappedEta(0);
// edm::LogError("L1TCaloSummary") << "outputMappedPhi: " << outputMappedPhi;
// edm::LogError("L1TCaloSummary") << "outputMappedPhi(0): " << outputMappedPhi(0);


  float predictedPhi = outputMappedPhi(0);  // First value in the result (phi)
  float predictedEta = outputMappedEta(0);  // Second value in the result (eta)
// edm::LogError("L1TCaloSummary") << "outputMappedPHI: " << predictedPhi;
// edm::LogError("L1TCaloSummary") << "outputMappedETA: " << predictedEta;


  std::array<double, 2> lowerBound = {{6.508498, 6.5006247}};
  std::array<double, 2> upperBound = {{7.4892583, 7.476416}};


  double threshold = 0.2;

  bool zeroEnergy = (predictedPhi <= lowerBound[0] + threshold && predictedPhi >= upperBound[0] - threshold &&
                    predictedEta <= lowerBound[1] + threshold && predictedEta >= upperBound[1] - threshold); 

  if(zeroEnergy){
    predictedEta=-99;
    predictedPhi=-99;
  }
                    
  if (!zeroEnergy) {
      for (int dPhi = -1; dPhi <= 1; ++dPhi) {
          for (int dEta = -1; dEta <= 1; ++dEta) {
              int phiIndex = predictedPhi + dPhi;
              int etaIndex = predictedEta + dEta;
                     
              if (phiIndex >= 0 && phiIndex < 18 && etaIndex >= 0 && etaIndex < 14) {
                  // flattened index in modelInputB
                  int flatIndex = phiIndex * 14 + etaIndex;
                  
                  // Accumulate the energy
                  pt += static_cast<double>(modelInput[flatIndex]);
              }
          }
      }
  }
  uint32_t iPredictedEta = static_cast<int>(round(predictedEta));
  uint32_t iPredictedPhi = static_cast<int>(round(predictedPhi));

  for (const L1CaloRegion& i : *regionCollection) {
      if ((i.gctEta() - 4 == iPredictedEta) && (i.gctPhi() == iPredictedPhi)) {
  	  //std::cout<<"HERE!!!!!! "<<i.gctEta()<<"\t"<<i.gctPhi()<<std::endl;
  
          UCTRegionIndex r = g.getUCTRegionIndexFromL1CaloRegion(i.gctEta(), i.gctPhi());
          UCTTowerIndex t = g.getUCTTowerIndexFromL1CaloRegion(r, i.raw());
          
          uint32_t absCaloEta = std::abs(t.first);
          uint32_t absCaloPhi = std::abs(t.second);
          bool negativeEta = (t.first < 0);
  
          eta = g.getUCTTowerEta(absCaloEta);
          if (negativeEta)
              eta = -eta;
          phi = g.getUCTTowerPhi(absCaloPhi);
  
          bJetCands->push_back(L1JetParticle(math::PtEtaPhiMLorentzVector(pt, eta, phi, mass), L1JetParticle::kCentral));
  
          break;
      }
  }

  iEvent.put(std::move(bJetCands), "Boosted");
  //Write out anomaly score
  iEvent.put(std::move(CICADAScore), "CICADAScore");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template <class INPUT, class OUTPUT>
void L1TCaloSummary<INPUT, OUTPUT>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// Initial version, X.0.0, input/output typing
typedef L1TCaloSummary<ap_ufixed<10, 10>, ap_fixed<11, 5>> L1TCaloSummary_CICADA_v1p0p0;
typedef L1TCaloSummary<ap_uint<10>, ap_ufixed<16, 8>> L1TCaloSummary_CICADA_v2p0p0;
DEFINE_FWK_MODULE(L1TCaloSummary_CICADA_v1p0p0);
DEFINE_FWK_MODULE(L1TCaloSummary_CICADA_v2p0p0);
// X.1.0 version input.output typing
typedef L1TCaloSummary<ap_ufixed<10, 10>, ap_fixed<11, 5>> L1TCaloSummary_CICADA_v1p1p0;
typedef L1TCaloSummary<ap_uint<10>, ap_ufixed<16, 8>> L1TCaloSummary_CICADA_v2p1p0;
DEFINE_FWK_MODULE(L1TCaloSummary_CICADA_v1p1p0);
DEFINE_FWK_MODULE(L1TCaloSummary_CICADA_v2p1p0);
// X.1.1 version input/output typing
typedef L1TCaloSummary<ap_uint<10>, ap_ufixed<16, 8, AP_RND, AP_SAT, AP_SAT>> L1TCaloSummary_CICADA_vXp1p1;
DEFINE_FWK_MODULE(L1TCaloSummary_CICADA_vXp1p1);
// X.1.2 version input/output typing
typedef L1TCaloSummary<ap_uint<10>, ap_ufixed<16, 8, AP_RND_CONV, AP_SAT>> L1TCaloSummary_CICADA_vXp1p2;
DEFINE_FWK_MODULE(L1TCaloSummary_CICADA_vXp1p2);
