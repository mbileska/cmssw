#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "ap_fixed.h"
#include "hls4ml/emulator.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "DataFormats/L1CaloTrigger/interface/CICADA.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1TGlobal/interface/AXOL1TLScore.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace {
  constexpr unsigned int kAXOInputs = 57;
  constexpr unsigned int kCICADAInputs = 252;

  template <typename ResultType, typename LossType>
  LossType readAXOResult(hls4mlEmulator::Model& model) {
    std::pair<ResultType, LossType> result;
    model.read_result(&result);
    return result.second;
  }

  template <typename Collection>
  float scoreAtBxZero(const Collection& scores, const char* product) {
    if (scores.getFirstBX() > 0 || scores.getLastBX() < 0 || scores.isEmpty(0)) {
      throw cms::Exception("ProductNotFound") << product << " has no score at BX=0";
    }
    return scores.at(0, 0).getAXOScore();
  }

  template <>
  float scoreAtBxZero<l1t::CICADABxCollection>(const l1t::CICADABxCollection& scores, const char* product) {
    if (scores.getFirstBX() > 0 || scores.getLastBX() < 0 || scores.isEmpty(0)) {
      throw cms::Exception("ProductNotFound") << product << " has no score at BX=0";
    }
    return scores.at(0, 0);
  }
}  // namespace

class L1AnomalyScoreTableProducer : public edm::stream::EDProducer<edm::stream::WatchRuns> {
public:
  explicit L1AnomalyScoreTableProducer(const edm::ParameterSet& config)
      : axoScoreToken_(consumes<AXOL1TLScoreBxCollection>(config.getParameter<edm::InputTag>("axoScore"))),
        cicadaScoreToken_(consumes<l1t::CICADABxCollection>(config.getParameter<edm::InputTag>("cicadaScore"))),
        muonToken_(consumes<l1t::MuonBxCollection>(config.getParameter<edm::InputTag>("muons"))),
        egammaToken_(consumes<l1t::EGammaBxCollection>(config.getParameter<edm::InputTag>("egammas"))),
        jetToken_(consumes<l1t::JetBxCollection>(config.getParameter<edm::InputTag>("jets"))),
        etsumToken_(consumes<l1t::EtSumBxCollection>(config.getParameter<edm::InputTag>("etsums"))),
        regionToken_(consumes<L1CaloRegionCollection>(config.getParameter<edm::InputTag>("regions"))),
        menuToken_(esConsumes<edm::Transition::BeginRun>()),
        cicadaModelName_(config.getParameter<std::string>("cicadaModel")),
        cicadaLoader_(cicadaModelName_),
        maxScoreDifference_(config.getParameter<double>("maxScoreDifference")) {
    cicadaModel_ = cicadaLoader_.load_model();
    produces<nanoaod::FlatTable>("AXO");
    produces<nanoaod::FlatTable>("CICADA");
    produces<nanoaod::FlatTable>("AXOInput");
    produces<nanoaod::FlatTable>("CICADAInput");
  }

  void beginRun(edm::Run const&, edm::EventSetup const& setup) override {
    const auto& menu = setup.getData(menuToken_);
    std::set<std::string> versions;
    std::set<int> relativeBxs;
    for (const auto& conditionEntry : menu.getConditionMap()) {
      for (const auto& object : conditionEntry.second.getObjects()) {
        if (object.getType() != tmeventsetup::Axol1tl)
          continue;
        relativeBxs.insert(object.getBxOffset());
        for (const auto& cut : object.getCuts()) {
          if (cut.getCutType() == tmeventsetup::Model && !cut.getData().empty())
            versions.insert(cut.getData());
        }
      }
    }
    if (versions.size() != 1) {
      throw cms::Exception("Configuration") << "Expected exactly one AXOL1TL model in L1 menu '" << menu.getName()
                                             << "', found " << versions.size();
    }
    if (relativeBxs.size() != 1 || *relativeBxs.begin() != 0) {
      throw cms::Exception("Configuration")
          << "Independent AXOL1TL closure currently requires every menu condition at relative BX=0";
    }

    axoModelVersion_ = *versions.begin();
    axoModelName_ = "GTADModel_" + axoModelVersion_;
    if (axoModelName_ != "GTADModel_v3" && axoModelName_ != "GTADModel_v4" && axoModelName_ != "GTADModel_v5") {
      throw cms::Exception("Configuration") << "Unsupported AXOL1TL model output type for '" << axoModelName_ << "'";
    }
    axoLoader_.reset(axoModelName_);
    axoModel_ = axoLoader_.load_model();
    edm::LogInfo("L1AnomalyScoreTableProducer") << "Independent inference uses " << axoModelName_ << " and "
                                                 << cicadaModelName_;
  }

  void produce(edm::Event& event, const edm::EventSetup&) override {
    using AXOInput = ap_fixed<18, 13>;
    using AXOLoss = ap_ufixed<18, 14>;
    using CICADAInput = ap_uint<10>;
    using CICADAOutput = ap_ufixed<16, 8, AP_RND_CONV, AP_SAT>;

    const auto& officialAXOScores = event.get(axoScoreToken_);
    const auto& officialCICADAScores = event.get(cicadaScoreToken_);
    const float officialAXO = scoreAtBxZero(officialAXOScores, "AXOL1TL score collection");
    const float officialCICADA = scoreAtBxZero(officialCICADAScores, "CICADA score collection");

    std::array<AXOInput, kAXOInputs> axoInput{};
    std::vector<int16_t> axoValue2(kAXOInputs, 0);
    fillAXOInput(event, axoInput, axoValue2);

    axoModel_->prepare_input(axoInput.data());
    axoModel_->predict();
    AXOLoss loss;
    if (axoModelName_ == "GTADModel_v3" || axoModelName_ == "GTADModel_v4") {
      using AXOResult = std::array<ap_fixed<10, 7, AP_RND_CONV, AP_SAT>, 8>;
      loss = readAXOResult<AXOResult, AXOLoss>(*axoModel_);
    } else {
      using AXOResult = ap_fixed<18, 14, AP_RND_CONV, AP_SAT>;
      loss = readAXOResult<AXOResult, AXOLoss>(*axoModel_);
    }
    const float recomputedAXO = loss.to_float() * 16.f;

    std::array<CICADAInput, kCICADAInputs> cicadaInput{};
    std::vector<uint16_t> cicadaEt(kCICADAInputs, 0);
    std::vector<uint8_t> cicadaPresent(kCICADAInputs, 0);
    std::vector<uint8_t> cicadaGctEta(kCICADAInputs, 0);
    std::vector<uint8_t> cicadaGctPhi(kCICADAInputs, 0);
    fillCICADAInput(event, cicadaInput, cicadaEt, cicadaPresent, cicadaGctEta, cicadaGctPhi);

    cicadaModel_->prepare_input(cicadaInput.data());
    cicadaModel_->predict();
    std::array<CICADAOutput, 1> cicadaResult{};
    cicadaModel_->read_result(cicadaResult.data());
    const float recomputedCICADA = cicadaResult[0].to_float();

    checkDifference("AXOL1TL", officialAXO, recomputedAXO);
    checkDifference("CICADA", officialCICADA, recomputedCICADA);

    auto axoTable = std::make_unique<nanoaod::FlatTable>(1, "L1AXO", true);
    axoTable->addColumnValue<float>("score", officialAXO, "AXOL1TL score produced inside the uGT emulator at BX=0");
    axoTable->addColumnValue<float>(
        "scoreRecomputed", recomputedAXO, "AXOL1TL score independently recomputed from the saved 57-value tensor");
    axoTable->addColumnValue<float>(
        "scoreDifference", recomputedAXO - officialAXO, "recomputed minus uGT AXOL1TL score");

    auto cicadaTable = std::make_unique<nanoaod::FlatTable>(1, "L1CICADA", true);
    cicadaTable->addColumnValue<float>(
        "score", officialCICADA, "CICADA score produced inside the Layer-1 summary emulator at BX=0");
    cicadaTable->addColumnValue<float>(
        "scoreRecomputed", recomputedCICADA, "CICADA score independently recomputed from the saved 252-region tensor");
    cicadaTable->addColumnValue<float>(
        "scoreDifference", recomputedCICADA - officialCICADA, "recomputed minus Layer-1 CICADA score");

    auto axoInputTable = std::make_unique<nanoaod::FlatTable>(kAXOInputs, "L1AXOInput", false);
    axoInputTable->addColumn<int16_t>(
        "value2", axoValue2, "twice the exact AXOL1TL input; divide by two before inference");

    auto cicadaInputTable = std::make_unique<nanoaod::FlatTable>(kCICADAInputs, "L1CICADAInput", false);
    cicadaInputTable->addColumn<uint16_t>("et", cicadaEt, "10-bit CICADA region ET, ordered as 14*gctPhi+gctEta-4");
    cicadaInputTable->addColumn<uint8_t>("gctEta", cicadaGctEta, "GCT eta index (4..17)");
    cicadaInputTable->addColumn<uint8_t>("gctPhi", cicadaGctPhi, "GCT phi index (0..17)");
    cicadaInputTable->addColumn<uint8_t>(
        "present", cicadaPresent, "one when the input region collection contained this eta-phi cell");

    event.put(std::move(axoTable), "AXO");
    event.put(std::move(cicadaTable), "CICADA");
    event.put(std::move(axoInputTable), "AXOInput");
    event.put(std::move(cicadaInputTable), "CICADAInput");
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("axoScore", edm::InputTag("simGtStage2Digis", "AXOScore"));
    desc.add<edm::InputTag>("cicadaScore", edm::InputTag("simCaloStage2Layer1Summary", "CICADAScore"));
    desc.add<edm::InputTag>("muons", edm::InputTag("simGmtStage2Digis"));
    desc.add<edm::InputTag>("egammas", edm::InputTag("simCaloStage2Digis"));
    desc.add<edm::InputTag>("jets", edm::InputTag("simCaloStage2Digis"));
    desc.add<edm::InputTag>("etsums", edm::InputTag("simCaloStage2Digis"));
    desc.add<edm::InputTag>("regions", edm::InputTag("simCaloStage2Layer1Digis"));
    desc.add<std::string>("cicadaModel", "CICADAModel_v2p2p0");
    desc.add<double>("maxScoreDifference", 1.e-6)
        ->setComment("negative disables the event-level agreement check");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  template <typename AXOInput>
  void fillAXOInput(const edm::Event& event,
                    std::array<AXOInput, kAXOInputs>& input,
                    std::vector<int16_t>& value2) const {
    const auto& etsums = event.get(etsumToken_);
    if (etsums.getFirstBX() <= 0 && etsums.getLastBX() >= 0) {
      for (auto it = etsums.begin(0); it != etsums.end(0); ++it) {
        if (it->getType() == l1t::EtSum::EtSumType::kMissingEt) {
          input[0] = it->hwPt() * .5;
          input[2] = it->hwPhi();
          value2[0] = it->hwPt();
          value2[2] = 2 * it->hwPhi();
        }
      }
    }

    const auto& egammas = event.get(egammaToken_);
    if (egammas.getFirstBX() <= 0 && egammas.getLastBX() >= 0) {
      unsigned int object = 0;
      for (auto it = egammas.begin(0); it != egammas.end(0) && object < 4; ++it, ++object) {
        const unsigned int offset = 3 + 3 * object;
        input[offset] = it->hwPt() * .5;
        input[offset + 1] = it->hwEta();
        input[offset + 2] = it->hwPhi();
        value2[offset] = it->hwPt();
        value2[offset + 1] = 2 * it->hwEta();
        value2[offset + 2] = 2 * it->hwPhi();
      }
    }

    const auto& muons = event.get(muonToken_);
    if (muons.getFirstBX() <= 0 && muons.getLastBX() >= 0) {
      unsigned int object = 0;
      for (auto it = muons.begin(0); it != muons.end(0) && object < 4; ++it, ++object) {
        const unsigned int offset = 15 + 3 * object;
        input[offset] = it->hwPt() * .5;
        input[offset + 1] = it->hwEtaAtVtx();
        input[offset + 2] = it->hwPhiAtVtx();
        value2[offset] = it->hwPt();
        value2[offset + 1] = 2 * it->hwEtaAtVtx();
        value2[offset + 2] = 2 * it->hwPhiAtVtx();
      }
    }

    const auto& jets = event.get(jetToken_);
    if (jets.getFirstBX() <= 0 && jets.getLastBX() >= 0) {
      unsigned int object = 0;
      for (auto it = jets.begin(0); it != jets.end(0) && object < 10; ++it, ++object) {
        const unsigned int offset = 27 + 3 * object;
        input[offset] = it->hwPt() * .5;
        input[offset + 1] = it->hwEta();
        input[offset + 2] = it->hwPhi();
        value2[offset] = it->hwPt();
        value2[offset + 1] = 2 * it->hwEta();
        value2[offset + 2] = 2 * it->hwPhi();
      }
    }
  }

  template <typename CICADAInput>
  void fillCICADAInput(const edm::Event& event,
                       std::array<CICADAInput, kCICADAInputs>& input,
                       std::vector<uint16_t>& et,
                       std::vector<uint8_t>& present,
                       std::vector<uint8_t>& gctEta,
                       std::vector<uint8_t>& gctPhi) const {
    for (unsigned int phi = 0; phi < 18; ++phi) {
      for (unsigned int eta = 4; eta < 18; ++eta) {
        const unsigned int index = 14 * phi + eta - 4;
        gctEta[index] = eta;
        gctPhi[index] = phi;
      }
    }

    const auto& regions = event.get(regionToken_);
    for (const auto& region : regions) {
      if (region.gctPhi() >= 18 || region.gctEta() < 4 || region.gctEta() >= 18) {
        throw cms::Exception("LogicError") << "CICADA region outside expected grid: gctEta=" << region.gctEta()
                                           << ", gctPhi=" << region.gctPhi();
      }
      if (region.et() > 1023) {
        throw cms::Exception("LogicError") << "CICADA region ET exceeds its 10-bit model input: " << region.et();
      }
      const unsigned int index = 14 * region.gctPhi() + region.gctEta() - 4;
      if (present[index]) {
        throw cms::Exception("LogicError") << "Duplicate CICADA region at gctEta=" << region.gctEta()
                                           << ", gctPhi=" << region.gctPhi();
      }
      input[index] = region.et();
      et[index] = region.et();
      present[index] = 1;
    }
  }

  void checkDifference(const char* name, float official, float recomputed) const {
    if (maxScoreDifference_ >= 0. && std::abs(recomputed - official) > maxScoreDifference_) {
      throw cms::Exception("LogicError") << name << " independent score does not match official emulation: official="
                                         << official << ", recomputed=" << recomputed
                                         << ", allowed absolute difference=" << maxScoreDifference_;
    }
  }

  const edm::EDGetTokenT<AXOL1TLScoreBxCollection> axoScoreToken_;
  const edm::EDGetTokenT<l1t::CICADABxCollection> cicadaScoreToken_;
  const edm::EDGetTokenT<l1t::MuonBxCollection> muonToken_;
  const edm::EDGetTokenT<l1t::EGammaBxCollection> egammaToken_;
  const edm::EDGetTokenT<l1t::JetBxCollection> jetToken_;
  const edm::EDGetTokenT<l1t::EtSumBxCollection> etsumToken_;
  const edm::EDGetTokenT<L1CaloRegionCollection> regionToken_;
  const edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> menuToken_;
  std::string axoModelVersion_;
  std::string axoModelName_;
  const std::string cicadaModelName_;
  hls4mlEmulator::ModelLoader axoLoader_;
  hls4mlEmulator::ModelLoader cicadaLoader_;
  std::shared_ptr<hls4mlEmulator::Model> axoModel_;
  std::shared_ptr<hls4mlEmulator::Model> cicadaModel_;
  const double maxScoreDifference_;
};

DEFINE_FWK_MODULE(L1AnomalyScoreTableProducer);
