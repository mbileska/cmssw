/*
 * Phase-2 HF/HGCAL (MHH) mixer emulator.
 */

#include <ap_int.h>

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "L1Trigger/L1CaloTrigger/interface/MHH_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_MHH_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_MHH_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/MHH_cpp.h"

class Phase2L1MHHEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1MHHEmulator(const edm::ParameterSet&);
  ~Phase2L1MHHEmulator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  std::array<edm::EDGetTokenT<std::vector<uint64_t>>, p2mhh::MHH_N_INPUT_LINKS> inputLinkTokens_;
  bool debug_;
};

Phase2L1MHHEmulator::Phase2L1MHHEmulator(const edm::ParameterSet& config)
    : debug_(config.getParameter<bool>("debug")) {
  const auto inputLinks = config.getParameter<std::vector<edm::InputTag>>("inputLinks");
  if (inputLinks.size() != p2mhh::MHH_N_INPUT_LINKS) {
    throw cms::Exception("Phase2L1MHHEmulator")
        << "Expected exactly " << p2mhh::MHH_N_INPUT_LINKS << " MHH input links";
  }

  for (unsigned int i = 0; i < inputLinkTokens_.size(); ++i) {
    inputLinkTokens_[i] = consumes<std::vector<uint64_t>>(inputLinks[i]);
  }
  for (int i = 0; i < p2mhh::MHH_N_OUTPUT_LINKS; ++i) {
    produces<std::vector<uint64_t>>(std::string("LinkOut") + std::to_string(i));
  }
}

void Phase2L1MHHEmulator::produce(edm::Event& event, const edm::EventSetup&) {
  std::array<ap_uint<576>, p2mhh::MHH_N_INPUT_LINKS> linkIn{};
  std::array<ap_uint<576>, p2mhh::MHH_N_OUTPUT_LINKS> linkOut{};

  for (unsigned int i = 0; i < inputLinkTokens_.size(); ++i) {
    edm::Handle<std::vector<uint64_t>> words;
    event.getByToken(inputLinkTokens_[i], words);
    if (!words.isValid()) {
      throw cms::Exception("Phase2L1MHHEmulator") << "MHH input link " << i << " is missing";
    }
    if (words->size() != 9) {
      throw cms::Exception("Phase2L1MHHEmulator")
          << "MHH input link " << i << " has " << words->size() << " words; expected 9";
    }
    for (unsigned int word = 0; word < 9; ++word) {
      linkIn[i].range(word * 64 + 63, word * 64) = words->at(word);
    }
  }

  p2mhh::mhh_algo_top(linkIn.data(), linkOut.data());

  for (int i = 0; i < p2mhh::MHH_N_OUTPUT_LINKS; ++i) {
    auto words = std::make_unique<std::vector<uint64_t>>();
    words->reserve(9);
    for (unsigned int word = 0; word < 9; ++word) {
      words->push_back(linkOut[i].range(word * 64 + 63, word * 64).to_uint64());
    }
    if (debug_) {
      edm::LogVerbatim("Phase2L1MHHEmulator") << "Output link " << i << " contains 9 words";
    }
    event.put(std::move(words), std::string("LinkOut") + std::to_string(i));
  }
}

void Phase2L1MHHEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<std::vector<edm::InputTag>>("inputLinks", {});
  description.add<bool>("debug", false);
  descriptions.add("phase2L1MHHEmulator", description);
}

DEFINE_FWK_MODULE(Phase2L1MHHEmulator);
