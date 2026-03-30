/* 
 * Description: Phase 2 GCT SumCard emulator
 * Author: Mila Bileska
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <cstdint>

// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1CaloTrigger/interface/GCTSum_h.h"
#include "L1Trigger/L1CaloTrigger/interface/GCTSum_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/GCTSumToGT_h.h"
#include "L1Trigger/L1CaloTrigger/interface/GCTSumToGT_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_GCT_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_GCT_cpp.h"

////////////////////////////////////////////////////////////////////////////////

// Declare the Phase2L1GCTSumEmulator class and its methods

class Phase2L1GCTSumEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1GCTSumEmulator(const edm::ParameterSet&);
  ~Phase2L1GCTSumEmulator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  std::array<edm::EDGetTokenT<std::vector<uint64_t> >, 24> inputLinkTokens_;
  bool debug_;
};

//////////////////////////////////////////////////////////////////

Phase2L1GCTSumEmulator::Phase2L1GCTSumEmulator(const edm::ParameterSet& iConfig)
    : debug_(iConfig.getParameter<bool>("debug")) {
  const std::vector<edm::InputTag> inputLinks = iConfig.getParameter<std::vector<edm::InputTag> >("inputLinks");
  if (inputLinks.size() != 24) {
    throw cms::Exception("Phase2L1GCTSumEmulator") << "Expected exactly 24 input links (12 positive eta + 12 negative eta)";
  }

  for (unsigned int i = 0; i < 24; ++i) {
    inputLinkTokens_[i] = consumes<std::vector<uint64_t> >(inputLinks[i]);
  }

  produces<std::vector<uint64_t> >("LinkOut0");
  produces<std::vector<uint64_t> >("LinkOut1");
  produces<std::vector<uint64_t> >("LinkOut2");
  produces<std::vector<uint64_t> >("LinkOut3");
  produces<std::vector<uint64_t> >("LinkOut4");
  produces<std::vector<uint64_t> >("LinkOut5");
}

void Phase2L1GCTSumEmulator::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  (void)iSetup;

  std::array<ap_uint<576>, 12> link_in_pos;
  std::array<ap_uint<576>, 12> link_in_neg;
  std::array<ap_uint<576>, 3> link_out_pos;
  std::array<ap_uint<576>, 3> link_out_neg;
  std::array<ap_uint<576>, 6> link_in_gt;
  std::array<ap_uint<576>, 6> link_out_gt;

  for (unsigned int i = 0; i < 12; ++i) {
    edm::Handle<std::vector<uint64_t> > handle;
    iEvent.getByToken(inputLinkTokens_[i], handle);

    ap_uint<576> packed = 0;
    if (handle.isValid()) {
      for (unsigned int word = 0; word < handle->size() && word < 9; ++word) {
        packed.range((word * 64) + 63, word * 64) = (*handle)[word];
      }
    }
    link_in_pos[i] = packed;
  }

  for (unsigned int i = 0; i < 12; ++i) {
    edm::Handle<std::vector<uint64_t> > handle;
    iEvent.getByToken(inputLinkTokens_[12 + i], handle);

    ap_uint<576> packed = 0;
    if (handle.isValid()) {
      for (unsigned int word = 0; word < handle->size() && word < 9; ++word) {
        packed.range((word * 64) + 63, word * 64) = (*handle)[word];
      }
    }
    link_in_neg[i] = packed;
  }

  p2gctsum::algo_top(link_in_pos.data(), link_out_pos.data());
  p2gctsum::algo_top(link_in_neg.data(), link_out_neg.data());

  link_in_gt[0] = link_out_pos[0];
  link_in_gt[1] = link_out_pos[1];
  link_in_gt[2] = link_out_pos[2];
  link_in_gt[3] = link_out_neg[0];
  link_in_gt[4] = link_out_neg[1];
  link_in_gt[5] = link_out_neg[2];

  p2gctsumGT::algo_top_GT(link_in_gt.data(), link_out_gt.data());

  for (unsigned int i = 0; i < 6; ++i) {
    std::unique_ptr<std::vector<uint64_t> > outWords = std::make_unique<std::vector<uint64_t> >();
    outWords->reserve(9);

    for (unsigned int word = 0; word < 9; ++word) {
      outWords->push_back((uint64_t)link_out_gt[i].range((word * 64) + 63, word * 64));
    }

    if (debug_) {
      edm::LogVerbatim("Phase2L1GCTSumEmulator") << "Output Link " << i << " has " << outWords->size() << " 64b words";
    }

    iEvent.put(std::move(outWords), std::string("LinkOut") + std::to_string(i));
  }
}

void Phase2L1GCTSumEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<edm::InputTag> >("inputLinks", std::vector<edm::InputTag>());
  desc.add<bool>("debug", false);
  descriptions.add("phase2L1GCTSumEmulator", desc);
}

DEFINE_FWK_MODULE(Phase2L1GCTSumEmulator);
