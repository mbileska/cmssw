/* 
 * Description: Test-vector producer for Phase 2 GCT SumCard emulator
 * Provides deterministic link-level patterns, including Alexander Savin's
 * pT-ordering vector for one IP1 side.
 Author: Mila Bileska
 */

// system include files
#include <ap_int.h>
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

class GCTSumTestVectorProducer : public edm::stream::EDProducer<> {
public:
  explicit GCTSumTestVectorProducer(const edm::ParameterSet&);
  ~GCTSumTestVectorProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  void fillPtSort(std::array<ap_uint<576>, 24>& links) const;
  void fillCyclicPattern(std::array<ap_uint<576>, 24>& links, unsigned long long eventNumber) const;
  std::vector<uint64_t> unpack576ToWords(ap_uint<576> data) const;

  std::string patternMode_;
  bool debug_;
};

GCTSumTestVectorProducer::GCTSumTestVectorProducer(const edm::ParameterSet& iConfig)
    : patternMode_(iConfig.getParameter<std::string>("patternMode")),
      debug_(iConfig.getParameter<bool>("debug")) {
  for (unsigned int i = 0; i < 24; ++i) {
    produces<std::vector<uint64_t> >(std::string("LinkIn") + std::to_string(i));
  }
}

std::vector<uint64_t> GCTSumTestVectorProducer::unpack576ToWords(ap_uint<576> data) const {
  std::vector<uint64_t> out;
  out.reserve(9);
  for (unsigned int word = 0; word < 9; ++word) {
    out.push_back((uint64_t)data.range((word * 64) + 63, word * 64));
  }
  return out;
}

void GCTSumTestVectorProducer::fillPtSort(std::array<ap_uint<576>, 24>& links) const {
  for (auto& link : links)
    link = 0;

  // This reproduces Alexander Savin's standalone TV at the IP1 input:
  // - populate only one IP1 side (positive eta here: links 0..11)
  // - use simple, well defined 48-bit payload values
  // - test pT ordering only
  //
  // Per source:
  //   link A : EGs (slots 0..5), EGIs left zero
  //   link B : Jets (slots 0..5), Taus left zero
  //   link C : sums left zero

  for (unsigned int i = 0; i < 6; ++i) {
    ap_uint<10> start = i * 48;
    ap_uint<10> end = start + 47;

    // Source 0
    links[0].range(end, start) = (ap_uint<48>)(2 * i + 1);
    links[1].range(end, start) = (ap_uint<48>)(8 * i + 2);

    // Source 1
    links[3].range(end, start) = (ap_uint<48>)(32 * i + 3);
    links[4].range(end, start) = (ap_uint<48>)(32 * i + 4);

    // Source 2
    links[6].range(end, start) = (ap_uint<48>)(16 * i + 5);
    links[7].range(end, start) = (ap_uint<48>)(16 * i + 6);

    // Source 3
    links[9].range(end, start) = (ap_uint<48>)(32 * i + 7);
    links[10].range(end, start) = (ap_uint<48>)(64 * i + 8);
  }

  // links[2], links[5], links[8], links[11] remain zero (no sums)
  // links[12..23] remain zero (negative eta side empty)
}

void GCTSumTestVectorProducer::fillCyclicPattern(std::array<ap_uint<576>, 24>& links,
                                                 unsigned long long eventNumber) const {
  for (auto& link : links)
    link = 0;

  unsigned int pattern = (unsigned int)((eventNumber - 1ULL) % 6ULL);

  switch (pattern) {
    case 0:
      // all zero
      break;

    case 1:
      // single positive-side EG-like object
      links[0].range(47, 0) = (ap_uint<48>)0x321;
      break;

    case 2:
      // positive-side EG and positive-side sums
      links[0].range(47, 0) = (ap_uint<48>)0x6E1;
      links[2].range(47, 0) = (ap_uint<48>)0x141;
      links[2].range(95, 48) = (ap_uint<48>)0x041;
      break;

    case 3:
      // dense positive-side pattern
      for (unsigned int i = 0; i < 6; ++i) {
        ap_uint<10> start = i * 48;
        ap_uint<10> end = start + 47;
        links[0].range(end, start) = (ap_uint<48>)(0x4C1 + 0x20 * i);
        links[3].range(end, start) = (ap_uint<48>)(0x741 + 0x20 * i);
        links[6].range(end, start) = (ap_uint<48>)(0x9C1 + 0x20 * i);
        links[9].range(end, start) = (ap_uint<48>)(0xC41 + 0x20 * i);
      }
      for (unsigned int i = 0; i < 4; ++i) {
        ap_uint<10> start = i * 48;
        ap_uint<10> end = start + 47;
        links[2].range(end, start) = (ap_uint<48>)(0x281 + 0x20 * i);
      }
      break;

    case 4:
      // positive and negative eta both populated
      links[0].range(47, 0) = (ap_uint<48>)0x8C1;
      links[6].range(47, 0) = (ap_uint<48>)0xB41;
      links[14].range(47, 0) = (ap_uint<48>)0xBE1;
      links[18].range(47, 0) = (ap_uint<48>)0xF01;
      links[22].range(95, 48) = (ap_uint<48>)0x7A1;
      break;

    case 5:
      // sparse hadron/tau-like pattern and a sum
      links[7].range(47, 0) = (ap_uint<48>)0x961;
      links[10].range(47, 0) = (ap_uint<48>)0xF01;
      links[22].range(95, 48) = (ap_uint<48>)0x15C;
      break;
  }
}

void GCTSumTestVectorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  (void)iSetup;

  std::array<ap_uint<576>, 24> links;
  for (auto& link : links)
    link = 0;

  unsigned long long eventNumber = iEvent.id().event();

  if (patternMode_ == "ptsort") {
    fillPtSort(links);
  } else {
    fillCyclicPattern(links, eventNumber);
  }

  for (unsigned int i = 0; i < 24; ++i) {
    auto outWords = std::make_unique<std::vector<uint64_t> >(unpack576ToWords(links[i]));
    iEvent.put(std::move(outWords), std::string("LinkIn") + std::to_string(i));
  }
}

void GCTSumTestVectorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("patternMode", "cyclic");
  desc.add<bool>("debug", false);
  descriptions.add("gctSumTestVectorProducer", desc);
}

DEFINE_FWK_MODULE(GCTSumTestVectorProducer);
