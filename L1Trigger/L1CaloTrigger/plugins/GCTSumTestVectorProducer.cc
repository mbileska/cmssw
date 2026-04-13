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

  std::vector<uint64_t> unpack576ToWords(ap_uint<576> data) const;

  ap_uint<48> packEG(ap_uint<12> et, ap_uint<7> eta, ap_uint<9> phi, ap_uint<1> isBarrel) const;
  ap_uint<48> packHadron(ap_uint<12> et, ap_uint<6> eta, ap_uint<9> phi, ap_uint<4> seed, ap_uint<1> isBarrel) const;
  ap_uint<48> packSum(ap_uint<12> ex, ap_uint<12> ey, ap_uint<12> ht, ap_uint<12> spare = 0) const;

  void fillValidationPattern(std::array<ap_uint<576>, 24>& links, unsigned long long eventNumber) const;

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

ap_uint<48> GCTSumTestVectorProducer::packEG(ap_uint<12> et, ap_uint<7> eta, ap_uint<9> phi, ap_uint<1> isBarrel) const {
  ap_uint<48> out = 0;
  out.range(11, 0) = et;
  out.range(18, 12) = eta;
  out.range(27, 19) = phi;
  out.range(47, 47) = isBarrel;
  return out;
}

ap_uint<48> GCTSumTestVectorProducer::packHadron(ap_uint<12> et,
                                                 ap_uint<6> eta,
                                                 ap_uint<9> phi,
                                                 ap_uint<4> seed,
                                                 ap_uint<1> isBarrel) const {
  ap_uint<48> out = 0;
  out.range(11, 0) = et;
  out.range(17, 12) = eta;
  out.range(26, 18) = phi;
  out.range(30, 27) = seed;
  out.range(47, 47) = isBarrel;
  return out;
}

ap_uint<48> GCTSumTestVectorProducer::packSum(ap_uint<12> ex,
                                              ap_uint<12> ey,
                                              ap_uint<12> ht,
                                              ap_uint<12> spare) const {
  ap_uint<48> out = 0;
  out.range(11, 0) = ex;
  out.range(23, 12) = ey;
  out.range(35, 24) = ht;
  out.range(47, 36) = spare;
  return out;
}

void GCTSumTestVectorProducer::fillValidationPattern(std::array<ap_uint<576>, 24>& links,
                                                     unsigned long long eventNumber) const {
  for (auto& link : links)
    link = 0;

  // 10-pattern validation cycle, repeated by cmsRun over 20 events
  const unsigned int pattern = (unsigned int)((eventNumber - 1ULL) % 10ULL);

  auto putWord48 = [&](unsigned int linkIdx, unsigned int slot, ap_uint<48> word) {
    const ap_uint<10> start = slot * 48;
    const ap_uint<10> end = start + 47;
    links[linkIdx].range(end, start) = word;
  };

  switch (pattern) {
    case 0:
      // all-zero input
      break;

    case 1:
      // single positive-side EG-like object
      putWord48(0, 0, packEG(50, 10, 20, 0));
      break;

    case 2:
      // positive-side EG plus positive-side sum
      putWord48(0, 0, packEG(40, 11, 30, 0));
      putWord48(2, 0, packSum(5, 7, 9, 0));
      break;

    case 3:
      // dense positive-side ordering for EG and EGiso
      putWord48(0, 0, packEG(10, 10, 15, 0));
      putWord48(0, 1, packEG(60, 10, 16, 0));
      putWord48(3, 0, packEG(30, 10, 17, 1));
      putWord48(6, 0, packEG(50, 10, 18, 1));
      putWord48(9, 0, packEG(20, 10, 19, 1));
      putWord48(9, 1, packEG(40, 10, 20, 1));

      putWord48(0, 6, packEG(12, 10, 25, 0));
      putWord48(3, 6, packEG(32, 10, 26, 1));
      putWord48(6, 6, packEG(52, 10, 27, 1));
      putWord48(9, 6, packEG(22, 10, 28, 1));
      putWord48(9, 7, packEG(42, 10, 29, 1));
      break;

    case 4:
      // positive and negative eta both populated
      putWord48(0, 0, packEG(45, 10, 40, 0));                 // positive EG
      putWord48(1, 0, packHadron(55, 12, 45, 3, 1));          // positive jet
      putWord48(12, 0, packEG(35, 10, 140, 0));               // negative EG
      putWord48(13, 6, packHadron(25, 12, 145, 2, 1));        // negative tau
      break;

    case 5:
      // sparse hadron/tau plus sum
      putWord48(1, 0, packHadron(55, 9, 60, 4, 1));           // positive jet
      putWord48(1, 6, packHadron(35, 9, 61, 2, 1));           // positive tau
      putWord48(2, 0, packSum(4, 6, 8, 0));
      break;

    case 6:
      // stitching: endcap + barrel same phi should merge
      putWord48(0, 0, packEG(20, 31, 70, 0));                 // endcap-like source 0
      putWord48(3, 0, packEG(12, 0, 70, 1));                  // barrel-like source 1
      break;

    case 7:
      // stitching: dphi = +1 should merge
      putWord48(0, 0, packEG(18, 31, 80, 0));
      putWord48(3, 0, packEG(11, 0, 81, 1));
      break;

    case 8:
      // no stitch: dphi > 1
      putWord48(0, 0, packEG(20, 31, 90, 0));
      putWord48(3, 0, packEG(12, 0, 94, 1));
      break;

    case 9:
      // sum aggregation: positive + negative
      putWord48(2, 0, packSum(5, 7, 9, 0));                   // positive side source 0 sum
      putWord48(14, 0, packSum(11, 13, 17, 0));               // negative side source 0 sum
      break;
  }
}

void GCTSumTestVectorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  (void)iSetup;

  std::array<ap_uint<576>, 24> links;
  for (auto& link : links)
    link = 0;

  const unsigned long long eventNumber = iEvent.id().event();

  if (patternMode_ == "validation") {
    fillValidationPattern(links, eventNumber);
  } else {
    // default to validation
    fillValidationPattern(links, eventNumber);
  }

  for (unsigned int i = 0; i < 24; ++i) {
    auto outWords = std::make_unique<std::vector<uint64_t> >(unpack576ToWords(links[i]));
    iEvent.put(std::move(outWords), std::string("LinkIn") + std::to_string(i));
  }
}

void GCTSumTestVectorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("patternMode", "validation");
  desc.add<bool>("debug", false);
  descriptions.add("gctSumTestVectorProducer", desc);
}

DEFINE_FWK_MODULE(GCTSumTestVectorProducer);