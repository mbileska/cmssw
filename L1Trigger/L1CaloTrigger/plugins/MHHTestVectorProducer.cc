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
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

namespace {
constexpr unsigned int kInputLinks = 16;

ap_uint<48> packGamma(unsigned int et, unsigned int eta, unsigned int phi) {
  ap_uint<48> word = 0;
  word.range(11, 0) = et;
  word.range(18, 12) = eta;
  word.range(27, 19) = phi;
  return word;
}

ap_uint<48> packHadron(unsigned int et, unsigned int eta, unsigned int phi, unsigned int seed) {
  ap_uint<48> word = 0;
  word.range(11, 0) = et;
  word.range(17, 12) = eta;
  word.range(26, 18) = phi;
  word.range(30, 27) = seed;
  return word;
}

ap_uint<48> packSignedSums(int first, int second, int third) {
  ap_uint<48> word = 0;
  word.range(11, 0) = static_cast<ap_uint<12>>(first);
  word.range(23, 12) = static_cast<ap_uint<12>>(second);
  word.range(35, 24) = static_cast<ap_uint<12>>(third);
  return word;
}

ap_uint<48> packUnsignedSum(unsigned int value) {
  ap_uint<48> word = 0;
  word.range(11, 0) = value;
  return word;
}

void putObject(std::array<ap_uint<576>, kInputLinks>& links,
               unsigned int link,
               unsigned int slot,
               ap_uint<48> object) {
  links[link].range(slot * 48 + 47, slot * 48) = object;
}
}  // namespace

class MHHTestVectorProducer : public edm::stream::EDProducer<> {
public:
  explicit MHHTestVectorProducer(const edm::ParameterSet&);
  ~MHHTestVectorProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;
  void fillPattern(std::array<ap_uint<576>, kInputLinks>&, unsigned long long) const;
};

MHHTestVectorProducer::MHHTestVectorProducer(const edm::ParameterSet&) {
  for (unsigned int link = 0; link < kInputLinks; ++link) {
    produces<std::vector<uint64_t>>(std::string("LinkIn") + std::to_string(link));
  }
}

void MHHTestVectorProducer::fillPattern(std::array<ap_uint<576>, kInputLinks>& links,
                                       unsigned long long eventNumber) const {
  for (auto& link : links) {
    link = 0;
  }

  switch ((eventNumber - 1ULL) % 8ULL) {
    case 0:
      break;
    case 1:
      putObject(links, 0, 0, packGamma(40, 70, 20));
      break;
    case 2:
      putObject(links, 8, 0, packGamma(55, 60, 100));
      break;
    case 3:
      putObject(links, 0, 0, packGamma(30, 85, 44));
      putObject(links, 8, 0, packGamma(20, 84, 44));
      break;
    case 4:
      putObject(links, 0, 0, packGamma(30, 85, 44));
      putObject(links, 8, 0, packGamma(20, 84, 49));
      break;
    case 5:
      putObject(links, 1, 0, packHadron(35, 6, 7, 3));
      putObject(links, 9, 0, packHadron(25, 5, 7, 2));
      break;
    case 6:
      putObject(links, 2, 0, packSignedSums(10, -20, 30));
      putObject(links, 2, 1, packSignedSums(-5, 15, -25));
      putObject(links, 2, 2, packUnsignedSum(100));
      putObject(links, 2, 3, packUnsignedSum(200));
      putObject(links, 2, 4, packUnsignedSum(3));
      putObject(links, 10, 0, packSignedSums(7, 8, 9));
      putObject(links, 10, 1, packSignedSums(4, 5, 6));
      putObject(links, 10, 2, packUnsignedSum(50));
      putObject(links, 10, 3, packUnsignedSum(75));
      putObject(links, 10, 4, packUnsignedSum(4));
      break;
    case 7:
      for (unsigned int link : {3U, 4U, 5U, 6U, 7U, 11U, 12U, 13U, 14U, 15U}) {
        links[link].range(63, 0) = 0xfeed000000000000ULL + link;
      }
      break;
  }
}

void MHHTestVectorProducer::produce(edm::Event& event, const edm::EventSetup&) {
  std::array<ap_uint<576>, kInputLinks> links{};
  fillPattern(links, event.id().event());

  for (unsigned int link = 0; link < kInputLinks; ++link) {
    auto words = std::make_unique<std::vector<uint64_t>>();
    words->reserve(9);
    for (unsigned int word = 0; word < 9; ++word) {
      words->push_back(links[link].range(word * 64 + 63, word * 64).to_uint64());
    }
    event.put(std::move(words), std::string("LinkIn") + std::to_string(link));
  }
}

void MHHTestVectorProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  descriptions.add("mhhTestVectorProducer", description);
}

DEFINE_FWK_MODULE(MHHTestVectorProducer);
