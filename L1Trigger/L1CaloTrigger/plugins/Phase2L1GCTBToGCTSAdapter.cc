// Private proof-of-concept bridge from the current GCTB IP2 frames to the
// GCTS source ABI v2. The derived sums and empty isolated-EG slots are not
// claimed to be GCTB firmware outputs.
#include <ap_int.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/L1TCalorimeterPhase2/interface/GCT_output.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace {
  constexpr unsigned int kWords = 9;
  constexpr unsigned int kObjects = 6;
  constexpr unsigned int kSTEta = 6;
  constexpr unsigned int kSTPhiPerCard = 8;
  constexpr double kPi = 3.14159265358979323846;
  using Frame = ap_uint<576>;
  using Collection = l1tp2::gctOutputLinkCollection;

  unsigned int clampUnsigned(long value) { return static_cast<unsigned int>(std::clamp(value, 0L, 4095L)); }
  unsigned int clampSigned(long value) {
    return static_cast<unsigned int>(std::clamp(value, -2048L, 2047L)) & 0xFFFU;
  }
  double phiAt(unsigned int card, unsigned int localPhi, bool positive, double phase) {
    const unsigned int physicalLocalPhi = positive ? localPhi : kSTPhiPerCard - 1 - localPhi;
    return phase + 2.0 * kPi * (card * kSTPhiPerCard + physicalLocalPhi + 0.5) / 24.0;
  }
  std::vector<uint64_t> toWords(const Frame& frame) {
    std::vector<uint64_t> words;
    words.reserve(kWords);
    for (unsigned int i = 0; i < kWords; ++i)
      words.push_back(frame.range(64 * i + 63, 64 * i).to_uint64());
    return words;
  }
  Frame cleanEG(const Frame& source, bool positive) {
    Frame result = 0;
    // IP2's literal selected indices are 0,1,8,9,15,22. They came from
    // four RCT-region links, but the source retains only 0..29 local phi.
    // RCT_IO.h reverses local phi in negative-eta cards.
    constexpr std::array<unsigned int, 6> positivePair{{0, 0, 1, 1, 1, 2}};
    constexpr std::array<unsigned int, 6> negativePair{{3, 3, 2, 2, 2, 1}};
    for (unsigned int i = 0; i < kObjects; ++i) {
      const ap_uint<48> object = source.range(48 * i + 47, 48 * i);
      const unsigned int pt = object.range(11, 0).to_uint();
      if (pt == 0)
        continue;
      const unsigned int eta = object.range(18, 12).to_uint();
      const unsigned int rctPhi = object.range(25, 19).to_uint();
      if (eta >= 85 || rctPhi >= 30)
        throw cms::Exception("GCTBToGCTS") << "GCTB EG coordinate outside the RCT card range";
      const unsigned int pair = positive ? positivePair[i] : negativePair[i];
      const unsigned int phi = 30 * pair + (positive ? rctPhi : 29 - rctPhi);
      ap_uint<48> packed = pt;
      packed |= ap_uint<48>(eta) << 12;
      packed |= ap_uint<48>(phi) << 19;
      result.range(48 * i + 47, 48 * i) = packed;
    }
    // The upper six records are deliberately empty: current IP2 has no iso-EG.
    return result;
  }
  Frame cleanJetTau(const Frame& source, unsigned int card, bool positive, double phase,
                    double& htx, double& hty, unsigned int& ht) {
    Frame result = 0;
    for (unsigned int i = 0; i < 2 * kObjects; ++i) {
      const ap_uint<48> object = source.range(48 * i + 47, 48 * i);
      const unsigned int pt = object.range(11, 0).to_uint();
      if (pt == 0)
        continue;
      const unsigned int paddedEta = object.range(17, 12).to_uint();
      const unsigned int paddedPhi = object.range(24, 18).to_uint();
      if (paddedEta < 2 || paddedEta >= 2 + kSTEta || paddedPhi < 2 ||
          paddedPhi >= 2 + kSTPhiPerCard)
        throw cms::Exception("GCTBToGCTS") << "Nonzero GCTB jet/tau has invalid padded coordinates";
      const unsigned int eta = paddedEta - 2;
      const unsigned int localPhi = paddedPhi - 2;
      const unsigned int phi = positive ? localPhi : kSTPhiPerCard - 1 - localPhi;
      const unsigned int ratio = object.range(30, 27).to_uint();
      ap_uint<48> packed = pt;
      packed |= ap_uint<48>(eta) << 12;
      packed |= ap_uint<48>(phi) << 18;
      packed |= ap_uint<48>(ratio) << 27;
      result.range(48 * i + 47, 48 * i) = packed;
      if (i < kObjects) {
        const double angle = phiAt(card, localPhi, positive, phase);
        ht += pt;
        htx += pt * std::cos(angle);
        hty += pt * std::sin(angle);
      }
    }
    return result;
  }
  Frame derivedSums(const Collection& pre, unsigned int card, bool positive, double phase,
                    double htx, double hty, unsigned int ht) {
    // IP2's six ST inputs per side include two overlap columns on each edge.
    // Only input positions 1..4 are unique to this card's 8-column core.
    constexpr std::array<unsigned int, 4> positiveST{{9, 10, 17, 18}};
    constexpr std::array<unsigned int, 4> negativeST{{21, 22, 13, 14}};
    const auto& positions = positive ? positiveST : negativeST;
    double ex = 0.0, ey = 0.0;
    unsigned int et = 0, nObj = 0;
    for (unsigned int link = 0; link < positions.size(); ++link) {
      const Frame frame = pre[positions[link]].data();
      for (unsigned int column = 0; column < 2; ++column) {
        const unsigned int localPhi = 2 * link + column;
        const double angle = phiAt(card, localPhi, positive, phase);
        for (unsigned int eta = 0; eta < kSTEta; ++eta) {
          const unsigned int slot = column * kSTEta + eta;
          const unsigned int pt = frame.range(48 * slot + 11, 48 * slot).to_uint();
          et += pt;
          nObj += pt != 0;
          ex += pt * std::cos(angle);
          ey += pt * std::sin(angle);
        }
      }
    }
    Frame result = 0;
    result.range(11, 0) = clampSigned(std::lround(ex));
    result.range(23, 12) = clampSigned(std::lround(ey));
    result.range(35, 24) = clampSigned(std::lround(htx));
    result.range(47, 36) = clampSigned(std::lround(hty));
    result.range(75, 64) = clampUnsigned(et);
    result.range(87, 76) = clampUnsigned(ht);
    result.range(99, 88) = clampUnsigned(nObj);
    return result;
  }
}  // namespace

class Phase2L1GCTBToGCTSAdapter : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1GCTBToGCTSAdapter(const edm::ParameterSet& config)
      : pre_(consumes<Collection>(config.getParameter<edm::InputTag>("preIP2"))),
        post_(consumes<Collection>(config.getParameter<edm::InputTag>("postIP2"))),
        card_(config.getParameter<unsigned int>("cardIndex")),
        phiPhase_(config.getParameter<double>("phiPhaseRadians")) {
    if (card_ >= 3)
      throw cms::Exception("GCTBToGCTS") << "cardIndex must be 0, 1, or 2";
    if (!std::isfinite(phiPhase_))
      throw cms::Exception("GCTBToGCTS") << "phiPhaseRadians must be finite";
    for (const char* side : {"Positive", "Negative"})
      for (unsigned int i = 0; i < 3; ++i)
        produces<std::vector<uint64_t>>(std::string("LinkOut") + side + std::to_string(i));
    produces<std::vector<uint64_t>>("AssumptionFlags");
  }
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription description;
    description.add<edm::InputTag>("preIP2", edm::InputTag());
    description.add<edm::InputTag>("postIP2", edm::InputTag());
    description.add<unsigned int>("cardIndex", 0);
    description.add<double>("phiPhaseRadians", -kPi);
    descriptions.add("l1tPhase2L1GCTBToGCTSAdapter", description);
  }
private:
  void produce(edm::Event& event, const edm::EventSetup&) override {
    const auto pre = event.getHandle(pre_);
    const auto post = event.getHandle(post_);
    if (!pre.isValid() || pre->size() != 24 || !post.isValid() || post->size() != 6)
      throw cms::Exception("GCTBToGCTS") << "Expected exactly 24 PreIP2 and 6 PostIP2 frames";
    for (unsigned int side = 0; side < 2; ++side) {
      const bool positive = side == 0;
      const std::string label = positive ? "Positive" : "Negative";
      double htx = 0.0, hty = 0.0;
      unsigned int ht = 0;
      const Frame eg = cleanEG((*post)[positive ? 0 : 3].data(), positive);
      const Frame jets = cleanJetTau((*post)[positive ? 1 : 4].data(), card_, positive, phiPhase_, htx, hty, ht);
      const Frame sums = derivedSums(*pre, card_, positive, phiPhase_, htx, hty, ht);
      event.put(std::make_unique<std::vector<uint64_t>>(toWords(eg)), "LinkOut" + label + "0");
      event.put(std::make_unique<std::vector<uint64_t>>(toWords(jets)), "LinkOut" + label + "1");
      event.put(std::make_unique<std::vector<uint64_t>>(toWords(sums)), "LinkOut" + label + "2");
    }
    // bit 0: iso-EG absent; bit 1: sums derived from ST/jets;
    // bit 2: padded jet/tau indices normalized; bit 3: phi phase assumed.
    event.put(std::make_unique<std::vector<uint64_t>>(std::vector<uint64_t>{15}), "AssumptionFlags");
  }
  edm::EDGetTokenT<Collection> pre_, post_;
  unsigned int card_;
  double phiPhase_;
};

DEFINE_FWK_MODULE(Phase2L1GCTBToGCTSAdapter);
