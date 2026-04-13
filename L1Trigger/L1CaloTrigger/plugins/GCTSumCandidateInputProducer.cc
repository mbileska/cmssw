/*
 * Description: MC-backed input producer for the Phase 2 GCT SumCard emulator.
 *
 * This producer only packs the parts of the GCT Sum input contract that are
 * exposed by CMSSW collections in this repo:
 * - barrel EG / isoEG candidates from GCT EGammas
 * - barrel and HGCal hadronic candidates from GCT jets / taus
 *
 * The remaining upstream contracts are intentionally left explicit:
 * - endcap EM links are left empty
 * - partial sum words are zero-filled
 * - HF hadronic objects are not injected because this SumCard interface only
 *   documents one non-barrel source per eta side in the available code
 */

#include <ap_int.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/L1TCalorimeterPhase2/interface/Phase2L1CaloJet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

namespace {

constexpr unsigned int kSides = 2;
constexpr unsigned int kSourcesPerSide = 4;
constexpr unsigned int kInputLinks = 24;
constexpr unsigned int kLinksPerSource = 3;
constexpr unsigned int kWordsPerLink = 9;
constexpr unsigned int kObjectsPerCategory = 6;

constexpr int kPositiveGammaEtaOffset = 85;
constexpr int kNegativeGammaEtaMax = 84;
constexpr int kGammaPhiRotation = 100;

constexpr int kBarrelTowerPhiRotation = 20;
constexpr int kPositiveBarrelTowerEtaOffset = 17;
constexpr int kNegativeBarrelTowerEtaMax = 16;

constexpr int kPositiveHgcalTowerEtaOffset = 18;
constexpr int kNegativeHgcalTowerEtaMax = 17;

constexpr float kBarrelEtaMax = 1.5f;
constexpr float kHgcalEtaMin = 1.479f;
constexpr float kHgcalEtaMax = 3.0f;

struct PackedWord {
  ap_uint<48> word = 0;
  unsigned int et = 0;
};

using SourceWordLists = std::array<std::vector<PackedWord>, kSides * kSourcesPerSide>;

int wrapModulo(int value, int modulo) {
  int wrapped = value % modulo;
  if (wrapped < 0) {
    wrapped += modulo;
  }
  return wrapped;
}

unsigned int quantizeEt(float et) {
  if (et <= 0.f) {
    return 0;
  }
  const long rounded = std::lround(et);
  if (rounded <= 0L) {
    return 0;
  }
  if (rounded > 0xFFF) {
    return 0xFFF;
  }
  return static_cast<unsigned int>(rounded);
}

unsigned int ratioFromSeed(float totalEt, float seedEt) {
  if (totalEt <= 0.f || seedEt <= 0.f) {
    return 0;
  }
  const long rounded = std::lround(totalEt / seedEt);
  if (rounded <= 0L) {
    return 0;
  }
  if (rounded > 0xF) {
    return 0xF;
  }
  return static_cast<unsigned int>(rounded);
}

ap_uint<48> packGammaWord(unsigned int et, unsigned int eta, unsigned int phi, bool isBarrel) {
  ap_uint<48> word = 0;
  word.range(11, 0) = et;
  word.range(18, 12) = eta;
  word.range(27, 19) = phi;
  word.range(47, 47) = isBarrel ? 1 : 0;
  return word;
}

ap_uint<48> packHadWord(unsigned int et, unsigned int eta, unsigned int phi, unsigned int ratio, bool isBarrel) {
  ap_uint<48> word = 0;
  word.range(11, 0) = et;
  word.range(17, 12) = eta;
  word.range(26, 18) = phi;
  word.range(30, 27) = ratio;
  word.range(47, 47) = isBarrel ? 1 : 0;
  return word;
}

void sortAndTrim(std::vector<PackedWord>& words) {
  std::stable_sort(words.begin(), words.end(), [](const PackedWord& lhs, const PackedWord& rhs) {
    return lhs.et > rhs.et;
  });
  if (words.size() > kObjectsPerCategory) {
    words.resize(kObjectsPerCategory);
  }
}

unsigned int sideSourceIndex(bool isPositiveEta, unsigned int source) {
  return (isPositiveEta ? 0U : kSourcesPerSide) + source;
}

void writeWords(ap_uint<576>& link, unsigned int slotOffset, const std::vector<PackedWord>& words) {
  for (unsigned int i = 0; i < words.size() && i < kObjectsPerCategory; ++i) {
    const unsigned int slot = slotOffset + i;
    link.range((slot * 48) + 47, slot * 48) = words[i].word;
  }
}

std::vector<uint64_t> unpackToWords(ap_uint<576> link) {
  std::vector<uint64_t> words;
  words.reserve(kWordsPerLink);
  for (unsigned int i = 0; i < kWordsPerLink; ++i) {
    words.push_back(link.range((i * 64) + 63, i * 64).to_uint64());
  }
  return words;
}

}  // namespace

class GCTSumCandidateInputProducer : public edm::stream::EDProducer<> {
public:
  explicit GCTSumCandidateInputProducer(const edm::ParameterSet&);
  ~GCTSumCandidateInputProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<l1t::EGammaBxCollection> gctEGammasToken_;
  edm::EDGetTokenT<l1tp2::Phase2L1CaloJetCollection> gctJetsToken_;
  bool debug_;
};

GCTSumCandidateInputProducer::GCTSumCandidateInputProducer(const edm::ParameterSet& iConfig)
    : gctEGammasToken_(consumes<l1t::EGammaBxCollection>(iConfig.getParameter<edm::InputTag>("gctEGammas"))),
      gctJetsToken_(consumes<l1tp2::Phase2L1CaloJetCollection>(iConfig.getParameter<edm::InputTag>("gctJets"))),
      debug_(iConfig.getParameter<bool>("debug")) {
  for (unsigned int i = 0; i < kInputLinks; ++i) {
    produces<std::vector<uint64_t> >(std::string("LinkIn") + std::to_string(i));
  }
}

void GCTSumCandidateInputProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  (void)iSetup;

  SourceWordLists egWords;
  SourceWordLists egiWords;
  SourceWordLists jetWords;
  SourceWordLists tauWords;

  edm::Handle<l1t::EGammaBxCollection> gctEGammas;
  iEvent.getByToken(gctEGammasToken_, gctEGammas);
  if (gctEGammas.isValid()) {
    for (auto it = gctEGammas->begin(0); it != gctEGammas->end(0); ++it) {
      const l1t::EGamma& eg = *it;
      const unsigned int hwPt = std::min<unsigned int>(eg.hwPt(), 0xFFF);
      if (hwPt == 0U) {
        continue;
      }

      const int hwEta = eg.hwEta();
      const int hwPhi = eg.hwPhi();
      const bool isPositiveEta = (hwEta >= kPositiveGammaEtaOffset);
      const int localEta = isPositiveEta ? (hwEta - kPositiveGammaEtaOffset) : (kNegativeGammaEtaMax - hwEta);
      if (localEta < 0 || localEta > 84) {
        continue;
      }

      const int sumGlobalPhi = wrapModulo(hwPhi - kGammaPhiRotation, 360);
      const unsigned int source = 1U + static_cast<unsigned int>(sumGlobalPhi / 120);
      if (source >= kSourcesPerSide) {
        continue;
      }
      const unsigned int localPhi = static_cast<unsigned int>(sumGlobalPhi % 120);

      const PackedWord word{packGammaWord(hwPt, static_cast<unsigned int>(localEta), localPhi, true), hwPt};
      const unsigned int sourceIndex = sideSourceIndex(isPositiveEta, source);
      egWords[sourceIndex].push_back(word);
      if ((eg.hwQual() & 0x1) != 0) {
        egiWords[sourceIndex].push_back(word);
      }
    }
  }

  edm::Handle<l1tp2::Phase2L1CaloJetCollection> gctJets;
  iEvent.getByToken(gctJetsToken_, gctJets);
  if (gctJets.isValid()) {
    for (const auto& jet : *gctJets) {
      const float absEta = std::abs(jet.jetEta());
      const bool isBarrel = (absEta <= kBarrelEtaMax);
      const bool isHgcal = (absEta > kHgcalEtaMin && absEta <= kHgcalEtaMax);
      if (!isBarrel && !isHgcal) {
        continue;
      }

      const int rawEta = jet.jetIEta();
      const int rawPhi = jet.jetIPhi();
      const unsigned int sumGlobalPhi = static_cast<unsigned int>(wrapModulo(rawPhi - kBarrelTowerPhiRotation, 72) / 3);

      const bool isPositiveEta = isBarrel ? (rawEta >= kPositiveBarrelTowerEtaOffset) : (rawEta >= kPositiveHgcalTowerEtaOffset);
      int localEta = 0;
      unsigned int source = 0;
      unsigned int localPhi = 0;

      if (isBarrel) {
        source = 1U + static_cast<unsigned int>(sumGlobalPhi / 8);
        if (source >= kSourcesPerSide) {
          continue;
        }
        localPhi = sumGlobalPhi % 8;
        localEta = isPositiveEta ? ((rawEta - kPositiveBarrelTowerEtaOffset) / 3) : ((kNegativeBarrelTowerEtaMax - rawEta) / 3);
      } else {
        source = 0U;
        localPhi = sumGlobalPhi;
        localEta = isPositiveEta ? ((rawEta - kPositiveHgcalTowerEtaOffset) / 3) : ((kNegativeHgcalTowerEtaMax - rawEta) / 3);
      }

      if (localEta < 0 || localEta > 5) {
        continue;
      }

      const unsigned int jetEt = quantizeEt(jet.jetEt());
      const unsigned int tauEt = quantizeEt(jet.tauEt());
      const unsigned int jetRatio = ratioFromSeed(jet.jetEt(), jet.towerEt());
      const unsigned int tauRatio = ratioFromSeed(jet.tauEt(), jet.towerEt());
      const unsigned int sourceIndex = sideSourceIndex(isPositiveEta, source);

      if (jetEt > 0U) {
        jetWords[sourceIndex].push_back(
            PackedWord{packHadWord(jetEt, static_cast<unsigned int>(localEta), localPhi, jetRatio, isBarrel), jetEt});
      }
      if (tauEt > 0U) {
        tauWords[sourceIndex].push_back(
            PackedWord{packHadWord(tauEt, static_cast<unsigned int>(localEta), localPhi, tauRatio, isBarrel), tauEt});
      }
    }
  }

  for (auto& words : egWords)
    sortAndTrim(words);
  for (auto& words : egiWords)
    sortAndTrim(words);
  for (auto& words : jetWords)
    sortAndTrim(words);
  for (auto& words : tauWords)
    sortAndTrim(words);

  std::array<ap_uint<576>, kInputLinks> links;
  links.fill(0);

  for (unsigned int side = 0; side < kSides; ++side) {
    const bool isPositiveEta = (side == 0U);
    const unsigned int sideLinkOffset = isPositiveEta ? 0U : 12U;
    for (unsigned int source = 0; source < kSourcesPerSide; ++source) {
      const unsigned int sourceIndex = sideSourceIndex(isPositiveEta, source);
      const unsigned int linkBase = sideLinkOffset + (source * kLinksPerSource);

      writeWords(links[linkBase + 0], 0, egWords[sourceIndex]);
      writeWords(links[linkBase + 0], 6, egiWords[sourceIndex]);
      writeWords(links[linkBase + 1], 0, jetWords[sourceIndex]);
      writeWords(links[linkBase + 1], 6, tauWords[sourceIndex]);
    }
  }

  if (debug_) {
    unsigned int nEg = 0;
    unsigned int nEgi = 0;
    unsigned int nJet = 0;
    unsigned int nTau = 0;
    for (const auto& words : egWords)
      nEg += words.size();
    for (const auto& words : egiWords)
      nEgi += words.size();
    for (const auto& words : jetWords)
      nJet += words.size();
    for (const auto& words : tauWords)
      nTau += words.size();

    edm::LogVerbatim("GCTSumCandidateInputProducer") << "event=" << iEvent.id().event()
                                                     << " packedEG=" << nEg
                                                     << " packedEGI=" << nEgi
                                                     << " packedJet=" << nJet
                                                     << " packedTau=" << nTau
                                                     << " sums=0";
  }

  for (unsigned int link = 0; link < kInputLinks; ++link) {
    auto outWords = std::make_unique<std::vector<uint64_t> >(unpackToWords(links[link]));
    iEvent.put(std::move(outWords), std::string("LinkIn") + std::to_string(link));
  }
}

void GCTSumCandidateInputProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("gctEGammas", edm::InputTag("l1tPhase2L1CaloEGammaEmulator", "GCTEGammas"));
  desc.add<edm::InputTag>("gctJets", edm::InputTag("l1tPhase2CaloJetEmulator", "GCTJet"));
  desc.add<bool>("debug", false);
  descriptions.add("gctSumCandidateInputProducer", desc);
}

DEFINE_FWK_MODULE(GCTSumCandidateInputProducer);
