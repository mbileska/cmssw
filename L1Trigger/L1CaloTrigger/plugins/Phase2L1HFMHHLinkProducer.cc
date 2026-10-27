/*
 * Additive HF source-link producer for the MHH emulator.
 *
 * The existing HF emulator is intentionally not modified.  Its published
 * mixer jet/tau records contain region-local phi after candidates from three
 * regions have been sorted, so global phi cannot be recovered downstream.
 * This producer instead consumes the same HF trigger primitives, builds a
 * global integration-level MHH jet/tau source, and emits source sums ABI v2.
 * It does not replace or claim bit-exact equivalence to the HF IP1 firmware;
 * that requires an upstream HF interface retaining global phi.
 *
 * Link 0 (EG/isolated EG) is zero because the current HF algorithm defines no
 * such collections.  Link 1 contains six 48-bit jets followed by six 48-bit
 * taus.  Link 2 contains Ex/Ey/HTx/HTy, total ET/HT, and NObj.  Components are
 * global energy-flow sums and are not negated missing-energy vectors.
 */

#include <ap_int.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "CalibFormats/CaloTPG/interface/CaloTPGRecord.h"
#include "CalibFormats/CaloTPG/interface/CaloTPGTranscoder.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1CaloJetEmulator.h"
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1MHHLinkPacker.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

namespace {

  constexpr unsigned int kWordsPerLink = 9;
  constexpr unsigned int kOutputLinks = 3;
  constexpr unsigned int kObjectsPerCategory = 6;
  constexpr unsigned int kHFHadEtaOffset = 6;
  constexpr unsigned int kHadPhiBins = 24;
  constexpr unsigned int kFineEtaBins = nHfEta;
  constexpr unsigned int kFinePhiBins = nHfPhi;
  constexpr double kPi = 3.14159265358979323846;

  struct Candidate {
    unsigned int et = 0;
    unsigned int eta = 0;
    unsigned int phi = 0;
    unsigned int seedRatio = 0;
  };

  unsigned int quantizeUnsigned(double value, double lsb, unsigned int maximum) {
    if (value <= 0.) {
      return 0;
    }
    const long quantized = std::lround(value / lsb);
    return static_cast<unsigned int>(std::clamp<long>(quantized, 0, maximum));
  }

  int quantizeSigned(double value, double lsb) {
    const long quantized = std::lround(value / lsb);
    return static_cast<int>(std::clamp<long>(quantized, -2048, 2047));
  }

  unsigned int seedRatio(double totalEt, double seedEt) {
    if (totalEt <= 0. || seedEt <= 0.) {
      return 0;
    }
    return static_cast<unsigned int>(std::clamp<long>(std::lround(totalEt / seedEt), 0, 0xF));
  }

  void sortAndTrim(std::vector<Candidate>& candidates) {
    std::stable_sort(candidates.begin(), candidates.end(), [](const Candidate& lhs, const Candidate& rhs) {
      return lhs.et > rhs.et;
    });
    if (candidates.size() > kObjectsPerCategory) {
      candidates.resize(kObjectsPerCategory);
    }
  }

  void putCandidates(ap_uint<576>& link, unsigned int slotOffset, const std::vector<Candidate>& candidates) {
    for (unsigned int i = 0; i < candidates.size(); ++i) {
      const unsigned int slot = slotOffset + i;
      link.range(slot * 48 + 47, slot * 48) = p2mhhlink::packHadron(
          candidates[i].et, candidates[i].eta, candidates[i].phi, candidates[i].seedRatio);
    }
  }

  std::vector<uint64_t> toWords(const ap_uint<576>& link) {
    std::vector<uint64_t> words;
    words.reserve(kWordsPerLink);
    for (unsigned int i = 0; i < kWordsPerLink; ++i) {
      words.push_back(link.range(i * 64 + 63, i * 64).to_uint64());
    }
    return words;
  }

}  // namespace

class Phase2L1HFMHHLinkProducer : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1HFMHHLinkProducer(const edm::ParameterSet&);
  ~Phase2L1HFMHHLinkProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<HcalTrigPrimDigiCollection> hfToken_;
  edm::ESGetToken<CaloTPGTranscoder, CaloTPGRecord> decoderToken_;
  int zSide_;
  double etLsb_;
  double minTowerEt_;
  double jetSeedThreshold_;
  bool debug_;
};

Phase2L1HFMHHLinkProducer::Phase2L1HFMHHLinkProducer(const edm::ParameterSet& config)
    : hfToken_(consumes<HcalTrigPrimDigiCollection>(config.getParameter<edm::InputTag>("hcalDigis"))),
      decoderToken_(esConsumes<CaloTPGTranscoder, CaloTPGRecord>()),
      zSide_(config.getParameter<int>("zSide")),
      etLsb_(config.getParameter<double>("etLsb")),
      minTowerEt_(config.getParameter<double>("minTowerEt")),
      jetSeedThreshold_(config.getParameter<double>("jetSeedThreshold")),
      debug_(config.getParameter<bool>("debug")) {
  if (zSide_ != -1 && zSide_ != 1) {
    throw cms::Exception("Phase2L1HFMHHLinkProducer") << "zSide must be either -1 or +1";
  }
  if (etLsb_ <= 0.) {
    throw cms::Exception("Phase2L1HFMHHLinkProducer") << "etLsb must be positive";
  }
  if (minTowerEt_ < 0. || jetSeedThreshold_ < 0.) {
    throw cms::Exception("Phase2L1HFMHHLinkProducer") << "energy thresholds must be non-negative";
  }

  for (unsigned int link = 0; link < kOutputLinks; ++link) {
    produces<std::vector<uint64_t>>(std::string("LinkIn") + std::to_string(link));
  }
}

void Phase2L1HFMHHLinkProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  const auto hfHandle = event.getHandle(hfToken_);
  if (!hfHandle.isValid()) {
    throw cms::Exception("Phase2L1HFMHHLinkProducer") << "Configured HF trigger-primitive collection is missing";
  }
  const auto& decoder = setup.getData(decoderToken_);

  std::array<ap_uint<576>, kOutputLinks> links{};
  float towerGrid[kFineEtaBins][kFinePhiBins] = {};
  double sumEt = 0.;
  double ex = 0.;
  double ey = 0.;
  unsigned int nObjects = 0;

  for (const auto& hit : *hfHandle) {
    const int iEta = hit.id().ieta();
    const int absEta = std::abs(iEta);
    if ((iEta * zSide_) <= 0 || absEta < 30 || absEta > 41 || hit.id().version() != 1) {
      continue;
    }

    const double et = std::max(0., static_cast<double>(decoder.hcaletValue(hit.id(), hit.t0())));
    if (et <= minTowerEt_) {
      continue;
    }

    const int localEta = absEta - 30;
    const int iPhi = hit.id().iphi() - 1;
    if (localEta < 0 || localEta >= 12 || iPhi < 0 || iPhi >= static_cast<int>(kFinePhiBins)) {
      continue;
    }

    // Match the existing Phase-2 HF tower splitting: each physical eta bin is
    // split in two; 10-degree and 20-degree towers are also split into the
    // corresponding 5-degree phi bins.  The split conserves physical ET.
    towerGrid[2 * localEta][iPhi] = static_cast<float>(et / 2.);
    towerGrid[2 * localEta + 1][iPhi] = static_cast<float>(et / 2.);
    if ((localEta < 2 || localEta >= 10) && (iPhi % 4 == 2)) {
      for (int phiStep = 0; phiStep < 4; ++phiStep) {
        const int splitPhi = (iPhi + phiStep) % static_cast<int>(kFinePhiBins);
        towerGrid[2 * localEta][splitPhi] = static_cast<float>(et / 8.);
        towerGrid[2 * localEta + 1][splitPhi] = static_cast<float>(et / 8.);
      }
    } else if (localEta >= 2 && localEta < 10 && (iPhi % 2 == 0)) {
      for (int phiStep = 0; phiStep < 2; ++phiStep) {
        const int splitPhi = (iPhi + phiStep) % static_cast<int>(kFinePhiBins);
        towerGrid[2 * localEta][splitPhi] = static_cast<float>(et / 4.);
        towerGrid[2 * localEta + 1][splitPhi] = static_cast<float>(et / 4.);
      }
    }

    const double phi = l1t::CaloTools::towerPhi(iEta, hit.id().iphi());
    sumEt += et;
    ex += et * std::cos(phi);
    ey += et * std::sin(phi);
    ++nObjects;
  }

  gctobj::GCTsupertower_t superTowers[nSTEta][nSTPhi];
  gctobj::makeST_hf(towerGrid, superTowers);

  std::vector<Candidate> jetCandidates;
  std::vector<Candidate> tauCandidates;
  double ht = 0.;
  double htx = 0.;
  double hty = 0.;
  for (int i = 0; i < nJets; ++i) {
    const gctobj::jetInfo jet = gctobj::getRegion(superTowers, static_cast<float>(jetSeedThreshold_));
    if (jet.energy <= 0.) {
      continue;
    }

    Candidate jetCandidate;
    jetCandidate.et = quantizeUnsigned(jet.energy, etLsb_, 0xFFF);
    jetCandidate.eta = kHFHadEtaOffset + std::min(static_cast<unsigned int>(jet.etaCenter / 3), 7U);
    jetCandidate.phi = static_cast<unsigned int>(jet.phiCenter / 3) % kHadPhiBins;
    jetCandidate.seedRatio = seedRatio(jet.energy, jet.energyMax);
    if (jetCandidate.et != 0) {
      jetCandidates.push_back(jetCandidate);
      ht += jet.energy;
      const double jetPhi =
          (static_cast<double>(jet.phiCenter) + 0.5) * 2. * kPi / static_cast<double>(kFinePhiBins);
      htx += jet.energy * std::cos(jetPhi);
      hty += jet.energy * std::sin(jetPhi);
    }

    Candidate tauCandidate = jetCandidate;
    tauCandidate.et = quantizeUnsigned(jet.tauEt, etLsb_, 0xFFF);
    tauCandidate.seedRatio = seedRatio(jet.tauEt, jet.energyMax);
    if (tauCandidate.et != 0) {
      tauCandidates.push_back(tauCandidate);
    }
  }
  sortAndTrim(jetCandidates);
  sortAndTrim(tauCandidates);

  // Link 0 stays zero: the current HF algorithm has no EG or isolated-EG
  // source collection for the MHH object contract.
  putCandidates(links[1], 0, jetCandidates);
  putCandidates(links[1], kObjectsPerCategory, tauCandidates);
  links[2] = p2mhhlink::packSumsV2(quantizeSigned(ex, etLsb_),
                                    quantizeSigned(ey, etLsb_),
                                    quantizeSigned(htx, etLsb_),
                                    quantizeSigned(hty, etLsb_),
                                    quantizeUnsigned(sumEt, etLsb_, 0xFFF),
                                    quantizeUnsigned(ht, etLsb_, 0xFFF),
                                    std::min(nObjects, 0xFFFU));

  if (debug_) {
    edm::LogVerbatim("Phase2L1HFMHHLinkProducer")
        << "event=" << event.id().event() << " zSide=" << zSide_ << " jets=" << jetCandidates.size()
        << " taus=" << tauCandidates.size() << " towers=" << nObjects << " Ex=" << ex << " Ey=" << ey
        << " HTx=" << htx << " HTy=" << hty << " ET=" << sumEt << " HT=" << ht;
  }

  for (unsigned int link = 0; link < kOutputLinks; ++link) {
    event.put(std::make_unique<std::vector<uint64_t>>(toWords(links[link])),
              std::string("LinkIn") + std::to_string(link));
  }
}

void Phase2L1HFMHHLinkProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("hcalDigis", edm::InputTag("simHcalTriggerPrimitiveDigis"));
  description.add<int>("zSide", 1);
  description.add<double>("etLsb", 0.5);
  description.add<double>("minTowerEt", 0.0);
  description.add<double>("jetSeedThreshold", 3.0);
  description.add<bool>("debug", false);
  descriptions.add("l1tPhase2L1HFMHHLinkProducer", description);
}

DEFINE_FWK_MODULE(Phase2L1HFMHHLinkProducer);
