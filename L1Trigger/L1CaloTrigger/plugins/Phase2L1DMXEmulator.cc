/*
 * Provisional Phase-2 HGCAL-to-MHH link producer.
 *
 * This module provides the CMSSW plumbing needed for an alpha integration of
 * the DMX/MHH chain.  It converts BX=0 HGCAL trigger towers and multiclusters
 * into the 16 x 576-bit input-link contract consumed by Phase2L1MHHEmulator.
 * Only HGCAL links 8, 9, and 10 are populated; all other links are zero-filled.
 *
 * The link serialization matches MHH sums ABI v2.  Ex/Ey and total ET are
 * formed from the accepted HGCAL towers; NObj counts those towers.  HTx/HTy
 * and total HT are formed from the accepted jets.  Cartesian components are
 * global energy-flow sums, not negated missing-energy vectors.  The
 * object-building and isolation choices remain provisional until an
 * authoritative DMX specification exists.
 */

#include <ap_int.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/L1THGCal/interface/HGCalTower.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "L1Trigger/L1CaloTrigger/interface/MHH_h.h"
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1CaloJetEmulator.h"
#include "L1Trigger/L1CaloTrigger/interface/Phase2L1MHHLinkPacker.h"

namespace {

constexpr unsigned int kWordsPerLink = 9;
constexpr unsigned int kObjectsPerCategory = 6;
constexpr unsigned int kHgcalLinkBase = 8;
constexpr unsigned int kGammaEtaBins = 85;
constexpr unsigned int kGammaPhiBins = 360;
constexpr unsigned int kHadEtaBins = 6;
constexpr unsigned int kHadPhiBins = 24;
constexpr unsigned int kTowerEtaBins = nHgcalEta / 2;
constexpr unsigned int kTowerPhiBins = nHgcalPhi;
constexpr double kPi = 3.14159265358979323846;

struct Candidate {
  unsigned int et = 0;
  unsigned int eta = 0;
  unsigned int phi = 0;
  unsigned int seedRatio = 0;
};

double normalizedPhi(double phi, double offset) {
  double result = std::fmod(phi + offset, 2. * kPi);
  if (result < 0.) {
    result += 2. * kPi;
  }
  return result;
}

unsigned int phiBin(double phi, unsigned int bins, double offset) {
  const double scaled = normalizedPhi(phi, offset) * static_cast<double>(bins) / (2. * kPi);
  return std::min(static_cast<unsigned int>(scaled), bins - 1);
}

unsigned int etaBin(double absEta, double etaMin, double etaMax, unsigned int bins) {
  const double fraction = (absEta - etaMin) / (etaMax - etaMin);
  const int bin = static_cast<int>(std::floor(fraction * static_cast<double>(bins)));
  return static_cast<unsigned int>(std::clamp(bin, 0, static_cast<int>(bins) - 1));
}

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

ap_uint<48> packGamma(const Candidate& candidate) {
  return p2mhhlink::packGamma(candidate.et, candidate.eta, candidate.phi);
}

ap_uint<48> packHadron(const Candidate& candidate) {
  return p2mhhlink::packHadron(candidate.et, candidate.eta, candidate.phi, candidate.seedRatio);
}

template <typename Pack>
void putCandidates(ap_uint<576>& link, unsigned int slotOffset, const std::vector<Candidate>& candidates, Pack pack) {
  for (unsigned int i = 0; i < candidates.size(); ++i) {
    const unsigned int slot = slotOffset + i;
    link.range(slot * 48 + 47, slot * 48) = pack(candidates[i]);
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

class Phase2L1DMXEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1DMXEmulator(const edm::ParameterSet&);
  ~Phase2L1DMXEmulator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<l1t::HGCalTowerBxCollection> hgcalTowersToken_;
  edm::EDGetTokenT<l1t::HGCalMulticlusterBxCollection> hgcalMulticlustersToken_;
  int zSide_;
  double etaMin_;
  double etaMax_;
  double etLsb_;
  double phiOffset_;
  double minTowerEt_;
  double minEgEt_;
  double maxEgHoverE_;
  double jetSeedThreshold_;
  bool duplicateEGAsIso_;
  bool debug_;
};

Phase2L1DMXEmulator::Phase2L1DMXEmulator(const edm::ParameterSet& config)
    : hgcalTowersToken_(
          consumes<l1t::HGCalTowerBxCollection>(config.getParameter<edm::InputTag>("hgcalTowers"))),
      hgcalMulticlustersToken_(consumes<l1t::HGCalMulticlusterBxCollection>(
          config.getParameter<edm::InputTag>("hgcalMulticlusters"))),
      zSide_(config.getParameter<int>("zSide")),
      etaMin_(config.getParameter<double>("etaMin")),
      etaMax_(config.getParameter<double>("etaMax")),
      etLsb_(config.getParameter<double>("etLsb")),
      phiOffset_(config.getParameter<double>("phiOffset")),
      minTowerEt_(config.getParameter<double>("minTowerEt")),
      minEgEt_(config.getParameter<double>("minEgEt")),
      maxEgHoverE_(config.getParameter<double>("maxEgHoverE")),
      jetSeedThreshold_(config.getParameter<double>("jetSeedThreshold")),
      duplicateEGAsIso_(config.getParameter<bool>("duplicateEGAsIso")),
      debug_(config.getParameter<bool>("debug")) {
  if (zSide_ != -1 && zSide_ != 1) {
    throw cms::Exception("Phase2L1DMXEmulator") << "zSide must be either -1 or +1";
  }
  if (etaMin_ < 0. || etaMax_ <= etaMin_) {
    throw cms::Exception("Phase2L1DMXEmulator") << "Require 0 <= etaMin < etaMax";
  }
  if (etLsb_ <= 0.) {
    throw cms::Exception("Phase2L1DMXEmulator") << "etLsb must be positive";
  }

  for (int link = 0; link < p2mhh::MHH_N_INPUT_LINKS; ++link) {
    produces<std::vector<uint64_t> >(std::string("LinkIn") + std::to_string(link));
  }
}

void Phase2L1DMXEmulator::produce(edm::Event& event, const edm::EventSetup&) {
  edm::Handle<l1t::HGCalTowerBxCollection> towers;
  event.getByToken(hgcalTowersToken_, towers);
  if (!towers.isValid()) {
    throw cms::Exception("Phase2L1DMXEmulator") << "Configured HGCAL tower collection is missing";
  }

  edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters;
  event.getByToken(hgcalMulticlustersToken_, multiclusters);
  if (!multiclusters.isValid()) {
    throw cms::Exception("Phase2L1DMXEmulator") << "Configured HGCAL multicluster collection is missing";
  }

  std::array<ap_uint<576>, p2mhh::MHH_N_INPUT_LINKS> links{};
  std::vector<Candidate> egCandidates;
  std::vector<Candidate> isoEgCandidates;
  std::vector<Candidate> jetCandidates;
  std::vector<Candidate> tauCandidates;

  for (auto cluster = multiclusters->begin(0); cluster != multiclusters->end(0); ++cluster) {
    const double absEta = std::abs(cluster->eta());
    if ((cluster->eta() * zSide_) <= 0. || absEta < etaMin_ || absEta >= etaMax_) {
      continue;
    }
    if (cluster->pt() < minEgEt_ || cluster->hOverE() > maxEgHoverE_) {
      continue;
    }

    Candidate candidate;
    candidate.et = quantizeUnsigned(cluster->pt(), etLsb_, 0xFFF);
    candidate.eta = etaBin(absEta, etaMin_, etaMax_, kGammaEtaBins);
    candidate.phi = phiBin(cluster->phi(), kGammaPhiBins, phiOffset_);
    if (candidate.et != 0) {
      egCandidates.push_back(candidate);
    }
  }
  sortAndTrim(egCandidates);
  if (duplicateEGAsIso_) {
    // Integration-only fallback: no isolation observable is supplied by the
    // multicluster input contract used here.
    isoEgCandidates = egCandidates;
  }

  float towerGrid[kTowerEtaBins][kTowerPhiBins] = {};
  double sumEt = 0.;
  double ex = 0.;
  double ey = 0.;
  unsigned int nObjects = 0;

  for (auto tower = towers->begin(0); tower != towers->end(0); ++tower) {
    const double absEta = std::abs(tower->eta());
    if (tower->zside() != zSide_ || absEta < etaMin_ || absEta >= etaMax_) {
      continue;
    }

    const double et = std::max(0., tower->etEm() + tower->etHad());
    if (et < minTowerEt_) {
      continue;
    }

    const unsigned int iEta = etaBin(absEta, etaMin_, etaMax_, kTowerEtaBins);
    const unsigned int iPhi = phiBin(tower->phi(), kTowerPhiBins, phiOffset_);
    towerGrid[iEta][iPhi] += static_cast<float>(et);
    sumEt += et;
    // The grid offset is only a binning convention.  Source-sum components
    // must remain in the detector-global Cartesian frame.
    ex += et * std::cos(tower->phi());
    ey += et * std::sin(tower->phi());
    ++nObjects;
  }

  gctobj::GCTsupertower_t superTowers[nSTEta][nSTPhi];
  gctobj::makeST_hgcal(towerGrid, superTowers);

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
    jetCandidate.eta = std::min(static_cast<unsigned int>(jet.etaCenter / 3), kHadEtaBins - 1);
    jetCandidate.phi = static_cast<unsigned int>(jet.phiCenter / 3) % kHadPhiBins;
    jetCandidate.seedRatio = seedRatio(jet.energy, jet.energyMax);
    if (jetCandidate.et != 0) {
      jetCandidates.push_back(jetCandidate);
      ht += jet.energy;

      // phiCenter is the 5-degree tower-bin index in the offset grid.  Use
      // its bin centre and undo the grid rotation to recover global phi.
      const double jetPhi = normalizedPhi(
          (static_cast<double>(jet.phiCenter) + 0.5) * 2. * kPi / static_cast<double>(kTowerPhiBins),
          -phiOffset_);
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

  putCandidates(links[kHgcalLinkBase], 0, egCandidates, packGamma);
  putCandidates(links[kHgcalLinkBase], kObjectsPerCategory, isoEgCandidates, packGamma);
  putCandidates(links[kHgcalLinkBase + 1], 0, jetCandidates, packHadron);
  putCandidates(links[kHgcalLinkBase + 1], kObjectsPerCategory, tauCandidates, packHadron);

  links[kHgcalLinkBase + 2] = p2mhhlink::packSumsV2(quantizeSigned(ex, etLsb_),
                                                    quantizeSigned(ey, etLsb_),
                                                    quantizeSigned(htx, etLsb_),
                                                    quantizeSigned(hty, etLsb_),
                                                    quantizeUnsigned(sumEt, etLsb_, 0xFFF),
                                                    quantizeUnsigned(ht, etLsb_, 0xFFF),
                                                    std::min(nObjects, 0xFFFU));

  if (debug_) {
    edm::LogVerbatim("Phase2L1DMXEmulator")
        << "event=" << event.id().event() << " zSide=" << zSide_ << " EG=" << egCandidates.size()
        << " isoEG=" << isoEgCandidates.size() << " jet=" << jetCandidates.size()
        << " tau=" << tauCandidates.size() << " towers=" << nObjects << " Ex=" << ex << " Ey=" << ey
        << " HTx=" << htx << " HTy=" << hty << " ET=" << sumEt << " HT=" << ht;
  }

  for (int link = 0; link < p2mhh::MHH_N_INPUT_LINKS; ++link) {
    auto words = std::make_unique<std::vector<uint64_t> >(toWords(links[link]));
    event.put(std::move(words), std::string("LinkIn") + std::to_string(link));
  }
}

void Phase2L1DMXEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("hgcalTowers", edm::InputTag("l1tHGCalTowerProducer", "HGCalTowerProcessor"));
  description.add<edm::InputTag>(
      "hgcalMulticlusters",
      edm::InputTag("l1tHGCalBackEndLayer2Producer", "HGCalBackendLayer2Processor3DClustering"));
  description.add<int>("zSide", 1);
  description.add<double>("etaMin", 1.479);
  description.add<double>("etaMax", 3.045);
  description.add<double>("etLsb", 0.5);
  description.add<double>("phiOffset", 0.0);
  description.add<double>("minTowerEt", 1.0);
  description.add<double>("minEgEt", 0.0);
  description.add<double>("maxEgHoverE", 0.3);
  description.add<double>("jetSeedThreshold", 3.0);
  description.add<bool>("duplicateEGAsIso", false);
  description.add<bool>("debug", false);
  descriptions.add("l1tPhase2L1DMXEmulator", description);
}

DEFINE_FWK_MODULE(Phase2L1DMXEmulator);
