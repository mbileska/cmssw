#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

namespace {
  enum class Process : uint8_t { qcd = 0, wjets = 1, zjets = 2, ttbar = 3 };
  enum class WDecay : uint8_t { electron = 1, muon = 2, tau = 3 };
  enum class ZDecay : uint8_t { electron = 1, muon = 2, tau = 3, neutrino = 4 };
  enum class TopWDecay : uint8_t { hadronic = 0, electron = 1, muon = 2, tau = 3 };

  Process parseProcess(const std::string& value) {
    if (value == "qcd")
      return Process::qcd;
    if (value == "wjets")
      return Process::wjets;
    if (value == "zjets")
      return Process::zjets;
    if (value == "ttbar")
      return Process::ttbar;
    throw cms::Exception("Configuration") << "Unknown background process '" << value
                                           << "'; choose qcd, wjets, zjets, or ttbar";
  }

  bool hasDaughterWithPdgId(const reco::Candidate& particle, int pdgId) {
    for (size_t index = 0; index < particle.numberOfDaughters(); ++index) {
      if (particle.daughter(index)->pdgId() == pdgId)
        return true;
    }
    return false;
  }

  const reco::Candidate* finalCopy(const reco::Candidate* particle) {
    const reco::Candidate* current = particle;
    while (current != nullptr) {
      const reco::Candidate* next = nullptr;
      for (size_t index = 0; index < current->numberOfDaughters(); ++index) {
        const auto* daughter = current->daughter(index);
        if (daughter->pdgId() == current->pdgId()) {
          next = daughter;
          break;
        }
      }
      if (next == nullptr)
        return current;
      current = next;
    }
    return particle;
  }

  bool isHardProcessCandidate(const reco::GenParticle& particle) {
    return particle.isHardProcess() || particle.fromHardProcessBeforeFSR() || particle.statusFlags().fromHardProcess();
  }

  struct LHELeptons {
    std::array<std::vector<int>, 3> charged;
    std::array<std::vector<int>, 3> neutrinos;
  };

  LHELeptons lheLeptons(const LHEEventProduct& lhe) {
    LHELeptons result;
    const auto& particles = lhe.hepeup();
    for (int index = 0; index < particles.NUP; ++index) {
      if (particles.ISTUP[index] != 1)
        continue;
      const int id = particles.IDUP[index];
      const int absId = std::abs(id);
      if (absId == 11 || absId == 13 || absId == 15)
        result.charged[(absId - 11) / 2].push_back(id);
      else if (absId == 12 || absId == 14 || absId == 16)
        result.neutrinos[(absId - 12) / 2].push_back(id);
    }
    return result;
  }

  WDecay lheWDecay(const LHEEventProduct& lhe) {
    const auto leptons = lheLeptons(lhe);
    unsigned int mode = 3;
    for (unsigned int index = 0; index < 3; ++index) {
      if (leptons.charged[index].size() == 1 && leptons.neutrinos[index].size() == 1 &&
          leptons.charged[index][0] * leptons.neutrinos[index][0] < 0) {
        if (mode != 3)
          throw cms::Exception("LogicError") << "More than one direct W decay mode in outgoing LHE particles";
        mode = index;
      } else if (!leptons.charged[index].empty() || !leptons.neutrinos[index].empty()) {
        throw cms::Exception("LogicError")
            << "Outgoing LHE leptons are not exactly one charge-compatible l-nu pair from a direct W decay";
      }
    }
    if (mode == 3)
      throw cms::Exception("LogicError") << "No direct W decay found in outgoing LHE particles";
    return static_cast<WDecay>(mode + 1);
  }

  ZDecay lheZDecay(const LHEEventProduct& lhe) {
    const auto leptons = lheLeptons(lhe);
    unsigned int mode = 4;
    for (unsigned int index = 0; index < 3; ++index) {
      const bool chargedPair =
          leptons.charged[index].size() == 2 && leptons.charged[index][0] == -leptons.charged[index][1];
      const bool neutrinoPair =
          leptons.neutrinos[index].size() == 2 && leptons.neutrinos[index][0] == -leptons.neutrinos[index][1];
      if (chargedPair || neutrinoPair) {
        if (mode != 4 || (chargedPair && neutrinoPair))
          throw cms::Exception("LogicError") << "More than one direct Z decay mode in outgoing LHE particles";
        mode = chargedPair ? index : 3;
      } else if (!leptons.charged[index].empty() || !leptons.neutrinos[index].empty()) {
        throw cms::Exception("LogicError")
            << "Outgoing LHE leptons are not exactly one opposite-sign same-flavor pair from a direct Z decay";
      }
    }
    if (mode == 4)
      throw cms::Exception("LogicError") << "No direct Z decay found in outgoing LHE particles";
    return static_cast<ZDecay>(mode + 1);
  }

  TopWDecay topWDecay(const reco::Candidate& inputW) {
    const auto& w = *finalCopy(&inputW);
    bool quark = false;
    for (size_t index = 0; index < w.numberOfDaughters(); ++index) {
      const int absId = std::abs(w.daughter(index)->pdgId());
      if (absId == 11)
        return TopWDecay::electron;
      if (absId == 13)
        return TopWDecay::muon;
      if (absId == 15)
        return TopWDecay::tau;
      quark |= absId >= 1 && absId <= 6;
    }
    if (quark)
      return TopWDecay::hadronic;
    throw cms::Exception("LogicError") << "Could not assign the direct decay of a top-quark W";
  }

  unsigned int binIndex(double value, const std::vector<double>& upperEdges) {
    if (!std::isfinite(value) || value < 0.)
      throw cms::Exception("LogicError") << "Invalid non-negative kinematic value " << value;
    return std::distance(upperEdges.begin(), std::upper_bound(upperEdges.begin(), upperEdges.end(), value));
  }

  double lheHT(const LHEEventProduct& lhe) {
    const auto& particles = lhe.hepeup();
    double result = 0.;
    for (int index = 0; index < particles.NUP; ++index) {
      const int absId = std::abs(particles.IDUP[index]);
      if (particles.ISTUP[index] == 1 && (absId <= 6 || absId == 21))
        result += std::hypot(particles.PUP[index][0], particles.PUP[index][1]);
    }
    return result;
  }

  std::array<TopWDecay, 2> topDecays(const reco::GenParticleCollection& particles) {
    std::vector<const reco::Candidate*> tops;
    for (const auto& particle : particles) {
      if (std::abs(particle.pdgId()) == 6 && isHardProcessCandidate(particle) &&
          !hasDaughterWithPdgId(particle, particle.pdgId())) {
        tops.push_back(&particle);
      }
    }
    if (tops.size() != 2)
      throw cms::Exception("LogicError") << "Expected exactly two final hard-process top quarks, found "
                                         << tops.size();
    if (tops[0]->pdgId() != -tops[1]->pdgId())
      throw cms::Exception("LogicError") << "The two final hard-process top quarks are not a top-antitop pair";

    std::array<TopWDecay, 2> result;
    for (unsigned int topIndex = 0; topIndex < tops.size(); ++topIndex) {
      const reco::Candidate* w = nullptr;
      for (size_t daughterIndex = 0; daughterIndex < tops[topIndex]->numberOfDaughters(); ++daughterIndex) {
        const auto* daughter = tops[topIndex]->daughter(daughterIndex);
        if (std::abs(daughter->pdgId()) == 24) {
          w = daughter;
          break;
        }
      }
      if (w == nullptr)
        throw cms::Exception("LogicError") << "Final top quark has no direct W daughter";
      result[topIndex] = topWDecay(*w);
    }
    std::sort(result.begin(), result.end(), [](TopWDecay left, TopWDecay right) {
      return static_cast<unsigned int>(left) < static_cast<unsigned int>(right);
    });
    return result;
  }

  uint8_t topPairIndex(const std::array<TopWDecay, 2>& decays) {
    // HH, HE, HM, HT, EE, EM, ET, MM, MT, TT
    static constexpr uint8_t lookup[4][4] = {{0, 1, 2, 3}, {1, 4, 5, 6}, {2, 5, 7, 8}, {3, 6, 8, 9}};
    return lookup[static_cast<unsigned int>(decays[0])][static_cast<unsigned int>(decays[1])];
  }
}  // namespace

class BackgroundCategoryTableProducer : public edm::global::EDProducer<> {
public:
  explicit BackgroundCategoryTableProducer(const edm::ParameterSet& config)
      : process_(parseProcess(config.getParameter<std::string>("process"))),
        sample_(config.getParameter<unsigned int>("sample")),
        genParticlesToken_(
            consumes<reco::GenParticleCollection>(config.getParameter<edm::InputTag>("genParticles"))),
        generatorToken_(consumes<GenEventInfoProduct>(config.getParameter<edm::InputTag>("generator"))),
        lheToken_(consumes<LHEEventProduct>(config.getParameter<edm::InputTag>("lheEvent"))) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& event, const edm::EventSetup&) const override {
    uint8_t decay = 0;
    uint8_t kinematicBin = 0;
    uint8_t category = 0;
    float kinematicValue = -1.f;

    if (process_ == Process::qcd) {
      const auto& generator = event.get(generatorToken_);
      if (!generator.hasBinningValues())
        throw cms::Exception("LogicError")
            << "QCD category requires generator.binningValues()[0] as the generated hard-pT value";
      kinematicValue = generator.binningValues().front();
      if (kinematicValue < 15.f)
        throw cms::Exception("LogicError") << "QCD generated hard pT is below the declared 15 GeV phase-space floor";
      kinematicBin = binIndex(kinematicValue, {50., 100., 200., 400., 800.});
      category = kinematicBin;
    } else if (process_ == Process::wjets) {
      const auto lhe = event.getHandle(lheToken_);
      if (!lhe.isValid())
        throw cms::Exception("ProductNotFound") << "W+jets category requires externalLHEProducer for LHE HT";
      kinematicValue = lheHT(*lhe);
      kinematicBin = binIndex(kinematicValue, {100., 400.});
      decay = static_cast<uint8_t>(lheWDecay(*lhe));
      category = 6 + (decay - 1) * 3 + kinematicBin;
    } else if (process_ == Process::zjets) {
      const auto lhe = event.getHandle(lheToken_);
      if (!lhe.isValid())
        throw cms::Exception("ProductNotFound") << "Z/DY+jets category requires externalLHEProducer for LHE HT";
      kinematicValue = lheHT(*lhe);
      kinematicBin = binIndex(kinematicValue, {100., 400.});
      decay = static_cast<uint8_t>(lheZDecay(*lhe));
      category = 15 + (decay - 1) * 3 + kinematicBin;
    } else {
      decay = topPairIndex(topDecays(event.get(genParticlesToken_)));
      category = 27 + decay;
    }

    auto table = std::make_unique<nanoaod::FlatTable>(1, "Background", true);
    table->addColumnValue<uint8_t>("process", static_cast<uint8_t>(process_), "0=QCD, 1=W+jets, 2=Z/DY+jets, 3=ttbar");
    table->addColumnValue<uint32_t>("sample", sample_, "parent sample ID from background_samples_2025.json");
    table->addColumnValue<uint8_t>("decay", decay, "process-specific exclusive direct-decay index");
    table->addColumnValue<uint8_t>("kinematicBin", kinematicBin, "process-specific exclusive kinematic bin");
    table->addColumnValue<uint8_t>("category", category, "exclusive category ID from background_categories.json");
    table->addColumnValue<float>(
        "kinematicValue", kinematicValue, "QCD generated hard pT or W/Z LHE HT in GeV; -1 for ttbar");
    event.put(std::move(table));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<std::string>("process")->setComment("qcd, wjets, zjets, or ttbar; never inferred from event content");
    desc.add<unsigned int>("sample")->setComment("parent sample ID used to keep normalization provenance");
    desc.add<edm::InputTag>("genParticles", edm::InputTag("genParticles"));
    desc.add<edm::InputTag>("generator", edm::InputTag("generator"));
    desc.add<edm::InputTag>("lheEvent", edm::InputTag("externalLHEProducer"));
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const Process process_;
  const uint32_t sample_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BackgroundCategoryTableProducer);
