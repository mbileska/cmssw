#include <ap_int.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/GCT_output.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/RCT_output.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TTree.h"

// Companion to MHHGCTSTrainingNtuplizer. Each event contains all RCT and
// GCTB Pre/Post frames plus the simple decoded energy fields below. The
// frame ordering is documented in DIST_TRIGGER_POC.md.
class Phase2L1FullCaloStageNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Phase2L1FullCaloStageNtuplizer(const edm::ParameterSet& config) {
    usesResource("TFileService");
    rct_ = tokens<l1tp2::rctOutputLinkCollection, 4>(config, "rctLinks");
    preIP1_ = tokens<l1tp2::rctOutputLinkCollection, 6>(config, "preIP1Links");
    postIP1_ = tokens<l1tp2::gctOutputLinkCollection, 6>(config, "postIP1Links");
    preIP2_ = tokens<l1tp2::gctOutputLinkCollection, 3>(config, "preIP2Links");
    postIP2_ = tokens<l1tp2::gctOutputLinkCollection, 3>(config, "postIP2Links");
    adapterFlags_ = tokens<std::vector<uint64_t>, 3>(config, "adapterFlags");

    edm::Service<TFileService> file;
    tree_ = file->make<TTree>("stages", "RCT and GCTB pipeline frames and decoded ET");
    tree_->Branch("run", &run_, "run/i");
    tree_->Branch("lumi", &lumi_, "lumi/i");
    tree_->Branch("event", &event_, "event/l");
    tree_->Branch("rctWords", &rctWords_);
    tree_->Branch("gctbPreIP1Words", &preIP1Words_);
    tree_->Branch("gctbPostIP1Words", &postIP1Words_);
    tree_->Branch("gctbPreIP2Words", &preIP2Words_);
    tree_->Branch("gctbPostIP2Words", &postIP2Words_);
    tree_->Branch("rctTowerEtGeV", &rctTowerEtGeV_);
    tree_->Branch("gctbSTEtGeV", &gctbSTEtGeV_);
    tree_->Branch("gctbEGPtGeV", &gctbEGPtGeV_);
    tree_->Branch("gctbJetPtGeV", &gctbJetPtGeV_);
    tree_->Branch("gctbTauPtGeV", &gctbTauPtGeV_);
    tree_->Branch("gctbAdapterFlags", &flags_);
  }

private:
  template <typename Collection, std::size_t N>
  std::array<edm::EDGetTokenT<Collection>, N> tokens(const edm::ParameterSet& config, const char* name) {
    const auto tags = config.getParameter<std::vector<edm::InputTag>>(name);
    if (tags.size() != N)
      throw cms::Exception("FullCaloStageNtuplizer") << name << " must have " << N << " tags";
    std::array<edm::EDGetTokenT<Collection>, N> result;
    for (std::size_t i = 0; i < N; ++i)
      result[i] = consumes<Collection>(tags[i]);
    return result;
  }

  template <typename Collection, std::size_t N>
  void read(edm::Event const& event, const std::array<edm::EDGetTokenT<Collection>, N>& input,
            std::size_t framesPerCollection, std::vector<uint64_t>& destination) {
    destination.clear();
    destination.reserve(N * framesPerCollection * 9);
    for (const auto& token : input) {
      const auto collection = event.getHandle(token);
      if (!collection.isValid() || collection->size() != framesPerCollection)
        throw cms::Exception("FullCaloStageNtuplizer") << "Missing or wrong-sized stage collection";
      for (const auto& frame : *collection) {
        const ap_uint<576> data = frame.data();
        for (unsigned word = 0; word < 9; ++word)
          destination.push_back(data.range(64 * word + 63, 64 * word).to_uint64());
      }
    }
  }

  void analyze(edm::Event const& event, edm::EventSetup const&) override {
    run_ = event.id().run();
    lumi_ = event.id().luminosityBlock();
    event_ = event.id().event();
    read(event, rct_, 24, rctWords_);
    read(event, preIP1_, 32, preIP1Words_);
    read(event, postIP1_, 21, postIP1Words_);
    read(event, preIP2_, 24, preIP2Words_);
    read(event, postIP2_, 6, postIP2Words_);

    flags_.clear();
    for (const auto& token : adapterFlags_) {
      const auto flags = event.getHandle(token);
      if (!flags.isValid() || flags->size() != 1)
        throw cms::Exception("FullCaloStageNtuplizer") << "Missing GCTB adapter assumption flags";
      flags_.push_back((*flags)[0]);
    }

    // RCT link 0 contains EG clusters. Links 1..3 contain 34 tower records
    // each (10-bit ET in the low bits of each 16-bit record).
    rctTowerEtGeV_.clear();
    rctTowerEtGeV_.reserve(24 * 102);
    for (unsigned card = 0; card < 24; ++card)
      for (unsigned link = 1; link < 4; ++link)
        for (unsigned slot = 0; slot < 34; ++slot) {
          const unsigned bit = 16 * slot;
          const unsigned index = (link * 24 + card) * 9 + bit / 64;
          const unsigned shift = bit % 64;
          // 16-bit records never straddle 64-bit transport words.
          const uint64_t record = (rctWords_[index] >> shift) & 0xFFFFU;
          rctTowerEtGeV_.push_back(0.5f * static_cast<float>(record & 0x3FFU));
        }

    // Core GCTB STs occupy the four central PreIP2 links on each eta side.
    constexpr unsigned positive[4] = {9, 10, 17, 18};
    constexpr unsigned negative[4] = {21, 22, 13, 14};
    gctbSTEtGeV_.clear();
    gctbSTEtGeV_.reserve(3 * 2 * 48);
    for (unsigned card = 0; card < 3; ++card)
      for (unsigned side = 0; side < 2; ++side)
        for (unsigned eta = 0; eta < 6; ++eta)
          for (unsigned phi = 0; phi < 8; ++phi) {
            const unsigned link = side == 0 ? positive[phi / 2] : negative[phi / 2];
            const unsigned slot = eta + 6 * (phi % 2);
            const unsigned frame = card * 24 + link;
            gctbSTEtGeV_.push_back(0.5f * static_cast<float>(recordPt(preIP2Words_, frame, slot)));
          }

    gctbEGPtGeV_.clear();
    gctbJetPtGeV_.clear();
    gctbTauPtGeV_.clear();
    for (unsigned card = 0; card < 3; ++card)
      for (unsigned side = 0; side < 2; ++side) {
        const unsigned egFrame = card * 6 + (side == 0 ? 0 : 3);
        const unsigned hadFrame = card * 6 + (side == 0 ? 1 : 4);
        for (unsigned slot = 0; slot < 6; ++slot) {
          gctbEGPtGeV_.push_back(0.5f * recordPt(postIP2Words_, egFrame, slot));
          gctbJetPtGeV_.push_back(0.5f * recordPt(postIP2Words_, hadFrame, slot));
          gctbTauPtGeV_.push_back(0.5f * recordPt(postIP2Words_, hadFrame, slot + 6));
        }
      }
    tree_->Fill();
  }

  static unsigned recordPt(const std::vector<uint64_t>& words, unsigned frame, unsigned slot) {
    const unsigned bit = slot * 48;
    const unsigned index = frame * 9 + bit / 64;
    const unsigned shift = bit % 64;
    if (shift <= 52)
      return (words[index] >> shift) & 0xFFFU;
    return ((words[index] >> shift) | (words[index + 1] << (64 - shift))) & 0xFFFU;
  }

  std::array<edm::EDGetTokenT<l1tp2::rctOutputLinkCollection>, 4> rct_;
  std::array<edm::EDGetTokenT<l1tp2::rctOutputLinkCollection>, 6> preIP1_;
  std::array<edm::EDGetTokenT<l1tp2::gctOutputLinkCollection>, 6> postIP1_;
  std::array<edm::EDGetTokenT<l1tp2::gctOutputLinkCollection>, 3> preIP2_, postIP2_;
  std::array<edm::EDGetTokenT<std::vector<uint64_t>>, 3> adapterFlags_;
  TTree* tree_{};
  unsigned int run_{}, lumi_{};
  unsigned long long event_{};
  std::vector<uint64_t> rctWords_, preIP1Words_, postIP1Words_, preIP2Words_, postIP2Words_, flags_;
  std::vector<float> rctTowerEtGeV_, gctbSTEtGeV_, gctbEGPtGeV_, gctbJetPtGeV_, gctbTauPtGeV_;
};

DEFINE_FWK_MODULE(Phase2L1FullCaloStageNtuplizer);
