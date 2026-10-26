/*
 * Description:
 *   CMSSW EDProducer for the complete Phase-2 GCT IP1 -> IP2 emulator chain.
 *
 *   Event flow:
 *     1. Consume the four RCT output-link collections.
 *     2. Build six 32-link pre-IP1 SLR bundles in the Figure-6 RCT-pair
 *        label space.
 *     3. Run the existing IP1 firmware independently for all six SLR bundles.
 *     4. Route each SLR1/SLR3 pair into one 24-link pre-IP2 GCT-card bundle.
 *     5. Run the supplied IP2 firmware for GCT1, GCT2, and GCT3.
 *     6. Save every boundary: PreIP1, PostIP1, PreIP2, and PostIP2.
 *
 *   The RCT-card ordering and the IP1-to-IP2 router are deliberately explicit.
 *   IP1 regions 0..3 are the four positive-eta RCT collection entries in the
 *   Figure-6 pair order.  Regions 4..7 are the matching negative-eta entries in
 *   reverse order, as required by IP1's R[i] <-> R[7-i] boundary logic.
 *
 *   IMPORTANT: the current RCT producer partitions phi from -pi in 30-degree
 *   blocks, while Figure 6 places the physical RCT-card boundaries on a
 *   different 30-degree grid.  Therefore the identity between a collection
 *   pair index and a physical Figure-6 card-pair label is still unverified.
 *   See SOURCE_BASIS_AND_LIMITS.md in the delivered overlay.
 *
 *   EG inputs use IP1's dedicated links 9..12.  Supertower inputs retain the
 *   left overlap card, four central cards, and right overlap card in the PDF
 *   pair-label order.  All 24 attached-firmware IP2 input positions are
 *   documented in buildIP2Input().
 *
 *   One adapter is required between the checked-in IP1 and attached IP2: IP1
 *   serializes each 12-ST word as [eta0 phi0, eta0 phi1, eta1 phi0, ...], while
 *   IP2 fillInputGCTst() reads [eta0..eta5 phi0, eta0..eta5 phi1].  The helper
 *   transposeIP1STForIP2() performs only that slot permutation.  PostIP1 is
 *   saved before the permutation and PreIP2 is saved after it, so both sides of
 *   the interface can be compared directly.
 */

#include <ap_int.h>

#include <array>
#include <cstddef>
#include <memory>
#include <string>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/L1TCalorimeterPhase2/interface/GCT_output.h"
#include "DataFormats/L1TCalorimeterPhase2/interface/RCT_output.h"

// Existing IP1 firmware implementation.
#include "L1Trigger/L1CaloTrigger/interface/algo_topIP1gct_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort16_h.h"
#include "L1Trigger/L1CaloTrigger/interface/pfcluster_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/algo_topIP1gct_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort16_cpp.h"

// Namespaced IP2 firmware implementation added by this integration.
#include "L1Trigger/L1CaloTrigger/interface/algo_topIP2gct_h.h"
#include "L1Trigger/L1CaloTrigger/interface/gctforgtIP2_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/algo_topIP2gct_cpp.h"

class Phase2L1CaloL1GCTEmulator : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1CaloL1GCTEmulator(const edm::ParameterSet&);
  ~Phase2L1CaloL1GCTEmulator() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  using RCTCollection = l1tp2::rctOutputLinkCollection;
  using GCTCollection = l1tp2::gctOutputLinkCollection;
  using LinkWord = ap_uint<576>;

  static constexpr int kNRCTCards = 24;
  static constexpr int kNGCTSLRs = 6;
  static constexpr int kNGCTCards = 3;
  static constexpr int kLinksPerRCTCard = 4;
  static constexpr int kLinksPerIP1Input = 32;
  static constexpr int kIP1OutputLinks = N_OUTPUT_LINKS;

  using IP1Input = std::array<LinkWord, kLinksPerIP1Input>;
  using IP1Output = std::array<LinkWord, kIP1OutputLinks>;
  using IP2Input = std::array<LinkWord, gctip2::kInputLinks>;
  using IP2Output = std::array<LinkWord, gctip2::kOutputLinks>;

  edm::EDGetTokenT<RCTCollection> link0Src_;
  edm::EDGetTokenT<RCTCollection> link1Src_;
  edm::EDGetTokenT<RCTCollection> link2Src_;
  edm::EDGetTokenT<RCTCollection> link3Src_;

  LinkWord getRCTLink(const RCTCollection& link0,
                      const RCTCollection& link1,
                      const RCTCollection& link2,
                      const RCTCollection& link3,
                      int cardIndex,
                      int localLink) const;
  void appendRCTCard(IP1Input& output,
                     int& writeIndex,
                     const RCTCollection& link0,
                     const RCTCollection& link1,
                     const RCTCollection& link2,
                     const RCTCollection& link3,
                     int cardIndex) const;
  IP1Input buildIP1Input(const RCTCollection& link0,
                         const RCTCollection& link1,
                         const RCTCollection& link2,
                         const RCTCollection& link3,
                         int slrIndex) const;
  IP1Output runIP1(const IP1Input& input) const;
  static LinkWord transposeIP1STForIP2(const LinkWord& input);
  IP2Input buildIP2Input(const IP1Output& slr1, const IP1Output& slr3) const;
  IP2Output runIP2(const IP2Input& input) const;

  static int pairToNegativeCard(int pairIndex) { return 2 * pairIndex; }
  static int pairToPositiveCard(int pairIndex) { return 2 * pairIndex + 1; }
};

Phase2L1CaloL1GCTEmulator::Phase2L1CaloL1GCTEmulator(const edm::ParameterSet& config)
    : link0Src_(consumes<RCTCollection>(config.getParameter<edm::InputTag>("LinkOut0"))),
      link1Src_(consumes<RCTCollection>(config.getParameter<edm::InputTag>("LinkOut1"))),
      link2Src_(consumes<RCTCollection>(config.getParameter<edm::InputTag>("LinkOut2"))),
      link3Src_(consumes<RCTCollection>(config.getParameter<edm::InputTag>("LinkOut3"))) {
  static_assert(kLinksPerIP1Input == N_INPUT_LINKS,
                "The CMSSW IP1 router must provide exactly the firmware N_INPUT_LINKS.");
  static_assert(kIP1OutputLinks >= 21, "IP2 routing requires all 21 IP1 output links.");
  static_assert(gctip2::kInputLinks == 24, "The supplied IP2 firmware expects 24 input links.");
  static_assert(gctip2::kOutputLinks == 6, "The supplied IP2 firmware emits 6 output links.");

  const std::array<std::string, kNGCTSLRs> preIP1Names{{
      "GCT1SLR3PreIP1", "GCT1SLR1PreIP1", "GCT2SLR3PreIP1",
      "GCT2SLR1PreIP1", "GCT3SLR3PreIP1", "GCT3SLR1PreIP1"}};
  const std::array<std::string, kNGCTSLRs> postIP1Names{{
      "GCT1SLR3PostIP1", "GCT1SLR1PostIP1", "GCT2SLR3PostIP1",
      "GCT2SLR1PostIP1", "GCT3SLR3PostIP1", "GCT3SLR1PostIP1"}};
  const std::array<std::string, kNGCTCards> preIP2Names{{"GCT1PreIP2", "GCT2PreIP2", "GCT3PreIP2"}};
  const std::array<std::string, kNGCTCards> postIP2Names{{"GCT1PostIP2", "GCT2PostIP2", "GCT3PostIP2"}};

  for (const auto& name : preIP1Names) {
    produces<RCTCollection>(name);
  }
  for (const auto& name : postIP1Names) {
    produces<GCTCollection>(name);
  }
  for (const auto& name : preIP2Names) {
    produces<GCTCollection>(name);
  }
  for (const auto& name : postIP2Names) {
    produces<GCTCollection>(name);
  }
}

Phase2L1CaloL1GCTEmulator::LinkWord Phase2L1CaloL1GCTEmulator::getRCTLink(
    const RCTCollection& link0,
    const RCTCollection& link1,
    const RCTCollection& link2,
    const RCTCollection& link3,
    int cardIndex,
    int localLink) const {
  if (cardIndex < 0 || cardIndex >= kNRCTCards) {
    throw cms::Exception("Phase2L1CaloL1GCTEmulator") << "Invalid RCT card index " << cardIndex;
  }
  if (static_cast<int>(link0.size()) <= cardIndex || static_cast<int>(link1.size()) <= cardIndex ||
      static_cast<int>(link2.size()) <= cardIndex || static_cast<int>(link3.size()) <= cardIndex) {
    throw cms::Exception("Phase2L1CaloL1GCTEmulator")
        << "RCT output collections are too short for card " << cardIndex << ". Sizes are "
        << link0.size() << ", " << link1.size() << ", " << link2.size() << ", " << link3.size() << '.';
  }

  switch (localLink) {
    case 0:
      return link0[cardIndex].data();
    case 1:
      return link1[cardIndex].data();
    case 2:
      return link2[cardIndex].data();
    case 3:
      return link3[cardIndex].data();
    default:
      throw cms::Exception("Phase2L1CaloL1GCTEmulator") << "Invalid local RCT link " << localLink;
  }
}

void Phase2L1CaloL1GCTEmulator::appendRCTCard(IP1Input& output,
                                              int& writeIndex,
                                              const RCTCollection& link0,
                                              const RCTCollection& link1,
                                              const RCTCollection& link2,
                                              const RCTCollection& link3,
                                              int cardIndex) const {
  for (int localLink = 0; localLink < kLinksPerRCTCard; ++localLink) {
    output[writeIndex++] = getRCTLink(link0, link1, link2, link3, cardIndex, localLink);
  }
}

Phase2L1CaloL1GCTEmulator::IP1Input Phase2L1CaloL1GCTEmulator::buildIP1Input(
    const RCTCollection& link0,
    const RCTCollection& link1,
    const RCTCollection& link2,
    const RCTCollection& link3,
    int slrIndex) const {
  // RCT eta-pair labels for the six overlapping 120-degree SLR input
  // windows in Figure 6 of the L1 Calorimeter interface note.  Each row
  // contains the left overlap pair, two core pairs, and the right overlap pair
  // in the PDF's increasing-phi label order.
  //
  // The current event collections use pair p = entries (2*p, 2*p+1), with even
  // entries treated as negative eta and odd entries as positive eta.  Mapping
  // collection pair p directly to PDF physical pair label p is an explicit,
  // provisional assumption: the current RCT producer's phi boundaries do not
  // yet reproduce Figure 6.  This table therefore fixes the GCT combinatorics
  // but does not by itself establish end-to-end physical card alignment.
  static constexpr int kRegionPairMap[kNGCTSLRs][4] = {
      {11, 0, 1, 2},   // Region 0: GCT1 SLR3
      {1, 2, 3, 4},    // Region 1: GCT1 SLR1
      {3, 4, 5, 6},    // Region 2: GCT2 SLR3
      {5, 6, 7, 8},    // Region 3: GCT2 SLR1
      {7, 8, 9, 10},   // Region 4: GCT3 SLR3
      {9, 10, 11, 0}   // Region 5: GCT3 SLR1
  };

  if (slrIndex < 0 || slrIndex >= kNGCTSLRs) {
    throw cms::Exception("Phase2L1CaloL1GCTEmulator") << "Invalid GCT SLR index " << slrIndex;
  }

  IP1Input output{};
  int writeIndex = 0;

  // The IP1 firmware treats RCTRegions[0..3] as one eta side and
  // RCTRegions[4..7] as the other.  It joins adjacent regions within each
  // group and explicitly pairs RCTRegions[i] with RCTRegions[7-i] across eta.
  // Consequently the two eta sides must not be interleaved.
  //
  // Positive eta is stored first in the PDF pair-label order.  Negative eta is
  // stored second in reverse pair-label order so region i and region 7-i use
  // the same collection pair when IP1 mirrors the eta boundary.
  for (int phiSlot = 0; phiSlot < 4; ++phiSlot) {
    const int pair = kRegionPairMap[slrIndex][phiSlot];
    appendRCTCard(output, writeIndex, link0, link1, link2, link3, pairToPositiveCard(pair));
  }
  for (int phiSlot = 3; phiSlot >= 0; --phiSlot) {
    const int pair = kRegionPairMap[slrIndex][phiSlot];
    appendRCTCard(output, writeIndex, link0, link1, link2, link3, pairToNegativeCard(pair));
  }

  if (writeIndex != kLinksPerIP1Input) {
    throw cms::Exception("Phase2L1CaloL1GCTEmulator")
        << "Internal IP1 routing error: wrote " << writeIndex << " links instead of " << kLinksPerIP1Input;
  }
  return output;
}

Phase2L1CaloL1GCTEmulator::IP1Output Phase2L1CaloL1GCTEmulator::runIP1(const IP1Input& input) const {
  ap_uint<576> firmwareInput[N_INPUT_LINKS];
  ap_uint<576> firmwareOutput[N_OUTPUT_LINKS];
  for (int index = 0; index < N_INPUT_LINKS; ++index) {
    firmwareInput[index] = input[index];
  }
  for (int index = 0; index < N_OUTPUT_LINKS; ++index) {
    firmwareOutput[index] = 0;
  }

  algo_topIP1gct(firmwareInput, firmwareOutput);

  IP1Output output{};
  for (int index = 0; index < N_OUTPUT_LINKS; ++index) {
    output[index] = firmwareOutput[index];
  }
  return output;
}

Phase2L1CaloL1GCTEmulator::LinkWord Phase2L1CaloL1GCTEmulator::transposeIP1STForIP2(
    const LinkWord& input) {
  // Existing IP1 fillSTlink():
  //   source slot = 2*eta + phi, eta=0..5, phi=0..1.
  // Attached IP2 fillInputGCTst():
  //   destination slot = eta + 6*phi.
  //
  // Every slot is copied as an opaque 48-bit object.  No energy, coordinate,
  // flag, or precision field is changed here.
  LinkWord output = 0;
  for (int eta = 0; eta < 6; ++eta) {
    for (int phi = 0; phi < 2; ++phi) {
      const int sourceSlot = 2 * eta + phi;
      const int destinationSlot = eta + 6 * phi;
      output.range(48 * destinationSlot + 47, 48 * destinationSlot) =
          input.range(48 * sourceSlot + 47, 48 * sourceSlot);
    }
  }
  return output;
}

Phase2L1CaloL1GCTEmulator::IP2Input Phase2L1CaloL1GCTEmulator::buildIP2Input(
    const IP1Output& slr1,
    const IP1Output& slr3) const {
  IP2Input output{};

  // --------------------------------------------------------------------------
  // EG routing
  // --------------------------------------------------------------------------
  // IP1 links 9..12 are the four dedicated, full 9-object EG links:
  //   link  9 = positive-eta local RCT region 1
  //   link 10 = positive-eta local RCT region 2
  //   link 11 = negative-eta local RCT region 5
  //   link 12 = negative-eta local RCT region 6
  //
  // The supplied IP2 top function reads positive eta in the order
  // inputs (0,1,4,5) and negative eta in the order (6,7,2,3).
  //
  // In the provisional PDF pair-label topology, SLR3 provides the lower-label
  // half of the four core RCT pairs and SLR1 provides the higher-label half.
  // This is not yet a proven detector-phi statement because the current RCT
  // producer partitions phi on a different boundary grid.  The positive
  // sequence is
  // therefore SLR3 regions 1,2 followed by SLR1 regions 1,2.  Because the
  // negative IP1 regions are stored in reverse phi order, the corresponding
  // negative sequence is SLR1 regions 5,6 followed by SLR3 regions 5,6.
  output[0] = slr3[9];
  output[1] = slr3[10];
  output[2] = slr3[11];
  output[3] = slr3[12];
  output[4] = slr1[9];
  output[5] = slr1[10];
  output[6] = slr1[11];
  output[7] = slr1[12];

  // --------------------------------------------------------------------------
  // Supertower routing
  // --------------------------------------------------------------------------
  // IP1 links 13..20 correspond directly to local RCT regions 0..7:
  //   links 13..16 = positive eta, increasing PDF pair-label slot;
  //   links 17..20 = negative eta, decreasing PDF pair-label slot.
  //
  // IP2 needs six RCT-card columns: left overlap + four central cards + right
  // overlap.  For positive eta this is SLR3 slots 0,1,2 followed by SLR1
  // slots 1,2,3.  The negative side is the exact reverse in this provisional
  // pair-label topology, as required by fillOverlap() in the supplied IP2 top
  // function.
  constexpr std::array<int, 3> kSLR3Positive{{13, 14, 15}};
  constexpr std::array<int, 3> kSLR1Positive{{14, 15, 16}};
  constexpr std::array<int, 3> kSLR1Negative{{17, 18, 19}};
  constexpr std::array<int, 3> kSLR3Negative{{18, 19, 20}};

  for (std::size_t index = 0; index < 3; ++index) {
    output[8 + index] = transposeIP1STForIP2(slr3[kSLR3Positive[index]]);
    output[17 + index] = transposeIP1STForIP2(slr1[kSLR1Positive[index]]);
    output[20 + index] = transposeIP1STForIP2(slr1[kSLR1Negative[index]]);
    output[13 + index] = transposeIP1STForIP2(slr3[kSLR3Negative[index]]);
  }

  // Attached-firmware inputs 11, 12, 16, and 23 are unused and remain zero.
  // PreIP2 stores this exact post-adapter 24-link array.
  return output;
}

Phase2L1CaloL1GCTEmulator::IP2Output Phase2L1CaloL1GCTEmulator::runIP2(const IP2Input& input) const {
  IP2Output output{};
  gctip2::algoTop(input, output);
  return output;
}

void Phase2L1CaloL1GCTEmulator::produce(edm::Event& event, const edm::EventSetup&) {
  const auto link0 = event.getHandle(link0Src_);
  const auto link1 = event.getHandle(link1Src_);
  const auto link2 = event.getHandle(link2Src_);
  const auto link3 = event.getHandle(link3Src_);
  if (!link0.isValid() || !link1.isValid() || !link2.isValid() || !link3.isValid()) {
    throw cms::Exception("Phase2L1CaloL1GCTEmulator") << "One or more RCT output-link products are missing.";
  }

  const std::array<std::string, kNGCTSLRs> preIP1Names{{
      "GCT1SLR3PreIP1", "GCT1SLR1PreIP1", "GCT2SLR3PreIP1",
      "GCT2SLR1PreIP1", "GCT3SLR3PreIP1", "GCT3SLR1PreIP1"}};
  const std::array<std::string, kNGCTSLRs> postIP1Names{{
      "GCT1SLR3PostIP1", "GCT1SLR1PostIP1", "GCT2SLR3PostIP1",
      "GCT2SLR1PostIP1", "GCT3SLR3PostIP1", "GCT3SLR1PostIP1"}};
  const std::array<std::string, kNGCTCards> preIP2Names{{"GCT1PreIP2", "GCT2PreIP2", "GCT3PreIP2"}};
  const std::array<std::string, kNGCTCards> postIP2Names{{"GCT1PostIP2", "GCT2PostIP2", "GCT3PostIP2"}};

  std::array<std::unique_ptr<RCTCollection>, kNGCTSLRs> preIP1Products;
  std::array<std::unique_ptr<GCTCollection>, kNGCTSLRs> postIP1Products;
  std::array<std::unique_ptr<GCTCollection>, kNGCTCards> preIP2Products;
  std::array<std::unique_ptr<GCTCollection>, kNGCTCards> postIP2Products;
  std::array<IP1Output, kNGCTSLRs> ip1Outputs{};

  for (int slr = 0; slr < kNGCTSLRs; ++slr) {
    preIP1Products[slr] = std::make_unique<RCTCollection>();
    postIP1Products[slr] = std::make_unique<GCTCollection>();

    const IP1Input ip1Input = buildIP1Input(*link0, *link1, *link2, *link3, slr);
    for (const LinkWord& word : ip1Input) {
      preIP1Products[slr]->emplace_back(word);
    }

    ip1Outputs[slr] = runIP1(ip1Input);
    for (const LinkWord& word : ip1Outputs[slr]) {
      postIP1Products[slr]->emplace_back(word);
    }
  }

  // Product ordering is [GCT1 SLR3, GCT1 SLR1, GCT2 SLR3, GCT2 SLR1, ...].
  // Therefore SLR3 = 2*gct and SLR1 = 2*gct+1.
  for (int gct = 0; gct < kNGCTCards; ++gct) {
    preIP2Products[gct] = std::make_unique<GCTCollection>();
    postIP2Products[gct] = std::make_unique<GCTCollection>();

    const IP1Output& slr3 = ip1Outputs[2 * gct];
    const IP1Output& slr1 = ip1Outputs[2 * gct + 1];
    const IP2Input ip2Input = buildIP2Input(slr1, slr3);
    for (const LinkWord& word : ip2Input) {
      preIP2Products[gct]->emplace_back(word);
    }

    const IP2Output ip2Output = runIP2(ip2Input);
    for (const LinkWord& word : ip2Output) {
      postIP2Products[gct]->emplace_back(word);
    }
  }

  edm::LogInfo("GCTEmulator")
      << "Produced per event: 6 x 32 PreIP1 links, 6 x " << kIP1OutputLinks
      << " PostIP1 links, 3 x " << gctip2::kInputLinks << " PreIP2 links, and 3 x "
      << gctip2::kOutputLinks << " PostIP2 links.";

  for (int slr = 0; slr < kNGCTSLRs; ++slr) {
    event.put(std::move(preIP1Products[slr]), preIP1Names[slr]);
    event.put(std::move(postIP1Products[slr]), postIP1Names[slr]);
  }
  for (int gct = 0; gct < kNGCTCards; ++gct) {
    event.put(std::move(preIP2Products[gct]), preIP2Names[gct]);
    event.put(std::move(postIP2Products[gct]), postIP2Names[gct]);
  }
}

void Phase2L1CaloL1GCTEmulator::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription description;
  description.add<edm::InputTag>("LinkOut0", edm::InputTag("l1tPhase2RCTEmulatorProducer", "LinkOut0"));
  description.add<edm::InputTag>("LinkOut1", edm::InputTag("l1tPhase2RCTEmulatorProducer", "LinkOut1"));
  description.add<edm::InputTag>("LinkOut2", edm::InputTag("l1tPhase2RCTEmulatorProducer", "LinkOut2"));
  description.add<edm::InputTag>("LinkOut3", edm::InputTag("l1tPhase2RCTEmulatorProducer", "LinkOut3"));
  descriptions.add("Phase2L1CaloL1GCTEmulator", description);
}

DEFINE_FWK_MODULE(Phase2L1CaloL1GCTEmulator);
