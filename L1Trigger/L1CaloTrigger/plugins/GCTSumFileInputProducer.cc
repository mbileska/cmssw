/*
 * Description: File-backed input producer for Phase 2 GCT SumCard emulator
 */

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace {

constexpr unsigned int kWordsPerEvent = 9;
constexpr unsigned int kLinksPerSide = 12;
constexpr unsigned int kInputLinks = 24;

using LinkWords = std::array<uint64_t, kWordsPerEvent>;
using EventWords = std::array<LinkWords, kInputLinks>;

std::string trimCopy(const std::string& input) {
  std::size_t begin = 0;
  while (begin < input.size() && std::isspace(static_cast<unsigned char>(input[begin])) != 0) {
    ++begin;
  }

  std::size_t end = input.size();
  while (end > begin && std::isspace(static_cast<unsigned char>(input[end - 1])) != 0) {
    --end;
  }

  return input.substr(begin, end - begin);
}

bool isDataLine(const std::string& line) {
  const std::string trimmed = trimCopy(line);
  if (trimmed.empty() || trimmed[0] == '#') {
    return false;
  }

  std::istringstream iss(trimmed);
  std::string firstToken;
  iss >> firstToken;
  if (firstToken.empty()) {
    return false;
  }

  std::size_t begin = 0;
  if (firstToken.size() > 2 && firstToken[0] == '0' && (firstToken[1] == 'x' || firstToken[1] == 'X')) {
    begin = 2;
  }
  if (begin >= firstToken.size()) {
    return false;
  }

  for (std::size_t i = begin; i < firstToken.size(); ++i) {
    if (std::isxdigit(static_cast<unsigned char>(firstToken[i])) == 0) {
      return false;
    }
  }
  return true;
}

uint64_t parseHexWord64(const std::string& token) {
  return std::stoull(token, nullptr, 0);
}

std::vector<EventWords> readEvents(const std::string& path, unsigned int posOffset, unsigned int negOffset) {
  std::ifstream input(path.c_str());
  if (!input.is_open()) {
    throw cms::Exception("GCTSumFileInputProducer") << "Could not open input vector file: " << path;
  }

  std::vector<std::array<uint64_t, kInputLinks> > rows;
  std::string line;
  const unsigned int minColumns = std::max(posOffset + kLinksPerSide, negOffset + kLinksPerSide);

  while (std::getline(input, line)) {
    if (!isDataLine(line)) {
      continue;
    }

    std::istringstream iss(line);
    std::string rowToken;
    iss >> rowToken;

    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token) {
      tokens.push_back(token);
    }

    if (tokens.size() < minColumns) {
      throw cms::Exception("GCTSumFileInputProducer")
          << "Vector row in " << path << " has " << tokens.size() << " payload columns, expected at least " << minColumns;
    }

    std::array<uint64_t, kInputLinks> row{};
    for (unsigned int i = 0; i < kLinksPerSide; ++i) {
      row[i] = parseHexWord64(tokens[posOffset + i]);
      row[kLinksPerSide + i] = parseHexWord64(tokens[negOffset + i]);
    }
    rows.push_back(row);
  }

  if (rows.empty()) {
    throw cms::Exception("GCTSumFileInputProducer") << "No data rows were found in input vector file: " << path;
  }

  if ((rows.size() % kWordsPerEvent) != 0) {
    throw cms::Exception("GCTSumFileInputProducer")
        << "Input vector rows (" << rows.size() << ") are not a multiple of " << kWordsPerEvent;
  }

  const std::size_t eventCount = rows.size() / kWordsPerEvent;
  std::vector<EventWords> events(eventCount);
  for (std::size_t event = 0; event < eventCount; ++event) {
    for (unsigned int link = 0; link < kInputLinks; ++link) {
      events[event][link].fill(0);
    }

    for (unsigned int row = 0; row < kWordsPerEvent; ++row) {
      const std::size_t rowIndex = event * kWordsPerEvent + row;
      for (unsigned int link = 0; link < kInputLinks; ++link) {
        events[event][link][row] = rows[rowIndex][link];
      }
    }
  }

  return events;
}

}  // namespace

class GCTSumFileInputProducer : public edm::stream::EDProducer<> {
public:
  explicit GCTSumFileInputProducer(const edm::ParameterSet&);
  ~GCTSumFileInputProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  std::vector<EventWords> events_;
  std::size_t eventOffset_;
  bool wrapAround_;
  bool debug_;
};

GCTSumFileInputProducer::GCTSumFileInputProducer(const edm::ParameterSet& iConfig)
    : events_(readEvents(iConfig.getParameter<std::string>("inputFile"),
                         iConfig.getParameter<unsigned int>("posOffset"),
                         iConfig.getParameter<unsigned int>("negOffset"))),
      eventOffset_(iConfig.getParameter<unsigned int>("eventOffset")),
      wrapAround_(iConfig.getParameter<bool>("wrapAround")),
      debug_(iConfig.getParameter<bool>("debug")) {
  for (unsigned int i = 0; i < kInputLinks; ++i) {
    produces<std::vector<uint64_t> >(std::string("LinkIn") + std::to_string(i));
  }
}

void GCTSumFileInputProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  (void)iSetup;

  if (events_.empty()) {
    throw cms::Exception("GCTSumFileInputProducer") << "No events were loaded from the configured input file";
  }

  const unsigned long long eventNumber = iEvent.id().event();
  if (eventNumber == 0) {
    throw cms::Exception("GCTSumFileInputProducer") << "Expected event numbering to start at 1";
  }

  std::size_t eventIndex = eventOffset_ + static_cast<std::size_t>(eventNumber - 1ULL);
  if (wrapAround_) {
    eventIndex %= events_.size();
  } else if (eventIndex >= events_.size()) {
    throw cms::Exception("GCTSumFileInputProducer")
        << "Requested event index " << eventIndex << " but only " << events_.size()
        << " events are available in the input file. Set maxEvents accordingly or enable wrapAround.";
  }

  if (debug_) {
    edm::LogVerbatim("GCTSumFileInputProducer")
        << "Producing GCT Sum file event " << eventIndex << " for edm event " << eventNumber;
  }

  for (unsigned int link = 0; link < kInputLinks; ++link) {
    auto out = std::make_unique<std::vector<uint64_t> >(events_[eventIndex][link].begin(), events_[eventIndex][link].end());
    iEvent.put(std::move(out), std::string("LinkIn") + std::to_string(link));
  }
}

void GCTSumFileInputProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("inputFile", "");
  desc.add<unsigned int>("posOffset", 0);
  desc.add<unsigned int>("negOffset", 12);
  desc.add<unsigned int>("eventOffset", 0);
  desc.add<bool>("wrapAround", false);
  desc.add<bool>("debug", false);
  descriptions.add("gctSumFileInputProducer", desc);
}

DEFINE_FWK_MODULE(GCTSumFileInputProducer);
