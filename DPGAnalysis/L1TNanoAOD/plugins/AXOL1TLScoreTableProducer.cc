#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/L1TGlobal/interface/AXOL1TLScore.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

class AXOL1TLScoreTableProducer : public edm::global::EDProducer<> {
public:
  explicit AXOL1TLScoreTableProducer(const edm::ParameterSet& iConfig)
      : axoToken_(consumes<AXOL1TLScoreBxCollection>(
            iConfig.getParameter<edm::InputTag>("src"))),
        name_(iConfig.getParameter<std::string>("name")) {
    produces<nanoaod::FlatTable>();
  }

  void produce(edm::StreamID, edm::Event& iEvent,
	       const edm::EventSetup&) const override {

    const auto& axoScores = iEvent.get(axoToken_);

      // DEBUG
    edm::LogWarning("AXOL1TLTableProducer") 
      << "firstBX=" << axoScores.getFirstBX() 
      << " lastBX=" << axoScores.getLastBX()
      << " size=" << axoScores.size();
    for (int bx = axoScores.getFirstBX(); bx <= axoScores.getLastBX(); ++bx) {
      edm::LogWarning("AXOL1TLTableProducer")
        << "BX=" << bx << " isEmpty=" << axoScores.isEmpty(bx);
      if (!axoScores.isEmpty(bx)) {
	edm::LogWarning("AXOL1TLTableProducer")
          << "  score=" << axoScores.at(bx, 0).getAXOScore();
      }
    }
    // END DEBUG
    
    float score = -999.f;
    if (axoScores.getFirstBX() <= 0 && axoScores.getLastBX() >= 0 &&
	!axoScores.isEmpty(0)) {
      score = axoScores.at(0, 0).getAXOScore();
    }

    auto table = std::make_unique<nanoaod::FlatTable>(1, name_, true);
    table->addColumnValue<float>(
				 "score", score, "AXOL1TL anomaly detection score at BX=0");

    iEvent.put(std::move(table));
  }
    
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src",
        edm::InputTag("simGtStage2Digis", "AXOScore"));
    desc.add<std::string>("name", "L1AXO");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<AXOL1TLScoreBxCollection> axoToken_;
  const std::string name_;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(AXOL1TLScoreTableProducer);
