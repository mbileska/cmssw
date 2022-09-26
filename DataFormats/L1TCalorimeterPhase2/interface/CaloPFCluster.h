#ifndef DataFormats_L1TCalorimeterPhase2_CaloPFCluster_h
#define DataFormats_L1TCalorimeterPhase2_CaloPFCluster_h

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace l1tp2 {

  class CaloPFCluster : public l1t::L1Candidate {
  public:
    CaloPFCluster() : l1t::L1Candidate(), clusterEt_(0.), clusterIEta_(-99), clusterIPhi_(-99){};

    CaloPFCluster(const PolarLorentzVector& p4, float clusterEt, int clusterIEta, int clusterIPhi)
        : l1t::L1Candidate(p4), clusterEt_(clusterEt), clusterIEta_(clusterIEta), clusterIPhi_(clusterIPhi){};

    inline float clusterEt() const { return clusterEt_; };
    inline int clusterIEta() const { return clusterIEta_; };
    inline int clusterIPhi() const { return clusterIPhi_; };
    void setClusterEt(float clusterEtIn) { clusterEt_ = clusterEtIn; };
    void setClusterIEta(float clusterIEtaIn) { clusterIEta_ = clusterIEtaIn; };
    void setClusterIPhi(float clusterIPhiIn) { clusterIPhi_ = clusterIPhiIn; };

  private:
    // ET
    float clusterEt_ = 0.;
    // GCT ieta
    int clusterIEta_ = -99;
    // GCT iphi
    int clusterIPhi_ = -99;
  };

  // Concrete collection of output objects (with extra tuning information)
  typedef std::vector<l1tp2::CaloPFCluster> CaloPFClusterCollection;
}  // namespace l1tp2
#endif
