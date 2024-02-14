#ifndef L1Trigger_L1TCaloLayer1_UCTLogging_hh
#define L1Trigger_L1TCaloLayer1_UCTLogging_hh

#define CMSSW

#ifdef CMSSW
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#define LOG_ERROR edm::LogError("L1TCaloLayer1")
#else
#define LOG_ERROR std::cerr
#endif

#endif
