/*
 * Description:
 *   Deterministic source-faithfulness test for the attached GCT IP2 firmware.
 *
 *   The expected values intentionally test what the attachment actually does,
 *   not what the interface PDF ultimately wants the subsystem to do:
 *     - fixed EG source indices 0,1,8,9,15,22;
 *     - only two populated jet slots;
 *     - only two populated tau slots;
 *     - padded seed eta/phi serialized unchanged;
 *     - Ex=Ey=0 and HT equal to the interior ST-energy sum;
 *     - output link 5 remains zero.
 */

#include <ap_int.h>

#include <array>
#include <iostream>
#include <string>

#include "L1Trigger/L1CaloTrigger/interface/algo_topIP2gct_h.h"
#include "L1Trigger/L1CaloTrigger/interface/gctforgtIP2_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/algo_topIP2gct_cpp.h"

namespace {

bool expectEqual(const std::string& label, unsigned int observed, unsigned int expected) {
  if (observed == expected) {
    return true;
  }
  std::cerr << "[FAIL] " << label << ": observed " << observed << ", expected " << expected << '\n';
  return false;
}

void putIP1EG(gctip2::LinkWord& link, int slot, unsigned int energy, unsigned int eta, unsigned int phi) {
  ap_uint<64> object = 0;
  object.range(11, 0) = energy;
  object.range(18, 12) = eta;
  object.range(25, 19) = phi;
  link.range(64 * slot + 63, 64 * slot) = object;
}

void putIP2OrderedST(gctip2::LinkWord& link, int slot, unsigned int energy) {
  gctip2::STower tower;
  tower.energy = energy;
  link.range(48 * slot + 47, 48 * slot) = tower.pack();
}

}  // namespace

int main() {
  std::array<gctip2::LinkWord, gctip2::kInputLinks> input{};
  std::array<gctip2::LinkWord, gctip2::kOutputLinks> output{};

  // Positive EG inputs are read by the top as 0,1,4,5.
  for (int slot = 0; slot < 8; ++slot) {
    putIP1EG(input[0], slot, 100 + slot, 10 + slot, 20 + slot);
    putIP1EG(input[1], slot, 200 + slot, 30 + slot, 40 + slot);
    putIP1EG(input[4], slot, 300 + slot, 50 + slot, 60 + slot);
    putIP1EG(input[5], slot, 400 + slot, 70 + slot, 80 + slot);
  }

  // input[9] becomes padded phi columns 2 and 3.  The attached ST contract is
  // slots0..5 = eta0..5 in the first phi column and slots6..11 = eta0..5 in
  // the second phi column.
  putIP2OrderedST(input[9], 0, 10);  // padded (eta,phi)=(2,2)
  putIP2OrderedST(input[9], 6, 5);   // padded (eta,phi)=(2,3)

  gctip2::algoTop(input, output);

  bool ok = true;

  // Fixed source indices from createEGs(): 0,1,8,9,15,22.
  constexpr std::array<unsigned int, 6> expectedEG{{100, 101, 200, 201, 207, 306}};
  for (int slot = 0; slot < 6; ++slot) {
    gctip2::GctEcalCluster eg;
    eg.unpack(output[0].range(48 * slot + 47, 48 * slot));
    ok &= expectEqual("positive EG energy slot " + std::to_string(slot),
                      static_cast<unsigned int>(eg.energy),
                      expectedEG[slot]);
  }

  gctip2::Jet jet0;
  jet0.unpack(output[1].range(47, 0));
  ok &= expectEqual("jet0 energy", static_cast<unsigned int>(jet0.energy), 15);
  ok &= expectEqual("jet0 padded eta", static_cast<unsigned int>(static_cast<ap_uint<6>>(jet0.eta)), 2);
  ok &= expectEqual("jet0 padded phi", static_cast<unsigned int>(static_cast<ap_uint<7>>(jet0.phi)), 2);

  for (int slot = 1; slot < 6; ++slot) {
    gctip2::Jet jet;
    jet.unpack(output[1].range(48 * slot + 47, 48 * slot));
    ok &= expectEqual("jet energy slot " + std::to_string(slot), static_cast<unsigned int>(jet.energy), 0);
  }

  gctip2::Tau tau0;
  gctip2::Tau tau1;
  tau0.unpack(output[1].range(48 * 6 + 47, 48 * 6));
  tau1.unpack(output[1].range(48 * 7 + 47, 48 * 7));
  ok &= expectEqual("tau0 energy", static_cast<unsigned int>(tau0.energy), 10);
  ok &= expectEqual("tau0 padded eta", static_cast<unsigned int>(static_cast<ap_uint<6>>(tau0.eta)), 2);
  ok &= expectEqual("tau0 padded phi", static_cast<unsigned int>(static_cast<ap_uint<7>>(tau0.phi)), 2);
  ok &= expectEqual("tau1 energy", static_cast<unsigned int>(tau1.energy), 5);
  ok &= expectEqual("tau1 padded eta", static_cast<unsigned int>(static_cast<ap_uint<6>>(tau1.eta)), 2);
  ok &= expectEqual("tau1 padded phi", static_cast<unsigned int>(static_cast<ap_uint<7>>(tau1.phi)), 3);

  for (int slot = 8; slot < 12; ++slot) {
    gctip2::Tau tau;
    tau.unpack(output[1].range(48 * slot + 47, 48 * slot));
    ok &= expectEqual("tau energy slot " + std::to_string(slot - 6), static_cast<unsigned int>(tau.energy), 0);
  }

  gctip2::Sums sums;
  sums.unpack(output[2].range(47, 0));
  ok &= expectEqual("sum Ex", static_cast<unsigned int>(sums.Ex), 0);
  ok &= expectEqual("sum Ey", static_cast<unsigned int>(sums.Ey), 0);
  ok &= expectEqual("sum HT", static_cast<unsigned int>(sums.HT), 15);
  ok &= expectEqual("spare output link", static_cast<unsigned int>(output[5] != 0), 0);

  if (!ok) {
    return 1;
  }
  std::cout << "[PASS] attached-IP2 source-faithful deterministic vector\n";
  return 0;
}
