#include <ap_int.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

#include "L1Trigger/L1CaloTrigger/interface/MHH_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_MHH_h.h"
#include "L1Trigger/L1CaloTrigger/interface/bitonicSort32_MHH_cpp.h"
#include "L1Trigger/L1CaloTrigger/interface/MHH_cpp.h"

namespace {

void require(bool condition, const std::string& message) {
  if (!condition)
    throw std::runtime_error(message);
}

void setObject(ap_uint<576>& link, unsigned int slot, ap_uint<48> object) {
  link.range(slot * 48 + 47, slot * 48) = object;
}

ap_uint<48> objectAt(const ap_uint<576>& link, unsigned int slot) {
  return link.range(slot * 48 + 47, slot * 48);
}

ap_uint<48> gamma(unsigned int et, unsigned int eta, unsigned int phi, bool hgcal) {
  ap_uint<48> out = 0;
  out.range(11, 0) = et;
  out.range(18, 12) = eta;
  out.range(27, 19) = phi;
  out[47] = hgcal;
  return out;
}

ap_uint<48> hadron(unsigned int et, unsigned int eta, unsigned int phi, unsigned int ratio, bool hgcal) {
  ap_uint<48> out = 0;
  out.range(11, 0) = et;
  out.range(17, 12) = eta;
  out.range(26, 18) = phi;
  out.range(30, 27) = ratio;
  out[47] = hgcal;
  return out;
}

ap_uint<576> sourceSums(
    int ex, int ey, int htx, int hty, unsigned int sumEt, unsigned int ht, unsigned int count) {
  ap_uint<576> out = 0;
  out.range(11, 0) = ap_int<12>(ex);
  out.range(23, 12) = ap_int<12>(ey);
  out.range(35, 24) = ap_int<12>(htx);
  out.range(47, 36) = ap_int<12>(hty);
  out.range(75, 64) = ap_uint<12>(sumEt);
  out.range(87, 76) = ap_uint<12>(ht);
  out.range(99, 88) = ap_uint<12>(count);
  return out;
}

std::array<ap_uint<576>, p2mhh::MHH_N_OUTPUT_LINKS> run(
    std::array<ap_uint<576>, p2mhh::MHH_N_INPUT_LINKS> input) {
  std::array<ap_uint<576>, p2mhh::MHH_N_OUTPUT_LINKS> output{};
  p2mhh::mhh_algo_top(input.data(), output.data());
  return output;
}

void testSumsAbiV2AndReservedLinks() {
  std::array<ap_uint<576>, p2mhh::MHH_N_INPUT_LINKS> input{};

  // Exercise signed and unsigned 12-bit modular reduction exactly as in the
  // frozen mhh-sums-v2 firmware vectors.
  input[2] = sourceSums(2047, -2048, 2040, -2040, 4090, 4091, 4095);
  input[10] = sourceSums(1, -1, 100, -100, 20, 10, 2);

  // The remaining ten input links are ignored by the current MHH algorithm.
  for (unsigned int link : {3U, 4U, 5U, 6U, 7U, 11U, 12U, 13U, 14U, 15U})
    input[link] = ~ap_uint<576>(0);

  const auto output = run(input);
  const ap_uint<576>& sums = output[2];

  require((ap_int<12>)sums.range(11, 0) == -2048, "MHH Ex modular sum differs from firmware");
  require((ap_int<12>)sums.range(23, 12) == 2047, "MHH Ey modular sum differs from firmware");
  require((ap_int<12>)sums.range(35, 24) == -1956, "MHH HTx modular sum differs from firmware");
  require((ap_int<12>)sums.range(47, 36) == 1956, "MHH HTy modular sum differs from firmware");
  require((ap_uint<12>)sums.range(75, 64) == 14, "MHH total-ET modular sum differs from firmware");
  require((ap_uint<12>)sums.range(87, 76) == 5, "MHH total-HT modular sum differs from firmware");
  require((ap_uint<12>)sums.range(99, 88) == 1, "MHH object-count modular sum differs from firmware");
  require(sums.range(63, 48) == 0, "MHH sums word 0 reserved bits are nonzero");
  require(sums.range(127, 100) == 0, "MHH sums word 1 reserved bits are nonzero");
  require(sums.range(575, 128) == 0, "MHH sums words 2--8 are nonzero");
  require(output[0] == 0 && output[1] == 0, "empty MHH object inputs must produce zero object links");
  require(output[3] == 0 && output[4] == 0 && output[5] == 0, "MHH outputs 3--5 must remain reserved zero");
}

void testBoundaryStitchingAndPacking() {
  std::array<ap_uint<576>, p2mhh::MHH_N_INPUT_LINKS> input{};

  // Gamma match across the 359/0 wrap: HF wins and retains its coordinates.
  setObject(input[0], 0, gamma(100, 85, 359, false));
  setObject(input[8], 0, gamma(20, 84, 0, true));

  // Equal ET selects HGCAL, matching the firmware's explicit tie behavior.
  setObject(input[0], 1, gamma(50, 85, 100, false));
  setObject(input[8], 1, gamma(50, 84, 100, true));

  // Isolated-EG stitching saturates at the 12-bit maximum.
  setObject(input[0], 6, gamma(3000, 85, 40, false));
  setObject(input[8], 6, gamma(2000, 84, 40, true));

  // Hadronic stitching keeps eta/phi/ratio from the higher-ET HGCAL object.
  setObject(input[1], 0, hadron(30, 6, 23, 2, false));
  setObject(input[9], 0, hadron(70, 5, 0, 11, true));

  const auto output = run(input);

  const ap_uint<48> leadingEg = objectAt(output[0], 0);
  require(leadingEg.range(11, 0) == 120 && leadingEg.range(18, 12) == 85 &&
              leadingEg.range(27, 19) == 359 && leadingEg[47] == 0,
          "wrapped MHH gamma stitching differs from firmware");

  const ap_uint<48> tiedEg = objectAt(output[0], 1);
  require(tiedEg.range(11, 0) == 100 && tiedEg.range(18, 12) == 84 &&
              tiedEg.range(27, 19) == 100 && tiedEg[47] == 1,
          "MHH equal-ET gamma ownership differs from firmware");

  const ap_uint<48> leadingIso = objectAt(output[0], 6);
  require(leadingIso.range(11, 0) == 0xFFF, "MHH isolated-EG stitched ET did not saturate");

  const ap_uint<48> leadingJet = objectAt(output[1], 0);
  require(leadingJet.range(11, 0) == 100 && leadingJet.range(17, 12) == 5 &&
              leadingJet.range(26, 18) == 0 && leadingJet.range(30, 27) == 11 && leadingJet[47] == 1,
          "MHH hadronic stitching or 48-bit packing differs from firmware");
}

}  // namespace

#ifdef GCTSUM_SYSTEMC_AP_INT
int sc_main(int, char**) {
#else
int main() {
#endif
  try {
    testSumsAbiV2AndReservedLinks();
    testBoundaryStitchingAndPacking();
  } catch (const std::exception& error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }

  std::cout << "PASS: MHH clean emulator matches frozen sums ABI v2 and object packing\n";
  return 0;
}
