#include "algo_topIP2_h.h"
#include <iostream>
#include <iomanip>

using namespace l1ct;

int main() {

    ap_uint<576> link_in[N_HF_REGIONS];
    ap_uint<576> link_out[N_HF_REGIONS];

//    ap_uint<64> input_hex[N_HF_REGIONS][NCALO] = {
//
//        {0x040307B, 0x00000000004050BD, 0x000000000030217A, 0x00000000004071C5, 0x00000000007051EB, 0x0000000000305266},
//
//        {0x0902000, 0x00000000005070FC, 0x0000000000404135, 0x0000000000702137, 0x0000000000305138, 0x00000000003031EC},
//
//		{0x0902000, 0x0000000000702137, 0x000000000030217A, 0x0000000000404192, 0x00000000005071E7, 0x0000000000305266},
//
//		{0x0402008, 0x000000000030707C, 0x00000000004060D8, 0x000000000040311C, 0x000000000070317B, 0x000000000030317B},
//
//		{0x04020AA, 0x00000000005070CE, 0x0000000000302125, 0x0000000000705134, 0x000000000050415E, 0x0000000000307182},
//
//		{0x0902000, 0x00000000007020F7, 0x0000000000302182, 0x0000000000404192, 0x00000000003051B4, 0x00000000005071E7}
//    };



//     Pack PFClusters into links
//    for(int r=0;r<N_HF_REGIONS;r++){
//
//        link_in[r] = 0;
//
//        for(int c=0;c<NCALO;c++){
//            link_in[r].range((c+1)*64-1,c*64) = input_hex[r][c];
//        }
//    }

    link_in[0] = ap_uint<576>("00000000000000000000000000000000040307B00000000004050BD000000000030217A00000000004071C500000000007051EB0000000000305266", 16);
	link_in[1] = ap_uint<576>("0000000000000000000000000000000009020000000000005070FC00000000004041350000000000702137000000000030513800000000003031EC", 16);
	link_in[2] = ap_uint<576>("000000000000000000000000000000000902000000000000702137000000000030217A000000000040419200000000005071E70000000000305266", 16);
	link_in[3] = ap_uint<576>("00000000000000000000000000000000040200800000000030707C00000000004060D8000000000040311C000000000070317B000000000030317B", 16);
	link_in[4] = ap_uint<576>("0000000000000000000000000000000004020AA0000000005070CE00000000003021250000000000705134000000000050415E0000000000307182", 16);
	link_in[5] = ap_uint<576>("0000000000000000000000000000000009020000000000007020F70000000000302182000000000040419200000000003051B400000000005071E7", 16);

    // Print input PFClusters
//     std::cout << "\n===== INPUT PFCLUSTERS =====\n";
//
//    for(int r=0;r<N_HF_REGIONS;r++){
//
//        std::cout<<"Region "<<r<<"\n";
//
//        for(int c=0;c<NCALO;c++){
//
//            std::cout<<"  PF["<<c<<"] = "
//                     <<std::hex<<input_hex[r][c]
//                     <<std::dec<<"\n";
//        }
//    }


    // Print raw input links
//    std::cout<<"\n===== RAW INPUT LINKS =====\n";
//
//    for(int r=0;r<N_HF_REGIONS;r++){
//
//        std::cout<<"link_in["<<r<<"] = "
//                 <<std::hex<<link_in[r]
//                 <<std::dec<<"\n";
//    }

    // Print raw input links
    std::cout << "\n===== RAW INPUT LINKS =====\n";

    for(int r=0; r<N_HF_REGIONS; r++) {

        std::cout << "link_in[" << r << "] = 0x"
                  << link_in[r].to_string(16)
                  << "\n";
    }


    // Run IP
    algo_topIP2(link_in,link_out);

    // Print hex output links
    std::cout<<"\n===== RAW OUTPUT LINKS =====\n";

//    for(int r=0;r<N_HF_REGIONS;r++){
//
//        std::cout<<"link_out["<<r<<"] = "
//                 <<std::hex<<link_out[r]
//                 <<std::dec<<"\n";
//    }
    for(int r=0; r<N_HF_REGIONS; r++) {

        std::cout << "link_out[" << r << "] = 0x"
                  << link_out[r].to_string(16)
                  << "\n";
    }


    // Decode Puppi Objects
    std::cout<<"\n===== DECODED PUPPI OBJECTS =====\n";

    for(int r=0;r<N_HF_REGIONS;r++){

        std::cout<<"Region "<<r<<"\n";

        for(int i=0;i<NNEUTRALS;i++){

            PuppiObj obj;
            obj.initFromBits(link_out[r].range((i+1)*64-1,i*64));

            std::cout<<"  puppi["<<i<<"] "
                     <<"Pt="<<obj.hwPt
                     <<" Eta="<<obj.hwEta
                     <<" Phi="<<obj.hwPhi
                     <<" Id="<<obj.hwId.rawId()
                     <<"\n";
        }
    }

    return 0;
}
