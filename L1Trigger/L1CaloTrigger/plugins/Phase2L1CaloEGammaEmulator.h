//------------------------------------
// Helper functions for Phase2L1CaloEGammaEmulator.h
//------------------------------------

#include <ap_int.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>

#ifndef EGAMMA_CRYSTALS_PRODUCER_H
#define EGAMMA_CRYSTALS_PRODUCER_H

static constexpr int dimPhi = 4;
static constexpr int dimEta = 17;
static constexpr int dimCard = 36;
static constexpr int nClustersPerLink = 2;



//--------------------------------------------------------// 

// Helper function to print a 4x17x36 int array to a ostream f.                                                             // oneCard: optional argument. If set to some value, only prints
// one card.   

inline void printL1ArrayInt(ofstream& f, int array[dimPhi][dimEta][dimCard],
			    string desc = "",
			    bool printOnlyOneCard = false,
			    int oneCard = -1) {

  
  for (int kk = 0; kk < dimCard; kk++) {
    
    if ((!printOnlyOneCard) || (printOnlyOneCard && (kk == oneCard))) {
      
      f << "[CARD " << kk << "]: " << desc << std::endl;
      for (int ii = 0; ii < dimPhi; ii++) {
	for (int jj = 0; jj < dimEta; jj++) {
	  f << array[ii][jj][kk] << "\t";
	}
	f << std::endl;
      }
    }
  } 

}

//--------------------------------------------------------//                                                                  

// Helper function to print a 4x17x36 float array to a ostream f.                                                                    
          
inline void printL1ArrayFloat(ofstream& f, float array[dimPhi][dimEta][dimCard],
			      string desc = "", bool printOnlyOneCard = false, 
			      int oneCard = -1,
			      unsigned int precision = 3) {

  for (int kk = 0; kk < dimCard; kk++) {

    if ((!printOnlyOneCard) || (printOnlyOneCard && (kk == oneCard))) {
      f << "[CARD " << kk << "]: " << desc << std::endl;
      for (int ii = 0; ii < dimPhi; ii++) {
	for (int jj = 0; jj < dimEta; jj++) {
	  f << std::setprecision(precision) << array[ii][jj][kk] << "\t";
	}
	f << std::endl;
      }
    }
  }
}

//--------------------------------------------------------//

// Helper function to print a 4x17x36 ap_uint<12> array to a ostream f,
// where the ap_uint<12> are encoded floats (divide by 8).

// float et = hit.encodedEt() / 8.;

inline void printL1ArrayEncodedEt(ofstream& f, ap_uint<12> array[dimPhi][dimEta][dimCard],
				  string desc = "", bool printOnlyOneCard = false, 
				  int oneCard = -1, 
				  unsigned int precision = 3 ) {

  for (int kk = 0; kk < dimCard; kk++) {

    if ((!printOnlyOneCard) || (printOnlyOneCard && (kk == oneCard))) {
      f << "[CARD " << kk << "]: " << desc << std::endl;
      for (int ii = 0; ii < dimPhi; ii++) {
	for (int jj = 0; jj < dimEta; jj++) {
	  float fVal = array[ii][jj][kk] / 8.;
	  ap_uint<12> uVal = array[ii][jj][kk];
	  
	  f << std::setprecision(precision) << fVal << "\t";
	}
	// Write newline after each row
	f << std::endl;
      }
    }
  }

}

//--------------------------------------------------------//                                                  

// Helper function to print a 4x17x36 ap_uint<12> array to a ostream f,                                                    
// where the ap_uint<12> are encoded floats (divide by 8).

// Use LSB 0.5
// https://github.com/cms-sw/cmssw/blob/master/L1Trigger/L1TCalorimeter/python/caloParams_cfi.py#L15

inline void printL1ArrayCompressedEt(ofstream& f, ap_uint<12> array[dimPhi][dimEta][dimCard],
				     string desc = "", bool printOnlyOneCard = false, 
				     int oneCard = -1,
				     unsigned int precision = 3 ) {
  
  
  float LSB = 0.5;

  for (int kk = 0; kk < dimCard; kk++) {

    if ((!printOnlyOneCard) || (printOnlyOneCard && (kk == oneCard))) {
      f << "[CARD " << kk << "]: " << desc << std::endl;
      for (int ii = 0; ii < dimPhi; ii++) {
	for (int jj = 0; jj < dimEta; jj++) {
	  float fVal = array[ii][jj][kk] * LSB;
	  ap_uint<12> uVal = array[ii][jj][kk];
	  
        f << std::setprecision(precision) << fVal << "\t";
	}
	// Write newline after each row                                                                                                    
	f << std::endl;
      }
    }
  }

}

//--------------------------------------------------------//                                                                 

// Helper function to print a 4x2x36 int array to a ostream f.                                                              

inline void printL1Array_ClusterInt(ofstream& f, int array[dimPhi][nClustersPerLink][dimCard],
				    string desc = "", bool printOnlyOneCard = false, int oneCard = -1) {

  for (int kk = 0; kk < dimCard; kk++) {
    if ((!printOnlyOneCard) || (printOnlyOneCard && (kk == oneCard))) {

      f << "[CARD " << kk << "]: " << desc << std::endl;
      for (int ii = 0; ii < dimPhi; ii++) {
	for (int jj = 0; jj < nClustersPerLink; jj++) {
	  f << array[ii][jj][kk] << "\t";
	}
	f << std::endl;
      }
    }
  }
}

//--------------------------------------------------------//

// Helper function to print a 4x2x36 ap uint array to a ostream f.

inline void printL1Array_ClusterUint(ofstream& f, ap_uint<12> array[dimPhi][nClustersPerLink][dimCard],
				     string desc = "", bool printOnlyOneCard = false, int oneCard = -1) {

  for (int kk = 0; kk < dimCard; kk++) {
    if ((!printOnlyOneCard) || (printOnlyOneCard && (kk == oneCard))) {

      f << "[CARD " << kk << "]: " << desc << std::endl;
      for (int ii = 0; ii < dimPhi; ii++) {
        for (int jj = 0; jj < nClustersPerLink; jj++) {
	  // Convert to float
          f << (array[ii][jj][kk])/8. << "\t";
        }
        f << std::endl;
      }
    }
  } 
  
}

//--------------------------------------------------------//                                                                  

// Helper function to print a 4x2x36 float array to a ostream f.                                                            

inline void printL1Array_ClusterFloat(ofstream& f, float array[dimPhi][nClustersPerLink][dimCard],
                                     string desc = "", bool printOnlyOneCard = false, int oneCard = -1) {

  for (int kk = 0; kk < dimCard; kk++) {
    if ((!printOnlyOneCard) || (printOnlyOneCard && (kk == oneCard))) {

      f << "[CARD " << kk << "]: " << desc << std::endl;
      for (int ii = 0; ii < dimPhi; ii++) {
        for (int jj = 0; jj < nClustersPerLink; jj++) {
          f << array[ii][jj][kk] << "\t";
        }
        f << std::endl;
      }
    }
  }

}




#endif

