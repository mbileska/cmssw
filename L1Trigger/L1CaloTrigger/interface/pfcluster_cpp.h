#include "algo_topIP1gct_h.h"
#include "bitonicSort16_h.h"

    
     loop  bestOf2e(const ap_uint<12>& e0, const ap_uint<12>& e1) {
        loop x = (e0 > e1) ? 0:1 ;
        return x;
    }

    rcttower  bestOf2(const rcttower& ecaltp0, const rcttower& ecaltp1) {
        rcttower x;
        x = (ecaltp0.energy > ecaltp1.energy)?ecaltp0:ecaltp1;
        return x;
    }

void findmaxoftenAe(rcttower towers[TEN], rcttower &MaxTEN){


	rcttower Step1[5] ;

	rcttower Step2[2] ;


	   for(loop k=0; k<TEN; k=k+2){
        Step1[k/2]  = towers[k+bestOf2e(towers[k].energy, towers[k+1].energy)] ;
        }

	   for(loop k=0; k<4; k=k+2){
        Step2[k/2]  = Step1[k+bestOf2e(Step1[k].energy, Step1[k+1].energy)] ;
        }

        rcttower x1  = Step2[bestOf2e(Step2[0].energy, Step2[1].energy)] ;

	MaxTEN  = bestOf2(x1,Step1[4]) ;
}

void findmaxoftenA(rcttower towers[TEN], rcttower &MaxTEN){


	rcttower Step1[5] ;

	rcttower Step2[2] ;


	   for(loop k=0; k<TEN; k=k+2){

        Step1[k/2]  = bestOf2(towers[k], towers[k+1]) ;
        }

	   for(loop k=0; k<4; k=k+2){

        Step2[k/2]  = bestOf2(Step1[k], Step1[k+1]) ;
        }


        rcttower x1 = bestOf2(Step2[0], Step2[1]) ;

	MaxTEN  = bestOf2(x1,Step1[4]) ;

}

void findmaxoftenB(rcttower towers[TEN], rcttower &MaxTEN){


	rcttower Step1[5] ;

	rcttower Step2[2] ;


	   for(loop k=0; k<TEN; k=k+2){
        Step1[k/2]  = bestOf2(towers[k], towers[k+1]) ;
        }

	   for(loop k=0; k<4; k=k+2){
        Step2[k/2]  = bestOf2(Step1[k], Step1[k+1]) ;
        }


        rcttower x1 = bestOf2(Step2[0], Step2[1]) ;

	MaxTEN  = bestOf2(x1,Step1[4]) ;

}

void getseedposition(rcttower towers[TEN][TEN], rcttower &Seed){



        rcttower EtaSlicesOUT[TEN] ;


        int delay_var=0;
        for(loop i=0; i<1; i++){
        rcttower towers1D[TEN] ;

                        for(loop j=0; j<TEN; j++){
						#pragma HLS PIPELINE

                        rcttower tmp = towers[i][j] ;
                        towers1D[j] = tmp ;
                        }
         delay_var += 1;
        findmaxoftenAe(towers1D, EtaSlicesOUT[i]) ;
        }

        for(loop i=delay_var; i<TEN; i++){
                rcttower towers1D[TEN] ;
  
                                for(loop j=0; j<TEN; j++){
								//#pragma HLS UNROLL
                                rcttower tmp = towers[i][j] ;
                                towers1D[j] = tmp ;
                                }
                findmaxoftenAe(towers1D, EtaSlicesOUT[i]) ;
                }

        findmaxoftenB(EtaSlicesOUT, Seed) ;

    }



ap_uint<12> computeMaskedSum(rcttower towers_seq[TEN][TEN], towermask mask[TEN+2][TEN+2]) {


	ap_uint<12> tmpSums[10];



  for (int i = 0; i < 10; i++) {
          tmpSums[i] = 0; // Reset tmpSums
    }



  ap_uint<12> tmpValue = 0;
  ap_uint<12> energytmp =0;
  ap_uint<12> energy =0;

  for (int i = 0; i < TEN; i++) {
    for (int j = 0; j < TEN; j++) {
      energytmp = towers_seq[i][j].energy;
      energy = energytmp * mask[i + 1][j + 1].energy;
      tmpSums[i] +=  energy;
    }

  }


  for (int i = 0; i < 10; i++) {
      tmpValue += tmpSums[i];
    }

  return tmpValue;
}




void getcluster(rcttower towers[TEN][TEN], rcttower Seed, PFcluster &output) {



	rcttower towers_seq[TEN][TEN];

	for (int blk = 0; blk < 5; blk++){
	  for (int i = 0; i < 2; i++) {
	    int row = blk * 2 + i;
	    if (row < TEN) {
	      for (int j = 0; j < TEN; j++) {
	        towers_seq[row][j].energy = towers[row][j].energy;
	      }
	    }
	  }
	}

	towermask mask[TEN+2][TEN+2] ;


	ap_uint<8> eta = Seed.eta ;
	ap_uint<7> phi = Seed.phi ;

	for (int i = 0; i < 3; i++) {
		  if (i >= eta && i <= eta + 2) {
		    for (int j = 0; j < TEN + 2; j++) {
		      if (j >= phi && j <= phi + 2)
		        mask[i][j].energy = 1;
		    }
		  }
		}

		// Remaining rows, slow access in blocks of 3
		for (int blk = 1; blk < 4; blk++) {
		  for (int i = 0; i < 3; i++) {
		    int row = blk * 3 + i ;
		    if (row >= TEN + 2) continue;
		    if (row >= eta && row <= eta + 2) {
		      for (int j = 0; j < TEN + 2; j++) {
		        if (j >= phi && j <= phi + 2)
		          mask[row][j].energy = 1;
		      }
		    }
		  }
		}

	//loadTowerSeq(towers, towers_seq);
	//generateMask(mask, Seed.eta, Seed.phi);
	ap_uint<12> tmpValue = computeMaskedSum(towers_seq, mask);


      /* ap_uint<12> tmpValue = 0  ;

          for(loop i=0; i<TEN; i++){
          	for(loop j=0; j<TEN; j++){
       	ap_uint<12> energytmp = towers_seq[i][j].energy ;
       	ap_uint<12> energy = energytmp *  mask[i+1][j+1].energy ;
       	ap_uint<12> tmp = tmpValue + energy ;
       	tmpValue = tmp  ;
		}
     }*/

	if(Seed.eta < 2 || Seed.eta > 7 || Seed.phi < 2 || Seed.phi > 7) tmpValue  = 0 ;
	output.energy = tmpValue ;
	output.eta = Seed.eta ;
	output.phi = Seed.phi ;
	output.getHoE() ; 
}



void zerrowtowers(rcttower towers[TEN][TEN], rcttower Seed) {





	towermask mask[TEN+2][TEN+2] ;
	towermask maskI[TEN+2][TEN+2] ;


	ap_uint<8> eta = Seed.eta ;
	ap_uint<7> phi = Seed.phi ;


        for(loop i=0; i<TEN+2; i++){

	if(i >= eta && i <= eta+2){
        for(loop j=0; j<TEN+2; j++){

	if(j >= phi && j <= phi+2 ) mask[i][j].energy =1 ;
	}}}

        	for(loop i=0; i<TEN+2; i++){

        	for(loop j=0; j<TEN+2; j++){

        	maskI[i][j].energy = ((ap_uint<1>)1 - mask[i][j].energy) ;
         	}}

         for(loop i=0; i<TEN; i++){

           for(loop j=0; j<TEN; j++){

		ap_uint<12> tmp = towers[i][j].energy ;
		ap_uint<12> energy = tmp * maskI[i+1][j+1].energy  ;
        	towers[i][j].energy = energy ;
		}}

}

void getPFClusters1(rcttower  subregiontowers[TEN][TEN], PFcluster  pfcluster[FIVE]){



		rcttower  Seed ;

		int delay_var = 0;
			for(loop k=0; k<1; k++){

				getseedposition(subregiontowers, Seed) ;
			getcluster(subregiontowers,Seed,pfcluster[k]) ;
			zerrowtowers(subregiontowers, Seed) ;
			delay_var += 1;
			}

			for(loop k= delay_var; k<FIVE; k++){

					getseedposition(subregiontowers, Seed) ;
				getcluster(subregiontowers,Seed,pfcluster[k]) ;
				zerrowtowers(subregiontowers, Seed) ;
				}

}


void createPFClusters(rcttower regiongcttowers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4], PFcluster pfclusters[N_GCT_PFCLUSTERS_REGION]) {





rcttower regiongcttowers_seq[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4];
rcttower  subregiontowers1[TEN][TEN] ;
rcttower  subregiontowers2[TEN][TEN] ;
rcttower  subregiontowers3[TEN][TEN] ;
PFcluster  pfclusters1[FIVE] ;
PFcluster  pfclusters2[FIVE] ;
PFcluster  pfclusters3[FIVE] ;
PFcluster  localIN[Nbclusters] ;
PFcluster  localOUT[Nbclusters] ;
rcttower subregiontowers[TEN][TEN];
PFcluster pfclustersTemp[Nbclusters];
PFcluster pfclusters_all[Nbclusters];






	for (int r = 0; r < 3; r++) {
		    for (int i = 0; i < TEN; i++) {
		        for (int k = 0; k < TEN; k++) {
	            subregiontowers[i][k] = regiongcttowers[i + r*6][k];
		        }
		    }

		    // Call once per region
		    getPFClusters1(subregiontowers, pfclustersTemp);


		        // Write to the region-specific slot
		        for (int i = 0; i < FIVE; i++) {
		            pfclusters_all[r * FIVE + i] = pfclustersTemp[i];
		        }


		}




	// Copy result to correct pfclustersX array
	for (int i = 0; i < FIVE; i++) {

			pfclusters1[i] = pfclusters_all[i];

			pfclusters2[i] = pfclusters_all[i+5];

			pfclusters3[i] = pfclusters_all[i+10];
		}



	for(loop i=0; i<FIVE; i++){

          localIN[i+0] = pfclusters1[i] ;
          localIN[i+5] = pfclusters2[i] ;
          localIN[i+10] = pfclusters3[i] ;
	}


	 bitonicSort16AB(localIN, localOUT);


	for (int i = 0; i < 12; i++) {
	  pfclusters[i] = localOUT[15 - i];
	}


}
