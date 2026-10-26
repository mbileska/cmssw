

#include "algo_topIP1gct_h.h"

void fillInputClusters(ap_uint<576> input, rctecalcluster rctecalclusters[N_RCT_CLUSTERS_REGION], towermask Grid[N_RCT_CLUSTERS_REGION]){


     for(loop i=0; i<N_RCT_CLUSTERS_REGION; i++) {

    	ap_uint<10> start   = 64*i;
    	ap_uint<10> end = start + 63;
	rctecalclusters[i].fillrctecalcluster(input.range(end, start)) ;
	ap_uint<5> posEta = rctecalclusters[i].eta/5+2 ;
	ap_uint<4> posPhi = rctecalclusters[i].phi/5+2 ;
	for(loop k=2; k<N_GCT_TOWERS_REGION_ETA+2; k++) {

		for(loop m=2; m<N_GCT_TOWERS_REGION_PHI+2; m++) {

		if (posEta == k && posPhi == m &&  rctecalclusters[i].energy > 0) { Grid[i].energy = 1 ; Grid[i].eta = k ; Grid[i].phi = m ; }
		else { }
		}
	}}} 


void fillInputTowers1(ap_uint<576> input, rcttower rcttowers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4]){


     for(loop i=0; i<N_GCT_TOWERS_REGION_ETA; i++) 
	{ 

    	ap_uint<10> start   = 16*i;
    	ap_uint<10> end = start + 15;
	rcttowers[i+2][2].fillrcttower(input.range(end, start)) ;
	}
     for(loop i=0; i<N_GCT_TOWERS_REGION_ETA; i++) 
	{ 

    	ap_uint<10> start1   = 16*(17+i);
    	ap_uint<10> end1 = start1 + 15;
	rcttowers[i+2][3].fillrcttower(input.range(end1, start1)) ;
	}
	} 


void fillInputTowers2(ap_uint<576> input, rcttower rcttowers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4]){


     for(loop i=0; i<N_GCT_TOWERS_REGION_ETA; i++) 
	{ 

    	ap_uint<10> start   = 16*i;
    	ap_uint<10> end = start + 15;
	rcttowers[i+2][4].fillrcttower(input.range(end, start)) ;
	}
     for(loop i=0; i<N_GCT_TOWERS_REGION_ETA; i++) 
	{ 

    	ap_uint<10> start1   = 16*(17+i);
    	ap_uint<10> end1 = start1 + 15;
	rcttowers[i+2][5].fillrcttower(input.range(end1, start1)) ;
	}
	} 


void fillInputTowers3(ap_uint<576> input, rcttower rcttowers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4]){


     for(loop i=0; i<N_GCT_TOWERS_REGION_ETA; i++) 
	{ 

    	ap_uint<10> start   = 16*i;
    	ap_uint<10> end = start + 15;
	rcttowers[i+2][6].fillrcttower(input.range(end, start)) ;
	}
     for(loop i=0; i<N_GCT_TOWERS_REGION_ETA; i++) 
	{ 

    	ap_uint<10> start1   = 16*(17+i);
    	ap_uint<10> end1 = start1 + 15;
	rcttowers[i+2][7].fillrcttower(input.range(end1, start1)) ;
	}
	} 

void buildfulltowers(rcttower rcttowers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4], rctecalcluster rctecalclusters[N_RCT_CLUSTERS_REGION], towermask Grid[N_RCT_CLUSTERS_REGION]){

	  rcttower rcttowers_seq[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4];


	  // Load
	 for (int k = 2; k < N_GCT_TOWERS_REGION_ETA + 2; k++) {
	   for (int m = 2; m < N_GCT_TOWERS_REGION_PHI + 2; m++) {
	     rcttowers_seq[k][m] = rcttowers[k][m];
	   }
	 }


	  for (loop i = 0; i < N_RCT_CLUSTERS_REGION; i++) {
	    ap_uint<12> tmp1 = rctecalclusters[i].energy;
	    for (loop m = 2; m < N_GCT_TOWERS_REGION_PHI+2; m++) {
	      for (loop k = 2; k < N_GCT_TOWERS_REGION_ETA+2; k++) {
	        if (Grid[i].eta == k && Grid[i].phi == m) {
	          rcttowers_seq[k][m].energy += tmp1;
	        }
	      }
	    }
	  }

	  // Write
	  for (int k = 2; k < N_GCT_TOWERS_REGION_ETA + 2; k++) {
	    for (int m = 2; m < N_GCT_TOWERS_REGION_PHI + 2; m++) {
	      rcttowers[k][m] = rcttowers_seq[k][m];  // Write
	    }
	  }

	}








void getSums(rctecalcluster low, rctecalcluster high, ap_uint<12> &EnergySum, ap_uint<10> &EnergySum5x5){

ap_uint<12> EnergyH = high.energy ;
ap_uint<10> EnergyH5x5 = high.et5x5 ;

ap_uint<12> EnergyL = low.energy ;
ap_uint<10> EnergyL5x5 = low.et5x5 ;

EnergySum = EnergyH + EnergyL ;
EnergySum5x5 = EnergyH5x5 + EnergyL5x5 ;


}

void stitchClustersOnEta(rctecalcluster rctecalclusters1[N_RCT_CLUSTERS_REGION], rctecalcluster rctecalclusters2[N_RCT_CLUSTERS_REGION]){


/* there are 9 clusters in each region, we stitch the energy and 0 one of two that has
are at eta=24 Low with eta=0 High if the dPhi<=2 */

for(loop j=0; j<9; j++){  /* R1 */

for(loop i=0; i<9; i++){  /* R2 */


bool phi1 = (rctecalclusters1[j].phi == 29) ? true : false ;
bool phi2 = (rctecalclusters2[i].phi == 0) ? true : false ;


ap_uint<5> eta1 = rctecalclusters1[j].eta ;
ap_uint<5> eta2 = rctecalclusters2[i].eta ;

ap_uint<5> dEta = (eta1 > eta2) ? (eta1 - eta2) : (eta2 - eta1) ;
bool deta = (dEta <= 1) ? true : false ;

bool stitch = phi1 && phi2 && deta ;

// cout << " Stitch 1 " << rctecalclusters1[j].phi << " 2 " << rctecalclusters2[i].phi << " ok? " << stitch << endl ;

ap_uint<12> Energy1 = rctecalclusters1[j].energy ;
ap_uint<12> Energy2 = rctecalclusters2[i].energy ;

ap_uint<12> EnergySum ;
ap_uint<10> EnergySum5x5 ;

getSums(rctecalclusters1[j], rctecalclusters2[i], EnergySum, EnergySum5x5) ;

if(stitch)
{
	if(Energy1 > Energy2){
	rctecalclusters1[j].energy = EnergySum ;
	rctecalclusters1[j].et5x5 = EnergySum5x5 ;
	rctecalclusters2[i].energy = 0 ;
	}
	else {
	rctecalclusters1[j].energy = 0 ;
	rctecalclusters2[i].energy = EnergySum ;
	rctecalclusters2[i].et5x5 = EnergySum5x5 ;
	}
}

}}

}

void stitchClustersOnPhi(rctecalcluster rctecalclusters1[N_RCT_CLUSTERS_REGION], rctecalcluster rctecalclusters2[N_RCT_CLUSTERS_REGION]){



for(loop j=0; j<9; j++){  /* R1 */
    
for(loop i=0; i<9; i++){  /* R2 */


bool eta1 = (rctecalclusters1[j].eta == 0) ? true : false ;
bool eta2 = (rctecalclusters2[i].eta == 0) ? true : false ;

ap_uint<5> phi1 = rctecalclusters1[j].phi ;
ap_uint<5> phi2 = rctecalclusters2[i].phi ;
ap_uint<5> phiSum = phi1 + phi2 ;
ap_uint<5> phiSum1 = phi1 + phi2 - 29;
ap_uint<5> phiSum2 = 29 - phi1 - phi2 ;

ap_uint<5> dPhi = (phiSum > 29) ? phiSum1 : phiSum2 ;

bool dphi = (dPhi <= 2) ? true : false ;

bool stitch = eta1 && eta2 && dphi ;

ap_uint<12> Energy1 = rctecalclusters1[j].energy ;
ap_uint<12> Energy2 = rctecalclusters2[i].energy ;

ap_uint<12> EnergySum ;
ap_uint<10> EnergySum5x5 ;

getSums(rctecalclusters1[j], rctecalclusters2[i], EnergySum, EnergySum5x5) ;

if(stitch)
{
	if(Energy1 > Energy2){
	rctecalclusters1[j].energy = EnergySum ;
	rctecalclusters1[j].et5x5 = EnergySum5x5 ;
	rctecalclusters2[i].energy = 0 ;
	}
	else {
	rctecalclusters1[j].energy = 0 ;
	rctecalclusters2[i].energy = EnergySum ;
	rctecalclusters2[i].et5x5 = EnergySum5x5 ;
	}
}

}}

}


void fillRegionsEta(rcttower R0Towers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4], rcttower R1Towers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4]){


	for(loop ieta=2; ieta<20; ieta++){

		R0Towers[ieta][8] = R1Towers[ieta][2] ;
		R0Towers[ieta][9] = R1Towers[ieta][3] ;
		R1Towers[ieta][0] = R0Towers[ieta][6] ;
		R1Towers[ieta][1] = R0Towers[ieta][7] ;
	}
}

void fillRegionsPhi(rcttower R0Towers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4], rcttower R1Towers[N_GCT_TOWERS_REGION_ETA+5][N_GCT_TOWERS_REGION_PHI+4]){


	for(loop iphi=0; iphi<10; iphi++){

		R0Towers[1][iphi] = R1Towers[2][9-iphi] ;
		R0Towers[0][iphi] = R1Towers[3][9-iphi] ;

		R1Towers[1][iphi] = R0Towers[2][9-iphi] ;
		R1Towers[0][iphi] = R0Towers[3][9-iphi] ;
	}
}


void copy_to_output(ap_uint<576> link_tmp[N_OUTPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]){



    for(loop i=0; i<N_OUTPUT_LINKS; i++) link_out[i]=link_tmp[i] ;
}

void algo_topIP1gct(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]){


    rctregion RCTRegions[N_RCT_REGIONS];


    ap_uint<576> link_tmp[N_OUTPUT_LINKS];


    for(loop i=0; i<N_RCT_REGIONS; i++){

        fillInputClusters(link_in[4*i].range(575, 0), RCTRegions[i].rctecalclusters, RCTRegions[i].Grid);
        fillInputTowers1(link_in[4*i+1].range(575, 0), RCTRegions[i].rcttowers);
        fillInputTowers2(link_in[4*i+2].range(575, 0), RCTRegions[i].rcttowers);
        fillInputTowers3(link_in[4*i+3].range(575, 0), RCTRegions[i].rcttowers);
        buildfulltowers(RCTRegions[i].rcttowers, RCTRegions[i].rctecalclusters, RCTRegions[i].Grid);
    }

    for(loop i=0; i<N_RCT_REGIONS; i=i+2){

        stitchClustersOnEta(RCTRegions[i].rctecalclusters, RCTRegions[i+1].rctecalclusters);
    }

    stitchClustersOnEta(RCTRegions[1].rctecalclusters, RCTRegions[2].rctecalclusters);
    stitchClustersOnEta(RCTRegions[5].rctecalclusters, RCTRegions[6].rctecalclusters);

    for(loop i=0; i<N_RCT_REGIONS/2; i++){

        stitchClustersOnPhi(RCTRegions[i].rctecalclusters, RCTRegions[7-i].rctecalclusters);
    }

    for(loop i=0; i<3; i++){

        fillRegionsEta(RCTRegions[i].rcttowers, RCTRegions[i+1].rcttowers);
    }
    for(loop i=0; i<3; i++){

        fillRegionsEta(RCTRegions[i+4].rcttowers, RCTRegions[i+5].rcttowers);
    }
    for(loop i=0; i<N_RCT_REGIONS/2; i++){

        fillRegionsPhi(RCTRegions[i].rcttowers, RCTRegions[7-i].rcttowers);
    }

    for(loop i=0; i<N_RCT_REGIONS; i++){

        for(loop k=0; k<N_GCT_TOWERS_REGION_ETA+5; k++){

            for(loop m=0; m<N_GCT_TOWERS_REGION_PHI+4; m++){

                RCTRegions[i].rcttowers[k][m].eta = k;
                RCTRegions[i].rcttowers[k][m].phi = m;
            }
        }
    }

    int delay_var = 0;
    for(loop i=0; i<1; i++){
        createPFClusters(RCTRegions[i].rcttowers, RCTRegions[i].pfclusters);
        delay_var += 1;
    }
    for(loop i=delay_var; i<N_RCT_REGIONS; i++){
        createPFClusters(RCTRegions[i].rcttowers, RCTRegions[i].pfclusters);
    }

    for(loop i=0; i<N_RCT_REGIONS; i++){

        RCTRegions[i].completeEG();
    }

    for(loop i=0; i<N_RCT_REGIONS; i++){

        RCTRegions[i].calculateST();
    }

    for(loop i=0; i<N_OUTPUT_LINKS; i++){

        link_tmp[i] = 0;
    }

    for(loop i=0; i<1; i++){

        RCTRegions[4*i+1].fillEGlink1n8(link_tmp[0]);
        RCTRegions[4*i+2].fillEGlink8n1(link_tmp[1]);
        RCTRegions[4*i+5].fillEGlink17(link_tmp[1], link_tmp[2]);
        RCTRegions[4*i+6].fillEGlink26(link_tmp[2], link_tmp[3]);
        RCTRegions[4*i+1].fillPFlink39(link_tmp[3], link_tmp[4]);
        RCTRegions[4*i+2].fillPFlink93(link_tmp[5], link_tmp[6]);
        RCTRegions[4*i+5].fillPFlink66(link_tmp[6], link_tmp[7]);
        RCTRegions[4*i+6].fillPFlink39(link_tmp[7], link_tmp[8]);
    }

    for(loop i=1; i<3; i++){

        RCTRegions[i].fillEGlink(link_tmp[i+8]);
        RCTRegions[i+4].fillEGlink(link_tmp[i+10]);
    }

    for(loop i=0; i<N_RCT_REGIONS; i++){

        RCTRegions[i].fillSTlink(link_tmp[i+13]);
    }

    copy_to_output(link_tmp, link_out); //ap_vld sync
}
