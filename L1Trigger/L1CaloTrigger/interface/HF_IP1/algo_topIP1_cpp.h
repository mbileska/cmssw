#include "algo_topIP1_h.h"
#include "bubl_sorter_h.h"

void processInputLinks(ap_uint<576> link_in[N_INPUT_LINKS], hfregion HFRegion[N_HF_REGIONS]){

// -2 because last 2 are 4x towers 
// we have initally 10 A and 10 B towers that we need to split in 2 each
// and after that 11th A and 11th B we need to split in 4
//
	for(loop j=0; j<N_INPUT_LINKS/3; j++){
		for(loop i=0; i<8; i++){
			ap_uint<8> startA = i*10 ;
			ap_uint<8> startB = startA+110 ;
			ap_uint<8> endA = startA+7 ;
			ap_uint<8> endB = startB+7 ;
			ap_uint<8> startAfb = startA+8 ;
			ap_uint<8> startBfb = startB+8 ;
			ap_uint<8> endAfb = startAfb+1 ;
			ap_uint<8> endBfb = startBfb+1 ;

			ap_uint<8> valueA1[3]; 
			ap_uint<8> valueA1s[3]; 
			ap_uint<2> valueA1fb[3]; 
			ap_uint<8> valueB1[3]; 
			ap_uint<8> valueB1s[3]; 
			ap_uint<2> valueB1fb[3]; 

			for(loop k=0; k<3; k++){
			ap_uint<8> valueA = link_in[3*j+k].range(endA,startA) ; 
			ap_uint<8> valueB = link_in[3*j+k].range(endB,startB) ; 
			ap_uint<8> valueAfb = link_in[3*j+k].range(endAfb,startAfb) ; 
			ap_uint<8> valueBfb = link_in[3*j+k].range(endBfb,startBfb) ; 
			valueA1[k] = valueA>>1 ; valueA1s[k] = valueA - valueA1[k] ; valueA1fb[k] = valueAfb ;
			valueB1[k] = valueB>>1 ; valueB1s[k] = valueB - valueB1[k] ; valueB1fb[k] = valueBfb ;
			}

			loop cntr = 2*j ; loop cntr1 = i+2 ;

			HFRegion[cntr].hftowers[cntr1][2].energy = valueA1[0] ;
			HFRegion[cntr].hftowers[cntr1][3].energy = valueA1s[0] ;
			HFRegion[cntr].hftowers[cntr1][4].energy = valueB1[0] ;
			HFRegion[cntr].hftowers[cntr1][5].energy = valueB1s[0] ;
			HFRegion[cntr].hftowers[cntr1][6].energy = valueA1[1] ;
			HFRegion[cntr].hftowers[cntr1][7].energy = valueA1s[1] ;

			HFRegion[cntr+1].hftowers[cntr1][2].energy = valueB1[1] ;
			HFRegion[cntr+1].hftowers[cntr1][3].energy = valueB1s[1] ;
			HFRegion[cntr+1].hftowers[cntr1][4].energy = valueA1[2] ;
			HFRegion[cntr+1].hftowers[cntr1][5].energy = valueA1s[2] ;
			HFRegion[cntr+1].hftowers[cntr1][6].energy = valueB1[2] ;
			HFRegion[cntr+1].hftowers[cntr1][7].energy = valueB1s[2] ;

			HFRegion[cntr].hftowers[cntr1][2].fb = valueA1fb[0] ;
			HFRegion[cntr].hftowers[cntr1][3].fb = valueA1fb[0] ;
			HFRegion[cntr].hftowers[cntr1][4].fb = valueB1fb[0] ;
			HFRegion[cntr].hftowers[cntr1][5].fb = valueB1fb[0] ;
			HFRegion[cntr].hftowers[cntr1][6].fb = valueA1fb[1] ;
			HFRegion[cntr].hftowers[cntr1][7].fb = valueA1fb[1] ;

			HFRegion[cntr+1].hftowers[cntr1][2].fb = valueB1fb[1] ;
			HFRegion[cntr+1].hftowers[cntr1][3].fb = valueB1fb[1] ;
			HFRegion[cntr+1].hftowers[cntr1][4].fb = valueA1fb[2] ;
			HFRegion[cntr+1].hftowers[cntr1][5].fb = valueA1fb[2] ;
			HFRegion[cntr+1].hftowers[cntr1][6].fb = valueB1fb[2] ;
			HFRegion[cntr+1].hftowers[cntr1][7].fb = valueB1fb[2] ;
	
		}	

		for(loop i=0; i<6; i++){
			ap_uint<8> startA = (i+4)*10 ;
			ap_uint<8> startB = startA+110 ;
			ap_uint<8> endA = startA+7 ;
			ap_uint<8> endB = startB+7 ;
			ap_uint<8> startAfb = startA+8 ;
			ap_uint<8> startBfb = startB+8 ;
			ap_uint<8> endAfb = startAfb+1 ;
			ap_uint<8> endBfb = startBfb+1 ;

			ap_uint<8> valueA1[3]; 
			ap_uint<8> valueA1s[3]; 
			ap_uint<2> valueA1fb[3]; 
			ap_uint<8> valueB1[3]; 
			ap_uint<8> valueB1s[3]; 
			ap_uint<2> valueB1fb[3]; 

			for(loop k=0; k<3; k++){
			ap_uint<8> valueA = link_in[3*j+k].range(endA,startA) ; 
			ap_uint<8> valueB = link_in[3*j+k].range(endB,startB) ; 
			ap_uint<8> valueAfb = link_in[3*j+k].range(endAfb,startAfb) ; 
			ap_uint<8> valueBfb = link_in[3*j+k].range(endBfb,startBfb) ; 
			valueA1[k] = valueA>>1 ; valueA1s[k] = valueA - valueA1[k] ; valueA1fb[k] = valueAfb ;
			valueB1[k] = valueB>>1 ; valueB1s[k] = valueB - valueB1[k] ; valueB1fb[k] = valueBfb ;
			}
			loop cntr = 2*j+12;

			HFRegion[cntr].hftowers[i][2].energy = valueA1[0] ;
			HFRegion[cntr].hftowers[i][3].energy = valueA1s[0] ;
			HFRegion[cntr].hftowers[i][4].energy = valueB1[0] ;
			HFRegion[cntr].hftowers[i][5].energy = valueB1s[0] ;
			HFRegion[cntr].hftowers[i][6].energy = valueA1[1] ;
			HFRegion[cntr].hftowers[i][7].energy = valueA1s[1] ;

			HFRegion[cntr+1].hftowers[i][2].energy = valueB1[1] ;
			HFRegion[cntr+1].hftowers[i][3].energy = valueB1s[1] ;
			HFRegion[cntr+1].hftowers[i][4].energy = valueA1[2] ;
			HFRegion[cntr+1].hftowers[i][5].energy = valueA1s[2] ;
			HFRegion[cntr+1].hftowers[i][6].energy = valueB1[2] ;
			HFRegion[cntr+1].hftowers[i][7].energy = valueB1s[2] ;

			HFRegion[cntr].hftowers[i][2].fb = valueA1fb[0] ;
			HFRegion[cntr].hftowers[i][3].fb = valueA1fb[0] ;
			HFRegion[cntr].hftowers[i][4].fb = valueB1fb[0] ;
			HFRegion[cntr].hftowers[i][5].fb = valueB1fb[0] ;
			HFRegion[cntr].hftowers[i][6].fb = valueA1fb[1] ;
			HFRegion[cntr].hftowers[i][7].fb = valueA1fb[1] ;

			HFRegion[cntr+1].hftowers[i][2].fb = valueB1fb[1] ;
			HFRegion[cntr+1].hftowers[i][3].fb = valueB1fb[1] ;
			HFRegion[cntr+1].hftowers[i][4].fb = valueA1fb[2] ;
			HFRegion[cntr+1].hftowers[i][5].fb = valueA1fb[2] ;
			HFRegion[cntr+1].hftowers[i][6].fb = valueB1fb[2] ;
			HFRegion[cntr+1].hftowers[i][7].fb = valueB1fb[2] ;
			}

		for(loop i=0; i<2; i++){
			ap_uint<8> startA = 100+i*110 ; 
			ap_uint<8> endA = startA+7 ;
			ap_uint<8> startAfb = startA+8 ;
			ap_uint<8> endAfb = startAfb+1 ;

			ap_uint<8> valueA1[3]; 
			ap_uint<8> valueA1s[3]; 
			ap_uint<8> valueA1ss[3]; 
			ap_uint<2> valueA1fb[3]; 

			for(loop k=0; k<3; k++){
			ap_uint<8> valueA = link_in[3*j+k].range(endA,startA) ; 
			ap_uint<8> valueAfb = link_in[3*j+k].range(endAfb,startAfb) ; 
			valueA1[k] = valueA>>2 ; valueA1s[k] = (valueA - 2*valueA1[k])>>1 ; 
			valueA1ss[k] = valueA - valueA1s[k] - 2*valueA1[k] ; valueA1fb[k] = valueAfb ;
			}
			
			loop cntr12 = 2*j+12 ; loop cntr13 = cntr12+1 ; 
			loop cntr1 = i+6 ;

			HFRegion[cntr12].hftowers[cntr1][2].energy = valueA1[0] ;
			HFRegion[cntr12].hftowers[cntr1][3].energy = valueA1s[0] ;
			HFRegion[cntr12].hftowers[cntr1][4].energy = valueA1ss[0] ;
			HFRegion[cntr12].hftowers[cntr1][5].energy = valueA1[0] ;
			HFRegion[cntr12].hftowers[cntr1][6].energy = valueA1[1] ;
			HFRegion[cntr12].hftowers[cntr1][7].energy = valueA1s[1] ;

			HFRegion[cntr13].hftowers[cntr1][2].energy = valueA1ss[1] ;
			HFRegion[cntr13].hftowers[cntr1][3].energy = valueA1[1] ;
			HFRegion[cntr13].hftowers[cntr1][4].energy = valueA1[2] ;
			HFRegion[cntr13].hftowers[cntr1][5].energy = valueA1s[2] ;
			HFRegion[cntr13].hftowers[cntr1][6].energy = valueA1ss[2] ;
			HFRegion[cntr13].hftowers[cntr1][7].energy = valueA1[2] ;

			HFRegion[cntr12].hftowers[cntr1][2].fb = valueA1fb[0] ;
			HFRegion[cntr12].hftowers[cntr1][3].fb = valueA1fb[0] ;
			HFRegion[cntr12].hftowers[cntr1][4].fb = valueA1fb[0] ;
			HFRegion[cntr12].hftowers[cntr1][5].fb = valueA1fb[0] ;
			HFRegion[cntr12].hftowers[cntr1][6].fb = valueA1fb[1] ;
			HFRegion[cntr12].hftowers[cntr1][7].fb = valueA1fb[1] ;

			HFRegion[cntr13].hftowers[cntr1][2].fb = valueA1fb[1] ;
			HFRegion[cntr13].hftowers[cntr1][3].fb = valueA1fb[1] ;
			HFRegion[cntr13].hftowers[cntr1][4].fb = valueA1fb[2] ;
			HFRegion[cntr13].hftowers[cntr1][5].fb = valueA1fb[2] ;
			HFRegion[cntr13].hftowers[cntr1][6].fb = valueA1fb[2] ;
			HFRegion[cntr13].hftowers[cntr1][7].fb = valueA1fb[2] ;
	
		}}	


}

void calculateST(hfregion HFRegionTMP, stower stowers[2][2]){

	for(loop i=0; i<2; i++){
		for(loop k=0; k<2; k++){
			ap_uint<12> tmp0;
			ap_uint<12> tmp1;
			ap_uint<12> tmp2;
			ap_uint<12> tmp3;
			ap_uint<12> tmp4;
			ap_uint<12> tmp5;
			ap_uint<4> a = 3*i+2 ;
			ap_uint<4> b = a+1 ;
			ap_uint<4> c = a+2 ;
			ap_uint<4> d = 3*k+2 ;
			ap_uint<4> e = d+1 ;
			ap_uint<4> f = d+2 ;

	tmp0 = HFRegionTMP.hftowers[a][d].energy + HFRegionTMP.hftowers[a][e].energy ;
	tmp1 = HFRegionTMP.hftowers[a][f].energy + HFRegionTMP.hftowers[b][d].energy ;
	tmp2 = HFRegionTMP.hftowers[b][e].energy + HFRegionTMP.hftowers[b][f].energy ;
	tmp3 = HFRegionTMP.hftowers[c][d].energy + HFRegionTMP.hftowers[c][e].energy ;
	tmp4 = tmp0 + tmp1;
	tmp5 = tmp2 + tmp3;
			stowers[i][k].energy = tmp4 +tmp5 + HFRegionTMP.hftowers[c][f].energy ;
			stowers[i][k].SeedEnergy = HFRegionTMP.hftowers[b][e].energy;
//			stowers[i][k].eta = i ;
//			stowers[i][k].phi = k ;
	}}

}

void calculateSTregions(stregion STRegion[24], stower HFSTowers1[STOWERS_ETA_R][STOWERS_PHI_R+4], stower  HFSTowers2[STOWERS_ETA_R][STOWERS_PHI_R+4], stower HFSTowers3[STOWERS_ETA_R][STOWERS_PHI_R+4]){

HFSTowers1[0][0] = STRegion[11].stowers[0][0] ;
HFSTowers1[1][0] = STRegion[11].stowers[1][0] ;
HFSTowers1[2][0] = STRegion[23].stowers[0][0] ;
HFSTowers1[3][0] = STRegion[23].stowers[1][0] ;

HFSTowers1[0][1] = STRegion[11].stowers[0][1] ;
HFSTowers1[1][1] = STRegion[11].stowers[1][1] ;
HFSTowers1[2][1] = STRegion[23].stowers[0][1] ;
HFSTowers1[3][1] = STRegion[23].stowers[1][1] ;

for(loop i=0; i<5; i++){
for(loop k=0; k<2; k++){
HFSTowers1[0][2+k+2*i] = STRegion[i].stowers[0][k] ;
HFSTowers1[1][2+k+2*i] = STRegion[i].stowers[1][k] ;
HFSTowers1[2][2+k+2*i] = STRegion[12+i].stowers[0][k] ;
HFSTowers1[3][2+k+2*i] = STRegion[12+i].stowers[1][k] ;
}}
for(loop i=0; i<5; i++){
for(loop k=0; k<2; k++){
HFSTowers2[0][2+k+2*i] = STRegion[i+3].stowers[0][k] ;
HFSTowers2[1][2+k+2*i] = STRegion[i+3].stowers[1][k] ;
HFSTowers2[2][2+k+2*i] = STRegion[15+i].stowers[0][k] ;
HFSTowers2[3][2+k+2*i] = STRegion[15+i].stowers[1][k] ;
}}
for(loop i=0; i<5; i++){
for(loop k=0; k<2; k++){
HFSTowers3[0][k+2*i] = STRegion[i+7].stowers[0][k] ;
HFSTowers3[1][k+2*i] = STRegion[i+7].stowers[1][k] ;
HFSTowers3[2][k+2*i] = STRegion[19+i].stowers[0][k] ;
HFSTowers3[3][k+2*i] = STRegion[19+i].stowers[1][k] ;
}}

HFSTowers3[0][10] = STRegion[0].stowers[0][0] ;
HFSTowers3[1][10] = STRegion[0].stowers[1][0] ;
HFSTowers3[2][10] = STRegion[12].stowers[0][0] ;
HFSTowers3[3][10] = STRegion[12].stowers[1][0] ;

HFSTowers3[0][11] = STRegion[0].stowers[0][1] ;
HFSTowers3[1][11] = STRegion[0].stowers[1][1] ;
HFSTowers3[2][11] = STRegion[12].stowers[0][1] ;
HFSTowers3[3][11] = STRegion[12].stowers[1][1] ;

	for (loop i = 0; i < STOWERS_ETA_R; i++) {
	        for (loop j = 0; j < STOWERS_PHI_R+4; j++) {
		HFSTowers1[i][j].eta = i ; HFSTowers1[i][j].phi = j ;
		HFSTowers2[i][j].eta = i ; HFSTowers2[i][j].phi = j ;
		HFSTowers3[i][j].eta = i ; HFSTowers3[i][j].phi = j ;
	}}


}

void fillOutputLink(PFcluster pfclusters[16], ap_uint<576> &linkOut){
	linkOut = 0;

	ap_uint<10> start ;
	ap_uint<10> end ;

	for(loop j=0; j<6; j++){
		start=j*64 ; end=start+63;
		linkOut.range(end, start) = pfclusters[j].getPFcluster() ;
	    }
}


void copy_to_output(ap_uint<576> link_tmp[N_OUTPUT_LINKS_CL1+N_OUTPUT_LINKS_MIX], ap_uint<576> link_out_copy[N_OUTPUT_LINKS_CL1+N_OUTPUT_LINKS_MIX]){
	for(loop i=0; i<N_OUTPUT_LINKS_CL1+N_OUTPUT_LINKS_MIX; i++)
		link_out_copy[i]=link_tmp[i];
	}


void algo_topIP1(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS_CL1+N_OUTPUT_LINKS_MIX]){

	// we define the towers 12x72 eta phi and create 24 regions
	// 6x6 after adding 2 overlap becomes 10x10, where eta additions
	// and phi are from correspondent towers

	ap_uint<576> link_tmp[N_OUTPUT_LINKS_CL1+N_OUTPUT_LINKS_MIX] ;
	ap_uint<576> link_out_copy[N_OUTPUT_LINKS_CL1+N_OUTPUT_LINKS_MIX];

	for(loop i=0; i<N_OUTPUT_LINKS_CL1+N_OUTPUT_LINKS_MIX; i++){
		link_tmp[i] = 0 ;
		link_out_copy[i] = 0 ;
		}

	hfregion HFRegion[N_HF_REGIONS] ; // for PF clusters
	processInputLinks(link_in, HFRegion) ;


	for(loop i=0; i<8; i++){
	HFRegion[0].hftowers[i+2][0] = HFRegion[11].hftowers[i+2][6] ;
	HFRegion[0].hftowers[i+2][1] = HFRegion[11].hftowers[i+2][7] ;
	HFRegion[12].hftowers[i][0] = HFRegion[23].hftowers[2][6] ;
	HFRegion[12].hftowers[i][1] = HFRegion[23].hftowers[2][7] ;
	for(loop k=0; k<11 ; k++){
	HFRegion[k+1].hftowers[i+2][0] = HFRegion[k].hftowers[i+2][6] ;
	HFRegion[k+1].hftowers[i+2][1] = HFRegion[k].hftowers[i+2][7] ;
	HFRegion[k+13].hftowers[i][0] = HFRegion[k+12].hftowers[2][6] ;
	HFRegion[k+13].hftowers[i][1] = HFRegion[k+12].hftowers[2][7] ;
	HFRegion[k].hftowers[i+2][8] = HFRegion[k+1].hftowers[i+2][0] ;
	HFRegion[k].hftowers[i+2][9] = HFRegion[k+1].hftowers[i+2][1] ;
	HFRegion[k+12].hftowers[i][8] = HFRegion[k+13].hftowers[2][0] ;
	HFRegion[k+12].hftowers[i][9] = HFRegion[k+13].hftowers[2][1] ;
	}
	HFRegion[11].hftowers[i+2][8] = HFRegion[0].hftowers[i+2][0] ;
	HFRegion[11].hftowers[i+2][9] = HFRegion[0].hftowers[i+2][1] ;
	HFRegion[23].hftowers[i][8] = HFRegion[12].hftowers[2][0] ;
	HFRegion[23].hftowers[i][9] = HFRegion[12].hftowers[2][1] ;
	}

// we create arrays with 2 extra phi's on each side for later overlap
// for eta we have 2 first or 2 last etas ==0 , dont need to fill
// for phi we need to use pattern ...70 71 0 1 therefore we 
// create extra loops for this to include wedge 17 or 0 that are
// out of sequence , HFRegions 0...11 low eta 12...23 high eta
// we use HFRegions to calculate 2x2 ST in each region 
//
	stregion STRegion[N_HF_REGIONS] ; // for ST 
	pfregion PFRegion[N_HF_REGIONS] ; // for pf clusters

	for(loop i=0; i<N_HF_REGIONS; i++){
		stower stowers[2][2];
		calculateST(HFRegion[i], stowers);
		STRegion[i].stowers[0][0]=stowers[0][0];
		STRegion[i].stowers[0][1]=stowers[0][1];
		STRegion[i].stowers[1][0]=stowers[1][0];
		STRegion[i].stowers[1][1]=stowers[1][1];
		}

	sums HFSums ; // only one variable for now
	createSums(STRegion, HFSums) ;


	// PF clusters from towers
	for(loop k=0; k<N_HF_REGIONS; k++){
		for(loop i=0; i<N_HF_TOWERS_REGION_ETA+4; i++){
			for(loop j=0; j<N_HF_TOWERS_REGION_PHI+4; j++){
			HFRegion[k].hftowers[i][j].eta = i ; 
			HFRegion[k].hftowers[i][j].phi = j ; 
			}}
		getPFClusters(HFRegion[k].hftowers,PFRegion[k].pfclusters) ;
		}

	for(loop k=0; k<N_HF_REGIONS/4; k++){
		PFcluster pfclustersIn[16] ;
		for(loop i=0; i<4; i++){
			pfclustersIn[i] = PFRegion[2*k].pfclusters[i] ;
			pfclustersIn[i+4] = PFRegion[2*k+1].pfclusters[i] ;
			pfclustersIn[i+8] = PFRegion[2*k+12].pfclusters[i] ;
			pfclustersIn[i+12] = PFRegion[2*k+13].pfclusters[i] ;
			}
		PFcluster pfclustersOut[16] ;
		bubl_sorter(pfclustersIn, pfclustersOut);
		fillOutputLink(pfclustersOut, link_tmp[k]) ;
	for(int i=0;i<15;i++){cout<<"cluster "<< i << " " << pfclustersIn[i].energy <<endl;}
	for(int i=0;i<15;i++){cout<<"cluster sorted "<< i << " " << pfclustersOut[i].energy <<endl;}
		}

	stower HFSTowers1[STOWERS_ETA_R][STOWERS_PHI_R+4] ; // for jets and taus
	stower HFSTowers2[STOWERS_ETA_R][STOWERS_PHI_R+4] ; // for jets and taus
	stower HFSTowers3[STOWERS_ETA_R][STOWERS_PHI_R+4] ; // for jets and taus
	stower HFSTowersJets1[STOWERS_ETA_R][STOWERS_PHI_R+4] ; // for jets and taus
	stower HFSTowersJets2[STOWERS_ETA_R][STOWERS_PHI_R+4] ; // for jets and taus
	stower HFSTowersJets3[STOWERS_ETA_R][STOWERS_PHI_R+4] ; // for jets and taus
	

	tau HFTaus1[5];
	jet HFJets1[5];
	tau HFTaus2[5];
	jet HFJets2[5];
	tau HFTaus3[5];
	jet HFJets3[5];


	calculateSTregions(STRegion, HFSTowers1, HFSTowers2, HFSTowers3) ;

	for(loop i=0; i<STOWERS_ETA_R; i++){
		for(loop k=0; k<STOWERS_PHI_R+4; k++){
			HFSTowersJets1[i][k] = HFSTowers1[i][k];
			HFSTowersJets2[i][k] = HFSTowers2[i][k];
			HFSTowersJets3[i][k] = HFSTowers3[i][k];
		}
	}

	createTaus(HFSTowers1, HFTaus1) ;
	createJets(HFSTowersJets1, HFJets1) ;
	createTaus(HFSTowers2, HFTaus2) ;
	createJets(HFSTowersJets2, HFJets2) ;
	createTaus(HFSTowers3, HFTaus3) ;
	createJets(HFSTowersJets3, HFJets3) ;

	tau HFTausin[16];
	jet HFJetsin[16];
	tau HFTausout[16];
	jet HFJetsout[16];
	for (int i = 0; i < 5; i++) {
	    HFTausin[i]     = HFTaus1[i];
	    HFTausin[i+5]   = HFTaus2[i];
	    HFTausin[i+10]  = HFTaus3[i];
	}

	for (int i = 0; i < 5; i++) {
		    HFJetsin[i]     = HFJets1[i];
		    HFJetsin[i+5]   = HFJets2[i];
		    HFJetsin[i+10]  = HFJets3[i];
		}


	bubl_sorter(HFTausin, HFTausout);
	bubl_sorter(HFJetsin, HFJetsout);

#ifndef __SYNTHESIS__
	for(int i=0;i<15;i++){cout<<"tau "<< i << " " << HFTausin[i].energy <<endl;}
	for(int i=0;i<15;i++){cout<<"sorted tau "<< i << " " << HFTausout[i].energy<<endl;}
	for(int i=0;i<15;i++){cout<<"jets "<< i << " "<< HFJetsin[i].energy<<endl;}
	for(int i=0;i<15;i++){cout<<"sorted jets " << i << " " << HFJetsout[i].energy<<endl;}
#endif

	ap_uint<10> start ;
	ap_uint<10> end ;

	for(loop j=0; j<6; j++){
		start=j*64 ; end=start+63; 
		link_tmp[6].range(end, start) = HFJetsout[j].getJet() ;
		link_tmp[7].range(end, start) = HFTausout[j].getTau() ;
		link_tmp[8].range(end, start) = HFSums.getSum() ;

		}

	copy_to_output(link_tmp, link_out) ;
//#ifndef __SYNTHESIS__
//for (int j = 0; j < N_SECTORS_PF; j++) {
//    std::cout << "LINK_out1[" << j << "]: "
//              << link_out1[j].to_string(16)
//              << std::endl;
//}
//#endif
//sync_output_links(link_out_copy, link_out);
//#ifndef __SYNTHESIS__
//for (int j = 0; j < N_SECTORS_PF; j++) {
//    std::cout << "LINK_out[" << j << "]: "
//              << link_out[j].to_string(16)
//              << std::endl;
//}
//#endif

}
