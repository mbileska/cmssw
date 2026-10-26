#include "algo_topIP1_h.h"


void getseedpositionPF(hftower towers[N_HF_TOWERS_REGION_ETA+4][N_HF_TOWERS_REGION_PHI+4], hftower &Seed){
	hftower Seed1 ; hftower Seed2; // this helps very little but some

	for(loop i=0; i<N_HF_TOWERS_REGION_ETA+4; i++){
		for (loop j=0;j<5 ;j++){
		if (towers[i][j].energy >= Seed1.energy) Seed1 = towers[i][j] ; }
		for (loop j=5; j<N_HF_TOWERS_REGION_PHI+4;j++){
		if (towers[i][j].energy >= Seed2.energy) Seed2 = towers[i][j] ; }
	}
		if (Seed1.energy >= Seed2.energy) Seed = Seed1; else Seed = Seed2 ;
}


void getseedpositionTau(stower stowers[STOWERS_ETA_R][STOWERS_PHI_R+4], stower &Seed){

	for(loop i=0; i<STOWERS_ETA_R; i++){
		for(loop j=2; j<STOWERS_PHI_R+2; j++){
		ap_uint<12> tmp = stowers[i][j].energy ;
			if (tmp >= Seed.energy) Seed = stowers[i][j] ;
	}}
}

void getseedpositionJet(stower stowers[STOWERS_ETA_R+4][STOWERS_PHI_R+4], stower &Seed){

	for(loop i=0; i<STOWERS_ETA_R+4; i++){
		for(loop j=0; j<STOWERS_PHI_R+4; j++){
		ap_uint<12> tmp = stowers[i][j].energy ;
			if (tmp >= Seed.energy) Seed = stowers[i][j] ;
		}}
}

void getclusterPF(hftower towers[N_HF_TOWERS_REGION_ETA+4][N_HF_TOWERS_REGION_PHI+4], hftower Seed, PFcluster &output) {


	ap_uint<5> eta = Seed.eta ;
	ap_uint<7> phi = Seed.phi ;

	ap_uint<5> etamin; if( eta == 0) etamin = 0; else etamin = eta-1 ;
	ap_uint<5> etamax = eta+1 ;
	ap_uint<7> phimin; if( phi == 0) phimin = 0; else phimin = phi-1 ;
	ap_uint<7> phimax = phi+1 ;

	ap_uint<12> tmpValue = 0  ;

	for(loop i=0; i<N_HF_TOWERS_REGION_ETA+4; i++){
		for(loop j=0; j<N_HF_TOWERS_REGION_PHI+4; j++){
	if(towers[i][j].eta >= etamin && towers[i][j].eta <= etamax && towers[i][j].phi >= phimin && towers[i][j].phi <= phimax )
			{ap_uint<12> energy = towers[i][j].energy ;
			tmpValue +=  energy ;}
//			tmpValue +=  energy ;
//			towers[i][j].energy = 0 ;}
	}}


//cout << " tmp " << tmpValue << " eta " << Seed.eta << " phi " <<  Seed.phi <<  endl ;
	if(Seed.eta < 2 || Seed.eta > 7 || Seed.phi < 2 || Seed.phi > 7) tmpValue  = 0 ;
//cout << " tmp1 " << tmpValue << endl ;
	output.energy = tmpValue ;
	output.eta = eta ;
	output.phi = phi ;
}

void zerrotowersPF(hftower towers[N_HF_TOWERS_REGION_ETA+4][N_HF_TOWERS_REGION_PHI+4], hftower Seed) {
	ap_uint<5> eta = Seed.eta ;
	ap_uint<7> phi = Seed.phi ;

	ap_uint<5> etamin; if( eta == 0) etamin = 0; else etamin = eta-1 ;
	ap_uint<5> etamax = eta+1 ;
	ap_uint<7> phimin; if( phi == 0) phimin = 0; else phimin = phi-1 ;
	ap_uint<7> phimax = phi+1 ;

	for(loop i=0; i<N_HF_TOWERS_REGION_ETA+4; i++){
		for(loop j=0; j<N_HF_TOWERS_REGION_PHI+4; j++){
	if(towers[i][j].eta >= etamin && towers[i][j].eta <= etamax && towers[i][j].phi >= phimin && towers[i][j].phi <= phimax )
			towers[i][j].energy = 0 ;
	}}
}

void getPFClusters(hftower towers[N_HF_TOWERS_REGION_ETA+4][N_HF_TOWERS_REGION_PHI+4], PFcluster pfclusters[N_HF_PFCLUSTERS_REGION]) {

	for (loop i=0;i<N_HF_PFCLUSTERS_REGION;i++){
		hftower  Seed ;
		getseedpositionPF(towers, Seed) ;
		getclusterPF(towers, Seed, pfclusters[i]);
		zerrotowersPF(towers, Seed);
	}
}


void createSums(stregion region[N_HF_REGIONS] , sums& sum) {
	ap_uint<16> energySum = 0;

	for (loop i = 0; i < N_HF_REGIONS; i++) {
		ap_uint<16> tmp = region[i].stowers[0][0].energy + region[i].stowers[0][1].energy;
		ap_uint<16> tmp2 = region[i].stowers[1][0].energy + region[i].stowers[1][1].energy;
		ap_uint<16> tmp3 = tmp +tmp2;
	energySum += tmp3;
	}

	if (energySum >= 4095)
		energySum = 4095;
		sum.energy = energySum;
}

void gettau(stower stowers[STOWERS_ETA_R][STOWERS_PHI_R+4], stower Seed, tau &output) {
	ap_uint<4> eta = Seed.eta ;
	ap_uint<5> phi = Seed.phi ;

	ap_uint<12> tmpValue = 0 ;

	for (loop i = 0; i < STOWERS_ETA_R; i++) {
	        for (loop j = 2; j < STOWERS_PHI_R+2; j++) {
	if( i == eta && j == phi ){
	            tmpValue = stowers[i][j].energy;
	 //           stowers[i][j].energy = 0;
	        }
	    }}

	    output.energy = tmpValue;
	    output.eta    = eta;
	    output.phi    = phi;
	}

void zerrotau(stower stowers[STOWERS_ETA_R][STOWERS_PHI_R+4], stower Seed) {
	ap_uint<4> eta = Seed.eta ;
	ap_uint<5> phi = Seed.phi ;

	ap_uint<12> tmpValue = 0 ;

	for (loop i = 0; i < STOWERS_ETA_R; i++) {
	        for (loop j = 2; j < STOWERS_PHI_R+2; j++) {
	if( i == eta && j == phi ){
	            stowers[i][j].energy = 0;
	        }
	    }}
	}


void createTaus(stower stowers[STOWERS_ETA_R][STOWERS_PHI_R+4], tau taus[5]) {

	for (int i=0; i<5; i++){
		stower  Seed ;
		getseedpositionTau(stowers, Seed);
		gettau(stowers,Seed,taus[i]);
		zerrotau(stowers,Seed);
		}
}

void getJet(stower stowers[STOWERS_ETA_R+4][STOWERS_PHI_R+4], stower Seed, jet &output) {

	ap_uint<4> eta = Seed.eta ;
	ap_uint<5> phi = Seed.phi ;

	ap_uint<12> tmpValue = 0  ;

	ap_uint<4> etamin; if( eta == 0) etamin = 0; else etamin = eta-1 ;
	ap_uint<4> etamax = eta+1 ;
	ap_uint<5> phimin; if( phi == 0) phimin = 0; else phimin = phi-1 ;
	ap_uint<5> phimax = phi+1 ;

	for(loop i=0; i<STOWERS_ETA_R+4; i++){
		for(loop j=0; j<STOWERS_PHI_R+4; j++){
//		   ap_uint<12> energy = stowers[i][j].energy ;
	if(stowers[i][j].eta >= etamin && stowers[i][j].eta <= etamax && stowers[i][j].phi >= phimin && stowers[i][j].phi <= phimax )
		{ap_uint<12> energy = stowers[i][j].energy ;
		tmpValue += energy ;}
//	stowers[i][j].energy = 0 ;}

		}}

	if(eta < 2 || eta > 5 || phi < 2 || phi > 9) tmpValue  = 0 ;
	output.energy = tmpValue ;
	output.eta = eta ;
	output.phi = phi ;
 }

void zerrostowersJet(stower stowers[STOWERS_ETA_R+4][STOWERS_PHI_R+4], stower Seed) {

	ap_uint<4> eta = Seed.eta ;
	ap_uint<5> phi = Seed.phi ;

	ap_uint<12> tmpValue = 0  ;

	ap_uint<4> etamin; if( eta == 0) etamin = 0; else etamin = eta-1 ;
	ap_uint<4> etamax = eta+1 ;
	ap_uint<5> phimin; if( phi == 0) phimin = 0; else phimin = phi-1 ;
	ap_uint<5> phimax = phi+1 ;

	for(loop i=0; i<STOWERS_ETA_R+4; i++){
		for(loop j=0; j<STOWERS_PHI_R+4; j++){
	if(stowers[i][j].eta >= etamin && stowers[i][j].eta <= etamax && stowers[i][j].phi >= phimin && stowers[i][j].phi <= phimax )
		stowers[i][j].energy = 0 ;

		}}

}

void createJets(stower stowers[STOWERS_ETA_R][STOWERS_PHI_R+4], jet jets[5]) {
	stower stowersExtended[STOWERS_ETA_R+4][STOWERS_PHI_R+4] ;

	for(loop i=0; i<STOWERS_ETA_R; i++){
		for(loop j=0; j<STOWERS_PHI_R+4; j++){
			stowersExtended[i+2][j] = stowers[i][j] ;
			}}

	for(loop i=0; i<STOWERS_ETA_R+4; i++){
		for(loop j=0; j<STOWERS_PHI_R+4; j++){
			stowersExtended[i][j].eta = i ; stowersExtended[i][j].phi = j ;
			}}


	for (int i=0;i<5;i++){
		stower  Seed ;
		getseedpositionJet(stowersExtended, Seed) ;
		getJet(stowersExtended,Seed,jets[i]) ;
		zerrostowersJet(stowersExtended, Seed) ;
	}
}
	
