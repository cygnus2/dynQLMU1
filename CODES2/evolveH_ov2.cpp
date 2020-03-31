/* compute spectral weights in each flippable plaquette sector */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

extern void initState(int, int, int*);
extern void flippedHist(int, int, int*);

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void evolveH_ov2(int sector){
    extern void printconf(bool*);
    std::vector<double> initC;
    int i,ix,iy,p,q;
    int m,k,j,r;
    int nFlip[VOL+1]; 
    double specWT[VOL+1];
    double t;
    double betaR,betaI,betaM,betaTot;
    double fprof[VOL];
    int sizet, nTot;
    FILE *fptr,*fptr1; 

    /* find the relevant cartoon states with the specified number of flippable plaquettes */
    if((LX == 2) && (LY == 2)){
      std::cout<<"This routine does not include the case LX=2 and LY=2"<<std::endl; exit(0);}
    if(LY > 2){ std::cout<<"This routine does not work for LY>2"<<std::endl; exit(0); }
    
    sizet = Wind[sector].nBasis;

    /* initialize */
    //for(k=0;k<=VOL;k++) nFlip[k]=0; 

    flippedHist(sector, VOL+1, nFlip);
    
    /* compute #-basis states with k-flippable plaquettes */
    //for(k=0; k<sizet; k++){
    //	  p = Wind[sector].nflip[k];
    //	  nFlip[p]++;
    //}

    /* obtain the initial starting state */
    q=0;
    initState(sector,INIT,&q);
    if(q<0){ 
      std::cout<<"Initial state not found!"<<std::endl; exit(0); }
    else{
      std::cout<<"Starting state is basis state = "<<q<<std::endl;
      /* store the overlap of the initial state with the eigenvectors */
      for(p=0; p<sizet; p++){
        initC.push_back(Wind[sector].evecs[p*sizet+q]);
      }
    }

    /* now compute the overlap in each sector */
    fptr = fopen("spectral_WT.dat","w");
    fptr1= fopen("state_Prof.dat","w");
    for(t=Ti; t<Tf; t=t+dT){

      /* initialize spectral weights */
      for(k=VOL;k>=0;k--) specWT[k]=0.0;

      /* initialize flux profile */
      for(k=0;k<VOL;k++) fprof[k]=0.0; 

      for(k=0; k<sizet; k++){
          p = Wind[sector].nflip[k];
          betaR = 0.0; betaI = 0.0;
          for(m=0; m<sizet; m++){
	     betaR += Wind[sector].evecs[m*sizet+k]*initC[m]*cos(-Wind[sector].evals[m]*t);	  
	     betaI += Wind[sector].evecs[m*sizet+k]*initC[m]*sin(-Wind[sector].evals[m]*t);	  
	  }
          betaM = betaR*betaR + betaI*betaI;
          specWT[p] = specWT[p] + betaM;
          
          /* get the flippability profile at time t */
          for(r=0;r<VOL;r++){
	  	  if(Wind[sector].xflip[k][r]) fprof[r] += betaM; 	  
	  }
      } // close loop over basis states
      betaTot = 0.0;
      /* print the spectral weights, starting from maximal flippable plaquettes,
       * increasing with the number of non-flippable plaquettes */
      fprintf(fptr,"%lf ",t);
      for(k=VOL; k>=0; k--){
	  fprintf(fptr,"%lf ",specWT[k]);
          betaTot += specWT[k];	  
      }
      fprintf(fptr,"%lf \n",betaTot);
      /* print the flippability profile at each times */
      fprintf(fptr1,"%lf ",t);
      for(k=0;k<VOL;k++){ fprintf(fptr1,"%lf ",fprof[k]); }
       fprintf(fptr1,"\n");
    }
    fclose(fptr);
    fclose(fptr1);

    /* free memory */
    initC.clear();
}

void flippedHist(int sector, int T, int *nFlip){ 
  int k,p,sizet,tot; 
  FILE *flipH;
  if(T!=(VOL+1)) std::cout<<"Mismatch in the routine flipped_hist. Check!"<<std::endl;
  sizet = Wind[sector].nBasis;
  /* initialize */
  for(k=0;k<=VOL;k++) nFlip[k]=0;
  /* compute #-basis states with k-flippable plaquettes */
  for(k=0; k<sizet; k++){
	  p = Wind[sector].nflip[k];
	  nFlip[p]++;
  }
  /* print histogram */
  tot=0;
  flipH = fopen("flipH.dat","w");
  for(k=0;k<=VOL;k++){
    tot = tot + nFlip[k];
    fprintf(flipH,"%d %d\n",k,nFlip[k]);
  }
  fclose(flipH);
  if(tot != sizet) std::cout<<"Size mismatch: sizet="<<sizet<<" tot="<<tot<<std::endl;	
}
