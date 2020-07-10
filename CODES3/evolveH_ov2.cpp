/* compute spectral weights in each flippable plaquette sector */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

//extern void flippedHist(int, int, int*);
// eigenstate: |w_k> = \sum_l B_l |b_l>, |b_l> bag state
void evolveH_ov2(int sector){
  MKL_INT p,q,q1,r,k,m;
  int sizet,tsect;
  double t;
  std::vector<double> fprof(VOL);
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi, alphaPi0, alpha0Pi;
  std::vector<std::vector<double>> flipx(VOL, std::vector<double>(tsect, 0.0));
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double phiRE0Pi, phiIM0Pi, phiREPi0, phiIMPi0;
  double sum1, sum2;
  FILE *fptr;

  /* initialize the starting state */
  if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
  q1 = INITq;
  std::cout<<"In function evolveH_ov2. Starting state is basis state = "<<q1<<std::endl;


  // calculate flipx(x) in the bag basis: true for all momenta sectors
  for(r=0; r<VOL; r++){
    for(p=0; p<sizet; p++){
     flipx[r][Wind[sector].Tflag[p]-1] += Wind[sector].xflip[p][r]*Wind[sector].Tdgen[p]/((double)VOL);
  }}

  std::cout<<"Test example"<<std::endl;
  for(r=0; r<VOL; r++){
    for(p=0; p<tsect; p++){
      std::cout<<flipx[r][p]<<" ";
    }
    std::cout<<std::endl;
  }

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  if(INIT==0){
    for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]*inorm);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]*inorm);
    }
  }
  else if(INIT==4){
    for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]*inorm);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]*inorm);
       alpha0Pi.push_back(Wind[sector].evecs_K0Pi[p*tsect+INITbag]*inorm);
       alphaPi0.push_back(Wind[sector].evecs_KPi0[p*tsect+INITbag]*inorm);
    }
  }
  std::cout<<"Overlap with initial states recorded"<<std::endl;
  // compute the state profile
  fptr= fopen("state_Prof.dat","w");
  if(INIT==0){ // for initial condition = 0 with fully flippable plaquettes
    for(t=Ti; t<Tf; t=t+dT){
      // initialize the flux profile
      fprof.assign(VOL, 0.0);
      for(k=0; k<tsect; k++){
        phiRE00 = 0.0; phiIM00 = 0.0; phiREPiPi=0.0; phiIMPiPi=0.0;
        for(m=0; m<tsect; m++){
          /* contributions from sector (0,0) */
          phiRE00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*cos(Wind[sector].evals_K00[m]*t);
          phiIM00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*sin(Wind[sector].evals_K00[m]*t);
          /* contributions from sector (pi,pi) */
          phiREPiPi += alphaPiPi[m]*INITphasePiPi*Wind[sector].evecs_KPiPi[tsect*m+k]*cos(Wind[sector].evals_KPiPi[m]*t);
          phiIMPiPi += alphaPiPi[m]*INITphasePiPi*Wind[sector].evecs_KPiPi[tsect*m+k]*sin(Wind[sector].evals_KPiPi[m]*t);
        }
        sum1 = phiRE00*phiRE00 +  phiIM00*phiIM00 + phiREPiPi*phiREPiPi + phiIMPiPi*phiIMPiPi;
        for(r=0;r<VOL;r++) fprof[r] += flipx[r][k]*sum1;
      } // close the loop over bag states
      /* print the flippability profile at each times */
      fprintf(fptr,"%lf ",t);
      for(r=0;r<VOL;r++){ fprintf(fptr,"% lf ",fprof[r]); }
      fprintf(fptr,"\n");
    }
  }
  else if(INIT==4){ // for initial condition = 4
    std::cout<<"Come back later! "<<std::endl;
  }
  fclose(fptr);
  /* free memory */
  alpha00.clear(); alphaPiPi.clear(); alpha0Pi.clear(); alphaPi0.clear();
  flipx.clear(); fprof.clear();
}

/*
void flippedHist(int sector, int T, int *nFlip){
  int k,p,sizet,tot;
  FILE *flipH;
  if(T!=(VOL+1)) std::cout<<"Mismatch in the routine flipped_hist. Check!"<<std::endl;
  sizet = Wind[sector].nBasis;
  // initialize
  for(k=0;k<=VOL;k++) nFlip[k]=0;
  // compute #-basis states with k-flippable plaquettes
  for(k=0; k<sizet; k++){
	  p = Wind[sector].nflip[k];
	  nFlip[p]++;
  }
  // print histogram
  tot=0;
  flipH = fopen("flipH.dat","w");
  for(k=0;k<=VOL;k++){
    tot = tot + nFlip[k];
    fprintf(flipH,"%d %d\n",k,nFlip[k]);
  }
  fclose(flipH);
  if(tot != sizet) std::cout<<"Size mismatch: sizet="<<sizet<<" tot="<<tot<<std::endl;
}*/
