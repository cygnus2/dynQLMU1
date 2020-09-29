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
void evolveH_ov2_INIT0(int sector){
  MKL_INT p,q,q1,r,k,m;
  int sizet,tsect;
  double t;
  std::vector<double> fprof(VOL);
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi;
  std::vector<double> cos00(tsect), cosPiPi(tsect), sin00(tsect), sinPiPi(tsect);
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double sum1, sum2, norm;
  FILE *fptr;

  /* initialize the starting state */
  q1 = INITq;
  std::cout<<"In function evolveH_ov2. Starting state is basis state = "<<q1<<std::endl;

  // depending on the initial state, the flipx needs to be expressed in the
  // appropriate basis. For INITq=0 [where only (0,0) and (pi,pi) contribute]
  // we have three independent bag "matrix" elements (00, 11 and 01=10).
  //00 --> sector (0,0) connecting (0,0);
  //11 --> sector (pi,pi) connecting (pi,pi)
  //01 --> sector (0,0) connecting (pi,pi)
  // Also, for the moment we note that the diagonal matrix elements will always be
  // the same as 00. This follows since the actual phase factors are +1/-1. Thus,
  // for the diagonal elements, there is always a square of the phase factor, which
  // is equal to unity. We use to save memory.
  std::vector<std::vector<double>> flipx00(VOL, std::vector<double>(tsect, 0.0));
  //std::vector<std::vector<double>> flipx11(VOL, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> flipx01(VOL, std::vector<double>(tsect, 0.0));

  // calculate matrix elements of flipx(x) in the bag basis. We need this only for
  // the momenta which have overlap with the initial state
  for(r=0; r<VOL; r++){
    for(p=0; p<sizet; p++){
     q=Wind[sector].Tflag[p]-1;
     norm=Wind[sector].Tdgen[p]/((double)VOL);
     flipx00[r][q] += Wind[sector].xflip[p][r]*norm;
     flipx01[r][q] += Wind[sector].xflip[p][r]*Wind[sector].FPiPi[p]*norm;
  }}

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]);
  }

  std::cout<<"Overlap with initial states recorded"<<std::endl;
  // compute the state profile
  fptr= fopen("state_Prof.dat","w");
  for(t=Ti; t<Tf; t=t+dT){
      // initialize the flux profile
      fprof.assign(VOL, 0.0);
      for(std::size_t ii = 0; ii < tsect; ii++){
         cos00[ii]   = cos(Wind[sector].evals_K00[ii]*t);
         cosPiPi[ii] = cos(Wind[sector].evals_KPiPi[ii]*t);
         sin00[ii]   = sin(Wind[sector].evals_K00[ii]*t);
         sinPiPi[ii] = sin(Wind[sector].evals_KPiPi[ii]*t);
      }
      for(k=0; k<tsect; k++){
        phiRE00 = 0.0; phiIM00 = 0.0; phiREPiPi=0.0; phiIMPiPi=0.0;
        for(m=0; m<tsect; m++){
          /* contributions from sector (0,0) */
          phiRE00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*cos00[m];
          phiIM00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*sin00[m];
          //phiRE00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*cos(Wind[sector].evals_K00[m]*t);
          //phiIM00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*sin(Wind[sector].evals_K00[m]*t);
          /* contributions from sector (pi,pi) */
          phiREPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*cosPiPi[m];
          phiIMPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*sinPiPi[m];
          //phiREPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*cos(Wind[sector].evals_KPiPi[m]*t);
          //phiIMPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*sin(Wind[sector].evals_KPiPi[m]*t);
        }
        sum1 = (phiRE00*phiRE00 +  phiIM00*phiIM00 + phiREPiPi*phiREPiPi + phiIMPiPi*phiIMPiPi)*inorm*inorm;
        sum2 = 2*(phiRE00*phiREPiPi + phiIM00*phiIMPiPi)*inorm*inorm*INITphasePiPi;
        //std::cout<<sum1<<" "<<sum2<<std::endl;
        for(r=0;r<VOL;r++) fprof[r] += (flipx00[r][k]*sum1 + flipx01[r][k]*sum2);
      } // close the loop over bag states
      /* print the flippability profile at each times */
      fprintf(fptr,"%lf ",t);
      for(r=0;r<VOL;r++){ fprintf(fptr,"% lf ",fprof[r]); }
      fprintf(fptr,"\n");
  }
  fclose(fptr);
  /* free memory */
  alpha00.clear(); alphaPiPi.clear(); fprof.clear();
  flipx00.clear(); flipx01.clear();
  cos00.clear(); sin00.clear(); cosPiPi.clear(); sinPiPi.clear();
}


void evolveH_ov2_INIT4(int sector){
  MKL_INT p,q,q1,r,k,m;
  int sizet,tsect;
  double t;
  std::vector<double> fprof(VOL);
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi, alphaPi0, alpha0Pi;
  std::vector<double> cos00(tsect), cosPiPi(tsect), cosPi0(tsect), cos0Pi(tsect);
  std::vector<double> sin00(tsect), sinPiPi(tsect), sinPi0(tsect), sin0Pi(tsect);
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double phiRE0Pi, phiIM0Pi, phiREPi0, phiIMPi0;
  double sum1, sum2, sum3, sum4, sum5, sum6, sum7, norm;
  FILE *fptr;

  /* initialize the starting state */
  q1 = INITq;
  std::cout<<"In function evolveH_ov2. Starting state is basis state = "<<q1<<std::endl;
  //For INIT=4, the matrix elements 00, 01, 02, 03, 12, 13, 23
  // 0 --> sector (0,0); 1 --> sector (pi,pi); 2 --> sector (pi,0); 3 --> sector (0,pi)
  //00 --> sector (0,0) connecting (0,0) [= 11 = 22 = 33]
  std::vector<std::vector<double>> flipx00(VOL, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> flipx01(VOL, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> flipx02(VOL, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> flipx03(VOL, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> flipx12(VOL, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> flipx13(VOL, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> flipx23(VOL, std::vector<double>(tsect, 0.0));

  // calculate matrix elements of flipx(x) in the bag basis. We need this only for
  // the momenta which have overlap with the initial state
  for(r=0; r<VOL; r++){
    for(p=0; p<sizet; p++){
     q=Wind[sector].Tflag[p]-1;
     norm=Wind[sector].Tdgen[p]/((double)VOL);
     flipx00[r][q] += Wind[sector].xflip[p][r]*norm;
     flipx01[r][q] += Wind[sector].xflip[p][r]*Wind[sector].FPiPi[p]*norm;
     flipx02[r][q] += Wind[sector].xflip[p][r]*Wind[sector].FPi0[p]*norm;
     flipx03[r][q] += Wind[sector].xflip[p][r]*Wind[sector].F0Pi[p]*norm;
     flipx12[r][q] += Wind[sector].xflip[p][r]*Wind[sector].FPiPi[p]*Wind[sector].FPi0[p]*norm;
     flipx13[r][q] += Wind[sector].xflip[p][r]*Wind[sector].FPiPi[p]*Wind[sector].F0Pi[p]*norm;
     flipx23[r][q] += Wind[sector].xflip[p][r]*Wind[sector].FPi0[p]*Wind[sector].F0Pi[p]*norm;
  }}

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  for(p=0; p<tsect; p++){
   alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]);
   alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]);
   alphaPi0.push_back(Wind[sector].evecs_KPi0[p*tsect+INITbag]);
   alpha0Pi.push_back(Wind[sector].evecs_K0Pi[p*tsect+INITbag]);
  }

  // compute the state profile
  fptr= fopen("state_Prof.dat","w");
  for(t=Ti; t<Tf; t=t+dT){
      // initialize the flux profile
      fprof.assign(VOL, 0.0);
      for(std::size_t ii = 0; ii < tsect; ii++){
         cos00[ii]   = cos(Wind[sector].evals_K00[ii]*t);
         cosPiPi[ii] = cos(Wind[sector].evals_KPiPi[ii]*t);
         cosPi0[ii]  = cos(Wind[sector].evals_KPi0[ii]*t);
         cos0Pi[ii]  = cos(Wind[sector].evals_K0Pi[ii]*t);
         sin00[ii]   = sin(Wind[sector].evals_K00[ii]*t);
         sinPiPi[ii] = sin(Wind[sector].evals_KPiPi[ii]*t);
         sinPi0[ii]  = sin(Wind[sector].evals_KPi0[ii]*t);
         sin0Pi[ii]  = sin(Wind[sector].evals_K0Pi[ii]*t);
      }
      for(k=0; k<tsect; k++){
        phiRE00=0.0; phiIM00=0.0; phiREPiPi=0.0; phiIMPiPi=0.0;
        phiREPi0=0.0; phiIMPi0=0.0; phiRE0Pi=0.0; phiIM0Pi=0.0;
        for(m=0; m<tsect; m++){
          /* contributions from sector (0,0) */
          phiRE00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*cos00[m];
          phiIM00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*sin00[m];
          //phiRE00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*cos(Wind[sector].evals_K00[m]*t);
          //phiIM00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*sin(Wind[sector].evals_K00[m]*t);
          /* contributions from sector (pi,pi) */
          phiREPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*cosPiPi[m];
          phiIMPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*sinPiPi[m];
          //phiREPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*cos(Wind[sector].evals_KPiPi[m]*t);
          //phiIMPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*sin(Wind[sector].evals_KPiPi[m]*t);
          /* contributions from sector (pi,0) */
          phiREPi0 += alphaPi0[m]*Wind[sector].evecs_KPi0[tsect*m+k]*cosPi0[m];
          phiIMPi0 += alphaPi0[m]*Wind[sector].evecs_KPi0[tsect*m+k]*sinPi0[m];
          //phiREPi0 += alphaPi0[m]*Wind[sector].evecs_KPi0[tsect*m+k]*cos(Wind[sector].evals_KPi0[m]*t);
          //phiIMPi0 += alphaPi0[m]*Wind[sector].evecs_KPi0[tsect*m+k]*sin(Wind[sector].evals_KPi0[m]*t);
          /* contributions from sector (0,pi) */
          phiRE0Pi += alpha0Pi[m]*Wind[sector].evecs_K0Pi[tsect*m+k]*cos0Pi[m];
          phiIM0Pi += alpha0Pi[m]*Wind[sector].evecs_K0Pi[tsect*m+k]*sin0Pi[m];
          //phiRE0Pi += alpha0Pi[m]*Wind[sector].evecs_K0Pi[tsect*m+k]*cos(Wind[sector].evals_K0Pi[m]*t);
          //phiIM0Pi += alpha0Pi[m]*Wind[sector].evecs_K0Pi[tsect*m+k]*sin(Wind[sector].evals_K0Pi[m]*t);
        }
        sum1 = (phiRE00*phiRE00 +  phiIM00*phiIM00 + phiREPiPi*phiREPiPi + phiIMPiPi*phiIMPiPi
             +  phiREPi0*phiREPi0 + phiIMPi0*phiIMPi0 + phiRE0Pi*phiRE0Pi + phiIM0Pi*phiIM0Pi)*inorm*inorm;
        sum2 = 2*(phiRE00*phiREPiPi + phiIM00*phiIMPiPi)*inorm*inorm*INITphasePiPi;
        sum3 = 2*(phiRE00*phiREPi0 + phiIM00*phiIMPi0)*inorm*inorm*INITphasePi0;
        sum4 = 2*(phiRE00*phiRE0Pi + phiIM00*phiIM0Pi)*inorm*inorm*INITphase0Pi;
        sum5 = 2*(phiREPiPi*phiREPi0 + phiIMPiPi*phiIMPi0)*inorm*inorm*INITphasePiPi*INITphasePi0;
        sum6 = 2*(phiREPiPi*phiRE0Pi + phiIMPiPi*phiIM0Pi)*inorm*inorm*INITphasePiPi*INITphase0Pi;
        sum7 = 2*(phiREPi0*phiRE0Pi + phiIMPi0*phiIM0Pi)*inorm*inorm*INITphasePi0*INITphase0Pi;
        //std::cout<<sum1<<" "<<sum2<<std::endl;
        for(r=0;r<VOL;r++) fprof[r] += (flipx00[r][k]*sum1 + flipx01[r][k]*sum2 + flipx02[r][k]*sum3 +
                      + flipx03[r][k]*sum4 + flipx12[r][k]*sum5 + flipx13[r][k]*sum6 + flipx23[r][k]*sum7);
      } // close the loop over bag states
      /* print the flippability profile at each times */
      fprintf(fptr,"%lf ",t);
      for(r=0;r<VOL;r++){ fprintf(fptr,"% lf ",fprof[r]); }
      fprintf(fptr,"\n");
  }
  fclose(fptr);
  /* free memory */
  fprof.clear();
  alpha00.clear(); alphaPiPi.clear(); alphaPi0.clear(); alpha0Pi.clear();
  flipx00.clear(); flipx01.clear(); flipx02.clear(); flipx03.clear();
  flipx12.clear(); flipx13.clear(); flipx23.clear();
  cos00.clear(); cosPiPi.clear(); cosPi0.clear(); cos0Pi.clear();
  sin00.clear(); sinPiPi.clear(); sinPi0.clear(); sin0Pi.clear();
}
