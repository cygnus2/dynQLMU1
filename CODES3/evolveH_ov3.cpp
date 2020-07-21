/* compute the expectation values of Ey(x) if INIT=4, dEy(x) if INIT=0 */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// The notation and conventions followed in this routine are the same as evolveH_ov2.cpp
// INIT=0 ==> dEy is non-zero, Ey is zero; ONLY calculate dEy here!
void evolveH_ov3_INIT0(int sector){
  MKL_INT p,q,q1,r,k,m;
  int sizet,tsect;
  double t;
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi;
  std::vector<double> cos00(tsect), cosPiPi(tsect), sin00(tsect), sinPiPi(tsect);
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double sum1, sum2, norm;
  std::vector<double> dEyProf(LX);
  FILE *fptr1;

  /* initialize the starting state */
  q1 = INITq;
  std::cout<<"In function evolveH_ov3, computing dEy. Starting state is basis state = "<<q1<<std::endl;

  std::vector<std::vector<double>> dEy00(LX, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> dEy01(LX, std::vector<double>(tsect, 0.0));
  // calculate matrix elements of dEy(x) in the bag basis.
  for(r=0; r<LX; r++){
    for(p=0; p<sizet; p++){
     q=Wind[sector].Tflag[p]-1;
     norm=Wind[sector].Tdgen[p]/((double)VOL);
     dEy00[r][q] += Wind[sector].dEy[p][r]*norm;
     dEy01[r][q] += Wind[sector].dEy[p][r]*Wind[sector].FPiPi[p]*norm;
  }}

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]);
  }

  // compute the state profile
  fptr1 = fopen("dEyProf.dat","w");
  for(t=Ti; t<Tf; t=t+dT){
      // initialize the flux profile
      dEyProf.assign(LX, 0.0);
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
          /* contributions from sector (pi,pi) */
          phiREPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*cosPiPi[m];
          phiIMPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*sinPiPi[m];
        }
        sum1 = (phiRE00*phiRE00 +  phiIM00*phiIM00 + phiREPiPi*phiREPiPi + phiIMPiPi*phiIMPiPi)*inorm*inorm;
        sum2 = 2*(phiRE00*phiREPiPi + phiIM00*phiIMPiPi)*inorm*inorm*INITphasePiPi;
        //std::cout<<sum1<<" "<<sum2<<std::endl;
        for(r=0;r<LX;r++){
           dEyProf[r] += (dEy00[r][k]*sum1  + dEy01[r][k]*sum2);
        }
      } // close the loop over bag states
      /* print dEy profile at each times */
      fprintf(fptr1,"%lf ",t);
      for(r=0;r<LX;r++){ fprintf(fptr1,"% lf ",dEyProf[r]); }
      fprintf(fptr1,"\n"); 
  }
  fclose(fptr1); 
  /* free memory */
  dEyProf.clear(); 
  alpha00.clear(); alphaPiPi.clear();
  dEy00.clear(); dEy01.clear(); 
  cos00.clear(); sin00.clear(); cosPiPi.clear(); sinPiPi.clear();
}

// INIT=4 ==> Ey is non-zero, dEy is zero; ONLY calculate Ey here!
void evolveH_ov3_INIT4(int sector){
  MKL_INT p,q,q1,r,k,m;
  int sizet,tsect;
  double t;
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi, alphaPi0, alpha0Pi;
  std::vector<double> cos00(tsect), cosPiPi(tsect), cosPi0(tsect), cos0Pi(tsect);
  std::vector<double> sin00(tsect), sinPiPi(tsect), sinPi0(tsect), sin0Pi(tsect);
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double phiRE0Pi, phiIM0Pi, phiREPi0, phiIMPi0;
  double sum1, sum2, sum3, sum4, sum5, sum6, sum7, norm;
  std::vector<double> EyProf(LX);
  FILE *fptr1;

  /* initialize the starting state */
  q1 = INITq;
  std::cout<<"In function evolveH_ov3, computing Ey. Starting state is basis state = "<<q1<<std::endl;

  std::vector<std::vector<double>> Ey00(LX, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> Ey01(LX, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> Ey02(LX, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> Ey03(LX, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> Ey12(LX, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> Ey13(LX, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> Ey23(LX, std::vector<double>(tsect, 0.0));
  // calculate matrix elements of Ey(x) in the bag basis.
  for(r=0; r<LX; r++){
    for(p=0; p<sizet; p++){
     q=Wind[sector].Tflag[p]-1;
     norm=Wind[sector].Tdgen[p]/((double)VOL);
     Ey00[r][q] += Wind[sector].Ey[p][r]*norm;
     Ey01[r][q] += Wind[sector].Ey[p][r]*Wind[sector].FPiPi[p]*norm;
     Ey02[r][q] += Wind[sector].Ey[p][r]*Wind[sector].FPi0[p]*norm;
     Ey03[r][q] += Wind[sector].Ey[p][r]*Wind[sector].F0Pi[p]*norm;
     Ey12[r][q] += Wind[sector].Ey[p][r]*Wind[sector].FPiPi[p]*Wind[sector].FPi0[p]*norm;
     Ey13[r][q] += Wind[sector].Ey[p][r]*Wind[sector].FPiPi[p]*Wind[sector].F0Pi[p]*norm;
     Ey23[r][q] += Wind[sector].Ey[p][r]*Wind[sector].FPi0[p]*Wind[sector].F0Pi[p]*norm;
  }}

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]);
       alphaPi0.push_back(Wind[sector].evecs_KPi0[p*tsect+INITbag]);
       alpha0Pi.push_back(Wind[sector].evecs_K0Pi[p*tsect+INITbag]);
  }

  // compute the state profile
  fptr1 = fopen("EyProf.dat","w");
  for(t=Ti; t<Tf; t=t+dT){
      // initialize the flux profile
      EyProf.assign(LX, 0.0);
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
        phiRE00 = 0.0; phiIM00 = 0.0; phiREPiPi=0.0; phiIMPiPi=0.0;
        phiREPi0=0.0; phiIMPi0=0.0; phiRE0Pi=0.0; phiIM0Pi=0.0;
        for(m=0; m<tsect; m++){
          /* contributions from sector (0,0) */
          phiRE00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*cos00[m];
          phiIM00 += alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*sin00[m];
          /* contributions from sector (pi,pi) */
          phiREPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*cosPiPi[m];
          phiIMPiPi += alphaPiPi[m]*Wind[sector].evecs_KPiPi[tsect*m+k]*sinPiPi[m];
          /* contributions from sector (pi,0) */
          phiREPi0 += alphaPi0[m]*Wind[sector].evecs_KPi0[tsect*m+k]*cosPi0[m];
          phiIMPi0 += alphaPi0[m]*Wind[sector].evecs_KPi0[tsect*m+k]*sinPi0[m];
          /* contributions from sector (0,pi) */
          phiRE0Pi += alpha0Pi[m]*Wind[sector].evecs_K0Pi[tsect*m+k]*cos0Pi[m];
          phiIM0Pi += alpha0Pi[m]*Wind[sector].evecs_K0Pi[tsect*m+k]*sin0Pi[m];
        }
        sum1 = (phiRE00*phiRE00 +  phiIM00*phiIM00 + phiREPiPi*phiREPiPi + phiIMPiPi*phiIMPiPi
             +  phiREPi0*phiREPi0 + phiIMPi0*phiIMPi0 + phiRE0Pi*phiRE0Pi + phiIM0Pi*phiIM0Pi)*inorm*inorm;
        sum2 = 2*(phiRE00*phiREPiPi + phiIM00*phiIMPiPi)*inorm*inorm*INITphasePiPi;
        sum3 = 2*(phiRE00*phiREPi0 + phiIM00*phiIMPi0)*inorm*inorm*INITphasePi0;
        sum4 = 2*(phiRE00*phiRE0Pi + phiIM00*phiIM0Pi)*inorm*inorm*INITphase0Pi;
        sum5 = 2*(phiREPiPi*phiREPi0 + phiIMPiPi*phiIMPi0)*inorm*inorm*INITphasePiPi*INITphasePi0;
        sum6 = 2*(phiREPiPi*phiRE0Pi + phiIMPiPi*phiIM0Pi)*inorm*inorm*INITphasePiPi*INITphase0Pi;
        sum7 = 2*(phiREPi0*phiRE0Pi + phiIMPi0*phiIM0Pi)*inorm*inorm*INITphasePi0*INITphase0Pi;
        for(r=0;r<LX;r++){
           EyProf[r]  += (Ey00[r][k]*sum1  + Ey01[r][k]*sum2 + Ey02[r][k]*sum3 + Ey03[r][k]*sum4
                         + Ey12[r][k]*sum5 + Ey13[r][k]*sum6 + Ey23[r][k]*sum7);
        }
      } // close the loop over bag states
      /* print the Ey profile at each times */
      fprintf(fptr1,"%lf ",t);
      for(r=0;r<LX;r++){ fprintf(fptr1,"% lf ",EyProf[r]); }
      fprintf(fptr1,"\n"); 
  }
  fclose(fptr1); 
  /* free memory */
  EyProf.clear(); 
  alpha00.clear(); alphaPiPi.clear(); alphaPi0.clear(); alpha0Pi.clear();
  Ey00.clear(); Ey01.clear(); Ey02.clear(); Ey03.clear(); Ey12.clear(); Ey13.clear(); Ey23.clear();
  cos00.clear(); cosPiPi.clear(); cosPi0.clear(); cos0Pi.clear();
  sin00.clear(); sinPiPi.clear(); sinPi0.clear(); sin0Pi.clear();
}
