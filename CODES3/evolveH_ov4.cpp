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
// INIT=0 ==> compute the Ey and the oflip correlators
void evolveH_ov4_INIT0(int sector){
  MKL_INT p,q,q1,r,k,m;
  int sizet,tsect;
  double t;
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi;
  std::vector<double> cos00(tsect), cosPiPi(tsect), sin00(tsect), sinPiPi(tsect);
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double sum1, sum2, norm;
  double expCf0, expCf1;
  std::vector<double> Cf0sec00(tsect, 0.0); //overlap with identical sectors
  std::vector<double> Cf0sec01(tsect, 0.0); //mixing between two sectors
  std::vector<double> Cf1sec00(tsect, 0.0);
  std::vector<double> Cf1sec01(tsect, 0.0);
  //std::vector<double> avgd1, avgd2, avgh1, avgh2, avgv1, avgv2;
  FILE *fptr1;

  /* initialize the starting state */
  q1 = INITq;
  std::cout<<"In function evolveH_ov4, computing correlators. Starting state is basis state = "<<q1<<std::endl;

  // calculate Ey correlation functions
  for(p=0; p<sizet; p++){
     q=Wind[sector].Tflag[p]-1;
     norm=Wind[sector].Tdgen[p]/((double)VOL);
     Cf0sec00[q] += Wind[sector].CEy0[p]*norm;
     Cf0sec01[q] += Wind[sector].CEy0[p]*Wind[sector].FPiPi[p]*norm;
     Cf1sec00[q] += Wind[sector].CEy1[p]*norm;
     Cf1sec01[q] += Wind[sector].CEy1[p]*Wind[sector].FPiPi[p]*norm;
  }

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]);
  }

  // compute the state profile
  fptr1 = fopen("CorrF.dat","w");
  for(t=Ti; t<Tf; t=t+dT){
      // initialize the variables
      expCf0=0.0;  expCf1=0.0;
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
        // sum1 connects the same sectors, while sum2 has the mixing
        expCf0 = Cf0sec00[k]*sum1 + Cf0sec01[k]*sum2;
        expCf1 = Cf1sec00[k]*sum1 + Cf1sec01[k]*sum2;
      } // close the loop over bag states
      /* print the correlation functions */
      fprintf(fptr1,"%lf % .12lf % .12lf \n",t,expCf0,expCf1);
  }
  fclose(fptr1);
  /* free memory */
  alpha00.clear(); alphaPiPi.clear();
  Cf0sec00.clear(); Cf0sec01.clear(); Cf1sec00.clear(); Cf1sec01.clear();
  cos00.clear(); sin00.clear(); cosPiPi.clear(); sinPiPi.clear();
}

// INIT=4 ==> calculate correlators
void evolveH_ov4_INIT4(int sector){
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
  double expCf0, expCf1;
  std::vector<double> Cf0sec00(tsect, 0.0);  std::vector<double> Cf1sec00(tsect, 0.0);
  std::vector<double> Cf0sec01(tsect, 0.0);  std::vector<double> Cf1sec01(tsect, 0.0);
  std::vector<double> Cf0sec02(tsect, 0.0);  std::vector<double> Cf1sec02(tsect, 0.0);
  std::vector<double> Cf0sec03(tsect, 0.0);  std::vector<double> Cf1sec03(tsect, 0.0);
  std::vector<double> Cf0sec12(tsect, 0.0);  std::vector<double> Cf1sec12(tsect, 0.0);
  std::vector<double> Cf0sec13(tsect, 0.0);  std::vector<double> Cf1sec13(tsect, 0.0);
  std::vector<double> Cf0sec23(tsect, 0.0);  std::vector<double> Cf1sec23(tsect, 0.0);

  FILE *fptr1;

  /* initialize the starting state */
  q1 = INITq;
  std::cout<<"In function evolveH_ov4, computing correlators. Starting state is basis state = "<<q1<<std::endl;

  // calculate matrix elements of correlators in the bag basis.
  for(p=0; p<sizet; p++){
     q=Wind[sector].Tflag[p]-1;
     norm=Wind[sector].Tdgen[p]/((double)VOL);
     Cf0sec00[q] += Wind[sector].CEy0[p]*norm;
     Cf0sec01[q] += Wind[sector].CEy0[p]*Wind[sector].FPiPi[p]*norm;
     Cf0sec02[q] += Wind[sector].CEy0[p]*Wind[sector].FPi0[p]*norm;
     Cf0sec03[q] += Wind[sector].CEy0[p]*Wind[sector].F0Pi[p]*norm;
     Cf0sec12[q] += Wind[sector].CEy0[p]*Wind[sector].FPiPi[p]*Wind[sector].FPi0[p]*norm;
     Cf0sec13[q] += Wind[sector].CEy0[p]*Wind[sector].FPiPi[p]*Wind[sector].F0Pi[p]*norm;
     Cf0sec23[q] += Wind[sector].CEy0[p]*Wind[sector].FPi0[p]*Wind[sector].F0Pi[p]*norm;
     Cf1sec00[q] += Wind[sector].CEy1[p]*norm;
     Cf1sec01[q] += Wind[sector].CEy1[p]*Wind[sector].FPiPi[p]*norm;
     Cf1sec02[q] += Wind[sector].CEy1[p]*Wind[sector].FPi0[p]*norm;
     Cf1sec03[q] += Wind[sector].CEy1[p]*Wind[sector].F0Pi[p]*norm;
     Cf1sec12[q] += Wind[sector].CEy1[p]*Wind[sector].FPiPi[p]*Wind[sector].FPi0[p]*norm;
     Cf1sec13[q] += Wind[sector].CEy1[p]*Wind[sector].FPiPi[p]*Wind[sector].F0Pi[p]*norm;
     Cf1sec23[q] += Wind[sector].CEy1[p]*Wind[sector].FPi0[p]*Wind[sector].F0Pi[p]*norm;
  }

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]);
       alphaPi0.push_back(Wind[sector].evecs_KPi0[p*tsect+INITbag]);
       alpha0Pi.push_back(Wind[sector].evecs_K0Pi[p*tsect+INITbag]);
  }

  // compute the state profile
  fptr1 = fopen("CorrF.dat","w");
  for(t=Ti; t<Tf; t=t+dT){
      // initialize the expectation values
      expCf0=0.0; expCf1=0.0;
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
        // expectation value
        expCf0 += Cf0sec00[k]*sum1 + Cf0sec01[k]*sum2 + Cf0sec02[k]*sum3 + Cf0sec03[k]*sum4
                 + Cf0sec12[k]*sum5 + Cf0sec13[k]*sum6 + Cf0sec23[k]*sum7;
        expCf1 += Cf1sec00[k]*sum1 + Cf1sec01[k]*sum2 + Cf1sec02[k]*sum3 + Cf1sec03[k]*sum4
                 + Cf1sec12[k]*sum5 + Cf1sec13[k]*sum6 + Cf1sec23[k]*sum7;
      } // close the loop over bag states
      /* print the correlators at each times */
      fprintf(fptr1,"%lf % .12lf % .12lf \n",t, expCf0, expCf1);
  }
  fclose(fptr1);
  /* free memory */
  alpha00.clear(); alphaPiPi.clear(); alphaPi0.clear(); alpha0Pi.clear();
  cos00.clear(); cosPiPi.clear(); cosPi0.clear(); cos0Pi.clear();
  sin00.clear(); sinPiPi.clear(); sinPi0.clear(); sin0Pi.clear();
  Cf0sec00.clear(); Cf0sec01.clear(); Cf0sec02.clear(); Cf0sec03.clear();
  Cf0sec12.clear(); Cf0sec13.clear(); Cf0sec23.clear();
  Cf1sec00.clear(); Cf1sec01.clear(); Cf1sec02.clear(); Cf1sec03.clear();
  Cf1sec12.clear(); Cf1sec13.clear(); Cf1sec23.clear();
}
