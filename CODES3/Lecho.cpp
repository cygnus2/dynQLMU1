/*  Evolve states in real-time with the Hamiltonian. This routine computes the *
 *  overlap probability to the initial state, and the Loschmidt Echo. All the  *
 *  initial states can be handled with this routine.                           *
 */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void Lecho_INIT0(int sector){
  MKL_INT p,q,r,sizet;
  std::size_t k,m;
  int q1,tsect;
  double t;
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi;
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double fidelity;
  FILE *outf;

  /* initialize the starting state */
  if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
  q1 = INITq;
  std::cout<<"In function Lecho. Starting state is basis state = "<<q1<<std::endl;

  for(p=0; p<tsect; p++){
    alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]*inorm);
    alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]*inorm);
  }

  outf = fopen("overlap.dat","w");
  for(t=Ti;t<Tf;t=t+dT){
    phiRE00=0.0;    phiIM00=0.0;
    phiREPiPi=0.0;  phiIMPiPi=0.0;
    for(k=0; k<tsect; k++){
      /* contributions from sector (0,0) */
      phiRE00 += alpha00[k]*alpha00[k]*cos(Wind[sector].evals_K00[k]*t);
      phiIM00 += alpha00[k]*alpha00[k]*sin(Wind[sector].evals_K00[k]*t);
      /* contributions from sector (pi,pi) */
      phiREPiPi += alphaPiPi[k]*alphaPiPi[k]*cos(Wind[sector].evals_KPiPi[k]*t);
      phiIMPiPi += alphaPiPi[k]*alphaPiPi[k]*sin(Wind[sector].evals_KPiPi[k]*t);
      }
    fidelity = (phiRE00 + phiREPiPi)*(phiRE00 + phiREPiPi) + (phiIM00 + phiIMPiPi)*(phiIM00 + phiIMPiPi);
    fprintf(outf,"%lf %.12lf\n",t,fidelity);
  }
  fclose(outf);
  // deallocate
  alpha00.clear();  alphaPiPi.clear();
}

void Lecho_INIT4(int sector){
  MKL_INT p,q,r,sizet;
  std::size_t k,m;
  int q1,tsect;
  double t;
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi,alphaPi0,alpha0Pi;
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double phiREPi0, phiIMPi0, phiRE0Pi, phiIM0Pi;
  double fidelity;
  FILE *outf;

  /* initialize the starting state */
  if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
  q1 = INITq;
  std::cout<<"In function Lecho. Starting state is basis state = "<<q1<<std::endl;

  for(p=0; p<tsect; p++){
    alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]*inorm);
    alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]*inorm);
    alphaPi0.push_back(Wind[sector].evecs_KPi0[p*tsect+INITbag]*inorm);
    alpha0Pi.push_back(Wind[sector].evecs_K0Pi[p*tsect+INITbag]*inorm);
  }

  outf = fopen("overlap.dat","w");
  for(t=Ti;t<Tf;t=t+dT){
    phiRE00=0.0;    phiIM00=0.0;  phiREPiPi=0.0;  phiIMPiPi=0.0;
    phiREPi0=0.0;   phiIMPi0=0.0; phiRE0Pi=0.0;   phiIM0Pi=0.0;
    for(k=0; k<tsect; k++){
      /* contributions from sector (0,0) */
      phiRE00 += alpha00[k]*alpha00[k]*cos(Wind[sector].evals_K00[k]*t);
      phiIM00 += alpha00[k]*alpha00[k]*sin(Wind[sector].evals_K00[k]*t);
      /* contributions from sector (pi,pi) */
      phiREPiPi += alphaPiPi[k]*alphaPiPi[k]*cos(Wind[sector].evals_KPiPi[k]*t);
      phiIMPiPi += alphaPiPi[k]*alphaPiPi[k]*sin(Wind[sector].evals_KPiPi[k]*t);
      /* contributions from sector (pi,0) */
      phiREPi0 += alphaPi0[k]*alphaPi0[k]*cos(Wind[sector].evals_KPi0[k]*t);
      phiIMPi0 += alphaPi0[k]*alphaPi0[k]*sin(Wind[sector].evals_KPi0[k]*t);
      /* contributions from sector (0,pi) */
      phiRE0Pi += alpha0Pi[k]*alpha0Pi[k]*cos(Wind[sector].evals_K0Pi[k]*t);
      phiIM0Pi += alpha0Pi[k]*alpha0Pi[k]*sin(Wind[sector].evals_K0Pi[k]*t);
    }
    fidelity = (phiRE00 + phiREPiPi + phiREPi0 + phiRE0Pi)*(phiRE00 + phiREPiPi + phiREPi0 + phiRE0Pi)
      + (phiIM00 + phiIMPiPi + phiIMPi0 + phiIM0Pi)*(phiIM00 + phiIMPiPi + phiIMPi0 + phiIM0Pi);
    //fprintf(outf,"%lf %.12lf %.12lf %.12lf %.12lf %.12lf\n",t,phiRE00,phiIM00,phiREPiPi,phiIMPiPi,fidelity);
    fprintf(outf,"%lf %.12lf\n",t,fidelity);
  }
  fclose(outf);
  // deallocate
  alpha00.clear();  alphaPiPi.clear(); alphaPi0.clear(); alpha0Pi.clear();
}
