/* calculate the time evolution of Oflip on an initial state */
/* we need information from different momentum sectors depending on the intial state */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include<chrono>
#include "define.h"

extern void calc_Oflip(int, std::vector<double>&);

// eigenstate: |w_k> = \sum_l B_l |b_l>, |b_l> bag state
void calc_Oflipt(int sector){
  MKL_INT p,q,r,sizet;
  std::size_t k,m;
  int q1,tsect;
  double t;
  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<double> alpha00, alphaPiPi, alphaPi0, alpha0Pi;
  std::vector<double> oflip(tsect, 0.0);
  std::vector<double> cos00(tsect);
  std::vector<double> cosPiPi(tsect); 
  std::vector<double> cosPi0(tsect); 
  std::vector<double> cos0Pi(tsect);
  double phiRE00, phiIM00, phiREPiPi, phiIMPiPi;
  double phiRE0Pi, phiIM0Pi, phiREPi0, phiIMPi0;
  double oflipt;
  FILE *outf;

  /* initialize the starting state */
  if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
  q1 = INITq;
  std::cout<<"In function calc_Oflipt. Starting state is basis state = "<<q1<<std::endl;

  // Of in the bag basis: true for all momenta sectors
  for(p=0; p<sizet; p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }

  // < w_k | IN >; k-th eigenvector; IN=initial state
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

 // Get starting timepoint
 auto start = std::chrono::high_resolution_clock::now();

 // compute the static expectation values
 // value in the Diagonal Ensemble NOT computed yet!
 //calc_Oflip(sector, oflip);
 // compute the time evolution
 outf = fopen("OflipT.dat","w");
 if(INIT==0){ // for initial condition = 0 with fully flippable plaquettes
  for(t=Ti;t<Tf;t=t+dT){
     oflipt = 0.0;
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
       oflipt += (phiRE00*phiRE00 +  phiIM00*phiIM00 + phiREPiPi*phiREPiPi + phiIMPiPi*phiIMPiPi)*oflip[k];
     }
     oflipt = oflipt/((double)VOL);
     fprintf(outf,"%.4lf %.12lf ",t,oflipt);
     fprintf(outf,"\n");
  } //close loop over time
} //close for the case with INIT=0
else if(INIT==4){
  for(t=Ti;t<Tf;t=t+dT){
     oflipt = 0.0;
     for(std::size_t ii = 0; ii < tsect; ii++){
        cos00[ii]   = cos(Wind[sector].evals_K00[ii]*t);
        cosPiPi[ii] = cos(Wind[sector].evals_KPiPi[ii]*t);
        cosPi0[ii]  = cos(Wind[sector].evals_KPi0[ii]*t);
        cos0Pi[ii]  = cos(Wind[sector].evals_K0Pi[ii]*t);
     }
     for(k=0; k<tsect; k++){
       phiRE00 = 0.0; phiIM00 = 0.0; phiREPiPi=0.0; phiIMPiPi=0.0;
       phiRE0Pi = 0.0; phiIM0Pi = 0.0; phiREPi0=0.0; phiIMPi0=0.0;
       // off diagonal terms
       for(m=0; m<tsect; m++){
         /* contributions from sector (0,0) */
         phiRE00 += 2.0*alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*cos00[m];
         phiIM00 += 2.0*alpha00[m]*Wind[sector].evecs_K00[tsect*m+k]*sqrt(1.0-cos00[m]*cos00[m]);
         /* contributions from sector (pi,pi) */
         phiREPiPi += 2.0*alphaPiPi[m]*INITphasePiPi*Wind[sector].evecs_KPiPi[tsect*m+k]*cosPiPi[m];
         phiIMPiPi += 2.0*alphaPiPi[m]*INITphasePiPi*Wind[sector].evecs_KPiPi[tsect*m+k]*sqrt(1.0-cosPiPi[m]*cosPiPi[m]);
         /* contributions from sector (pi,0) + (0,pi) */
         phiREPi0 += 2.0*alphaPi0[m]*INITphasePi0*Wind[sector].evecs_KPi0[tsect*m+k]*cosPi0[m];
         phiIMPi0 += 2.0*alphaPi0[m]*INITphasePi0*Wind[sector].evecs_KPi0[tsect*m+k]*sqrt(1.0-cosPi0[m]*cosPi0[m]);
         phiRE0Pi += 2.0*alpha0Pi[m]*INITphase0Pi*Wind[sector].evecs_K0Pi[tsect*m+k]*cos0Pi[m];
         phiIM0Pi += 2.0*alpha0Pi[m]*INITphase0Pi*Wind[sector].evecs_K0Pi[tsect*m+k]*sqrt(1.0-cos0Pi[m]*cos0Pi[m]);
       }
       // diagonal terms
       /* contributions from sector (0,0) */
       //phiRE00 += alpha00[k]*Wind[sector].evecs_K00[tsect*k+k]*cos00[k];
       //phiIM00 += alpha00[k]*Wind[sector].evecs_K00[tsect*k+k]*sqrt(1.0-cos00[k]*cos00[k]);
       /* contributions from sector (pi,pi) */
       //phiREPiPi += alphaPiPi[k]*INITphasePiPi*Wind[sector].evecs_KPiPi[tsect*k+k]*cosPiPi[k];
       //phiIMPiPi += alphaPiPi[k]*INITphasePiPi*Wind[sector].evecs_KPiPi[tsect*k+k]*sqrt(1.0-cosPiPi[k]*cosPiPi[k]);
       /* contributions from sector (pi,0) + (0,pi) */
       //phiREPi0 += alphaPi0[k]*INITphasePi0*Wind[sector].evecs_KPi0[tsect*k+k]*cosPi0[k];
       //phiIMPi0 += alphaPi0[k]*INITphasePi0*Wind[sector].evecs_KPi0[tsect*k+k]*sqrt(1.0-cosPi0[k]*cosPi0[k]);
       //phiRE0Pi += alpha0Pi[k]*INITphase0Pi*Wind[sector].evecs_K0Pi[tsect*k+k]*cos0Pi[k];
       //phiIM0Pi += alpha0Pi[k]*INITphase0Pi*Wind[sector].evecs_K0Pi[tsect*k+k]*sqrt(1.0-cos0Pi[k]*cos0Pi[k]);

       oflipt += (phiRE00*phiRE00 +  phiIM00*phiIM00 + phiREPiPi*phiREPiPi + phiIMPiPi*phiIMPiPi +
         phiREPi0*phiREPi0 + phiIMPi0*phiIMPi0 + phiRE0Pi*phiRE0Pi + phiIM0Pi*phiIM0Pi)*oflip[k];
     }
     oflipt = oflipt/((double)VOL);
     fprintf(outf,"%.4lf %.12lf ",t,oflipt);
     fprintf(outf,"\n");
  } //close loop over time
 }
 fclose(outf);

 // Get ending timepoint
 auto stop = std::chrono::high_resolution_clock::now();
 auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
 std::cout<<"Time taken for calc_Oflipt="<<duration.count()<< " secs"<<std::endl;


 /* clear memory */
 alpha00.clear(); alphaPiPi.clear(); alpha0Pi.clear(); alphaPi0.clear();
 oflip.clear();
}

/* calculate Oflip in the chosen winding sector for the different momenta  */
/* Oflip = (1/VOL) \sum_{xy, plaq} (1=flippable plaq, 0 otherwise)         */
// Notation: The eigenvectors in of the translation matrix is denoted as |w_k>.
// |w_k> = \sum_l B_l |b_l>, where b_l are the "bag states", ie, states which
// are labelled with the translation flag. We calculate:
// <w_k| O_flip |w_k> = \sum_l \sum_m B_l B_m <b_l| O_flip |b_m>
// |b_l> = \sum_alpha c_alpha |alpha>; alpha are the cartoon-states (in the given W-no sector)
// <b_l| O_flip |b_m> = \sum_{beta,alpha} c_alpha c_beta <beta_l| O_flip |alpha_m>
//                    = \sum_{alpha} |c_alpha|^2 <alpha_l| O_flip |alpha_l> delta_{l,m}
// Hence, we finall have, for each k,
// <w_k| O_flip |w_k> = \sum_l |B_l|^2 <b_l| O_flip |b_l>
void calc_Oflip(int sector, std::vector<double>& of){
  int k,l,p;
  unsigned int sizet = Wind[sector].trans_sectors;
  double Oflip_avg;
  FILE *outf;

  outf = fopen("Oflip.dat","w");
  fprintf(outf,"# Results for (kx,ky)=(0,0) \n");
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs_K00[k*sizet+l]*Wind[sector].evecs_K00[k*sizet+l]*of[l];
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%lf %lf\n",Wind[sector].evals_K00[k],Oflip_avg);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(Pi,Pi) \n");
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs_KPiPi[k*sizet+l]* Wind[sector].evecs_KPiPi[k*sizet+l]*of[l];
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%lf %lf\n",Wind[sector].evals_KPiPi[k],Oflip_avg);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(Pi,0) \n");
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs_KPi0[k*sizet+l]* Wind[sector].evecs_KPi0[k*sizet+l]*of[l];
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%lf %lf\n",Wind[sector].evals_KPi0[k],Oflip_avg);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(0,Pi) \n");
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs_K0Pi[k*sizet+l]* Wind[sector].evecs_K0Pi[k*sizet+l]*of[l];
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%lf %lf\n",Wind[sector].evals_K0Pi[k],Oflip_avg);
  }
 fclose(outf);
}
