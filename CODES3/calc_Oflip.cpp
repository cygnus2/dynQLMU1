/* calculate Oflip in the chosen winding sector for the different momenta  */
/* Oflip = (1/VOL) \sum_{xy, plaq} (1=flippable plaq, 0 otherwise)         */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Note that this routine is already included in "calc_oflipt.cpp"
// There is no point in running this code here.

// Notation: The eigenvectors in of the translation matrix is denoted as |w_k>.
// |w_k> = \sum_l B_l |b_l>, where b_l are the "bag states", ie, states which
// are labelled with the translation flag. We calculate:
// <w_k| O_flip |w_k> = \sum_l \sum_m B_l B_m <b_l| O_flip |b_m>
// |b_l> = \sum_alpha c_alpha |alpha>; alpha are the cartoon-states (in the given W-no sector)
// <b_l| O_flip |b_m> = \sum_{beta,alpha} c_alpha c_beta <beta_l| O_flip |alpha_m>
//                    = \sum_{alpha} |c_alpha|^2 <alpha_l| O_flip |alpha_l> delta_{l,m}
// Hence, we finall have, for each k,
// <w_k| O_flip |w_k> = \sum_l |B_l|^2 <b_l| O_flip |b_l>

void calc_Oflip(int sector){
  int k,l,p;
  unsigned int sizet = Wind[sector].trans_sectors;
  std::vector<double> oflip(sizet);
  double Oflip_avg;
  FILE *outf;

  // initialize
  for(k=0;k<sizet;k++) oflip.push_back(0.0);
  // Calculate: < b_l | O_flip | b_l >; The expression below is correct for the momenta:
  // (pi,0), (0,pi) and (pi,pi); since the phase factors cancel out.
  for(p=0;p<Wind[sector].nBasis;p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }

  outf = fopen("Oflip.dat","w");
  fprintf(outf,"# Results for (kx,ky)=(0,0) \n");
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs_K00[k*sizet+l]* Wind[sector].evecs_K00[k*sizet+l]*oflip[l];
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%.12lf %.12lf\n",Wind[sector].evals_K00[k],Oflip_avg);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(Pi,Pi) \n");
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs_KPiPi[k*sizet+l]* Wind[sector].evecs_KPiPi[k*sizet+l]*oflip[l];
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%.12lf %.12lf\n",Wind[sector].evals_KPiPi[k],Oflip_avg);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(Pi,0) \n");
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs_KPi0[k*sizet+l]* Wind[sector].evecs_KPi0[k*sizet+l]*oflip[l];
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%.12lf %.12lf\n",Wind[sector].evals_KPi0[k],Oflip_avg);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(0,Pi) \n");
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs_K0Pi[k*sizet+l]* Wind[sector].evecs_K0Pi[k*sizet+l]*oflip[l];
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%.12lf %.12lf\n",Wind[sector].evals_K0Pi[k],Oflip_avg);
  }
 fclose(outf);
}
