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
void calc_Oflip(int sector){
  int k,l,p,kk,ll;
  int chkindex;
  unsigned int sizet = Wind[sector].trans_sectors;
  std::vector<double> oflip(sizet);
  std::vector<double> cEy(sizet);
  double Oflip_avg, cEy_avg;
  FILE *outf;
  double amp, shannonE, IPR;

  // initialize
  for(k=0;k<sizet;k++) { oflip.push_back(0.0); cEy.push_back(0.0); }

  // Calculate: < b_l | O_flip | b_l >; The expression below is correct for the momenta:
  // (pi,0), (0,pi) and (pi,pi); since the phase factors cancel out.
  for(p=0;p<Wind[sector].nBasis;p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
     cEy[Wind[sector].Tflag[p]-1]   += Wind[sector].corrE1[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }

  outf = fopen("Oflip.dat","w");
  fprintf(outf,"# Results for (kx,ky)=(0,0) \n");
  for(k=0;k<sizet;k++){
    shannonE=0.0; IPR=0.0;
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0; cEy_avg = 0.0;
    for(l=0;l<sizet;l++){
     amp        = Wind[sector].evecs_K00[k*sizet+l];
     Oflip_avg += amp*amp*oflip[l];
     cEy_avg   += amp*amp*cEy[l];
     if(fabs(amp) > 1e-10) shannonE  -= (amp*amp)*log(amp*amp);
     IPR       += (amp*amp*amp*amp);
    }
    Oflip_avg /= ((double)VOL);   cEy_avg /= ((double)VOL);
    fprintf(outf,"%.12lf %.12lf %.12lf %.12lf %.12lf\n",Wind[sector].evals_K00[k],Oflip_avg,cEy_avg,shannonE,IPR);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(Pi,Pi) \n");
  for(k=0;k<sizet;k++){
    shannonE=0.0; IPR=0.0;
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0; cEy_avg = 0.0;
    for(l=0;l<sizet;l++){
     amp        = Wind[sector].evecs_KPiPi[k*sizet+l];
     Oflip_avg += amp*amp*oflip[l];
     cEy_avg   += amp*amp*cEy[l];
     if(fabs(amp) > 1e-10) shannonE  -= (amp*amp)*log(amp*amp);
     IPR       += (amp*amp*amp*amp);
    }
    Oflip_avg /= ((double)VOL);  cEy_avg /= ((double)VOL);
    fprintf(outf,"%.12lf %.12lf %.12lf %.12lf %.12lf\n",Wind[sector].evals_KPiPi[k],Oflip_avg,cEy_avg,shannonE,IPR);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(Pi,0) \n");
  // this needs to look for correct labels to the translation basis
  for(k=0;k<sizet;k++){
    shannonE=0.0; IPR=0.0;
    Oflip_avg = 0.0; cEy_avg = 0.0;
    // if this is not a physical bag then skip
    kk = labelPi0[k]; if(kk==-997) continue;
    for(l=0;l<sizet;l++){
      ll = labelPi0[l]; if(ll==-997) continue;
      amp        = Wind[sector].evecs_KPi0[kk*nDimPi0+ll];
      Oflip_avg += amp*amp*oflip[l];
      cEy_avg   += amp*amp*cEy[l];
      if(fabs(amp) > 1e-10) shannonE  -= (amp*amp)*log(amp*amp);
      IPR       += (amp*amp*amp*amp);
    }
    Oflip_avg /= ((double)VOL); cEy_avg /= ((double)VOL);
    fprintf(outf,"%.12lf %.12lf %.12lf %.12lf %.12lf\n",Wind[sector].evals_KPi0[kk],Oflip_avg,cEy_avg,shannonE,IPR);
  }
  fprintf(outf,"\n\n# Results for (kx,ky)=(0,Pi) \n");
  // this needs to look for correct labels to the translation basis
  for(k=0;k<sizet;k++){
    shannonE=0.0; IPR=0.0;
    Oflip_avg = 0.0; cEy_avg = 0.0;
    // if this is not a physical bag then skip
    kk = label0Pi[k]; if(kk==-997) continue;
    for(l=0;l<sizet;l++){
      ll = label0Pi[l]; if(ll==-997) continue;
      amp        = Wind[sector].evecs_K0Pi[kk*nDim0Pi+ll];
      Oflip_avg += amp*amp*oflip[l];
      cEy_avg   += amp*amp*cEy[l];
      if(fabs(amp) > 1e-10) shannonE  -= (amp*amp)*log(amp*amp);
      IPR       += (amp*amp*amp*amp);
    }
    Oflip_avg /= ((double)VOL); cEy_avg /= ((double)VOL);
    fprintf(outf,"%.12lf %.12lf %.12lf %.12lf %.12lf\n",Wind[sector].evals_K0Pi[kk],Oflip_avg,cEy_avg,shannonE,IPR);
  }
 fclose(outf);
}
