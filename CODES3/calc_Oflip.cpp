/* calculate the local operator Oflip in a given winding number sector */
/* Oflip calculates the expectation value of the number of flippable
   plaquettes in a given eigenstate, normalized with the total number 
   of plaquettes   */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: The eigenvectors in of the translation matrix is denoted as |w_k>. 
// |w_k> = \sum_l B_l |b_l>, where b_l are the "bag states", ie, states which
// are labelled with the translation flag. We calculate:
// <w_k| O_flip |w_k> = \sum_l \sum_m B_l B_m <b_l| O_flip |b_m>
// |b_l> = \sum_alpha c_alpha |alpha>; alpha are the ice-states (in the given W-no sector)
// <b_l| O_flip |b_m> = \sum_{l,m} c_alpha c_beta <beta_l| O_flip |alpha_m> 
//                    = \sum_l |c_alpha|^2 <alpha_l| O_flip |alpha_l> delta_{l,m}
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

  // calculate the matrix element: < b_l | O_flip | b_l >
  for(p=0;p<Wind[sector].nBasis;p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }

  outf = fopen("Oflip.dat","w");
  fprintf(outf,"# Results in (Wx,Wy)=(%d,%d), and (kx,ky)=(0,0) \n",Wind[sector].Wx,Wind[sector].Wy);
  for(k=0;k<sizet;k++){
    // calculate the expectation value in each eigenstate in translation basis
    Oflip_avg = 0.0;
    for(l=0;l<sizet;l++){
     Oflip_avg += Wind[sector].evecs[k*sizet+l]* Wind[sector].evecs[k*sizet+l]*oflip[l]; 
    }
    Oflip_avg /= ((double)VOL);
    fprintf(outf,"%lf %lf\n",Wind[sector].evals[k],Oflip_avg);
  }
 fclose(outf);
}
