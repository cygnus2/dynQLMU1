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

// Notation: eigenstate |psi_n> = \sum_k \alpha_k |k>, |k> is a basis state in
//           specified winding number (wx,wy) sector.
void calc_Oflip(int sector){
  MKL_INT p,q,sizet;
  double v_q, O_q, E_q;
  double Oflip_avg, CorrEy1_avg;
  FILE *outf;

  sizet = Wind[sector].nBasis;
  outf = fopen("LocalOps.dat","w");
  fprintf(outf,"# Results of winding number sector (%d,%d) \n",Wind[sector].Wx,Wind[sector].Wy);
  // scan through all the eigenvalues
  for(p=0;p<sizet;p++){
    // calculate the expectation value in each eigenstate
    Oflip_avg = 0.0; CorrEy1_avg = 0.0;
    for(q=0;q<sizet;q++){
      v_q = Wind[sector].evecs[p*sizet+q];
      O_q = Wind[sector].nflip[q];
      E_q = Wind[sector].CEy1[q];
      Oflip_avg   += v_q*v_q*O_q;
      CorrEy1_avg += v_q*v_q*E_q;
    }
    Oflip_avg   = Oflip_avg/((double)VOL);
    fprintf(outf,"%.12lf %.12lf %.12lf\n",Wind[sector].evals[p],Oflip_avg, CorrEy1_avg);
  }
  fclose(outf);
}
