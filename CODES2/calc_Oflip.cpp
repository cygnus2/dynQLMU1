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

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a basis state in 
//           specified winding number (wx,wy) sector.
void calc_Oflip(int sector){
  int p,q;
  double Oflip_avg;
  FILE *outf;
  outf = fopen("Oflip.dat","w");
  fprintf(outf,"# Results of winding number sector (%d,%d) \n",Wind[sector].Wx,Wind[sector].Wy);
  // scan through all the eigenvalues
  for(p=0;p<Wind[sector].nBasis;p++){
    // calculate the expectation value in each eigenstate
    Oflip_avg = 0.0;
    for(q=0;q<Wind[sector].nBasis;q++){
      Oflip_avg += Wind[sector].evecs[p][q]*Wind[sector].evecs[p][q]*Wind[sector].nflip[q];
    }
    Oflip_avg = Oflip_avg/((double)VOL);
    fprintf(outf,"%lf %lf\n",Wind[sector].evals[p],Oflip_avg);
  }
 fclose(outf);
}
