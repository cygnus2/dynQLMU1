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
  unsigned int sizet;
  double Oflip_avg;
  FILE *outf;
  outf = fopen("Oflip.dat","w");
  fprintf(outf,"# Results of winding number sector (%d,%d) \n",Wind[sector].Wx,Wind[sector].Wy);

  // scan through all the eigenvalues
  sizet = Wind[sector].nBasis;  
  for(p=0;p<sizet;p++){
    //std::cout<<" state = "<<p<<" flippable plaq ="<<Wind[sector].nflip[p]<<std::endl;
    // calculate the expectation value in each eigenstate
    Oflip_avg = 0.0;
    for(q=0;q<sizet;q++){
      Oflip_avg += Wind[sector].evecs[p*sizet+q]*Wind[sector].evecs[p*sizet+q]*Wind[sector].nflip[q];
    }
    Oflip_avg = Oflip_avg/((double)VOL);
    fprintf(outf,"%lf %lf\n",Wind[sector].evals[p],Oflip_avg);
  }
 fclose(outf);
}
