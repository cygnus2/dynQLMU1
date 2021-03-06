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
  MKL_INT p, q, sizet,nB;
  double v_q, v_l, O_q, O_l; 
  double Oflip_avg;
  FILE *outf;
  outf = fopen("Oflip.dat","w");
  //fprintf(outf,"# Results of winding number sector (%d,%d) in C=+ sector\n",Wind[sector].Wx,Wind[sector].Wy);

  // scan through all the eigenvalues
  nB    = Wind[sector].nBasis;
  sizet = nB/2;  
  for(p=0;p<sizet;p++){
    //std::cout<<" state = "<<p<<" flippable plaq ="<<Wind[sector].nflip[p]<<std::endl;
    // calculate the expectation value in each eigenstate
    Oflip_avg = 0.0;
    for(q=0;q<sizet;q++){
      v_q = Wind[sector].evecs[p*sizet+q]; 
      O_q = Wind[sector].nflip[q]; 
      O_l = Wind[sector].nflip[nB-1-q];
      Oflip_avg += 0.5*v_q*v_q*(O_q + O_l);
    }
    Oflip_avg = Oflip_avg/((double)VOL);
    fprintf(outf,"%.12lf %.12lf\n",Wind[sector].evals[p],Oflip_avg);
  }
 fclose(outf);
}

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a basis state in 
//           specified winding number (wx,wy) sector.
// this routine calculates the Oflip for the ED done in the whole (Wx,Wy)=(0,0) 
// sector without using the CC symmetry. Useful to see what to expect when a
// discrete symmetry is not taken into account
void calc_Oflip_all(int sector){
  unsigned int p, q, sizet;
  double Oflip_avg;
  FILE *outf;
  outf = fopen("Oflip_full.dat","w");
  fprintf(outf,"# Results of winding number sector (%d,%d)\n",Wind[sector].Wx,Wind[sector].Wy);

  // scan through all the eigenvalues
  sizet  = Wind[sector].nBasis;
  for(p=0;p<sizet;p++){
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
