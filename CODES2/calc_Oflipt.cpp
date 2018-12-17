/* calculate the time evolution of Oflip on an initial state */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

extern void cartoonState(int, int, std::vector<bool>& );

// Notation: eigenstate |psi_n> = \sum_k \alpha_k |k>, |k> is a basis state in 
//           specified winding number (wx,wy) sector.
void calc_Oflipt(int sector, int wx, int wy){
  MKL_INT p,q,sizet,q1,k,m;
  double t;
  sizet =  Wind[sector].nBasis;
  std::vector<bool> cart1(2*VOL);
  std::vector<double> alpha;
  double phiRE, phiIM, oflipt; 
  FILE *outf;

  /* construct cartoon state */
  cartoonState(wx, wy, cart1);
  q1=Wind[sector].binscan(cart1);
  for(p=0; p<sizet; p++){
     alpha.push_back(Wind[sector].evecs[p*sizet+q1]);
  }
  outf = fopen("Evol_Oflip.dat","w");

  for(t=Ti;t<Tf;t=t+dT){
     oflipt = 0.0;
     for(k=0; k<sizet; k++){
       phiRE = 0.0; phiIM = 0.0;
       for(m=0; m<sizet; m++){
         phiRE += alpha[m]*Wind[sector].evecs[sizet*m+k]*cos(Wind[sector].evals[m]*t);
         phiIM += alpha[m]*Wind[sector].evecs[sizet*m+k]*sin(Wind[sector].evals[m]*t);
       }  
       oflipt += (phiRE*phiRE +  phiIM*phiIM)*Wind[sector].nflip[k];
     }
     fprintf(outf,"%.4lf %.12lf\n",t,oflipt);
  }
  fclose(outf);
  /* clear memory */
  cart1.clear();
  alpha.clear();



}
