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
  double Oflip_diag;
  double v_q, O_q;
  double Oflip_avg;
  FILE *outf;
  
  /* construct cartoon state */
  cartoonState(wx, wy, cart1);
  q1=Wind[sector].binscan(cart1);
  for(p=0; p<sizet; p++){
     alpha.push_back(Wind[sector].evecs[p*sizet+q1]);
  }

  // Compute <PSI_M| O_flip|PSI_M> for every eigenstate PSI_M for the diag observable
  outf = fopen("Oflip.dat","w");
  fprintf(outf,"# Results of winding number sector (%d,%d) \n",Wind[sector].Wx,Wind[sector].Wy);
  Oflip_diag = 0.0;
  // scan through all the eigenvalues
  for(p=0;p<sizet;p++){
    // calculate the expectation value in each eigenstate
    Oflip_avg = 0.0;
    for(q=0;q<sizet;q++){
      v_q = Wind[sector].evecs[p*sizet+q];
      O_q = Wind[sector].nflip[q];
      Oflip_avg += v_q*v_q*O_q;
    }
    Oflip_diag += Oflip_avg*alpha[p]*alpha[p];
    //Oflip_avg = Oflip_avg/((double)VOL);
    fprintf(outf,"%.12lf %.12lf\n",Wind[sector].evals[p],Oflip_avg);
  }
 fclose(outf);

 // the time-evolution bit
 outf = fopen("Evol_Oflip.dat","w");
 fprintf(outf,"# value in the diagonal observable = %f \n",Oflip_diag);
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

// This routine serves as a check of the calculation of the real-time
// evolution done in the first routine. It takes much longer (as explained
// in the notes) and also checked in actual run-time.
void calc_Oflipt2(int sector, int wx, int wy){
  MKL_INT p,q,sizet,q1,k,l,m,n;
  double t;
  sizet =  Wind[sector].nBasis;
  std::vector<bool> cart1(2*VOL);
  std::vector<double> alpha;
  double Omn,oflipt; 
  FILE *outf;

  /* construct cartoon state */
  cartoonState(wx, wy, cart1);
  q1=Wind[sector].binscan(cart1);
  for(p=0; p<sizet; p++){
     alpha.push_back(Wind[sector].evecs[p*sizet+q1]);
  }
  outf = fopen("Evol_Oflip2.dat","w");

  for(t=Ti;t<Tf;t=t+dT){
     oflipt = 0.0;
     for(m=0; m<sizet; m++){
     for(n=0; n<sizet; n++){
       Omn = 0.0;
       for(k=0; k<sizet; k++){
         Omn += Wind[sector].evecs[sizet*m+k]*Wind[sector].evecs[sizet*n+k]*Wind[sector].nflip[k];
       }  
       oflipt += alpha[m]*alpha[n]*Omn*cos((Wind[sector].evals[m] - Wind[sector].evals[n])*t);
     } }
     fprintf(outf,"%.4lf %.12lf\n",t,oflipt);
  }
  fclose(outf);
  /* clear memory */
  cart1.clear();
  alpha.clear();
}

