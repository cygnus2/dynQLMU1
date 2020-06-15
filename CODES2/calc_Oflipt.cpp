/* calculate the time evolution of Oflip on an initial state */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

//extern void initState(int, int, int*);

// Notation: eigenstate |psi_n> = \sum_k \alpha_k |k>, |k> is a basis state in
//           specified winding number (wx,wy) sector.
void calc_Oflipt(int sector){
  MKL_INT p,q,r,sizet,k,m;
  MKL_INT sizesp;
  int q1;
  double t;
  sizet =  Wind[sector].nBasis;
  std::vector<bool> cart1(2*VOL);
  std::vector<double> alpha;
  //std::vector<double> CR(LX/2), CR_diag(LX/2);
  //std::vector<double> temp(LX/2);
  double phiRE, phiIM, oflipt;
  double Oflip_diag;
  double v_q, O_q;
  double Oflip_avg;
  FILE *outf;

  /* initialize the starting state */
  if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
  q1 = INITq;
  std::cout<<"In function calc_Oflipt. Starting state is basis state = "<<q1<<std::endl;
  for(p=0; p<sizet; p++){
     alpha.push_back(Wind[sector].evecs[p*sizet+q1]);
  }

  // Compute <PSI_M| O_flip|PSI_M> & <PSI_M| C_flip |PSI_M> for all PSI_M
  outf = fopen("Oflip.dat","w");
  fprintf(outf,"# Results of winding number sector (%d,%d) \n",Wind[sector].Wx,Wind[sector].Wy);
  Oflip_diag = 0.0;
  //for(r=0;r<(LX/2);r++) CR_diag[r]=0.0;
  // scan through all the eigenvalues
  for(p=0;p<sizet;p++){
    // calculate the expectation value in each eigenstate
    Oflip_avg = 0.0;
    //for(r=0;r<(LX/2);r++) temp[r]=0.0;
    for(q=0;q<sizet;q++){
      v_q = Wind[sector].evecs[p*sizet+q];
      O_q = Wind[sector].nflip[q];
      Oflip_avg += v_q*v_q*O_q;
      //for(r=0;r<(LX/2);r++) temp[r] += v_q*v_q*Wind[sector].cflip[q][r];
    }
    Oflip_diag += Oflip_avg*alpha[p]*alpha[p];
    Oflip_avg = Oflip_avg/((double)VOL);
    //for(r=0;r<(LX/2);r++) CR_diag[r] += temp[r]*alpha[p]*alpha[p];
    fprintf(outf,"%.12lf %.12lf\n",Wind[sector].evals[p],Oflip_avg);
  }
 fclose(outf);

 // diagonal ensemble results
 //for(r=0;r<(LX/2);r++) CR_diag[r] = CR_diag[r]/((double)VOL);
 Oflip_diag = Oflip_diag/((double)VOL);
 outf = fopen("OflipT.dat","w");
 fprintf(outf,"# value in the diagonal observable = %.12f ",Oflip_diag);
 //for(r=0;r<(LX/2);r++) fprintf(outf, " %lf ",CR_diag[r]);
 fprintf(outf,"\n");

 // the time-evolution
 for(t=Ti;t<Tf;t=t+dT){
     oflipt = 0.0;
     //for(r=0;r<(LX/2);r++) CR[r]=0.0;
     for(k=0; k<sizet; k++){
       phiRE = 0.0; phiIM = 0.0;
       for(m=0; m<sizet; m++){
         phiRE += alpha[m]*Wind[sector].evecs[sizet*m+k]*cos(-Wind[sector].evals[m]*t);
         phiIM += alpha[m]*Wind[sector].evecs[sizet*m+k]*sin(-Wind[sector].evals[m]*t);
       }
       oflipt += (phiRE*phiRE +  phiIM*phiIM)*Wind[sector].nflip[k];
       //for(r=0;r<(LX/2);r++) CR[r] += (phiRE*phiRE +  phiIM*phiIM)*Wind[sector].cflip[k][r];
     }
     oflipt = oflipt/((double)VOL);
     //for(r=0;r<(LX/2);r++) CR[r] = CR[r]/((double)VOL);
     fprintf(outf,"%.4lf %.12lf ",t,oflipt);
     //for(r=0;r<(LX/2);r++) fprintf(outf," %.8lf ",CR[r]);
     fprintf(outf,"\n");
 }
 fclose(outf);
 /* clear memory */
 alpha.clear();
 //CR.clear(); CR_diag.clear(); temp.clear();

}

// This routine serves as a check of the calculation of the real-time
// evolution done in the first routine. It takes much longer (as explained
// in the notes) and also checked in actual run-time.
void calc_Oflipt2(int sector){
  MKL_INT p,q,sizet,k,l,m,n;
  int q1;
  double t;
  sizet =  Wind[sector].nBasis;
  std::vector<double> alpha;
  double Omn,oflipt;
  FILE *outf;

  /* construct cartoon state */
  //initState(sector, INIT, &q1);
  q1 = INITq;
  std::cout<<"Starting state is basis state = "<<q1<<std::endl;
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
  alpha.clear();
}
