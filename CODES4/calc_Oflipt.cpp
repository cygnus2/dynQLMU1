/* calculate the time evolution of Oflip on an initial state */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: eigenstate |psi_n> = \sum_k \alpha_k |k>
void calc_Oflipt(){
  MKL_INT p,q,r,sizet,q1,k,m;
  MKL_INT sizesp;
  double t;
  sizet =  Wind.nBasis;
  double phiRE, phiIM, phiT, oflipt;
  double v_q, O_q;
  FILE *fptr,*fptr1;
  std::vector<double> initC;
  double fprof[VOL];

  // choosing different values of ch and cx allows choice of different
  // initial states with the same number of flippable plaquettes
  int ch, cx;

  /* initialize state in maximal flippable state */
  cx = 0; ch = 0; //choose the first state which has nflipMax plaquettes
  for(k=0; k<sizet; k++){
	    if( Wind.nflip[k] == Wind.nflipMax ){
         if(cx == ch){ q = k; break; }
		     cx++;
	    }
  }
  std::cout<<"Starting state is basis state = "<<q<<std::endl;
  /* store the overlap of the initial state with the eigenvectors */
  for(p=0; p<sizet; p++){
   initC.push_back(Wind.evecs[p*sizet+q]);
  }

  // open file to write
  fptr = fopen("OflipT.dat","w");
  fptr1= fopen("Fprof.dat","w");

  // the time-evolution
  for(t=Ti;t<Tf;t=t+dT){

     /* initialize flux profile and total flux */
     for(k=0;k<VOL;k++) fprof[k]=0.0;
     oflipt = 0.0;

     for(k=0; k<sizet; k++){
       phiRE = 0.0; phiIM = 0.0;
       for(m=0; m<sizet; m++){
         phiRE += initC[m]*Wind.evecs[sizet*m+k]*cos(Wind.evals[m]*t);
         phiIM += initC[m]*Wind.evecs[sizet*m+k]*sin(Wind.evals[m]*t);
       }
       phiT = phiRE*phiRE + phiIM*phiIM;
       oflipt += phiT*Wind.nflip[k];
       /* get the flippability profile at time t */
       for(r=0;r<VOL;r++){
         if(Wind.xflip[k][r]) fprof[r] += phiT;
       }
     } // close the loop over basis states

     oflipt = oflipt/((double)VOL);
     fprintf(fptr,"%.4lf %.12lf \n",t,oflipt);

     // print the profile
     fprintf(fptr1,"%.4lf ",t);
     for(r=0; r<VOL; r++) fprintf(fptr1, "  %.6lf  ", fprof[r]);
     fprintf(fptr1,"\n");

 }
 fclose(fptr);
 fclose(fptr1);

 /* clear memory */
 initC.clear();

}
