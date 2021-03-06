/*  Evolve states in real-time with the Hamiltonian. This routine computes the *
 *  overlap probability to the initial state, and the Loschmidt Echo. All the  *
 *  initial states can be handled with this routine.                           *
 */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void Lecho(int sector){
   int ix,iy,parity,p,q,q1,q2;
   int sizet;
   double t;
   double temp;
   double ampl_RE,ampl_IM;
   double lEcho,fidelity;
   std::vector<double> initC;
   FILE *outf;

   sizet = Wind[sector].nBasis;
   /* obtain the initial starting state */
   if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
   q = INITq;
   std::cout<<"Starting state is basis state = "<<q<<std::endl;
   /* store the overlap of the initial state with the eigenvectors */
   for(p=0; p<sizet; p++){
       initC.push_back(Wind[sector].evecs[p*sizet+q]);
   }

   outf = fopen("overlap.dat","w");
   for(t=Ti;t<Tf;t=t+dT){
     ampl_RE=0.0; ampl_IM=0.0;
     for(p=0; p<sizet; p++){
          ampl_RE = ampl_RE + initC[p]*initC[p]*cos(Wind[sector].evals[p]*t);
          ampl_IM = ampl_IM + initC[p]*initC[p]*sin(Wind[sector].evals[p]*t);
     } // close for-loop
     fidelity = ampl_RE*ampl_RE + ampl_IM*ampl_IM;
     fprintf(outf,"%lf %.12lf %.12lf %.12lf\n",t,ampl_RE,ampl_IM,fidelity);
   }
   fclose(outf);

   initC.clear();
}
