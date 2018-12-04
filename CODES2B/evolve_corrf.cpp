/* real-time correlation function */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

void evolve_corrf1(int sector){
   int ix,iy,parity,p,q1,q2,r;
   int sizet;
   double t;
   std::vector<bool> cart1(2*VOL), cart2(2*VOL);
   std::vector<double> alpha1,alpha2;
   FILE *outf1;

   std::vector<double> state1RE(LX/2, 0.0), state1IM(LX/2, 0.0);
   std::vector<double> temp1(LX/2, 0.0);
   sizet = Wind[sector].nBasis;

   outf1 = fopen("corrf1.dat","w");
   for(p=0; p<sizet; p++){
    for(r=0; r<(LX/2); r++){
    fprintf(outf1,"% lf ",Wind[sector].oflip[p][r]); 
    }
    fprintf(outf1,"\n");
   }
   fclose(outf1);

   /* construct cartoon state */
   for(iy=0;iy<LY;iy++){
   for(ix=0;ix<LX;ix++){
    parity=(ix+iy)%2;
    p = 2*(iy*LX+ix);
    if(parity){
       cart1[p]=false; cart1[p+1]=true; cart2[p]=true;  cart2[p+1]=false;
    }
    else{
       cart1[p]=true; cart1[p+1]=false; cart2[p]=false;  cart2[p+1]=true;
    }
   }}

   q1=Wind[sector].binscan(cart1); 
   q2=Wind[sector].binscan(cart2);
   for(p=0; p<sizet; p++){
     alpha1.push_back(Wind[sector].evecs[p*sizet+q1]);
     alpha2.push_back(Wind[sector].evecs[p*sizet+q2]);
   }

   outf1 = fopen("CORRF1_STATE1.dat","w");
   // compute the real-time evolution
   for(p=0; p<sizet; p++){
      if(fabs(Wind[sector].evals[p]) < 1e-10){
      for(r=0;r<(LX/2);r++){
        temp1[r] = temp1[r] + alpha1[p]*alpha1[p]*Wind[sector].oflip[p][r];
      }
   }}

   for(t=Ti;t<Tf;t=t+dT){
     for(r=0;r<(LX/2);r++){
        state1RE[r]=0.0; state1IM[r]=0.0; 
     }
     for(p=0; p<sizet; p++){
        if(fabs(Wind[sector].evals[p]) > 1e-10){
          for(r=0;r<(LX/2);r++){
             state1RE[r] = state1RE[r] + alpha1[p]*alpha1[p]*Wind[sector].oflip[p][r]*cos(Wind[sector].evals[p]*t);
             state1IM[r] = state1IM[r] + alpha1[p]*alpha1[p]*Wind[sector].oflip[p][r]*sin(Wind[sector].evals[p]*t);
          }
        } // close if-loop
     } // close for-loop over eigenstates
     fprintf(outf1,"%.4f ",t);
     for(r=0;r<(LX/2);r++){
        state1RE[r]+=temp1[r]; 
        fprintf(outf1,"% lf % lf ",state1RE[r],state1IM[r]); 
     }
     fprintf(outf1,"\n");
   }
   fclose(outf1);  
   alpha1.clear();
   alpha2.clear();
}
