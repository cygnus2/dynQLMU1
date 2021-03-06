/* evolve states in real-time with the Hamiltonian */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void evolve_cartoons(std::vector<double> &evals, std::vector<std::vector<double>> &evecs){
   extern int scan(std::vector<bool>&);
   int ix,iy,parity,p,q1,q2;
   double t;
   double temp1,temp2;
   double ampl1,ampl2,lam1,lam2;
   std::vector<bool> cart1(2*VOL), cart2(2*VOL);
   std::vector<double> alpha1,alpha2;
   FILE *outf;
   outf = fopen("overlap.dat","w");
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
   for(p=0;p<NH;p++){
     q1=scan(cart1); 
     q2=scan(cart2);
     alpha1.push_back(evecs[p][q1]);
     alpha2.push_back(evecs[p][q2]);
   }
   // check that the coefficients are correctly set
   //for(p=0;p<NH;p++) std::cout<<alpha1[p]<<" "<<alpha2[p]<<std::endl;
   // compute the real-time evolution
   temp1 = 0.0; temp2 = 0.0;
   for(p=0;p<NH;p++){
    if(fabs(evals[p]) < 1e-10){
      temp1 = temp1 + alpha1[p]*alpha1[p];
      temp2 = temp2 + alpha2[p]*alpha1[p];
    }
   }
   //std::cout<<ampl1<<" "<<ampl2<<std::endl;
   for(t=Ti;t<Tf;t=t+dT){
     ampl1=0.0; ampl2=0.0;
     for(p=0;p<NH;p++){
        if(fabs(evals[p]) > 1e-10){
          ampl1 = ampl1 + alpha1[p]*alpha1[p]*cos(evals[p]*t);
          ampl2 = ampl2 + alpha2[p]*alpha1[p]*cos(evals[p]*t);
        } // close if-loop
     } // close for-loop
     ampl1 = ampl1 + temp1; ampl2 = ampl2 + temp2;
     //std::cout<<t<<" "<<ampl1<<std::endl;
     lam1 = -log(ampl1*ampl1)/VOL;
     lam2 = -log(ampl2*ampl2)/VOL;
     fprintf(outf,"%.2f % lf % lf % lf % lf\n",t,ampl1,ampl2,lam1,lam2);
   }
   fclose(outf);
}
