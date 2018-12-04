/* construct the Hamiltonian matrix from the basis vectors */
// This constructs the Hamiltonian in a chosen winding number basis
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"
#include "mkl.h"
#include "mkl_solvers_ee.h"
#include "mkl_types.h"

void constH(int sector){

   //extern void eigcheck(std::vector<double>&, std::vector<std::vector<double>>&,int); 
   //extern void diag_LAPACK_RRR(MKL_INT, std::vector<double>&, std::vector<double>&, std::vector<double>&);
   //extern void diag_LAPACK(MKL_INT, std::vector<double>&, std::vector<double>&, std::vector<double>&);
   extern void printmatrix(std::vector<MKL_INT>&,std::vector<MKL_INT>&,std::vector<double>&);
   extern void actH(int, std::vector<double>&, std::vector<double>&);
   MKL_INT i,j,k;
   int nK = Wind[sector].nKry;
   std::vector<double> alpha(nK, 0.0);
   std::vector<double> beta(nK, 0.0);
   std::cout<<"Krylov space = "<<nK<<std::endl;
   // workspace variables to construct the tridiagonal matrix in Krylov space
   MKL_INT nB = Wind[sector].nBasis;
   // to store the coefficients
   std::vector<double> init(nB,0.0);
   std::vector<double> v0(nB,0.0);
   std::vector<double> v1(nB,0.0);
   std::vector<double> v2(nB,0.0);

   /* initialize the starting vector */
   init[2]=1.0;
   std::cout<<"Starting state"<<std::endl;
   for(i=0; i<nB; i++){ std::cout<<init[i]<<" "; }
   std::cout<<std::endl;

   v0 = init; 
   // act the Hamiltonian on the initial state
   actH(sector, v0, v1);
   std::cout<<"Vector after action with the Hamiltonian"<<std::endl;
   for(i=0; i<nB; i++){
     std::cout<<v1[i]<<" ";
   }
   std::cout<<std::endl;
   alpha[0] = std::inner_product(v0.begin(), v0.end(), v1.begin(), 0.0);
   //std::cout<<alpha[0]<<std::endl;
   for(i=0; i<nB; i++){
      v2[i]=v1[i]-alpha[0]*v0[i];
   }
   
   // start the Lanczos loop
   for(j=1; j<nK; j++){
      v0 = v1;
      beta[j] = std::inner_product(v2.begin(), v2.end(), v2.begin(), 0.0);
      if(beta[j] < 1e-12) std::cout<<"Norm too small"<<std::endl;
      else{  beta[j] = std::sqrt(beta[j]);
             for(k=0; k<nB; k++){
                v1[k] = v2[k]/beta[j];
             }
      }
      actH(sector, v1, v2);
      alpha[j] = std::inner_product(v2.begin(), v2.end(), v1.begin(), 0.0);
      for(k=0; k<nB; k++){
         v2[k] = v2[k] - alpha[j]*v1[k] - beta[j]*v0[k];
      }
      std::cout<<j<<"  "<<beta[j]<<" "<<alpha[j-1]<<" "<<std::endl;
   }

   
}

int WindNo::binscan(std::vector<bool> &newstate){
     unsigned int m;
     // binary search of the sorted array  
     std::vector<std::vector<bool>>::iterator it;
     it = std::lower_bound(basisVec.begin(),basisVec.end(),newstate);
     m  = std::distance(basisVec.begin(),it);
     if(it == basisVec.end()){
       std::cout<<"Element not found here! "<<std::endl;
       return -100;
     }
     return m;
}

void actH(int sector, std::vector<double> &v0, std::vector<double> &v1){
   int i,p,q;
   int p1,p2,p3,p4;
   bool pxy,pyz,pzw,pwx;
   // to store the individual basis states
   std::vector<bool> newstate(2*VOL);

   //assign v1 to zero
   v1.assign(v1.size(), 0.0);

   MKL_INT nB = Wind[sector].nBasis;
   for(i=0; i<nB; i++){
      if( fabs(v0[i]) < 1e-6 ) continue;
      newstate = Wind[sector].basisVec[i]; 
      /* act on the basis state with the Hamiltonian */
      /* a single plaquette is arranged as 
                pzw
             o-------o
             |       |
        pwx  |   p   |  pyz
             |       |
             o-------o
                pxy
      */
      for(p=0;p<VOL;p++){
         // check if the plaquette p is flippable 
         p1=2*p; p2=2*next[DIM+1][p]+1; p3=2*next[DIM+2][p]; p4=2*p+1;
         pxy=newstate[p1]; pyz=newstate[p2]; pzw=newstate[p3]; pwx=newstate[p4];
         if((pxy==pyz)&&(pzw==pwx)&&(pwx!=pxy)){
             // If flippable, act with the Hamiltonian 
             newstate[p1]=!newstate[p1]; newstate[p2]=!newstate[p2];
             newstate[p3]=!newstate[p3]; newstate[p4]=!newstate[p4];
             // check which state it is by scanning other states in the same sector
             q=Wind[sector].binscan(newstate);
             // construct the new vector 
             v1[q]=v1[q]-1.0;
             //flip back the plq
             newstate[p1]=!newstate[p1]; newstate[p2]=!newstate[p2];
             newstate[p3]=!newstate[p3]; newstate[p4]=!newstate[p4];
         }
      }
      v1[i] += lam*Wind[sector].nflip[i]; 
   }

}
