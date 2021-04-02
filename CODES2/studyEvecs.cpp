#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<iterator>
#include "define.h"

void print2file(int, int, FILE*);
void diag_LAPACK_RRR2(int, std::vector<double>&, std::vector<double>&, std::vector<double>&);
// This routine prints out diagnostics of eigenstates which are infinite temperature
// states in the spectrum and having low entropy. The aim is to study if "scars" exist
void studyEvecs(int sector){
   double targetEN,NF;
   double cutoff, amp, prob;
   double check;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,sizet;
   FILE *fptr1,*fptr2;

   cutoff = 0.001;
   sizet = Wind[sector].nBasis;
   // scan for states with the same energy density as INIT=4; NF=Nplaq/2
   targetEN = lam*VOL/2.0;  NF=VOL/2.0;
   // the second possibility is only for the Ly=4 ladders
   //NF=3.0*VOL/8.0;  targetEN = lam*NF;
   printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
   num_Eigst=0;
   for(p=0;p<sizet;p++){
     if(fabs(Wind[sector].evals[p]-targetEN) < 1e-10){
        num_Eigst++;
        ev_list.push_back(p);
     }
   }
   printf("#-of eigenstates found = %d \n",num_Eigst);
   //printf("#-of ice states above cutoffProb = %lf in each eigenstate: \n",cutoff);
   fptr1 = fopen("EvecList.dat","w");
   fptr2 = fopen("BasisList.dat","w");
   for(p=0; p<num_Eigst; p++){
     totbasisState=0;
     check=0.0;
     fprintf(fptr2,"Important Basis States of eigenvector %d\n",ev_list[p]);
     for(q=0; q<sizet; q++){
       amp    = Wind[sector].evecs[ev_list[p]*sizet + q];
       prob   = amp*amp;
       check += prob;
       if(prob > cutoff){
          totbasisState++;
	        //printf("#-of flippable states in basis state=%d is %d\n",q,Wind[sector].nflip[q]);
          print2file(sector, q, fptr2);
       }
       fprintf(fptr1,"%.6le ",prob);
     }
     if( fabs(check-1.0) > 1e-10) printf("Normalization of evector %d is %.12lf\n",ev_list[p],check);
     fprintf(fptr1,"\n");
     printf("Eigenstate = %d, #-of ice states= %d\n",ev_list[p],totbasisState);
   }
   fclose(fptr1);
   fclose(fptr2);
   // check if these eigenstates are eigenstates of the Okin and Opot separately
   for(p=0; p<num_Eigst; p++){
     check=0.0;
     for(q=0; q<sizet; q++){
        amp    = Wind[sector].evecs[ev_list[p]*sizet + q];
        amp    = amp*(Wind[sector].nflip[q] - NF);
        check  += amp*amp;
     }
     printf("norm || Opot|psi> - NF|psi> || = %.12le \n",check);
   }
}

void studyEvecsLy4(int sector){
   double targetEN,NF;
   double cutoff, amp, prob;
   double check;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,sizet;
   FILE *fptr1,*fptr2;

   cutoff = 0.001;
   sizet = Wind[sector].nBasis;
   // the second possibility is only for the Ly=4 ladders
   if(LX==4){
     NF=3.0*VOL/8.0;  targetEN = lam*NF;
   }
   else if(LX==6){
     NF=5.0*VOL/12.0; targetEN = lam*NF;
   }
   else{
     printf("No scars for this lattice at this energy \n");
     return;
   }
   printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
   num_Eigst=0;
   for(p=0;p<sizet;p++){
     if(fabs(Wind[sector].evals[p]-targetEN) < 1e-10){
        num_Eigst++;
        ev_list.push_back(p);
     }
   }
   printf("#-of eigenstates found = %d \n",num_Eigst);
   //printf("#-of ice states above cutoffProb = %lf in each eigenstate: \n",cutoff);
   fptr1 = fopen("EvecList2.dat","w");
   fptr2 = fopen("BasisList2.dat","w");
   for(p=0; p<num_Eigst; p++){
     totbasisState=0;
     check=0.0;
     fprintf(fptr2,"Important Basis States of eigenvector %d\n",ev_list[p]);
     for(q=0; q<sizet; q++){
       amp    = Wind[sector].evecs[ev_list[p]*sizet + q];
       prob   = amp*amp;
       check += prob;
       if(prob > cutoff){
          totbasisState++;
	        //printf("#-of flippable states in basis state=%d is %d\n",q,Wind[sector].nflip[q]);
          print2file(sector, q, fptr2);
       }
       fprintf(fptr1,"%.6le ",prob);
     }
     if( fabs(check-1.0) > 1e-10) printf("Normalization of evector %d is %.12lf\n",ev_list[p],check);
     fprintf(fptr1,"\n");
     printf("Eigenstate = %d, #-of ice states= %d\n",ev_list[p],totbasisState);
   }
   fclose(fptr1);
   fclose(fptr2);
   // check if these eigenstates are eigenstates of the Okin and Opot separately
   for(p=0; p<num_Eigst; p++){
     check=0.0;
     for(q=0; q<sizet; q++){
        amp    = Wind[sector].evecs[ev_list[p]*sizet + q];
        amp    = amp*(Wind[sector].nflip[q] - NF);
        check  += amp*amp;
     }
     printf("norm || Opot|psi> - NF|psi> || = %.12le \n",check);
   }
}

void studyEvecs2(int sector){
  int i,j,p;
  double cutoff;
  double cI, cJ;
  // nZero counts the zero modes of the oKin;
  int nZero, nZero2, sizet;
  std::vector<int> ev_list;
  std::vector<double> Opot;
  std::vector<double> evalPot, evecPot;

  sizet = Wind[sector].nBasis;
  cutoff = 1e-10;
  // store the eigenvector labels whose eigenvalues are zero
  nZero=0;
  for(i=0; i<sizet; i++){
    if(fabs(Wind[sector].evals[i]) < cutoff){
       nZero++; ev_list.push_back(i);
     }
  }
  std::cout<<"#-of zero modes ="<<nZero<<std::endl;
  //for(i=0; i<nZero; i++) std::cout<<"identified state = "<<ev_list[i]<<std::endl;

  // construct the Opot( nZero x nZero ) matrix
  nZero2 = nZero*nZero;
  for(i=0; i<nZero2; i++) Opot.push_back(0.0);

  for(i=0; i<nZero; i++){
    for(j=0; j<nZero; j++){
      for(p=0; p<sizet; p++){
        cI = Wind[sector].evecs[ev_list[i]*sizet + p];
        cJ = Wind[sector].evecs[ev_list[j]*sizet + p];
        Opot[i+j*nZero] += cI*cJ*Wind[sector].nflip[p];
      }
    }
  }
  // diagonalize Opot operator in the zero mode subspace
  diag_LAPACK_RRR2(nZero, Opot, evalPot, evecPot);

  // clear memory
  Opot.clear(); evalPot.clear(); evecPot.clear();
  ev_list.clear();
}
