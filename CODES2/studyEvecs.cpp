#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<iterator>
#include "define.h"

void print2file(int, int, FILE*);
// This routine prints out diagnostics of eigenstates which are infinite energy
// states in the spectrum and having low entropy. The aim is to study if "scars"
// exist
void studyEvecs(int sector){
   double targetEN;
   double cutoff, amp, prob;
   double check;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,sizet;
   FILE *fptr1,*fptr2;

   cutoff = 0.07;
   sizet = Wind[sector].nBasis;
   // scan for states with the same energy density as INIT=4
   targetEN = lam*VOL/2.0;
   printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
   num_Eigst=0;
   for(p=0;p<sizet;p++){
     if(fabs(Wind[sector].evals[p]-targetEN) < 1e-10){
        num_Eigst++;
        ev_list.push_back(p);
     }
   }
   printf("#-of eigenstates found = %d \n",num_Eigst);
   printf("#-of ice states above cutoffProb = %lf in each eigenstate: \n",cutoff);
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
}
