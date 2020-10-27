#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<iterator>
#include "define.h"

void print2file(int, int, FILE*);
// This routine prints out diagnostics of eigenstates which are infinite temperature
// states in the spectrum and having low entropy. The aim is to study if "scars" exist
void studyEvecs00(int sector, double cutoff){
   double targetEN, amp, prob;
   double check;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,sizet;
   FILE *fptr1,*fptr2;

   sizet = Wind[sector].trans_sectors;
   // scan for states with the same energy density as INIT=4
   targetEN = lam*VOL/2.0;
   printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
   num_Eigst=0;
   for(p=0; p<sizet; p++){
     if(fabs(Wind[sector].evals_K00[p]-targetEN) < 1e-10){
        num_Eigst++;
        ev_list.push_back(p);
     }
   }
   printf("#-of eigenstates found in momentum (0,0)= %d \n",num_Eigst);
   printf("#-of ice states above cutoffProb = %lf in each eigenstate: \n",cutoff);
   fptr1 = fopen("EvecList00.dat","w");
   fptr2 = fopen("BasisList00.dat","w");
   for(p=0; p<num_Eigst; p++){
     totbasisState=0;
     check=0.0;
     for(q=0; q<sizet; q++){
       amp    = Wind[sector].evecs_K00[ev_list[p]*sizet + q];
       prob   = amp*amp;
       if(prob > cutoff) { totbasisState++; print2file(sector, q, fptr2);  }
       fprintf(fptr1,"% .6le ",amp);
       check += prob;
     }
     if( fabs(check-1.0) > 1e-10) printf("Warning! Normalization of evector %d is %.12lf\n",ev_list[p],check);
     fprintf(fptr1,"\n");
     printf("Eigenstate = %d, #-of ice states= %d\n",ev_list[p],totbasisState);
   }
   fclose(fptr1); fclose(fptr2);
}

// Same routine, but deals with the momentum (Pi,Pi) sector
void studyEvecsPiPi(int sector, double cutoff){
   double targetEN, amp, prob;
   double check;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,sizet;
   FILE *fptr1,*fptr2;

   sizet = Wind[sector].trans_sectors;
   // scan for states with the same energy density as INIT=4
   targetEN = lam*VOL/2.0;
   printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
   num_Eigst=0;
   for(p=0; p<sizet; p++){
     if(fabs(Wind[sector].evals_KPiPi[p]-targetEN) < 1e-10){
        num_Eigst++;
        ev_list.push_back(p);
     }
   }
   printf("#-of eigenstates found in momentum (Pi,Pi) = %d \n",num_Eigst);
   printf("#-of ice states above cutoffProb = %lf in each eigenstate: \n",cutoff);
   fptr1 = fopen("EvecListPiPi.dat","w");
   fptr2 = fopen("BasisListPiPi.dat","w");
   for(p=0; p<num_Eigst; p++){
     totbasisState=0;
     check=0.0;
     for(q=0; q<sizet; q++){
       amp    = Wind[sector].evecs_KPiPi[ev_list[p]*sizet + q];
       prob   = amp*amp;
       if(prob > cutoff){ totbasisState++; print2file(sector, q, fptr2);  }
       fprintf(fptr1,"% .6le ",amp);
       check += prob;
     }
     if( fabs(check-1.0) > 1e-10) printf("Warning! Normalization of evector %d is %.12lf\n",ev_list[p],check);
     fprintf(fptr1,"\n");
     printf("Eigenstate = %d, #-of ice states= %d\n",ev_list[p],totbasisState);
   }
   fclose(fptr1); fclose(fptr2);
}

// Same routine, but deals with the momentum (Pi,0) sector
void studyEvecsPi0(int sector, double cutoff){
   double targetEN, amp, prob;
   double check, chkPi0;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,sizet;
   FILE *fptr1, *fptr2;

   sizet = Wind[sector].trans_sectors;
   // scan for states with the same energy density as INIT=4
   targetEN = lam*VOL/2.0;
   printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
   num_Eigst=0; chkPi0  = 0;
   for(p=0; p<sizet; p++){
     if(spurPi0[chkPi0] == p) chkPi0++;
     else{   // skip spurious states
        if(fabs(Wind[sector].evals_KPi0[p]-targetEN) < 1e-10){
          num_Eigst++;
          ev_list.push_back(p);
        }
     }
   }
   printf("#-of eigenstates found in momentum (Pi,0) = %d \n",num_Eigst);
   printf("#-of ice states above cutoffProb = %lf in each eigenstate: \n",cutoff);
   fptr1 = fopen("EvecListPi0.dat","w");
   fptr2 = fopen("BasisListPi0.dat","w");
   for(p=0; p<num_Eigst; p++){
     totbasisState=0;
     check=0.0;
     for(q=0; q<sizet; q++){
       amp    = Wind[sector].evecs_KPi0[ev_list[p]*sizet + q];
       prob   = amp*amp;
       if(prob > cutoff) { totbasisState++; print2file(sector, q, fptr2);  }
       fprintf(fptr1,"% .6le ",amp);
       check += prob;
     }
     if( fabs(check-1.0) > 1e-10) printf("Warning! Normalization of evector %d is %.12lf\n",ev_list[p],check);
     fprintf(fptr1,"\n");
     printf("Eigenstate = %d, #-of ice states= %d\n",ev_list[p],totbasisState);
   }
   fclose(fptr1); fclose(fptr2);
}

// Same routine, but deals with the momentum (0,Pi) sector
void studyEvecs0Pi(int sector, double cutoff){
   double targetEN, amp, prob;
   double check, chk0Pi;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,sizet;
   FILE *fptr1,*fptr2;

   sizet = Wind[sector].trans_sectors;
   // scan for states with the same energy density as INIT=4
   targetEN = lam*VOL/2.0;
   printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
   num_Eigst=0; chk0Pi  = 0;
   for(p=0; p<sizet; p++){
     if(spur0Pi[chk0Pi] == p) chk0Pi++;
     else{   // skip spurious states
        if(fabs(Wind[sector].evals_K0Pi[p]-targetEN) < 1e-10){
          num_Eigst++;
          ev_list.push_back(p);
        }
     }
   }
   printf("#-of eigenstates found in momentum (0,Pi) = %d \n",num_Eigst);
   printf("#-of ice states above cutoffProb = %lf in each eigenstate: \n",cutoff);
   fptr1 = fopen("EvecList0Pi.dat","w");
   fptr2 = fopen("BasisList0Pi.dat","w");
   for(p=0; p<num_Eigst; p++){
     totbasisState=0;
     check=0.0;
     for(q=0; q<sizet; q++){
       amp    = Wind[sector].evecs_K0Pi[ev_list[p]*sizet + q];
       prob   = amp*amp;
       if(prob > cutoff){ totbasisState++; print2file(sector, q, fptr2);  }
       fprintf(fptr1,"% .6le ",amp);
       check += prob;
     }
     if( fabs(check-1.0) > 1e-10) printf("Warning! Normalization of evector %d is %.12lf\n",ev_list[p],check);
     fprintf(fptr1,"\n");
     printf("Eigenstate = %d, #-of ice states= %d\n",ev_list[p],totbasisState);
   }
   fclose(fptr1); fclose(fptr2);
}
