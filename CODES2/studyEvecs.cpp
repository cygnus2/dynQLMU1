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
   double targetEN, NF;
   double cutoff, amp, prob;
   double check,ele,vio;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,r,sizet;
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
          if(Wind[sector].nflip[q] != NF) printf("(%d, %d) \n",q,Wind[sector].nflip[q]);
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
   // check if these eigenstates are eigenstates of the Opot separately
   for(p=0; p<num_Eigst; p++){
     printf("Scar state %d = ",p);
     check=0.0;
     for(q=0; q<sizet; q++){
        amp    = Wind[sector].evecs[ev_list[p]*sizet + q];
        amp    = amp*(Wind[sector].nflip[q] - NF);
        check  += amp*amp;
     }
     printf("norm || Opot|psi> - NF|psi> || = %.12le \n",check);
     // check if they are eigenstates of Okin:
     // [Okin]_rq [psi^i]_q = amp_r;  vio = sum_r (amp_r^2)
     vio=0.0;
     for(r=0; r<sizet; r++){
       check=0.0;
       for(q=0; q<sizet; q++){
         if(r==q) continue;
         amp = Wind[sector].evecs[ev_list[p]*sizet + q];
         ele = Wind[sector].getH(r,q);
         //if(ev_list[p] == q) continue;
         check += amp*ele;
       }
       vio += check*check;
     }
     printf("norm || Okin |psi> || = %.12le\n", vio);
   }
}

void studyEvecsLy4(int sector){
   double targetEN,NF;
   double cutoff, amp, prob;
   double check,ele,vio;
   int num_Eigst;
   std::vector<int> ev_list;
   int totbasisState;
   int p,q,r,sizet;
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
          if(Wind[sector].nflip[q] != NF) printf("(%d, %d) \n",q,Wind[sector].nflip[q]);
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
   // check if these eigenstates are eigenstates of the Opot
   for(p=0; p<num_Eigst; p++){
     check=0.0;
     for(q=0; q<sizet; q++){
        amp    = Wind[sector].evecs[ev_list[p]*sizet + q];
        amp    = amp*(Wind[sector].nflip[q] - NF);
        check += amp*amp;
     }
     printf("norm || Opot|psi> - NF|psi> || = %.12le \n",check);
     // check if they are eigenstates of Okin
     // [Okin]_rq [psi^i]_q = amp^i_r;  vio[i] = sum_r (amp^i_r^2)
     vio=0.0;
     for(r=0; r<sizet; r++){
       check=0.0;
       for(q=0; q<sizet; q++){
         if(r==q) continue;
         amp = Wind[sector].evecs[ev_list[p]*sizet + q];
         ele = Wind[sector].getH(r,q);
         check += amp*ele;
       }
       vio += check*check;
     }
     printf("norm || Okin |psi> || = %.12le\n", vio);
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

// this routine does all the big bad checks needed to recognize scars as
// eigenvectors in the Okin = 0 manifold
void studyEvecs2_Ly4(int sector){
  int i,j,k,p;
  double cutoff;
  double cI, cJ;
  // (nZero, nTwo) counts the (zero,+/-2) modes of the oKin;
  int nZero, nTwo, sizet;
  int nTot, nTot2;
  double NF1, NF2;
  std::vector<int> ev_list;
  std::vector<double> Opot;
  std::vector<double> evalPot, evecPot;
  nTwo=0;

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
  // store the eigenvector labels whose eigenvalues are +2 and -2
  //nTwo=0;
  //for(i=0; i<sizet; i++){
  //  if( fabs((fabs(Wind[sector].evals[i])-2.0000)) < cutoff){
  //     nTwo++; ev_list.push_back(i);
  //  }
  //}
  //std::cout<<"#-of two modes ="<<nTwo<<std::endl;
  //for(i=0; i<nZero; i++) std::cout<<"identified state = "<<ev_list[i]<<std::endl;

  // construct the Opot( nTot x nTot ) matrix
  nTot  = nZero + nTwo;
  nTot2 = nTot*nTot;
  for(i=0; i<nTot2; i++) Opot.push_back(0.0);

  for(i=0; i<nTot; i++){
    for(j=0; j<nTot; j++){
      for(p=0; p<sizet; p++){
        cI = Wind[sector].evecs[ev_list[i]*sizet + p];
        cJ = Wind[sector].evecs[ev_list[j]*sizet + p];
        Opot[i+j*nTot] += cI*cJ*Wind[sector].nflip[p];
      }
    }
  }
  // diagonalize Opot operator in the zero+two mode subspace
  diag_LAPACK_RRR2(nTot, Opot, evalPot, evecPot);

  // Some basic sanity checks: normalization, eigenvalue equations
  double norm, ovl, vio, avgOp;
  std::vector<double> tvec(nTot, 0.0);
  int m,n;
  vio=0.0;
  // check orthornormality
  for(i=0; i<nTot; i++){
    norm=0.0;
    for(j=0; j<nTot; j++) norm += evecPot[i*nTot + j]*evecPot[i*nTot + j];
    if(fabs(norm-1.0) > 1e-10) printf("norm = %.12lf\n",norm);
  }
  for(i=0; i<nTot; i++){
    for(j=i+1; j<nTot; j++){
      ovl=0.0;
      for(k=0; k<nTot; k++){
        ovl += evecPot[i*nTot+k]*evecPot[j*nTot+k];
      }
      if(fabs(ovl) > 1e-10 ) printf("%d %d %lf\n",i,j,ovl);
    }
  }
  // eigenvalue equation
  for(i=0; i<nTot; i++){
    vio=0.0;
    for(m=0; m<nTot; m++){
      tvec[m]=0.0;
      for(n=0; n<nTot; n++) tvec[m] += Opot[m + n*nTot]*evecPot[i*nTot + n];
      ovl = tvec[m] - evalPot[i]*evecPot[i*nTot + m];
      vio += ovl*ovl;
    }
    if(vio > 1e-6) printf("eigenvector %d =%lf\n", i, vio);
  }

  // Express the scars in the electric flux basis; Nz=zero mode basis
  /* |Scar_i> = \sum_{k=1,..,Nz} c^i_k |z_k>
              = \sum_{k=i,..,Nz} c^i_k (\sum_{j=1,..,N} a^k_j |E_j>)
              = \sum_{j=1,..,N}  cJ |E_j>,
               where c_J = \sum_{i=1,..,Nz} c^i_k * a^k_j
   */
  NF1=VOL/2;
  int nS=0;
  for(i=0; i<nTot; i++){
    if(fabs(evalPot[i]-NF1) > cutoff) continue;
    printf("Scar state=%d with (Okin, Opot)=(0,%lf) \n",i, evalPot[i]);
    nS++;
    //check the contributing ice basis states in each scar
    norm=0.0; avgOp=0.0;
    for(j=0; j<sizet; j++){
       cJ=0.0;
       for(k=0; k<nTot; k++){
         cJ += evecPot[i*nTot + k]*Wind[sector].evecs[ev_list[k]*sizet + j];
       }
       norm += cJ*cJ;
       avgOp += cJ*cJ*Wind[sector].nflip[j];
       //if( (cJ*cJ) > 0.001) printf("basis state=%d, nFlip=%d, cJ=%lf\n",j,Wind[sector].nflip[j],cJ);
    }
    printf("Norm=%.12lf, <Opot>=%.12lf\n",norm,avgOp);
  }
  printf("#-of-scars=%d\n",nS);

  // clear memory
  Opot.clear(); evalPot.clear(); evecPot.clear();
  ev_list.clear();
}
