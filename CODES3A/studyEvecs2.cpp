#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<iterator>
#include "define.h"

void print2file(int, int, FILE*);

void studyEvecs2_K00(int sector, double cutoff){
  double targetEN, amp, prob, check;
  int num_Eigst, totbasisState;
  int i,j,p,q;
  double cI, cJ;
  FILE *fptr, *fptr1, *fptr2, *fptr3;
  // nZero counts the zero modes of the oKin in the momenta (0,0) basis;
  int nZero, nZero2, tsect, sizet;
  sizet = Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<int> ev_list;
  std::vector<int> scarList;
  // oflip stores the Opot (=Oflip) in the bag basis
  std::vector<double> oflip(tsect, 0.0);
  std::vector<double> evecBag(tsect, 0.0);
  std::vector<std::vector<double>> Opot;
  std::vector<double> init;
  std::vector<double> evalPot, evecPot;

  // Of in the bag basis: true for all momenta sectors
  for(p=0; p<sizet; p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }
  // store the eigenvector labels whose eigenvalues are zero
  nZero=0;
  for(i=0; i<tsect; i++){
    if(fabs(Wind[sector].evals_K00[i]) < 1e-10){
       nZero++; ev_list.push_back(i);
     }
  }
  std::cout<<"#-of zero modes ="<<nZero<<std::endl;
  // construct the Opot( nZero x nZero ) matrix
  nZero2 = nZero*nZero;
  for(i=0; i<nZero; i++) init.push_back(0.0);
  for(i=0; i<nZero; i++) Opot.push_back(init);
  std::cout<<"Size allocated for Opot ="<<Opot.size()<<std::endl;

  for(i=0; i<nZero; i++){
    for(j=0; j<nZero; j++){
      for(p=0; p<tsect; p++){
        cI = Wind[sector].evecs_K00[ev_list[i]*tsect + p];
        cJ = Wind[sector].evecs_K00[ev_list[j]*tsect + p];
        Opot[i][j] += cI*cJ*oflip[p];
      }
    }
  }
  // diagonalize Opot operator in the zero mode subspace
  diag_LAPACK_RRR(nZero, Opot, evalPot, evecPot);

  // print eigenvectors
  //fptr = fopen("ZeroMode00EigenF.dat","w");
  //for(i=0; i<nZero; i++){
  //  for(j=0; j<nZero; j++){
  //    fprintf(fptr,"%.12lf ",Wind[sector].evecs_K00[i*tsect + j]);
  //  }
  //  fprintf(fptr,"\n");
  //}
  //fclose(fptr);

  // fileprint eigenvalues; and locate scar states
  targetEN = VOL/2.0; num_Eigst=0;
  printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
  fptr = fopen("OpotEvals00.dat","w");
  for(i=0; i<nZero; i++){
    fprintf(fptr, "%.12lf \n",evalPot[i]);
    if(fabs(evalPot[i]-targetEN) < 1e-10){
       num_Eigst++;
       scarList.push_back(i);
    }
  }
  fclose(fptr);

  // print the eigenvector of the scars
  fptr3 = fopen("ZeroModeSuper00.dat","w");
  for(i=0; i<nZero; i++){
    if(fabs(evalPot[i]-targetEN) < 1e-10){ //this is a scar state
      for(j=0; j<nZero; j++){
        fprintf(fptr3, "% .6le ",evecPot[i*nZero + j]);
      }
      fprintf(fptr3, "\n");
    }
  }
  fclose(fptr3);

  // obtain and print the scar eigenstates
  printf("#-of eigenstates found in momentum (0,0)= %d \n",num_Eigst);
  printf("#-of ice states above cutoff Prob = %lf in each eigenstate: \n",cutoff);
  fptr1 = fopen("EvecList00.dat","w");
  fptr2 = fopen("BasisList00.dat","w");
  for(i=0; i<num_Eigst; i++){
    // express the scar eigenvector in the translation bag basis
    for(p=0; p<tsect; p++){
      evecBag[p] = 0.0;
      for(j=0; j<nZero; j++){
         evecBag[p] += evecPot[scarList[i]*nZero + j]*Wind[sector].evecs_K00[ev_list[j]*tsect + p];
      }
    }
    // check normalization, and for the special bag states
    totbasisState=0; check=0.0;
    for(q=0; q<tsect; q++){
      amp    = evecBag[q];
      prob   = amp*amp;
      if(prob > cutoff) { totbasisState++; print2file(sector, q, fptr2);  }
      fprintf(fptr1,"% .6le ",amp);
      check += prob;
    }
    if( fabs(check-1.0) > 1e-10) printf("Warning! Normalization of evector %d is %.12lf\n",ev_list[p],check);
    fprintf(fptr1,"\n");
    printf("Eigenstate of Opot =%d, #-of ice states= %d\n",ev_list[i],totbasisState);
  }
  fclose(fptr1); fclose(fptr2);

  // clear memory
  Opot.clear(); evalPot.clear(); evecPot.clear();
  ev_list.clear(); oflip.clear(); init.clear();
  evecBag.clear(); scarList.clear();
}

void studyEvecs2_KPiPi(int sector, double cutoff){
  double targetEN, amp, prob, check;
  int num_Eigst, totbasisState;
  int i,j,p,q;
  double cI, cJ;
  FILE *fptr, *fptr1, *fptr2;
  // nZero counts the zero modes of the oKin in the momenta (pi,pi) basis;
  int nZero, nZero2, tsect, sizet;
  sizet = Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<int> ev_list;
  std::vector<int> scarList;
  // oflip stores the Opot (=Oflip) in the bag basis
  std::vector<double> oflip(tsect, 0.0);
  std::vector<double> evecBag(tsect, 0.0);
  std::vector<std::vector<double>> Opot;
  std::vector<double> init;
  std::vector<double> evalPot, evecPot;

  // Of in the bag basis: true for all momenta sectors
  for(p=0; p<sizet; p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }
  // store the eigenvector labels whose eigenvalues are zero
  nZero=0;
  for(i=0; i<tsect; i++){
    if(fabs(Wind[sector].evals_KPiPi[i]) < 1e-10){
       nZero++; ev_list.push_back(i);
     }
  }
  std::cout<<"#-of zero modes ="<<nZero<<std::endl;
  // construct the Opot( nZero x nZero ) matrix
  nZero2 = nZero*nZero;
  for(i=0; i<nZero; i++) init.push_back(0.0);
  for(i=0; i<nZero; i++) Opot.push_back(init);
  std::cout<<"Size allocated for Opot ="<<Opot.size()<<std::endl;

  for(i=0; i<nZero; i++){
    for(j=0; j<nZero; j++){
      for(p=0; p<tsect; p++){
        cI = Wind[sector].evecs_KPiPi[ev_list[i]*tsect + p];
        cJ = Wind[sector].evecs_KPiPi[ev_list[j]*tsect + p];
        Opot[i][j] += cI*cJ*oflip[p];
      }
    }
  }
  // diagonalize Opot operator in the zero mode subspace
  diag_LAPACK_RRR(nZero, Opot, evalPot, evecPot);

  // fileprint eigenvalues; and locate scar states
  targetEN = VOL/2.0; num_Eigst=0.0;
  printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
  fptr = fopen("OpotEvalsPiPi.dat","w");
  for(i=0; i<nZero; i++){
    fprintf(fptr, "%.12lf \n",evalPot[i]);
    if(fabs(evalPot[i]-targetEN) < 1e-10){
       num_Eigst++;
       scarList.push_back(i);
    }
  }
  fclose(fptr);

  // obtain and print the scar eigenstates
  printf("#-of eigenstates found in momentum (Pi,Pi)= %d \n",num_Eigst);
  printf("#-of ice states above cutoff Prob = %lf in each eigenstate: \n",cutoff);
  fptr1 = fopen("EvecListPiPi.dat","w");
  fptr2 = fopen("BasisListPiPi.dat","w");
  for(i=0; i<num_Eigst; i++){
    // express the scar eigenvector in the translation bag basis
    for(p=0; p<tsect; p++){
      evecBag[p] = 0.0;
      for(j=0; j<nZero; j++){
         evecBag[p] += evecPot[scarList[i]*nZero + j]*Wind[sector].evecs_KPiPi[ev_list[j]*tsect + p];
      }
    }
    // check normalization, and for the special bag states
    totbasisState=0; check=0.0;
    for(q=0; q<tsect; q++){
      amp    = evecBag[q];
      prob   = amp*amp;
      if(prob > cutoff) { totbasisState++; print2file(sector, q, fptr2);  }
      fprintf(fptr1,"% .6le ",amp);
      check += prob;
    }
    if( fabs(check-1.0) > 1e-10) printf("Warning! Normalization of evector %d is %.12lf\n",ev_list[p],check);
    fprintf(fptr1,"\n");
    printf("Eigenstate of Opot =%d, #-of ice states= %d\n",ev_list[i],totbasisState);
  }
  fclose(fptr1); fclose(fptr2);

  // clear memory
  Opot.clear(); evalPot.clear(); evecPot.clear();
  ev_list.clear(); oflip.clear(); init.clear();
  evecBag.clear(); scarList.clear();
}

void studyEvecs2_KPi0(int sector, double cutoff){
  double targetEN, amp, prob, check;
  int num_Eigst, totbasisState;
  int i,j,p,q,pp;
  double cI, cJ;
  FILE *fptr, *fptr1, *fptr2;
  // nZero counts the zero modes of the oKin in the momenta (0,0) basis;
  int nZero, nZero2, tsect, sizet;
  sizet = Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<int> ev_list;
  std::vector<int> scarList;
  // oflip stores the Opot (=Oflip) in the bag basis
  std::vector<double> oflip(tsect, 0.0);
  std::vector<double> evecBag(tsect, 0.0);
  std::vector<std::vector<double>> Opot;
  std::vector<double> init;
  std::vector<double> evalPot, evecPot;

  // Of in the bag basis: true for all momenta sectors
  for(p=0; p<sizet; p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }
  // store the eigenvector labels whose eigenvalues are zero
  nZero=0;
  for(i=0; i<nDimPi0; i++){
    if(fabs(Wind[sector].evals_KPi0[i]) < 1e-10){
       nZero++; ev_list.push_back(i);
    }
  }
  std::cout<<"#-of zero modes ="<<nZero<<std::endl;
  // construct the Opot( nZero x nZero ) matrix
  nZero2 = nZero*nZero;
  for(i=0; i<nZero; i++) init.push_back(0.0);
  for(i=0; i<nZero; i++) Opot.push_back(init);
  std::cout<<"Size allocated for Opot ="<<Opot.size()<<std::endl;

  for(i=0; i<nZero; i++){
    for(j=0; j<nZero; j++){
      for(p=0; p<tsect; p++){
        pp = labelPi0[p];
        if(pp == -997) continue;
        cI = Wind[sector].evecs_KPi0[ev_list[i]*nDimPi0 + pp];
        cJ = Wind[sector].evecs_KPi0[ev_list[j]*nDimPi0 + pp];
        Opot[i][j] += cI*cJ*oflip[p];
      }
    }
  }
  // diagonalize Opot operator in the zero mode subspace
  diag_LAPACK_RRR(nZero, Opot, evalPot, evecPot);

  // fileprint non-spurious eigenvalues; and locate scar states
  targetEN = VOL/2.0; num_Eigst=0.0;
  printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
  fptr = fopen("OpotEvalsPi0.dat","w");
  for(i=0; i<nZero; i++){
     fprintf(fptr, "%.12lf\n",evalPot[i]);
     if(fabs(evalPot[i]-targetEN) < 1e-10){
       num_Eigst++;
       scarList.push_back(i);
    }
  }
  fclose(fptr);

  // obtain and print the scar eigenstates
  printf("#-of eigenstates found in momentum (Pi,0)= %d \n",num_Eigst);
  printf("#-of ice states above cutoff Prob = %lf in each eigenstate: \n",cutoff);
  fptr1 = fopen("EvecListPi0.dat","w");
  fptr2 = fopen("BasisListPi0.dat","w");
  for(i=0; i<num_Eigst; i++){
    // express the scar eigenvector in the translation bag basis
    for(p=0; p<tsect; p++){
      evecBag[p] = 0.0;
      pp = labelPi0[p];
      if(pp == -997) continue;
      for(j=0; j<nZero; j++){
         evecBag[p] += evecPot[scarList[i]*nZero + j]*Wind[sector].evecs_KPi0[ev_list[j]*nDimPi0 + pp];
      }
    }
    // check normalization, and for the special bag states
    totbasisState=0; check=0.0;
    for(q=0; q<tsect; q++){
      amp    = evecBag[q];
      prob   = amp*amp;
      if(prob > cutoff) { totbasisState++; print2file(sector, q, fptr2);  }
      fprintf(fptr1,"% .6le ",amp);
      check += prob;
    }
    if( fabs(check-1.0) > 1e-10) printf("Warning! Normalization of evector %d is %.12lf\n",ev_list[p],check);
    fprintf(fptr1,"\n");
    printf("Eigenstate of Opot =%d, #-of ice states= %d\n",ev_list[i],totbasisState);
  }
  fclose(fptr1); fclose(fptr2);

  // clear memory
  Opot.clear(); evalPot.clear(); evecPot.clear();
  oflip.clear(); init.clear();
  evecBag.clear(); scarList.clear();
}

void studyEvecs2_K0Pi(int sector, double cutoff){
  double targetEN, amp, prob, check;
  int num_Eigst, totbasisState;
  int i,j,p,q,pp;
  double cI, cJ;
  FILE *fptr, *fptr1, *fptr2;
  // nZero counts the zero modes of the oKin in the momenta (0,0) basis;
  int nZero, nZero2, tsect, sizet;
  sizet = Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<int> ev_list;
  std::vector<int> scarList;
  // oflip stores the Opot (=Oflip) in the bag basis
  std::vector<double> oflip(tsect, 0.0);
  std::vector<double> evecBag(tsect, 0.0);
  std::vector<std::vector<double>> Opot;
  std::vector<double> init;
  std::vector<double> evalPot, evecPot;

  // Of in the bag basis: true for all momenta sectors
  for(p=0; p<sizet; p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }
  // store the eigenvector labels whose eigenvalues are zero
  nZero=0;
  for(i=0; i<nDim0Pi; i++){
    if(fabs(Wind[sector].evals_K0Pi[i]) < 1e-10){
       nZero++; ev_list.push_back(i);
     }
  }
  std::cout<<"#-of zero modes ="<<nZero<<std::endl;
  // construct the Opot( nZero x nZero ) matrix
  nZero2 = nZero*nZero;
  for(i=0; i<nZero; i++) init.push_back(0.0);
  for(i=0; i<nZero; i++) Opot.push_back(init);
  std::cout<<"Size allocated for Opot ="<<Opot.size()<<std::endl;

  for(i=0; i<nZero; i++){
    for(j=0; j<nZero; j++){
      for(p=0; p<tsect; p++){
        pp = label0Pi[p];
        if(pp == -997) continue;
        cI = Wind[sector].evecs_K0Pi[ev_list[i]*nDim0Pi + pp];
        cJ = Wind[sector].evecs_K0Pi[ev_list[j]*nDim0Pi + pp];
        Opot[i][j] += cI*cJ*oflip[p];
      }
    }
  }
  // diagonalize Opot operator in the zero mode subspace
  diag_LAPACK_RRR(nZero, Opot, evalPot, evecPot);

  // fileprint non-spurious eigenvalues; and locate scar states
  targetEN = VOL/2.0; num_Eigst=0.0;
  printf("Looking for eigenstates with energy = %.12lf\n",targetEN);
  fptr = fopen("OpotEvals0Pi.dat","w");
  for(i=0; i<nZero; i++){
    fprintf(fptr, "%.12lf\n",evalPot[i]);
    if(fabs(evalPot[i]-targetEN) < 1e-10){
       num_Eigst++;
       scarList.push_back(i);
    }
  }
  fclose(fptr);

  // obtain and print the scar eigenstates
  printf("#-of eigenstates found in momentum (0,Pi)= %d \n",num_Eigst);
  printf("#-of ice states above cutoff Prob = %lf in each eigenstate: \n",cutoff);
  fptr1 = fopen("EvecList0Pi.dat","w");
  fptr2 = fopen("BasisList0Pi.dat","w");
  for(i=0; i<num_Eigst; i++){
    // express the scar eigenvector in the translation bag basis
    for(p=0; p<tsect; p++){
      evecBag[p] = 0.0;
      pp = label0Pi[p];
      if(pp == -997) continue;
      for(j=0; j<nZero; j++){
         evecBag[p] += evecPot[scarList[i]*nZero + j]*Wind[sector].evecs_K0Pi[ev_list[j]*nDim0Pi + pp];
      }
    }
    // check normalization, and for the special bag states
    totbasisState=0; check=0.0;
    for(q=0; q<tsect; q++){
      amp    = evecBag[q];
      prob   = amp*amp;
      if(prob > cutoff) { totbasisState++; print2file(sector, q, fptr2);  }
      fprintf(fptr1,"% .6le ",amp);
      check += prob;
    }
    if( fabs(check-1.0) > 1e-10) printf("Warning! Normalization of evector %d is %.12lf\n",ev_list[p],check);
    fprintf(fptr1,"\n");
    printf("Eigenstate of Opot =%d, #-of ice states= %d\n",ev_list[i],totbasisState);
  }
  fclose(fptr1); fclose(fptr2);

  // clear memory
  Opot.clear(); evalPot.clear(); evecPot.clear();
  ev_list.clear(); oflip.clear(); init.clear();
  evecBag.clear(); scarList.clear();
}
