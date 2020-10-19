#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<iterator>
#include "define.h"

void studyEvecs2_K00(int sector, double cutoff){
  int i,j,p;
  double cI, cJ;
  FILE *fptr;
  // nZero counts the zero modes of the oKin in the momenta (0,0) basis;
  int nZero, nZero2, tsect, sizet;
  sizet = Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<int> ev_list;
  // oflip stores the Opot (=Oflip) in the bag basis
  std::vector<double> oflip(tsect, 0.0);
  std::vector<std::vector<double>> Opot;
  std::vector<double> init;
  std::vector<double> evalPot, evecPot;

  // Of in the bag basis: true for all momenta sectors
  for(p=0; p<sizet; p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }
  // store the eigenvector labels whose eigenvalues are zero
  cutoff = 1e-10; nZero=0;
  for(i=0; i<tsect; i++){
    if(fabs(Wind[sector].evals_K00[i]) < cutoff){
       nZero++; ev_list.push_back(i);
     }
  }
  std::cout<<"#-of zero modes ="<<nZero<<std::endl;
  //for(i=0; i<nZero; i++) std::cout<<"identified state = "<<ev_list[i]<<std::endl;
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

  // fileprint eigenvalues
  fptr = fopen("OpotEvals00.dat","w");
  for(i=0; i<nZero; i++){
    fprintf(fptr, "%.12lf \n",evalPot[i]);
  }
  fclose(fptr);

  // clear memory
  Opot.clear(); evalPot.clear(); evecPot.clear();
  ev_list.clear(); oflip.clear(); init.clear();
}

void studyEvecs2_KPiPi(int sector, double cutoff){
  int i,j,p;
  double cI, cJ;
  FILE *fptr;
  // nZero counts the zero modes of the oKin in the momenta (0,0) basis;
  int nZero, nZero2, tsect, sizet;
  sizet = Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  std::vector<int> ev_list;
  // oflip stores the Opot (=Oflip) in the bag basis
  std::vector<double> oflip(tsect, 0.0);
  std::vector<std::vector<double>> Opot;
  std::vector<double> init;
  std::vector<double> evalPot, evecPot;

  // Of in the bag basis: true for all momenta sectors
  for(p=0; p<sizet; p++){
     oflip[Wind[sector].Tflag[p]-1] += Wind[sector].nflip[p]*Wind[sector].Tdgen[p]/((double)VOL);
  }
  // store the eigenvector labels whose eigenvalues are zero
  cutoff = 1e-10; nZero=0;
  for(i=0; i<tsect; i++){
    if(fabs(Wind[sector].evals_KPiPi[i]) < cutoff){
       nZero++; ev_list.push_back(i);
     }
  }
  std::cout<<"#-of zero modes ="<<nZero<<std::endl;
  //for(i=0; i<nZero; i++) std::cout<<"identified state = "<<ev_list[i]<<std::endl;
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

  // fileprint eigenvalues
  fptr = fopen("OpotEvalsPiPi.dat","w");
  for(i=0; i<nZero; i++){
    fprintf(fptr, "%.12lf \n",evalPot[i]);
  }
  fclose(fptr);

  // clear memory
  Opot.clear(); evalPot.clear(); evecPot.clear();
  ev_list.clear(); oflip.clear(); init.clear();
}
