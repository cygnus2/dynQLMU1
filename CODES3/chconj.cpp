/* compute the charge conjugation parity of a given eigenstate */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// the zero winding sector has the additional discrete charge conjugation symmetry
void checkCCpartners(int sector){
  int i,q,flag,r,s;
  int size;

  size = Wind[sector].nBasis;
  std::vector<bool> flipstate(2*VOL);
  for(i=0; i<size; i++){
    flipstate = Wind[sector].basisVec[i];
    flipstate.flip();
    q=Wind[sector].binscan(flipstate);
    r=Wind[sector].Tflag[i]; s=Wind[sector].Tflag[q];
    //std::cout<<"(i,q) = "<<i<<" , "<<q<<";  Tflags ="<<r<<" , "<<s<<std::endl;
    if(q != (size-1-i)){
      std::cout<<"CC partner of "<<i<<" is not detected at "<<size-1-i<<std::endl;
      exit(0);
    }
  }
}

void ChargeConjEval1(int sector, std::vector<double> &CCj00, std::vector<double> &CCjPiPi){
  int p,q,r,s,k;
  int sizet,tsect;
  double ele;
  double CC00,CCPiPi;

  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  // the charge conjugation matrix in the bag basis; initalize first
  std::vector<std::vector<double>> CCconj00(tsect, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> CCconjPiPi(tsect, std::vector<double>(tsect, 0.0));

  // run a loop over nBasis to construct the C matrix in bag basis, then compute
  // charge conjugation expectation values the eigenvectors
  for(p=0; p<sizet; p++){
    q=sizet-1-p; // CC partner for state q
    r=Wind[sector].Tflag[p]-1; s=Wind[sector].Tflag[q]-1;
    ele=sqrt(Wind[sector].Tdgen[p]*Wind[sector].Tdgen[q])/((double)VOL);
    CCconj00[r][s]   += ele; // phases are all +1 for mom (0,0)
    CCconjPiPi[r][s] += ele*Wind[sector].FPiPi[p]*Wind[sector].FPiPi[q];
  }
  // now compute the expectation values over eigenstates
  for(k=0; k<tsect; k++){
    CC00 = 0.0; CCPiPi = 0.0;
    for(p=0; p<tsect; p++){
    for(q=0; q<tsect; q++){
      CC00   += CCconj00[p][q]*Wind[sector].evecs_K00[k*tsect+p]*Wind[sector].evecs_K00[k*tsect+q];
      CCPiPi += CCconjPiPi[p][q]*Wind[sector].evecs_KPiPi[k*tsect+p]*Wind[sector].evecs_KPiPi[k*tsect+q];
    }}
    //std::cout<<"Eigenstate ="<<k<<"; CConj evals = "<<CC00<<" , "<<CCPiPi<<std::endl;
    if(fabs((fabs(CC00)-1.0)) > 1e-10) printf("Not a good C-eigenstate of (0,0) \n");
    else CCj00.push_back(CC00);
    if(fabs((fabs(CCPiPi)-1.0)) > 1e-10) printf("Not a good C-eigenstate of (Pi,Pi) \n");
    else CCjPiPi.push_back(CCPiPi);
  }
  // clear memory
  CCconj00.clear(); CCconjPiPi.clear();
}

void ChargeConjEval2(int sector, std::vector<double> &CCjPi0, std::vector<double> &CCj0Pi){
  int p,q,r,s,k;
  int sizet,tsect;
  double ele;
  double CCPi0,CC0Pi;

  sizet =  Wind[sector].nBasis;
  tsect =  Wind[sector].trans_sectors;
  // the charge conjugation matrix in the bag basis; initalize first
  std::vector<std::vector<double>> CCconjPi0(tsect, std::vector<double>(tsect, 0.0));
  std::vector<std::vector<double>> CCconj0Pi(tsect, std::vector<double>(tsect, 0.0));

  // run a loop over nBasis to construct the C matrix in bag basis, then compute
  // charge conjugation expectation values the eigenvectors
  for(p=0; p<sizet; p++){
    q=sizet-1-p; // CC partner for state q
    r=Wind[sector].Tflag[p]-1; s=Wind[sector].Tflag[q]-1;
    ele=sqrt(Wind[sector].Tdgen[p]*Wind[sector].Tdgen[q])/((double)VOL);
    CCconjPi0[r][s] += ele*Wind[sector].FPi0[p]*Wind[sector].FPi0[q];
    CCconj0Pi[r][s] += ele*Wind[sector].F0Pi[p]*Wind[sector].F0Pi[q];
  }
  // now compute the expectation values over eigenstates
  for(k=0; k<tsect; k++){
    CCPi0 = 0.0; CC0Pi = 0.0;
    for(p=0; p<tsect; p++){
    for(q=0; q<tsect; q++){
      CCPi0 += CCconjPi0[p][q]*Wind[sector].evecs_KPi0[k*tsect+p]*Wind[sector].evecs_KPi0[k*tsect+q];
      CC0Pi += CCconj0Pi[p][q]*Wind[sector].evecs_K0Pi[k*tsect+p]*Wind[sector].evecs_K0Pi[k*tsect+q];
    }}
    std::cout<<"Eigenstate ="<<k<<"; CConj evals = "<<CCPi0<<" , "<<CC0Pi<<std::endl;
    if(fabs((fabs(CCPi0)-1.0)) > 1e-10) printf("Not a good C-eigenstate of (Pi,0) \n");
    else CCjPi0.push_back(CCPi0);
    if(fabs((fabs(CC0Pi)-1.0)) > 1e-10) printf("Not a good C-eigenstate of (0,Pi) \n");
    else CCj0Pi.push_back(CC0Pi);
  }
  // clear memory
  CCconjPi0.clear(); CCconj0Pi.clear();
}

void calcCCvalues(int sector){
  int i,j,tsect;
  int chkPi0, chk0Pi;
  FILE *fptr1,*fptr2;
  // vectors to store the charge conjugate eigenvalue
  std::vector<double> CCj00, CCjPiPi;
  std::vector<double> CCjPi0, CCj0Pi;
  tsect = Wind[sector].trans_sectors;
  // compute the charge conjugation eigenvalues in each sector
  if(INIT == 0 || INIT == 4){
    ChargeConjEval1(sector, CCj00, CCjPiPi);
    fptr1 = fopen("CCvals00.dat","w");
    fptr2 = fopen("CCvalsPiPi.dat","w");
    for(i=0; i<tsect; i++){
      fprintf(fptr1,"%.12lf % .4lf\n",Wind[sector].evals_K00[i],CCj00[i]);
      fprintf(fptr2,"%.12lf % .4lf\n",Wind[sector].evals_KPiPi[i],CCjPiPi[i]);
    }
    fclose(fptr1); fclose(fptr2);
  }
  if(INIT == 4){
    ChargeConjEval2(sector, CCjPi0, CCj0Pi);
    fptr1 = fopen("CCvalsPi0.dat","w");
    fptr2 = fopen("CCvals0Pi.dat","w");
    // initialize variables needed to track the spurious eigenstates
    chkPi0 = 0; chk0Pi = 0;
    for(i=0; i<tsect; i++){
       if(spurPi0[chkPi0] == i) chkPi0++;
       else{
         fprintf(fptr1,"%.12lf % .4lf\n",Wind[sector].evals_KPi0[i],CCjPi0[i]);
       }
    }
    for(i=0; i<tsect; i++){
       if(spur0Pi[chk0Pi] == i) chk0Pi++;
       else{
         fprintf(fptr2,"%.12lf % .4lf\n",Wind[sector].evals_K0Pi[i],CCj0Pi[i]);
       }
    }
    fclose(fptr1); fclose(fptr2);
  }
  // free other vectors
  CCj00.clear(); CCjPiPi.clear(); CCjPi0.clear(); CCj0Pi.clear();
}
