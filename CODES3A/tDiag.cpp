#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>
#include<chrono>

// basically the same subroutine as trans_Hamil_INIT0 but
void trans_diag(int sector){
  extern void diag_LAPACK(int,std::vector<std::vector<double>>&,std::vector<double>&,
      std::vector<double>&);
  //int i,j,
  int k,l;
  int tsect, sizet;
  double ele, norm;
  FILE *fptr;

  sizet = Wind[sector].nBasis;
  tsect = Wind[sector].trans_sectors;
  std::cout<<"=================================================================================="<<std::endl;
  std::cout<<"Constructing the Hamiltonian in the (kx,ky)=(0,0); (Pi,Pi) sectors."<<std::endl;
  std::cout<<"Dimension of matrix = "<<tsect<<std::endl;

  // allocate space for hamil_Kxy
  Wind[sector].allocate_Kxy(4);

  // Get starting timepoint
  auto start = std::chrono::high_resolution_clock::now();

  // construct the hamil_kxy matrix for sectors (0,0) and (pi,pi)
  for(std::size_t i=0;i<sizet;i++){
     k=Wind[sector].Tflag[i]-1;
     // off-diagonal elements
     for(std::size_t j=i+1;j<sizet;j++){
       ele = Wind[sector].getH(i,j);
       if(ele == 0) continue;
       l   = Wind[sector].Tflag[j]-1;
       norm= sqrt(Wind[sector].Tdgen[i]*Wind[sector].Tdgen[j])/((double)VOL);
       Wind[sector].hamil_K00[k][l] +=  2*ele*norm;
       if((Wind[sector].momPiPi[k]) && (Wind[sector].momPiPi[l]))
        Wind[sector].hamil_KPiPi[k][l] += 2*ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[j]*norm;
     }
     // diagonal elements
     ele = Wind[sector].getH(i,i);
     if(ele == 0) continue;
     norm= Wind[sector].Tdgen[i]/((double)VOL);
     Wind[sector].hamil_K00[k][k] +=  ele*norm;
     if(Wind[sector].momPiPi[k])
       Wind[sector].hamil_KPiPi[k][k] += ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[i]*norm;
  }

  std::cout<<"===================================================================="<<std::endl;
  std::cout<<"Constructing the Hamiltonian in the (kx,ky)= (Pi,0); (0,Pi) sectors."<<std::endl;
  dimPi0 = Wind[sector].trans_sectors - listPi0.size();
  dim0Pi = Wind[sector].trans_sectors - list0Pi.size();
  std::cout<<"Dimension of matrix (Pi,0)= "<<dimPi0<<std::endl;
  std::cout<<"Dimension of matrix (0,Pi)= "<<dim0Pi<<std::endl;

  // Get ending timepoint
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time taken for Hamiltonian construction in t-basis="<<duration.count()<< " secs"<<std::endl;
  // construct the hamil_kxy matrix
  for(std::size_t i=0;i<Wind[sector].nBasis;i++){
     k=Wind[sector].Tflag[i]-1;
     // though the flags start from 1; the indices start from 0
     //for(j=0;j<Wind[sector].nBasis;j++){
     // off-diagonal elements
     for(std::size_t j=i+1;j<Wind[sector].nBasis;j++){
       ele = Wind[sector].getH(i,j);
       if(ele == 0) continue;
       //if(ele==0.0) continue;
       l   = Wind[sector].Tflag[j]-1;
       norm= sqrt(Wind[sector].Tdgen[i]*Wind[sector].Tdgen[j])/((double)VOL);
       Wind[sector].hamil_K00[k][l] +=  2*ele*norm;
       if((Wind[sector].momPi0[k]) && (Wind[sector].momPi0[l]))
        Wind[sector].hamil_KPi0[k][l] += 2*ele*Wind[sector].FPi0[i]*Wind[sector].FPi0[j]*norm;
       if((Wind[sector].mom0Pi[k]) && (Wind[sector].mom0Pi[l]))
        Wind[sector].hamil_K0Pi[k][l] += 2*ele*Wind[sector].F0Pi[i]*Wind[sector].F0Pi[j]*norm;
       if((Wind[sector].momPiPi[k]) && (Wind[sector].momPiPi[l]))
        Wind[sector].hamil_KPiPi[k][l] += 2*ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[j]*norm;
   }
   // diagonal elements
   ele = Wind[sector].getH(i,i);
   if(ele == 0) continue;
   norm= Wind[sector].Tdgen[i]/((double)VOL);
   Wind[sector].hamil_K00[k][k] +=  ele*norm;
   if(Wind[sector].momPi0[k])
    Wind[sector].hamil_KPi0[k][k] += ele*Wind[sector].FPi0[i]*Wind[sector].FPi0[i]*norm;
   if(Wind[sector].mom0Pi[k])
    Wind[sector].hamil_K0Pi[k][k] += ele*Wind[sector].F0Pi[i]*Wind[sector].F0Pi[i]*norm;
   if(Wind[sector].momPiPi[k])
    Wind[sector].hamil_KPiPi[k][k] += ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[i]*norm;
  }

  // print matrix in translation basis
  //print_matrixTbasis( "Hamiltonian for (0,0)   sector", sector, 1 );
  //print_matrixTbasis( "Hamiltonian for (Pi,Pi) sector", sector, 4 );

  // clear file
  if(CHKDIAG){
    fptr=fopen("eigencheck.dat","w");
    fclose(fptr);
  }

  // remove memory for the original Hamiltonian which is not needed any more
  Wind[sector].hamil.clear();
  Wind[sector].rows.clear();
  Wind[sector].cols.clear();

  // diagonalize the matrixes with a LAPACK routine
  diag_LAPACK_RRR(Wind[sector].trans_sectors,Wind[sector].hamil_K00,
    Wind[sector].evals_K00,Wind[sector].evecs_K00);

  diag_LAPACK_RRR(Wind[sector].trans_sectors,Wind[sector].hamil_KPiPi,
    Wind[sector].evals_KPiPi,Wind[sector].evecs_KPiPi);

  // deallocate space for hamil_Kxy
  Wind[sector].deallocate_Kxy(INIT);
}

void trans_Hamil_INIT4(int sector){


  // construct the hamil_kxy matrix
  for(std::size_t i=0;i<Wind[sector].nBasis;i++){
     k=Wind[sector].Tflag[i]-1;
     // though the flags start from 1; the indices start from 0
     //for(j=0;j<Wind[sector].nBasis;j++){
     // off-diagonal elements
     for(std::size_t j=i+1;j<Wind[sector].nBasis;j++){
       ele = Wind[sector].getH(i,j);
       if(ele == 0) continue;
       //if(ele==0.0) continue;
       l   = Wind[sector].Tflag[j]-1;
       norm= sqrt(Wind[sector].Tdgen[i]*Wind[sector].Tdgen[j])/((double)VOL);
       Wind[sector].hamil_K00[k][l] +=  2*ele*norm;
       if((Wind[sector].momPi0[k]) && (Wind[sector].momPi0[l]))
        Wind[sector].hamil_KPi0[k][l] += 2*ele*Wind[sector].FPi0[i]*Wind[sector].FPi0[j]*norm;
       if((Wind[sector].mom0Pi[k]) && (Wind[sector].mom0Pi[l]))
        Wind[sector].hamil_K0Pi[k][l] += 2*ele*Wind[sector].F0Pi[i]*Wind[sector].F0Pi[j]*norm;
       if((Wind[sector].momPiPi[k]) && (Wind[sector].momPiPi[l]))
        Wind[sector].hamil_KPiPi[k][l] += 2*ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[j]*norm;
   }
   // diagonal elements
   ele = Wind[sector].getH(i,i);
   if(ele == 0) continue;
   norm= Wind[sector].Tdgen[i]/((double)VOL);
   Wind[sector].hamil_K00[k][k] +=  ele*norm;
   if(Wind[sector].momPi0[k])
    Wind[sector].hamil_KPi0[k][k] += ele*Wind[sector].FPi0[i]*Wind[sector].FPi0[i]*norm;
   if(Wind[sector].mom0Pi[k])
    Wind[sector].hamil_K0Pi[k][k] += ele*Wind[sector].F0Pi[i]*Wind[sector].F0Pi[i]*norm;
   if(Wind[sector].momPiPi[k])
    Wind[sector].hamil_KPiPi[k][k] += ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[i]*norm;
  }

  // Get ending timepoint
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time taken for Hamiltonian construction in t-basis="<<duration.count()<< " secs"<<std::endl;

  // print matrix in translation basis
  //print_matrixTbasis( "Hamiltonian for (0,0)   sector", sector, 1 );
  //print_matrixTbasis( "Hamiltonian for (Pi,0)  sector", sector, 2 );
  //print_matrixTbasis( "Hamiltonian for (0,Pi)  sector", sector, 3 );
  //print_matrixTbasis( "Hamiltonian for (Pi,Pi) sector", sector, 4 );

  // clear file
  if(CHKDIAG){
    fptr=fopen("eigencheck.dat","w");
    fclose(fptr);
  }

  // remove memory for the original Hamiltonian which is not needed any more
  Wind[sector].hamil.clear();
  Wind[sector].rows.clear();
  Wind[sector].cols.clear();

  // diagonalize the matrixes with a LAPACK routine
  // measure and print the times needed to do ED

  start = std::chrono::high_resolution_clock::now();
  diag_LAPACK_RRR(Wind[sector].trans_sectors,Wind[sector].hamil_K00,
    Wind[sector].evals_K00,Wind[sector].evecs_K00);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time to diag H in (0,0) sector="<<duration.count()<< " secs"<<std::endl;

  start = std::chrono::high_resolution_clock::now();
  diag_LAPACK_RRR(Wind[sector].trans_sectors,Wind[sector].hamil_KPi0,
    Wind[sector].evals_KPi0,Wind[sector].evecs_KPi0);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time to diag H in (Pi,0) sector="<<duration.count()<< " secs"<<std::endl;

  start = std::chrono::high_resolution_clock::now();
  diag_LAPACK_RRR(Wind[sector].trans_sectors,Wind[sector].hamil_K0Pi,
    Wind[sector].evals_K0Pi,Wind[sector].evecs_K0Pi);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time to diag H in (0,Pi) sector="<<duration.count()<< " secs"<<std::endl;

  start = std::chrono::high_resolution_clock::now();
  diag_LAPACK_RRR(Wind[sector].trans_sectors,Wind[sector].hamil_KPiPi,
    Wind[sector].evals_KPiPi,Wind[sector].evecs_KPiPi);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time to diag H in (Pi,Pi) sector="<<duration.count()<< " secs"<<std::endl;

  // deallocate space for hamil_Kxy
  Wind[sector].deallocate_Kxy(INIT);

  // routine to print eigenvalues and eigenvectors
  //printf("========= Sector (0,0) =======================================================\n");
  //print_eigsys(Wind[sector].trans_sectors, Wind[sector].evals_K00, Wind[sector].evecs_K00);
  //printf("========= Sector (Pi,Pi) =======================================================\n");
  //print_eigsys(Wind[sector].trans_sectors, Wind[sector].evals_KPiPi, Wind[sector].evecs_KPiPi);
  //printf("========= Sector (Pi,0) =======================================================\n");
  //print_eigsys(Wind[sector].trans_sectors, Wind[sector].evals_KPi0, Wind[sector].evecs_KPi0);
  //printf("========= Sector (0,Pi) =======================================================\n");
  //print_eigsys(Wind[sector].trans_sectors, Wind[sector].evals_K0Pi, Wind[sector].evecs_K0Pi);
  //check_eigsys(sector);
}
