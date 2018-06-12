/* construct the Hamiltonian matrix in the charge conjugate basis from the ice states */
// only works in the (Wx,Wy)=(0,0) sector
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "mkl_lapacke.h"
#include "define.h"

// About the charge conjugation operation: because of the binary
// sort, the basis states are arranged in a perfectly symmetric
// order, such that the first and the last state in the list are
// charge conjugation partners. Same for the second and the
// next-to-last state etc. This obviates the need for a separate
// flag; and an analytic formula can be written to compute the
// matrix elements in terms of the Hamiltonian in the ice basis
extern void eig_print(std::vector<double>&,std::vector<double>&,int);
void chconj(int sector){
  void diag_LAPACK1(MKL_INT, MKL_INT, MKL_INT, double*, double*);
  // check the correct sector
  if((Wind[sector].Wx!=0)&&(Wind[sector].Wy!=0)){
    printf("This routine does not work in this sector. \n");
    exit(0);
  }

  unsigned int i,j,q,count,N2;
  MKL_INT N, LDA, info;
  double *W,*A; 
  double ele;
  
  // putting flags is not explicitly needed due to sorting
  //std::vector<bool> newstate(2*VOL);
  //printf("No of basis states here = %ld \n",Wind[sector].nBasis);
  //Wind[sector].initCflag();
  // construct and label the charge conjugation partners
  //count=1;
  //for(i=0;i<Wind[sector].nBasis;i++){
  //    if(Wind[sector].Cflag[i]) continue;
  //    newstate = Wind[sector].basisVec[i];
  //    newstate.flip();
      //search for the charge conjugated spin with the fast binary search 
  //    q=Wind[sector].binscan(newstate);
  //    Wind[sector].Cflag[i]=count; 
  //    Wind[sector].Cflag[q]=count;
  //    count++;
  //}
  //count--;
  //if(count != Wind[sector].nBasis/2) std::cout<<"Hilbert space not halved!"<<std::endl;
  // print charge conjugation flags
  //for(i=0;i<Wind[sector].nBasis;i++) std::cout<<"ice states = "<<i<<"; flag ="<<Wind[sector].Cflag[i]<<std::endl;

  // construct the matrix in the charge conjugation basis
  N2= Wind[sector].nBasis;
  N = Wind[sector].nBasis/2;
  LDA = N;
  if(N2%2) std::cout<<"Hilbert space in (0,0) sector not even!"<<std::endl;
  W = (double*)malloc(N*sizeof(double));
  A = (double*)malloc(N*LDA*sizeof(double));
 
  for(i=0;i<N;i++){
  for(j=0;j<N;j++){
     A[i*N+j] = 0.5*(Wind[sector].getH(i,j) + Wind[sector].getH(i,N2-1-j) + Wind[sector].getH(N2-1-i,j) + Wind[sector].getH(N2-1-i,N2-1-j));
  }}
  
  // diagonalize matrix
  diag_LAPACK1(N,LDA,info,W,A);

  // copy eigenvalues and eigenvectors to main code
  Wind[sector].evals.insert(Wind[sector].evals.begin(), W, W+N);
  Wind[sector].evecs.insert(Wind[sector].evecs.begin(), A, A+N*N);

  // although these routines seem trivial, they are extremely important to
  // check the allocation of the eigenvectors and eigenvalues (on small lattices)
  // print copied evals and evecs
  //eig_print(Wind[sector].evals,Wind[sector].evecs,N);

  free(A); free(W);  
}

 // label the states which are C-partners with a common flag. 
 // this routine initializes the flag vector with -1.
 //void WindNo::initCflag(){
 //   for(unsigned i=0;i<nBasis;i++) Cflag.push_back(0); 
 //}

