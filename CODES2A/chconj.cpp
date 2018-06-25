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
extern void print2file(MKL_INT,MKL_INT,std::vector<double>&);
extern void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

void chconj(int sector){
  //void diag_LAPACK_RRR(MKL_INT, MKL_INT, std::vector<double>&, std::vector<double>&, std::vector<double>&);
  
  // check the correct sector
  if((Wind[sector].Wx!=0)&&(Wind[sector].Wy!=0)){
    printf("This routine does not work in this sector. \n");
    exit(0);
  }

  MKL_INT i,j,q,count,N2;
  double ele;
  MKL_INT pos, N, NSQ;
  
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
  std::cout<<"Going to construct the matrix in the CC-basis"<<std::endl;
  N2= Wind[sector].nBasis;
  N = Wind[sector].nBasis/2;
  NSQ = N*N;
  std::vector<double> a(NSQ,0.0);
  if(N2%2) std::cout<<"Hilbert space in (0,0) sector not even!"<<std::endl;
  else std::cout<<"Matrix size in C=+ and (0,0) ="<<N<<std::endl;

  //Wind[sector].check_getH();
  /* initialize the matrix */
  for(i=0;i<N;i++){
  for(j=0;j<N;j++){
     ele = 0.5*(Wind[sector].getH(i,j) + Wind[sector].getH(i,N2-1-j) 
         + Wind[sector].getH(N2-1-i,j) + Wind[sector].getH(N2-1-i,N2-1-j));
     pos = i*N+j;
     a.at(pos) = ele;
  }}
  
  // space for eigenvalues and eigenvectors
  //Wind[sector].evals.resize(N,0.0);
  //Wind[sector].evecs.resize(NSQ,0.0); 

  MKL_INT LDZ, LDA, NSELECT, info;    
  MKL_INT il, iu, m;
  MKL_INT LWORK=-1;
  double abstol, vl, vu;
  //std::vector<double> w(1);
  std::vector<double> w(N,0.0);
  std::vector<double> z(NSQ,0.0);
  std::vector<MKL_INT> isuppz(2*N,0); 
   /* set array sizes */
   LDZ = N;
   NSELECT = N;
   LDA = N;
   abstol = -1;
   double* A = &a[0];
   MKL_INT* ISUPPZ = &isuppz[0];
   double* W = &w[0];
   double* Z = &z[0];
   //double* W = &Wind[sector].evals[0];
   //double* Z = &Wind[sector].evecs[0];
   // estimate workspace
   //info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'A', 'U', N, A, LDA,
   //                    0., 0., 0, 0, abstol, &m, W, Z, LDZ, ISUPPZ );
   //LWORK = static_cast<int>(work.front());
   //std::cout<<"Workspace required = "<<LWORK<<std::endl;
   //try {
   //   w.resize(LWORK);
   //} catch (const std::exception& e) {
   //   std::cerr << "Failed to allocate work space." << std::endl;
   //   exit(1);
   //}
   printf("Going to RRR routine. \n");
   //info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'A', 'U', N, A, LDA,
   //                    vl, vu, il, iu, abstol, &m, W, Z, LDZ, ISUPPZ );
   info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'A', 'U', N, A, LDA,
                       0., 0., 0, 0, abstol, &m, W, Z, LDZ, ISUPPZ );
   // check eigenvectors
   //check_eigvecs(N,a,w,z);

   //printf("Going to Divide-and-Conquer routine. \n");
   //info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', N, A, LDA, W );
   /* Check for convergence */
   if( info > 0 ) {
       printf( "The algorithm failed to compute eigenvalues.\n" );
       exit( 1 );
   }
   printf("Back from the routine and working properly. Info = %ld. Eval[0]=%e\n",info,W[0]);
   //fileprint_matrix( "Eigenvalues.dat", 1, N, W, 1 );
   isuppz.clear();

  // although these routines seem trivial, they are extremely important to
  // check the allocation of the eigenvectors and eigenvalues (on small lattices)
  // print copied evals and evecs
  //eig_print(Wind[sector].evals,Wind[sector].evecs,N);

  // calculate the O_flip on the eigenvectors
  unsigned int p;
  double Oflip_avg;
  double v_q,O_q,O_l;
  FILE *outf;
  outf = fopen("Oflip.dat","w");
  for(p=0;p<N;p++){
    Oflip_avg = 0.0;
    for(q=0;q<N;q++){
      v_q = z[p*N+q];
      O_q = Wind[sector].nflip[q];
      O_l = Wind[sector].nflip[2*N-1-q];
      Oflip_avg += 0.5*(v_q*v_q*O_q + v_q*v_q*O_l);
    }
    Oflip_avg = Oflip_avg/((double)VOL);
    fprintf(outf,"%.12lf %.12lf\n",w[p],Oflip_avg);
  }
 fclose(outf);

 w.clear();
 z.clear();
}

 // label the states which are C-partners with a common flag. 
 // this routine initializes the flag vector with -1.
 //void WindNo::initCflag(){
 //   for(unsigned i=0;i<nBasis;i++) Cflag.push_back(0); 
 //}

extern void print2file(MKL_INT N, MKL_INT NSQ, std::vector<double> &hamil){
  std::ofstream Outfile;
  Outfile.open("sparse_mat.bin",std::ios::out);
  Outfile.write((char*)&N, sizeof(MKL_INT));
  Outfile.write((char*)&NSQ, sizeof(MKL_INT));
  Outfile.write((char*)&hamil[0], hamil.size()*sizeof(double));
  Outfile.close();
}

//extern void check_eigvecs(MKL_INT, std::vector<double>& matrix, std::vector<double>& evals, std::vector<double>& evecs){
//  MKL_INT i,j;
//  double tryeval;
//  for(i=0;i<N;i++){
//   tryeval=0.0; 
//   for(j=0;j<N;j++){
//  }
//  }
//}

