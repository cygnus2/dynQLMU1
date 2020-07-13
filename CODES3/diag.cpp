#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<vector>
#include<iterator>
#include<cmath>
#include "mkl_lapacke.h"
#include "define.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void eig_print(std::vector<double>&,std::vector<double>&,int);
extern void check_eigvecs(MKL_INT, std::vector<double>&, std::vector<double>&, std::vector<double>&);

void diag_LAPACK(int size, std::vector<std::vector<double>>& matrix, std::vector<double>& evals,
      std::vector<double>& evecs){
    unsigned i,j;
    std::vector<double> acopy(size*size);
    // variabes to pass to LAPACK routine
    MKL_INT N, LDA, info;
    double *W, *A;
    N=size; LDA=N;
    W = (double*)malloc(N*sizeof(double));
    A = (double*)malloc(N*LDA*sizeof(double));
    for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A[i+j*N]=matrix[i][j];
      if(CHKDIAG) acopy[i+j*N]=matrix[i][j];
    }}
    /* print full matrix */
    //print_matrix("Full matrix", N, N, A, LDA);

    /* Executable statements */
    printf( "LAPACKE_dsyev (column-major, high-level) Example Program Results\n" );

    /* Solve eigenproblem */
    info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', N, A, LDA, W );
    printf("Back from the routine and working properly. Info = %ld. Eval[0]=%f\n",info,W[0]);
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }
    // copy eigenvalues and eigenvectors to main code
    evals.insert(evals.begin(), W, W+size);
    evecs.insert(evecs.begin(), A, A+size*size);

    // although these routines seem trivial, they are extremely important to
    // check the allocation of the eigenvectors and eigenvalues (on small lattices)
    // print copied evals and evecs
    //eig_print(evals,evecs,size);

    // check the eigenvectors. Be careful this requires an O(N^3) time
    if(CHKDIAG) check_eigvecs(N, acopy, evals, evecs);

    /* Print eigenvalues */
    //fileprint_matrix( "Eigenvalues_6x2.dat", 1, N, W, 1 );
    //print_matrix( "Eigenvalues", 1, N, W, 1 );
    /* Print eigenvectors */
    //print_matrix( "Eigenvectors (stored columnwise)", N, N, A, LDA );

    free(A); free(W); acopy.clear();
}

void diag_LAPACK_RRR(int size, std::vector<std::vector<double>>& matrix, std::vector<double>& evals,
  std::vector<double>& evecs){

  MKL_INT i, j;
  MKL_INT LDZ, LDA, NSELECT, info;
  MKL_INT il, iu, m;
  MKL_INT N,NSQ;
  double abstol, vl, vu;
  std::vector<double> acopy(size*size);
  double *W, *Z, *A;
  MKL_INT *ISUPPZ;

  N = size; LDA = N; LDZ = N;
  NSELECT = N;
  NSQ = N*N;
  W = (double*)malloc(N*sizeof(double));
  Z = (double*)malloc(N*LDA*sizeof(double));
  A = (double*)malloc(N*LDA*sizeof(double));
  ISUPPZ = (MKL_INT*)malloc(2*N*sizeof(MKL_INT));
  for(i=0;i<N;i++){
  for(j=0;j<N;j++){
    A[i+j*N]=matrix[i][j];
    if(CHKDIAG) acopy[i+j*N]=matrix[i][j];
  }}

  abstol = -1;
  printf("Going to RRR routine. \n");
  info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'A', 'U', N, A, LDA,
                       0., 0., 0, 0, abstol, &m, W, Z, LDZ, ISUPPZ );
  /* Check for convergence */
  if(info > 0){
     printf( "The algorithm failed to compute eigenvalues.\n" );
     exit(1);
  }
  printf("Back from the routine and working properly. Info = %ld. Eval[0]=%e\n",info,W[0]);
  // copy eigenvalues and eigenvectors to main code
  evals.insert(evals.begin(), W, W+N);
  evecs.insert(evecs.begin(), Z, Z+N*N);

  // check the eigenvectors. Be careful this requires an O(N^3) time
  if(CHKDIAG) check_eigvecs(N, acopy, evals, evecs);

  free(A); free(W); free(Z); free(ISUPPZ);
  acopy.clear();
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " % f", a[i+j*lda] );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a matrix to a file */
void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        FILE *outfile;
        outfile = fopen(desc,"w");
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) fprintf(outfile, " %f \n", a[i+j*lda] );
                fprintf(outfile, "\n" );
        }
}

void eig_print(std::vector<double> &evals,std::vector<double> &evecs,int size){
   unsigned i,j;
   printf("Eigenvalues \n");
   for(i=0;i<size;i++) printf("% lf \n",evals[i]);
   printf("\n");
   printf("Eigenvectors \n");
   for(i=0;i<size;i++){ for(j=0;j<size;j++){
      printf("% lf \n",evecs[i*size+j]);
    }
    printf("\n");
   }
 }

void check_eigvecs(MKL_INT size, std::vector<double> &matrix, std::vector<double> &evals,
  std::vector<double> &evecs){
  MKL_INT i,j,k;
  MKL_INT tryvec;
  FILE *outfile;
  std::vector<double> vec1(size,0.0);
  std::vector<double> vec2(size,0.0);
  double *v1, *v2;
  CBLAS_LAYOUT layout;
  CBLAS_UPLO uplo;
  double tryeval;
  double* aa = &matrix[0];
  layout = CblasColMajor;
  uplo   = CblasLower;

  // recalculate the eigenvalues by acting the Hamiltonian on the eigenvectors
  // check for any deviation
  if (CHKDIAG==1){
     outfile = fopen("eigencheck.dat","a");
     for(i=0;i<size;i++){
       for(j=0;j<size;j++){ k = i*size + j; vec1.at(j) = evecs[k]; }
       // do matrix multiplication y = alpha*A*x + beta*y; incx=incy=1
       v1 = &vec1[0]; v2 = &vec2[0];
       cblas_dsymv(layout, uplo, size, 1.0, aa, size, v1, 1, 0.0, v2, 1);
       tryeval = cblas_ddot(size, v2, 1, v1, 1) - evals[i];
       //if(tryeval > 1e-10) fprintf(outfile,"%ld %.12lf\n",i,tryeval);
       fprintf(outfile,"%ld %.12lf\n",i,tryeval);
     }
     fclose(outfile);
  }
  if(CHKDIAG==2){
     outfile = fopen("eigencheck.dat","a");
     for(i=0;i<10;i++){
       tryvec = std::floor((rand()/((double)RAND_MAX))*size);
       for(j=0;j<size;j++){ k = tryvec*size + j; vec1.at(j) = evecs[k]; }
       // do matrix multiplication y = alpha*A*x + beta*y; incx=incy=1
       v1 = &vec1[0]; v2 = &vec2[0];
       cblas_dsymv(layout, uplo, size, 1.0, aa, size, v1, 1, 0.0, v2, 1);
       tryeval = cblas_ddot(size, v2, 1, v1, 1) - evals[tryvec];
       //if(tryeval > 1e-10) fprintf(outfile,"%ld %.12lf\n",i,tryeval);
       fprintf(outfile,"%ld %.12lf\n",tryvec,tryeval);
     }
     fclose(outfile);
  }
  vec1.clear();
  vec2.clear();
}
