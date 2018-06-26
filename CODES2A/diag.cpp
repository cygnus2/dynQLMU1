#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<vector>
#include<iterator>
#include "mkl_lapacke.h"
#include "define.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void eig_print(std::vector<double>&,std::vector<double>&,int);
extern void check_eigvecs(MKL_INT, std::vector<double>&, std::vector<double>&, std::vector<double>&);

void diag_LAPACK(int sector, std::vector<double>& matrix, std::vector<double>& evals, 
 std::vector<double>& evecs){
    unsigned int i,j;
    // variabes to pass to LAPACK routine
    MKL_INT N, LDA, info;
    double *W, *A;
    N=Wind[sector].nBasis; 
    LDA=N;
    W = (double*)malloc(N*sizeof(double));
    A = (double*)malloc(N*LDA*sizeof(double));
    // convert the matrix to a long row
    for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      A[i*N+j]=Wind[sector].getH(i,j);
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
    evals.insert(evals.begin(), W, W+N);
    evecs.insert(evecs.begin(), A, A+N*N);

    // although these routines seem trivial, they are extremely important to
    // check the allocation of the eigenvectors and eigenvalues (on small lattices)
    // print copied evals and evecs
    //eig_print(evals,evecs,N);
    
    /* Print eigenvalues */
    //fileprint_matrix( "Eigenvalues.dat", 1, N, W, 1 );
    //print_matrix( "Eigenvalues", 1, N, W, 1 );
    /* Print eigenvectors */
    //print_matrix( "Eigenvectors (stored columnwise)", N, N, A, LDA );

    free(A); free(W);
}

void diag_LAPACK_RRR(MKL_INT N, MKL_INT NSQ, int sector, std::vector<double>& matrix, std::vector<double>& eval, std::vector<double>& evec){

  MKL_INT LDZ, LDA, NSELECT, info;    
  MKL_INT il, iu, m;
  double abstol, vl, vu;
  std::vector<double> acopy;
  std::vector<double> w(N,0.0);
  std::vector<double> z(NSQ,0.0);
  std::vector<MKL_INT> isuppz(2*N,0); 

  if(CHKDIAG) acopy = matrix; 
  /* set array sizes */
  LDZ = N;
  NSELECT = N;
  LDA = N;
  abstol = -1;
  double* A = &matrix[0];
  MKL_INT* ISUPPZ = &isuppz[0];
  double* W = &w[0];
  double* Z = &z[0];
  printf("Going to RRR routine. \n");
  info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'A', 'U', N, A, LDA,
                       0., 0., 0, 0, abstol, &m, W, Z, LDZ, ISUPPZ );

  printf("Back from the routine and working properly. Info = %ld. Eval[0]=%e\n",info,W[0]);
  //fileprint_matrix( "Eigenvalues.dat", 1, N, W, 1 );
  //eig_print(w,z,N);

  // check the eigenvectors. Be careful this requires an O(N^3) time
  check_eigvecs(N, acopy, w, z); 
  
  // clear memory
  w.clear();
  z.clear();
  isuppz.clear();

 }

// Notes on printing the eigenvalues and eigenvectors: Note that we have chosen a column-major
// layout, which means that the eigenvectors are arranged column-wsie. This means that the routine
// to print the matrix must have the form = i + j*N, where i is the index over row, and j is the
// index over columns. While printing the symmetric matrix, this however does not matter.
// Also, I think 'U' in the routine specifications means that the upper diagonal is free to be
// used as a workspace for the routine. Thus, only specifying the lower-diagonal is enough for
// the lapack routine.
/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* aa, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " % e", aa[i+j*lda] );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a matrix to a file */
void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* aa, MKL_INT lda ) {
        MKL_INT i, j;
        FILE *outfile;
        outfile = fopen(desc,"w");
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) fprintf(outfile, " % .12lf \n", aa[i+j*lda] );
                fprintf(outfile, "\n" );
        }
}

 void eig_print(std::vector<double> &evals,std::vector<double> &evecs,int size){
   unsigned i,j;
   printf("Eigenvalues \n");
   for(i=0;i<size;i++) printf("% le \n",evals[i]);
   printf("\n");
   printf("Eigenvectors \n");
   for(i=0;i<size;i++){ for(j=0;j<size;j++){
      printf("% le \n",evecs[i*size+j]);
    }
    printf("\n");
   }
 };

extern void check_eigvecs(MKL_INT size, std::vector<double> &matrix, std::vector<double> &evals, std::vector<double> &evecs){
  MKL_INT i,j,k;
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
  outfile = fopen("eigencheck.dat","w");
  for(i=0;i<size;i++){ 
   for(j=0;j<size;j++){ k = i*size + j; vec1.at(j) = evecs[k]; }
   // do matrix multiplication y = alpha*A*x + beta*y; incx=incy=1
   v1 = &vec1[0]; v2 = &vec2[0];
   cblas_dsymv(layout, uplo, size, 1.0, aa, size, v1, 1, 0.0, v2, 1);
   tryeval = cblas_ddot(size, v2, 1, v1, 1) - evals[i];
   if(tryeval > 1e-10) fprintf(outfile,"%ld %.12lf\n",i,tryeval);
  }
  fclose(outfile);
  vec1.clear();
  vec2.clear();  

}

