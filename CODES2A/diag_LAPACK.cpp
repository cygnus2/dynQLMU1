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

void diag_LAPACK(int sector, std::vector<double>& matrix, std::vector<double>& evals, std::vector<double>& evecs){
    unsigned i,j;
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
    printf("Back from the routine and working properly. Info = %d. Eval[0]=%f\n",info,W[0]);
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

void diag_LAPACK1(MKL_INT N, MKL_INT LDA, MKL_INT info, double *W, double *A){
    unsigned i,j;
    
    /* print full matrix */
    print_matrix("Full matrix", N, N, A, LDA);

    /* Executable statements */
    printf( "LAPACKE_dsyev (column-major, high-level) Example Program Results\n" );

    /* Solve eigenproblem */
    info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', N, A, LDA, W );
    printf("Back from the routine and working properly. Info = %d. Eval[0]=%f\n",info,W[0]);
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }

    /* Print eigenvalues */
    //fileprint_matrix( "Eigenvalues.dat", 1, N, W, 1 );
    print_matrix( "Eigenvalues", 1, N, W, 1 );
    /* Print eigenvectors */
    print_matrix( "Eigenvectors (stored columnwise)", N, N, A, LDA );
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
 };
