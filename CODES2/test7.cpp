/*******************************************************************************
*   Copyright(C) 2010-2015 Intel Corporation. All Rights Reserved.
*   
*   The source code, information  and  material ("Material") contained herein is
*   owned  by Intel Corporation or its suppliers or licensors, and title to such
*   Material remains  with Intel Corporation  or its suppliers or licensors. The
*   Material  contains proprietary information  of  Intel or  its  suppliers and
*   licensors. The  Material is protected by worldwide copyright laws and treaty
*   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
*   modified, published, uploaded, posted, transmitted, distributed or disclosed
*   in any way  without Intel's  prior  express written  permission. No  license
*   under  any patent, copyright  or  other intellectual property rights  in the
*   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
*   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
*   intellectual  property  rights must  be express  and  approved  by  Intel in
*   writing.
*   
*   *Third Party trademarks are the property of their respective owners.
*   
*   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
*   this  notice or  any other notice embedded  in Materials by Intel or Intel's
*   suppliers or licensors in any way.
*
********************************************************************************
*/
/*
   LAPACKE_dsyev Example.
   ======================

   Program computes all eigenvalues and eigenvectors of a real symmetric
   matrix A:

     1.96  -6.49  -0.47  -7.20  -0.65
    -6.49   3.80  -6.39   1.50  -6.34
    -0.47  -6.39   4.17  -1.51   2.67
    -7.20   1.50  -1.51   5.70   1.80
    -0.65  -6.34   2.67   1.80  -7.10

   Description.
   ============

   The routine computes all eigenvalues and, optionally, eigenvectors of an
   n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.

   Example Program Results.
   ========================

 LAPACKE_dsyev (column-major, high-level) Example Program Results

 Eigenvalues
 -11.07  -6.23   0.86   8.87  16.09

 Eigenvectors (stored columnwise)
  -0.30  -0.61   0.40  -0.37   0.49
  -0.51  -0.29  -0.41  -0.36  -0.61
  -0.08  -0.38  -0.66   0.50   0.40
   0.00  -0.45   0.46   0.62  -0.46
  -0.80   0.45   0.17   0.31   0.16
*/
// this uses relatively robust representations 
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<vector>
#include<iterator>
#include<math.h>
#include "mkl_lapacke.h"
#include "mkl.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Main program */
int main() {
    MKL_INT k,pos;
    MKL_INT i,j;
    MKL_INT N, NSQ, LDA, LDZ, NSELECT, info;    
    MKL_INT il, iu, m;
    long double abstol, vl, vu;

    N=80000;
    NSQ=N*N;
    LDA=N;
    LDZ=N;

    /* allocate space for the arrays */
    std::vector<MKL_INT> isuppz(2*N,0); 
    std::vector<double> w(N,0.0);
    std::vector<double> z(NSQ,0.0);
    std::vector<double> a(NSQ,0.0);
   
    std::cout<<"size of vector "<<a.size()<<std::endl;
    std::cout<<"maximum capacity of vector "<<a.max_size()<<std::endl;
    // make the matrix
    for(i=0; i < N; i++){
      j=i; pos=i*N+j;
      a.at(pos) = 1.0;
      j=(i+1)%N; pos = i*N+j;
      a.at(pos) = -0.5; 
      j=(i-1+N)%N; pos = i*N+j;
      a.at(pos) = -0.5; 
     }

    double* A = &a[0];
    /* print full matrix */
    //print_matrix("Full matrix", N, N, A, LDA);

    /* Executable statements */
    printf( "LAPACKE_dsyevr (column-major, high-level) RRR Example Program Results\n" );

    /* Solve eigenproblem */
    abstol=-1;
    MKL_INT* ISUPPZ = &isuppz[0];
    double* W = &w[0];
    double* Z = &z[0];
    info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'A', 'U', N, A, LDA,
                        vl, vu, il, iu, abstol, &m, W, Z, LDZ, ISUPPZ );
    //info = LAPACKE_dsyevd( LAPACK_COL_MAJOR, 'V', 'U', N, A, LDA, W );

    printf("Back from the routine and working properly. Info = %d. Eval[0]=%e\n",info,W[0]);
    printf("No of eigenvalues = %d\n",m);
 
    /* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm failed to compute eigenvalues.\n" );
        exit( 1 );
    }
    /* Print eigenvalues */
    //print_matrix( "Eigenvalues", 1, N, W, 1 );
    fileprint_matrix( "EigenvaluesRRR_test.dat", 1, N, W, 1 );
    /* Print eigenvectors */
    //print_matrix( "Eigenvectors (stored columnwise)", N, N, A, LDA );

    isuppz.clear(); 
    w.clear();
    z.clear();
    a.clear();
    
    exit( 0 );
} /* End of LAPACKE_dsyevr Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* aa, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " % f", aa[i+j*lda] );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a matrix */
void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* aa, MKL_INT lda ) {
	MKL_INT i, j;
        FILE *outfile;
        outfile = fopen(desc,"w");
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) fprintf(outfile, " %le \n", aa[i+j*lda] );
		fprintf(outfile, "\n" );
	}
}
