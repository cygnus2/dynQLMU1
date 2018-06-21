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
#include "mkl_lapacke.h"
#include "mkl.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Main program */
int main() {
       
    unsigned int k;
    unsigned int i,j;
    unsigned int curr_row, curr_col, n_entries, col_index;
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> vals;

    /* read matrix A from a binary file */
    std::ifstream inFile ("sparse_mat.bin", std::ios::in | std::ios::binary); 
    int nrows, ncols;
    inFile.read((char*)&nrows,sizeof(int));
    inFile.read((char*)&ncols,sizeof(int));
    std::cout<<"No of rows = "<<(nrows-1)<<std::endl;
    std::cout<<"No of cols = "<<ncols<<std::endl;
    rows.resize(nrows);
    cols.resize(ncols);
    vals.resize(ncols);
    inFile.read((char*)&rows[0],rows.size()*sizeof(int));
    inFile.read((char*)&cols[0],cols.size()*sizeof(int));
    inFile.read((char*)&vals[0],vals.size()*sizeof(double));
    inFile.close();

    /* print to check if matrix is read in correctly */
    //std::cout<<"rows = ";
    //for(k=0;k<rows.size();k++) std::cout<<rows[k]<<" ";
    //std::cout<<" "<<std::endl;
    //std::cout<<"cols = ";
    //for(k=0;k<cols.size();k++) std::cout<<cols[k]<<" ";
    //std::cout<<" "<<std::endl;
    //std::cout<<"Hamiltonian = ";
    //for(k=0;k<vals.size();k++) std::cout<<vals[k]<<" ";
    //std::cout<<" "<<std::endl;

    MKL_INT N, NSQ, LDA, LDZ, NSELECT, info;    
    MKL_INT il, iu, m;
    double abstol, vl, vu;

    /* set array sizes */
    N = nrows-1;
    NSQ = N*N;
    LDA = N;
    LDZ = N;
    NSELECT = N;

    /* allocate space for the arrays */
    std::vector<MKL_INT> isuppz(2*N); 
    std::vector<double> w(N);
    std::vector<double> z(NSQ);
    std::vector<double> a(NSQ);
   
    /* initialize the matrix */
    for(i=0;i<N*LDA;i++) a.push_back(i);

    /* reconstruct the full matrix */
    col_index = 0;
    for(i=0;i<N;i++){
       curr_row  = i;
       n_entries = rows[i+1]-rows[i];
       for(j=0;j<n_entries;j++){
         curr_col = cols[col_index]-1;
         a[curr_row*N + curr_col] = vals[col_index];
         col_index++;
       }
    }
    std::cout<<"Full matrix constructed"<<std::endl;
    /* print full matrix */
    //print_matrix("Full matrix", N, N, A, LDA);

    /* Executable statements */
    printf( "LAPACKE_dsyevr (column-major, high-level) RRR Example Program Results\n" );

    /* Solve eigenproblem */
    MKL_INT* ISUPPZ = &isuppz[0];
    double* A = &a[0];
    double* W = &w[0];
    double* Z = &z[0];
    info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'A', 'U', N, A, LDA,
                        vl, vu, il, iu, abstol, &m, W, Z, LDZ, ISUPPZ );


    printf("Back from the routine and working properly. Info = %d. Eval[0]=%f\n",info,W[0]);
 
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
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " % f", a[i+j*lda] );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a matrix */
void fileprint_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
        FILE *outfile;
        outfile = fopen(desc,"w");
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) fprintf(outfile, " %f \n", a[i+j*lda] );
		fprintf(outfile, "\n" );
	}
}
