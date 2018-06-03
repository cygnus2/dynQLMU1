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
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<stdio.h>
#include<vector>
#include<iterator>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
//#define N 5
//#define LDA N

/* Main program */
int main() {
       
    unsigned int k;
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> vals;
    /* read matrix A from a binary file */
    std::ifstream inFile ("sparse_mat.bin", std::ios::in | std::ios::binary); 
    int nrows, ncols;
    inFile.read((char*)&nrows,sizeof(int));
    inFile.read((char*)&ncols,sizeof(int));
    std::cout<<"No of rows = "<<nrows<<std::endl;
    std::cout<<"No of cols = "<<ncols<<std::endl;
    rows.resize(nrows);
    cols.resize(ncols);
    vals.resize(ncols);
    inFile.read((char*)&rows[0],rows.size()*sizeof(int));
    inFile.read((char*)&cols[0],cols.size()*sizeof(int));
    inFile.read((char*)&vals[0],vals.size()*sizeof(double));
    inFile.close();
    /* print to check if matrix is read in correctly */
    std::cout<<"rows = ";
    for(k=0;k<rows.size();k++) std::cout<<rows[k]<<" ";
    std::cout<<" "<<std::endl;
    std::cout<<"cols = ";
    for(k=0;k<cols.size();k++) std::cout<<cols[k]<<" ";
    std::cout<<" "<<std::endl;
    std::cout<<"Hamiltonian = ";
    for(k=0;k<vals.size();k++) std::cout<<vals[k]<<" ";
    std::cout<<" "<<std::endl;

   
    /* construct the full matrix */
    
	/* Locals */
	//MKL_INT n = N, lda = LDA, info;
	/* Local arrays */
	//double w[N];
	//double a[LDA*N] = {
	//    1.96,  0.00,  0.00,  0.00,  0.00,
	//   -6.49,  3.80,  0.00,  0.00,  0.00,
	//   -0.47, -6.39,  4.17,  0.00,  0.00,
	//  -7.20,  1.50, -1.51,  5.70,  0.00,
	//   -0.65, -6.34,  2.67,  1.80, -7.10
	//};
	/* Executable statements */
	//printf( "LAPACKE_dsyev (column-major, high-level) Example Program Results\n" );
	/* Solve eigenproblem */
	//info = LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V', 'U', n, a, lda, w );
	/* Check for convergence */
	//if( info > 0 ) {
	//	printf( "The algorithm failed to compute eigenvalues.\n" );
	//	exit( 1 );
	//}
	/* Print eigenvalues */
	//print_matrix( "Eigenvalues", 1, n, w, 1 );
	/* Print eigenvectors */
	//print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
	exit( 0 );
} /* End of LAPACKE_dsyev Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}
