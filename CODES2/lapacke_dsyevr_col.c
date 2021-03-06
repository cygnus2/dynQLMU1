/*******************************************************************************
* Copyright 2010-2016 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
   LAPACKE_dsyevr Example.
   =======================

   Program computes the smallest eigenvalues and the corresponding
   eigenvectors of a real symmetric matrix A using the Relatively Robust
   Representations, where A is:

     0.67  -0.20   0.19  -1.06   0.46
    -0.20   3.82  -0.13   1.06  -0.48
     0.19  -0.13   3.27   0.11   1.10
    -1.06   1.06   0.11   5.86  -0.98
     0.46  -0.48   1.10  -0.98   3.54

   Description.
   ============

   The routine computes selected eigenvalues and, optionally, eigenvectors of
   an n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.
   Eigenvalues and eigenvectors can be selected by specifying either a range
   of values or a range of indices for the desired eigenvalues.

   Example Program Results.
   ========================

 LAPACKE_dsyevr (column-major, high-level) Example Program Results

 The total number of eigenvalues found: 3

 Selected eigenvalues
   0.43   2.14   3.37

 Selected eigenvectors (stored columnwise)
  -0.98  -0.01  -0.08
   0.01   0.02  -0.93
   0.04  -0.69  -0.07
  -0.18   0.19   0.31
   0.07   0.69  -0.13
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
#define N 5
#define NSELECT 5
#define LDA N
#define LDZ N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, il, iu, m, lda = LDA, ldz = LDZ, info;
	double abstol, vl, vu;
	/* Local arrays */
	MKL_INT isuppz[2*N];
	double w[N], z[LDZ*NSELECT];
	double a[LDA*N] = {
	    0.67,  0.00,  0.00,  0.00,  0.00,
	   -0.20,  3.82,  0.00,  0.00,  0.00,
	    0.19, -0.13,  3.27,  0.00,  0.00,
	   -1.06,  1.06,  0.11,  5.86,  0.00,
	    0.46, -0.48,  1.10, -0.98,  3.54
	};
	/* Executable statements */
	printf( "LAPACKE_dsyevr (column-major, high-level) Example Program Results\n" );
	/* Negative abstol means using the default value */
	abstol = -1.0;
	/* Set il, iu to compute NSELECT smallest eigenvalues */
	il = 1;
	iu = NSELECT;
	/* Solve eigenproblem */
	info = LAPACKE_dsyevr( LAPACK_COL_MAJOR, 'V', 'I', 'U', n, a, lda,
			vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	/* Print the number of eigenvalues found */
	printf( "\n The total number of eigenvalues found:%2i\n", m );
	/* Print eigenvalues */
	print_matrix( "Selected eigenvalues", 1, m, w, 1 );
	/* Print eigenvectors */
	print_matrix( "Selected eigenvectors (stored columnwise)", n, m, z, ldz );
	exit( 0 );
} /* End of LAPACKE_dsyevr Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}
