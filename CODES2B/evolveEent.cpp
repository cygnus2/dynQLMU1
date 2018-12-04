/* evolve Entanglement Entropy of cartoon states in real-time with the Hamiltonian */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

extern void print_zmatrix( char*, MKL_INT, MKL_INT, MKL_Complex16*, MKL_INT ); 
extern void patch(std::vector<bool>&, std::vector<bool>&, std::vector<bool>&);
extern void print_zmatrix( char* , MKL_INT , MKL_INT , MKL_Complex16* , MKL_INT  );
extern void print_rmatrix( char* , MKL_INT , MKL_INT , double* , MKL_INT  );
extern int checkGL2(std::vector<bool>&);
extern void createLookupTable(int, MKL_INT, MKL_INT, std::vector<MKL_INT>&);

void evolve_Eent(int sector){
   
   int i,ix,iy,parity,p,q1,q2;
   int sizet,nchi,d,k,m;
   double t,entE;
   sizet = Wind[sector].nBasis;
   std::vector<bool> cart1(2*VOL);
   std::vector<double> alpha_real(sizet,0.0);
   std::vector<double> alpha_imag(sizet,0.0);
   // this stores which subsystem states map to states of the full system
   std::vector<MKL_INT> sub2main;
   FILE *outf;
   double overlap,norm;
  
   sizet = Wind[sector].nBasis; 
   /* construct cartoon state */
   for(iy=0;iy<LY;iy++){
   for(ix=0;ix<LX;ix++){
    parity=(ix+iy)%2;
    p = 2*(iy*LX+ix);
    if(parity){ cart1[p]=false; cart1[p+1]=true; }
    else{ cart1[p]=true; cart1[p+1]=false; }
   }}

   q1=Wind[sector].binscan(cart1);

   // Construct the spin basis for the sub-systems
   LEN_B = LX - LEN_A;
   VOL_A = LEN_A*LY; VOL_B = LEN_B*LY;
   createBasis(sector);
   sub2main.assign(DA*DB, -5); // negative initial value
   createLookupTable(sector,DA,DB,sub2main);

   outf = fopen("EENT_tevol.dat","w");
   // calculate alpha(t)
   for(t=Ti; t<Tf; t=t+dT){
      alpha_real.assign(alpha_real.size(), 0.0);
      alpha_imag.assign(alpha_imag.size(), 0.0);
      for(k=0; k<sizet; k++){
      for(m=0; m<sizet; m++){
          if(fabs(Wind[sector].evals[m]) < 1e-10){
               alpha_real[k] += Wind[sector].evecs[m*sizet+k]*Wind[sector].evecs[m*sizet+q1];
          }
          else{
               alpha_real[k] += Wind[sector].evecs[m*sizet+k]*Wind[sector].evecs[m*sizet+q1]*cos(-Wind[sector].evals[m]*t);
               alpha_imag[k] += Wind[sector].evecs[m*sizet+k]*Wind[sector].evecs[m*sizet+q1]*sin(-Wind[sector].evals[m]*t);

          }
      }}
      // check norm
      //norm=0.0;
      //for(k=0; k<sizet; k++){ 
      //  norm += alpha_real[k]*alpha_real[k] + alpha_imag[k]*alpha_imag[k];
      //};
      //if( fabs(norm-1.0) >  1e-6) std::cout<<"t = "<<t<<" Norm = "<<norm<<std::endl;
      entE = schmidtDecomRT(alpha_real,alpha_imag,eA,eB,sector,sub2main);
      fprintf(outf,"%lf %lf\n",t,entE);
   }
   fclose(outf);
   alpha_real.clear(); alpha_imag.clear(); sub2main.clear();
}

double schmidtDecomRT(std::vector<double> &alpha_real, std::vector<double> &alpha_imag, 
  std::vector<std::vector<bool>> &eA, std::vector<std::vector<bool>> &eB, int sector, std::vector<MKL_INT> &sub2main){

  int i,j,p,k;
  int flagGI,count;
  double norm,EE;
  std::vector<bool> cA(2*VOL_A),cB(2*VOL_B),conf(2*VOL);
  MKL_Complex16 *chi, *u, *vt;
  double *chi_svd, *superb;
  MKL_INT M, N, LDA, info;
  MKL_INT dmin;
  MKL_INT ldu, ldvt;
  M = DA; N = DB;
  ldu = M; ldvt = N;
  if(DA < DB) dmin = DA;
  else dmin = DB;
  chi     = (MKL_Complex16*)malloc((M*N)*sizeof(MKL_Complex16));
  //u       = (MKL_Complex16*)malloc((ldu*M)*sizeof(MKL_Complex16));
  //vt      = (MKL_Complex16*)malloc((ldvt*N)*sizeof(MKL_Complex16));
  chi_svd = (double*)malloc(dmin*sizeof(double)); 
  superb  = (double*)malloc((dmin-1)*sizeof(double));

  // note that this is ROW_MAJOR_REPRESENTATION in contrast to the representations
  // used in the diagonalization routines. 
  for(i=0; i<(DA*DB); i++){
     if(sub2main[i] == -5) chi[i] = {0.0 , 0.0};
     else chi[i] = {alpha_real[sub2main[i]], alpha_imag[sub2main[i]]};
  }
  // check norm
  //norm = 0.0;
  //for(i=0;i<(DA*DB);i++) norm += (chi[i].real*chi[i].real + chi[i].imag*chi[i].imag);
  //if( fabs(norm-1.0) >  1e-6) std::cout<<"Norm = "<<norm<<std::endl;

  //printf( "LAPACKE_dgesvd (row-major, high-level) Example Program Results\n" );
  /* Compute SVD */
  info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N', DA, DB, chi, DB,
                  chi_svd, u, ldu, vt, ldvt, superb );
  /* Check for convergence */
  if( info > 0 ) {
         printf( "The algorithm computing SVD failed to converge.\n" );
         exit( 1 );
  } 
  /* Print singular values */
  //print_rmatrix( "Singular values", 1, dmin, chi_svd, 1 );
  /* Print left singular vectors */
  //print_zmatrix( "Left singular vectors (stored columnwise)", M, M, u, ldu );
  /* Print right singular vectors */
  //print_zmatrix( "Right singular vectors (stored rowwise)", M, N, vt, ldvt );

  // check norm
  //norm = 0;
  //for(i=0;i<dmin; i++) norm += chi_svd[i]*chi_svd[i];
  //if( fabs(norm - 1.0) > 1e-6) std::cout<<"Norm of svd ="<<norm<<std::endl;

  EE=0.0;
  for(i=0; i<dmin; i++){
      if(chi_svd[i] < 1e-6) continue;
      EE -= chi_svd[i]*chi_svd[i]*log(chi_svd[i]*chi_svd[i]);
  }
  //std::cout<<"Entanglement Entropy = "<<EE<<std::endl;

  // clear memory
  free(chi); free(chi_svd); free(superb);
  //free(u); free(vt);
  return EE;
}

/* Auxiliary routine: printing a complex matrix */
void print_zmatrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", a[i*lda+j].real, a[i*lda+j].imag );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}

void createLookupTable(int sector, MKL_INT DA, MKL_INT DB, std::vector<MKL_INT> &sub2main){
   MKL_INT i,j,p;
   int flagGI;
   std::vector<bool> cA(2*VOL_A),cB(2*VOL_B),conf(2*VOL);

   for(i=0; i<DA; i++){
   for(j=0; j<DB; j++){
      cA = eA[i]; cB = eB[j];
      patch(conf,cA,cB);
      // if Gauss Law not satisfied, skip the patching
      flagGI = checkGL2(conf);
      if(flagGI==0){ sub2main[i*DB+j]=-1; continue; }
      // match with the corresponding basis state in the winding number sector
      p = Wind[sector].binscan2(conf);
      if(p == -100) continue;
      else sub2main[i*DB+j] = p;
   }}
}
