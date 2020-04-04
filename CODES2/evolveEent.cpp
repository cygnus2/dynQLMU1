/* evolve Entanglement Entropy of cartoon states in real-time with the Hamiltonian */
/* if so chosen, this also computes the results in the diagonal ensemble           */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include<iterator>
#include "define.h"

extern void printconf1(std::vector<bool>);
extern void printconfA(std::vector<std::vector<bool>>);
extern void printconfB(std::vector<std::vector<bool>>);
extern void print_zmatrix( char*, MKL_INT, MKL_INT, MKL_Complex16*, MKL_INT ); 
extern void print_rmatrix( char* , MKL_INT , MKL_INT , double* , MKL_INT  );
extern void patch(std::vector<bool>&, std::vector<bool>&, std::vector<bool>&);
extern int checkGL2(std::vector<bool>&);
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void findBasisA(std::vector<std::vector<bool>>&, std::vector<bool>&);
extern void findBasisB(std::vector<std::vector<bool>>&, std::vector<bool>&);
extern void createLookupTable(int, MKL_INT, MKL_INT, std::vector<MKL_INT>&);
extern double schmidtDecom(std::vector<double>&, int, std::vector<MKL_INT>&);
extern double schmidtDecomRT(std::vector<double>&, std::vector<double>&, int, std::vector<MKL_INT>&);
extern void initState(int, int, int*);
extern double entanglementEntropy(int, int, std::vector<MKL_INT>&);

// local variables used 
unsigned int LEN_A,LEN_B,VOL_A,VOL_B;
unsigned int DA,DB,NCHI;
// basis for subsystems A and B
std::vector<std::vector<bool>> eA;
std::vector<std::vector<bool>> eB;


void evolve_Eent(int sector){
   int i,ix,iy,parity,p,q1,q2;
   int sizet,nchi,d,k,m;
   int initC;
   double t,EE,entE,diagEE;
   double initE;
   sizet = Wind[sector].nBasis;
   std::vector<double> alpha_real(sizet,0.0);
   std::vector<double> alpha_imag(sizet,0.0);
   // this stores which subsystem states map to states of the full system
   std::vector<MKL_INT> sub2main;
   FILE *outf;
   double overlap,norm;
   sizet = Wind[sector].nBasis; 

   /* construct cartoon state */
   initState(sector, INIT, &initC);
   std::cout<<"In routine evolve_Eent. Starting state is basis state = "<<initC<<std::endl;
   // calculate the average energy of the initial state
   initE=0.0;
   for(k=0; k<sizet; k++){
     initE += Wind[sector].evals[k]*Wind[sector].evecs[k*sizet+initC]*Wind[sector].evecs[k*sizet+initC]; 
   }

   // Construct the spin basis for the sub-systems
   LEN_B = LX - LEN_A;
   VOL_A = LEN_A*LY; VOL_B = LEN_B*LY;
   createBasis(sector);
   sub2main.assign(DA*DB, -5); // negative initial value
   createLookupTable(sector,DA,DB,sub2main);

   // Calculate diagonal ensemble
   diagEE = -100000;
   diagEE = entanglementEntropy(initC,sector,sub2main);

   outf = fopen("EENT_tevol.dat","w");
   fprintf(outf,"# avg energy of initial state = %.8lf\n",initE);
   fprintf(outf,"# diagonal entropy for initial state = %.8lf\n",diagEE);
   // calculate alpha(t)
   for(t=Ti; t<Tf; t=t+dT){
      alpha_real.assign(alpha_real.size(), 0.0);
      alpha_imag.assign(alpha_imag.size(), 0.0);
      for(k=0; k<sizet; k++){
      for(m=0; m<sizet; m++){
          if(fabs(Wind[sector].evals[m]) < 1e-10){
               alpha_real[k] += Wind[sector].evecs[m*sizet+k]*Wind[sector].evecs[m*sizet+initC];
          }
          else{
               alpha_real[k] += Wind[sector].evecs[m*sizet+k]*Wind[sector].evecs[m*sizet+initC]*cos(-Wind[sector].evals[m]*t);
               alpha_imag[k] += Wind[sector].evecs[m*sizet+k]*Wind[sector].evecs[m*sizet+initC]*sin(-Wind[sector].evals[m]*t);
          }
      }}
      // check norm
      //norm=0.0;
      //for(k=0; k<sizet; k++){ 
      //  norm += alpha_real[k]*alpha_real[k] + alpha_imag[k]*alpha_imag[k];
      //};
      //if( fabs(norm-1.0) >  1e-6) std::cout<<"t = "<<t<<" Norm = "<<norm<<std::endl;
      entE = schmidtDecomRT(alpha_real,alpha_imag,sector,sub2main);
      fprintf(outf,"%lf %lf\n",t,entE);
   }
   fclose(outf);
   alpha_real.clear(); alpha_imag.clear(); sub2main.clear();
   eA.clear(); eB.clear();
}

double schmidtDecomRT(std::vector<double> &alpha_real, std::vector<double> &alpha_imag,int sector, std::vector<MKL_INT> &sub2main){
  int i,j,p,k;
  int flagGI,count;
  double norm,EE;
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
  norm = 0.0;
  for(i=0;i<(DA*DB);i++) norm += (chi[i].real*chi[i].real + chi[i].imag*chi[i].imag);
  if( fabs(norm-1.0) >  1e-6) std::cout<<"Norm = "<<norm<<std::endl;

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

double schmidtDecom(std::vector<double> &vec, int sector, std::vector<MKL_INT> &sub2main){
  int i,j,p,k;
  int flagGI;
  double norm,EE;
  std::vector<bool> cA(2*VOL_A),cB(2*VOL_B),conf(2*VOL);
  double *chi, *chi_svd, *u, *vt;
  double *superb;
  MKL_INT M, N, LDA, info;
  MKL_INT dmin;
  MKL_INT ldu, ldvt;
  M = DA; N = DB;
  ldu = M; ldvt = N;
  if(DA < DB) dmin = DA;
  else dmin = DB;
  chi     = (double*)malloc((M*N)*sizeof(double));
  chi_svd = (double*)malloc(dmin*sizeof(double)); 
  superb  = (double*)malloc((dmin-1)*sizeof(double));

  // a simple check; remove it later
  //if(NCHI != dmin) std::cout<<"why does NCHI change from one evec to another?"<<std::endl;

  // initialize chi
  for(i=0;i<(DA*DB);i++){
     if(sub2main[i] == -5) chi[i] = 0.0;
     else                  chi[i] = vec[sub2main[i]];
  } 

  // check norm
  norm = 0.0;
  for(i=0;i<(DA*DB);i++) norm += chi[i]*chi[i];
  if( fabs(norm-1.0) >  1e-6) std::cout<<"Norm = "<<norm<<std::endl;
  
  //print_matrix("Density matrix", DA, DB, chi, DA);
  //printf( "LAPACKE_dgesvd (row-major, high-level) Example Program Results\n" );
  /* Compute SVD */
  //info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'A', 'A', DA, DB, chi, DB,
  //                chi_svd, u, ldu, vt, ldvt, superb );
  info = LAPACKE_dgesvd( LAPACK_ROW_MAJOR, 'N', 'N', DA, DB, chi, DB,
                  chi_svd, u, ldu, vt, ldvt, superb );
  /* Check for convergence */
  if( info > 0 ) {
         printf( "The algorithm computing SVD failed to converge.\n" );
         exit( 1 );
  } 
  /* Print singular values */
  //print_matrix( "Singular values", 1, dmin, chi_svd, 1 );
  // check norm
  norm = 0;
  for(i=0;i<dmin; i++) norm += chi_svd[i]*chi_svd[i];
  if( fabs(norm - 1.0) > 1e-6) std::cout<<"Norm of svd ="<<norm<<std::endl;
  
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

void createLookupTable(int sector, MKL_INT DA, MKL_INT DB, std::vector<MKL_INT> &sub2main){
   MKL_INT i,j,p;
   int flagGI;
   std::vector<bool> cA(2*VOL_A),cB(2*VOL_B),conf(2*VOL);
   int count;

   count=0;
   for(i=0; i<DA; i++){
   for(j=0; j<DB; j++){
      cA = eA[i]; cB = eB[j];
      patch(conf,cA,cB);
      // if Gauss Law not satisfied, skip the patching
      flagGI = checkGL2(conf);
      if(flagGI==0){ sub2main[i*DB+j]=-1; continue; }
      // match with the corresponding basis state in the winding number sector
      // it is very important that binscan2() is used. binscan() seems to give
      // wrong results. Why? 
      //p = Wind[sector].binscan2(conf);
      p = Wind[sector].scan(conf);
      if(p == -100) continue;
      else sub2main[i*DB+j] = p;
      count++;
   }}
   if(count != Wind[sector].nBasis) std::cout<<"possible error in schmidtDecom "<<std::endl;
}

// calculate the EE for all the states in a chosen sector for a fixed LA
// returns the result of the diagonal ensemble
// q1 is the initial state
double entanglementEntropy(int initC, int sector,  std::vector<MKL_INT> &sub2main){
  int p,i;
  int q1,q2;
  FILE *outf;
  MKL_INT sizet;
  std::vector<double> alpha;
  double EE,EENT_diag;
  double norm;
  // eigenvector and eigenvalue to calculate EE
  double sel_eval;
  std::vector<double> sel_evec;


  sizet = Wind[sector].nBasis;
  /* construct cartoon state and calculate the overlap with the eigenstates */
  for(p=0; p<sizet; p++){
     alpha.push_back(Wind[sector].evecs[p*sizet+initC]);
  }

  outf = fopen("EntE.dat","w");
  fprintf(outf,"# Entropy of Entanglement as a function of the eigenvalues for LA = %d\n",LEN_A);
  fprintf(outf,"# Eigenvalues  Entanglement Entropy \n");
 
  // Calculate the EE for each of the eigenstates 
  EENT_diag = 0.0;
  for(p=0; p<sizet; p++){
    // initialize the state and the entropy
    sel_evec.clear(); EE = -100.0;
    sel_eval = Wind[sector].evals[p];
    // two ways of copying the vector; both give same result
    for(i=0; i<sizet; i++){
       sel_evec.push_back(Wind[sector].evecs[p*sizet+i]);
    }
    //q1 = p*sizet; q2 = (p+1)*sizet; 
    //sel_evec.insert(sel_evec.begin(), Wind[sector].evecs.begin()+q1, Wind[sector].evecs.begin()+q2);
    //std::cout<<"Going to do Schmidt decompose eigenvector = "<< p << std::endl;
    //printvec(sel_evec);
    EE = schmidtDecom(sel_evec,sector,sub2main);
    EENT_diag += alpha[p]*alpha[p]*EE;
    // write to file
    fprintf(outf,"%lf %lf\n",sel_eval,EE);
  }
  fclose(outf);
  // clear the allocated memory in this routine 
  alpha.clear();  

  // return the result of the diagonal ensemble
  return EENT_diag;
}

// calculate the (gauge-invariant) basis of sub-systems A and B
void createBasis(int sector){
  int i,source,sink;
  int ix,iy;
  int p,p1;
  std::vector<bool> confA(2*VOL_A),confB(2*VOL_B),tconf(VOL2);

  // clear the subsytem basis before using them
  eA.clear(); eB.clear();

  for(i=0; i<Wind[sector].nBasis; i++){
    tconf = Wind[sector].basisVec[i];
    // collect A vectors    
    p1 = 0; 
    for(iy=0; iy<LY;    iy++){
    for(ix=0; ix<LEN_A; ix++){
      p = iy*LX + ix;
      confA[2*p1] = tconf[2*p]; confA[2*p1+1] = tconf[2*p+1];
      p1++;
    }}
    if((int)confA.size() != 2*VOL_A) std::cout<<"Size mismatch in A! "<<std::endl;
    eA.push_back(confA);
    // collect B vectors     
    p1 = 0;
    for(iy=0; iy<LY;    iy++){
    for(ix=0; ix<LEN_B; ix++){
      p = iy*LX + (ix+LEN_A);
      confB[2*p1] = tconf[2*p]; confB[2*p1+1] = tconf[2*p+1];
      p1++;
    }}
    if((int)confB.size() != 2*VOL_B) std::cout<<"Size mismatch in B! "<<std::endl;
    eB.push_back(confB);
  }
  // sort the sub-system basis vectors
  std::sort(eA.begin(),eA.end());
  std::sort(eB.begin(),eB.end());
  // remove duplicates
  eA.erase( unique(eA.begin(),eA.end() ), eA.end());
  eB.erase( unique(eB.begin(),eB.end() ), eB.end());
  DA = eA.size(); DB = eB.size();
  std::cout<<"Basis of subsystem A after sort+delete = "<<DA<<std::endl;
  std::cout<<"Basis of subsystem B after sort+delete = "<<DB<<std::endl;
  // initialize the number of SVD vectors
  if(DA < DB) NCHI = DA;
  else NCHI = DB;

  // print subsystem configuration
  //printconfA(eA);
  //printconfB(eB);
}

// Glue together the basis states of the subsystems
void patch(std::vector<bool> &conf, std::vector<bool> &cA, std::vector<bool> &cB){
   int ix,iy;
   int p,p1;
   //std::cout<<"New patch"<<std::endl;
   for(iy=0;iy<LY;iy++){
      for(ix=0;ix<LEN_A;ix++){
          p  = iy*LEN_A + ix;
          p1 = iy*LX    + ix;
          //std::cout<<ix<<" "<<iy<<" "<<p<<" "<<p1<<std::endl;
          conf[2*p1] = cA[2*p];  conf[2*p1+1] = cA[2*p+1];
      }
      for(ix=0;ix<LEN_B;ix++){
          p  = iy*LEN_B + ix;
          p1 = iy*LX    + (ix+LEN_A);
          //std::cout<<ix<<" "<<iy<<" "<<p<<" "<<p1<<std::endl;
          conf[2*p1] = cB[2*p];  conf[2*p1+1] = cB[2*p+1];
      }
   }
}

// Gauss Law routine to check which of the patched configs are valid
int checkGL2(std::vector<bool> &conf){
  int e1,e2,e3,e4;
  int p,pp,qx,qy;
  int ix,iy;
  for(p=0; p<VOL; p++){
  //ix=0;
  //for(iy=0; iy<LY; iy++){
  //  p = iy*LX+ix;
    qx=2*next[DIM-1][p]; 
    qy=2*next[DIM-2][p]+1; 
    pp=2*p;
    e1 = (conf[pp])? 1:-1; e2 = (conf[pp+1])? 1:-1; e3 = (conf[qx])? 1:-1; e4 = (conf[qy])? 1:-1;
    // if Q=e1+e2-e3-e4=0 return 1 (Gauss' Law OK) else return 0
    if(e1+e2-e3-e4) return 0; 
  }
  //ix=LEN_A;
  //for(iy=0; iy<LY; iy++){
  //  p = iy*LX+ix;
  //  qx=2*next[DIM-1][p]; 
  //  qy=2*next[DIM-2][p]+1; 
  //  pp=2*p;
  //  e1 = (conf[pp])? 1:-1; e2 = (conf[pp+1])? 1:-1; e3 = (conf[qx])? 1:-1; e4 = (conf[qy])? 1:-1;
    // if Q=e1+e2-e3-e4=0 return 1 (Gauss' Law OK) else return 0
  //  if(e1+e2-e3-e4) return 0; 
  //}
  return 1;
}

// old routines; not used in the code
// check routines. Finds a given basis state for the subsystem A
void findBasisA(std::vector<std::vector<bool>> &eA, std::vector<bool> &newstate){
     unsigned int m;
     // binary search of the sorted array  
     std::vector<std::vector<bool>>::iterator it; 
     it = std::find(eA.begin(),eA.end(),newstate);
     m  = std::distance(eA.begin(),it);
     if( m == eA.size() ){
       std::cout<<"Element not found!"<<std::endl;
     }
     else std::cout<<"Test state matches with Basis A found at "<<m<<std::endl;
}

// check routines. Finds a given basis state for the subsystem B
void findBasisB(std::vector<std::vector<bool>> &eB, std::vector<bool> &newstate){
     unsigned int m;
     // binary search of the sorted array  
     std::vector<std::vector<bool>>::iterator it; 
     it = std::find(eB.begin(),eB.end(),newstate);
     m  = std::distance(eB.begin(),it);
     if( m == eB.size() ){
       std::cout<<"Element not found!"<<std::endl;
     }
     else std::cout<<"Test state matches with Basis B found at "<<m<<std::endl;
}

