/* compute the entanglement entropy */
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
extern void createLookupTable(int, MKL_INT, MKL_INT, std::vector<MKL_INT>&);
extern double schmidtDecom(std::vector<double>&, int, std::vector<MKL_INT>&);

// variables used in this set of functions
unsigned int LEN_A,LEN_B,VOL_A,VOL_B;
unsigned int DA,DB,NCHI;
// basis for subsystems A and B
std::vector<std::vector<bool>> eA;
std::vector<std::vector<bool>> eB;
// entanglement Entropy for the eigenvector
double EE;

void entanglementEnt(int sector){
  int p,i;
  int sizet;
  FILE *outf;
  // eigenvector and eigenvalue to calculate EE
  double sel_eval;
  std::vector<double> sel_evec;


  sizet = Wind[sector].trans_sectors;
  // this stores which subsystem states map to states of the full system
  std::vector<MKL_INT> sub2main;

  // Construct the spin basis for the sub-systems
  LEN_B = LX - LEN_A;
  VOL_A = LEN_A*LY; VOL_B = LEN_B*LY;
  createBasis(sector);
  sub2main.assign(DA*DB, -5); // negative initial value
  createLookupTable(sector,DA,DB,sub2main);

  outf = fopen("EntE.dat","w");
  fprintf(outf,"# Entropy of Entanglement as a function of the eigenvalues for LA = %d\n",LEN_A);
  fprintf(outf,"# Eigenvalues  Entanglement Entropy \n");

  // Calculate the EE for each of the eigenstates
  for(p=0; p<sizet; p++){
    // initialize the state and the entropy
    sel_evec.clear(); EE = -100.0;
    sel_eval = Wind[sector].evals_K00[p];
    //sel_evec.assign(Wind[sector].evecs_K00[p*sizet], Wind[sector].evecs_K00[(p+1)*sizet-1]);
    for(i=0; i < sizet; i++){
       sel_evec.push_back(Wind[sector].evecs_K00[p*sizet + i]);
    }
    std::cout<<"Going to do Schmidt decompose eigenvector = "<< p << std::endl;
    std::cout<<"Size of the vector = "<<int(sel_evec.size())<<std::endl;
    //printvec(sel_evec);
    EE = schmidtDecom(sel_evec,sector,sub2main);
    // write to file
    fprintf(outf,"%lf %lf\n",sel_eval,EE);
  }
  fclose(outf);
  // free memory from the spin basis
  eA.clear();
  eB.clear();
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
      if(flagGI==0) continue;
      // match with the corresponding basis state in the winding number sector
      // it is very important that binscan2() is used. binscan() seems to give
      // wrong results. Why?
      p = Wind[sector].binscan2(conf);
      //p = Wind[sector].scan(conf);
      if(p == -100) continue;
      else sub2main[i*DB+j] = p;
      count++;
   }}
   if(count != Wind[sector].nBasis) std::cout<<"possible error in schmidtDecom "<<std::endl;
}


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


double schmidtDecom(std::vector<double> &vec, int sector, std::vector<MKL_INT> &sub2main){
  int i,j,p,k,l;
  int flagGI,count,sizet;
  double norm;
  sizet = Wind[sector].trans_sectors;
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
  //u       = (double*)malloc((M*M)*sizeof(double));
  //vt      = (double*)malloc((N*N)*sizeof(double));

  // count the frequency
  //for(i=0; i<sizet; i++){
  //   freq[i] = std::count(Wind[sector].Tflag.begin(), Wind[sector].Tflag.end(), i+1);
  //   //std::cout<<"Flag "<<i+1<<" occurs "<<freq[i]<<" times"<<std::endl;
  //}
  // this is counted and stored in Tbag[q]

  // initialize chi
  for(i=0; i<(DA*DB); i++){
    if(sub2main[i]==-5) chi[i]=0.0;
    else{
      k     = sub2main[i];                // get the ice state corresponding to the patch (eA,eB)
      l     = Wind[sector].Tflag[k]-1;    // get the translation bag for the ice state
      norm  = sqrt(Wind[sector].Tbag[l]); // get the frequency
      chi[i]= vec[l]/norm;
    }
  }

  // initialize chi
  //for(i=0;i<(DA*DB);i++) chi[i]=0.0;
  // calculate chi
  //count=0;
  //for(i=0; i<DA; i++){
  //for(j=0; j<DB; j++){
  //    cA = eA[i]; cB = eB[j];
  //    patch(conf,cA,cB);
      // if Gauss Law not satisfied, skip
  //    flagGI = checkGL2(conf);
  //    if(flagGI==0) continue;
  //    count++;
      // match with the corresponding basis state in the winding number sector
  //    p = Wind[sector].binscan2(conf);
  //    if(p == -100) continue;
      //std::cout<<i<<" "<<j<<" "<<p<<" "<<Wind[sector].Tflag[p]-1<<std::endl;
      // All ice states with the same Tflag get the same amplitude
  //    k = Wind[sector].Tflag[p]-1;
  //    chi[i*DB+j] = vec[k]/sqrt(freq[k]);
  //}}
  //std::cout<<"total no of GI states by patching "<<count<<std::endl;

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

void printvec(std::vector<double> &vec){
 int i;
 double norm;
 //std::cout<<"size ="<<vec.size()<<std::endl;
 norm = 0.0;
 for(i=0;i<vec.size();i++){
    norm += vec[i]*vec[i];
 //   std::cout<<vec[i]<<" ";
 }
 //std::cout<<std::endl;
 std::cout<<"Norm = "<<norm<<std::endl;
}
