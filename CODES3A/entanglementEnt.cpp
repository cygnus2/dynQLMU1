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
#include<chrono>

extern void printconf1(std::vector<bool>);
extern void printconfA(std::vector<std::vector<bool>>);
extern void printconfB(std::vector<std::vector<bool>>);
extern void print_zmatrix( char*, MKL_INT, MKL_INT, MKL_Complex16*, MKL_INT );
extern void print_rmatrix( char* , MKL_INT , MKL_INT , double* , MKL_INT  );
extern void patch(std::vector<bool>&, std::vector<bool>&, std::vector<bool>&);
extern int checkGL2(std::vector<bool>&);
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void createLookupTable(int, MKL_INT, MKL_INT, std::vector<MKL_INT>&);
extern double schmidtDecom(std::vector<double>&, int, std::vector<MKL_INT>&, size_t);

// variables used in this set of functions
unsigned int LEN_A,LEN_B,VOL_A,VOL_B;
unsigned int DA,DB,NCHI;
// basis for subsystems A and B
std::vector<std::vector<bool>> eA;
std::vector<std::vector<bool>> eB;
// entanglement Entropy for the eigenvector
double EE;

void entanglementEnt_INIT0(int sector){
  int p,i;
  int tsect;
  FILE *outf1, *outf2;
  // eigenvector and eigenvalue to calculate EE
  double sel_eval00, sel_evalPiPi;
  std::vector<double> sel_evec00, sel_evecPiPi;
  // overlap with the initial state
  std::vector<double> alpha00, alphaPiPi;
  tsect = Wind[sector].trans_sectors;
  // this stores which subsystem states map to states of the full system
  std::vector<MKL_INT> sub2main;

  // Construct the spin basis for the sub-systems
  LEN_B = LX - LEN_A;
  VOL_A = LEN_A*LY; VOL_B = LEN_B*LY;
  createBasis(sector);
  sub2main.assign(DA*DB, -5); // negative initial value
  // Get starting timepoint
  auto start = std::chrono::high_resolution_clock::now();
  createLookupTable(sector,DA,DB,sub2main);
  // Get ending timepoint
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time taken for lookup table="<<duration.count()<< " secs"<<std::endl;

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]*inorm);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]*inorm);
  }

  outf1 = fopen("EntE00.dat","w");
  outf2 = fopen("EntEPiPi.dat","w");
  fprintf(outf1,"# Entropy of Entanglement as a function of the eigenvalues for LA = %d\n",LEN_A);
  fprintf(outf1,"# Eigenvalues  Entanglement Entropy  Overlap_w_InitC\n");
  fprintf(outf2,"# Entropy of Entanglement as a function of the eigenvalues for LA = %d\n",LEN_A);
  fprintf(outf2,"# Eigenvalues  Entanglement Entropy  Overlap_w_InitC\n");

  // Calculate the EE for each of the eigenstates
  for(p=0; p<tsect; p++){
    // initialize the state and the entropy
    sel_evec00.clear(); sel_evecPiPi.clear(); EE = -100.0;
    sel_eval00  = Wind[sector].evals_K00[p];
    sel_evalPiPi= Wind[sector].evals_KPiPi[p];
    //sel_evec.assign(Wind[sector].evecs_K00[p*sizet], Wind[sector].evecs_K00[(p+1)*sizet-1]);
    for(i=0; i < tsect; i++){
       sel_evec00.push_back(Wind[sector].evecs_K00[p*tsect + i]);
       sel_evecPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect + i]);
    }
    //std::cout<<"Going to do Schmidt decompose eigenvector = "<< p << std::endl;
    //std::cout<<"Size of the vector = "<<int(sel_evec.size())<<std::endl;
    //printvec(sel_evec);
    EE = schmidtDecom(sel_evec00,sector,sub2main,0);
    fprintf(outf1,"%.12lf %.12lf %.12lf\n",sel_eval00,EE,alpha00[p]*alpha00[p]);
    EE = schmidtDecom(sel_evecPiPi,sector,sub2main,1);
    fprintf(outf2,"%.12lf %.12lf %.12lf\n",sel_evalPiPi,EE,alphaPiPi[p]*alphaPiPi[p]);
  }
  fclose(outf1);
  fclose(outf2);
  // free memory from the spin basis
  eA.clear(); eB.clear();
}

void entanglementEnt_INIT4(int sector){
  int p,i,tsect;
  int chkPi0, chk0Pi;
  FILE *outf1, *outf2, *outf3, *outf4;
  // eigenvector and eigenvalue to calculate EE
  double sel_eval00, sel_evalPiPi, sel_evalPi0, sel_eval0Pi;
  std::vector<double> sel_evec00, sel_evecPiPi;
  std::vector<double> sel_evecPi0, sel_evec0Pi;
  // overlap with the initial state
  std::vector<double> alpha00, alphaPiPi;
  std::vector<double> alphaPi0, alpha0Pi;
  tsect = Wind[sector].trans_sectors;
  // this stores which subsystem states map to states of the full system
  std::vector<MKL_INT> sub2main;

  // Construct the spin basis for the sub-systems
  LEN_B = LX - LEN_A;
  VOL_A = LEN_A*LY; VOL_B = LEN_B*LY;
  createBasis(sector);
  sub2main.assign(DA*DB, -5); // negative initial value
  // Get starting timepoint
  auto start = std::chrono::high_resolution_clock::now();
  createLookupTable(sector,DA,DB,sub2main);
  // Get ending timepoint
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time taken for lookup table="<<duration.count()<< " secs"<<std::endl;

  // < w_k | IN >; k-th eigenvector; IN=initial state; details about initial state
  for(p=0; p<tsect; p++){
       alpha00.push_back(Wind[sector].evecs_K00[p*tsect+INITbag]*inorm);
       alphaPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect+INITbag]*inorm);
       alphaPi0.push_back(Wind[sector].evecs_KPi0[p*tsect+INITbag]*inorm);
       alpha0Pi.push_back(Wind[sector].evecs_K0Pi[p*tsect+INITbag]*inorm);
  }

  outf1 = fopen("EntE00.dat","w");
  outf2 = fopen("EntEPiPi.dat","w");
  outf3 = fopen("EntEPi0.dat","w");
  outf4 = fopen("EntE0Pi.dat","w");
  fprintf(outf1,"# Entropy of Entanglement as a function of the eigenvalues for LA = %d\n",LEN_A);
  fprintf(outf1,"# Eigenvalues  Entanglement Entropy  Overlap_w_InitC\n");
  fprintf(outf2,"# Entropy of Entanglement as a function of the eigenvalues for LA = %d\n",LEN_A);
  fprintf(outf2,"# Eigenvalues  Entanglement Entropy  Overlap_w_InitC\n");
  fprintf(outf3,"# Entropy of Entanglement as a function of the eigenvalues for LA = %d\n",LEN_A);
  fprintf(outf3,"# Eigenvalues  Entanglement Entropy  Overlap_w_InitC\n");
  fprintf(outf4,"# Entropy of Entanglement as a function of the eigenvalues for LA = %d\n",LEN_A);
  fprintf(outf4,"# Eigenvalues  Entanglement Entropy  Overlap_w_InitC\n");

  // initialize variables needed to track the spurious eigenstates
  chkPi0 = 0; chk0Pi = 0;
  // Get starting timepoint
  start = std::chrono::high_resolution_clock::now();
  // Calculate the EE for each of the eigenstates
  for(p=0; p<tsect; p++){
    // initialize the state and the entropy
    sel_evec00.clear();  sel_evecPiPi.clear(); EE = -100.0;
    sel_evecPi0.clear(); sel_evec0Pi.clear();
    sel_eval00  = Wind[sector].evals_K00[p];
    sel_evalPiPi= Wind[sector].evals_KPiPi[p];
    sel_evalPi0 = Wind[sector].evals_KPi0[p];
    sel_eval0Pi = Wind[sector].evals_K0Pi[p];
    //sel_evec.assign(Wind[sector].evecs_K00[p*sizet], Wind[sector].evecs_K00[(p+1)*sizet-1]);
    for(i=0; i < tsect; i++){
       sel_evec00.push_back(Wind[sector].evecs_K00[p*tsect + i]);
       sel_evecPiPi.push_back(Wind[sector].evecs_KPiPi[p*tsect + i]);
       sel_evecPi0.push_back(Wind[sector].evecs_KPi0[p*tsect + i]);
       sel_evec0Pi.push_back(Wind[sector].evecs_K0Pi[p*tsect + i]);
    }
    //std::cout<<"Going to do Schmidt decompose eigenvector = "<< p << std::endl;
    //std::cout<<"Size of the vector = "<<int(sel_evec.size())<<std::endl;
    //printvec(sel_evec);
    EE = schmidtDecom(sel_evec00,sector,sub2main,0);
    fprintf(outf1,"%.12lf %.12lf %.12lf\n",sel_eval00,EE,alpha00[p]*alpha00[p]);
    EE = schmidtDecom(sel_evecPiPi,sector,sub2main,1);
    fprintf(outf2,"%.12lf %.12lf %.12lf\n",sel_evalPiPi,EE,alphaPiPi[p]*alphaPiPi[p]);
    if(spurPi0[chkPi0] == p) chkPi0++;
    else{
     EE = schmidtDecom(sel_evecPi0,sector,sub2main,2);
     fprintf(outf3,"%.12lf %.12lf %.12lf\n",sel_evalPi0,EE,alphaPi0[p]*alphaPi0[p]);
    }
    if(spur0Pi[chk0Pi] == p) chk0Pi++;
    else{
      EE = schmidtDecom(sel_evec0Pi,sector,sub2main,3);
      fprintf(outf4,"%.12lf %.12lf %.12lf\n",sel_eval0Pi,EE,alpha0Pi[p]*alpha0Pi[p]);
    }
  }
  // Get ending timepoint
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time taken for EE calculation="<<duration.count()<< " secs"<<std::endl;

  fclose(outf1);
  fclose(outf2);
  fclose(outf3);
  fclose(outf4);
  // free memory from the spin basis
  eA.clear();  eB.clear();
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


double schmidtDecom(std::vector<double> &vec, int sector, std::vector<MKL_INT> &sub2main, size_t cs){
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

  if(cs==0){
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
  }
  else if(cs==1){
    // initialize chi
    for(i=0; i<(DA*DB); i++){
      if(sub2main[i]==-5) chi[i]=0.0;
      else{
        k     = sub2main[i];                // get the ice state corresponding to the patch (eA,eB)
        l     = Wind[sector].Tflag[k]-1;    // get the translation bag for the ice state
        norm  = sqrt(Wind[sector].Tbag[l]); // get the frequency
        chi[i]= Wind[sector].FPiPi[k]*vec[l]/norm;
      }
    }
  }
  else if(cs==2){
    // initialize chi
    for(i=0; i<(DA*DB); i++){
      if(sub2main[i]==-5) chi[i]=0.0;
      else{
        k     = sub2main[i];                // get the ice state corresponding to the patch (eA,eB)
        l     = Wind[sector].Tflag[k]-1;    // get the translation bag for the ice state
        norm  = sqrt(Wind[sector].Tbag[l]); // get the frequency
        chi[i]= Wind[sector].FPi0[k]*vec[l]/norm;
      }
    }
  }
  else if(cs==3){
    // initialize chi
    for(i=0; i<(DA*DB); i++){
      if(sub2main[i]==-5) chi[i]=0.0;
      else{
        k     = sub2main[i];                // get the ice state corresponding to the patch (eA,eB)
        l     = Wind[sector].Tflag[k]-1;    // get the translation bag for the ice state
        norm  = sqrt(Wind[sector].Tbag[l]); // get the frequency
        chi[i]= Wind[sector].F0Pi[k]*vec[l]/norm;
      }
    }
  }

  // check norm
  //norm = 0.0;
  //for(i=0;i<(DA*DB);i++) norm += chi[i]*chi[i];
  //if( fabs(norm-1.0) >  1e-6) std::cout<<"Norm = "<<norm<<std::endl;

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
