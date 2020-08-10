#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>


/* according to the rules of cpp, the variables are declared here
 * and also in the header file as extern such that they are avl to
 * the other functions.
 */
int *next[2*DIM+1];
int *nextCHK[2*DIM+1];
int *chk2lin,*lin2chk;
int LX,LY,VOL,VOL2;
// Labels the winding number sectors
// lookup[LX/2+wx][LY/2+wy] refers to the (wx,wy) winding number sector
int nWindSector;
int **lookup;
std::vector<WindNo> Wind;
unsigned int nWind;
double lam,Ti,Tf,dT;
int NTOT,NH;
std::vector<std::vector<bool>> basis;
std::vector<std::vector<bool>> basis_nonflip;
std::vector<std::vector<bool>> basis_flip;
int CHKDIAG, STORE_SVD;
int INIT, INITq, INITbag;
int INITphasePi0, INITphase0Pi, INITphasePiPi;
double inorm;
int CALC;
// if CALC=0 calculate everything
// CALC==1,2,3,4 compute OflipT, state_Prof, Ey(dEy), EENT respectively

int main(){
  FILE *fptr;
  char string[50];
  int i,d,p,q;
  int x,y;
  int wx,wy;
  int sector;
  //WindNo SectorZero;
  extern void initneighbor(void);
  extern void conststates(void);
  extern void printbasis(void);
  extern int** allocateint2d(int, int);
  extern void deallocateint2d(int**,int,int);
  extern int calc_WindNo(int,int);

  fptr = fopen("QUEUE","r");
  if(fptr == NULL){
      printf("could not open QUEUE FILE to open\n");
      exit(1);
  }
  fscanf(fptr,"%s %d\n",string,&LX);
  fscanf(fptr,"%s %d\n",string,&LY);
  fscanf(fptr,"%s %lf\n",string,&lam);
  fscanf(fptr,"%s %lf\n",string,&Ti);
  fscanf(fptr,"%s %lf\n",string,&Tf);
  fscanf(fptr,"%s %lf\n",string,&dT);
  fscanf(fptr,"%s %d\n",string,&LEN_A);
  fscanf(fptr,"%s %d\n",string,&INIT);
  fscanf(fptr,"%s %d\n",string,&CALC);
  fclose(fptr);
  if(( LX%2 != 0 )||( LY%2 !=0 )) { printf("Code does not work with odd LX and/or LY. \n"); exit(0); }
  if(LX<LY) printf("Please make sure LX >= LY. Unforseen errors can occur otherwise. \n");
  VOL = LX*LY;
  VOL2 = VOL/2;

  std::cout<<"Initial state ="<<INIT<<std::endl;
  std::cout<<"lambda ="<<lam<<std::endl;

  // decide whether to check the results of the diagonalization
  CHKDIAG=0;

  /* Initialize nearest neighbours */
  for(i=0;i<=2*DIM;i++){
    next[i] = (int *)malloc(VOL*sizeof(int));
    nextCHK[i] = (int *)malloc(VOL*sizeof(int));
  }

  /* Initialize checkerboard co-ordinates */
  lin2chk = (int *)malloc(VOL*sizeof(int));
  chk2lin = (int *)malloc(VOL*sizeof(int));
  initneighbor();

  /* Winding number sectors */
  lookup = allocateint2d(LX+1,LY+1);

  /* build basis states satisfying Gauss' Law */
  conststates();

  /* get number of winding number sectors */
  nWind = calc_WindNo(LX,LY);
  Wind.reserve(nWind);
  winding_no_decompose();

  // get the winding number sector (wx,wy)
  wx = 0; wy = 0;
  sector = lookup[LX/2+wx][LY/2+wy];
  constH(sector);

  /* breakup into translation sectors */
  trans_decompose(sector);

  /* If INIT==0, Hamiltonian in the (kx,ky)=(0,0) & (Pi,Pi) sectors are constructed */
  /* If INIT==4, H in sectors (0,0), (Pi,Pi), (Pi,0), (0,Pi) are constructed        */
  if(INIT==0)      trans_Hamil_INIT0(sector);
  else if(INIT==4) trans_Hamil_INIT4(sector);

  INITq = -1;
  // initialize the starting state (once and for all the routines)
  initState(sector, INIT, &INITq);
  // initalize the bag details for the initial state
  INITbag = Wind[sector].Tflag[INITq]-1;
  inorm = sqrt(Wind[sector].Tdgen[INITq]/(double)VOL);
  INITphasePi0 = Wind[sector].FPi0[INITq];
  INITphase0Pi = Wind[sector].F0Pi[INITq];
  INITphasePiPi= Wind[sector].FPiPi[INITq];
  std::cout<<"Initial state belongs to bag ="<<INITbag<<std::endl;
  std::cout<<"Normalization ="<<inorm<<std::endl;
  std::cout<<"Phases: (Pi,0):"<<INITphasePi0<<";  (0,Pi):"<<INITphase0Pi<<"; (Pi,Pi):"<<INITphasePiPi<<std::endl;

  // real time evolution of <PHI(t)| O_flip |PHI(t)>
  // starting from specified initial states in each sector (see notes)
  if(CALC==0 || CALC==1)  calc_Oflipt(sector);

  // real-time evolution state_Prof
  if(CALC==0 || CALC==2){
     if(INIT==0)      evolveH_ov2_INIT0(sector);
     else if(INIT==4) evolveH_ov2_INIT4(sector);
  }

  // real-time evolution dEy (for INIT=0) or Ey (for INIT=4)
  if(CALC==0 || CALC==3){
     if(INIT==0)      evolveH_ov3_INIT0(sector);
     else if(INIT==4) evolveH_ov3_INIT4(sector);
  }

  // calculate the Entanglement Entropy for the states
  if(CALC==0 || CALC==4){
     if(INIT==0)      entanglementEnt_INIT0(sector);
     else if(INIT==4) entanglementEnt_INIT4(sector);
  }

  // calculate the Fidelity
  if(CALC==0 || CALC==5){
     if(INIT==0)      Lecho_INIT0(sector);
     else if(INIT==4) Lecho_INIT4(sector);
  }

  /* Clear memory */
  for(i=0;i<=2*DIM;i++){  free(next[i]); free(nextCHK[i]); }
  free(chk2lin); free(lin2chk);
  deallocateint2d(lookup,LX+1,LY+1);
  Wind.clear();

  return 0;
}
