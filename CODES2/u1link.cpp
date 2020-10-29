#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>
#include<chrono>
#include<ctime>

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
int STORE_SVD;
int CHKDIAG;
// initial state is specified by W and LR and INIT
// W is the amount of y-flux=Wy.
// LR=0 flux to the left (subsystem LA), LR=1 flux to the right (subsystem LB).
// INIT=0, symm broken states with all flippable plaquettes
// INIT=1, C1-C2-C1 domain wall states with fused domain walls
// INIT=2,3 domain wall states with inter-dw distance ~ LX/2
// INIT=4, domain walls and anti-domain walls arranged as alternating manner
// INIT=5, domain walls and anti-domain walls are phase separated.
// For winding number states INIT has a value 10+WX.
int WX, WY, LR;
int INIT, INITq;

int main(){
  FILE *fptr;
  char string[50];
  int i,d,p,q;
  int x,y;
  int sector;
  WindNo SectorZero;
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
  fscanf(fptr,"%s %d\n",string,&WX);
  fscanf(fptr,"%s %d\n",string,&WY);
  fscanf(fptr,"%s %d\n",string,&LR);
  fscanf(fptr,"%s %d\n",string,&INIT);
  fclose(fptr);
  if(( LX%2 != 0 )||( LY%2 !=0 )) { printf("Code does not work with odd LX and/or LY. \n"); exit(0); }
  if(LX<LY) printf("Please make sure LX >= LY. Unforseen errors can occur otherwise. \n");
  VOL = LX*LY;
  VOL2 = VOL/2;

  if((WX!=0)&&((INIT)!=(10+std::abs(WX)))){
     std::cout<<"Please specify INIT = 10 + abs(WX)"<<std::endl; exit(0);
  }

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
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  conststates();
  std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
  std::cout << "Time needed " << elapsed.count() << " milliseconds." << std::endl;

  /* get number of winding number sectors */
  nWind = calc_WindNo(LX,LY);
  Wind.reserve(nWind);
  winding_no_decompose();

  printf("Chosen (Wx,Wy) sector = (%d,%d)\n",WX,WY);
  sector = lookup[LX/2+WX][LY/2+WY];
  constH(sector);

  INITq = -1;
  //set the starting state once (and for all!)
  initState(sector, INIT, &INITq);

  // calculate the <psi_n| O_flip |psi_n>, <psi_n|CEy1|psi_n> for every eigenstate psi_n
  calc_Oflip(sector);

  // real time evolution of <PHI(t)| O_flip |PHI(t)> and <PHI(t)| Cflip |PHI(t)>
  // starting from specified initial states in each sector (see notes)
  // note that recalculates the same as the previous routine, so don't use both!
  //calc_Oflipt(sector);

  // real time evolution of <PHI(t)| O_kin |PHI(t)>,
  //calc_Okint(sector);

  //FilePrintBasis(sector);

  // Lochschmidt Echo
  //Lecho(sector);

  // real-time evolution of initial states (evolveH_ov2
  // has all the functionalities of evolveH_ov1 built in!)
  //// evolveH_ov1(sector);
  //evolveH_ov2(sector);
  //evolveH_ov3(sector);

  // Entanglement Entropy calculations
  evolve_Eent(sector);

  // Investigate potential "scar" states
  if(lam == 0.0) studyEvecs2(sector);
  else           studyEvecs(sector);

  // Set up ED in the basis of states with flippable plaq = VOL/2
  //scarDiag(sector);

  /* Clear memory */
  for(i=0;i<=2*DIM;i++){  free(next[i]); free(nextCHK[i]); }
  free(chk2lin); free(lin2chk);
  deallocateint2d(lookup,LX+1,LY+1);
  Wind.clear();

  return 0;
}
