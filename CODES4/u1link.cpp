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
// QTOT = total charge either 1, or 2
// +ve charge is placed at (LXP,LYP), -ve charge at (LXM,LYM)
int QTOT,LXP,LYP,LXM,LYM;
int locQP, locQM;
double lam,Ti,Tf,dT;
int NTOT,NH;
WindNo Wind;
std::vector<std::vector<bool>> basis;
std::vector<std::vector<bool>> basis_nonflip;
std::vector<std::vector<bool>> basis_flip;
int CHKDIAG;

int main(){
  FILE *fptr;
  char string[50];
  int i,d,p,q;
  int x,y;
  int sector;
  extern void initneighbor(void);
  extern void conststates(void);
  extern void printbasis(void);
  extern int** allocateint2d(int, int);
  extern void deallocateint2d(int**,int,int);

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
  fscanf(fptr,"%s %d\n",string,&QTOT);
  fscanf(fptr,"%s %d\n",string,&LXP);
  fscanf(fptr,"%s %d\n",string,&LYP);
  fscanf(fptr,"%s %d\n",string,&LXM);
  fscanf(fptr,"%s %d\n",string,&LYM);
  fclose(fptr);
  if(( LX%2 != 0 )||( LY%2 !=0 )) { printf("Code does not work with odd LX and/or LY. \n"); exit(0); }
  if(LX<LY) printf("Please make sure LX >= LY. Unforseen errors can occur otherwise. \n");
  VOL  = LX*LY;
  VOL2 = VOL/2;

  // decide whether to check the results of the diagonalization
  CHKDIAG=1;

  /* Initialize nearest neighbours */
  for(i=0;i<=2*DIM;i++){
    next[i]    = (int *)malloc(VOL*sizeof(int));
    nextCHK[i] = (int *)malloc(VOL*sizeof(int));
  }

  /* Initialize checkerboard co-ordinates */
  lin2chk = (int *)malloc(VOL*sizeof(int));
  chk2lin = (int *)malloc(VOL*sizeof(int));
  initneighbor();

  // get the lattice co-ordinates for positive and negative charges
  locQP = lin2chk[LYP*LX + LXP];
  locQM = lin2chk[LYM*LX + LXM];
  if(locQP>=VOL2){ printf("+Q is not on the right sublattice \n"); exit(0); }
  if(locQM>=VOL2){ printf("-Q is not on the right sublattice \n"); exit(0); }
  printf("Q = %d, placed at (%d,%d); checkerboard site = %d\n", QTOT,LXP,LYP,locQP);
  printf("Q = %d, placed at (%d,%d); checkerboard site = %d\n",-QTOT,LXM,LYM,locQM);

  /* build basis states satisfying Gauss' Law */
  if(QTOT==1)       conststates_Q1();
  else if(QTOT==2)  conststates_Q2();
  else{  std::cout<<"Wrong choice \n"<<std::endl; exit(0); }

  //select the winding number (0,0)
  select_winding();

  constH();

  // calculate the <psi_n| O_flip |psi_n>, for every eigenstate psi_n
  calc_Oflip();

  // real time evolution of <PHI(t)| O_flip |PHI(t)> with initial states
  // with maximum number of flippable plaquettes
  calc_Oflipt();

  FilePrintPlaq();

  /* Clear memory */
  for(i=0;i<=2*DIM;i++){  free(next[i]); free(nextCHK[i]); }
  free(chk2lin); free(lin2chk);

  return 0;
}

void FilePrintPlaq(){
 FILE *fptr;
 int i,p;
 fptr=fopen("BASIS_QQBAR.dat","w");
 for(i=0;i<Wind.nBasis;i++){
   for(p=0;p<VOL;p++){
    if(Wind.xflip[i][p]) fprintf(fptr, " 1 ");
    else                 fprintf(fptr, " 0 ");
   }
   fprintf(fptr,"\n");
  }
 fclose(fptr);
}
