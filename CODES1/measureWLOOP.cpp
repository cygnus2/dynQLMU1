#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#include<armadillo>

#include "define.hxx"
#include "routines.h"

double getloopval(int m,int n, int ***state,int tot){
  int i,j,p,nsp,flag,pp;
  int countx,county,sx;
  double **loop;
  double genergy,gstate[N1],loopavg;
  double dumen,dumstate[N1];
  int newstate[DIM][VOL];
  FILE *EVL,*EVC;

  /* Read in the ground state */
  EVL=fopen("evals.dat","r");
  EVC=fopen("evecs.dat","r");
  for(i=0;i<N1;i++){
    fscanf(EVL,"%d %le",&sx,&dumen);
    for(j=0;j<N1;j++){
     fscanf(EVC,"%le",&dumstate[j]);
    }
    if(i==0) {
      genergy=dumen;
      for(j=0;j<N1;j++)
        gstate[j]=dumstate[j];
    }
  }
  fclose(EVL);
  fclose(EVC);


  loopavg=0.0;

  /* construct the Wilson loop operator */
  
  loop=allocatedouble2d(tot,tot);
  for(i=0;i<tot;i++){
  for(p=0;p<VOL;p++){
    countx=0;county=0;
    for(j=0;j<VOL;j++){
       newstate[0][j] = state[i][0][j];
       newstate[1][j] = state[i][1][j];
    }
    /* Act the Wilson loop operator */
    flag=1;
    nsp=p;
    /* m-units in the x-direction */
    while(countx<m){
     if(newstate[0][nsp]==-1) newstate[0][nsp]=1;
     else {newstate[0][nsp]=0; flag=0; }
     countx++; 
     nsp  = next[DIM+1][nsp];
    }
    /* n-units in the y-direction */
    while(county<n){
     if(newstate[1][nsp]==-1) newstate[1][nsp]=1;
     else {newstate[1][nsp]=0; flag=0; }
     county++;
     nsp = next[DIM+2][nsp]; 
    }
    /* m-units backward in x-direction */
    while(countx){
     if(newstate[0][next[DIM-1][nsp]]==1) newstate[0][next[DIM-1][nsp]]=-1;
     else { newstate[0][next[DIM-1][nsp]]=0; flag=0; }
     countx--;
     nsp  = next[DIM-1][nsp]; 
    }
    /* n-units bakward in the y-direction */
    while(county){
     if(newstate[1][next[DIM-2][nsp]]==1) newstate[1][next[DIM-2][nsp]]=-1;
     else { newstate[1][next[DIM-2][nsp]]=0; flag=0; }
     county--;
     nsp = next[DIM-2][nsp];
    }
    if(p!=nsp) { printf("Error\n"); exit(0);}
    if (flag==1){
    /* compare which state it is */
    pp=scan(newstate,state,tot);
    loop[i][pp] += 1.0;
    }
  }
  }
  /* Calculate the expectation value in the ground state. 
   * Note that can again be restricted to be in the subspace of the
   * non-trivial states since the ground state is in this subspace
   */
  for(i=0;i<N1;i++){
  for(j=0;j<N1;j++){
    loop[i][j]/=((double)VOL);
    loopavg += gstate[i]*gstate[j]*loop[i][j];
  }
  }
 
  deallocatedouble2d(loop,tot,tot);
  return loopavg;
}
 
void measureWLOOP(int N1,int N2){
  int N = N1 + N2;
  int ***state;
  double LOOP[LX-1][LY-1];
  FILE *spin;
  int sx,sy,i,j;

  /* Get all states */
  state=allocate3d(N,DIM,VOL);
  spin = fopen("SPINSTATES","r");
  for(i=0;i<N1;i++){
  for(j=0;j<VOL;j++){
    fscanf(spin,"%d %d",&sx,&sy);
    state[i][0][j]=sx;
    state[i][1][j]=sy;
  }
  }
  fclose(spin);
  /* Since the ground state lies in the subspace of non-trivial flux states, we need not worry
   * about the trivial part
   */
  spin = fopen("TSPINSTATES","r");
  for(i=N1;i<N;i++){
  for(j=0;j<VOL;j++){
    fscanf(spin,"%d %d",&sx,&sy);
    state[i][0][j]=sx;
    state[i][1][j]=sy;
  }
  }
  fclose(spin);
  /* calc loop */
  for(i=0;i<LX-1;i++){
  for(j=0;j<LY-1;j++){
     LOOP[i][j] = getloopval(i+1,j+1,state,N);
     printf("LOOP: (Nx=%d,Ny=%d) = % .5e\n",i+1,j+1,LOOP[i][j]);
  }
  }

  deallocate3d(state,N,DIM,VOL);
}
