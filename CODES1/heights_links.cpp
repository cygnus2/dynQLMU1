// In this code, we count all the configurations of the 4x4 lattice
// which satisfy the Gauss' Law constraint for Q=0, +2, -2 (with a
// zero charge neutrality condition. The counting is done both 
// with the link and the height configurations

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>
#define DIM 2

/* according to the rules of cpp, the variables are declared here 
 * and also in the header file as extern such that they are avl to
 * the other functions.
 */
int *next[2*DIM+1];
int *nextCHK[2*DIM+1];
int *chk2lin,*lin2chk;
int LX,LY,VOL,VOL2;
double lam,Ti,Tf,dT;
long int NTOT;

int main(){
  int i,d,p,q;
  int x,y;
  extern void initneighbor(void);
  extern long int constQ02(void);
  
  /* Lattice parameters */ 
  LX=4; LY=4;
  VOL  = LX*LY; 
  VOL2 = 2*VOL;

  /* Initialize nearest neighbours */
  for(i=0;i<=2*DIM;i++){
    next[i] = (int *)malloc(VOL*sizeof(int));
    nextCHK[i] = (int *)malloc(VOL*sizeof(int));
  }

  /* Initialize checkerboard co-ordinates */
  lin2chk = (int *)malloc(VOL*sizeof(int));
  chk2lin = (int *)malloc(VOL*sizeof(int));
  initneighbor();

  NTOT = constQ02();
  std::cout<<"Total states returned ="<<NTOT<<std::endl;

  /* Clear memory */
  for(i=0;i<=2*DIM;i++){  free(next[i]); free(nextCHK[i]); }
  free(chk2lin); free(lin2chk);
  return 0;
}

long int constQ02(){
   extern void Dec2Bin(long int, std::vector<bool>&);
   long int NST,count;
   std::vector<bool> conf;
   conf.resize(VOL2);

   NST=pow(2,VOL2);
   count=0;
   while(count<NST){
      Dec2Bin(count,conf);
      count++;
   }
   return count;
}


// Function to convert a given decimal number to a binary array
void Dec2Bin(long int inputNum, std::vector<bool> &state){
    int index,i;
    index=0;
    //printf("%d  ",inputNum);
    while (inputNum > 0){
        if(inputNum%2) state[index]=true; 
        else state[index]=false;
        inputNum = inputNum/2;
        index++;
    }
    while(index < 2*VOL){
       state[index]=false; index++;
    }
    //printf("index = %d;  ",index);
    //for(unsigned int k=2*VOL;k>0;k--) std::cout<<state[k-1];
    //printf("\n");
}


void initneighbor(void){
  int p,x,y,q;
  int lin;

  /* lattice linear co-ordinates */
  for(p=0;p<VOL;p++){
      x = (p%LX);
      y = (p/LX)%LY;

      next[DIM][p] = p;
      next[DIM+1][p] = y*LX + ((x+1)%LX);
      next[DIM+2][p] = ((y+1)%LY)*LX + x;
      next[DIM-1][p] = y*LX + ((x-1+LX)%LX);
      next[DIM-2][p] = ((y-1+LY)%LY)*LX + x;
  }

  /* checkerboard to linear and vice-versa */
  p=0;
  for(y=0;y<LY;y++){
  for(x=0;x<LX;x++){
    if((x+y)%2==0){
       q = y*LX + x;
       chk2lin[p] = q; lin2chk[q] = p;
       printf("checkerboard site = %d(%d,%d); Lin site = %d\n",p,x,y,q);
       p++;
    }
  }}
  for(y=0;y<LY;y++){
  for(x=0;x<LX;x++){
    if((x+y)%2==1){
       q = y*LX + x;
       chk2lin[p]= q; lin2chk[q] = p;
       printf("checkerboard site = %d(%d,%d); Lin site = %d\n",p,x,y,q);
       p++;
    }
  }}
  // nextCHK targets the neighboring sites in the string 
  // representing the GL implementation on the even sites
  for(p=0;p<VOL;p++){
      lin = chk2lin[p];
      nextCHK[DIM][p]   = p;
      nextCHK[DIM+1][p] = lin2chk[next[DIM+1][lin]];
      nextCHK[DIM-1][p] = lin2chk[next[DIM-1][lin]];
      nextCHK[DIM+2][p] = lin2chk[next[DIM+2][lin]];
      nextCHK[DIM-2][p] = lin2chk[next[DIM-2][lin]];
      //printf("checkerboard NN. site = %d (in chkrboard = %d)\n",lin,p);
      //printf("+x=%d,   -x=%d,   +y=%d,   -y=%d\n",nextCHK[DIM+1][p],
      //nextCHK[DIM-1][p],nextCHK[DIM+2][p],nextCHK[DIM-2][p]); 
   }
}

