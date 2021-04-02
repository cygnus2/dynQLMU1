#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<iterator>
#include "define.h"

void print2file(int, int, FILE*);
void storeState(FILE*, std::vector<bool>&);
void collectBasis(int, int*, std::vector<std::vector<bool>>&);
void diag_LAPACK_RRR2(int, std::vector<double>&, std::vector<double>&, std::vector<double>&);
// Translation in X,Y
void initTflag(int sector, std::vector<int>&, std::vector<int>&);
void TransX(std::vector<bool>&, std::vector<bool>&, int);
void TransY(std::vector<bool>&, std::vector<bool>&, int);
int whichIceState(std::vector<bool>&, std::vector<std::vector<bool>> &);
int whichIceState2(std::vector<bool>&, std::vector<std::vector<bool>> &);
int actOkin(int, std::vector<std::vector<bool>>&, std::vector<int>&);
// This routine set up the oKin matrix in the basis of ice states which have
// #-of-flippable states = VOL/2, and diagonalizes it. This is a fast "recipe" to
// construct scar states
void scarDiag(int sector){
    int nbScar, flag, trans_sectors;
    int nScarBags;
    int i, q, ix, iy;
    std::vector<std::vector<bool>> basisScar;
    // labelling the Translation flags and the degeneracy, and NFLIP
    std::vector<int> Tflag;
    std::vector<int> Tdgen;
    // arrays for translating and collecting translated states
    std::vector<bool> init(2*VOL);
    std::vector<bool> new1(2*VOL);
    std::vector<bool> new2(2*VOL);

    // obtain the basis states for the scar state
    nbScar = 0;
    collectBasis(sector, &nbScar, basisScar);
    std::cout<<"#-of basis states for the scars ="<<nbScar<<std::endl;
    std::cout<<"#-of basis states for the scars ="<<basisScar.size()<<std::endl;

    // initialize the translation flag
    initTflag(nbScar, Tflag, Tdgen);
    // calculate the translational bag basis; the procedure we follow is exactly
    // the same as in CODES3.
    flag=1;
    for(i=0; i<nbScar; i++){
      // when the flag is non-zero, this has already been considered
      if(Tflag[i]) continue;
      init = basisScar[i];
      // set flag for a unassigned state
      Tflag[i]=flag;
      // go through all possible translations, in x and in y
      for(iy=0;iy<LY;iy++){
      for(ix=0;ix<LX;ix++){
         TransX(init,new1,ix);
         TransY(new1,new2,iy);
         // check which basis state it corresponds to
         q = whichIceState2(new2, basisScar);
         // set flag of the new state
         if(q==-100){ printf("Error. Translated state not found. \n"); exit(0); }
         Tflag[q]=flag;
         // increase degeneracy
         Tdgen[q]++;
      }}
      flag++;
   }
   flag--;
   trans_sectors = flag;
   std::cout<<"Total number of translational sectors = "<<trans_sectors<<std::endl;
   //for(i=0; i<nbScar; i++) printf("state=%d, Tflag=%d\n",i,Tflag[i]);
   // at this point we have assigned bag indexes to all the states. We need to
   // check, for a representative state in a bag, if it gives rise to states with
   // unchanged number of flippable plaquettes. This function returns the result.
   nScarBags = actOkin(trans_sectors, basisScar, Tflag);

   // clear memory
   Tflag.clear(); Tdgen.clear(); basisScar.clear();
   init.clear();  new1.clear();  new2.clear();
}

void collectBasis(int sector, int *nbScar, std::vector<std::vector<bool>> &basisScar){
  int k,sizet;
  std::vector<bool> conf(2*VOL);

  sizet = Wind[sector].nBasis;
  /* compute #-basis states with k-flippable plaquettes */
  for(k=0; k<sizet; k++){
   if(Wind[sector].nflip[k] == VOL2){
     conf = Wind[sector].basisVec[k];
     basisScar.push_back(conf);
     (*nbScar)++;
   }
  }
  //std::sort(basisScar.begin(),basisScar.end());
}

void initTflag(int nbScar, std::vector<int> &Tflag, std::vector<int> &Tdgen){
    int i;
    for(i=0; i<nbScar; i++){
      Tflag.push_back(0);
      Tdgen.push_back(0);
    }
}

// Translate a basis state in the x-direction: state --> stateTx;  by lattice translations (lx,0)
// return the translated state
void TransX(std::vector<bool> &state,std::vector<bool> &stateT,int lx){
    int ix,iy,p,q;
    for(iy=0;iy<LY;iy++){
    for(ix=0;ix<LX;ix++){
       p = iy*LX + ix;         // initial point
       q = iy*LX + (ix+lx)%LX; // shifted by lx lattice units in +x dir
       stateT[2*q]=state[2*p]; stateT[2*q+1]=state[2*p+1];
    }}
}

// Translate a basis state in the y-direction: state --> stateTx;  by lattice translations (0,ly)
// return the translated state
void TransY(std::vector<bool> &state,std::vector<bool> &stateT,int ly){
    int ix,iy,p,q;
    for(iy=0;iy<LY;iy++){
    for(ix=0;ix<LX;ix++){
       p = iy*LX + ix;            // initial point
       q = ((iy+ly)%LY)*LX + ix;  // shifted by ly lattice units in +y dir
       stateT[2*q]=state[2*p]; stateT[2*q+1]=state[2*p+1];
    }}
}

int whichIceState(std::vector<bool> &newstate, std::vector<std::vector<bool>> &basisScar){
     unsigned int m;
     // binary search of the sorted array
     std::vector<std::vector<bool>>::iterator it;
     it = std::lower_bound(basisScar.begin(),basisScar.end(),newstate);
     m  = std::distance(basisScar.begin(),it);
     if(it == basisScar.end()){
       //std::cout<<"Element not found here! "<<std::endl;
       return -100;
     }
     return m;
}

int whichIceState2(std::vector<bool> &newstate, std::vector<std::vector<bool>> &basisScar){
   int nmax;
   nmax = basisScar.size();
   for(int m=0;m<nmax;m++){
     if(newstate == basisScar[m]) return m;
   }
   //std::cout<<"Element not found here!"<<std::endl;
   return -100;
 }

int actOkin(int trans_sectors, std::vector<std::vector<bool>> &basisScar, std::vector<int> &Tflag){
  int i,j,p,q,count;
  int p1,p2,p3,p4;
  bool pxy, pyz, pzw, pwx;
  FILE *fptr;
  // array to store states
  std::vector<bool> conf(2*VOL);
  std::vector<std::vector<bool>> repStates;

  // store representative states in an array
  for(i=0; i<trans_sectors; i++){
    for(j=0; j<basisScar.size(); j++){
      conf = basisScar[j];
      if(Tflag[j]==(i+1)){ repStates.push_back(conf); break; }
    }
  }
  printf("#-of-states stored=%ld\n",repStates.size());

  fptr = fopen("NflipVOL2.dat","w");
  count=0;
  for(i=0; i<repStates.size(); i++){
    // check if this state can be flipped without changing the #-flippable plaquettes
    conf = basisScar[i];
    /* act on the basis state with the Hamiltonian */
    /* a single plaquette is arranged as
              pzw
           o-------o
           |       |
      pwx  |   p   |  pyz
           |       |
           o-------o
              pxy
    */
    for(p=0;p<VOL;p++){
      // Find if a single plaquette is flippable
      p1=2*p; p2=2*next[DIM+1][p]+1; p3=2*next[DIM+2][p]; p4=2*p+1;
      pxy=conf[p1]; pyz=conf[p2]; pzw=conf[p3]; pwx=conf[p4];
      if((pxy==pyz)&&(pzw==pwx)&&(pwx!=pxy)){
       // If flippable, act with Okin
       conf[p1]=!conf[p1]; conf[p2]=!conf[p2];
       conf[p3]=!conf[p3]; conf[p4]=!conf[p4];
       // check for the new state
       q = whichIceState2(conf, basisScar);
       // flip back
       conf[p1]=!conf[p1]; conf[p2]=!conf[p2];
       conf[p3]=!conf[p3]; conf[p4]=!conf[p4];
       if(q >= 0){
         //printf("Match found for state = %d, at point = %d; q=%d\n",i,p,q);
         storeState(fptr, conf);
         count++;  break; }
      } // close if-statement
    } // close for loop over p VOL
    printf("Current value of count = %d\n",count);
  }
  fclose(fptr);
  printf("#-bags with at least one flippable plaquette that does not change nFlip= %d\n",count);
  // clear memory
  repStates.clear(); conf.clear();
  return count;
}

void storeState(FILE* fptr, std::vector<bool> &conf){
  int p;
  for(p=0;p<2*VOL;p++){
    fprintf(fptr,"%d ",(int)conf[p]);
  }
  fprintf(fptr,"\n");
}
