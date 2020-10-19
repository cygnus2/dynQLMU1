#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<iterator>
#include "define.h"

void print2file(int, int, FILE*);
void collectBasis(int, int*, std::vector<std::vector<bool>>&);
void diag_LAPACK_RRR2(int, std::vector<double>&, std::vector<double>&, std::vector<double>&);
// Translation in X,Y
void initTflag(int sector, std::vector<int>&, std::vector<int>&);
void TransX(std::vector<bool>&, std::vector<bool>&, int);
void TransY(std::vector<bool>&, std::vector<bool>&, int);
int whichIceState(std::vector<bool>&, std::vector<std::vector<bool>> &);
// This routine set up the oKin matrix in the basis of ice states which have
// #-of-flippable states = VOL/2, and diagonalizes it. This is a fast "recipe" to
// construct scar states
void scarDiag(int sector){
    int nbScar, flag, trans_sectors;
    int i, q, ix, iy;
    std::vector<std::vector<bool>> basisScar;
    // labelling the Translation flags and the degeneracy
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
         q = whichIceState(new2, basisScar);
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
       std::cout<<"Element not found here! "<<std::endl;
       return -100;
     }
     return m;
}
