#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

 
 /* Decompose the Winding Number sector into states connected by translation flags */
 void trans_decompose(int sector){
   int i,flag,count;
   int ix,iy,q;
   std::vector<bool> init(2*VOL);
   std::vector<bool> new1(2*VOL);
   std::vector<bool> new2(2*VOL);

   // initialize the translation flag
   Wind[sector].initTflag();
   std::cout<<"No of basis states in this sector = "<<Wind[sector].nBasis<<std::endl;

   // Note that with the calculation of the degeneracy one has to be careful. 
   // Either start the loop from (Tx,Ty)=(0,0) and then only increase the 
   // degeneracy in the loop. Otherwise, increase the Tgen every time you find
   // an unmarked state, and then do translations starting from tx=1,.. and ty=1,... 
   flag=1;
   for(i=0;i<Wind[sector].nBasis;i++){
      // when the flag is non-zero, this has already been considered
      if(Wind[sector].Tflag[i]) continue;
      init = Wind[sector].basisVec[i];
      // set flag for a unassigned state
      Wind[sector].Tflag[i]=flag;
      // increment the degeneracy
      Wind[sector].Tdgen[i]++;
      // go through all possible translations, in x and in y
      for(iy=1;iy<LY;iy++){
      for(ix=1;ix<LX;ix++){
         Wind[sector].WindNo::TransX(init,new1,ix);   
         Wind[sector].WindNo::TransY(new1,new2,iy);   
         // check which basis state it corresponds to
         //q = Wind[sector].scan(new2);
         q = Wind[sector].binscan(new2);
         // set flag of the new state
         if(q==-100){ printf("Error. Translated state not in Winding no sector. \n"); exit(0); }
         Wind[sector].Tflag[q]=flag;
         // increase degeneracy
         Wind[sector].Tdgen[q]++;
      }}
      flag++;
   }
   flag--;
   std::cout<<"current flag ="<<flag<<std::endl;
   Wind[sector].disp_Tprop();
 }

 // the inidividual basis (ice) states in the Winding number sectors are initially labelled with 0
 // the first bag of states labelled by lattice translations is labelled as 1, the second one as 2, and so on. 
 // the total number of bags is the highest counter. 
 // the degeneracy is initialized with zero.
 void WindNo::initTflag(){
    int i;
    for(i=0;i<nBasis;i++){
      Tflag.push_back(0); 
      Tdgen.push_back(0);
    }
 }
 
 // Translate a basis state in a given winding number in the x-direction 
 // state --> stateTx;  by lattice translations (lx,0) and return the translated state
 void WindNo::TransX(std::vector<bool> &state,std::vector<bool> &stateT,int lx){
    int ix,iy,p,q;
    for(iy=0;iy<LY;iy++){
    for(ix=0;ix<LX;ix++){
       p = iy*LX + ix;         // initial point
       q = iy*LX + (ix+lx)%LX; // shifted by lx lattice units in +x dir
       stateT[2*q]=state[2*p]; stateT[2*q+1]=state[2*p+1];  
    }}
 }

 // Translate a basis state in a given winding number in the y-direction 
 // state --> stateTx;  by lattice translations (0,ly) and return the translated state
 void WindNo::TransY(std::vector<bool> &state,std::vector<bool> &stateT,int ly){
    int ix,iy,p,q;
    for(iy=0;iy<LY;iy++){
    for(ix=0;ix<LX;ix++){
       p = iy*LX + ix;            // initial point
       q = ((iy+ly)%LY)*LX + ix;  // shifted by ly lattice units in +y dir
       stateT[2*q]=state[2*p]; stateT[2*q+1]=state[2*p+1];  
    }}
 }

 int WindNo::binscan(std::vector<bool> &newstate){
     unsigned int m;
     // binary search of the sorted array  
     std::vector<std::vector<bool>>::iterator it; 
     it = std::lower_bound(basisVec.begin(),basisVec.end(),newstate);
     m  = std::distance(basisVec.begin(),it);
     //if(it == basisVec.end()){
     //  std::cout<<"Element not found here! "<<std::endl; 
     //  return -100;
     //}
     return m;
 }

 void WindNo::disp_Tprop(){
   int i;
   for(i=0;i<nBasis;i++){
     std::cout<<"flag for basis state "<<i<<" = "<<Tflag[i]<<" Degeneracy ="<<Tdgen[i]<<std::endl;
   }
 }
