 /* routine to initialize the starting state */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>
#include<time.h>

void printconf(std::vector<bool>);
bool k00, kPiPi, kPi0, k0Pi;

/* The construction of the initial states follows the idea that
 * the cartoon states are identified with the already stored
 * basis states. The q-th basis (handled as a pointer) is sent
 * back into the main routine.  */
void initState(int sector, int INIT, int *q){
  int ix,iy,parity,k,p,r,py;
  int W,sign,wx,wy,sizet;
  int ch1,ch2; // allow the choice of a different initial state,
  int cx1,cx2; // i.e, with domain walls placed elsewhere
  int flag;
  int p1,p2,p3,p4,s1,s2,s3,s4;
  std::vector<int> posX(LX); // co-ordinates to set Ey[x]=+1,-1
  std::vector<bool> cart(2*VOL);

  int chk1,chk2; //to check the routines scan, binscan and binscan2.

  wx = Wind[sector].Wx;
  wy = Wind[sector].Wy;
  W = std::abs(wx);
  if(wy!=0){
	  printf("The option with Wy != 0 is not yet supported. \n"); exit(0);
  }
  std::cout<<"Routine to initialize states. INIT="<<INIT<<", wx="<<wx<<", wy="<<wy<<std::endl;
  // k00, kPiPi, kPi0, k0Pi indicate if the initial states have overlap in these sectors
  k00=false; kPiPi=false; kPi0=false; k0Pi=false;
  if(INIT==0){ k00=true; kPiPi=true; }
  else if(INIT==4){ k00=true; kPiPi=true; kPi0=true; k0Pi=true; }
  else{
    printf("No other initial state is coded in yet. \n"); exit(0);
  }

  /* initialization */
  *q = -100;
  sizet = Wind[sector].nBasis;
  /* initial states in (Wx,Wy)=(0,0) */
  if((wx==0)&&(wy==0)&&(INIT<=10)){
    if(INIT==0){ // symmetry broken cartoon state; all flippable plaqs
      for(iy=0;iy<LY;iy++){
      for(ix=0;ix<LX;ix++){
         parity=(ix+iy)%2;
         p = 2*(iy*LX+ix);
         if(parity){ cart[p]=false; cart[p+1]=true;  }
         else      { cart[p]=true;  cart[p+1]=false; }
      }}
      // locate the cartoon state in the list of basis states
      (*q) = Wind[sector].binscan(cart);
      std::cout<<"In routine initState. Starting state ="<<(*q)<<std::endl;
      std::cout<<"#-of-flippable plaqs ="<<Wind[sector].nflip[(*q)]<<std::endl;
      std::cout<<"Printing the flippability profile of initial state "<<std::endl;
      for(r=0; r<VOL; r++) std::cout<<r<<"  "<<Wind[sector].xflip[(*q)][r]<<std::endl;
      std::cout<<"Printing the Ey profile of the initial state:"<<std::endl;
      for(r=0;r<LX;r++) std::cout<<r<<"  "<<Wind[sector].Ey[*q][r]<<" "<<Wind[sector].dEy[*q][r]<<std::endl;
      //std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
      return;
    } // close INIT=0
    else if(INIT==1){ // 1-plaquette flipped cartoon state
      std::cout<<"The initial state is not compatible with the momentum sectors studied"<<std::endl;
      exit(0);
    } // close INIT=1
    else if(INIT==2){ // 2 domain wall cartoon state, even separation
      std::cout<<"The initial state is not compatible with the momentum sectors studied"<<std::endl;
      exit(0);
    } // close INIT=2
    else if(INIT==3){ // 2 domain wall cartoon state, odd separation
      std::cout<<"The initial state is not compatible with the momentum sectors studied"<<std::endl;
      exit(0);
    } // close INIT=3
    else if(INIT==4){ // domain walls cartoon state, staggered along x-direction
       for(ix=0;ix<LX;ix++){
	         parity = ix%2;
           for(iy=0;iy<LY;iy++){
              p = 2*(iy*LX+ix);
	            if(iy%2)   cart[p]=false;
	            else       cart[p]=true;
              py= p+1;
	            if(parity) cart[py]=false;
	            else       cart[py]=true;
         }
         sign=-sign;
       }
       /* find the relevant basis state in the hilbert space */
       //*q = Wind[sector].binscan2(cart);
       *q = Wind[sector].scan(cart);
       std::cout<<"The initial state is "<<(*q)<<std::endl;
       std::cout<<"#-of-flippable plaqs ="<<Wind[sector].nflip[(*q)]<<std::endl;
       std::cout<<"Printing the flippability profile of initial state "<<std::endl;
       for(r=0; r<VOL; r++) std::cout<<r<<"  "<<Wind[sector].xflip[(*q)][r]<<std::endl;
       std::cout<<"Printing the Ey profile of the initial state:"<<std::endl;
       for(r=0;r<LX;r++) std::cout<<r<<"  "<<Wind[sector].Ey[*q][r]<<" "<<Wind[sector].dEy[*q][r]<<std::endl;
       //std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
       return;
    } // close INIT=4
    else if(INIT==5){ // domain walls together, and anti-domain walls stacked together
      std::cout<<"The initial state is not compatible with the momentum sectors studied"<<std::endl;
      exit(0);
    } // close INIT=5
    else if(INIT==6){ // Ey =+1,-1 scattered randomly in x
      std::cout<<"The initial state is not compatible with the momentum sectors studied"<<std::endl;
      exit(0);
    } // close INIT=6
  } // close if((wx==0)&&(wy==0))
  else if((W>0)&&(wy==0)&&(INIT>10)){
      std::cout<<"Supply initial state in a given momentum sector "<<std::endl;
      exit(0);
  } // close else if((W>0)&&(wy==0)&&(INIT>10))
  // control comes here if not encountered the cases before
  if(*q==-100){
     std::cout<<"Correct initial state not selected. Exiting. "<<std::endl;
     exit(0);
  }
}
