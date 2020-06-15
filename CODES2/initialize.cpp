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

extern void printconf(std::vector<bool>);

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
      std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
      return;
    } // close INIT=0
    else if(INIT==1){ // 1-plaquette flipped cartoon state
       ch1 = 5; cx1 = 0;
       for(k=0; k<sizet; k++){
            if(Wind[sector].nflip[k]==(VOL-3)){
                if(cx1 == ch1){ *q = k;
			std::cout<<"In routine initState. Starting state="<<(*q)<<std::endl;
      std::cout<<"#-of-flippable plaqs ="<<Wind[sector].nflip[(*q)]<<std::endl;
      std::cout<<"Printing the flippability profile of initial state "<<std::endl;
      for(r=0; r<VOL; r++) std::cout<<r<<"  "<<Wind[sector].xflip[(*q)][r]<<std::endl;
      std::cout<<"Printing the Ey profile of the initial state:"<<std::endl;
      for(r=0;r<LX;r++) std::cout<<r<<"  "<<Wind[sector].Ey[*q][r]<<" "<<Wind[sector].dEy[*q][r]<<std::endl;
      std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
			return; }
                cx1++;
            }
       }
    } // close INIT=1
    else if(INIT==2){ // 2 domain wall cartoon state, even separation
       r = (LX/2);
       if(r%2==1) r=r-1;
       for(k=0; k<sizet; k++){
            if(Wind[sector].nflip[k]==(VOL-4)){
                    ix=0;   iy=0; p1=iy*LX+ix;  s1=Wind[sector].xflip[k][p1];
                    ix=1;   iy=1; p2=iy*LX+ix;  s2=Wind[sector].xflip[k][p2];
                    ix=r;   iy=0; p3=iy*LX+ix;  s3=Wind[sector].xflip[k][p3];
                    ix=r+1; iy=1; p4=iy*LX+ix;  s4=Wind[sector].xflip[k][p4];
                    if((!s1) && (!s2) && (!s3) && (!s4)){ *q = k;
                        std::cout<<"In routine initState. Starting state ="<<(*q)<<std::endl;
      std::cout<<"#-of-flippable plaqs ="<<Wind[sector].nflip[(*q)]<<std::endl;
      std::cout<<"Printing the flippability profile of initial state "<<std::endl;
      for(r=0; r<VOL; r++) std::cout<<r<<"  "<<Wind[sector].xflip[(*q)][r]<<std::endl;
      std::cout<<"Printing the Ey profile of the initial state:"<<std::endl;
      for(r=0;r<LX;r++) std::cout<<r<<"  "<<Wind[sector].Ey[*q][r]<<" "<<Wind[sector].dEy[*q][r]<<std::endl;
      std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
			return;}
            }
       }
    } // close INIT=2
    else if(INIT==3){ // 2 domain wall cartoon state, odd separation
       r = (LX/2);
       if(r%2==1) r=r+1;
       for(k=0; k<sizet; k++){
            if(Wind[sector].nflip[k]==(VOL-4)){
                    ix=0;     iy=0; p1=iy*LX+ix;  s1=Wind[sector].xflip[k][p1];
                    ix=1;     iy=1; p2=iy*LX+ix;  s2=Wind[sector].xflip[k][p2];
                    ix=r-1;   iy=1; p3=iy*LX+ix;  s3=Wind[sector].xflip[k][p3];
                    ix=r;     iy=0; p4=iy*LX+ix;  s4=Wind[sector].xflip[k][p4];
                    if((!s1) && (!s2) && (!s3) && (!s4)){ *q = k;
		  std::cout<<"In routine initState. Starting state = "<<(*q)<<std::endl;
      std::cout<<"#-of-flippable plaqs ="<<Wind[sector].nflip[(*q)]<<std::endl;
      std::cout<<"Printing the flippability profile of initial state "<<std::endl;
      for(r=0; r<VOL; r++) std::cout<<r<<"  "<<Wind[sector].xflip[(*q)][r]<<std::endl;
      std::cout<<"Printing the Ey profile of the initial state:"<<std::endl;
      for(r=0;r<LX;r++) std::cout<<r<<"  "<<Wind[sector].Ey[*q][r]<<" "<<Wind[sector].dEy[*q][r]<<std::endl;
      std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
		       return;}
            }
       }
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
       std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
       return;
    } // close INIT=4
    else if(INIT==5){ // domain walls together, and anti-domain walls stacked together
       for(ix=0;ix<(LX/2);ix++){
           for(iy=0;iy<LY;iy++){
              p = 2*(iy*LX+ix);
	            if(iy%2)   cart[p]=false;
	            else       cart[p]=true;
              py= p+1;
	            cart[py]=true;
           }
       }
       for(ix=(LX/2);ix<LX;ix++){
           for(iy=0;iy<LY;iy++){
              p = 2*(iy*LX+ix);
	            if(iy%2)   cart[p]=false;
	            else       cart[p]=true;
              py= p+1;
	            cart[py]=false;
           }
       }
       /* find the relevant basis state in the hilbert space */
       //*q = Wind[sector].binscan2(cart);
       *q = Wind[sector].scan(cart);
       std::cout<<"The initial state is "<<(*q)<<std::endl;
       std::cout<<"#-of-flippable plaqs ="<<Wind[sector].nflip[(*q)]<<std::endl;
       std::cout<<"Printing the Ey profile of the initial state:"<<std::endl;
       for(r=0;r<LX;r++) std::cout<<r<<"  "<<Wind[sector].Ey[*q][r]<<" "<<Wind[sector].dEy[*q][r]<<std::endl;
       std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
       std::cout<<"Printing the flippability profile of initial state "<<std::endl;
       for(r=0; r<VOL; r++) std::cout<<r<<"  "<<Wind[sector].xflip[(*q)][r]<<std::endl;
       return;
    } // close INIT=5
    else if(INIT==6){ // Ey =+1,-1 scattered randomly in x
       /* initialize random seed: */
       srand (time(NULL));
       // there is a (low) probability that the chosen Ey fluxes are all next to each other
       // need to insert a condition to stop this. Check with the #-of-flippable plaquettes
       do{
          for(r=0;r<LX;r++) posX[r]=-1; // initialize everything to -1
          k=0;  // counter to check how many distinct values are initialized
          while(k<(LX/2)){
               r=rand()%LX;
               if(posX[r]==1) continue;
               posX[r]=1;  k++;
          };
          printf("Randomly selected profile\n");
          for(k=0;k<LX;k++) std::cout<<"Ey["<<k<<"]= "<<posX[k]<<std::endl;
          // initialize the x-links
          for(ix=0;ix<LX;ix++){
          for(iy=0;iy<LY;iy++){
                p=2*(iy*LX+ix);
                if(iy%2)   cart[p]=false;
 	              else       cart[p]=true;
          }}
          // initialize the y-links
          for(ix=0;ix<LX;ix++){
          for(iy=0;iy<LY;iy++){
                p=2*(iy*LX+ix)+1;
                if(posX[ix]==1) cart[p]=true;
                else            cart[p]=false;
          }}
          /* find the relevant basis state in the hilbert space */
          *q = Wind[sector].scan(cart);
       }while(Wind[sector].nflip[*q]<=2);
       std::cout<<"The initial state is "<<(*q)<<std::endl;
       std::cout<<"#-of-flippable plaqs ="<<Wind[sector].nflip[(*q)]<<std::endl;
       std::cout<<"Printing the Ey profile of the initial state:"<<std::endl;
       for(r=0;r<LX;r++) std::cout<<r<<"  "<<Wind[sector].Ey[*q][r]<<" "<<Wind[sector].dEy[*q][r]<<std::endl;
       std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;
       std::cout<<"Printing the flippability profile of initial state "<<std::endl;
       for(r=0; r<VOL; r++) std::cout<<r<<"  "<<Wind[sector].xflip[(*q)][r]<<std::endl;
       return;
    } // close INIT=6
  } // close if((wx==0)&&(wy==0))
  else if((W>0)&&(wy==0)&&(INIT>10)){
      std::cout<<"Initial state in the winding sector "<<std::endl;
      if(wx>0) sign=1;
      else if(wx<0) sign=-1;
      else{ printf("Logic error. \n"); exit(0);}
      /* cartoon state in (Wx,Wy)=(WX,0) */
      if(LR==0){
        // the flux is in the right end (LB)
        //std::cout<<"Putting strings on the right end"<<std::endl;
        for(iy=0;iy<LY;iy++){
           for(ix=0;ix<(LX-2*W);ix++){
             parity=(ix+iy)%2;
             p = 2*(iy*LX+ix);
             if(parity){ cart[p]=true;  cart[p+1]=false; }
             else{       cart[p]=false; cart[p+1]=true; }
           }
           for(ix=(LX-2*W);ix<LX;ix++){
             parity=(iy)%2;
             p = 2*(iy*LX+ix);
             if(sign==1){
                 if(parity){ cart[p]=false; cart[p+1]=true; }
                 else{       cart[p]=true;  cart[p+1]=true; }
             }
             else if(sign==-1){
                 if(parity){ cart[p]=false; cart[p+1]=false; }
                 else{       cart[p]=true;  cart[p+1]=false; }
             }
           } // close loop for ix
        }    // close loop for iy
      } // close if(LR==0)
      else if(LR==1){
        // flux is in the left end (LA)
        //std::cout<<"Putting strings on the left end"<<std::endl;
        for(iy=0;iy<LY;iy++){
           for(ix=2*W;ix<LX;ix++){
              parity=(ix+iy)%2;
              p = 2*(iy*LX+ix);
              if(parity){ cart[p]=true;  cart[p+1]=false; }
              else{       cart[p]=false; cart[p+1]=true; }
           }
           for(ix=0;ix<2*W;ix++){
              parity=(iy)%2;
              p = 2*(iy*LX+ix);
              if(sign==1){
                    if(parity){ cart[p]=false; cart[p+1]=true; }
                    else{       cart[p]=true;  cart[p+1]=true; }
              }
              else if(sign==-1){
                    if(parity){ cart[p]=false; cart[p+1]=false; }
                    else{       cart[p]=true;  cart[p+1]=false; }
              }
           } // close loop ix
        }    // close loop iy
      }// close if(LR==1)
  /* print the cartoon state for explicit checking */
  //printconf(cart);
  /* find the relevant basis state in the hilbert space */
  //*q = Wind[sector].binscan2(cart);
  *q = Wind[sector].scan(cart);
  std::cout<<"The initial state is "<<(*q)<<std::endl;
  std::cout<<"#-of-flippable plaqs ="<<Wind[sector].nflip[(*q)]<<std::endl;
  std::cout<<"Printing the flippability profile of initial state "<<std::endl;
  for(r=0; r<VOL; r++) std::cout<<r<<"  "<<Wind[sector].xflip[(*q)][r]<<std::endl;
  std::cout<<"Printing the Ey profile of the initial state:"<<std::endl;
  for(r=0;r<LX;r++) std::cout<<r<<"  "<<Wind[sector].Ey[*q][r]<<" "<<Wind[sector].dEy[*q][r]<<std::endl;
  std::cout<<"CEy1 ="<<Wind[sector].CEy1[*q]<<"; CEy2 ="<<Wind[sector].CEy2[*q]<<";"<<std::endl;


  // sanity checks DELETE LATER!
  chk1=Wind[sector].scan(cart);
  chk2=Wind[sector].binscan2(cart);
  if(chk1!=(*q)) std::cout<<"Mismatch binscan vs scan. q="<<(*q)<<" chk1="<<chk1<<std::endl;
  if(chk2!=(*q)) std::cout<<"Mismatch binscan vs binscan2. q="<<(*q)<<" chk2="<<chk2<<std::endl;

  /* print the state in the stored list of basis states */
  //testc = Wind[sector].basisVec[(*q)];
  //printconf(testc);
  return;
  } // close else if((W>0)&&(wy==0)&&(INIT>10))
  // control comes here if not encountered the cases before
  if(*q==-100){
     std::cout<<"Correct initial state not selected. Exiting. "<<std::endl;
     exit(0);
  }
}
