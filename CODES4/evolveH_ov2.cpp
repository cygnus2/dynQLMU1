/* compute spectral weights in each flippable plaquette sector */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void evolveH_ov2(int sector){
    extern void printconf(bool*);
    std::vector<double> initC;
    std::vector<bool> cart(2*VOL);    
    int i,ix,iy,p,q,parity;
    int m,k,j,r;
    int nFlip[VOL+1]; 
    int ch1, ch2; // allow the choice of a different cartoon initial
    int cx1, cx2; // state, i.e., with domain walls placed elsewhere
    int p1,p2,p3,p4;
    bool s1,s2,s3,s4;
    double specWT[VOL+1];
    double t;
    double betaR,betaI,betaM,betaTot;
    double fprof[VOL];
    int sizet, nTot;
    FILE *fptr,*fptr1; 

    /* find the relevant cartoon states with the specified number of flippable plaquettes */
    if((LX == 2) && (LY == 2)){
      std::cout<<"This routine does not include the case LX=2 and LY=2"<<std::endl; exit(0);}
    if(LY > 2){ std::cout<<"This routine does not work for LY>2"<<std::endl; exit(0); }
    
    sizet = Wind[sector].nBasis;

    /* initialize */
    for(k=0;k<=VOL;k++) nFlip[k]=0; 
    
    /* compute #-basis states with k-flippable plaquettes */
    for(k=0; k<sizet; k++){
	  p = Wind[sector].nflip[k];
	  nFlip[p]++;
    }

    q=0;
    /* define the initial cartoon state */
    if(INIT==0){   // symmetry broken cartoon state
      for(iy=0;iy<LY;iy++){
      for(ix=0;ix<LX;ix++){
         parity=(ix+iy)%2;
         p = 2*(iy*LX+ix);
         if(parity){ cart[p]=false; cart[p+1]=true;  }
         else      { cart[p]=true;  cart[p+1]=false; }
      }}
      q=Wind[sector].binscan(cart); 
      std::cout<<"Starting state is basis state = "<<q<<std::endl;
      /* store the overlap of the initial state with the eigenvectors */
      for(p=0; p<sizet; p++){
        initC.push_back(Wind[sector].evecs[p*sizet+q]);
      }
    }
    else if(INIT==1){ // 1-plaquette flipped cartoon state
       ch1 = 2; cx1 = 0;
       for(k=0; k<sizet; k++){
	    if(Wind[sector].nflip[k]==(VOL-3)){
                if(cx1 == ch1){ q = k; break; }
		cx1++;
	    }   
       }
       std::cout<<"Starting state is basis state = "<<q<<std::endl;
      /* store the overlap of the initial state with the eigenvectors */
      for(p=0; p<sizet; p++){
        initC.push_back(Wind[sector].evecs[p*sizet+q]);
      }
    }
    else if(INIT==2){ // 2 domain wall cartoon state, even separation
       r = (LX/2);
       if(r%2==1) r=r-1;
       for(k=0; k<sizet; k++){
	    if(Wind[sector].nflip[k]==(VOL-4)){
		    ix=0;   iy=0; p1=iy*LX+ix;  s1=Wind[sector].xflip[k][p1];
		    ix=1;   iy=1; p2=iy*LX+ix;  s2=Wind[sector].xflip[k][p2];
		    ix=r;   iy=0; p3=iy*LX+ix;  s3=Wind[sector].xflip[k][p3];
		    ix=r+1; iy=1; p4=iy*LX+ix;  s4=Wind[sector].xflip[k][p4];
		    if((!s1) && (!s2) && (!s3) && (!s4)){ q = k; break;}
	    }  
       }
       std::cout<<"Starting state is basis state = "<<q<<std::endl;
      /* store the overlap of the initial state with the eigenvectors */
      for(p=0; p<sizet; p++){
        initC.push_back(Wind[sector].evecs[p*sizet+q]);
      }
    }
    else if(INIT==3){ // 2 domain wall cartoon state, odd separation
       r = (LX/2);
       if(r%2==1) r=r+1;       
       for(k=0; k<sizet; k++){
	    if(Wind[sector].nflip[k]==(VOL-4)){
		    ix=0;     iy=0; p1=iy*LX+ix;  s1=Wind[sector].xflip[k][p1];
		    ix=1;     iy=1; p2=iy*LX+ix;  s2=Wind[sector].xflip[k][p2];
		    ix=r-1;   iy=1; p3=iy*LX+ix;  s3=Wind[sector].xflip[k][p3];
		    ix=r;     iy=0; p4=iy*LX+ix;  s4=Wind[sector].xflip[k][p4];
		    if( (!s1) && (!s2) && (!s3) && (!s4)){ q = k; break;}
	    }   
       }
       if(k==sizet) std::cout<<"maximal state reached "<<std::endl;
       std::cout<<"Starting state is basis state = "<<q<<std::endl;
      /* store the overlap of the initial state with the eigenvectors */
      for(p=0; p<sizet; p++){
        initC.push_back(Wind[sector].evecs[p*sizet+q]);
      }
    }

    /* now compute the overlap in each sector */
    fptr = fopen("spectral_WT.dat","w");
    fptr1= fopen("state_Prof.dat","w");
    for(t=Ti; t<Tf; t=t+dT){

      /* initialize spectral weights */
      for(k=VOL;k>=0;k--) specWT[k]=0.0;

      /* initialize flux profile */
      for(k=0;k<VOL;k++) fprof[k]=0.0; 

      for(k=0; k<sizet; k++){
          p = Wind[sector].nflip[k];
          betaR = 0.0; betaI = 0.0;
          for(m=0; m<sizet; m++){
	     betaR += Wind[sector].evecs[m*sizet+k]*initC[m]*cos(-Wind[sector].evals[m]*t);	  
	     betaI += Wind[sector].evecs[m*sizet+k]*initC[m]*sin(-Wind[sector].evals[m]*t);	  
	  }
          betaM = betaR*betaR + betaI*betaI;
          specWT[p] = specWT[p] + betaM;
          
          /* get the flippability profile at time t */
          for(r=0;r<VOL;r++){
	  	  if(Wind[sector].xflip[k][r]) fprof[r] += betaM; 	  
	  }
      } // close loop over basis states
      betaTot = 0.0;
      /* print the spectral weights, starting from maximal flippable plaquettes,
       * increasing with the number of non-flippable plaquettes */
      fprintf(fptr,"%lf ",t);
      for(k=VOL; k>=0; k--){
	  fprintf(fptr,"%lf ",specWT[k]);
          betaTot += specWT[k];	  
      }
      fprintf(fptr,"%lf \n",betaTot);
      /* print the flippability profile at each times */
      fprintf(fptr1,"%lf ",t);
      for(k=0;k<VOL;k++){ fprintf(fptr1,"%lf ",fprof[k]); }
       fprintf(fptr1,"\n");
    }
    fclose(fptr);
    fclose(fptr1);

    /* free memory */
    initC.clear();
}

