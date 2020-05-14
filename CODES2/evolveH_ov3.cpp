/* compute the expectation values of Ey(x) */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

extern void initState(int, int, int*);

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void evolveH_ov3(int sector){
    extern void printconf(bool*);
    std::vector<double> initC;
    int i,ix,iy,p,q;
    int m,k,j,r;
    double t;
    double betaR,betaI,betaM,betaTot;
    double EyProf[LX];
    int sizet, nTot;
    FILE *fptr;

    /* find the relevant cartoon states with the specified number of flippable plaquettes */
    if((LX == 2) && (LY == 2)){
      std::cout<<"This routine does not include the case LX=2 and LY=2"<<std::endl; exit(0);}
    if(LY > 2){ std::cout<<"This routine does not work for LY>2"<<std::endl; exit(0); }

    sizet = Wind[sector].nBasis;

    /* obtain the initial starting state */
    q=0;
    initState(sector,INIT,&q);
    if(q<0){
      std::cout<<"Initial state not found!"<<std::endl; exit(0); }
    else{
      std::cout<<"Starting state is basis state = "<<q<<std::endl;
      /* store the overlap of the initial state with the eigenvectors */
      for(p=0; p<sizet; p++){
        initC.push_back(Wind[sector].evecs[p*sizet+q]);
      }
    }

    /* now compute the overlap in each sector */
    fptr = fopen("EyProf.dat","w");
    for(t=Ti; t<Tf; t=t+dT){

      /* initialize Ey profile */
      for(k=0;k<LX;k++) EyProf[k]=0.0;

      for(k=0; k<sizet; k++){
          betaR = 0.0; betaI = 0.0;
          for(m=0; m<sizet; m++){
       betaR += Wind[sector].evecs[m*sizet+k]*initC[m]*cos(-Wind[sector].evals[m]*t);
	     betaI += Wind[sector].evecs[m*sizet+k]*initC[m]*sin(-Wind[sector].evals[m]*t);
	    }
          betaM = betaR*betaR + betaI*betaI;
          /* get the Ey profile at time t */
          for(r=0;r<LX;r++){
	  	      EyProf[r] += Wind[sector].Ey[k][r]*betaM;
	        }
      } // close loop over basis states
      /* print the Ey profile at each times */
      fprintf(fptr,"%lf ",t);
      for(k=0;k<LX;k++){ fprintf(fptr,"%lf ",EyProf[k]); }
      fprintf(fptr,"\n");
    }
    fclose(fptr);

    /* free memory */
    initC.clear();
}
