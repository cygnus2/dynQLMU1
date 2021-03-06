/* compute the expectation values of Ey(x) */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

//extern void initState(int, int, int*);

// Notation: eigenstate |n> = \sum_k \alpha_k |k>, |k> is a cartoon state
void evolveH_ov3(int sector){
    extern void printconf(bool*);
    std::vector<double> initC;
    int i,ix,iy,p,q;
    int m,k,j,r;
    double t;
    double betaR,betaI,betaM,betaTot;
    double EyProf[LX],dEyProf[LX];
    double Cf0, Cf1;
    double avgd1, avgd2, avgh1, avgh2, avgv1, avgv2;
    int sizet, nTot;
    FILE *fptr,*fptr1,*fptr2;

    /* find the relevant cartoon states with the specified number of flippable plaquettes */
    if((LX == 2) && (LY == 2)){
      std::cout<<"This routine does not include the case LX=2 and LY=2"<<std::endl; exit(0);}
    if(LY > 2){ std::cout<<"This routine does not work for LY>2"<<std::endl; exit(0); }

    sizet = Wind[sector].nBasis;

    /* obtain the initial starting state */
    if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
    q = INITq;
    std::cout<<"Starting state is basis state = "<<q<<std::endl;
    /* store the overlap of the initial state with the eigenvectors */
    for(p=0; p<sizet; p++){
        initC.push_back(Wind[sector].evecs[p*sizet+q]);
    }

    /* now compute the overlap in each sector */
    fptr = fopen("EyProf.dat","w");
    fptr1= fopen("dEyProf.dat","w");
    fptr2= fopen("CorrF.dat","w");
    for(t=Ti; t<Tf; t=t+dT){

      /* initialize Ey profile */
      for(k=0;k<LX;k++) { EyProf[k]=0.0; dEyProf[k]=0.0; }
      /* initialize the correlation functions */
      avgd1=0.0; avgd2=0.0; avgh1=0.0; avgh2=0.0; avgv1=0.0; avgv2=0.0;

      for(k=0; k<sizet; k++){
          betaR = 0.0; betaI = 0.0;
          for(m=0; m<sizet; m++){
       betaR += Wind[sector].evecs[m*sizet+k]*initC[m]*cos(Wind[sector].evals[m]*t);
	     betaI += Wind[sector].evecs[m*sizet+k]*initC[m]*sin(Wind[sector].evals[m]*t);
	        }
          betaM = betaR*betaR + betaI*betaI;
          /* get the Ey profile at time t */
          for(r=0;r<LX;r++){
	  	      EyProf[r] += Wind[sector].Ey[k][r]*betaM;
            dEyProf[r]+= Wind[sector].dEy[k][r]*betaM;
	        }
          /* correlation functions */
          avgd1 += Wind[sector].OOd1[k]*betaM; avgd2 += Wind[sector].OOd2[k]*betaM;
          avgv1 += Wind[sector].OOv1[k]*betaM; avgv2 += Wind[sector].OOv2[k]*betaM;
          avgh1 += Wind[sector].OOh1[k]*betaM; avgh2 += Wind[sector].OOh2[k]*betaM;
      } // close loop over basis states
      /* print the Ey profile at each times */
      fprintf(fptr,"%lf ",t);
      fprintf(fptr1,"%lf ",t);
      for(r=0;r<LX;r++){ fprintf(fptr,"%lf ",EyProf[r]); fprintf(fptr1,"%lf ",dEyProf[r]); }
      fprintf(fptr,"\n"); fprintf(fptr1,"\n");
      fprintf(fptr2,"%.4lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf\n",t,
                avgd1,avgd2,avgv1,avgv2,avgh1,avgh2);
    }
    fclose(fptr);
    fclose(fptr1);
    fclose(fptr2);

    /* free memory */
    initC.clear();
}
