/* evolve states in real-time with the Hamiltonian */
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
// In this function, we compute the overlap of the evolved state with the
// ones of the ground state manifold (2 states) and the first excited state
// manifold (2*LX*LY) states. This is used to check the comparison with
// perturbation theory.
void evolveH_ov1(int sector){
    extern void printconf(bool*);
    int i,ix,iy,p,q;
    int m,k,j;
    double t;
    int sizet, nstates0, nstates1, nstates2, nstates3, nTot;
    int	FLIP1, FLIP2, FLIP3, FLIP4;
    /* list contains the basis number of the cartoon states with which
     * we want the overlap; alphaK is the init state; while alphaL are
     * the states of the first manifold with three non-flippable plaq,
     * and alphaM are the ones of the second manifold with four non-
     * flippable plaqs; Eq(11)(??) of the dynU1.pdf */
    std::vector<double> alphaK,alphaL,alphaM,alphaN;
    std::vector<double> initC;
    double betaR,betaI,betaM,betaTot;
    FILE *fptr;
    bool bconf[2*VOL];

    /* find the relevant cartoon states with the specified number of flippable plaquettes */
    if((LX == 2) && (LY == 2)){
       std::cout<<"This routine does not include the case LX=2 and LY=2"<<std::endl; exit(0);}
    if(LY > 2){ std::cout<<"This routine does not work for LY>2"<<std::endl; exit(0); }

    sizet = Wind[sector].nBasis;

    //initialize the starting state
    if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
    q = INITq;
    std::cout<<"Starting state is basis state = "<<q<<std::endl;
    /* store the overlap of the initial state with the eigenvectors */
    for(p=0; p<sizet; p++){
        initC.push_back(Wind[sector].evecs[p*sizet+q]);
    }

    FLIP1 = LX*LY;    nstates0=0;
    FLIP2 = LX*LY-3;  nstates1=0;
    FLIP3 = LX*LY-4;  nstates2=0;
    FLIP4 = LX*LY/2;  nstates3=0;

    for(i=0; i<sizet; i++){
        /* overlap of states with FLIP1 flippable plaqs with eigenvectors */
        if(Wind[sector].nflip[i] == FLIP1){ nstates0++;
        for(p=0; p<sizet; p++){
           alphaK.push_back(Wind[sector].evecs[p*sizet+i]);
         } }
        /* overlap of states with FLIP2 flippable plaqs with eigenvectors */
        if(Wind[sector].nflip[i] == FLIP2){ nstates1++;
        for(p=0; p<sizet; p++){
           alphaL.push_back(Wind[sector].evecs[p*sizet+i]);
        } }
        /* overlap of states with FLIP3 flippable plaqs with eigenvectors */
        if(Wind[sector].nflip[i] == FLIP3){ nstates2++;
        for(p=0; p<sizet; p++){
           alphaM.push_back(Wind[sector].evecs[p*sizet+i]);
        } }
        /* overlap of states with FLIP4 flippable plaqs with eigenvectors */
        if(Wind[sector].nflip[i] == FLIP4){ nstates3++;
        for(p=0; p<sizet; p++){
           alphaN.push_back(Wind[sector].evecs[p*sizet+i]);
        } }
    }
    std::cout<<"#-basis states with max flippable plaq ="<<nstates0<<std::endl;
    std::cout<<"#-basis states in manifold 1 ="<<nstates1<<std::endl;
    std::cout<<"#-basis states in manifold 2="<<nstates2<<std::endl;
    std::cout<<"#-basis states in manifold 3="<<nstates3<<std::endl;
    nTot = nstates0 + nstates1 + nstates2 + nstates3;
    std::cout<<"Computed total #-basis states ="<<nTot<<std::endl;

    /* overlap in the reduced "classical" basis */
    // Note that the spectral weight in each class of plaquettes is
    // summed over.
    fptr = fopen("overlap_RedBas.dat","w");
    for(t=Ti; t<Tf; t=t+dT){
      fprintf(fptr, "%lf ",t);
      betaTot = 0.0;
      for(k=0; k<nstates0; k++){
	betaR = 0.0; betaI = 0.0;
        for(m=0; m<sizet; m++){
             betaR += alphaK[k*sizet+m]*initC[m]*cos(Wind[sector].evals[m]*t);
             betaI += alphaK[k*sizet+m]*initC[m]*sin(Wind[sector].evals[m]*t);
        }
	betaM = betaR*betaR + betaI*betaI;
        //fprintf(fptr,"%lf ",betaM);
	betaTot += betaM;
      }
      fprintf(fptr,"%lf ",betaTot);
      for(k=0; k<nstates1; k++){
	betaR = 0.0; betaI = 0.0;
        for(m=0; m<sizet; m++){
             betaR += alphaL[k*sizet+m]*initC[m]*cos(Wind[sector].evals[m]*t);
             betaI += alphaL[k*sizet+m]*initC[m]*sin(Wind[sector].evals[m]*t);
        }
        betaM = betaR*betaR + betaI*betaI;
        //fprintf(fptr,"%lf ",betaM);
        betaTot += betaM;
      }
      fprintf(fptr,"%lf ",betaTot);
      for(k=0; k<nstates2; k++){
	betaR = 0.0; betaI = 0.0;
        for(m=0; m<sizet; m++){
             betaR += alphaM[k*sizet+m]*initC[m]*cos(Wind[sector].evals[m]*t);
             betaI += alphaM[k*sizet+m]*initC[m]*sin(Wind[sector].evals[m]*t);
        }
        betaM = betaR*betaR + betaI*betaI;
        //fprintf(fptr,"%lf ",betaM);
        betaTot += betaM;
      }
     fprintf(fptr,"%lf ",betaTot);
     for(k=0; k<nstates3; k++){
	betaR = 0.0; betaI = 0.0;
        for(m=0; m<sizet; m++){
             betaR += alphaN[k*sizet+m]*initC[m]*cos(Wind[sector].evals[m]*t);
             betaI += alphaN[k*sizet+m]*initC[m]*sin(Wind[sector].evals[m]*t);
        }
        betaM = betaR*betaR + betaI*betaI;
        //fprintf(fptr,"%lf ",betaM);
        betaTot += betaM;
      }
     fprintf(fptr,"%lf \n",betaTot);
   }
   fclose(fptr);

   /* free memory */
   alphaK.clear(); alphaL.clear(); alphaM.clear(); alphaN.clear();
   initC.clear();
}
