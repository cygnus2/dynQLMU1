/* evolve Entanglement Entropy of cartoon states in real-time with the Hamiltonian */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

std::vector<double> chiSVD_evec;

void evolve_Eent(int sector){
   
   int i,ix,iy,parity,p,q1,q2;
   int sizet,nchi,d;
   double t;
   std::vector<bool> cart1(2*VOL), cart2(2*VOL);
   std::vector<double> alpha1,alpha2;
   // chi1 and chi2 are the SVD time-dependent coefficients for 
   // the two cartoon states. Note, NCHI is not yet calculated!
   std::vector<double> chi1RE, chi1IM, chi2RE, chi2IM;
   std::vector<double> temp1, temp2;
   FILE *outf;
   double EE1, EE2, abschi1, abschi2;

   sizet = Wind[sector].nBasis;
   /* construct cartoon state */
   for(iy=0;iy<LY;iy++){
   for(ix=0;ix<LX;ix++){
    parity=(ix+iy)%2;
    p = 2*(iy*LX+ix);
    if(parity){
       cart1[p]=false; cart1[p+1]=true; cart2[p]=true;  cart2[p+1]=false;
    }
    else{
       cart1[p]=true; cart1[p+1]=false; cart2[p]=false;  cart2[p+1]=true;
    }
   }}

   for(p=0; p<sizet; p++){
     q1=Wind[sector].binscan(cart1); 
     q2=Wind[sector].binscan(cart2);
     alpha1.push_back(Wind[sector].evecs[p*sizet+q1]);
     alpha2.push_back(Wind[sector].evecs[p*sizet+q2]);
   }

  // Construct the spin basis for the sub-systems
  LEN_B = LX - LEN_A;
  VOL_A = LEN_A*LY; VOL_B = LEN_B*LY;
  createBasis(sector);

  // now reserve and initialize.
  chi1RE.reserve(NCHI); chi1IM.reserve(NCHI); chi2RE.reserve(NCHI); chi2IM.reserve(NCHI);
  temp1.reserve(NCHI); temp2.reserve(NCHI);

  // get the Schmidt coefficients for each eigenstates
  for(p=0; p<sizet; p++){
    sel_evec.clear(); 
    sel_eval = Wind[sector].evals[p];
    for(i=0; i<sizet; i++){
       sel_evec.push_back(Wind[sector].evecs[p*sizet+i]);
    }
    //std::cout<<"Vector "<<p<<std::endl;
    schmidtDecom(sel_evec,eA,eB,sector);
  }

  // the time independent part of the EE (from the zero eigenvalues)
  for(d=0;d<NCHI;d++){
    temp1[d]=0.0; temp2[d]=0.0;
    for(p=0; p<sizet; p++){
       if(chiSVD_evec[p*NCHI+d] < 1e-10) continue;
       if(fabs(Wind[sector].evals[p]) < 1e-10){
         temp1[d] = temp1[d] + alpha1[p]*chiSVD_evec[p*NCHI+d];
         temp2[d] = temp2[d] + alpha2[p]*chiSVD_evec[p*NCHI+d];
      }
    }
  }

  outf = fopen("EENT_tevol.dat","w");
  // real-time evolution of the entanglement entropy
  for(t=Ti;t<Tf;t=t+dT){
     EE1 = 0.0; EE2 = 0.0;
     for(d=0;d<NCHI;d++){
     chi1RE[d]=0.0; chi1IM[d]=0.0; chi2RE[d]=0.0; chi2IM[d]=0.0;   
     for(p=0; p<sizet; p++){
       if(chiSVD_evec[p*NCHI+d] < 1e-10) continue;
       if(fabs(Wind[sector].evals[p]) > 1e-10){
         chi1RE[d] = chi1RE[d] + alpha1[p]*chiSVD_evec[p*NCHI+d]*cos(Wind[sector].evals[p]*t);
         chi2RE[d] = chi2RE[d] + alpha2[p]*chiSVD_evec[p*NCHI+d]*cos(Wind[sector].evals[p]*t);
         chi1IM[d] = chi1IM[d] + alpha1[p]*chiSVD_evec[p*NCHI+d]*sin(Wind[sector].evals[p]*t);
         chi2IM[d] = chi2IM[d] + alpha2[p]*chiSVD_evec[p*NCHI+d]*sin(Wind[sector].evals[p]*t);
        }
     }} // close loops for p,d
     // get the Ent Entropy
     for(d=0;d<NCHI;d++){
        abschi1 = (chi1RE[d]+temp1[d])*(chi1RE[d]+temp1[d]) + chi1IM[d]*chi1IM[d]; 
        abschi2 = (chi2RE[d]+temp2[d])*(chi2RE[d]+temp2[d]) + chi2IM[d]*chi2IM[d];
        EE1 -= abschi1*log(abschi1);  
        EE2 -= abschi2*log(abschi2);  
     }
     fprintf(outf,"%.4f % lf % lf\n",t,EE1,EE2);
  }
  fclose(outf);

  chiSVD_evec.clear();
  chi1RE.clear(); chi2RE.clear(); chi1IM.clear(); chi2IM.clear();
  temp1.clear(); temp2.clear();
}

