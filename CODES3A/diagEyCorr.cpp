/* calculate the diagonal ensemble results of the EyEy correlator */
/* in the momentum (0,0) sector */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

void calc_diagEy(int sector, std::vector<int>& ScarBags){
  int k,l,m,p;
  int whichBag;
  int bag1, bag2, bag3, countBags;
  double t;
  double diagMax, diagMin, diagMid;
  double cEybag1, cEybag2, cEybag3;
  double initcEy;
  unsigned int sizet = Wind[sector].trans_sectors;
  std::vector<double> cEy(sizet, 0.0);
  std::vector<double> diagEy;
  std::vector<bool> flag(sizet, false);
  std::vector<double> sinf(sizet, 0.0), cosf(sizet, 0.0);
  double Bmk,Bml,diag_cEy;
  double phi1RE, phi1IM, phi2RE, phi2IM, phi3RE, phi3IM;
  double diff, diffmax;
  FILE *outf;

  // Calculate: < b_l | CEy | b_l >; The expression below is correct for the momenta:
  for(p=0;p<Wind[sector].nBasis;p++){
     cEy[Wind[sector].Tflag[p]-1] += Wind[sector].corrE1[p]*Wind[sector].Tdgen[p]/((double)VOL);
     if(Wind[sector].nflip[p] == VOL2){
       flag[Wind[sector].Tflag[p]-1] = true;
       //printf("%ld %d\n",Wind[sector].Tflag[p]-1,Wind[sector].nflip[p]);
     }
  }

  // bag1 and bag2 store the bags which have the maximum and minimum values of the
  // diagonal ensemble
  bag1 = 0; bag2 = 0; bag3 = 0;
  diagMax = 1000.0; diagMax = -1000.0;
  whichBag=0; countBags=0;
  outf = fopen("diag_CEy.dat","w");
  fprintf(outf,"# Results for (kx,ky)=(0,0) \n");
  for(k=0;k<sizet;k++){
    if(flag[k] == false) continue;
    // count bags which have the half the number of flippable plaquettes
    countBags++;
    diag_cEy=0.0;
    for(l=0;l<sizet;l++){
      for(m=0;m<sizet;m++){
          Bml = Wind[sector].evecs_K00[m*sizet+l];
          Bmk = Wind[sector].evecs_K00[m*sizet+k];
          diag_cEy += cEy[l]*Bmk*Bmk*Bml*Bml;
      }
    }
    diagEy.push_back(diag_cEy);
    if(ScarBags[whichBag] == k){
      fprintf(outf,"%d %d %.12lf %.12lf\n",1,ScarBags[whichBag],diag_cEy,cEy[k]);
      whichBag++;
      if(diagMin > diag_cEy) { diagMin = diag_cEy; bag1 = k; initcEy = cEy[k]; }
    }
    else{
       fprintf(outf,"%d %d %.12lf %.12lf\n",0,k,diag_cEy,cEy[k]);
       if(diagMax < diag_cEy){ diagMax = diag_cEy; bag2 = k; }
     }
  }
  fclose(outf);
  printf("size of bags which have half the flippable plaqs = %d\n",countBags);
  printf("size of diagEy   = %ld\n",diagEy.size());

  // find the third bag which has identical cEy but a diff diagonal ensemble
  whichBag=0;   // whichBag tracks the bags that contribute to the scar
  countBags=0;  // countBags track all those with the right energy density
  diffmax=0.0;
  for(k=0; k<sizet; k++){
    if(flag[k] == false) continue; // only bags with right init cEy
    if(ScarBags[whichBag]==k) whichBag++;
    else{ //interested in the non-scar bags
        if(initcEy == cEy[k]){
          diff = diagEy[countBags] - diagMin;
          if(diff > diffmax) { diffmax = diff; diagMid = diagEy[countBags]; bag3=k; }
        }
    }
    countBags++;
  }

  // bag1 is the tentative initia bag state, which contributes to the scar
  // bag2 is the initial state with the diag ensemble value farthest from the DE of bag1
  // bag3 is the initial state with same cEy as bag1, and DE farthest from bag1
  printf("Bag=%d, cEy_bag=%.12lf, cEy_DE=%.12lf\n",bag1,cEy[bag1],diagMin);
  printf("Bag=%d, cEy_bag=%.12lf, cEy_DE=%.12lf\n",bag2,cEy[bag2],diagMax);
  printf("Bag=%d, cEy_bag=%.12lf, cEy_DE=%.12lf\n",bag3,cEy[bag3],diagMid);

  // Do real time evolution with initial states in the diag min and max bags
  outf = fopen("tevol_cEy.dat","w");
  for(t=Ti; t<Tf; t=t+dT){
     cEybag1=0.0; cEybag2=0.0; cEybag3=0.0;
     for(std::size_t ii = 0; ii < sizet; ii++){
        cosf[ii]   = cos(Wind[sector].evals_K00[ii]*t);
        sinf[ii]   = sin(Wind[sector].evals_K00[ii]*t);
     }
     for(l=0; l<sizet; l++){
       phi1RE=0.0; phi1IM=0.0; phi2RE=0.0; phi2IM=0.0; phi3RE=0.0; phi3IM=0.0;
       for(k=0; k<sizet; k++){
         phi1RE += Wind[sector].evecs_K00[k*sizet+bag1]*Wind[sector].evecs_K00[k*sizet+l]*cosf[k];
         phi1IM += Wind[sector].evecs_K00[k*sizet+bag1]*Wind[sector].evecs_K00[k*sizet+l]*sinf[k];
         phi2RE += Wind[sector].evecs_K00[k*sizet+bag2]*Wind[sector].evecs_K00[k*sizet+l]*cosf[k];
         phi2IM += Wind[sector].evecs_K00[k*sizet+bag2]*Wind[sector].evecs_K00[k*sizet+l]*sinf[k];
         phi3RE += Wind[sector].evecs_K00[k*sizet+bag3]*Wind[sector].evecs_K00[k*sizet+l]*cosf[k];
         phi3IM += Wind[sector].evecs_K00[k*sizet+bag3]*Wind[sector].evecs_K00[k*sizet+l]*sinf[k];
       }
       cEybag1 += (phi1RE*phi1RE + phi1IM*phi1IM)*cEy[l];
       cEybag2 += (phi2RE*phi2RE + phi2IM*phi2IM)*cEy[l];
       cEybag3 += (phi3RE*phi3RE + phi3IM*phi3IM)*cEy[l];
     }
     fprintf(outf,"%.4f %.12le %.12le %.12le\n",t,cEybag1,cEybag2,cEybag3);
  }
  fclose(outf);
  cEy.clear(); flag.clear(); diagEy.clear();
  sinf.clear(); cosf.clear();
}
