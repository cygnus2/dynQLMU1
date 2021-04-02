/* calculate the local operator Oflip in a given winding number sector */
/* Oflip calculates the expectation value of the number of flippable
   plaquettes in a given eigenstate, normalized with the total number
   of plaquettes   */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

// Notation: eigenstate |psi_n> = \sum_k \alpha_k |k>, |k> is a basis state in
//           specified winding number (wx,wy) sector.
void calc_Oflip(int sector){
  MKL_INT p,q,sizet;
  double v_q, O_q, E_q;
  double Oflip_avg, CorrEy1_avg;
  double amp, shannonE, IPR;
  FILE *outf;

  sizet = Wind[sector].nBasis;
  outf = fopen("LocalOps.dat","w");
  fprintf(outf,"# Results of winding number sector (%d,%d) \n",Wind[sector].Wx,Wind[sector].Wy);
  // scan through all the eigenvalues
  for(p=0;p<sizet;p++){
    shannonE=0.0; IPR=0.0;
    // calculate the expectation value in each eigenstate
    Oflip_avg = 0.0; CorrEy1_avg = 0.0;
    for(q=0;q<sizet;q++){
      v_q = Wind[sector].evecs[p*sizet+q];
      O_q = Wind[sector].nflip[q];
      E_q = Wind[sector].CEy1[q];
      Oflip_avg   += v_q*v_q*O_q;
      CorrEy1_avg += v_q*v_q*E_q;
      if(fabs(v_q) > 1e-10) shannonE -= (v_q*v_q)*log(v_q*v_q);
      IPR += v_q*v_q*v_q*v_q;
    }
    //Oflip_avg   = Oflip_avg/((double)VOL);
    fprintf(outf,"%.12lf %.12lf %.12lf %.12lf %.12lf\n",Wind[sector].evals[p],Oflip_avg, CorrEy1_avg,shannonE,IPR);
  }
  fclose(outf);
}

void calc_CorrF(int sector){
  MKL_INT p,q,sizet;
  double v_q;
  double OOd1_avg, OOd2_avg, OOh1_avg, OOh2_avg, OOv1_avg, OOv2_avg;
  double OOd1_q, OOd2_q, OOv1_q, OOv2_q, OOh1_q, OOh2_q;
  FILE *outf;

  sizet = Wind[sector].nBasis;
  outf = fopen("CorrfGS.dat","w");
  fprintf(outf,"# Results of winding number sector (%d,%d) \n",Wind[sector].Wx,Wind[sector].Wy);
  // for the moment look only at the lowest 10 states
  for(p=0; p<10; p++){
    // calculate the expectation value in each eigenstate
    OOd1_avg = 0.0;  OOd2_avg = 0.0;
    OOh1_avg = 0.0;  OOh2_avg = 0.0;
    OOv1_avg = 0.0;  OOv2_avg = 0.0;
    for(q=0;q<sizet;q++){
      v_q = Wind[sector].evecs[p*sizet+q];
      OOd1_q = Wind[sector].OOd1[q]; OOd2_q = Wind[sector].OOd2[q];
      OOv1_q = Wind[sector].OOv1[q]; OOv2_q = Wind[sector].OOv2[q];
      OOh1_q = Wind[sector].OOh1[q]; OOh2_q = Wind[sector].OOh2[q];
      OOd1_avg += v_q*v_q*OOd1_q;
      OOd2_avg += v_q*v_q*OOd2_q;
      OOv1_avg += v_q*v_q*OOv1_q;
      OOv2_avg += v_q*v_q*OOv2_q;
      OOh1_avg += v_q*v_q*OOh1_q;
      OOh2_avg += v_q*v_q*OOh2_q;
    }
    fprintf(outf,"%.12lf %.12lf %.12lf %.12lf %.12lf %.12lf %.12lf\n",Wind[sector].evals[p], OOd1_avg,
      OOd2_avg, OOv1_avg, OOv2_avg, OOh1_avg,OOh2_avg);
  }
  fclose(outf);
}
