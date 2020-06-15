/* calculate the time evolution of Omag as evolution in real time */
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<vector>
#include<algorithm>
#include "define.h"

//extern void initState(int, int, int*);

// Notation: eigenstate |psi_n> = \sum_k \alpha_k |k>, |k> is a basis state in
//           specified winding number (wx,wy) sector.
void calc_Okint(int sector){
  MKL_INT p,q,sizet,k,l,m;
  MKL_INT row,col,row1,row2;
  int q1;
  double t;
  sizet =  Wind[sector].nBasis;
  std::vector<double> alpha;
  std::vector<double> phiRE(sizet),phiIM(sizet);
  double Okin_diag,Okin;
  //the imaginary part vanishes
  FILE *outf;

  /* construct cartoon state */
  if(INITq == -1){ std::cout<<"Error in initial state. Aborting. "<<std::endl; exit(0); }
  q1 = INITq;
  std::cout<<"In function calc_Okint. Starting state is basis state = "<<q1<<std::endl;
  for(p=0; p<sizet; p++){
     alpha.push_back(Wind[sector].evecs[p*sizet+q1]);
  }

  // Compute <PSI_M| Okin |PSI_M> for every eigenstate PSI_M for the diag observable
  Okin_diag = 0.0;
  // use the sparse storage matrix such that every operation is used profitably
  row = Wind[sector].rows.size()-1;
  for(k=0;k<row;k++){
     row1 = Wind[sector].rows[k]-1; row2 = Wind[sector].rows[k+1]-1;
     for(p=row1;p<row2;p++){
       l=Wind[sector].cols[p]-1;
       if(l==k) continue;
       for(m=0;m<sizet;m++){
          Okin_diag += alpha[m]*alpha[m]*Wind[sector].evecs[m*sizet+k]*Wind[sector].evecs[m*sizet+l];
       }
     }
  }
  Okin_diag = Okin_diag/((double)VOL);
  outf = fopen("OkinT.dat","w");
  fprintf(outf,"# value in the diagonal observable = %f \n",Okin_diag);
  // the time-evolution bit
  for(t=Ti;t<Tf;t=t+dT){
     Okin = 0.0;
     for(k=0;k<sizet;k++){
	phiRE[k]=0.0; phiIM[k]=0.0;
        for(m=0;m<sizet;m++){
	  phiRE[k] += alpha[m]*Wind[sector].evecs[m*sizet+k]*cos(Wind[sector].evals[m]*t);
	  phiIM[k] += alpha[m]*Wind[sector].evecs[m*sizet+k]*sin(Wind[sector].evals[m]*t);
        }
     }
     for(k=0;k<row;k++){
       row1 = Wind[sector].rows[k]-1; row2 = Wind[sector].rows[k+1]-1;
       for(p=row1;p<row2;p++){
         l=Wind[sector].cols[p]-1;
         if(l==k) continue;
         Okin += (phiRE[k]*phiRE[l] + phiIM[k]*phiIM[l]);
       }
     }
     Okin = Okin/((double)VOL);
     fprintf(outf,"%.4lf %.12lf\n",t,Okin);
  }
  fclose(outf);
  /* clear memory */
  alpha.clear();
  phiRE.clear(); phiIM.clear();
}
