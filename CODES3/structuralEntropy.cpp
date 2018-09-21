#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

double shanonE, IPR, structE;

void structuralEntropy(int sector){
   int p,i;
   double ampl;
   FILE *outf;

   outf = fopen("ShanonE.dat","w");
   fprintf(outf,"# Eigenvalues  Shanon_Entropy    Inverse_Participation_Ratio   Structural_Entropy\n");

   // Calculate the quantities for each of the eigenstates 
   for(p=0; p<Wind[sector].trans_sectors; p++){
      shanonE = 0.0; IPR = 0.0;
      for(i=0; i < Wind[sector].trans_sectors; i++){
         ampl = Wind[sector].evecs[p*Wind[sector].trans_sectors + i];
         shanonE -= (ampl*ampl)*log(ampl*ampl); 
         IPR     += (ampl*ampl*ampl*ampl);
      }
      IPR = 1.0/IPR;
      structE = shanonE - log(IPR);
      // write to file
      fprintf(outf,"%lf %lf %lf %lf\n",Wind[sector].evals[p],shanonE,IPR,structE);
  }
  fclose(outf);
}
