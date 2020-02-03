#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

/* Select the sector with zero winding */
void select_winding(){
  extern void fileprintBasis();
  unsigned int i,j,count,p,q,p1,p2,p3,p4;
  int ix,iy;
  bool pxy,pyz,pzw,pwx;

  // initialize winding number basis
  Wind.nBasis=0;

  for(i=0;i<basis.size();i++){
      Wind.nBasis++;
      Wind.basisVec.push_back(basis[i]);
  }// close loop over the i

  // sort the basis states
  Wind.sortbasis();

  // compute the no of flippable plaquettes in the basis (ice) states
  Wind.flip_plaq();

  // print the basis states, this is only for checking
  fileprintBasis();
}

// sorts the basis states in a given Winding number sector
void WindNo::sortbasis(){
  std::sort(basisVec.begin(),basisVec.end());
  std::cout<<"Basis sorted."<<std::endl;
}

// calculates the no of flippable plaquettes in each basis state
void WindNo::flip_plaq(){
    int i,j,p,p1,p2,p3,p4;
    int count_flip;
    int fmax;
    bool pxy,pyz,pzw,pwx;

    fmax = -1;
    //printf("In sector (% d, % d); basis states = %ld \n",Wx,Wy,nBasis);
    // error message for zero basis
    if(!nBasis){ std::cout<<"Total basis state not defined! "<<std::endl; exit(0); }
    for(i=0; i<nBasis; i++){
      // compute the diagonal term in the Hamiltonian
      /* a single plaquette is arranged as
                pzw
             o-------o
             |       |
        pwx  |   p   |  pyz
             |       |
             o-------o
                pxy
      */
      count_flip=0;
      for(j=0;j<VOL;j++){
          //p=chk2lin[j];
	        p=j;
          p1=2*p; p2=2*next[DIM+1][p]+1; p3=2*next[DIM+2][p]; p4=2*p+1;
          pxy=basisVec[i][p1]; pyz=basisVec[i][p2]; pzw=basisVec[i][p3]; pwx=basisVec[i][p4];
          if((pxy==pyz)&&(pzw==pwx)&&(pwx!=pxy)) count_flip++;
      }
      nflip.push_back(count_flip);
      if(count_flip > fmax) fmax = count_flip;
    }
    nflipMax = fmax;
    std::cout<<"Max possible flippable plaq in any basis state = "<<nflipMax<<std::endl;
    // printing the flippable states
    //for(i=0; i<nBasis; i++){
    //   std::cout<<" state ="<<i<<" #-flippable plaquettes="<<Wind.nflip[i]<<std::endl;
    //}
}

void fileprintBasis(){
 FILE *fptr;
 int i,p;
 fptr=fopen("BASISW00.dat","w");
 fprintf(fptr,"Printing the (flippable) basis states for (Lx,Ly)=(%d,%d) lattice.\n",LX,LY);
 for(i=0;i<Wind.nBasis;i++){
   fprintf(fptr,"basis %d ",i);
   for(p=0;p<2*VOL;p++){
    fprintf(fptr,"%d ",(int)Wind.basisVec[i][p]);
   }
   fprintf(fptr,"\n");
  }
 fclose(fptr);
}


 /*  Notes for selecting Zero Winding Number states with charges
  *  ======================================================
  *    Note that winding number cannot be defined in the presence
  *    of charges. To see this, consider the simple topological
  *    argument: The flux from a charge can either directly go
  *    the anti-charge using the shortest route. Equivalently,
  *    it can also wind around several times and then reach the
  *    anti-charge. This is equivalent to different winding
  *    number sectors.
  */
