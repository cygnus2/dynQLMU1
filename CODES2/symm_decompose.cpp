#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

 /* Decompose the full Hilbert space using symmetry transformations */
 /* Decompose according to Winding Numbers */
 void winding_no_decompose(){
  extern int** allocateint2d(int, int);
  extern void deallocateint2d(int**,int,int);
  extern void init_WindNo(std::vector<WindNo>&, int **, int, int);
  extern void initcheck_WindNo(std::vector<WindNo>&,int,int**);
  extern void display_basisNo(std::vector<WindNo>&,int);

  unsigned int i,j,count,p,p1,p2,p3,p4;
  int ix,iy,wx,wy;
  bool pxy,pyz,pzw,pwx;
  int count_flip;
  
  init_WindNo(Wind, lookup, LX+1, LY+1);
  // check the Winding number assignment
  //initcheck_WindNo(Wind,nWind,lookup);
  // scan basis states and pick out the Winding number
  for(i=0;i<basis.size();i++){
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
       p=chk2lin[j];
       p1=2*p; p2=2*next[DIM+1][p]+1; p3=2*next[DIM+2][p]; p4=2*p+1;
       pxy=basis[i][p1]; pyz=basis[i][p2]; pzw=basis[i][p3]; pwx=basis[i][p4];
       if((pxy==pyz)&&(pzw==pwx)&&(pwx!=pxy)) count_flip++;
    }
    // X-Winding (iy=0)
    wx=0;
    for(ix=0;ix<LX;ix++){
      p=2*ix+1; // linear co-ordinate of the basis state; y-link
      if(basis[i][p]) wx++;
      else wx--;
    }
    // Y-Winding (ix=0)
    wy=0;
    for(iy=0;iy<LY;iy++){
      p=2*iy*LX; // linear co-ordinate of the basis state; x-link
      if(basis[i][p]) wy++;
      else wy--;
    }
    wx=wx/2; wy=wy/2;
    /* copy to the relevant basis vector */
    count = lookup[(LX/2)+wx][(LY/2)+wy];
    Wind[count].nBasis++;
    Wind[count].basisVec.push_back(basis[i]);
  }
  display_basisNo(Wind,nWind);

  // sorting the basis states in each winding number sector
  //std::cout<<"Sorting the basis in each wind-no sector."<<std::endl;
  for(i=0; i<nWind; i++){
     Wind[i].sortbasis();
  }

  // compute the no of flippable plaquettes in the basis (ice) states
  // Note that this can only be calculated after sorting the basis!
  for(i=0; i<nWind; i++){
     Wind[i].flip_plaq();
  }

 }

/* Compute total number of sectors */
/* one sector for (0,0)
   Lx/2 and Ly/2 sectors for (Wx,0),(-Wx,0) and (0,Wy),(0,-Wy)
   Similarly (LX/2)*(LY/2) for (Wx,Wy),(-Wx,Wy),(Wx,-Wy),(-Wx,-Wy)
 */
int calc_WindNo(int LX,int LY){
  int nWindSector=1 + 2*(LX/2) + 2*(LY/2) + 4*(LX/2)*(LY/2); 
  printf("Total winding number sectors = %d\n",nWindSector);
  return nWindSector;
}

/* initialize Winding numbers */
void init_WindNo(std::vector<WindNo> &Wind,int **lookup, int row, int col){
  int count,ix,iy;
  WindNo iWind;

  count=0;
  /* (Wx,Wy) = (0,0) */
  lookup[LX/2][LY/2]=count;
  iWind.Wx = 0; iWind.Wy = 0; iWind.nBasis = 0; Wind.push_back(iWind); count++;

  /* (Wx,Wy) = (W,0) */
  for(ix=1;ix<=(LX/2);ix++){
    lookup[(LX/2)+ix][(LY/2)]=count;
    iWind.Wx = ix; iWind.Wy = 0; iWind.nBasis = 0; Wind.push_back(iWind); count++;

    lookup[(LX/2)-ix][(LY/2)]=count;
    iWind.Wx = -ix; iWind.Wy = 0; iWind.nBasis = 0; Wind.push_back(iWind); count++;
  }

  /* (Wx,Wy) = (0,W) */
  for(iy=1;iy<=LY/2;iy++){
    lookup[(LX/2)][(LY/2)+iy]=count;
    iWind.Wx = 0; iWind.Wy = iy; iWind.nBasis = 0; Wind.push_back(iWind); count++;

    lookup[(LX/2)][(LY/2)-iy]=count;
    iWind.Wx = 0; iWind.Wy = -iy; iWind.nBasis = 0; Wind.push_back(iWind); count++;
  }

  /* (Wx,Wy) */
  for(ix=1;ix<=LX/2;ix++){
  for(iy=1;iy<=LY/2;iy++){
    lookup[(LX/2)+ix][(LY/2)+iy]=count;
    iWind.Wx = ix; iWind.Wy = iy; iWind.nBasis = 0; Wind.push_back(iWind); count++;
    lookup[(LX/2)-ix][(LY/2)+iy]=count;
    iWind.Wx = -ix; iWind.Wy = iy; iWind.nBasis = 0; Wind.push_back(iWind); count++;
    lookup[(LX/2)+ix][(LY/2)-iy]=count;
    iWind.Wx = ix; iWind.Wy = -iy; iWind.nBasis = 0; Wind.push_back(iWind); count++;
    lookup[(LX/2)-ix][(LY/2)-iy]=count;
    iWind.Wx = -ix; iWind.Wy = -iy; iWind.nBasis = 0; Wind.push_back(iWind); count++;
  }}
}

void initcheck_WindNo(std::vector<WindNo> &Wind,int size,int **lookup){
  int i;
  std::cout<<"Checking the Winding Number initialization "<<std::endl;
  for(i=0;i<size;i++){
    printf(" count=%d;  (Wx,Wy)=(% d, % d); lookup[wx][wy]=%d\n",i,
    Wind[i].Wx,Wind[i].Wy,lookup[(LX/2)+Wind[i].Wx][(LY/2)+Wind[i].Wy]);
  }
}

void display_basisNo(std::vector<WindNo> &Wind,int size){
  int i,sum;
  sum=0;
  for(i=0;i<size;i++){
    //std::cout<<"count = "<<i<<std::endl;
    Wind[i].display();
    sum = sum + Wind[i].nBasis;
  }
  printf("Total states = %d\n",sum);
}

// sorts the basis states in a given Winding number sector
void WindNo::sortbasis(){
  std::sort(basisVec.begin(),basisVec.end());
}


// calculates the no of flippable plaquettes in each basis state
void WindNo::flip_plaq(){
    int i,j,p,p1,p2,p3,p4;
    int count_flip;
    bool pxy,pyz,pzw,pwx;

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
          p=chk2lin[j];
          p1=2*p; p2=2*next[DIM+1][p]+1; p3=2*next[DIM+2][p]; p4=2*p+1;
          pxy=basisVec[i][p1]; pyz=basisVec[i][p2]; pzw=basisVec[i][p3]; pwx=basisVec[i][p4];
          if((pxy==pyz)&&(pzw==pwx)&&(pwx!=pxy)) count_flip++;
      }
      nflip.push_back(count_flip);
      //std::cout<<"Basis state "<<i<<" = "<<count_flip<<std::endl;
    }
}

