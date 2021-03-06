#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>
#include<chrono>

void print_matrixTbasis(char *str, int, int);
void print_eigsys(int size, std::vector<double>&,  std::vector<double>&);
void check_eigsys(int);
/* Decompose the Winding Number sector into states connected by translation flags */
void trans_decompose(int sector){
   int i,flag,count;
   int ix,iy,q;
   std::vector<bool> init(2*VOL);
   std::vector<bool> new1(2*VOL);
   std::vector<bool> new2(2*VOL);
   std::vector<bool> phaseSET;

   // initialize the translation flag
   Wind[sector].initTflag();
   std::cout<<"No of basis states in this sector = "<<Wind[sector].nBasis<<std::endl;

   // Group states into "bags" related by translation and record degeneracies. For each
   // cartoon state consider translations of the form (Tx,Ty). Start with (0,0) which
   // corresponds to no translation. Every time a new state is generated, label it with
   // the same flag, and increase degeneracy in the loop. The inidividual basis (ice)
   // states in the Winding number sector are initialized with 0, the first bag of states
   // related by lattice translations is labelled as 1, the second one as 2, and so on.
   // The total number of bags is the highest counter. The degeneracy is initialized with zero.
   // #-of-states in bag[k] = VOL/(Tgen); Tgen = degeneracy of states in bag-k
   // For constructing the momenta in other sectors, we need phases of a basis state with
   // respect to a representative state in each bag. The vectors F0Pi, FPi0, FPiPi record
   // these phases of each basis state wrt the representative in each bag.
   flag=1;
   for(i=0;i<Wind[sector].nBasis;i++){
      // when the flag is non-zero, this has already been considered
      if(Wind[sector].Tflag[i]) continue;
      init = Wind[sector].basisVec[i];
      // set flag for a unassigned state
      Wind[sector].Tflag[i]=flag;
      // take this as the representative state of the bag-no=flag
      Wind[sector].isRep[i]=1;
      Wind[sector].F0Pi[i]=1;
      Wind[sector].FPi0[i]=1;
      Wind[sector].FPiPi[i]=1;
      // go through all possible translations, in x and in y
      for(iy=0;iy<LY;iy++){
      for(ix=0;ix<LX;ix++){
         Wind[sector].WindNo::TransX(init,new1,ix);
         Wind[sector].WindNo::TransY(new1,new2,iy);
         //check which basis state it corresponds to; function binscan is much faster!
         //q = Wind[sector].scan(new2);
         q = Wind[sector].binscan(new2);
         // set flag of the new state
         if(q==-100){ printf("Error. Translated state not in Winding no sector. \n"); exit(0); }
         Wind[sector].Tflag[q]=flag;
         // increase degeneracy
         Wind[sector].Tdgen[q]++;
         // set the phases of state q wrt the reference state of the bag
         Wind[sector].FPi0[q] = 1 - ((ix & 1) << 1);
         Wind[sector].F0Pi[q] = 1 - ((iy & 1) << 1);
         Wind[sector].FPiPi[q]= 1 - (((ix+iy) & 1)  << 1);
         //checks
         //if(Wind[sector].FPi0[q]!=1 && Wind[sector].FPi0[q]!=-1){
         //  printf("ERROR! q=%d, FPi0=%d\n",q,Wind[sector].FPi0[q]); exit(0); }
         //if(Wind[sector].F0Pi[q]!=1 && Wind[sector].F0Pi[q]!=-1){
         //  printf("ERROR! q=%d, FPi0=%d\n",q,Wind[sector].F0Pi[q]); exit(0); }
         //if(Wind[sector].FPiPi[q]!=1 && Wind[sector].FPiPi[q]!=-1){
         //  printf("ERROR! q=%d, FPiPi=%d\n",q,Wind[sector].FPiPi[q]); exit(0); }
      }}
      flag++;
   }
   flag--;
   Wind[sector].trans_sectors = flag;
   Wind[sector].tbag_count();

   // initialize the phase variables
   for(i=0;i<Wind[sector].trans_sectors;i++){
        Wind[sector].mom00.push_back(0);  Wind[sector].momPi0.push_back(0);
        Wind[sector].mom0Pi.push_back(0); Wind[sector].momPiPi.push_back(0);
        phaseSET.push_back(false);
   }

   // Compute the sum of phases for each translation bag which map the state to itself
   // F(a)=sum_{x=0,...,Lx-1; y=0,...,Ly-1} exp(-i*kx*x - i*ky*y).
   // A single state in a translation bag needs to be considered.
   for(i=0;i<Wind[sector].nBasis;i++){
     q=Wind[sector].Tflag[i];
     if(phaseSET[q-1]) continue;  // if sum of phases is calculated, skip
     else{                        // otherwise do all translations and compute
        init = Wind[sector].basisVec[i];
        for(iy=0;iy<LY;iy++){
        for(ix=0;ix<LX;ix++){
          Wind[sector].WindNo::TransX(init,new1,ix);
          Wind[sector].WindNo::TransY(new1,new2,iy);
          if(init == new2){
            // sum the phases for (pi,0), (0,pi), (pi,pi); sum for (0,0) is trivial
            // note that the operation (-1)^n = 1 - (n%2)*2 = 1 - ((n & 1) << 1)
            Wind[sector].mom00[q-1]   += 1;
            Wind[sector].momPi0[q-1]  += 1 - ((ix & 1) << 1);
            Wind[sector].mom0Pi[q-1]  += 1 - ((iy & 1) << 1);
            Wind[sector].momPiPi[q-1] += 1 - (((ix+iy) & 1) << 1);
          }
        }}
        phaseSET[q-1]=true;
        if(Wind[sector].momPiPi[q-1] == 0) listPiPi.push_back(q-1);
        if(Wind[sector].momPi0[q-1] == 0)  listPi0.push_back(q-1);
        if(Wind[sector].mom0Pi[q-1] == 0)  list0Pi.push_back(q-1);
     }
   }

   // relabel the translation bags in the Pi0 and 0Pi sectors
   if(listPiPi.size() > 0) {
     printf("Why are unphysical bags present in sector (Pi,Pi)? \n"); exit(0);
   }
   nDimPi0=0;   nDim0Pi=0;
   // note that this is an ordered list
   for(q=0; q<Wind[sector].trans_sectors; q++){
      if(Wind[sector].momPi0[q]){ labelPi0.push_back(nDimPi0); nDimPi0++; }
      else labelPi0.push_back(-997);
      if(Wind[sector].mom0Pi[q]){ label0Pi.push_back(nDim0Pi); nDim0Pi++; }
      else label0Pi.push_back(-997);
   }

   // print and display statement to check, debug etc
   std::cout<<"No of (kx,ky) sectors="<<Wind[sector].trans_sectors<<std::endl;
   Wind[sector].disp_Tprop();
   init.clear(); new1.clear(); new2.clear();
   phaseSET.clear();
}

void WindNo::initTflag(){
    int i;
    for(i=0;i<nBasis;i++){
      Tflag.push_back(0);
      Tdgen.push_back(0);
      FPi0.push_back(997);
      F0Pi.push_back(997);
      FPiPi.push_back(997);
      isRep.push_back(0);
    }
    // 997 is the largest prime smaller than 1000!
}

// Translate a basis state in a given winding number in the x-direction
// state --> stateTx;  by lattice translations (lx,0) and return the translated state
void WindNo::TransX(std::vector<bool> &state,std::vector<bool> &stateT,int lx){
    int ix,iy,p,q;
    for(iy=0;iy<LY;iy++){
    for(ix=0;ix<LX;ix++){
       p = iy*LX + ix;         // initial point
       q = iy*LX + (ix+lx)%LX; // shifted by lx lattice units in +x dir
       stateT[2*q]=state[2*p]; stateT[2*q+1]=state[2*p+1];
    }}
}

// Translate a basis state in a given winding number in the y-direction
// state --> stateTx;  by lattice translations (0,ly) and return the translated state
void WindNo::TransY(std::vector<bool> &state,std::vector<bool> &stateT,int ly){
    int ix,iy,p,q;
    for(iy=0;iy<LY;iy++){
    for(ix=0;ix<LX;ix++){
       p = iy*LX + ix;            // initial point
       q = ((iy+ly)%LY)*LX + ix;  // shifted by ly lattice units in +y dir
       stateT[2*q]=state[2*p]; stateT[2*q+1]=state[2*p+1];
    }}
}

int WindNo::binscan(std::vector<bool> &newstate){
     unsigned int m;
     // binary search of the sorted array
     std::vector<std::vector<bool>>::iterator it;
     it = std::lower_bound(basisVec.begin(),basisVec.end(),newstate);
     m  = std::distance(basisVec.begin(),it);
     if(it == basisVec.end()){
       std::cout<<"Element not found here! "<<std::endl;
       return -100;
     }
     return m;
}

// Alternate implementation of the search. Find where the required element
// exists in the vector with the iterator, then the index.
int WindNo::binscan2(std::vector<bool> &newstate){
     unsigned int m;
     std::vector<std::vector<bool>>::iterator it;
     it = std::find(basisVec.begin(),basisVec.end(),newstate);
     if(it == basisVec.end()){
       //std::cout<<"Element not found here! "<<std::endl;
       return -100;
     }
     m  = std::distance(basisVec.begin(),it);
     return m;
}

void WindNo::disp_Tprop(){
   int i;
   int count00, countPiPi, countPi0, count0Pi;
   //for(i=0;i<nBasis;i++){
   // printf("flag[%d]=%ld; Deg=%d; phasePi0=% d, phase0Pi=% d, phasePiPi=% d\n",i,Tflag[i],Tdgen[i],FPi0[i],F0Pi[i],FPiPi[i]);
   //}
   /*std::cout<<"Total number of translation-bags ="<<trans_sectors<<std::endl;
   for(i=0;i<trans_sectors;i++){
     std::cout<<"#-of-states in bag-"<<i+1<<"  ="<<Tbag[i]<<std::endl;
     std::cout<<"sum of phase factors for mom (0,0)   ="<< mom00[i]  <<std::endl;
     std::cout<<"sum of phase factors for mom (Pi,0)  ="<< momPi0[i] <<std::endl;
     std::cout<<"sum of phase factors for mom (0,Pi)  ="<< mom0Pi[i] <<std::endl;
     std::cout<<"sum of phase factors for mom (Pi,Pi) ="<< momPiPi[i]<<std::endl;
   }*/
   // Count and display which bags have zero form factors
   /*
   count00   = std::count(mom00.begin(), mom00.end(), 0);
   countPiPi = std::count(momPiPi.begin(), momPiPi.end(), 0);
   countPi0  = std::count(momPi0.begin(), momPi0.end(), 0);
   count0Pi  = std::count(mom0Pi.begin(), mom0Pi.end(), 0);
   std::cout<<"No of bags which do not contribute in sector (0,0)= "<<count00<<std::endl;
   std::cout<<"No of bags which do not contribute in sector (Pi,Pi)= "<<countPiPi<<std::endl;
   std::cout<<"No of bags which do not contribute in sector (Pi,0)= "<<countPi0<<std::endl;
   std::cout<<"No of bags which do not contribute in sector (0,Pi)= "<<count0Pi<<std::endl;
   */
   // display the same information via the lists
   std::cout<<"No of non-contributing bags in sector (Pi,Pi)="<<listPiPi.size()<<std::endl;
   for(i=0; i<listPiPi.size(); i++) std::cout<< listPiPi.at(i) << ' ';
   std::cout<<std::endl;
   std::cout<<"No of non-contributing bags in sector (Pi,0)="<<listPi0.size()<<std::endl;
   for(i=0; i<listPi0.size(); i++) std::cout<< listPi0.at(i) << ' ';
   std::cout<<std::endl;
   std::cout<<"No of non-contributing bags in sector (0,Pi)="<<list0Pi.size()<<std::endl;
   for(i=0; i<list0Pi.size(); i++) std::cout<< list0Pi.at(i) << ' ';
   std::cout<<std::endl;
   std::cout<<"Size of the Hilbert space in sector (Pi,0)="<<nDimPi0<<std::endl;
   std::cout<<"Size of the Hilbert space in sector (0,Pi)="<<nDim0Pi<<std::endl;
}

void trans_Hamil(int sector){
  extern void diag_LAPACK(int,std::vector<std::vector<double>>&,std::vector<double>&,
      std::vector<double>&);
  int k,l,kk,ll;
  int tsect, sizet;
  double ele, norm;
  FILE *fptr;

  sizet = Wind[sector].nBasis;
  tsect = Wind[sector].trans_sectors;
  // allocate space for hamil_Kxy
  Wind[sector].allocate_Kxy(INIT);
  std::cout<<"=================================================================================="<<std::endl;
  std::cout<<"Constructing the Hamiltonian in the (kx,ky)=(0,0); (Pi,Pi) sectors."<<std::endl;
  std::cout<<"Dimension of matrix = "<<tsect<<std::endl;
  std::cout<<"Value of INIT flag ="<<INIT<<std::endl;

  // Get starting timepoint
  auto start = std::chrono::high_resolution_clock::now();

  // It seems like the trick of used below only works for Ly=2
  // For Ly=4, we use the brute force approach
  // construct the hamil_kxy matrix in the (0,0) and (Pi,Pi) sectors
  if(LY==2){
    for(std::size_t i=0;i<sizet;i++){
       k=Wind[sector].Tflag[i]-1;
       // off-diagonal elements
       for(std::size_t j=i+1;j<sizet;j++){
         ele = Wind[sector].getH(i,j);
         if(ele == 0) continue;
         l   = Wind[sector].Tflag[j]-1;
         norm= sqrt(Wind[sector].Tdgen[i]*Wind[sector].Tdgen[j])/((double)VOL);

         Wind[sector].hamil_K00[k][l] +=  2*ele*norm;
         //if((Wind[sector].momPiPi[k]) && (Wind[sector].momPiPi[l]))  // none of these are zero
         Wind[sector].hamil_KPiPi[k][l] += 2*ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[j]*norm;
       }
       // diagonal elements
       ele = Wind[sector].getH(i,i);
       if(ele == 0) continue;
       norm= Wind[sector].Tdgen[i]/((double)VOL);
       Wind[sector].hamil_K00[k][k] +=  ele*norm;
       if(Wind[sector].momPiPi[k])
         Wind[sector].hamil_KPiPi[k][k] += ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[i]*norm;
    }
  } // close if (LY==2)
  else if(LY==4){
    // construct the hamil_kxy matrix brute force
    for(std::size_t i=0;i<sizet;i++){
      k=Wind[sector].Tflag[i]-1;
      for(std::size_t j=0;j<sizet;j++){
        ele = Wind[sector].getH(i,j);
        if(ele == 0) continue;
        l = Wind[sector].Tflag[j]-1;
        norm = sqrt(Wind[sector].Tdgen[i]*Wind[sector].Tdgen[j])/((double)VOL);
        Wind[sector].hamil_K00[k][l] += ele*norm;
        if((Wind[sector].momPiPi[k]) && (Wind[sector].momPiPi[l]))
            Wind[sector].hamil_KPiPi[k][l] += ele*Wind[sector].FPiPi[i]*Wind[sector].FPiPi[j]*norm;
      }
    }
  }// close if(LY==4)

  // construct the hamil_Kxy for (Pi,0)
  if(listPi0 == list0Pi){
    printf("Unphysical bags are the same in both sectors, constructing the matrices together \n");
  }
  else{ printf("Matrix needs to be separately constructed.\n"); exit(0); }
  if(LY==2){
    for(std::size_t i=0; i<sizet; i++){
       k=Wind[sector].Tflag[i]-1;
       // off-diagonal elements
       for(std::size_t j=i+1;j<sizet;j++){
         ele = Wind[sector].getH(i,j);
         if(ele == 0) continue;
         l   = Wind[sector].Tflag[j]-1;
         norm= sqrt(Wind[sector].Tdgen[i]*Wind[sector].Tdgen[j])/((double)VOL);
         kk  = labelPi0[k];  ll = labelPi0[l];
         if(kk == -997 || ll == -997) continue;
         Wind[sector].hamil_KPi0[kk][ll] += 2*ele*Wind[sector].FPi0[i]*Wind[sector].FPi0[j]*norm;
         Wind[sector].hamil_K0Pi[kk][ll] += 2*ele*Wind[sector].F0Pi[i]*Wind[sector].F0Pi[j]*norm;
       }
       // diagonal elements
       ele = Wind[sector].getH(i,i);
       if(ele == 0) continue;
       norm= Wind[sector].Tdgen[i]/((double)VOL);
       kk = labelPi0[k];
       if(kk == -997) continue;
       // phases appear as squared, and cancel out
       Wind[sector].hamil_KPi0[kk][kk] += ele*norm;
       kk = label0Pi[k];
       if(kk == -997) continue;
       Wind[sector].hamil_K0Pi[kk][kk] += ele*norm;
    }
  }//close if(LY==2)
  else if(LY==4){
    for(std::size_t i=0; i<sizet; i++){
       k=Wind[sector].Tflag[i]-1;
       // off-diagonal elements
       for(std::size_t j=0;j<sizet;j++){
         ele = Wind[sector].getH(i,j);
         if(ele == 0) continue;
         l   = Wind[sector].Tflag[j]-1;
         norm= sqrt(Wind[sector].Tdgen[i]*Wind[sector].Tdgen[j])/((double)VOL);
         kk  = labelPi0[k];  ll = labelPi0[l];
         if(kk == -997 || ll == -997) continue;
         Wind[sector].hamil_KPi0[kk][ll] += ele*Wind[sector].FPi0[i]*Wind[sector].FPi0[j]*norm;
         Wind[sector].hamil_K0Pi[kk][ll] += ele*Wind[sector].F0Pi[i]*Wind[sector].F0Pi[j]*norm;
       }
     }
  }// close if(LY==4)

  // Get ending timepoint
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time taken for Hamiltonian construction in t-basis="<<duration.count()<< " secs"<<std::endl;

  // clear file
  if(CHKDIAG){
    fptr=fopen("eigencheck.dat","w");
    fclose(fptr);
  }

  // remove memory for the original Hamiltonian which is not needed any more
  Wind[sector].hamil.clear();
  Wind[sector].rows.clear();
  Wind[sector].cols.clear();

  // diagonalize the matrixes with a LAPACK routine
  start = std::chrono::high_resolution_clock::now();
  diag_LAPACK_RRR(Wind[sector].trans_sectors,Wind[sector].hamil_K00,
    Wind[sector].evals_K00,Wind[sector].evecs_K00);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time to diag H in (0,0) sector="<<duration.count()<< " secs"<<std::endl;

  start = std::chrono::high_resolution_clock::now();
  diag_LAPACK_RRR(Wind[sector].trans_sectors,Wind[sector].hamil_KPiPi,
    Wind[sector].evals_KPiPi,Wind[sector].evecs_KPiPi);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time to diag H in (Pi,Pi) sector="<<duration.count()<< " secs"<<std::endl;

  start = std::chrono::high_resolution_clock::now();
  diag_LAPACK_RRR(nDimPi0,Wind[sector].hamil_KPi0,Wind[sector].evals_KPi0,Wind[sector].evecs_KPi0);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time to diag H in (Pi,0) sector="<<duration.count()<< " secs"<<std::endl;

  start = std::chrono::high_resolution_clock::now();
  diag_LAPACK_RRR(nDim0Pi,Wind[sector].hamil_K0Pi,Wind[sector].evals_K0Pi,Wind[sector].evecs_K0Pi);
  stop = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout<<"Time to diag H in (Pi,0) sector="<<duration.count()<< " secs"<<std::endl;

  // deallocate space for hamil_Kxy
  Wind[sector].deallocate_Kxy(INIT);
}

// print the Hamiltonian in the translation state basis
void print_matrixTbasis(char *str, int sector, int whichMom){
  int k,l;
  printf("===== %s =====\n",str);
  for(k=0;k<Wind[sector].trans_sectors;k++){
  for(l=0;l<Wind[sector].trans_sectors;l++){
    if(whichMom==1) printf("% .4lf ",Wind[sector].hamil_K00[k][l]);
    else if(whichMom==2) printf("% .4lf ",Wind[sector].hamil_KPi0[k][l]);
  }
  printf("\n");
  }
  if(whichMom==3 || whichMom==4){
     printf("Matrix print function not written. \n");
  }
}

double WindNo::getH(int p,int q){
   double ele;
   int row1,row2;
   //int c;

   ele=0.0;
   row1 = rows[p]-1; row2 = rows[p+1]-1;
   for(std::size_t c=row1;c<row2;c++){
     if((q+1)==cols[c]) { ele=hamil[c]; break; }
   }
   return ele;
}

void WindNo::check_getH(){
    int i,j;
    double ele;
    printf("==================\n");
    printf("Printing the full matrix for checking \n");
    for(i=0;i<nBasis;i++){
    for(j=0;j<nBasis;j++){
      ele = getH(i,j);
      printf("% .1lf ",ele);
    }
    printf("\n");
  }
}

void WindNo::allocate_Kxy(int initflag){
  std::size_t i,j;
  std::vector<double> myvec;
  if(trans_sectors <= 0){ printf("Error in allocation. \n"); exit(0); }
  // allocate for H00 and HPiPi
  if(initflag==0 || initflag==4){
    for(i=0; i<trans_sectors; i++) myvec.push_back(0.0);
    for(i=0; i<trans_sectors; i++){
      hamil_K00.push_back(myvec);  hamil_KPiPi.push_back(myvec);
    }
    myvec.clear();
  }
  if(initflag==4){
    // allocate for HPi0
    for(i=0; i<nDimPi0; i++) myvec.push_back(0.0);
    for(i=0; i<nDimPi0; i++) hamil_KPi0.push_back(myvec);
    myvec.clear();
    // allocate for HPi0
    for(i=0; i<nDim0Pi; i++) myvec.push_back(0.0);
    for(i=0; i<nDim0Pi; i++) hamil_K0Pi.push_back(myvec);
    myvec.clear();
  }
  else{
    std::cout<<"Wrong flag here. Aborting. "<<std::endl;
    exit(0);
  }
  printf("Successful allocation. \n");
  printf("======================\n");
}

void WindNo::deallocate_Kxy(int initflag){
  if(initflag==0 || initflag==4){
    hamil_K00.clear();
    hamil_KPiPi.clear();
  }
  if(initflag==4){
    hamil_KPi0.clear();
    hamil_K0Pi.clear();
  }
  else{
    std::cout<<"Wrong flag here. Aborting. "<<std::endl;
    exit(0);
  }
  printf("Successful de-allocation. \n");
  printf("======================\n");
}

void WindNo::tbag_count(){
   int i,q,ix,iy;
   int sum;
   sum = 0;
   for(i=0;i<trans_sectors;i++){
       // count the #-of-states with label i+1
       q = std::count(Tflag.begin(), Tflag.end(), i+1);
      Tbag.push_back(q);
      sum += q;
   }
   if(sum!=nBasis) std::cout<<"Basis mismatch: sum="<<sum<<";  nBasis="<<nBasis<<std::endl;
   sum = 0;
   sum = std::count(isRep.begin(), isRep.end(), 1);
   if(sum!=trans_sectors) std::cout<<"total #-of reference states don't match translation sectors"<<std::endl;
 }

void print_eigsys(int size, std::vector<double>& evals,  std::vector<double>& evecs){
  int i,j;
  double norm;
  for(i=0; i<size; i++){
     printf("Eigenvalue = %lf\n",evals[i]);
     norm = 0.0;
     for(j=0; j<size; j++){
       norm += evecs[size*i + j]*evecs[size*i + j];
       printf(" %lf ",evecs[size*i + j]);
     }
     printf("\n ");
     if(fabs(norm - 1.0) > 1e-8) printf("norm = %lf\n",norm);
  }
}

// check if the physical eigenvectors have any components on the bags with zero FF
void check_eigsys(int sector){
  printf("This function is redundant here; By construction sectors");
  printf("(Pi,0) and (0,Pi) have only physical bags.\n");
}
