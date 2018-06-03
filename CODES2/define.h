#ifndef __DEFINE_H_INCLUDED__
#define __DEFINE_H_INCLUDED__

#include <vector>
#include <iostream>
#include "mkl.h"
#include "mkl_types.h"
#define DIM 2
#define PI 2*asin(1.0)

extern int *next[2*DIM+1];
extern int *nextCHK[2*DIM+1];
extern int *chk2lin,*lin2chk;
extern int LX,LY,VOL,VOL2;
extern int nWindSector;
extern unsigned int nWind;
extern int **lookup;
extern double lam,Ti,Tf,dT;
/* NTOT = total no of basis states 
 * NH   = states not killed by H  */
extern int NTOT;

extern std::vector<std::vector<bool>> basis;
extern std::vector<std::vector<bool>> basis_nonflip;
extern std::vector<std::vector<bool>> basis_flip;

/* user defined classes */
class WindNo{
  public:
    // labeling the sectors
    int Wx,Wy;
    // basis states
    long int nBasis;
    std::vector<std::vector<bool>> basisVec;
    std::vector<int> nflip;
    // Hamiltonian in the sector
    std::vector<MKL_INT> rows,cols;
    std::vector<double> hamil;
    // eigenvalues and eigenvectors
    std::vector<double> evals;
    std::vector<std::vector<double>> evecs;

    // function to display variables
    void display(){
     printf("(Wx,Wy)=(% d,% d) with #-of-states= %ld \n",Wx,Wy,nBasis);
    } 

    // function to calculate the off-diagonal matrix elements of the Hamiltonian
    // in the winding number basis
    int scan(std::vector<bool>&);

    // Default constructor
    WindNo(){
      Wx = -LX; Wy = -LY;
      nBasis = 0;
    }
   
    // Destructor 
    ~WindNo(){
       basisVec.clear();
       rows.clear();
       cols.clear();
       hamil.clear();
       evals.clear();
       evecs.clear();
     }
};

extern std::vector<WindNo> Wind;

/* routines */
void initneighbor(void);
void conststates(void);
void constH(int);
void evolve_cartoons(std::vector<double>&, std::vector<std::vector<double>>&);
void winding_no_decompose(void);
void calc_Oflip(int);
#endif 
