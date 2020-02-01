#ifndef __DEFINE_H_INCLUDED__
#define __DEFINE_H_INCLUDED__

#include <vector>
#include <iostream>
#include "mkl.h"
#include "mkl_types.h"
#define DIM 2
#define PI 2*asin(1.0)

extern int CHKDIAG;
extern int *next[2*DIM+1];
extern int *nextCHK[2*DIM+1];
extern int *chk2lin,*lin2chk;
extern int LX,LY,VOL,VOL2;
extern int QTOT,LXP,LYP,LXM,LYM;
extern int locQP,locQM;
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
    // maximum number of flippable plaqs in a basis state
    int nflipMax; 
    // Hamiltonian in the sector stored as a sparse matrix
    // this is the only way of storing the matrix
    std::vector<MKL_INT> rows,cols;
    std::vector<double> hamil;
    double getH(int,int);
    // simple routine to check the direct way to extract 
    // elements from a sparse matrix in the CSC representation
    void check_getH();

    // eigenvalues and eigenvectors in the Winding number object
    std::vector<double> evals;
    std::vector<double> evecs;

    // information about correlation functions
    std::vector<std::vector<double>> cflip;

    // information about the flip(x)
    std::vector<std::vector<bool>> xflip;

    // function to display variables
    void display(){
     printf("(Wx,Wy)=(% d,% d) with #-of-states= %ld \n",Wx,Wy,nBasis);
    } 

    // function to count flippable plaquettes
    void flip_plaq(); 

    // function to sort the basis states 
    void sortbasis();

    // function to search a transformed state in the winding number basis
    // the binary search implemented on the sorted basis states clearly 
    // outperforms the linear search
    int scan(std::vector<bool>&);
    int binscan(std::vector<bool>&);
    int binscan2(std::vector<bool>&);

    // Default constructor
    WindNo(){
      Wx = 0; Wy = 0;
      nBasis = 0;
      nflipMax = -1;
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

extern WindNo Wind;

extern void printvec(std::vector<double>&);

/* routines */
void initneighbor(void);
void conststates_Q1(void);
void conststates_Q2(void);
void constH(void);
void select_winding(void);
void calc_Oflip();

//void calc_Oflipt(int,int,int);
//void entanglementEntropy(int,int,int);
//void evolve_cartoons(int);
//void evolveH_ov1(int);
//void evolveH_ov2(int);
//void flipped_hist(int);
//void calc_Oflipt2(int,int,int);
//void calc_Okint(int,int,int);
//void evolve_corrf1(int);
#endif 
