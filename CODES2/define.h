#ifndef __DEFINE_H_INCLUDED__
#define __DEFINE_H_INCLUDED__

#include <vector>
#include <iostream>
#include "mkl.h"
#include "mkl_types.h"
#define DIM 2
#define PI 2*asin(1.0)

extern int CHKDIAG,STORE_SVD;
extern int *next[2*DIM+1];
extern int *nextCHK[2*DIM+1];
extern int *chk2lin,*lin2chk;
extern int LX,LY,VOL,VOL2;
extern int nWindSector;
extern unsigned int nWind;
extern int **lookup;
extern double lam,Ti,Tf,dT;
extern int LOG;
extern int INIT,INITq;
/* NTOT = total no of basis states
 * NH   = states not killed by H  */
extern int NTOT;
extern int WX, WY, LR;
// initial state is specified by WX, WY and LR
// WX = (1/Lx)*SUM_(x,y=0) E_{x,y}.
// WY = (1/Ly)*SUM_(x=0,y) E_{x,y}.
// LR=0 flux to the left (subsystem LA), LR=1 flux to the right (subsystem LB).

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

    // flip(x)=1 (anti-clockwise), -1 (clockwise), 0 (not-flippable)
    std::vector<std::vector<int>> xflip;

    // information about the Ey(x): Ey = sum, dEy = diff,
    std::vector<std::vector<int>> Ey;
    std::vector<std::vector<int>> dEy;
    // CEy1 = corrf of Ey at r=1 = (1/Lx) \sum_{i=1}^Lx Ey(i) Ey(i+1)
    std::vector<double> CEy1;
    // OOd1, OOd2, OOv1, OOv2, OOh1, OOh2; see schematic diagram below
    std::vector<int> OOd1, OOd2;
    std::vector<int> OOv1, OOv2;
    std::vector<int> OOh1, OOh2;

    // function to display variables
    void display(){
     printf("(Wx,Wy)=(% d,% d) with #-of-states= %ld \n",Wx,Wy,nBasis);
    }

    // function to count flippable plaquettes
    void flip_plaq();

    // function to count Ey and related operators;
    void computeEy();

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
      Wx = -LX; Wy = -LY;
      nBasis = 0;
    }

    // Destructor
    ~WindNo(){
       basisVec.clear();
       rows.clear(); cols.clear(); hamil.clear();
       evals.clear(); evecs.clear();
       xflip.clear();
       nflip.clear();
       Ey.clear();   dEy.clear();
       CEy1.clear();
       OOd1.clear(); OOd2.clear();
       OOv1.clear(); OOv2.clear();
       OOh1.clear(); OOh2.clear();
     }
};

extern std::vector<WindNo> Wind;

// basis for subsystems A and B
extern std::vector<std::vector<bool>> eA;
extern std::vector<std::vector<bool>> eB;
// LEN_A, LEN_B are the sizes of subsystems A and B
// DA and DB are the respective Hilbert spaces
extern unsigned int LEN_A,LEN_B,VOL_A,VOL_B;
extern unsigned int DA, DB, NCHI;

// Schmidt matrix and decomposed SVD
extern double sel_eval;
extern std::vector<double> sel_evec;
extern double shanonE, IPR, structE;
extern std::vector<double> chiSVD_evec;
// functions to create the basis for subsystems A & B
extern void createBasis(int);
extern double schmidtDecom(std::vector<double>&, int, std::vector<MKL_INT>&);
extern double schmidtDecomRT(std::vector<double>&, std::vector<double>&,int, std::vector<MKL_INT>&);
extern void printvec(std::vector<double>&);

/* routines */
void initneighbor(void);
void conststates(void);
void initState(int, int, int*);
void constH(int);
void winding_no_decompose(void);
void calc_Oflip(int);
void calc_CorrF(int);
double entanglementEntropy(int,int,std::vector<MKL_INT>&);
void evolve_Eent(int);
void Lecho(int);
void evolveH_ov1(int);
void evolveH_ov2(int);
void evolveH_ov3(int);
void flipped_hist(int, int, int*);
void evolve_Eent(int);
void calc_Oflipt(int);
void calc_Oflipt2(int);
void calc_Okint(int);
void evolve_corrf1(int);
void studyEvecs(int);
void studyEvecsLy4(int);
void studyEvecs2(int);
void studyEvecs2_6x4(int);
void scarDiag(int);
void collectBasis(int, int*, std::vector<std::vector<bool>>&);
void FilePrintBasis(int);
void print2file(int, int, FILE*);
#endif

/*  Schematic set-up of the diagonal correlators
      |======|======|
      |  p4  |  p3  |
      |======|======|
      |  p1  |  p2  |
      |======|======|
      OOd1 = < flip(p1) * flip(p3) >; OOd2 = < flip(p2) * flip(p4) >;
      OOv1 = < flip(p1) * flip(p4) >; OOv2 = < flip(p2) * flip(p3) >;
      OOh1 = < flip(p1) * flip(p2) >; OOh2 = < flip(p3) * flip(p4) >;
*/
