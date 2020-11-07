#ifndef __DEFINE_H_INCLUDED__
#define __DEFINE_H_INCLUDED__

#include <vector>
#include <iostream>
#include "mkl.h"
#include "mkl_types.h"
#define DIM 2
#define PI 2*asin(1.0)

extern int CALC;
extern int CHKDIAG,STORE_SVD;
extern int *next[2*DIM+1];
extern int *nextCHK[2*DIM+1];
extern int *chk2lin,*lin2chk;
extern int LX,LY,VOL,VOL2;
extern bool k00, kPiPi, kPi0, k0Pi;
extern int nWindSector;
extern unsigned int nWind;
extern int **lookup;
extern double lam,Ti,Tf,dT;
extern int INIT,INITq;
extern int INITbag;
extern int INITphasePi0, INITphase0Pi, INITphasePiPi;
extern double inorm;
/* NTOT = total no of basis states
 * NH   = states not killed by H  */
extern int NTOT;
extern std::vector<std::vector<bool>> basis;
extern std::vector<std::vector<bool>> basis_nonflip;
extern std::vector<std::vector<bool>> basis_flip;
/* list to store the bags which have a zero form factor */
extern std::vector<int> listPiPi, listPi0, list0Pi;
/* list to store label physical bags */
extern size_t nDimPi0, nDim0Pi;
extern std::vector<int> labelPi0, label0Pi;
/* list the spurious eigenstates */
extern std::vector<int> spurPiPi, spurPi0, spur0Pi;
/* a flag to track the number of significant ice states */
extern double cutoff;

/* user defined classes */
class WindNo{
  public:
    // labeling the sectors
    int Wx,Wy;
    // numbering the total (kx,ky) sectors (#-of translation bags)
    int trans_sectors;
    // basis states
    long int nBasis;
    std::vector<std::vector<bool>> basisVec;
    std::vector<int> nflip;
    // labelling the Translation flags and the degeneracy
    std::vector<long int> Tflag;
    std::vector<int> Tdgen;
    // states in each of the "translation bags"
    std::vector<int> Tbag;
    // Note Tgen[1, ..., nBasis]; Tbag[1, ..., trans_sectors]
    // If state q is in Tbag[j], then Tgen[q]=VOL/Tbag[j]
    // relative phase of a basis vector wrt the reference state
    std::vector<int> FPi0, F0Pi, FPiPi;
    std::vector<int> isRep;

    // sum of the phase factors for the four momenta: (0,0), (pi,0), (0,pi), (pi,pi)
    std::vector<int> mom00, momPi0, mom0Pi, momPiPi;

    // flip(x)=1 (anti-clockwise), -1 (clockwise), 0 (not-flippable)
    std::vector<std::vector<int>> xflip;

    // information about the Ey(x): Ey = sum, dEy = diff,
    std::vector<std::vector<int>> Ey;
    std::vector<std::vector<int>> dEy;
    // CEy0,CEy1 = corrf of Ey in the y-direction
    std::vector<int> CEy0;
    std::vector<int> CEy1;
    // OOd1, OOd2, OOv1, OOv2, OOh1, OOh2; see schematic diagram below
    std::vector<int> OOd1, OOd2;
    std::vector<int> OOv1, OOv2;
    std::vector<int> OOh1, OOh2;
    // corrE1 = corrf of Ey at r=1 = (1/Lx) \sum_{i=1}^Lx Ey(i) Ey(i+1)
    // note that this is named differently in CODES2
    std::vector<double> corrE1;

    // full Hamiltonian in the WindNo sector stored as sparse matrix
    std::vector<MKL_INT> rows,cols;
    std::vector<double> hamil;
    double getH(int,int);
    // simple routine to check the direct way to extract
    // elements from a sparse matrix in the CSC representation
    void check_getH();
    // Hamiltonian in different (kx,ky) sectors
    std::vector<std::vector<double>> hamil_K00;
    std::vector<std::vector<double>> hamil_KPi0;
    std::vector<std::vector<double>> hamil_K0Pi;
    std::vector<std::vector<double>> hamil_KPiPi;

    void allocate_Kxy(int);
    void deallocate_Kxy(int);

    // eigenvalues and eigenvectors for the Hamiltonians in each of the four
    // momenta sectors
    std::vector<double> evals_K00, evals_KPi0, evals_K0Pi, evals_KPiPi;
    std::vector<double> evecs_K00, evecs_KPi0, evecs_K0Pi, evecs_KPiPi;

    // function to display variables
    void display(){
     printf("(Wx,Wy)=(% d,% d) with #-of-states= %ld \n",Wx,Wy,nBasis);
    }

    // Translation in X,Y
    void initTflag();
    void TransX(std::vector<bool>&, std::vector<bool>&, int);
    void TransY(std::vector<bool>&, std::vector<bool>&, int);
    void disp_Tprop();

    // count number of states in each "translation bag"
    void tbag_count();

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
       Tflag.clear(); Tdgen.clear(); Tbag.clear();
       hamil_K00.clear(); hamil_K0Pi.clear(); hamil_KPi0.clear(); hamil_KPiPi.clear();
       evals_K00.clear(); evals_K0Pi.clear(); evals_KPi0.clear(); evals_KPiPi.clear();
       evecs_K00.clear(); evecs_K0Pi.clear(); evecs_KPi0.clear(); evecs_KPiPi.clear();
       xflip.clear(); nflip.clear();  Ey.clear(); dEy.clear();
       FPi0.clear();  F0Pi.clear();   FPiPi.clear(); isRep.clear();
       mom00.clear(); momPi0.clear(); mom0Pi.clear(); momPiPi.clear();
    }
};
extern std::vector<WindNo> Wind;

// basis for subsystems A and B
//extern std::vector<std::vector<bool>> eA;
//extern std::vector<std::vector<bool>> eB;
// LEN_A, LEN_B are the sizes of subsystems A and B
// DA and DB are the respective Hilbert spaces
extern unsigned int LEN_A,LEN_B,VOL_A,VOL_B;
extern unsigned int DA, DB;

// functions to create the basis for subsystems A & B
extern void createBasis(int);
extern double schmidtDecom(std::vector<double>&, int, std::vector<MKL_INT>&);

/* routines */
void initneighbor(void);
void conststates(void);
void constH(int);
void winding_no_decompose(void);
void trans_decompose(int);
void trans_Hamil(int);
void checkCCpartners(int);
void ChargeConjEval1(int, std::vector<int>&, std::vector<int>&);
void ChargeConjEval2(int, std::vector<int>&, std::vector<int>&);
void calcCCvalues(int);
void diag_LAPACK(int, std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&);
void diag_LAPACK_RRR(int, std::vector<std::vector<double>>&, std::vector<double>&, std::vector<double>&);
void calc_Oflip(int);
void initState(int, int, int*);
void entanglementEnt_INIT0(int);
void entanglementEnt_INIT4(int);
void studyEvecs00(int, double);
void studyEvecsPiPi(int, double);
void studyEvecsPi0(int, double);
void studyEvecs0Pi(int, double);
void studyEvecs2_K00(int, double);
void studyEvecs2_KPiPi(int, double);
void studyEvecs2_KPi0(int, double);
void studyEvecs2_K0Pi(int, double);
void calc_diagEy(int, std::vector<int>&);
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
