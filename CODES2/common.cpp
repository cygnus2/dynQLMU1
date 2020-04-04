#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "define.h"

void initneighbor(void)
{
  int p,x,y,q;
  int lin;

  /* lattice linear co-ordinates */
  for(p=0;p<VOL;p++){
      x = (p%LX);
      y = (p/LX)%LY;

      next[DIM][p] = p;
      next[DIM+1][p] = y*LX + ((x+1)%LX);
      next[DIM+2][p] = ((y+1)%LY)*LX + x;
      next[DIM-1][p] = y*LX + ((x-1+LX)%LX);
      next[DIM-2][p] = ((y-1+LY)%LY)*LX + x;
  }
 
  /* checkerboard to linear and vice-versa */
  p=0;
  for(y=0;y<LY;y++){
  for(x=0;x<LX;x++){
    if((x+y)%2==0){
       q = y*LX + x;
       chk2lin[p] = q; lin2chk[q] = p;
       //printf("checkerboard site = %d(%d,%d); Lin site = %d\n",p,x,y,q);
       p++;
    }
  }}
  for(y=0;y<LY;y++){
  for(x=0;x<LX;x++){
    if((x+y)%2==1){
       q = y*LX + x;
       chk2lin[p]= q; lin2chk[q] = p;
       //printf("checkerboard site = %d(%d,%d); Lin site = %d\n",p,x,y,q);
       p++;
    }
  }} 

  // nextCHK targets the neighboring sites in the string 
  // representing the GL implementation on the even sites
  for(p=0;p<VOL;p++){
      lin = chk2lin[p];
      nextCHK[DIM][p]   = p;
      nextCHK[DIM+1][p] = lin2chk[next[DIM+1][lin]];
      nextCHK[DIM-1][p] = lin2chk[next[DIM-1][lin]]; 
      nextCHK[DIM+2][p] = lin2chk[next[DIM+2][lin]];
      nextCHK[DIM-2][p] = lin2chk[next[DIM-2][lin]];
      //printf("checkerboard NN. site = %d (in chkrboard = %d)\n",lin,p);
      //printf("+x=%d,   -x=%d,   +y=%d,   -y=%d\n",nextCHK[DIM+1][p],
      //nextCHK[DIM-1][p],nextCHK[DIM+2][p],nextCHK[DIM-2][p]); 
   }
}

/* configurations are stored such that i=2*p and i=2*p+1
   denotes the x and the y links of the site p=y*LX + x.
 */
//void printconf(bool *conf){
void printconf(std::vector<bool> conf){
  int p,d,x,y;
  for(p=0;p<VOL;p++){
   x=p%LX; y=p/LX;
   if(conf[2*p]==true)         printf("site (x,y)=(%d,%d); link-x=% d \n",x,y,1);
   else if(conf[2*p]==false)   printf("site (x,y)=(%d,%d); link-x=% d \n",x,y,-1);
   if(conf[2*p+1]==true)       printf("site (x,y)=(%d,%d); link-y=% d \n",x,y,1);
   else if(conf[2*p+1]==false) printf("site (x,y)=(%d,%d); link-y=% d \n",x,y,-1);
  }
}


/* Prints all the basis; use only for checking */
void printbasis(){
 FILE *fptr;
 int i,p;
 fptr=fopen("BASIS_FLIP.txt","w");
 fprintf(fptr,"Printing the (flippable) basis states for (Lx,Ly)=(%d,%d) lattice.\n",LX,LY);
 for(i=0;i<NTOT;i++){
   fprintf(fptr,"basis %d ",i);
   for(p=0;p<2*VOL;p++){
    fprintf(fptr,"%d ",(int)basis_flip[i][p]);
   }
   fprintf(fptr,"\n");
  }
 fclose(fptr);
}

void FilePrintBasis(int sector){
 FILE *fptr;
 int i,j,p,p1,p2,p3,p4;
 bool pxy,pyz,pzw,pwx;
 int sizet;

 sizet = Wind[sector].nBasis;
 fptr=fopen("BASISFL.dat","w");
 for(i=0;i<sizet;i++){
   for(j=0;j<VOL;j++){
      //p=chk2lin[j];
      p=j;
      p1=2*p; p2=2*next[DIM+1][p]+1; p3=2*next[DIM+2][p]; p4=2*p+1;
      pxy=Wind[sector].basisVec[i][p1]; pyz=Wind[sector].basisVec[i][p2]; 
      pzw=Wind[sector].basisVec[i][p3]; pwx=Wind[sector].basisVec[i][p4];
      if((pxy==pyz)&&(pzw==pwx)&&(pwx!=pxy)) fprintf(fptr,"%d ", 1);
      else                                   fprintf(fptr,"%d ", 0);
   }
   fprintf(fptr,"\n");
  }
 fclose(fptr);
}


int **allocateint2d(int row, int col){
  int i,j;
  int **mat;
  mat = (int **)malloc(row*sizeof(double*));
  if(mat==NULL) {printf("Out of memory\n"); exit(0);}

  for(i=0;i<row;i++){
   mat[i]=(int *)malloc(col*sizeof(double));
   if(mat[i]==NULL)  {printf("Out of memory\n"); exit(0);}
  }

 for(i=0;i<row;i++)
 for(j=0;j<col;j++)
  mat[i][j]=0;

 return mat;
}

void deallocateint2d(int **mat, int row, int col){
  int i;
  for(i=0;i<row;i++)
   free(mat[i]);

 free(mat);
}

double **allocatedouble2d(int row, int col){
  int i,j;
  double **mat;
  mat = (double **)malloc(row*sizeof(double*));
  if(mat==NULL) {printf("Out of memory\n"); exit(0);}

  for(i=0;i<row;i++){
   mat[i]=(double *)malloc(col*sizeof(double));
   if(mat[i]==NULL)  {printf("Out of memory\n"); exit(0);}
  }

 for(i=0;i<row;i++)
 for(j=0;j<col;j++)
  mat[i][j]=0.0;

 return mat;
}

void deallocatedouble2d(double **mat, int row, int col){
  int i;
  for(i=0;i<row;i++)
   free(mat[i]);

 free(mat);
}


/* Auxiliary routine: printing a complex matrix */
void print_zmatrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )
                        printf( " (%6.2f,%6.2f)", a[i*lda+j].real, a[i*lda+j].imag );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}

void printvec(std::vector<double> &vec){
 int i;
 double norm;
 norm = 0.0;
 for(i=0;i<vec.size();i++){
    norm += vec[i]*vec[i];
 }
 std::cout<<"Norm = "<<norm<<std::endl;
}


