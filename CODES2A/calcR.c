#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#define numstates 17906
//#define numstates 73789
#define numstates 1107
//#define numstates 141
int main()
{
  FILE *readthis;
  double eigval[(int)(numstates)];
  //double Oflip[(int)(numstates)];
  //double probstate[(int)(numstates)];
  int i,j;
  double MINGAP=1.0;
  //double beta;
  //double betastart;
  //double step_beta;
  //int numbeta;
  //double Zval;
  ///////////////////////////////////////////
  double ratio;
  double energydiff[(int)(numstates)-1];
  double mindelta;
  double maxdelta;
  double ratiodiff[(int)(numstates)-2];
  /////////////////
 
  /////////////////
  readthis=fopen("Eigenvalues.dat","r");

  for(i=0;i<(int)(numstates);i++)
    {
      fscanf(readthis,"%lf",&eigval[i]);
      //printf("%lf\n",eigval[i]);
    }
  fclose(readthis);
  /*************ratio calculation*******************/
  for(i=0;i<(int)(numstates)-1;i++)
    {
      energydiff[i]=eigval[i+1]-eigval[i];
      if(MINGAP > energydiff[i]) MINGAP = energydiff[i];
    }

  for(i=0;i<(int)(numstates)-2;i++)
    {
      if(energydiff[i+1]> energydiff[i])
    {
      mindelta=energydiff[i];
      maxdelta=energydiff[i+1];
    }
      else
    {
      mindelta=energydiff[i+1];
      maxdelta=energydiff[i];
    }
      ratiodiff[i]=mindelta/maxdelta;
      printf("%.12lf\n",ratiodiff[i]);
    }
  ratio=0.0;
  for(i=0;i<(int)(numstates)-2;i++)
    {
      ratio=ratio+ratiodiff[i];
    }
  ratio=ratio/(double)((int)(numstates)-2);
  printf("RATIO=%.12lf\n",ratio);
  printf("MinGAP = %.12lf\n",MINGAP);
  //getchar();
  /*************************************************/
 
 
  return 0;
}
