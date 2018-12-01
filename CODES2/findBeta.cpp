 /* code to find a beta for the average energy */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<algorithm>
#include<vector>
#include<iterator>

int main(){
 int Nstate = 282;
 double one,two,beta;
 double meanE,Z;
 int i;
 FILE *inp,*out;
 inp=fopen("try.dat","r");
 std::vector<double> ek;

 for(i=0; i<Nstate; i++){
   fscanf(inp,"%lf %lf",&one,&two);
   //std::cout<<one<<std::endl;
   ek.push_back(one);
 }
 fclose(inp);

 //for(i=0; i<Nstate; i++) printf("%d %lf\n",i,ek[i]);

 // scan in beta
 out=fopen("beta.dat","w");
 for(beta=-1.0; beta<=1.0; beta+=0.01){
   // calculate mean energy
   meanE=0.0; Z=0.0;
   for(i=0; i<Nstate; i++){
      Z += exp(-beta*ek[i]);
      meanE += ek[i]*exp(-beta*ek[i]);
   }
   meanE = meanE/Z;
   fprintf(out,"%lf %lf %lf\n",beta,Z,meanE);
 }
 fclose(out);
 
 return 0;
}
