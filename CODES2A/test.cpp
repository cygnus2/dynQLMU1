#include<stdio.h>
#include "ranlxd.h"

int main(){
  int SEED;
  printf("Enter an integer: \n");
  scanf("%d",&SEED);
  rlxd_init(1,SEED);
  printf("Random number initialized \n");

 }
