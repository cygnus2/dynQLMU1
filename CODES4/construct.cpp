/* construct gauge inequivalent states of the U(1) model  
 * on square lattice with external static charges Q=1/2 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<algorithm>
#include<iterator>
#include "define.h"
//#include <algorithm> //std::count
// For some notes see below

void conststates_Q1(){
  extern void add(char*);
  extern void addATpos(char*,int);
  extern int checkGL(bool*,int);
  extern int checkGLQP1(bool*,int);
  extern int checkGLQM1(bool*,int);
  extern void storeconf(bool*,int*);
  int count,count1;
  int flagFLUX;
  int base=6;
  long long int n,NST;
  char str[VOL2];
  int linkset[VOL];

  bool conf[2*VOL];
  int flagGI,flagJUMP;
  int i,j;
  int p,pp,l,cc;
  int q,xf,yf,xb,yb;
  int q0x,q0y,qx,qy;
  int Q,e1,e2,e3,e4;
  char cG;

  count = 0; 
  // step 1 in CHARGESTATES.txt 
  NST = pow(6,VOL2-2);
  std::cout<<"NST = "<<NST<<std::endl;
  for(p=0;p<VOL2;p++) str[p]='0';
  str[VOL2]='\0';

  for(i=0; i<4; i++){ // states for locQP
  for(j=0; j<4; j++){ // states for locQM
  n=0;
  printf("i=%d, j=%d, states=%d\n",i,j,count);
  while(n<NST){ // states over the Q=0 sites
    
    // track the links initialized  	  
    for(l=0;l<VOL;l++) linkset[l]=0; 

    // step 2+3 in CHARGESTATES.txt
    // initialize the links around the positive charge    
    q=chk2lin[locQP];   
    xf=next[DIM+1][q]; xb=next[DIM-1][q]; 
    yf=next[DIM+2][q]; yb=next[DIM-2][q];
    // convert to 1d co-ordinates on the lattice configuration
    q0x=2*q; q0y=2*q+1; qx=2*xb; qy=2*yb+1;
    if(i==0){      
      str[locQP]='A';      
      conf[q0x]=true;  conf[q0y]=true; conf[qx]=false;  conf[qy]=true; 
      linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; }
    else if(i==1){        
      str[locQP]='B';	    
      conf[q0x]=true;  conf[q0y]=false; conf[qx]=false;  conf[qy]=false; 
      linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; }
    else if(i==2){
       str[locQP]='C';	    
       conf[q0x]=true;  conf[q0y]=true; conf[qx]=true;   conf[qy]=false; 
       linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; }
    else if(i==3){ 
       str[locQP]='D';     
       conf[q0x]=false;  conf[q0y]=true; conf[qx]=false;   conf[qy]=false; 
       linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; }

    // initialize the links around the negative charge
    q=chk2lin[locQM]; 
    xf=next[DIM+1][q]; xb=next[DIM-1][q]; 
    yf=next[DIM+2][q]; yb=next[DIM-2][q];
    // convert to 1d co-ordinates on the lattice configuration
    q0x=2*q; q0y=2*q+1; qx=2*xb; qy=2*yb+1;
    if(j==0){         
      str[locQM]='E';	    
      conf[q0x]=false;  conf[q0y]=false; conf[qx]=true;  conf[qy]=false; 
      linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; }
    else if(j==1){     
      str[locQM]='F';	    
      conf[q0x]=false;  conf[q0y]=true; conf[qx]=true;   conf[qy]=true; 
      linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; }
    else if(j==2){        
      str[locQM]='G';	    
      conf[q0x]=false;  conf[q0y]=false; conf[qx]=false;   conf[qy]=true; 
      linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; }
    else if(j==3){        
      str[locQM]='H';	    
      conf[q0x]=true;  conf[q0y]=false; conf[qx]=true;   conf[qy]=true; 
      linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; }

    // Fill the other sites with Q=0 Gauss' Law      
    flagGI=1; flagJUMP=0;
    for(l=0;l<VOL2;l++){	    
      if((l==locQP)||(l==locQM)) continue;	    
      cG=str[l]; q=chk2lin[l]; 
      xf=next[DIM+1][q]; xb=next[DIM-1][q]; 
      yf=next[DIM+2][q]; yb=next[DIM-2][q];
      // convert to 1d co-ordinates on the lattice configuration
      q0x=2*q; q0y=2*q+1; qx=2*xb; qy=2*yb+1;
      // from Fig 1 (left to right) in https://arxiv.org/pdf/1311.2459.pdf
      switch(cG){
        case '0' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=true;  conf[qy]=false; 
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '1' : conf[q0x]=true;  conf[q0y]=true;  conf[qx]=true;  conf[qy]=true;  
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '2' : conf[q0x]=false; conf[q0y]=false; conf[qx]=false; conf[qy]=false;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '3' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=true;  conf[qy]=false;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '4' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=false; conf[qy]=true; 
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '5' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=false; conf[qy]=true; 
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
      } //close switch
      // check Gauss' Law at this point on the odd sublattices
      if(linkset[yb]==4){ flagGI = checkGL(conf,yb);
         if(flagGI==0){ addATpos(str,l); 
		  if((l<locQP)&&(l<locQM))      n=n+pow(6,VOL2-3-l);
		  else if((l<locQP)||(l<locQM)) n=n+pow(6,VOL2-2-l);
		  else                          n=n+pow(6,VOL2-1-l); 
		  flagJUMP=1; break; }
      }
      if(linkset[xb]==4){ flagGI = checkGL(conf,xb);    
         if(flagGI==0){ addATpos(str,l);  
		  if((l<locQP)&&(l<locQM))      n=n+pow(6,VOL2-3-l);
		  else if((l<locQP)||(l<locQM)) n=n+pow(6,VOL2-2-l);
		  else                          n=n+pow(6,VOL2-1-l); 
		  flagJUMP=1; break; }
      }
      if(linkset[xf]==4){ flagGI = checkGL(conf,xf);
         if(flagGI==0){ addATpos(str,l);  
		  if((l<locQP)&&(l<locQM))      n=n+pow(6,VOL2-3-l);
		  else if((l<locQP)||(l<locQM)) n=n+pow(6,VOL2-2-l);
		  else                          n=n+pow(6,VOL2-1-l); 
		  flagJUMP=1; break; }
      }
      if(linkset[yf]==4){ flagGI = checkGL(conf,yf);
         if(flagGI==0){ addATpos(str,l);  
		  if((l<locQP)&&(l<locQM))      n=n+pow(6,VOL2-3-l);
		  else if((l<locQP)||(l<locQM)) n=n+pow(6,VOL2-2-l);
		  else                          n=n+pow(6,VOL2-1-l); 
		  flagJUMP=1; break; }
      }
    }//close for loop over l
    //if(flagGI){
        // two versions of checks; the more restrictive one checks
        // over all lattice sites, choosing the routine for Q=0,+,-
        //for(p=0;p<VOL;p++){
        //  pp=chk2lin[p]; 
	//  if(p==locQP)      flagGI = checkGLQP1(conf,pp);      
	//  else if(p==locQM) flagGI = checkGLQM1(conf,pp);      
	//  else flagGI = checkGL(conf,pp);
        //  if(flagGI==0) break;
	//}
	// the faster, and the less restrictive one, checks only
	// the other checkerboard sites, where the charges are not
	// placed, and only have Q=0 Gauss Law
        //for(p=VOL2;p<VOL;p++){
        //  pp=chk2lin[p]; 
	//  flagGI = checkGL(conf,pp);
        //  if(flagGI==0) break;
        //}
    //}// if flagGI closure 

  // store the state if gauge invariant
  if(flagGI) { storeconf(conf,&count);
       	  //printf("state=%s\n",str); 
  }   
  if(flagJUMP==0){ add(str); n++; }

  } // close while loop over n
  }} // close for loops over i and j 
       
  printf("No of gauge invariant states: %d, %d\n",count,(int)basis.size());
  NTOT = count;

  //printf("Going to delete any duplicates\n");
  // delete duplicates
  //std::sort( basis.begin(),basis.end() );
  //basis.erase( unique(basis.begin(),basis.end() ), basis.end() );
  //count1 = basis.size();
  //printf("No of gauge invariant states after removal of duplicates: %d\n",count1);
}

void conststates_Q2(){
  extern void add(char*);
  extern void addATpos(char*,int);
  extern int checkGL(bool*,int);
  extern int checkGLQP2(bool*,int);
  extern int checkGLQM2(bool*,int);
  extern void storeconf(bool*,int*);
  int count,count1;
  int flagFLUX;
  int base=6;
  long long int n,NST;
  char str[VOL2];
  int linkset[VOL];

  bool conf[2*VOL];
  int flagGI,flagJUMP;
  int i,j;
  int p,pp,l,cc;
  int q,xf,yf,xb,yb;
  int q0x,q0y,qx,qy;
  int Q,e1,e2,e3,e4;
  char cG;

  count = 0; 
  // step 1 in CHARGESTATES.txt 
  NST = pow(6,VOL2-2);
  std::cout<<"NST = "<<NST<<std::endl;
  for(p=0;p<VOL2;p++) str[p]='0';
  str[VOL2]='\0';

  n=0;
  while(n<NST){ // states over the Q=0 sites
    
    // track the links initialized  	  
    for(l=0;l<VOL;l++) linkset[l]=0; 

    // initialize the links around the positive charge    
    q=chk2lin[locQP];   
    xf=next[DIM+1][q]; xb=next[DIM-1][q]; 
    yf=next[DIM+2][q]; yb=next[DIM-2][q];
    // convert to 1d co-ordinates on the lattice configuration
    q0x=2*q; q0y=2*q+1; qx=2*xb; qy=2*yb+1;
    str[locQP]='X';      
    conf[q0x]=true;  conf[q0y]=true; conf[qx]=false;  conf[qy]=false; 
    linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; 

    // initialize the links around the negative charge
    q=chk2lin[locQM]; 
    xf=next[DIM+1][q]; xb=next[DIM-1][q]; 
    yf=next[DIM+2][q]; yb=next[DIM-2][q];
    // convert to 1d co-ordinates on the lattice configuration
    q0x=2*q; q0y=2*q+1; qx=2*xb; qy=2*yb+1;
    str[locQM]='Y';	    
    conf[q0x]=false;  conf[q0y]=false; conf[qx]=true;  conf[qy]=true; 
    linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; 

    // Fill the other sites with Q=0 Gauss' Law      
    flagGI=1; flagJUMP=0;
    for(l=0;l<VOL2;l++){	    
      if((l==locQP)||(l==locQM)) continue;	    
      cG=str[l]; q=chk2lin[l]; 
      xf=next[DIM+1][q]; xb=next[DIM-1][q]; 
      yf=next[DIM+2][q]; yb=next[DIM-2][q];
      // convert to 1d co-ordinates on the lattice configuration
      q0x=2*q; q0y=2*q+1; qx=2*xb; qy=2*yb+1;
      // from Fig 1 (left to right) in https://arxiv.org/pdf/1311.2459.pdf
      switch(cG){
        case '0' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=true;  conf[qy]=false; 
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '1' : conf[q0x]=true;  conf[q0y]=true;  conf[qx]=true;  conf[qy]=true;  
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '2' : conf[q0x]=false; conf[q0y]=false; conf[qx]=false; conf[qy]=false;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '3' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=true;  conf[qy]=false;
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '4' : conf[q0x]=true;  conf[q0y]=false; conf[qx]=false; conf[qy]=true; 
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
        case '5' : conf[q0x]=false; conf[q0y]=true;  conf[qx]=false; conf[qy]=true; 
            linkset[q]=linkset[q]+4; linkset[xf]++; linkset[xb]++; linkset[yf]++; linkset[yb]++; break;
      } //close switch
      // check Gauss' Law at this point on the odd sublattices
      if(linkset[yb]==4){ flagGI = checkGL(conf,yb);
         if(flagGI==0){ addATpos(str,l); 
		  if((l<locQP)&&(l<locQM))      n=n+pow(6,VOL2-3-l);
		  else if((l<locQP)||(l<locQM)) n=n+pow(6,VOL2-2-l);
		  else                          n=n+pow(6,VOL2-1-l); 
		  flagJUMP=1; break; }
      }
      if(linkset[xb]==4){ flagGI = checkGL(conf,xb);    
         if(flagGI==0){ addATpos(str,l);  
		  if((l<locQP)&&(l<locQM))      n=n+pow(6,VOL2-3-l);
		  else if((l<locQP)||(l<locQM)) n=n+pow(6,VOL2-2-l);
		  else                          n=n+pow(6,VOL2-1-l); 
		  flagJUMP=1; break; }
      }
      if(linkset[xf]==4){ flagGI = checkGL(conf,xf);
         if(flagGI==0){ addATpos(str,l);  
		  if((l<locQP)&&(l<locQM))      n=n+pow(6,VOL2-3-l);
		  else if((l<locQP)||(l<locQM)) n=n+pow(6,VOL2-2-l);
		  else                          n=n+pow(6,VOL2-1-l); 
		  flagJUMP=1; break; }
      }
      if(linkset[yf]==4){ flagGI = checkGL(conf,yf);
         if(flagGI==0){ addATpos(str,l);  
		  if((l<locQP)&&(l<locQM))      n=n+pow(6,VOL2-3-l);
		  else if((l<locQP)||(l<locQM)) n=n+pow(6,VOL2-2-l);
		  else                          n=n+pow(6,VOL2-1-l); 
		  flagJUMP=1; break; }
      }
    }//close for loop over l
    if(flagGI){
        // two versions of checks; the more restrictive one checks
        // over all lattice sites, choosing the routine for Q=0,+,-
        //for(p=0;p<VOL;p++){
        //  pp=chk2lin[p]; 
	//  if(p==locQP)      flagGI = checkGLQP1(conf,pp);      
	//  else if(p==locQM) flagGI = checkGLQM1(conf,pp);      
	//  else flagGI = checkGL(conf,pp);
        //  if(flagGI==0) break;
	//}
	// the faster, and the less restrictive one, checks only
	// the other checkerboard sites, where the charges are not
	// placed, and only have Q=0 Gauss Law
        //for(p=VOL2;p<VOL;p++){
        //  pp=chk2lin[p]; 
	//  flagGI = checkGL(conf,pp);
        //  if(flagGI==0) break;
        //}
    }// if flagGI closure 

  // store the state if gauge invariant
  if(flagGI) { storeconf(conf,&count);
       	  //printf("state=%s\n",str); 
  }   
  if(flagJUMP==0){ add(str); n++; }

  } // close while loop over n
       
  printf("No of gauge invariant states: %d, %d\n",count,(int)basis.size());
  NTOT = count;

  //printf("Going to delete any duplicates\n");
  // delete duplicates
  //std::sort( basis.begin(),basis.end() );
  //basis.erase( unique(basis.begin(),basis.end() ), basis.end() );
  //count1 = basis.size();
  //printf("No of gauge invariant states after removal of duplicates: %d\n",count1);
}


int checkGL(bool *conf, int site){
  int e1,e2,e3,e4;
  int pp,qx,qy;
  qx=2*next[DIM-1][site]; 
  qy=2*next[DIM-2][site]+1; 
  pp=2*site;
  e1 = (conf[pp])? 1:-1; e2 = (conf[pp+1])? 1:-1; e3 = (conf[qx])? 1:-1; e4 = (conf[qy])? 1:-1;
  // if Q=e1+e2-e3-e4=0 return 1 (Gauss' Law OK) else return 0
  if(e1+e2-e3-e4) return 0;
  else return 1;
}

int checkGLQP1(bool *conf, int site){
  int e1,e2,e3,e4;
  int pp,qx,qy;
  qx=2*next[DIM-1][site]; 
  qy=2*next[DIM-2][site]+1; 
  pp=2*site;
  e1 = (conf[pp])? 1:-1; e2 = (conf[pp+1])? 1:-1; e3 = (conf[qx])? 1:-1; e4 = (conf[qy])? 1:-1;
  // if Q=e1+e2-e3-e4=1 return 1 (Gauss' Law OK) else return 0
  if((e1+e2-e3-e4)==1) return 0;
  else return 1;
}

int checkGLQM1(bool *conf, int site){
  int e1,e2,e3,e4;
  int pp,qx,qy;
  qx=2*next[DIM-1][site]; 
  qy=2*next[DIM-2][site]+1; 
  pp=2*site;
  e1 = (conf[pp])? 1:-1; e2 = (conf[pp+1])? 1:-1; e3 = (conf[qx])? 1:-1; e4 = (conf[qy])? 1:-1;
  // if Q=e1+e2-e3-e4=(-1) return 1 (Gauss' Law OK) else return 0
  if((e1+e2-e3-e4)==-1) return 0;
  else return 1;
}

int checkGLQP2(bool *conf, int site){
  int e1,e2,e3,e4;
  int pp,qx,qy;
  qx=2*next[DIM-1][site]; 
  qy=2*next[DIM-2][site]+1; 
  pp=2*site;
  e1 = (conf[pp])? 1:-1; e2 = (conf[pp+1])? 1:-1; e3 = (conf[qx])? 1:-1; e4 = (conf[qy])? 1:-1;
  // if Q=e1+e2-e3-e4=2 return 1 (Gauss' Law OK) else return 0
  if((e1+e2-e3-e4)==2) return 0;
  else return 1;
}

int checkGLQM2(bool *conf, int site){
  int e1,e2,e3,e4;
  int pp,qx,qy;
  qx=2*next[DIM-1][site]; 
  qy=2*next[DIM-2][site]+1; 
  pp=2*site;
  e1 = (conf[pp])? 1:-1; e2 = (conf[pp+1])? 1:-1; e3 = (conf[qx])? 1:-1; e4 = (conf[qy])? 1:-1;
  // if Q=e1+e2-e3-e4=0 return 1 (Gauss' Law OK) else return 0
  if((e1+e2-e3-e4)==(-2)) return 0;
  else return 1;
}


// Function to convert a given decimal number to a base 'base'
void fromDeci(long long int inputNum, int base, char* str){
    int index;
    index=VOL2;  // state looks reversed to fit neighbor indexing
    // (number)_10 --> (number)_(base) by repeatedly dividing it 
    // by base and taking remainder
    str[VOL2] = '\0';
    while (inputNum > 0){
        str[--index] = (char)(inputNum%base + '0');
        inputNum /= base;
    }
    while(index > 0){
        str[--index] = (char)('0');
    }
}

void add(char *str){
  int k,carry;
  carry=0; k=VOL2-1;
  if((k==locQP)||(k==locQM)) k--;
  if(str[k]<'5') { str[k]=str[k]+1; return; }
  else { str[k]='0'; carry=1; } 
  while(k>0){
    k--;
    if((k==locQP)||(k==locQM)) k--;
    if((carry==1)&&(str[k]=='5')) str[k]='0';
    else if(str[k]<'5') { str[k]=str[k]+1; carry=0; return; }
  }
}

void addATpos(char *str, int pos){
  int k,carry;
  carry=0;k=pos;
  if((k==locQP)||(k==locQM)) {printf("Error in state construction \n"); exit(0); }
  if(str[k]<'5') { str[k]=str[k]+1; return; }
  else { str[k]='0'; carry=1; } 
  while(k>0){
    k--;
    if((k==locQP)||(k==locQM)) k--;
    if((carry==1)&&(str[k]=='5')) str[k]='0';
    else if(str[k]<'5') { str[k]=str[k]+1; carry=0; return; }
  }
}

void storeconf(bool *conf,int* count){
  std::vector<bool> temp;
  temp.reserve(2*VOL);
  for(unsigned int i=0;i<2*VOL;i++) temp.push_back(conf[i]); 
  basis.push_back(temp);
  (*count)++;
  //printf("%d\n",(*count));
  //printf("current size = %d\n",(int)basis.size());
  temp.clear();
 }

/* 
 =======================================================================
 Implementation which help
 =======================================================================
 1. The state construction proceeds by application of the six possibilites
    of Gauss Law on all the even sites. These states can be represented in
    a base-6 notation. Eg: 00,01,02,...,05,10,11,...,55 can represent all
    the states of a 2x2 lattice. The sites are locally labelled using the
    checkerboard notation.
 2. The states use a bit-wise representation, ie, a true or false value 
    represents the {s^z} spin in the direction. 
 3. Normally one would have to have loop over 6^(VOL/2) possibilities and
    this will iterate over the possibilities represented as VOL/2 string-bits,
    with each bit between 0 and 1 (representing a particular Gauss' Law 
    realization on the even sites). The Gauss Law over the odd sites would
    now have to be separately checked, and only a config that satisfies all 
    the local Gauss' Laws would be accepted as a basis state.
 4. The algorithm improves the iterator (via a string addition method) 
 5. An important step in the algorithm is to identify the bottleneck in 
    iteration. For example consider a configuration on the 4x4 lattice, which 
    is represented as a string 23041533. The counting is as follows:
      string-array index:          0 1 2 3 4 5 6 7
      state-string      :          2 3 0 4 1 5 3 3
      checkerboard-label:          7 6 5 4 3 2 1 0

      (in the later versions, the checkerboard label is identical with the
       string-array index)

    It can happen that the forward x-neighbor of checkerboard site-label 4 
    (note checker-board site 4 = linear index 8, and (x,y) = (0,2)) does not 
    satisfy Gauss' Law. This is linear site 9 ==> (x,y) = (1,2). The standard 
    method would have rejected this configuration, then tried the next configuration, 
    represented as '2 3 0 4 1 5 3 4'. This however still has the GI violation
    at the aforementioned position (at string-index 3), and only after 
    6^(4) = 6^(VOL2-1-3) would it actually potentially be able to remove the 
    conflict. The current algorithm does it right away: it has a counter which
    measures how many of the links touching an odd site has been set, while
    it sets 4-links touching an even site using Gauss Law. When already
    4-links touch the odd neighbor (in +x, or -x, or +y or -y) then the
    string at the offending position is added, ie the configuration 
    '2 3 0 5 1 5 3 3' is checked. This already skips the 6^4 configurations
    which would have been checked without any effect. 
*/

/*  Note that this implementation does not avoid a "jam" that can happen
    at a different site. A more general algrithm is needed to untangle
    "jams" at multiple-vertices.
 */


