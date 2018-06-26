CC=icc
CXX=icpc
LIBS=-lgsl -lgslcblas -lm -mkl -lmkl_rt -lmkl_intel_lp64 -lmkl_core 
OBJS=ranlxd.o u1link.o common.o construct.o symm_decompose.o hamiltonian.o diag.o chconj.o calc_Oflip.o 
CFLAGS= -std=c99 -Wall
CXXFLAGS= -std=c++11 -Wall -mcmodel=large -shared-intel -DMKL_ILP64 

all: 2du1link

2du1link: ${OBJS}
	${CXX} -o 2du1link ${OBJS} ${LIBS}
clean:
	rm -rf 2du1link ${OBJS} *.optrpt 
