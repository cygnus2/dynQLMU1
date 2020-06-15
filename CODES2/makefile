CC=icc
CXX=icpc
LIBS=-lgsl -lgslcblas -lm -mkl -lmkl_rt -lmkl_intel_lp64 -lmkl_core 
OBJS=u1link.o common.o construct.o initialize.o symm_decompose.o diag.o hamiltonian.o Lecho.o evolveH_ov1.o evolveH_ov2.o evolveH_ov3.o evolveEent.o calc_Oflipt.o calc_Okint.o 
CFLAGS= -std=c99 -Wall
CXXFLAGS= -std=c++11 -Wall -mcmodel=large -shared-intel -DMKL_ILP64

all: 2du1link

2du1link: ${OBJS}
	${CXX} -o 2du1link ${OBJS} ${LIBS}
clean:
	rm -rf 2du1link ${OBJS} *.optrpt 
