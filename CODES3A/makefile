CC=icc
CXX=icpc
LIBS=-lgsl -lgslcblas -lm -mkl -lmkl_rt -lmkl_intel_lp64 -lmkl_core 
OBJS=u1link.o common.o construct.o symm_decompose.o symm_translation.o chconj.o initialize.o hamiltonian.o diag.o calc_Oflip.o entanglementEnt.o studyEvecs.o studyEvecs2.o
CFLAGS= -std=c99 -Wall
CXXFLAGS= -std=c++11 -Wall -mcmodel=large -shared-intel -DMKL_ILP64 
#CXXFLAGS= -std=c++11 -Wall -parallel -qopt-report-phase=par -qopt-report=5  -mcmodel=large -shared-intel -DMKL_ILP64 

all: 2du1link

2du1link: ${OBJS}
	${CXX} -o 2du1link ${OBJS} ${LIBS}
clean:
	rm -rf 2du1link ${OBJS} *.optrpt 
