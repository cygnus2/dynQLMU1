CC=icc
CXX=icpc
LIBS=-lgsl -lgslcblas -lm -mkl -lmkl_rt -lmkl_intel_lp64 -lmkl_core 
#LIB_PATH=-L/usr/local/armadillo/8.400.0/lib64
#HEADER_PATH=-I/usr/local/armadillo/8.400.0/include
#for openmp
#LIBS=-lgsl -lgslcblas -larmadillo -lsuperlu -lmkl_rt -lmkl_intel_lp64 -lmkl_core -lpthread -lm -ldl -lmkl_intel_thread -liomp5
#OBJS=u1link.o construct.o common.o hamiltonian.o teval.o chconj.o chconj1.o measurePOL.o measureWLOOP.o trans.o diagTrans.o diagTransTriv.o
OBJS=u1link.o common.o construct.o hamiltonian.o evolveH.o calc_Oflip.o 
#teval.o chconj.o chconj1.o measurePOL.o measureWLOOP.o trans.o diagTrans.o diagTransTriv.o calcOP.o
CFLAGS= -std=c99 -Wall
#CXXFLAGS= -std=c++11 -Wall -qopenmp -DARMA_DONT_USE_WRAPPER -DARMA_USE_HDF5
CXXFLAGS= -std=c++11 -Wall 

all: 2du1link

2du1link: ${OBJS}
	${CXX} -o 2du1link ${OBJS} ${LIBS}
clean:
	rm -rf 2du1link ${OBJS} *.optrpt 
