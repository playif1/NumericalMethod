CC = g++
SPARSE = hw5_sparse
MKL = hw5_intelMKL
CFLAGS = -O3 -DADD_
DEBUG = -g -Wall
ICC = icc
MKLFLAG = -mkl

all: sparse intelmkl

sparse: 
	$(CC) -c $(CFLAGS) sparse.cpp
	$(CC) -o $(SPARSE) sparse.o
   
debug:
	$(CC) -c $(CFLAGS) $(DEBUG) sparse.cpp
	$(CC) -o $(SPARSE) sparse.o

intelmkl:
	$(ICC) -o $(MKL) $(CFLAGS) intelMKL.c $(MKLFLAG) 

clean: 
	rm -f *.o $(SPARSE) $(MKL)

