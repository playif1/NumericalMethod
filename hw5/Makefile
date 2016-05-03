CC = g++
SPARSE = hw5_sparse
CFLAGS = -O3 -DADD_

all: sparse

sparse: 
	$(CC) -c $(CFLAGS) sparse.cpp
	$(CC) -o $(SPARSE)  sparse.o
   
clean: 
	rm -f *.o $(SPARSE) 

