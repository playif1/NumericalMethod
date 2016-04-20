dlvl = ./.
include $(dlvl)/CBLAS/Makefile.in
CC = gcc
ATLAS = hw4_atlas
CBLAS = hw4_cblas
MKL = hw4_intelMKL
goto = hw4_gotoblas
ATLASFLAG = -lcblas
ICC = icc
MKLFLAG = -mkl
CFLAGS = -O3 -DADD_

all: atlas cblas intelmkl gotoblas

atlas:
	$(CC) -o $(ATLAS) $(CFLAGS) atlas.c $(ATLASFLAG)

cblas: 
	$(CC) -c $(CFLAGS) -I./CBLAS/src cblas.c 
	$(LOADER) -o $(CBLAS) cblas.o $(CBLIB) $(BLLIB) 
   
intelmkl:
	$(ICC) -o $(MKL) $(CFLAGS) intelMKL.c $(MKLFLAG) 

gotoblas:
	sudo $(CC) -o $(goto) $(CFLAGS) gotoBLAS.c -I./GotoBLAS2 -pthread -lgoto -lgfortran

cleanall: 
	rm -f *.o $(CBLAS) $(ATLAS) $(MKL) $(goto)

