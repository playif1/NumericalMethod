dlvl = ./.
include $(dlvl)/CBLAS/Makefile.in
ATLAS=hw4_atlas
CBLAS=hw4_cblas
MKL=hw4_intelMKL
ATLASFLAG=-lcblas
ICC=icc
MKLFLAG=-mkl

all: atlas cblas intelmkl

atlas:
	$(CC) -o $(ATLAS) $(CFLAGS) atlas.c $(ATLASFLAG)

cblas: 
	$(CC) -c $(CFLAGS) -I./CBLAS/src cblas.c 
	$(LOADER) -o $(CBLAS) cblas.o $(CBLIB) $(BLLIB) 
   
intelmkl:
	$(ICC) -o $(MKL) intelMKL.c $(MKLFLAG) 

cleanall: 
	rm -f *.o $(CBLAS) $(ATLAS) $(MKL)

