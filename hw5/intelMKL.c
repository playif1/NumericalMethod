#define min(x,y) (((x) < (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "mkl.h"
#include "mkl_spblas.h"

void printAns(int m, int k, int n, double *A, double *B, double *C) {
    int i, o, j;
    printf("A:\n");
    for (i = 0; i<m; ++i) {
        printf("[");
        for (j = 0; j<k; ++j) {
            printf ("%g ", A[i]);
        }
        printf("]\n");
    }

    printf("B:\n");
    for (i = 0; i<k; ++i) {
        printf("[");
        for (j = 0; j<n; ++j) {
            printf ("%g ", B[i]);
        }
        printf("]\n");
    }

    printf("C:\n");
    for (i = 0; i<m; ++i) {
        printf("[");
        for (j = 0; j<n; ++j) {
            printf ("%g ", C[i]);
        }
        printf("]\n");
    }
}

void help(char *argv[], MKL_INT n, double sparsity) {
    printf("Usage: %s [-n {matrix size}] [-s {sparsity}]\n",argv[0]);
    printf("\tC = A * B\n");
    printf("\tA: sparse matrix (csr format here)\n");
    printf("\tB: dense matrix\n");
    printf("\tC: sparse matrix (csr format here)\n");
    printf("Arguments: \n");
    printf("\t-n: the matrix size of these three matrix, default value = %d\n", n);
    printf("\t-s: the sparsity (the ratio of zero elements) of A, default value = %f\n" ,sparsity);
}

int main(int argc, char *argv[]) {
    MKL_INT m = 8000, k = 8000, n = 8000;
    double sparsity = 0.99;
    int i, j;
    for(i = 1; i < argc; i++) {
        char *arg = argv[i];
        switch(arg[0]) {
            case '-':
                switch(arg[1]) {
                    case 'n':
                        ++i;
                        m = atoi(argv[i]);
                        k = atoi(argv[i]);
                        n = atoi(argv[i]);
                        break;
                    case 's':
                        sparsity = atof(argv[++i]);
                        break;
                    default:
                        help(argv, n, sparsity);
                        return 0;
                }
                break;
            default:
                help(argv, n, sparsity);
                return 0;
        }
    }
    if (argc <= 1)
        printf("using default value, matrix size = %d, sparsity = %f\n", m, sparsity);
    else
        printf("using value, matrix size = %d, sparsity = %f\n", m, sparsity);

    clock_t begin, end;
    double time_spent;

    struct timespec start, finish;
    double elapsed;

    /*"\n This example computes real matrix C=alpha*A*B+beta*C using \n"
            " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
            " alpha and beta are double precision scalars\n\n");*/

    double *A_dense, *B, *C;
    double alpha = 1.0;
    double beta = 0.0;
    const int allignment = 64;

    /*(" Allocating memory for matrices aligned on 64-byte boundary for better \n"
            " performance \n\n");*/
    A_dense = (double *)mkl_malloc( m*k*sizeof( double ), allignment );
    B = (double *)mkl_malloc( k*n*sizeof( double ), allignment );
    C = (double *)mkl_malloc( m*n*sizeof( double ), allignment );
    if (A_dense == NULL || B == NULL || C == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      mkl_free(A_dense);
      mkl_free(B);
      mkl_free(C);
      return 1;
    }

    /* (" Intializing matrix data \n\n");*/
    int nzmax = 0;
    srand(time(NULL));
    begin = clock();
    for (i = 0; i < (m*k); i++) {
        double val = rand() / (double)RAND_MAX;
        if ( val < sparsity ) {
            A_dense[i] = 0.0;
        } else {
            A_dense[i] = val;
            nzmax++;
        }
    }

    for (i = 0; i < (k*n); i++) {
        B[i] = rand();
    }

    for (i = 0; i < (m*n); i++) {
        C[i] = 0.0;
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The CPU time spent in dense matrix initialization is %g seconds\n", time_spent);

    // Convert A to a sparse matrix in CSR format.
    MKL_INT job[6];
    job[0] = 0;  // convert TO CSR.
    job[1] = 0;  // Zero-based indexing for input.
    job[2] = 0;  // Zero-based indexing for output.
    job[3] = 2;  // adns is a whole matrix A.
    job[4] = nzmax;  // Maximum number of non-zero elements allowed.
    job[5] = 3;  // all 3 arays are generated for output. 
   

    /* JOB: conversion parameters
    * m: number of rows of A.
    * k: number of columns of A.
    * adns: (input/output). Array containing non-zero elements of the matrix A.
    * lda: specifies the leading dimension of adns. must be at least max(1, m).
    * acsr: (input/output) array containing non-zero elements of the matrix A. 
    * ja: array containing the column indices.
    * ia length m+1,  rowIndex.
    * OUTPUT:
    * info: 0 if successful. i if interrupted at i-th row because of lack of space.
    */
    int info = -1;
    printf("nzmax:\t %d\n", nzmax);

    double *A_sparse = mkl_malloc(nzmax * sizeof(double), allignment);
    if (A_sparse == NULL) {
        printf("ERROR: Could not allocate enough space to A_sparse.\n");
        return 1;
    }
    MKL_INT *A_sparse_cols = mkl_malloc(nzmax * sizeof(MKL_INT), allignment);
    if (A_sparse_cols == NULL) {
        printf("ERROR: Could not allocate enough space to A_sparse_cols.\n");
        return 1;
    }
    MKL_INT *A_sparse_rowInd = mkl_malloc((m+1) * sizeof(MKL_INT), allignment);
    if (A_sparse_rowInd == NULL) {
        printf("ERROR: Could not allocate enough space to A_sparse_rowInd.\n");
        return 1;
    }

    begin = clock();
    clock_gettime(CLOCK_MONOTONIC, &start);

    mkl_ddnscsr(job, &m, &k, A_dense, &k, A_sparse, A_sparse_cols, A_sparse_rowInd, &info);
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The CPU time spent in sparse matrix transformation is %g seconds\n", time_spent);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("The wall-clock execution time spent in sparse matrix transformation is %g seconds\n\n", elapsed);

    if(info != 0) {
        printf("WARNING: info=%d, expected 0.\n", info);
    }
    assert(info == 0);

    char transa = 'n';
    MKL_INT ldb = n, ldc=n;
    char matdescra[6] = {'g', 'l', 'n', 'c', 'x', 'x'}; 

    /* (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");*/
    printf("Note that the following computation uses the same matrix A and B, but the store scheme of A is different.\n");
    begin = clock();
    clock_gettime(CLOCK_MONOTONIC, &start);

    mkl_dcsrmm(&transa, &m, &n, &k, &alpha, matdescra, A_sparse, A_sparse_cols, 
            A_sparse_rowInd, &(A_sparse_rowInd[1]), B, &ldb, &beta, C, &ldc);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The CPU time spent in matrix multiplication under intelMKL sparse format is %g seconds\n", time_spent);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("The wall-clock execution time spent in matrix multiplication under intelMKL sparse format is %g seconds\n\n", elapsed);

    // dense operation
    begin = clock();

    clock_gettime(CLOCK_MONOTONIC, &start);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha, A_dense, k, B, n, beta, C, n);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The CPU time spent in matrix multiplication under intelMKL dense format is %g seconds\n", time_spent);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("The wall-clock execution time spent in matrix multiplication under intelMKL dense format is %g seconds\n\n", elapsed);

    /*printf ("\n Computations completed.\n\n");*/

    /*    
    printf (" Top left corner of matrix A: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(k,6); j++) {
        printf ("%12.0f", A_dense[j+i*k]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix B: \n");
    for (i=0; i<min(k,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.0f", B[j+i*n]);
      }
      printf ("\n");
    }
    
    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", C[j+i*n]);
      }
      printf ("\n");
    }
    */
    /* ("\n Deallocating memory \n\n");*/
    mkl_free(A_dense);
    mkl_free(A_sparse);
    mkl_free(A_sparse_cols);
    mkl_free(A_sparse_rowInd);
    mkl_free(B);
    mkl_free(C);

    return 0;
}
