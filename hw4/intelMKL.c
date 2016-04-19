#define min(x,y) (((x) < (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mkl.h"

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


int main()
{
    double *A, *B, *C;
    int m, n, k, i, j;
    double alpha, beta;

    clock_t begin, end;
    double time_spent;


    /*"\n This example computes real matrix C=alpha*A*B+beta*C using \n"
            " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
            " alpha and beta are double precision scalars\n\n");*/

    m = 8000, k = 8000, n = 8000;
    /* (" Initializing data for matrix multiplication C=A*B for matrix \n"
            " A(%ix%i) and matrix B(%ix%i)\n\n", m, k, k, n);*/
    alpha = 1.0; beta = 0.0;

    /*(" Allocating memory for matrices aligned on 64-byte boundary for better \n"
            " performance \n\n");*/
    A = (double *)mkl_malloc( m*k*sizeof( double ), 64 );
    B = (double *)mkl_malloc( k*n*sizeof( double ), 64 );
    C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );
    if (A == NULL || B == NULL || C == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      mkl_free(A);
      mkl_free(B);
      mkl_free(C);
      return 1;
    }

    /* (" Intializing matrix data \n\n");*/
    begin = clock();
    for (i = 0; i < (m*k); i++) {
        A[i] = (double)(i+1);
    }

    for (i = 0; i < (k*n); i++) {
        B[i] = (double)(i+1);
    }

    for (i = 0; i < (m*n); i++) {
        C[i] = 0.0;
    }
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The time spent in matrix initialization is %g secends\n", time_spent);

    /* (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");*/
    begin = clock();

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha, A, k, B, n, beta, C, n);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The time spent in matrix multiplication is %g secends\n", time_spent);
    /*printf ("\n Computations completed.\n\n");*/

    printf (" Top left corner of matrix A: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(k,6); j++) {
        printf ("%12.0f", A[j+i*k]);
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
    /* ("\n Deallocating memory \n\n");*/
    mkl_free(A);
    mkl_free(B);
    mkl_free(C);

    return 0;
}
