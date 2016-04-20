#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cblas.h"

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

int main () {
    int m = 8000;
    int k = 8000;
    int n = 8000;

    double *A = malloc(m*k*sizeof(double));
    double *B = malloc(k*n*sizeof(double));
    double *C = malloc(m*n*sizeof(double));
    int i, o, j;

    clock_t begin, end;
    double time_spent;

    struct timespec start, finish;
    double elapsed;

    /* Initialize the three matrix */
    begin = clock();
    for(o = 0; o<m; o++) 
        for(i = 0; i<k; i++) 
            A[o*k + i] = (double)(i+1);    

    for(o = 0; o<k; o++) 
        for(i = 0; i<n; i++) 
            B[o*n + i] = (double)(i+1);    

    for(o = 0; o<m; o++) 
        for(i = 0; i<n; i++) 
            C[o*n + i] = 0.0;    

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The CPU time spent in matrix initialization is %g seconds\n", time_spent);

    /* Compute C = A B */
    begin = clock();

    clock_gettime(CLOCK_MONOTONIC, &start);

    cblas_dgemm (CblasRowMajor, 
                 CblasNoTrans, CblasNoTrans, m, n, k,
                 1.0, A, k, B, n, 0.0, C, n);
    
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The CPU time spent in matrix multiplication is %g seconds\n", time_spent);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("The wall-clock execution time spent in matrix multiplication is %g seconds\n", elapsed);

    /*print the answer if needed*/
    //printAns(m, k, n, A, B, C);
    free(A);
    free(B);
    free(C);

    return 0;  
}

