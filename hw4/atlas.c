#include <cblas.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

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

    /* Initialize the three matrix */
    begin = clock();
    for(o = 0; o<m; o++) 
        for(i = 0; i<k; i++) 
            A[o*k + i] = 1;    

    for(o = 0; o<k; o++) 
        for(i = 0; i<n; i++) 
            B[o*n + i] = 1;    

    for(o = 0; o<m; o++) 
        for(i = 0; i<n; i++) 
            C[o*n + i] = 0;    

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The time spent in matrix initialization is %g secends\n", time_spent);

    /* Compute C = A B */
    begin = clock();
    cblas_dgemm (CblasRowMajor, 
                 CblasNoTrans, CblasNoTrans, m, n, k,
                 1.0, A, k, B, n, 0.0, C, n);
    
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The time spent in matrix multiplication is %g secends\n", time_spent);

    /*print the answer if needed*/
    printAns(m, k, n, A, B, C);
    free(A);
    free(B);
    free(C);

    return 0;  
}

