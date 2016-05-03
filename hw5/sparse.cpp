#include "sparse.h"

csrMat::csrMat(int c, int r) {
    _col = c;
    _row = r;
    srand(time(NULL));
    int ptr = 0;
    _row_ptr.push_back(ptr);
    for(int i = 0; i < r; ++i) {
        for(int j = 0; j < c; ++j) {            
            // random initialize the sparse matrix with given sparsity
            if(rand() % 100 < 20) {
                ++ptr;
                _col_ind.push_back(j);
                _value.push_back((double) rand() / RAND_MAX);            
            }
        }
        _row_ptr.push_back(ptr);
    }
    _row_ptr.push_back(ptr+1);
}

void csrMat::printMat() {
    cout << "[" << endl;
    vector<int>::iterator it = _row_ptr.begin();
    for (int r = 0; r < _row; ++r) {
       int idx = *it;
       int num = *(it+1) - idx;
       for (int c = 0; c < _col; ++c) {
           if (c == _col_ind[idx] && num > 0) {
               cout << setw(4) << _value[idx] << " ";
               ++idx;
               --num;
           }
           else {
               cout << setw(4) << (double) 0 << " ";
           }
       }
       cout << endl;
       ++it;
    }
    cout << "]" << endl;
}

void csrMat::csr_mul_dense(csrMat& A, denseMat& B) {
    _row = A.get_row();
    _col = B.get_col();
    vector<int>::iterator rowPtr = A.get_row_ptr();
    vector<int>::iterator colInd = A.get_col_ind();
    vector<double>::iterator Avalue = A.get_value();
    double *Bvalue = B.get_value();
    int ptr = 0;
    for (int i = 0; i < _row; ++i) {
        _row_ptr.push_back(ptr);
        for (int j = 0; j < _col; ++j) {
            double v = 0;
            for (int k = rowPtr[i]; k < rowPtr[i+1]; ++k) {                
                v += Avalue[k]*Bvalue[colInd[k] * B.get_row() + j];
            }
            if (v <= 10e-9)
                continue;
            _col_ind.push_back(j);
            _value.push_back(v);
            ++ptr;
        }
    }
    _row_ptr.push_back(ptr);
}

int main () {
    int m = 8000;
    int k = 8000;
    int n = 8000;

    clock_t begin, end;
    double time_spent;

    struct timespec start, finish;
    double elapsed;

    /* Initialize the three matrix */
    begin = clock();

    csrMat A(m, k); 

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The CPU time spent in matrix initialization is %g seconds\n", time_spent);

    denseMat B(k, n);
    csrMat C;

    /* Compute C = A B */
    begin = clock();
    clock_gettime(CLOCK_MONOTONIC, &start);

    C.csr_mul_dense(A, B);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The CPU time spent in matrix multiplication is %g seconds\n", time_spent);

    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("The wall-clock execution time spent in matrix multiplication is %g seconds\n\n", elapsed);

    /*print the answer if needed*/
    /*
    cout << endl << "A = " << endl;
    A.printMat();
  
    cout << endl << "B = " << endl;
    B.printMat();
  
    cout << endl << "C = " << endl;
    C.printMat();
    */
    return 0;  
}

