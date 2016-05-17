#include <cmath>
#include "matrix.h"

using namespace::std;

// Solve the linear system Ax=b
class solver {
    public:
        solver(){};
        solver(matrix A, int iter) {
             _maxIter = iter;
             _A = A; 
             _epsilon = 10e-6;
             _dim = _A.get_col();
             _x = new double[_dim]();
             _b = new double[_dim]();
             set_b2Ae();
        }

        void set_x(double x) {
            for(int i = 0; i < _dim; ++i) { _x[i] = x; }
        }

        void set_b(double b) {
            for(int i = 0; i < _dim; ++i) {_b[i] = b; }
        }

        void set_b2Ae() {
            double *eye = new double[_dim]();
            for(int i = 0; i < _dim; ++i)
                eye[i] = 1.0;
            for(int i = 0; i < _dim; ++i)
                _b[i] = sigma_mul(i, 0, _dim, eye);
            delete[] eye;
        }

        double sigma_mul(int, int, int, double*);
        void print_x();
        void jacobi();
        void gauss_seidel();
        void conjugate_gradient();

    private:
        matrix _A;
        double *_x;
        double *_b;
        double _epsilon;
        int _dim;
        int _maxIter;
};

