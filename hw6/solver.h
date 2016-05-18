#include <vector>
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
             _epsilon = 10e-8;
             _dim = _A.get_col();
             _x.assign(_dim, 0);
             _b.assign(_dim, 1);
             //set_b2Ae();
        }

        void set_x(double x) {
            _x.assign(_x.size(), x);
        }

        void set_b(double b) {
            _b.assign(_b.size(), b);
        }

        void set_b2Ae() {
            vector<double> eye(_b.size(), 1.0);
            for(int i = 0; i < _b.size(); ++i)
                _b[i] = sigma_mul(i, 0, _dim, _A, eye);
        }

        double sigma_mul(int, int, int, matrix&, vector<double>&);
        double vec_norm(vector<double>&);
        double calculate_error();
        void diag_scaling();
        void print_x();
        void jacobi();
        void gauss_seidel();
        void conjugate_gradient(bool);

    private:
        matrix _A;
        matrix _C;
        vector<double> _x;
        vector<double> _b;
        double _epsilon;
        int _dim;
        int _maxIter;
};

