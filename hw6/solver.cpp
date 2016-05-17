#include "solver.h"

void solver::jacobi() {
    // Use the Jacobi method to solve the linear system.    
    double *next_x = new double[_dim]();
    double conv = 0.0;
    for(int i = 0; i < _dim; ++i) 
        next_x[i] = _x[i];
    
    for (int iter = 0; iter < _maxIter; ++iter) {        
        conv = 0.0;
        for(int i = 0; i < _dim; ++i) {
            next_x[i] = (_b[i] - sigma_mul(i, 0, i-1, _x) - sigma_mul(i, i+1, _dim, _x))/_A.get_data(i, i);
        }
        for(int i = 0; i < _dim; ++i) {
            conv += pow(next_x[i] - _x[i], 2);
            _x[i] = next_x[i];
        }
        conv = sqrt(conv);
        cout << "conv = " << conv << endl;
        if (conv < _epsilon) {
            //this method has converged.
            cout << "iter = " << iter << endl;
            break;
        }        
    }

    delete[] next_x;
}

void solver::gauss_seidel() {
    // Use the Gauss-Seidel method to solve the linear system.    
    double *next_x = new double[_dim]();
    double conv = 0.0;
    for(int i = 0; i < _dim; ++i) 
        next_x[i] = _x[i];

    for (int iter = 0; iter < _maxIter; ++iter) {        
        conv = 0.0;
        for(int i = 0; i < _dim; ++i) {
            next_x[i] = (_b[i] - sigma_mul(i, 0, i-1, next_x) - sigma_mul(i, i+1, _dim, _x))/_A.get_data(i, i);
        }
        for(int i = 0; i < _dim; ++i) {
            conv += pow(next_x[i] - _x[i], 2);
            _x[i] = next_x[i];
        }
        conv = sqrt(conv);
        cout << "conv = " << conv << endl;
        if (conv < _epsilon) {
            //this method has converged.
            cout << "iter = " << iter << endl;
            break;
        }        
    }

    delete[] next_x;
}

void solver::conjugate_gradient() {
    //TODO
}


double solver::sigma_mul(int i, int begin, int end, double *vec) {
    double result = 0.0;
    for(int j = begin; j < end; ++j) {
        result += _A.get_data(i, j)*vec[j];
    }
    return result;
} 

void solver::print_x() {
    for(int i = 0; i < _dim; ++i)
        cout << _x[i] << endl;//
}

void help(char *argv[], double maxIter, double default_x, string m) {
    cout << "Usage: " << argv[0] << "[-f {filename}] [-i {maxIter}] [-x {initial value of x}] [-m {ja | gs | cg}]" << endl;
    cout << "\tSolve the linear system Ax = b" << endl;
    cout << "\tA: matrix (given by -f argument)" << endl;
    cout << "\tx: the target vector" << endl;
    cout << "\tb: vector" << endl;
    cout << "Arguments: " << endl;
    cout << "\t-f: the filename of the matrix A" << endl;
    cout << "\t-i: the maximum number of iteration of the method, default value = " << maxIter << endl;
    cout << "\t-x: the initial value of x, default value = " << default_x << endl;
    cout << "\t-m: choose the method to solve, default method = " << m << endl;
    cout << "\t    supported method: ja | gs | cg for jacobi | gauss-seidel | conjugate gradient" << endl;
}

int main (int argc, char *argv[]) {
    string filename = "bcsstk14.mtx";
    string method = "ja";
    int maxIter = 50;
    double default_x = 50.0;    

    for(int i = 1; i < argc; i++) {
        char *arg = argv[i];
        switch(arg[0]) {
            case '-':
                switch(arg[1]) {
                    case 'f':
                        filename = argv[++i];
                        break;
                    case 'm':
                        method = argv[++i];
                        break;
                    case 'i':
                        maxIter = atoi(argv[++i]);
                        break;
                    case 'x':
                        default_x = atof(argv[++i]);
                        break;
                    default:
                        help(argv, maxIter, default_x, method);
                        return 0;
                }
                break;
            default:
                help(argv, maxIter, default_x, method);
                return 0;
        }
    }

    matrix A(filename);
    solver solMgr(A, maxIter);
    solMgr.set_x(500.0);
    if (method == "ja")
        solMgr.jacobi();
    else if (method == "gs")
        solMgr.gauss_seidel();
    else if (method == "cg")
        solMgr.conjugate_gradient();
    else
        cerr << "Flag error!" << endl;

    solMgr.print_x();

    return 0;
}
