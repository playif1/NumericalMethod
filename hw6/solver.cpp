#include "solver.h"

void solver::jacobi() {
    // Use the Jacobi method to solve the linear system.    
    vector<double> next_x(_x);
    double conv = 0.0;
    
    for (int iter = 0; iter < _maxIter; ++iter) {        
        conv = 0.0;
        for(int i = 0; i < _x.size(); ++i) {
            next_x[i] = (_b[i] - sigma_mul(i, 0, i-1, _A, _x) - sigma_mul(i, i+1, _dim, _A, _x))/_A.get_data(i, i);
        }
        for(int i = 0; i < _x.size(); ++i) {
            conv += pow(next_x[i] - _x[i], 2);
            _x[i] = next_x[i];
        }
        conv = sqrt(conv);
        cout << "Iteration: " << iter << ", conv = " << conv << endl;
        if (conv < _epsilon) {
            //this method has converged.
            cout << "total iteration times = " << iter << endl;
            break;
        }        
    }
}

void solver::gauss_seidel() {
    // Use the Gauss-Seidel method to solve the linear system.    
    vector<double> next_x(_x);
    double conv = 0.0;

    for (int iter = 0; iter < _maxIter; ++iter) {        
        conv = 0.0;
        for(int i = 0; i < _x.size(); ++i) {
            next_x[i] = (_b[i] - sigma_mul(i, 0, i-1, _A, next_x) - sigma_mul(i, i+1, _dim, _A, _x))/_A.get_data(i, i);
        }
        
        for(int i = 0; i < _x.size(); ++i) {
            conv += pow(next_x[i] - _x[i], 2);
            _x[i] = next_x[i];
        }
        conv = sqrt(conv);
        cout << "Iteration: " << iter << ", conv = " << conv << endl;
        if (conv < _epsilon) {
            //this method has converged.
            cout << "total iteration times = " << iter << endl;
            break;
        }        
    }
}

void solver::conjugate_gradient(bool diag = false) {
    // Use the Conjugate Gradient method to solve the linear system.    
    vector<double> tmp_b = _b;
    if (diag) {
        diag_scaling();        
        vector<double> eye(_b.size(), 1.0);
        for(int i = 0; i < _b.size(); ++i) {
            _b[i] = sigma_mul(i, 0, _dim, _C, eye);
        }
    }
    set_x(0.0);
    int iter = 0;
    vector<double> r(_b);
    vector<double> p(r.size(), 0.0);
    vector<double> w(r.size(), 0.0);
    double prev_rho = vec_dot(r, r);
    double rho = vec_dot(r, r);
    double alpha = 0.0;
    double beta = 0.0;
    double norm_b = sqrt(vec_dot(_b, _b));
    cout << "norm b = " << norm_b << endl;   

    for (int iter = 0; iter < _maxIter; ++iter) {
        //sqrt(rho) > _epsilon*norm_b
        cout << "Iteration: " << iter << ", sqrt(rho) = " << sqrt(rho) << endl;
        if (sqrt(rho) < _epsilon * norm_b) {
            break;
        }

        beta = rho/prev_rho;
        for(int i = 0; i < p.size(); ++i) {
            p[i] = r[i] + beta*p[i];
        }   
        
        // compute w
        if(!diag) {
            for(int i = 0; i < w.size(); ++i) {
                w[i] = sigma_mul(i, 0, p.size(), _A, p);
            }
        }
        else {
            vector<double> cp(w.size(), 0.0);
            vector<double> acp(w.size(), 0.0);
            for(int i = 0; i < p.size(); ++i) {
                cp[i] = sigma_mul(i, 0, p.size(), _C, p);
            }
            for(int i = 0; i < cp.size(); ++i) {
                acp[i] = sigma_mul(i, 0, cp.size(), _A, cp);
            }
            for(int i = 0; i < w.size(); ++i) {
                w[i] = sigma_mul(i, 0, acp.size(), _C, acp);
            }
        }

        // compute alpha, then update x and r.
        alpha = vec_dot(r, r) / vec_dot(p, w);
        for(int i = 0; i < _x.size(); ++i) {
            _x[i] += alpha * p[i];
        }
        for(int i = 0; i < r.size(); ++i) {
            r[i] = r[i] - alpha * w[i];
        }

        prev_rho = rho;
        rho = vec_dot(r, r);
    }

    if (diag) {
        vector<double> tmp(_x.size(), 0.0);
        for(int i = 0; i < _x.size(); ++i) {
            tmp[i] = sigma_mul(i, 0, _x.size(), _C, _x);
        }
        _x = tmp;
        _b = tmp_b;
    }
}


double solver::sigma_mul(int i, int begin, int end, matrix& mat, vector<double>& vec) {
    double result = 0.0;
    for(int j = begin; j < end; ++j) {
        result += mat.get_data(i, j)*vec[j];
    }
    return result;
} 

double solver::vec_dot(vector<double>& v1, vector<double>& v2) {
    double dot = 0.0;
    for(int i = 0; i < v1.size(); ++i) {
        dot += v1[i] * v2[i];
    }
    return dot;
}

void solver::print_x() {
    for(int i = 0; i < _dim; ++i)
        cout << _x[i] << endl;//
}

void solver::diag_scaling() {
    for(int i = 1; i <= _dim; ++i)
        _C.set_data(i, i, sqrt(_A.get_data(i, i)));
}

double solver::calculate_error() {
    double err = 0.0;
    vector<double> ax;
    for(int i = 0; i < _x.size(); ++i) {
        ax.push_back(sigma_mul(i, 0, _x.size(), _A, _x));
        cout << "ax = " << ax[i] << ", b = " << _b[i] << endl;
        err += pow(ax[i] - _b[i], 2);
    }
    err = sqrt(err);
    return err;
}

void help(char *argv[], double maxIter, double default_x, string m) {
    cout << "Usage: " << argv[0] << "[-f {filename}] [-i {maxIter}] [-x {initial value of x}] [-m {ja | gs | cg}] [-d]" << endl;
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
    cout << "\t-d: set the diagonal scaling flag, default is false." << endl;
}

int main (int argc, char *argv[]) {
    string filename = "bcsstk14.mtx";
    string method = "ja";
    int maxIter = 50;
    double default_x = 50.0;    
    bool diag = false;

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
                    case 'd':
                        diag = true;
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

    cout << "create the matrix from " << filename << endl;
    matrix A(filename);
    cout << "Initialize the solver with maxIter = " << maxIter << endl;
    solver solMgr(A, maxIter);
    cout << "Set initial x to " << default_x << endl;
    solMgr.set_x(default_x);
    cout << "Using -" << method << " to solve this linear system." << endl;
    if (method == "ja")
        solMgr.jacobi();
    else if (method == "gs")
        solMgr.gauss_seidel();
    else if (method == "cg")
        solMgr.conjugate_gradient(diag);
    else
        cerr << "Flag error!" << endl;

    cout << "The error is: " << solMgr.calculate_error() << endl;
    solMgr.print_x();

    return 0;
}
