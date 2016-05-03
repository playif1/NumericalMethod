#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <ctime>
#include <vector>

using namespace::std;

class denseMat {
    public:
        denseMat();
        denseMat(int c, int r) { 
            srand(time(NULL));
            _col = c;
            _row = r;
            _value = new double[r*c];
            for(int i = 0; i < r; ++i) {
                for(int j = 0; j < c; ++j) {
                    _value[i*r+j] = (double) rand() / RAND_MAX;
                    //_value[i*r+j] = 1.0+j%2;
                }
            }
        };

        ~denseMat(){ delete [] _value; };

        int get_row() {return _row;}
        int get_col() {return _col;}
        double* get_value() {return _value;}

        void printMat() {
            cout << "[" << endl;
            for(int i = 0; i < _row; ++i) {
                for(int j = 0; j < _col; ++j) {
                    cout << " " <<  _value[i*_row + j];
                }
                cout << endl;
            }
            cout << "]" << endl;
        }

    private:
        int _row;
        int _col;
        double *_value;
};

class csrMat {
    public:
        csrMat() {};
        csrMat(int, int);

        int get_row() {return _row;}
        int get_col() {return _col;}
        vector<double>::iterator get_value() {return _value.begin();}
        vector<int>::iterator get_col_ind() {return _col_ind.begin();}
        vector<int>::iterator get_row_ptr() {return _row_ptr.begin();}

        void printMat();
        void csr_mul_dense(csrMat&, denseMat&);

    private:
        int _row;
        int _col;
        vector<double> _value;
        vector<int> _col_ind;
        vector<int> _row_ptr;
};

