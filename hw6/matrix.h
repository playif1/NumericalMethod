#include <iostream>
#include <string>
#include <fstream>
#include <utility>
#include <map>

using namespace::std;

class matrix {
    public:
        matrix() {};
        matrix(int);
        matrix(string);
        int get_col() {return _cNum; }
        int get_row() {return _rNum; }
        double get_data(int r, int c) {
            map<pair<int, int>, double>::iterator it =  _data.find(pair<int, int>(r, c)); 
            return (it != _data.end())?it->second:0.0;
        }
        void set_data(int row, int col, double value) {
            _data.insert(pair<pair<int, int>, double>(pair<int, int>(row, col), value));
        }

    private:
        map<pair<int, int>, double> _data;
        int _rNum;
        int _cNum;
};
