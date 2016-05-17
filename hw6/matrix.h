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
        double get_data(int r, int c) {return _data[pair<int, int>(r, c)]; }

    private:
        map<pair<int, int>, double> _data;
        int _rNum;
        int _cNum;
};
