#include "matrix.h"

matrix::matrix(int n) {
    // Initialize the identity matrix with dimension n
    for(int i = 1; i <= n; ++i)
        _data.insert(pair<pair<int, int>, double>(pair<int, int>(i, i), 1));
}

matrix::matrix(string filename) {
    // read the matrix according to the filename
    ifstream ifs;
    ifs.open (filename, ios::in);
    char line[256];
    size_t found = 0;
    size_t found2 = 0;    
    string rowstr;
    string colstr;
    string datastr;
    int row = 0;
    int col = 0;   
    double data = 0.0;

    ifs.getline(line, 256); // skip the first line
    ifs.getline(line, 256);
    string s(line);
    found = s.find_first_of(" ");
    found2 = s.find_first_of(" ", found+1);
    rowstr = s.substr(0, found);
    _rNum = stoi(rowstr); // The row number
    colstr = s.substr(found+1, found2-found);
    _cNum = stoi(colstr); // The col number

    while(ifs.getline(line, 256)) {      
        string s(line);
        found = s.find_first_of(" ");
        found2 = s.find_first_of(" ", found+1);
        rowstr = s.substr(0, found);
        row = stoi(rowstr); // The row
        colstr = s.substr(found+1, found2-found);
        col = stoi(colstr);
        datastr = s.substr(found2+1);
        data = stod(datastr);
        // Note that the index start from zero
        _data.insert(pair<pair<int, int>, double>(pair<int, int>(row-1, col-1), data));
    }
    ifs.close();
}

