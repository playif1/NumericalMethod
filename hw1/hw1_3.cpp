#include <iostream>
#include <iomanip>
#include <cmath>

using namespace::std;

// Not used in this homework...
double round_cal(double val, double base = 10.0) {
   double p = int(log10(val));
   double shift = pow(base, p);
   return (round((val/shift)*100.0)/100.0)*shift;
}

int main() {
   double fi = 1;
   double A = -64;
   double B = 0;

   double cat_value = 1.0 - 1.0/(1+exp(A*fi+B));
   cout << "the catastrophic cancellation value = " << cat_value << endl;
   
   double reform_value  = exp(A*fi+B)/(1.0+exp(A*fi+B));
   cout << "the value after reforming the formula = " << reform_value << endl;

   return 0;
}

