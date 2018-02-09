//g++ -o LUdecompsolver.x  LUdecompsolver.cpp -larmadillo -llapack -lblas
#include <iostream>
#include <armadillo>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace arma;
ofstream ofile;


double real(double x){return 1.0-(1-exp(-10))*x - exp(-10*x);}
double RelativeError(double calc, double real){return fabs((real-calc)/real);}

int LUdecomp(int n, char* name){

  double h = 1/ (double) n;
  mat A = zeros<mat>(n-1,n-1);
  double* x = new double[n+1];

  for(int i = 0; i <= n; i++){
    x[i] = i*h;
  }

  A(0,0) = -2; A(0,1) = 1;
  A(n-2,n-2) = -2; A(n-2,n-3) = 1;
  for(int i = 1; i < n-2; i++){
    A(i,i) = -2;
    A(i,i-1) = 1;
    A(i,i+1) = 1;
  }

  vec b(n-1);
  for(int i = 0; i < n-1; i++){
    b(i) = h*h*100*exp(-10*(i+1)*h);
  }

  clock_t start, stop;
  start = clock();

  vec ans = solve(A,-b);

  stop = clock();
  double timeused = (double) (stop - start)/(CLOCKS_PER_SEC );

  cout << "Time used for LUdecom: " << timeused << endl;

  ofile.open(name);
  for(int i = 1; i < n-1; i++){
    double exact = real(x[i+1]);
    double Error = RelativeError(ans(i),exact);
    ofile << setw(15) << setprecision(8) << x[i+1];
    ofile << setw(15) << setprecision(8) << ans(i);
    ofile << setw(15) << setprecision(8) << exact;
    ofile << setw(15) << setprecision(8) << log10(Error) << endl;
  }
  ofile.close();
  delete [] x;
}

int main(int argc, char* argv[]){
  char* name = argv[1];
  int n = atoi(argv[2]);
  LUdecomp(n, name);
}
