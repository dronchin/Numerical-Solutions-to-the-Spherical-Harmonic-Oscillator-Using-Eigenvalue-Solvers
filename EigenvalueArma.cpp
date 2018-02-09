//g++ -o LUdecompsolver.x  LUdecompsolver.cpp -larmadillo -llapack -lblas
#include <iostream>
#include <armadillo>
#include <cmath>
// #include <ctime>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;
using namespace arma;
ofstream ofile;

mat makeTridiag(int n, mat A){
  double h = 1/ (double) n;
  double hh = h*h;
  A(0,0) = 2/hh; A(0,1) = -1/hh;
  A(n-1,n-1) = 2/hh; A(n-1,n-2) = -1/hh;
  for(int i = 1; i < n-1; i++){
    A(i,i) = 2/hh;
    A(i,i-1) = -1/hh;
    A(i,i+1) = -1/hh;
  }
  return A;
}

int main(int argc, char* argv[]){
  // char* name = argv[1];
  int n = atoi(argv[1]);
  double h = 1/ (double) n;
  double hh = h*h;
  const double PI = 3.141592653589793;

  mat A = zeros<mat>(n,n);
  A = makeTridiag(n,A);
  A.print();

  vec eigval;
  mat eigvec;

  eig_sym(eigval, A);
  eigval.print();
  cout << endl;

  // double *lam = new double[n];
  vec lam = zeros<vec>(n);
  for(int i = 0; i < n; i++){
    lam[i] = 2/hh + 2*(-1)/hh*cos(((i+1)*PI)/(n+1));
  }
  lam.print();



}
