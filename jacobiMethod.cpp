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

mat makeTridiagmat(int n, mat A){
  double h = 1.0/ (double) n;
  double hh = h*h;
  A(0,0) = 2.0/hh; A(0,1) = -1.0/hh;
  A(n-1,n-1) = 2.0/hh; A(n-1,n-2) = -1.0/hh;
  for(int i = 1; i < n-1; i++){
    A(i,i) = 2.0/hh;
    A(i,i-1) = -1.0/hh;
    A(i,i+1) = -1.0/hh;
  }
  return A;
}

double maxoffdiag(mat A, int n, int* k, int* l){
  double max = 0.0;
  for(int i=0; i<n;i++){
    for(int j=i+1; j<n; j++){
      // cout << i <<", " << j << ", " << A(i,j) << endl;
      if(fabs(A(i,j)) > max){
        max = fabs(A(i,j));
        cout << *k << endl;
        *k = i;
        *l = j;
      }
    }
  }
  return max;
}


int main(int argc, char* argv[]){
  // char* name = argv[1];
  int n = atoi(argv[1]);
  double h = 1/ (double) n;
  double hh = h*h;
  double max;

  mat A = zeros<mat>(n,n);
  A = makeTridiagmat(n, A);

  int k;
  int l;
  k = 0;
  l = 0;

  max = maxoffdiag(A,n,&k,&l);



}
