#include <iostream>
#include <armadillo>
#include <cmath>
// #include <ctime>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "jacobiMethod.h"

using namespace std;
using namespace arma;

int main(int argc, char* argv[]){
  // char* name = argv[1];
  int n = atoi(argv[1]);
  mat A = zeros<mat>(n,n);
  mat Acopy = zeros<mat>(n,n);
  A = makeTridiagmat(n);
  Acopy = makeTridiagmat(n);
  jacobi(n,A);
  A.print("A: ");
  vec B;
  B = getEigenvalues(A, n);
  B.print();

  vec eigval;
  eig_sym(eigval, Acopy);
  eigval.print("eigval: ");
}
