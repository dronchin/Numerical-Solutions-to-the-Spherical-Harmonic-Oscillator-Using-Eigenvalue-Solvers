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
  int n;
  cout << "Enter number of integration points n: ";
  cin >> n;
  double wr;
  cout << "Enter harmonic oscilator angular frequency omega: ";
  cin >> wr;
  bool interact;
  string yn;
  cout << "Do you want interacting particles? (y,n): ";
  cin >> yn;
  if(yn=="y" or yn=="Y"){
    interact = true;
  }else{
    interact = false;
  }

  mat A = zeros<mat>(n,n);
  mat Acopy = zeros<mat>(n,n);
  mat v = eye<mat>(n,n);

  double timer;
  A = makeTridiagmat(n,wr,true);
  Acopy = makeTridiagmat(n,wr,true);

  timer = jacobi(n,A,v);

  vec B;
  B = getEigenvalues(A, n);
  B.print();

  vec eigval;
  eig_sym(eigval, Acopy);
  eigval.print("eigval: ");
}
