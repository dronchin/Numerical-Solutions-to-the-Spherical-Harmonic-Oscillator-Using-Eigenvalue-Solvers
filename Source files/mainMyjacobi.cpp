#include <iostream>
#include <armadillo>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "jacobiMethod.h"

using namespace std;
using namespace arma;
ofstream ofile;


int main(int argc, char* argv[]){
  int n;
  cout << "Enter number of integration points n: ";
  cin >> n;
  double wr;
  cout << "Enter harmonic oscilator angular frequency omega: ";
  cin >> wr;
  double r_max;
  cout << "Eneter a Rho max: ";
  cin >> r_max;
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
  vec r = zeros<vec>(n);

  A = makeTridiagmat(n,r_max,wr,interact,r);
  Acopy = makeTridiagmat(n,r_max,wr,interact,r);
  clock_t start1, stop1;
  start1 = clock();
  jacobi(n,A,v);
  stop1 = clock();
  double timer1 = (double)(stop1 - start1)/(CLOCKS_PER_SEC);
  // A.print();
  vec B;
  B = getEigenvalues(A, n);
  cout << "My eigenvalues: " << endl;
  for(int i=0; i<5; i++){
    cout << B(i) << ", ";
  }
  cout << endl;

  vec eigval;
  mat eigenvec;
  clock_t start2,stop2;
  start2 = clock();
  eig_sym(eigval, eigenvec, Acopy);
  stop2 = clock();
  double timer2 = (double)(stop2 - start2)/(CLOCKS_PER_SEC );
  cout << "armadillo eigen values: " << endl;
  for(int i=0; i<5; i++){
    cout << eigval(i) << ", ";
  }
  cout << endl;

  cout << "My jacobi took: " << timer1 << "s" << endl;
  cout << "Armadillo took: " << timer2 << "s" << endl;
}
