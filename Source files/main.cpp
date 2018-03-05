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

  vec eigval;
  mat eigenvec;
  eig_sym(eigval, eigenvec, A);
  for(int i=0; i<5; i++){
    cout << eigval(i) << ", ";
  }
  cout << endl;
  if(interact){
    ofile.open("eigenvecs_interact.txt");
    for(int i = 0; i < n; i++){
      for(int j=0; j<n; j++){
        ofile << setprecision(8) << eigenvec(i,j) << " ";
      }
      ofile << endl;
    }
    ofile.close();
  }else{
    ofile.open("eigenvecs_noninteract.txt");
    for(int i = 0; i < n; i++){
      for(int j=0; j<n; j++){
        ofile << setprecision(8) << eigenvec(i,j) << " ";
      }
      ofile << endl;
    }
    ofile.close();
  }

  ofile.open("rhovals.txt");
  for(int i = 0; i < n; i++){
    ofile << setprecision(8) << r(i) << " ";
    ofile << endl;
  }
  ofile.close();








}
