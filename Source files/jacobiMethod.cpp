#include "jacobiMethod.h"


using namespace std;
using namespace arma;

mat makeTridiagmat(int n, double r_max, double wr, bool interact, vec &r){
  // double r_max = 8.0;
  mat A = zeros<mat>(n,n);
  double h = r_max/ (double) n;
  double hh = h*h;
  // vec r = zeros<vec>(n);
  double r_min;
  if(interact){
    r_min = 0.0000001;
  }else{
    r_min = 0.0;
  }
  for(int i=0; i<n; i++){
    r(i) = r_min + i*h;
  }
  // r.print();
  // beginning of matrix
  if(interact){
    A(0,0) = 2.0/hh + wr*wr*r(0)*r(0) + 1.0/r(0); A(0,1) = -1.0/hh;
  }else{
    A(0,0) = 2.0/hh + wr*wr*r(0)*r(0); A(0,1) = -1.0/hh;
  }
  // rest of matrix
  for(int i = 1; i < n-1; i++){
    if(interact){
      A(i,i) = 2.0/hh + wr*wr*r(i)*r(i) + 1.0/r(i);
      A(i,i-1) = -1.0/hh; A(i,i+1) = -1.0/hh;
    }else{
      A(i,i) = 2.0/hh + wr*wr*r(i)*r(i);
      A(i,i-1) = -1.0/hh; A(i,i+1) = -1.0/hh;
    }
  }
  if(interact){
    A(n-1,n-1) = 2.0/hh + wr*wr*r(n-1)*r(n-1) + 1.0/r(0); A(n-1,n-2) = -1.0/hh;
  }else{
    A(n-1,n-1) = 2.0/hh + wr*wr*r(n-1)*r(n-1); A(n-1,n-2) = -1.0/hh;
  }
  // A.print();
  return A;
}

double maxoffdiag(mat A, int n, int* k, int* l){
  double max = 0.0;
  for(int i=0; i<n;i++){
    for(int j=i+1; j<n; j++){
      // cout << i <<", " << j << ", " << A(i,j) << endl;
      if(fabs(A(i,j)) > max){
        max = fabs(A(i,j));
        *k = i;
        *l = j;
      }
    }
  }
  return max;
}

int rotate(mat &A, mat &v, int n, int k, int l){
  double s, c;
  if(A(k,l) != 0){
    double t, tau;
    tau = (A(l,l) - A(k,k))/(2.0*A(k,l));
    if ( tau > 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    } else {
      t = -1.0/( -tau + sqrt(1.0 + tau*tau));
    }

  c = 1.0/sqrt(1+t*t);
  s = c*t;

  } else{
    s = 0;
    c = 1;
  }
  double a_kk, a_ll, a_ik, a_il, v_ik, v_il;
  a_kk = A(k,k);
  a_ll = A(l,l);
  A(k,k) = a_kk*c*c + a_ll*s*s - 2.0*c*s*A(k,l);
  A(l,l) = a_kk*s*s + a_ll*c*c + 2.0*c*s*A(k,l);
  A(k,l) = 0;
  A(l,k) = 0;
  for(int i = 0; i < n; i++){
    if(i != k && i!= l){
      a_ik = A(i,k);
      a_il = A(i,l);
      A(i,k) = c*a_ik - s*a_il;
      A(k,i) = A(i,k);
      A(i,l) = c*a_il + s*a_ik;
      A(l,i) = A(i,l);
    }
    v_ik = v(i,k);
    v_il = v(i,l);
    v(i,k) = c*v_ik - s*v_il;
    v(i,l) = c*v_il + s*v_ik;
  }
  return 0;
}

vec getEigenvalues(mat A, int n){
  vec Eigenvalues = zeros<vec>(n);
  for(int i=0; i<n; i++){
    Eigenvalues[i] = A(i,i);
  }
  sort(Eigenvalues.begin(),Eigenvalues.end());
  return Eigenvalues;
}

int jacobi(int n, mat &A, mat& v){
  // double h = 1/ (double) n;
  // double hh = h*h;
  double max;
  double eps = 1e-8;
  double maxiterations = (double)n*(double)n*(double)n;
  int k = 0;
  int l = 0;
  // k = 0;
  // l = 0;


  max = maxoffdiag(A,n,&k,&l);
  int iterations = 0;
  while(fabs(max) > eps && (double) iterations < maxiterations){
    max = maxoffdiag(A,n,&k,&l);
    rotate(A,v,n,k,l);
    iterations++;
  }
  return 0;
}

double frobeniusNorm(mat A, int n){
  double frobn = 0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      frobn += A(i,j)*A(i,j);
    }
  }
  return frobn;
}
