#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "jacobiMethod.h"
#include <algorithm>
// #include "jacobiMethod.cpp"

TEST_CASE("Testing maxoffdiag of A(i,j)"){
  int n = 5;
  mat A;
  int x = 0; int y = 0;
  while(x == y){
    x = rand() % n;
    y = rand() % n;
  }
  float val = 1000000.0;
  vec r = zeros<vec>(n);
  A = makeTridiagmat(n,8,1,false,r);
  A(x,y) = val;
  A(y,x) = val;

  int k = 0; int l = 0;
  maxoffdiag(A, n, &k, &l);

  REQUIRE(k==x);
  REQUIRE(l==y);
  REQUIRE(A(k,l)==Approx(val));
}
TEST_CASE("Testing eigenvalues from jacobi"){
  int n = 5;
  mat A = zeros<mat>(n,n);
  mat Acopy = zeros<mat>(n,n);
  mat v = eye<mat>(n,n);
  vec r = zeros<vec>(n);
  A = makeTridiagmat(n,8,1,false,r);
  Acopy = makeTridiagmat(n,8,1,false,r);

  jacobi(n,A,v);

  vec B;
  B = getEigenvalues(Acopy, n);

  vec eigval;
  eig_sym(eigval, Acopy);

  vec Avals;
  Avals = getEigenvalues(A, n);

  // sort(Avals.begin(),Avals.end());

  for(int i=0; i < n; i++){
    REQUIRE(Avals(i) == Approx(eigval(i)));
  }
}
TEST_CASE("Rotate preserves frobenius norm"){
  int n = 5;
  mat A = zeros<mat>(n,n);
  mat v = eye<mat>(n,n);
  vec r = zeros<vec>(n);
  A = makeTridiagmat(n,8,1,false,r);
  double norm1;
  norm1 = frobeniusNorm(A,n);

  int x = 0; int y = 0;
  while(x == y){
    x = rand() % n;
    y = rand() % n;
  }

  rotate(A, v, n, x, y);
  double norm2;
  norm2 = frobeniusNorm(A,n);

  REQUIRE(norm1 == norm2);
}
TEST_CASE("Eignvector orthoganality"){
  int n = 4;
  mat A = zeros<mat>(n,n);
  mat v = eye<mat>(n,n);
  vec r = zeros<vec>(n);
  A = makeTridiagmat(n,8,1,false,r);

  jacobi(n, A, v);
  REQUIRE(dot(v.col(0),v.col(1)) == Approx(0));
  REQUIRE(dot(v.col(1),v.col(2)) == Approx(0));
  REQUIRE(dot(v.col(2),v.col(3)) == Approx(0));
  REQUIRE(dot(v.col(3),v.col(0)) == Approx(0));

  // REQUIRE(dot(v(0),v(1)) == 0);
}
