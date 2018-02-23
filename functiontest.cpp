#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "jacobiMethod.h"
#include <algorithm>
// #include "jacobiMethod.cpp"

TEST_CASE("Testing max A(i,j)"){
  int n = 3;
  mat A;
  int x = 0; int y = 0;
  while(x == y){
    x = rand() % n;
    y = rand() % n;
  }
  float val = 1000000.0;
  A = makeTridiagmat(n);
  A(x,y) = val;
  A(y,x) = val;

  int k = 0; int l = 0;
  maxoffdiag(A, n, &k, &l);

  REQUIRE(k==x);
  REQUIRE(l==y);
  REQUIRE(A(k,l)==Approx(val));
}
TEST_CASE("Testing eigen values"){
  int n = 5;
  mat A = zeros<mat>(n,n);
  mat Acopy = zeros<mat>(n,n);
  A = makeTridiagmat(n);
  Acopy = makeTridiagmat(n);
  jacobi(n,A);

  vec B;
  B = getEigenvalues(Acopy, n);

  vec eigval;
  eig_sym(eigval, Acopy);

  vec Avals;
  Avals = getEigenvalues(A, n);

  sort(Avals.begin(),Avals.end());

  for(int i=0; i < n; i++){
    REQUIRE(Avals(i) == Approx(eigval(i)));
  }
}
TEST_CASE("Rotate preserves frobenius norm"){
  int n = 5;
  mat A = zeros<mat>(n,n);
  A = makeTridiagmat(n);
  double norm1;
  norm1 = frobeniusNorm(A,n);

  int x = 0; int y = 0;
  while(x == y){
    x = rand() % n;
    y = rand() % n;
  }

  rotate(A, n, x, y);
  double norm2;
  norm2 = frobeniusNorm(A,n);

  REQUIRE(norm1 == norm2);
}
