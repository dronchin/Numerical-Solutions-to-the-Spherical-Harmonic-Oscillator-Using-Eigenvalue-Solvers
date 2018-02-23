#ifndef JACOBIMETHOD_H
#define	JACOBIMETHOD_H

#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<armadillo>
#include<time.h>
#include<vector>

using namespace std;
using namespace arma;

mat makeTridiagmat(int);
double maxoffdiag(mat,int,int*,int*);
int rotate(mat&,int,int,int);
vec getEigenvalues(mat,int);
int jacobi(int,mat&);
double frobeniusNorm(mat, int);

#endif /* JACOBIMETHOD_H */
