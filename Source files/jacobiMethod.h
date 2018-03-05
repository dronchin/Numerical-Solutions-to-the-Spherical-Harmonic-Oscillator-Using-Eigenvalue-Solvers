#ifndef JACOBIMETHOD_H
#define	JACOBIMETHOD_H

#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<armadillo>
#include<time.h>
#include<vector>
#include <algorithm>


using namespace std;
using namespace arma;

mat makeTridiagmat(int,double,double,bool,vec&);
double maxoffdiag(mat,int,int*,int*);
int rotate(mat&,mat&,int,int,int);
vec getEigenvalues(mat,int);
int jacobi(int,mat&,mat&);
double frobeniusNorm(mat,int);

#endif /* JACOBIMETHOD_H */
