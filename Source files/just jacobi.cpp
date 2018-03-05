int jacobi(int n, mat &A, mat& v){
  // double h = 1/ (double) n;
  // double hh = h*h;
  double max;
  double eps = 1e-8;
  double maxiterations = (double)n*(double)n*(double)n;
  int k;
  int l;
  k = 0;
  l = 0;


  max = maxoffdiag(A,n,&k,&l);
  int iterations = 0;
  while(fabs(max) > eps && (double) iterations < maxiterations){
    max = maxoffdiag(A,n,&k,&l);
    rotate(A,v,n,k,l);
    iterations++;
  }
  return 0;
}
