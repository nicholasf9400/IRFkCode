#include <Rcpp.h>
using namespace Rcpp;
using namespace sugar;

// [[Rcpp::export]]
double genCov(NumericVector h, NumericVector b){
  double a = sqrt(pow(h[0],2)+pow(h[1],2));
  double k = -b[0]*a;
  int blen = b.size();
  for(int p=1; p<blen; ++p){
    double c = pow(-1,p+1)*b[p]*pow(a,2*p+1);
    k += c;
  }
  return k;
}


// [[Rcpp::export]]
NumericMatrix GenCovMatrix(NumericMatrix points, NumericVector b){
  int n = points.nrow();
  NumericMatrix G = NumericMatrix(n,n);
  for(int i = 0; i<n; ++i){
    for(int j = i+1; j<n; ++j){
      NumericVector distance = points(j,_) - points(i,_);
      double g = genCov(distance, b);
      G(j,i) = g;
    }
  }
  return G;
}

// [[Rcpp::export]]
NumericVector GenCovVector(NumericMatrix points, NumericVector newpoint, NumericVector b){
  int n = points.nrow();
  NumericVector g = NumericMatrix(n,1);
  for(int i=0; i<n; ++i){
    NumericVector distance = newpoint - points(i,_);
    double q = genCov(distance, b);
    g[i] = q;
  }
  return g;
}
