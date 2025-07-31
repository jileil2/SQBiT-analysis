// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <algorithm> //std::for_each
#include <vector>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
Eigen::MatrixXd cppkp(Eigen::Map<Eigen::MatrixXd> a, Eigen::Map<Eigen::MatrixXd> b) {
  Eigen::MatrixXd res = kroneckerProduct(a,b);
  return res;
}


// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

#include <RcppArmadillo.h>
#include <RcppEigen.h>


// [[Rcpp::export]]
SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A,
                      Eigen::Map<Eigen::MatrixXd> B, 
                      int n_cores){
  
  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

using namespace Rcpp;

// [[Rcpp::export]]
double neg(double x) {
  // Test for x greater than zero
  if(x < 0) {
    // Return x
    return 1 ;
    // Otherwise
  } else {
    // Return negative x
    return 0;
  }
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double sum_cpp(NumericVector x) {
  // The size of x
  int n = x.size();
  // Initialize the result
  double result = 0;
  // Complete the loop specification
  for(int i = 0; i<n; i++) {
    // Add the next value
    result = result + x[i];
  }
  return result;
}


// [[Rcpp::export]]
double sum_rhotau(NumericVector x, const double tau) {
  // The size of x
  int n = x.size();
  // Initialize the result
  double result = 0;
  // Complete the loop specification
  for(int i = 0; i<n; i++) {
    // Add the next value
    result = result + x[i] * (tau - neg(x[i]));
  }
  return result;
}

// [[Rcpp::export]]
mat arma_dist(const mat& X){
  int n = X.n_rows;
  mat D(n, n, fill::zeros); // Allocate a matrix of dimension n x n
  for (int i = 0; i < n; i++) {
    for(int k = 0; k < i; k++){
      D(i, k) = sqrt(sum(pow(X.row(i) - X.row(k), 2)));
      D(k, i) = D(i, k);
    }
  }
  return D;
}