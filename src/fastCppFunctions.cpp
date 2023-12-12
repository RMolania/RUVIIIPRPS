#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
SEXP fastResidop(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd Y1 = B.transpose() * B ;
  Eigen::MatrixXd Y2 = Y1.inverse();
  Eigen::MatrixXd Y3 = B * Y2;
  Eigen::MatrixXd Y4 = B.transpose() * A;
  Eigen::MatrixXd Y5 = Y3 * Y4;
  Eigen::MatrixXd Y6 = A - Y5;
  return Rcpp::wrap(Y6);
}

// [[Rcpp::export]]
SEXP fastResidop2(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd Y1 = B.transpose() * B ;
  Eigen::MatrixXd Y2 = A - B * Y1.inverse() * B.transpose() * A;
  return Rcpp::wrap(Y2);
}

// [[Rcpp::export]]
SEXP matrixMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd Y = A * B ;
  return Rcpp::wrap(Y);
}

// [[Rcpp::export]]
SEXP matSubtraction(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd Y = A - B ;
  return Rcpp::wrap(Y);
}

// [[Rcpp::export]]
SEXP matTranspose(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::MatrixXd Y = A.transpose() ;
  return Rcpp::wrap(Y);
}

// [[Rcpp::export]]
SEXP matInverse(const Eigen::Map<Eigen::MatrixXd> A){
  Eigen::MatrixXd Y = A.inverse() ;
  return Rcpp::wrap(Y);
}

