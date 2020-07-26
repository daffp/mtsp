#include <RcppArmadillo.h>

#ifndef MTSP_UTILITY_H
#define MTSP_UTILITY_H

arma::irowvec rand_breaks(int minTour, int n, int nBreaks, arma::vec cumProb);

arma::irowvec randbreak( int n, int max_salesmen);

arma::imat CalcRange(arma::irowvec p_brk, int n);

double CalcTourLength(arma::irowvec Tour, arma::mat d, arma::mat d0);

double CalcTourLength2(arma::irowvec Tour, arma::mat d);

arma::imat flip(arma::imat myMat, int r, arma::uvec ij);

arma::imat swap(arma::imat myMat, int r, arma::uvec ij) ;

arma::imat slide(arma::imat myMat, int r, arma::uvec ij) ;

#endif
