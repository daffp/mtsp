#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;



arma::irowvec rand_breaks(int minTour, int n, int nBreaks, arma::vec cumProb) {
  
  arma::irowvec breaks;
  
  if (minTour==1)   {
    
    arma::irowvec tmpBreaks = as<arma::irowvec>(Rcpp::sample(n-1, n-1)); 
    breaks = arma::sort(tmpBreaks.head(nBreaks)) - 1;
    
  } else {
    
    double rand_no = as<double>(runif(1));
    arma::uvec nAdjust = arma::find(rand_no < cumProb, 1, "first");
    
    int nA = as_scalar(nAdjust);
    Rcpp::NumericVector spaces = ceiling(nBreaks* runif(nA));
    
    arma::irowvec adjust = arma::zeros<arma::irowvec>(nBreaks);
    for(int kk = 0; kk < nBreaks; ++kk){
      adjust(kk) = sum(spaces == (kk+1));
    }
    
    arma::irowvec temp = linspace<arma::irowvec>(1.0, nBreaks, nBreaks); 
    breaks = minTour*temp + cumsum(adjust) - 1; 
    
  }
  
  return breaks;
}


arma::irowvec randbreak( int n, int max_salesmen) {
  int num_breaks = max_salesmen - 1;
  arma::irowvec breaks = arma::sort(as<arma::irowvec>(Rcpp::sample(n-1, num_breaks)));
  return breaks;
}


arma::imat CalcRange(arma::irowvec p_brk, int n){
  
  n = n - 1; //for zero indexing
  int flag = 1;
  arma::uword l = p_brk.n_elem;
  arma::uword lm = l - 1;
  arma::imat rng = zeros<arma::imat>(l+1, 2);
  
  for(int i = 0; i < l; ++i){
    
    if( (flag == 1) & (p_brk(i) > 0)){
      arma::irowvec v(2); v(0) = 0; v(1) = p_brk(i);
      rng.row(i) = v;
      flag = 0;
    }
    else if(flag == 1) {
      arma::irowvec v = {0,-1};
      rng.row(i) = v;
    }
    else if (p_brk(i) <= p_brk(i-1)){
      arma::irowvec v(2); v(0) = p_brk(i-1); v(1) = p_brk(i);
      rng.row(i) = v;
    }
    else if (i < lm){
      arma::irowvec v(2); v(0) = p_brk(i-1)+1; v(1) = p_brk(i);
      rng.row(i) = v;
    }   
    else{
      arma::irowvec v(2); v(0) = p_brk(i-1)+1; v(1) = p_brk(i);
      rng.row(i) = v;
    }
  }
  
  if( (p_brk(lm) < n) & (p_brk(lm) !=1)){
    arma::irowvec v(2); v(0) = p_brk(lm)+1; v(1) = n;
    rng.row(l) = v;
  }  
  else if( (p_brk(lm) < n) & (p_brk(lm) == 1)){
    arma::irowvec v(2); v(0) = p_brk(lm); v(1) = n;
    rng.row(l) = v;
  }
  else{
    arma::irowvec v(2); v(0) = p_brk(lm); v(1) = n-1;
    rng.row(l) = v;
  }
  
  return rng;
}


double CalcTourLength(arma::irowvec Tour, arma::mat d, arma::mat d0){
  
  int indices = Tour.n_elem - 1;
  double VehicleTourLength = d0(Tour(0),Tour(1));
  
  for(int c = 1; c < (indices-1); ++c){
    VehicleTourLength +=  d(Tour(c+1), Tour(c));
  }
  VehicleTourLength += d0(Tour(indices),Tour(indices-1));
  
  return VehicleTourLength;
}


// no depot
double CalcTourLength2(arma::irowvec Tour, arma::mat d){
  
  int indices = Tour.n_elem - 1;
  double VehicleTourLength = 0.0;
  
  for(int c = 0; c < indices; ++c){
    VehicleTourLength +=  d(Tour(c+1), Tour(c));
  }
  
  return VehicleTourLength;
}


arma::imat flip(arma::imat myMat, int r, arma::uvec ij) {
  
  int i = ij(0); int j = ij(1);
  myMat(r, span(i, j)) = arma::reverse( myMat(r, span(i, j)), 1);
  return myMat;
}

arma::imat swap(arma::imat myMat, int r, arma::uvec ij) {
  
  arma::uvec row_idx = r + arma::zeros<arma::uvec>(1);
  myMat.submat(row_idx, ij) = arma::reverse( myMat.submat(row_idx, ij), 1);
  return myMat;
}


arma::imat slide(arma::imat myMat, int r, arma::uvec ij) {
  
  int i = ij(0); int j = ij(1);
  arma::uvec v = regspace<arma::uvec>(i, 1, j);
  arma::uvec sv = v;
  std::rotate( sv.begin(), sv.begin() + 1, sv.end() );
  
  arma::uvec row_idx = r + arma::zeros<arma::uvec>(1);
  myMat.submat(row_idx, v) = myMat.submat(row_idx, sv);
  return myMat ;
}


