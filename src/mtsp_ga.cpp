

#include <RcppArmadillo.h>
#include "utility.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;



//' Runs the MTSP solver
//'
//' @param xy A \code{n} by \code{2} matrix of coordinate positions. Each row has the \code{x} and \code{y} position of the n sites to be visited.
//' @param dmat A \code{n} by \code{n} matrix with the distances between the positions \code{xy}. 
//' @param nSalesmen An integer giving the number of salespersons (must be >= 2).
//' @param minTour An integer giving the minimum tour length for any of the salesmen.
//' @param CostType Possible values are the integers 1 or 2; 1 minimises the total distance travelled and 2 minimises the maximum distance travelled across the salepersons.
//' @param popSize The population size / search space of the GA (should be divisible by 8 but this is corrected within program).
//' @param numIter is the desired number of iterations for the algorithm to run.
//' @param Epsilon A double giving a multiplier on blancing the costs beteen total distance and minmax (only used for \code{costType} 2).
//' @param return_all Logical on whether ro return full ouput. Defaults to \code{FALSE}.
//' @return Returns the solution to the MTSP. Depending on \code{return_all} may return the inputs and \code{best_tour}: returns
//'  a list where each element gives the sites visited by each salesperson, \code{minDist}: is the best cost found by the algorithm,
//'  \code{minDistTotal} gives tht total distance travelled for all salespersons, and \code{minMaxTotal} gives the maximum distance travelled by any one salesperson.
//' @examples
//' set.seed(1)
//' n = 25
//' nSalesmen = 5
//' xy = matrix(rnorm(n*2), ncol=2)
//' dmat = as.matrix(dist(xy))
//' minTour = 3
//' numIter = 10
//' popSize = 8
//' CostType = 2
//' res = mtsp_ga(xy, dmat, nSalesmen, minTour, CostType, popSize, numIter, return_all = TRUE)
//' @export
//' @useDynLib mtsp
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
List mtsp_ga(arma::mat xy, arma::mat dmat, int nSalesmen, int minTour, int CostType, int popSize, int numIter, double Epsilon=1e-10, bool return_all=false) {
  
  // Initializations for Route Break Point Selection
  int n = xy.n_rows;
  int nBreaks = nSalesmen-1;
  int dof = n - minTour*nSalesmen ;
  arma::vec addto = arma::ones(dof+1);
  
  for(int k = 1; k < nBreaks; ++k){
    addto = arma::cumsum(addto);
  }
  arma::vec cumProb = arma::cumsum(addto)/sum(addto);
  
  // Initialize the Populations
  arma::imat popRoute = zeros<arma::imat>(popSize, n);
  arma::imat popBreak = zeros<arma::imat>(popSize, nBreaks);
  popRoute.row(0) = linspace<arma::irowvec>(0.0, n-1, n);
  popBreak.row(0)=  rand_breaks(minTour, n, nBreaks, cumProb);
  for(int k = 1; k < popSize; ++k){
    Rcpp::IntegerVector rs = Rcpp::seq(0,n-1);
    popRoute.row(k) = as<arma::irowvec>(Rcpp::sample(rs, n));
    popBreak.row(k) = rand_breaks(minTour=minTour, n=n, nBreaks=nBreaks, cumProb=cumProb);
  } 
  
  
  // Run the GA
  double globalMin = std::numeric_limits<double>::infinity();
  double globalTotalDist = 0.0;
  double globalMaxDist = 0.0;
  arma::vec totalDist = arma::zeros(popSize);
  arma::vec actual_total_dist = arma::zeros(popSize);
  arma::vec actual_max_dist = arma::zeros(popSize);
  arma::vec distHistory = arma::zeros(numIter);
  arma::imat tmpPopRoute = zeros<arma::imat>(8, n);
  arma::imat tmpPopBreak = zeros<arma::imat>(8, nBreaks);
  arma::imat newPopRoute = zeros<arma::imat>(popSize, n);
  arma::imat newPopBreak = zeros<arma::imat>(popSize, nBreaks);
  
  // extra
  arma::imat optRoute(1, n);
  arma::imat optBreak(1, nBreaks);
  double minDist;
  double actual_total_dist_out;
  double actual_max_dist_out;
  
  
  for(int iter = 0; iter < numIter; ++iter){
    
    // Evaluate Members of the Population
    for(int p = 0; p < popSize; ++p){
      arma::vec d = arma::zeros(nSalesmen);
      arma::irowvec pRoute = popRoute.row(p);
      arma::irowvec pBreak = popBreak.row(p);
      arma::imat a = join_rows(zeros<arma::irowvec>(1), pBreak+1);
      arma::imat b = join_rows(pBreak, arma::imat(1,1).fill(n-1)); 
      arma::imat rng = join_rows(conv_to<icolvec>::from(a), conv_to<icolvec>::from(b));
    
      for(int s = 0; s < nSalesmen; ++s){
        int rr = rng(s, 1);  
        int cc = rng(s, 0);
        d(s) = dmat(pRoute(rr), pRoute(cc)); // return journey
        
        for(int k = cc; k < rr; ++k){ 
          d(s) += dmat(pRoute(k), pRoute(k+1));
        }
      } // end s loop
      
      // Costs
      if(CostType == 1) {
        
        // keep these first two for output
        actual_total_dist(p) = sum(d);
        actual_max_dist(p) = max(d);
        totalDist(p) = sum(d);
        
      } else if (CostType == 2) {
        
        actual_total_dist(p) = sum(d);
        actual_max_dist(p) = max(d);
        totalDist(p) = actual_max_dist(p) + Epsilon*actual_total_dist(p);
        
      } 
      
    } // end p loop
    
    // Find the Best Route in the Population
    minDist = totalDist.min();
    int index = totalDist.index_min();
    actual_total_dist_out = actual_total_dist(index);
    actual_max_dist_out = actual_max_dist(index);
    
    distHistory(iter) = minDist;
    if(minDist < globalMin){
      globalMin = minDist;
      globalTotalDist = actual_total_dist_out;
      globalMaxDist = actual_max_dist_out;
      optRoute = popRoute.row(index);
      optBreak = popBreak.row(index);
    }
    
    
    // Genetic Algorithm Operators
    arma::uvec randomOrder = arma::randperm(popSize);
    arma::ivec pseq = regspace<arma::ivec>(8, 8, popSize);
    for(auto p : pseq){
      arma::uvec s = regspace<arma::uvec>(p-8, 1, p-1);
      arma::imat rtes = popRoute.rows(randomOrder.rows(s));
      arma::imat brks = popBreak.rows(randomOrder.rows(s));
      arma::vec dists = totalDist(randomOrder.rows(s));
      int idx = dists.index_min();
      arma::irowvec bestOf8Route = rtes.row(idx);
      arma::irowvec bestOf8Break = brks.row(idx);
      Rcpp::NumericVector vals = ceiling(n*runif(2)) - 1; // as zero indexing
      arma::uvec routeInsertionPoints = as<arma::uvec>(vals);
      std::sort(routeInsertionPoints.begin(), routeInsertionPoints.end());
      
      //Generate New Solutions
      for(int k = 0; k < 8; ++k){
        tmpPopRoute.row(k) = bestOf8Route;
        tmpPopBreak.row(k) = bestOf8Break;
        
        switch(k){
        case 1 :
          tmpPopRoute = flip(tmpPopRoute, k, routeInsertionPoints);
          break;
        case 2 :
          tmpPopRoute = swap(tmpPopRoute, k, routeInsertionPoints);
          break;
        case 3 :
          tmpPopRoute = slide(tmpPopRoute, k, routeInsertionPoints);
          break;
        case 4 :
          tmpPopBreak.row(k) = rand_breaks(minTour=minTour, n=n, nBreaks=nBreaks, cumProb=cumProb);
          break;
        case 5 :
          tmpPopRoute = flip(tmpPopRoute, k, routeInsertionPoints);
          tmpPopBreak.row(k) = rand_breaks(minTour=minTour, n=n, nBreaks=nBreaks, cumProb=cumProb);
          break;
        case 6 :
          tmpPopRoute = swap(tmpPopRoute, k, routeInsertionPoints);
          tmpPopBreak.row(k) = rand_breaks(minTour=minTour, n=n, nBreaks=nBreaks, cumProb=cumProb);
          break;
        case 7 :
          tmpPopRoute = slide(tmpPopRoute, k, routeInsertionPoints);
          tmpPopBreak.row(k) = rand_breaks(minTour=minTour, n=n, nBreaks=nBreaks, cumProb=cumProb);
          break;
        } // end switch
        
      } // end k loop -- new solutions
      
      newPopRoute.rows(s) = tmpPopRoute;
      newPopBreak.rows(s) = tmpPopBreak;
      
    } // end p loop -- population
    
    popRoute = newPopRoute;
    popBreak = newPopBreak;
    
  } // end iter loop -- iterations
  
  
  // calculate outputs
  Rcpp::List Tour(nSalesmen);
  arma::imat temp_brk = join_rows(optBreak, arma::imat(1,1).fill(n-1));
  Tour(0) = optRoute.cols(0, temp_brk(0)) + 1;
  for(int i = 1; i < nSalesmen; ++i){
    Tour(i) = optRoute.cols(temp_brk(i-1)+1, temp_brk(i)) + 1;
  }
  
  
  
  
  // Gather results
  List L;
  if(return_all){
    L = List::create(
      _["xy"] = xy,
      _["dmat"] = dmat,
      _["nSalesmen"] = nSalesmen,
      _["popSize"] = popSize,
      _["numIter"] = numIter,
      _["optRoute"] = optRoute,
      _["optBreak"] = optBreak,
      _["best_tour"] = Tour,
      _["minDist"] = globalMin,
      _["minDistTotal"] = globalTotalDist,
      _["minMaxTotal"] = globalMaxDist,
      _["distHistory"] = distHistory
    );
  } else {
    L = List::create(
      _["best_tour"] = Tour,
      _["minDistTotal"] = globalTotalDist,
      _["minMaxTotal"] = globalMaxDist
    );    
  }
  
  L.attr("class") = "mtsp_ga";
  
  return L;
  
} 





