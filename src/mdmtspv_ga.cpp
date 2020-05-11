

#include <RcppArmadillo.h>
#include "utility.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



//' Runs the MTSP solver
//'
//' @param xy A \code{n} by \code{2} matrix of coordinate positions. Each row has the \code{x} and \code{y} position of the n sites to be visited.
//' @param dmat A \code{n} by \code{n} matrix with the distances between the positions \code{xy}. 
//' @param depots A \code{nSales} by \code{2} matrix of coordinate positions for the starting positions / depots. This should have the same number of rows as there are salespersons.
//' @param nSalesmen An integer giving the number of salespersons (must be >= 2).
//' @param CostType Possible values are the integers 1 or 2; 1 minimises the total distance travelled and 2 minimises the maximum distance travelled across the salepersons.
//' @param popSize The population size / search space of the GA (should be divisible by 8 but this is corrected within program).
//' @param numIter is the number of desired iterations for the algorithm to run after a new best solution is found. 
//' @param Epsilon A double giving a multiplier on blancing the costs beteen total distance and minmax (only used for \code{costType} 2).
//' @param return_all Logical on whether ro return full ouput. Defaults to \code{FALSE}.
//' @return Returns the solution to the MTSP. Depending on \code{return_all} may return the inputs and \code{D0}: the distance from the depots to each site, \code{best_tour}: 
//' a matrix where each row gives the sites visited by each salesperson. The first and last values give the depot,  \code{minDist}: is the best cost found by the algorithm, 
//' and  \code{generation}: is the number of generations required by the algorithm to find the solution.
//' @examples
//' set.seed(1)
//' n = 25
//' nSalesmen = 5
//' xy = matrix(rnorm(n*2), ncol=2)
//' depots = matrix(rnorm(nSalesmen*2), ncol=2)
//' dmat = as.matrix(dist(xy))
//' mdmtspv_ga(xy, dmat, depots, nSalesmen, 1, 80, 10)
//' @export
//' @useDynLib mtsp
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
List mdmtspv_ga(arma::mat xy, arma::mat dmat, arma::mat depots,  int nSalesmen, int CostType, int popSize, int numIter, double Epsilon=1e-10, bool return_all=false) {

  // Initializations for Route Break Point Selection
  int n = xy.n_rows;
  int nBreaks = nSalesmen-1;

  arma::mat D0 = zeros<arma::mat>(nSalesmen, n);
  for(int i = 0; i < nSalesmen; ++i){
    for(int j = 0; j < n; ++j){
      D0(i,j) = arma::norm(depots.row(i) - xy.row(j)) ;
    }
  }

  // Initialize the Populations
  arma::imat popRoute = zeros<arma::imat>(popSize, n);
  arma::imat popBreak = zeros<arma::imat>(popSize, nBreaks);
  for(int k = 0; k < popSize; ++k){
    Rcpp::IntegerVector rs = Rcpp::seq(0,n-1);
    popRoute.row(k) = as<arma::irowvec>(Rcpp::sample(rs, n));
    popBreak.row(k) = randbreak( n, nSalesmen) ;
  }


  // Run the GA
  double globalMin = std::numeric_limits<double>::infinity();
  double globalTotalDist = 0.0;
  double globalMaxDist = 0.0;
  arma::vec totalDist = arma::zeros(popSize);
  arma::vec actual_total_dist = arma::zeros(popSize);
  arma::vec actual_max_dist = arma::zeros(popSize);
  arma::imat tmpPopRoute = zeros<arma::imat>(8, n);
  arma::imat tmpPopBreak = zeros<arma::imat>(8, nBreaks);
  arma::imat newPopRoute = zeros<arma::imat>(popSize, n);
  arma::imat newPopBreak = zeros<arma::imat>(popSize, nBreaks);

  // extra for scope
  arma::irowvec optRoute(1, n);
  arma::irowvec optBreak(1, nBreaks);
  arma::imat rng = zeros<arma::imat>(nSalesmen, 2);
  arma::imat rngKeep = zeros<arma::imat>(nSalesmen, 2);
  double minDist;
  double actual_total_dist_out;
  double actual_max_dist_out;
  
  int iter = -1, iter2go = -1, generation = 0;
  while(iter2go < (numIter-1)) {
    iter2go++; iter++;
    
    for(int p = 0; p < popSize; ++p){
      arma::vec d = arma::zeros(nSalesmen);
      arma::irowvec pRoute = popRoute.row(p);
      arma::irowvec pBreak = popBreak.row(p);
      rng = CalcRange(pBreak, n);

      for(int sa = 0; sa < nSalesmen; ++sa){
        
        int r0 = rng(sa, 0);
        int r1 = rng(sa, 1);
        if(r0 <= r1){
          
          // find a better way than this coercing
          arma::irowvec sap = arma::irowvec(1).fill(sa);
          arma::irowvec rng_range = pRoute.subvec(r0, r1);
          arma::irowvec Tour = join_rows(sap, rng_range, sap);  // returns to depot c(depot, breaks, depot)
          d(sa) = CalcTourLength(Tour, dmat, D0);
          
        } else {
          
          arma::irowvec sap = arma::irowvec(1).fill(sa);
          arma::irowvec Tour = join_rows(sap, sap); 
          d(sa) = 0.0;
          
        }
      } // sa; salesmen loop
      
      
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
      
    } // p: popsize loop
  
  
    // Find the Best Route in the Population
    minDist = totalDist.min();
    int index = totalDist.index_min();
    actual_total_dist_out = actual_total_dist(index);
    actual_max_dist_out = actual_max_dist(index);
    
    if(minDist < globalMin){
      iter2go = 0;
      generation = iter;
      globalMin = minDist;
      globalTotalDist = actual_total_dist_out;
      globalMaxDist = actual_max_dist_out;
      optRoute = popRoute.row(index);
      optBreak = popBreak.row(index);
      rngKeep = CalcRange(optBreak, n);
      }
  
  
    // Genetic Algorithm Operators
    //arma::uvec randomOrder = arma::randperm(popSize);
    Rcpp::IntegerVector s = Rcpp::sample(popSize, popSize)-1; 
    arma::uvec randomOrder = arma::sort(as<arma::uvec>(s));
    
    int ops = 8;
    arma::ivec pseq = regspace<arma::ivec>(ops, ops, popSize);
    for(auto p : pseq){
      arma::uvec s = regspace<arma::uvec>(p-ops, 1, p-1);
      arma::imat rtes = popRoute.rows(randomOrder.rows(s));
      arma::imat brks = popBreak.rows(randomOrder.rows(s));
      arma::vec dists = totalDist(randomOrder.rows(s));
      int idx = dists.index_min();
      arma::irowvec bestOf8Route = rtes.row(idx);
      arma::irowvec bestOf8Break = brks.row(idx);
      Rcpp::NumericVector vals = ceiling(n*runif(2)) - 1; // as zero indexing
      arma::uvec routeInsertionPoints = as<arma::uvec>(vals);
      std::sort(routeInsertionPoints.begin(), routeInsertionPoints.end()); // I,J

      //Generate New Solutions -- note that this doesn't include all of elad's
      for(int k = 0; k < ops; ++k){
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
            tmpPopBreak.row(k) = randbreak( n, nSalesmen);
            break;
          case 5 :
            tmpPopRoute = flip(tmpPopRoute, k, routeInsertionPoints);
            tmpPopBreak.row(k) = randbreak( n, nSalesmen);
            break;
          case 6 :
            tmpPopRoute = swap(tmpPopRoute, k, routeInsertionPoints);
            tmpPopBreak.row(k) = randbreak( n, nSalesmen);
            break;
          case 7 :
            tmpPopRoute = slide(tmpPopRoute, k, routeInsertionPoints);
            tmpPopBreak.row(k) = randbreak( n, nSalesmen);
            break;
          } // end switch

      } // end k loop -- new solutions
      
      newPopRoute.rows(s) = tmpPopRoute;
      newPopBreak.rows(s) = tmpPopBreak;

    } // end p loop -- population
  
     popRoute = newPopRoute;
     popBreak = newPopBreak;
  
  } // while loop

  
  // calculate outputs
  // no need to add the depot index as its givne by row/list index
  // if there are no values then it is NULL
  Rcpp::List Tour(nSalesmen);
  for(int i = 0; i < nSalesmen; ++i){
    int r0 = rngKeep(i, 0);
    int r1 = rngKeep(i, 1);
    if(r0 <= r1){
      Tour(i) = optRoute.subvec(r0, r1) + 1;
    } 
  }



  // Gather results
  List L;
  if(return_all){
                L = List::create(
                  _["xy"] = xy,
                  _["dmat"] = dmat,
                  _["depot"] = depots,
                  _["D0"] = D0,
                  _["nSalesmen"] = nSalesmen,
                  _["popSize"] = popSize,
                  _["numIter"] = numIter,
                  _["best_tour"] = Tour,
                  _["minDist"] = globalMin,
                  _["minDistTotal"] = globalTotalDist,
                  _["minMaxTotal"] = globalMaxDist,
                  _["generation"] = generation
                );
  } else {
    L = List::create(
      _["D0"] = D0,
      _["best_tour"] = Tour,
      _["minDistTotal"] = globalTotalDist,
      _["minMaxTotal"] = globalMaxDist
    );   
  }
  
  L.attr("class") = "mdmtspv_ga";
  
  return L;
  
}












