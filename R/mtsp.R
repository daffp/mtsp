#' Runs the MTSP solver
#'
#' @param xy A \code{n} by \code{2} matrix of coordinate positions. Each row has the \code{x} and \code{y} position of the n sites to be visited.
#' @param depots A \code{nSales} by \code{2} matrix of coordinate positions for the starting positions / depots. This should have the same number of rows as there are salespersons.
#' @param nSalesmen An integer giving the number of salespersons (must be >= 2).
#' @param minTour An integer giving the minimum tour length for any of the salesmen.
#' @param CostType CostType Possible values are the integers 1 or 2; 1 minimises the total distance travelled and 2 minimises the maximum distance travelled across the salepersons. 
#' @param algorithm Which algorithm to call. Possible choices are \code{"mdmtspv_ga"}, \code{"mdmtspv_ga2"}, or \code{"mtsp_ga"}. See details.
#' @param popSize The population size / search space of the GA. Defaults to 80 (should be divisible by 8 but this is corrected within program).
#' @param numIter is the number of desired iterations for the algorithm to run after a new best solution is found. 
#' @param Epsilon A double giving a multiplier on blancing the costs beteen total distance and minmax (only used for \code{costType} 2).
#' @param return_all Logical on whether ro return full ouput. Defaults to \code{FALSE}.
#' @return Returns a list with the solution to the MTSP. Depending on \code{return_all} it may include the inputs and the following:
#' \itemize{
#'   \item D0 The distance from the depots to each site.
#'   \item best_tour A matrix where each row gives the sites visited by each salesperson. The first and last values give the depot.
#'   \item minDist The best cost found by the algorithm. 
#'   \item minDistTotal The total distance across the salespersons at the interation where the best cost was found by the algorithm. 
#'   \item generation The number of generations required by the algorithm to find the solution.
#'   }
#' @details The algorithm \code{mdmtspv_ga} searches for the best solution with each salesperson starting and ending at a depot. The algorithms \code{"mdmtspv_ga2"}
#'    and \code{"mtsp_ga"} start and end at the first site on the tour - they do not include depots.
#' @importFrom stats dist
#' @examples
#' set.seed(1)
#' n = 25
#' xy = matrix(rnorm(n*2), ncol=2)
#' run = mtsp(xy, nSalesmen=5, CostType=2, popSize=80, 
#'            numIter=10, algorithm="mtsp_ga", return_all=TRUE)
#' run
#' summary(run)
#' depots = matrix(runif(10), nc=2)
#' run2 = mtsp(xy, depots=depots, nSalesmen=5, CostType=2, popSize=80,
#'             numIter=10, algorithm="mdmtspv_ga", return_all=TRUE)
#' @export
mtsp = function(xy, depots=NULL, nSalesmen=2, minTour=3, CostType=2, algorithm="mtsp_ga", popSize=80, numIter=10, Epsilon=1e-10, return_all=FALSE) {
  
  if (!is.null(depots) && nSalesmen != nrow(depots)) stop('A depot is required for each salesperson')
  if (nSalesmen < 2) stop(' nSalesmen should be >= 2')
  if (!CostType %in% 1:2) stop('Unknown CostType')
  if (ncol(xy) != 2) stop('Position matrix xy has incorrect dimensions')
  if (!algorithm %in% c("mdmtspv_ga", "mdmtspv_ga2", "mtsp_ga")) stop('Unknown algorithm')
  
  n = nrow(xy)
  nSalesmen = max(1, min(n, nSalesmen))
  minTour = max(1, min(floor(n/nSalesmen), minTour))
  popSize = max(8, 8*ceiling(popSize/8)) 
  numIter = max(1, numIter)
  dmat = as.matrix(dist(xy))
  
  output = switch(algorithm,
                  "mtsp_ga" = mtsp_ga(xy, dmat, nSalesmen, minTour, CostType, popSize, numIter, Epsilon, return_all),
                  "mdmtspv_ga" = mdmtspv_ga(xy, dmat, depots, nSalesmen, CostType, popSize, numIter, Epsilon, return_all),
                  "mdmtspv_ga2" = mdmtspv_ga2(xy, dmat, nSalesmen, CostType, popSize, numIter, Epsilon, return_all) 
                  )
  
  class(output) = c("mtsp", class(output))
  
  return(output)
}
 

#' @method print mtsp
#' @export
print.mtsp = function(x, ...){
  cat(paste("Total distance travelled for all salespersons =", round(x$minDistTotal, 2)), sep="\n")  
  cat(paste("Maximum distance travelled by a single salesperson =", round(x$minMaxTotal, 2)), sep="\n")  
}


#' Grabs the distances for a single salesperson
#'
#' @param x The list returned by \code{mtsp}.
#' @param ... Other arguments.
#' @export get_distances
get_distances <- function(x, ...) UseMethod("get_distances")


#' Grabs the distances for a single salesperson
#'
#' @param x The list returned by \code{mtsp}.
#' @param idx An index representing a single list element of \code{best_tour}.
#' @param paths The list with the \code{best_tour} returned by \code{mtsp}.
#' @param depot_mat The matrix \code{D0} giving the distances between the depots and positions, returned by \code{mtsp}.
#' @param dmat The distance matrix for the positions.
#' @param ... Not used.
#' @return Returns the distance for a single salesperson.
get_distances.mdmtspv_ga = function( x, idx, paths, depot_mat, dmat, ...){
  # site distances
  dmat_idx = cbind(paths[-length(paths)], paths[-1])
  dmat_distances = dmat[dmat_idx]
  
  # then distance to depot
  depot_idx = cbind(idx, c(paths[1], paths[length(paths)]))
  depot_distances = depot_mat[depot_idx]
  
  salesperson_distance = sum(dmat_distances, depot_distances)
  return(salesperson_distance)
}


#' Grabs the distances for a single salesperson
#'
#' @param x The list returned by \code{mtsp}.
#' @param idx An index representing a single list element of \code{best_tour}.
#' @param paths The list with the \code{best_tour} returned by \code{mtsp}.
#' @param dmat The distance matrix for the positions.
#' @param ... Not used.
#' @return Returns the distance for a single salesperson.
get_distances.mtsp_ga = function( x, idx, paths, dmat, ...){
  # site distances
  dmat_idx = rbind(cbind(paths[-length(paths)], paths[-1]),c(paths[length(paths)], paths[1]))
  dmat_distances = dmat[dmat_idx]
  salesperson_distance = sum(dmat_distances)
  return(salesperson_distance)
}





#' Calculates the distance for each saleperson
#'
#' @param object The list returned by \code{mtsp}.
#' @param xy A \code{n} by \code{2} matrix of coordinate positions. Each row has the \code{x} and \code{y} position of the n sites to be visited.
#' @param ... Not used.
#' @return Returns the distance for each salesperson.
#' @method summary mtsp
#' @export
summary.mtsp = function(object, xy=NULL, ...){

      if(!inherits(object, "mtsp")) stop("Cannot run summary on this object")
      if(is.null(object$dmat) & is.null(xy)) stop("Must either run mtsp with return_all = TRUE or else provide the positions matrix")

      if(is.null(object$dmat) & !is.null(xy)){
                  nxy = nrow(xy)
                  ncD0 = ncol(object$D0)
                  if(nxy != ncD0) stop("data xy have different number of coordinates than were used in the mtsp call")
                  object$dmat = as.matrix(dist(xy))
                }

      # Get distance for each salesperson
      distances = sapply(seq_along(object$best_tour), function(i) get_distances(object, i, paths=object$best_tour[[i]], depot_mat = object$D0, dmat=object$dmat))

      ans = list(distances=distances, tours=object$best_tour)
      class(ans) <- c("summary.mtsp", class(object))
      return(ans)
      }


#' @export
print.summary.mtsp = function(x, ...){
  
      cat(paste("Total distance travelled:", round(sum(x$distances),2)), sep="\n")
      cat("================================", sep="\n")
      
      if(inherits(x, "mdmtspv_ga")){
        for( i in seq_along(x$tours)) {
          dep = paste0("depot", i)
          print_tours = paste(c(dep, x$tours[[i]], dep), collapse=" -> ")
          dist = x$distances[i]
          cat(paste("Salesperson", i, ": distance travelled =", round(dist, 2)), sep="\n")
          cat(paste("Tour:", print_tours), sep="\n")
          cat("\n")
          }
      }
     if(inherits(x, "mtsp_ga")){
        for( i in seq_along(x$tours)) {
          print_tours = paste(c(x$tours[[i]], x$tours[[i]][1]), collapse=" -> ")
          dist = x$distances[i]
          cat(paste("Salesperson", i, ": distance travelled =", round(dist, 2)), sep="\n")
          cat(paste("Tour:", print_tours), sep="\n")
          cat("\n")
        }
     }
    }



#' Arranges the \code{mtsp} solution for plotting.
#'
#' @param x The list returned by \code{mtsp}.
#' @param dat A matrix of x,y positions passed to \code{mtsp}. Only needed when \code{return_all=FALSE} was used.
#' @param ... Other arguments.
#' @export 
fortify <- function(x, dat, ...) UseMethod("fortify")

#' @export 
fortify.mtsp_ga <- function(x, dat=NULL, ...){

                if(is.null(x$xy) & is.null(dat)) stop("The position matrix must be provided. Re-run mtsp with the option return_all=TRUE
                                                        or else provide the data.")

                plot_dat = do.call(rbind,
                                  lapply(seq_along(x$best_tour), function(tour) {
                                                                        route = x$best_tour[[tour]]
                                                                        from = x$xy[route,]
                                                                        colnames(from) = paste0("from_", c("x", "y"))
                                                                        to = x$xy[c(route[-1], route[1]),]
                                                                        colnames(to) = paste0("to_", c("x", "y"))
                                                                        pdat = data.frame(from, to, id=tour)
                                                                        pdat
                                                                        }
                                         ))
                return(plot_dat)
                }



#' Plot the MTSP solution
#' 
#' Plot the MTSP solution (currently only implemented for the non-depot algorithms).
#'
#' @param x The list returned by \code{mtsp}.
#' @param dat A matrix of x,y positions passed to \code{mtsp}. Only needed when \code{return_all=FALSE} was used.
#' @param ... Unused.
#' @examples
#' set.seed(1)
#' n = 25
#' xy = matrix(rnorm(n*2), ncol=2)
#' run = mtsp(xy, nSalesmen=5, CostType=2, popSize=80, 
#'            numIter=10, algorithm="mtsp_ga", return_all=TRUE)
#' plot(run) 
#' @importFrom graphics box segments
#' @export
plot.mtsp_ga  <- function(x, dat=NULL, ...){

  plot_dat = fortify(x, dat)

  plot(plot_dat[c("from_x","from_y")], axes=FALSE, xlab="", ylab="", pch=16, col=plot_dat$id+1)
  box()
  segments(plot_dat[["from_x"]], plot_dat[["from_y"]], plot_dat[["to_x"]], plot_dat[["to_y"]], col=plot_dat$id+1)

  } 


#' Plot
#'
#' @rdname plot.mtsp_ga 
#' @export 
ggplot <- function(x, dat, ...) UseMethod("ggplot")

#' @export 
ggplot.mtsp_ga  <- function(x, dat=NULL, ...){
  
    if (!requireNamespace("ggplot2", quietly = TRUE)) 
      stop("Package \"ggplot2\" needed for this function to work.")

  plot_dat = fortify(x, dat)

  ggplot2::ggplot(plot_dat, ggplot2::aes(x=from_x, y=from_y, xend=to_x, yend=to_y, col=factor(id))) +
    ggplot2::geom_point() + 
    ggplot2::geom_segment() +
    ggplot2::scale_colour_discrete(name = "Salesman") +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "top", panel.border = ggplot2::element_rect(fill=NA))
}




