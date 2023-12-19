
#' @title DiSCo_mixture
#' @description The alternative mixture of distributions approach in the paper
#' @param controls1 A list of controls
#' @param target The target unit
#' @param grid.min Minimal value of the grid on which the CDFs are evaluated.
#' @param grid.max Maximal value of the grid on which the CDFs are evaluated.
#' @param grid.rand Random grid on which the CDFs are evaluated.
#' @inheritParams DiSCo
#' @return A list containing the following elements:
#' \itemize{
#' \item \code{cdf } A matrix containing the CDFs of the target and control units evaluated on the grid.
#' \item \code{distance.opt } The optimal value of the Wasserstein distance.
#' \item \code{mean } The optimal value of the Wasserstein barycenter.
#' \item \code{target.order } The target unit, ordered.
#' \item \code{weights.opt } The optimal weights.
#' }
#' @keywords internal
DiSCo_mixture <- function(controls1, target, grid.min, grid.max, grid.rand, M) {

  ###### The mixture of distributions approach
  # we again only focus on the first half of the data
  # defining the grid on which we define the cumulative distribution functions
  # obtaining the minimal and maximal values among all supports
  # creating a list of controls with only the full data


  # Estimating the empirical CDFs
  CDF.control <- lapply(controls1,stats::ecdf)
  CDF.target <- stats::ecdf(target)


  # Evaluating the CDF on the random grid
  CDF.matrix <- matrix(0,nrow=length(grid.rand), ncol = (length(controls1)+1))
  CDF.matrix[,1] <- CDF.target(grid.rand)
  for (ii in 1:length(controls1)){
    CDF.matrix[,(ii+1)] <- CDF.control[[ii]](grid.rand)
  }

  # Solving the convex problem with CVXR
  # the variable we are after
  theweights <- CVXR::Variable(length(controls1))
  # the objective function
  objective <- CVXR::cvxr_norm((CDF.matrix[,2:ncol(CDF.matrix)] %*% theweights - CDF.matrix[,1]))

  # the constraints for the unit simplex
  constraints <- list(theweights>=0, CVXR::sum_entries(theweights) == 1)
  # the optimization problem
  problem <- CVXR::Problem(CVXR::Minimize(objective),constraints)
  # solving the optimization problem
  results <- CVXR::solve(problem, solver = "SCS")

  # returning the optimal weights and the value function which provides the
  # squared Wasserstein distance between the target and the corresponding barycenter
  theweights.opt <- results$getValue(theweights)
  thedistance.opt <- results$value*1/M*(grid.max - grid.min)


  themean <- CDF.matrix[,2:ncol(CDF.matrix)]%*%theweights.opt


  themean.order <- themean[order(grid.rand, decreasing=FALSE)]
  target.order <- CDF.matrix[order(grid.rand, decreasing=FALSE),1]


  return(list("weights.opt" = theweights.opt, "distance.opt" = thedistance.opt,
              "mean" = themean.order, "target.order" = target.order, "cdf" = CDF.matrix))

}
