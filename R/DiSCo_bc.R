#' Function for computing barycenters in the DiSCo method at every time period
#'
#' Compute barycenters in the DiSCo method at every time period, as in Definition 1,
#' Step 4 in \insertCite{gunsilius2023distributional}.
#'
#' @param controls.q List with matrices of control quantile functions
#' @param weights Vector of optimal synthetic control weights, computed using the DiSCo_weights_reg function.
#' @param evgrid Optional vector containing an evenly spaced grid on [0,1] on which the quantile function for the control units will be evaluated
#' By default, a grid of 100 points is used.
#' @return The quantile function of the barycenter associated with the "weights" evaluated at the vector "evgrid"
#' @references
#' \insertAllCited{}
DiSCo_bc <- function(controls.q, weights, evgrid = seq(from=0, to=1, length.out=101)){

  # Obtaining the Wasserstein barycenter as the average of the quantile functions
  # weighted by "weights" and evaluate it on the grid "evgrid"
  thebc <- controls.q%*%weights

  # returning the barycenter
  return(thebc)
}
