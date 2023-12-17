#' Function for computing barycenters in the DiSCo method at every time period
#'
#' Compute barycenters in the DiSCo method at every time period, as in Definition 1,
#' Step 4 in \insertCite{gunsilius2023distributional;textual}{DiSCos}.
#'
#' @param controls.q List with matrices of control quantile functions
#' @param weights Vector of optimal synthetic control weights, computed using the DiSCo_weights_reg function.
#' @inheritParams DiSCo
#' @return The quantile function of the barycenter associated with the "weights" evaluated at the vector "evgrid"
#' @references
#' \insertAllCited{}
#' @keywords internal
DiSCo_bc <- function(controls.q, weights, evgrid = seq(from=0, to=1, length.out=101)){

  # Obtaining the Wasserstein barycenter as the average of the quantile functions
  # weighted by "weights" and evaluate it on the grid "evgrid"
  thebc <- controls.q%*%weights

  # returning the barycenter
  return(thebc)
}
