
#' @title DiSCo_CI_iter
#'
#' @description Function for computing the confidence intervals in the DiSCo method in a single period
#' @param t Time period
#' @param controls_t List of control unit data for given period
#' @param target_t List of target unit data for given period
#' @inheritParams DiSCo_CI
#' @inheritParams DiSCo
#' @return The resampled counterfactual barycenter of the target unit
#' @keywords internal
DiSCo_CI_iter <- function(t, controls_t, target_t, grid, T0, M=1000,
                          evgrid = seq(from=0, to=1, length.out=1001), qmethod=NULL,
                          qtype=7,
                          mixture=FALSE, simplex=FALSE, replace=TRUE) {

  # resample target
  t_len <- length(target_t)
  mytar <- target_t[sample(1:t_len, floor(1*t_len), replace=replace)]
  mytar.q <- myQuant(mytar, evgrid, qmethod, qtype=qtype) # quantile
  mytar.cdf <- stats::ecdf(mytar)(grid)

  # resample controls
  mycon_list <- list()
  mycon.q <- matrix(0,nrow = length(evgrid), ncol=length(controls_t))
  mycon.cdf <- matrix(0,nrow = length(grid), ncol=length(controls_t)+1)
  mycon.cdf[,1] <- mytar.cdf


  for (ii in 1:length(controls_t)){
    controls_t_i <- controls_t[[ii]]
    c_len <- length(controls_t_i)
    mycon <- controls_t_i[sample(1:c_len, floor(1*c_len), replace=replace)] # resample
    mycon_list[[ii]] <- mycon
    mycon.q[,ii] <- myQuant(mycon, evgrid, qmethod, qtype=qtype) # resampled quantile
    mycon.cdf[,ii+1] <- stats::ecdf(mycon)(grid) # resampled cdf
  }

  if (t <= T0) { # if pre-treatment, calculate bootstrapped weights
    if (!mixture) {
      lambda <- DiSCo_weights_reg(mycon_list, mytar, M=M, qmethod=qmethod, simplex=simplex)
    } else {
      mixt <- DiSCo_mixture_solve(length(controls_t), mycon.cdf, min(grid), max(grid),
                                  grid, M, simplex)
      lambda <- mixt$weights.opt
    }
  } else { # if post-treatment, no weights are calculated, but we still want to resample the data for bootstrap
    lambda <- NULL
  }

  return(list("weights" = lambda, "target" = list("quantile" = mytar.q, "cdf" = mytar.cdf),
         "controls" = list("quantile" = mycon.q, "cdf" = mycon.cdf[,-1])))

}

#' @title bootCounterfactuals
#' @description Function for computing the bootstrapped counterfactuals in the DiSCo method
#' @param result_t A list containing the results of the DiSCo_CI_iter function
#' @param t The current time period
#' @inheritParams DiSCo_CI
#' @return A list containing the bootstrapped counterfactuals
#' @keywords internal
bootCounterfactuals <- function(result_t, t, mixture, weights, evgrid, grid) {
  if (mixture) {
    # calculate cdf and then back out quantile
    cdf_t <- result_t$controls$cdf %*%  weights # cdf
    q_t <- sapply(evgrid, function(y) grid[which(cdf_t >= (y-(1e-5)))[1]]) # tolerance accounts for inaccuracy (esp != 1)

  } else {
    # calculate quantile then back out cdf
    q_t <- DiSCo_bc(result_t$controls$quantile, weights, evgrid)
    cdf_t <- stats::ecdf(q_t)(grid)

  }
  # calculate bootstrapped counterfactual differences
  cdf_diff <- as.vector(result_t$target$cdf - cdf_t)
  q_diff <- as.vector(result_t$target$quantile - q_t)

  return(list("cdf" = cdf_t, "quantile" = q_t, "quantile_diff" = q_diff, "cdf_diff" = cdf_diff))
}



#' @title DiSCo_CI
#'
#' @description Function for computing the confidence intervals in the DiSCo method
#' using the bootstrap approach described in
#' @param redraw Integer indicating the current bootstrap redraw
#' @param controls A list containing the raw data for the control group
#' @param target A list containing the raw data for the target group
#' @param T_max Index of last time period
#' @param T0 Index of the last pre-treatment period
#' @param mc.cores Number of cores to use for parallelization
#' @param grid Grid to recompute the CDF on if `mixture` option is chosen
#' @inheritParams DiSCo
#' @return A list with the following components
#' \itemize{
#' \item \code{weights} The bootstrapped weights
#' \item \code{disco_boot} A list containing the bootstrapped counterfactuals,
#' with the following elements, each of which contains named elements called `upper` and `lower`
#' which are G x T matrices where G is the specified number of grid points and T is the number of time periods
#' }
#' @keywords internal
DiSCo_CI <- function(redraw, controls, target, T_max, T0, grid, mc.cores=1,
                     evgrid = seq(from=0, to=1, length.out=1001), qmethod=NULL, qtype=7,
                     M=1000,mixture=FALSE, simplex=FALSE, replace=TRUE) {


  boots.periods <- lapply(1:T_max, function(t) DiSCo_CI_iter(t, controls_t=controls[[t]],
                                            target_t=target[[t]], grid=grid[[t]], T0=T0, M=M,
                                             evgrid = evgrid, qmethod=qmethod, qtype=qtype,
                                            mixture=mixture, simplex=simplex, replace=replace)
         )

  # extract the weights
  weights <- 0
  for (t in 1:T0) {
    weights <- weights + boots.periods[[t]]$weights
  }
  weights <- (1/T0)*weights

  disco_boot <- list()

  # compute resampled counterfactuals
  disco_boot <- lapply(1:T_max, function(t)
    bootCounterfactuals(boots.periods[[t]], t, mixture, weights, evgrid, grid[[t]])
  )

  return(list("weights"=weights, "disco_boot"=disco_boot))

}
