
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




#' @title parseBoots
#' @description Function for parsing the bootstrapped counterfactuals in the DiSCo method
#' @param CI_temp A list containing the bootstrapped counterfactuals
#' @param cl The confidence level
#' @param q_disco The estimated quantiles around which to center
#' @param cdf_disco The estimated cdfs around which to center
#' @param q_obs The observed quantiles
#' @param cdf_obs The observed cdfs
#' @param uniform Whether to use uniform or pointwise confidence intervals
#' @return A list containing the confidence intervals for the quantiles and cdfs
#' @keywords internal
parseBoots <- function(CI_temp, cl, q_disco, cdf_disco, q_obs, cdf_obs, uniform=TRUE) {
  ## CIs for quantiles and cdfs

  # Extract and combine all data into a single call, and immediately convert to an array
  extract_and_combine <- function(data, attribute) {
    combined <- simplify2array(
      lapply(data, function(x) sapply(x$disco_boot, function(y) y[[attribute]]))
    )
  }
  # unfold list to matrix
  unfold <- function(x) x

  q_d <- sapply(q_disco, unfold)
  cdf_d <- sapply(cdf_disco, unfold)
  q_diff_d <- sapply(q_obs, unfold) - q_d
  cdf_diff_d <- sapply(cdf_obs, unfold) - cdf_d

  q_boot <- extract_and_combine(CI_temp, "quantile")
  cdf_boot <- extract_and_combine(CI_temp, "cdf")
  q_diff <- extract_and_combine(CI_temp, "quantile_diff")
  cdf_diff <- extract_and_combine(CI_temp, "cdf_diff")

  # calculate confidence intervals
  getCIs <- function(btmat, cl, og, uniform) {
    bt_diff <- sweep(btmat, c(1,2), og, "-")

    if (uniform) {
      bt_diff <- apply(abs(bt_diff), c(2,3), max) # outputs a G X T_max matrix
      t_all <- apply(bt_diff, c(1), function(x) stats::quantile(x, probs=cl))
      t_all <- t(replicate(nrow(og), t_all))
      upper <- og + t_all
      lower <- og - t_all
    } else {
      t_upper <- apply(bt_diff, c(1,2), function(x) stats::quantile(x, probs=(1-cl)/2)) # outputs a G X T_max matrix
      t_lower <- apply(bt_diff, c(1,2), function(x) stats::quantile(x, probs= cl+(1-cl)/2)) # outputs a G X T_max matrix
      upper <- og - t_upper
      lower <- og - t_lower
    }

    se <- apply(btmat, c(1,2), function(x) stats::sd(x)) # outputs a G X T_max matrix
    return(list("lower"=lower, "upper"=upper, "se" = se))
  }
  q_CI <- getCIs(q_boot, cl, q_d, uniform)
  cdf_CI <- getCIs(cdf_boot, cl, cdf_d, uniform)
  q_diff_CI <- getCIs(q_diff, cl, q_diff_d, uniform)
  cdf_diff_CI <- getCIs(cdf_diff, cl, cdf_diff_d, uniform)

  ## CIs for weights
  weights <- sapply(CI_temp, function(x) x$weights)
  weights_CI <- list(
    "upper" = apply(weights, c(1), function(x) stats::quantile(x, probs=cl+(1-cl)/2)),
    "lower"=  apply(weights, c(1), function(x) stats::quantile(x, probs=(1-cl)/2))
  )

  # create CI object
  CI_out <- list("quantile" = q_CI, "cdf" = cdf_CI, "quantile_diff" =
                   q_diff_CI, "cdf_diff" = cdf_diff_CI, "weights" = weights_CI,
                 bootmat = list("quantile" = q_boot, "cdf" = cdf_boot, "quantile_diff" = q_diff, "cdf_diff" = cdf_diff))
  return(CI_out)
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
