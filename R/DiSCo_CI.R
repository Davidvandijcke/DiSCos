
#' @title DiSCo_CI_iter
#'
#' @description Function for computing the confidence intervals in the DiSCo method in a single period
#' @param redraw Integer indicating the number of times the function has been called
#' @param controls List of control units
#' @param weights Vector of optimal synthetic control weights, computed using the DiSCo_weights_reg function.
#' @inheritParams DiSCo_CI
#' @inheritParams DiSCo
#' @return The resampled counterfactual barycenter of the target unit
#' @keywords internal
DiSCo_CI_iter <- function(redraw, controls, weights, grid, cl=0.95,
                          evgrid = seq(from=0, to=1, length.out=1001), qmethod=NULL,
                          mixture=FALSE) {
  # set.seed(redraw*1) # for reproducibility
  # drawing m = 100% of samples from controls

  mycon <- list()
  mycon.q <- matrix(0,nrow = length(evgrid), ncol=length(controls))
  for (ii in 1:length(controls)){
    sz <- length(controls[[ii]])
    m.c <- floor(1*sz)
    # sampling from controls
    idx <- sample(1:sz, m.c, replace=TRUE)
    mycon[[ii]] <- controls[[ii]][idx]
    # mycon.q[,ii] <- mapply(myquant, evgrid, MoreArgs = list(X=mycon[[ii]]))
    mycon.q[,ii] <- myQuant(mycon[[ii]], evgrid, qmethod)
  }


  if (!mixture) {
    out <- DiSCo_bc(mycon.q,weights,evgrid) # TODO: left off here
  } else {
    mycon.cdf <- apply(mycon.q, 2, function(x) stats::ecdf(x)(grid))
    cdf <- mycon.cdf %*% weights

    # this is slightly hacky cause we will have to recompute the cdf later but it allows me to keep the old DiSCoTEA function
    # the tolerance serves to account for the imperfect optimization of CVXR (esp fact that CDF can be <1 everywhere bcs of it)
    out <-  sapply(evgrid, function(y)
      grid[which(cdf >= (y-(1e-5)))[1]])
  }
  return(out)
}



#' @title DiSCo_CI
#'
#' @description Function for computing the confidence intervals in the DiSCo method
#' @param controls A list containing the raw data for the control group
#' @param weights A vector of optimal weights
#' @param mc.cores Number of cores to use for parallelization
#' @param num.redraws The number of bootstrap samples to draw
#' @param grid Grid to recompute the CDF on if `mixture` option is chosen
#' @inheritParams DiSCo
#' @return \code{DSC_CI} returns a list containing the following components:
#' \itemize{
#' \item \code{upper } A vector of the upper bound.
#' \item \code{lower } A vector of the lower bound.
#' \item \code{se } A vector of the standard errors of each counterfactual quantile estimate.
#' \item \code{bootmat } A matrix of the counterfactual quantile estimates for each bootstrap sample.
#' }
#' @keywords internal
DiSCo_CI <- function(controls, weights, grid, mc.cores=1, cl=0.95, num.redraws=500,
                     evgrid = seq(from=0, to=1, length.out=1001), qmethod=NULL,
                     mixture=FALSE){

  DSC_res2.CI <- mclapply.hack(1:num.redraws, DiSCo_CI_iter, controls=controls,
                               weights=weights, grid=grid, cl=cl, evgrid=evgrid, mc.cores=mc.cores,
                               mixture=mixture)

  DSC_res2.CI <- sapply(DSC_res2.CI, function(x) x)

  # DiSCo_CI_iter(1, controls=controls, bc=bc, weights=weights, cl=cl, evgrid = evgrid)
  # obtain the cl% confidence intervals
  if (!is.null(qmethod)){
    if (qmethod=="qkden") {
      # estimate bandwidth once
      bw <- stats::bw.nrd0(DSC_res2.CI[,1])
    }
  }
  # we do the quantiles in parallel cause the qmethod=qkden is slow
  CI.u <- unlist(mclapply.hack(1:dim(DSC_res2.CI)[1], function(x) myQuant(DSC_res2.CI[x,], q=cl+(1-cl)/2, qmethod=qmethod, bw=bw), mc.cores=mc.cores))
  CI.l <- unlist(mclapply.hack(1:dim(DSC_res2.CI)[1], function(x) myQuant(DSC_res2.CI[x,], q=(1-cl)/2, qmethod=qmethod, bw=bw), mc.cores=mc.cores))

  se <- apply(DSC_res2.CI,1,sd)


  return(list(upper=CI.u, lower=CI.l, se=se, bootmat=DSC_res2.CI))

}
