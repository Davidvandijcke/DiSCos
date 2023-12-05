

DiSCo_CI_iter <- function(redraw, controls, bc, weights, cl=0.95,
                          evgrid = seq(from=0, to=1, length.out=1001), qmethod=NULL){
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


  return(DiSCo_bc(mycon,mycon.q,weights,evgrid))
}




DiSCo_CI <- function(controls, bc, weights, mc.cores=1, cl=0.95, num.redraws=500, evgrid = seq(from=0, to=1, length.out=1001), qmethod=NULL){
  #' Confidence Intervals for the DiSCo quantile
  #'
  #' Function for computing the confidence intervals in the DiSCo method
  #' @param controls A list containing the raw data for the control group
  #' @param bc A list of barycenters for the control group
  #' @param weights A vector of optimal weights
  #' @param mc.cores Number of cores to use for parallelization
  #' @param cl The confidence level for the confidence interval
  #' @param num.redraws The number of bootstrap samples to draw
  #' @param evgrid The grid of quantiles to evaluate the counterfactual quantile function
  #' @return \code{DSC_CI} returns a list containing the following components:
  #' \itemize{
  #' \item{\code{upper} }{a vector of the upper bound.}
  #' \item{\code{lower} }{a vector of the lower bound.}
  #' \item{\code{se} }{A vector of the standard errors of each counterfactual quantile estimate.}
  #' \item{\code{bootmat} }{A matrix of the counterfactual quantile estimates for each bootstrap sample.}
  #' }


  DSC_res2.CI <- mclapply.hack(1:num.redraws, DiSCo_CI_iter, controls=controls, bc=bc,
                               weights=weights, cl=cl, evgrid=evgrid, mc.cores=mc.cores)

  DSC_res2.CI <- sapply(DSC_res2.CI, function(x) x)

  # DiSCo_CI_iter(1, controls=controls, bc=bc, weights=weights, cl=cl, evgrid = evgrid)
  # obtain the cl% confidence interval
  if (!is.null(qmethod)){
    if (qmethod=="qkden") {
      # estimate bandwidth once
      bw <- stats::bw.nrd0(DSC_res2.CI[,1])
    }
  }
  # we do the quantiles in parallel cause the qmethod=qkden is slow
  CI.u <- unlist(mclapply.hack(1:dim(DSC_res2.CI)[2], function(x) myQuant(DSC_res2.CI[,x], q=cl+(1-cl)/2, qmethod=qmethod, bw=bw), mc.cores=mc.cores))
  CI.l <- unlist(mclapply.hack(1:dim(DSC_res2.CI)[2], function(x) myQuant(DSC_res2.CI[,x], q=(1-cl)/2, qmethod=qmethod, bw=bw), mc.cores=mc.cores))

  se <- apply(DSC_res2.CI,1,sd)


  return(list(upper=CI.u, lower=CI.l, se=se, bootmat=DSC_res2.CI))

}
