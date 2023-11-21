controls=controls
bc=bc
weights=Weights_DiSCo_avg
mc.cores=num.cores
cl=cl
num.redraws=boots
evgrid = evgrid


DiSCo_CI_iter <- function(redraw, controls, bc, weights, cl=0.95, evgrid = seq(from=0, to=1, length.out=1001)){
  set.seed(redraw*1) # for reproducibility
  # drawing m = 100% of samples from controls

  mycon <- list()
  mycon.q <- matrix(0,nrow = length(evgrid), ncol=length(controls))
  for (ii in 1:length(controls)){
    sz <- length(controls[[ii]])
    m.c <- floor(1*sz)
    # sampling from controls
    idx <- sample(1:sz, m.c, replace=TRUE)
    mycon[[ii]] <- controls[[ii]][idx]
    mycon.q[,ii] <- mapply(myquant, evgrid, MoreArgs = list(X=mycon[[ii]]))
  }


  return(DiSCo_bc(mycon,mycon.q,weights,evgrid))
}




DiSCo_CI <- function(controls, bc, weights, mc.cores=1, cl=0.95, num.redraws=500, evgrid = seq(from=0, to=1, length.out=1001)){
  #' Confidence Intervals for the DiSCo quantile
  #'
  #' Function for computing the confidence interval (at a specific period) in the DSC method
  #' @param df_control a list or matrix:
  #' in list form, each list element contains a vector of observations for the given control unit;
  #' in matrix form, each column corresponds to one unit and each row is one observation.
  #' @param df_target a vector of observations.
  #' @param M an integer specifying the number of draws from the uniform distribution for approximating the integral.
  #' @param solver a solver for the optimization problem; see \code{\link[CVXR]{installed_solvers}} in CVXR for more options.
  #' @param w_t 0 or a vector. By default, i.e. 0, \code{DSC_CI} calculates the confidence interval of a pre-treatment period,
  #' weights of control units do not need to be specified. To calculate the confidence interval of a post-treatment period,
  #' a vector of optimal weights (pre-calculated using pre-treatment observations) is needed.
  #' @param cl a float specifying the confidence level.
  #' @param num.redraws an integer specifying the number of redraws used in the bootstrap approach.
  #' @param evgrid a vector of gridpoints used to evaluate quantile functions.
  #' @param graph \code{TRUE/FALSE}, indicating whether to output a plot of the confidence interval.
  #' @param y_name a string for the title of the y-axis.
  #' @param x_name a string for the title of the x-axis.
  #' @return \code{DSC_CI} returns a list containing the following components:
  #' \item{\code{CI.u}}{a vector of the upper bound.}
  #' \item{\code{CI.l}}{a vector of the lower bound.}
  #' @export
  #' @examples
  #' #simulated data from Gaussian mixture
  #' #ex_gmm() calls the simulated data
  #' #details can be found by ??ex_gmm
  #' DSC_CI(df_control=ex_gmm()$control, df_target=ex_gmm()$target, M=100, solver="SCS", w_t=0, cl=0.99, num.redraws=500, evgrid = seq(from=0, to=1, length.out=1001),graph=TRUE, y_name='y', x_name='x')


  DSC_res2.CI <- mclapply.hack(1:num.redraws, DiSCo_CI_iter, controls=controls, bc=bc,
                               weights=weights, cl=cl, evgrid=evgrid, mc.cores=mc.cores)

  DSC_res2.CI <- sapply(DSC_res2.CI, function(x) x$barycenter)

  # DiSCo_CI_iter(1, controls=controls, bc=bc, weights=weights, cl=cl, evgrid = evgrid)
  # obtain the cl% confidence interval
  CI.u <- apply(DSC_res2.CI,1,stats::quantile, probs=cl+(1-cl)/2)
  CI.l <- apply(DSC_res2.CI,1,stats::quantile, probs=(1-cl)/2)

  se <- apply(DSC_res2.CI,1,sd)


  return(list(upper=CI.u, lower=CI.l))

}
