#' Confidence Interval in the DSC Method
#'
#' Function for computing the confidence interval (at a specific period) in the DSC method
#' @param df_control a list or matrix:
#' in list form, each list element contains a vector of observations for the given control unit;
#' in matrix form, each column corresponds to one unit and each row is one observation.
#' @param df_target a vector of observations.
#' @param M an integer specifying the number of draws from the uniform distribution for approximating the integral.
#' @param solver a solver for the optimization problem; see \code{\link[CVXR]{installed_solvers}} in CVXR for more options.
#' @param w_t 0 or a vector. By default, i.e. 0, \code{DSC_CI} calculates the confidence interval of a pre-treatment period,
#' weights of control units are not needed to be specified. To calculate the confidence interval of a post-treatment period,
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
#' #detail can be found by ??ex_gmm
#' DSC_CI(df_control=ex_gmm()$control, df_target=ex_gmm()$target, M=100, solver="SCS", w_t=0, cl=0.99, num.redraws=500, evgrid = seq(from=0, to=1, length.out=1001),graph=TRUE, y_name='y', x_name='x')





DSC_CI=function(df_control, df_target, M=100, solver="SCS", w_t=0, cl=0.99, num.redraws=500, evgrid = seq(from=0, to=1, length.out=1001),
                graph=TRUE, y_name='y', x_name='x'){

  ## Creating a function for the empirical quantile function
  myquant <- function(X,q){
    # sort if unsorted
    if (is.unsorted(X)) X <- sort(X)
    # compute empirical CDF
    X.cdf <- 1:length(X) / length(X)
    # obtain the corresponding empirical quantile
    return(X[which(X.cdf >= q)[1]])
  }

  if (is.list(df_control)==FALSE) {
    controls.h <- list()
    for (ii in 1:ncol(df_control)){
      controls.h[[ii]] <- as.vector(df_control[,ii])
    }
    df_control <- controls.h
  }


  if (length(w_t)==1){
    #if pre-treatment, split the data to obtain weights

    target <- list(3)
    target[[3]] <- df_target
    v <- as.vector(c(rep(TRUE, floor(0.5*length(target[[3]]))),
                     rep(FALSE,ceiling(0.5*length(target[[3]])))))
    set.seed(1860)
    ind <- sample(v)
    target[[1]] <- target[[3]][ind]
    target[[2]] <- target[[3]][!ind]


    # obtaining all control states
    controls=list()
    for (ii in 1:length(df_control)){
      controls[[ii]] <- list(3)
      controls[[ii]][[3]] <- df_control[[ii]]
      v <- as.vector(c(rep(TRUE, floor(0.5*length(controls[[ii]][[3]]))),
                       rep(FALSE,ceiling(0.5*length(controls[[ii]][[3]])))))
      set.seed(1860)
      ind <- sample(v)
      controls[[ii]][[1]] <- controls[[ii]][[3]][ind]
      controls[[ii]][[2]] <- controls[[ii]][[3]][!ind]
    }


    # generating the target
    target.s <- mapply(myquant, evgrid, MoreArgs =
                         list(X=target[[1]]))

    # generating the controls
    controls1 <- list(length(controls))
    for (ii in 1:length(controls)){
      controls1[[ii]] <- controls[[ii]][[1]]
    }

    # obtaining the optimal weights
    DSC_res <- DSC_weights(controls1,target[[1]], M, choose_solver = solver)
    weights <- DSC_res[[1]]


    # use the second sample to obtain confidence intervals

    # generating the target
    target.s <- mapply(myquant, evgrid, MoreArgs = list(X=target[[2]]))


    controls2 <- list(length(controls))
    for (ii in 1:length(controls)){
      controls2[[ii]] <- controls[[ii]][[2]]
    }
  }else{
    #if post-treatment, use optimal weights
    weights=w_t
    controls2=df_control
  }



  # obtaining the barycenter of the control units
  DSC_res2 <- DSC_bc(controls2, weights, evgrid)


  # obtaining the barycenter of the control units
  DSC_res2.CI <- matrix(0, nrow = num.redraws, ncol = length(evgrid))
  DSC_res2.CI[1,] <- DSC_bc(controls2, weights, evgrid)[[2]]



  # Starting the loop for all redraws in the bootstrap approach

  pb <- txtProgressBar(min=0, max=100, style=3) #initiate progress bar

  for (redraws in 2:num.redraws){

    Sys.sleep(0.025)
    setTxtProgressBar(pb, round(redraws / num.redraws * 100)) # progress bar
    cat(if (redraws == num.redraws) '\n')

    set.seed(redraws*1) # for reproducibility
    # drawing m = 100% of samples from controls2

    mycon <- list(length(controls2))
    for (ii in 1:length(controls2)){
      sz <- length(controls2[[ii]])
      m.c <- floor(1*sz)
      # sampling from controls
      mycon[[ii]] <- sample(controls2[[ii]],m.c, replace=TRUE) # the sample
      # for the controls
    }

    # for the given weights lambda^* compute the barycenter
    DSC_res2.CI[redraws,] <- DSC_bc(mycon,weights,evgrid)[[2]]
  }

  # obtain the cl% confidence interval
  CI.u <- apply(DSC_res2.CI,2,quantile, probs=cl+(1-cl)/2)
  CI.l <- apply(DSC_res2.CI,2,quantile, probs=(1-cl)/2)

  if (graph==TRUE){

    # generating the target for the plot
    target.s <- mapply(myquant, evgrid, MoreArgs =
                         list(X=target[[2]]))

    length(which(weights>10^(-4))) # all "essentially non-zero" elements

    which(weights>10^(-4))

    # plotting
    plot(evgrid, target.s, xlab='',ylab='',
         type='l',lwd=1)
    lines(evgrid, DSC_res2[[2]],lwd=2, col='red')
    lines(evgrid, CI.u, lwd = 2, col='red', lty=2)
    lines(evgrid, CI.l, lwd = 2, col='red', lty=2)
    legend("topleft",legend = c("Target", "DSC", "CI"),
           col=c("black", "red", "red"),
           lty= c(1,1,2), lwd = c(2,2,2), cex = 1.5)
    title(ylab=y_name, line=2, cex.lab=1.5)
    title(xlab=x_name, line=2, cex.lab=1.5)
  }

  return(list(CI.u=CI.u, CI.l=CI.l))

}
