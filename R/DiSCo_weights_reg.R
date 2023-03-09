#' Function for obtaining the weights in the DiSCo method at every time period

#' Estimate the optimal weights for the distributional synthetic controls method by
#' solving the convex minimization problem in Eq. (2) in \insertref{gunsilius2020distributional}{DiSCo}
#' using a regression of the simulated target quantile on the simulated control quantiles, as in Eq. (3),
#' \eqn{\underset{\vec{\lambda} \in \Delta^J}{\operatorname{argmin}}\left\|\mathbb{Y}_t \vec{\lambda}_t-\vec{Y}_{1 t}\right\|_2^2}.
#' For the constrained optimization we rely on the package pracma
#' the control distributions can be given in list form, where each list element contains a
#' vector of observations for the given control unit, in matrix form;
#' in matrix- each column corresponds to one unit and each row is one observation.
#' The list-form is useful, because the number of draws for each control group can be different.
#' The target must be given as a vector.
#'
#' @param controls List with matrices of control distributions
#' @param target Matrix containing the target distribution
#' @param M Optional integer, number of draws from the uniform distribution for approximating the integral. See section 3.1 in the paper.
#'
#' @return
DiSCo_weights_reg <- function(controls,target, M = 500){

  if (!is.vector(target)){
    stop("Target needs to be given as a vector.")
  }
  if (!is.list(controls) && !is.matrix(controls)){
    stop ("Controls need to be given in either list-form or matrix form.")
  }
  # if the controls are given in matrix form we turn them into a list
  if (is.matrix(controls)) {
    controls.h <- list()
    for (ii in 1:ncol(controls)){
      controls.h[[ii]] <- as.vector(controls[,ii])
    }
    controls <- controls.h
  }
  # M is the number of draws from the uniform distribution for approximating the integral

  ## Creating a function for the empirical quantile function
  myquant <- function(X,q){
    # sort if unsorted
    if (is.unsorted(X)) X <- sort(X)
    # compute empirical CDF
    X.cdf <- 1:length(X) / length(X)
    # obtain the corresponding empirical quantile
    return(X[which(X.cdf >= q)[1]])
  }

  ## Sampling from this quantile function M times
  Mvec <- runif(M, min = 0, max = 1)
  controls.s <- matrix(0,nrow = M, ncol = length(controls))
  for (jj in 1:length(controls)){
    controls.s[,jj] <- mapply(myquant, Mvec, MoreArgs = list(X=controls[[jj]]))
  }

  target.s <- matrix(0, nrow = M, ncol=1)
  target.s[,1] <- mapply(myquant, Mvec, MoreArgs = list(X=target))

  ## Solving the optimization using constrained linear regression

    # the equality constraints
    Aequ <- matrix(rep(1,length(controls)),nrow = 1, ncol = length(controls))
    # if the values in controls.s and target.s are too large it can happen that we run into
    # overflow errors. For this reason we scale both the vector and the matrix by the Frobenius
    # norm of the matrix
    sc <- norm(controls.s,"2")

  weights.opt <- pracma::lsqlincon(controls.s/sc,target.s/sc, A=NULL,b=NULL,Aeq = Aequ, beq = 1, 0, 1)


  return(weights.opt)
}
