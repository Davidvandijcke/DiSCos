#' @title DiSCo_weights_reg
#' @description Function for obtaining the weights in the DiSCo method at every time period
#' @details Estimate the optimal weights for the distributional synthetic controls method.
#' solving the convex minimization problem in Eq. (2) in \insertCite{gunsilius2023distributional;textual}{DiSCos}..
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
#' @inheritParams DiSCo
#'
#' @return Vector of optimal synthetic control weights
#' @references
#' \insertAllCited{}
#' @keywords internal
DiSCo_weights_reg <- function(controls, target, M = 500, qmethod=NULL, qtype=7, simplex=FALSE, q_min=0, q_max=1){


  if (!is.null(qmethod)){
    if (qmethod=="extreme"){
      qmethod <- NULL # heavy tails of extreme quantiles mess up the regression, so just use the standard one
    }
  }

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

  ## Sampling from this quantile function M times
  Mvec <- stats::runif(M, min = q_min, max = q_max)
  controls.s <- matrix(0,nrow = M, ncol = length(controls))
  for (jj in 1:length(controls)){
    controls.s[,jj] <- myQuant(controls[[jj]], Mvec, qmethod, qtype=qtype)
  }

  target.s <- matrix(0, nrow = M, ncol=1)
  target.s[,1] <- myQuant(target, Mvec, qmethod, qtype=qtype)

  ## Solving the optimization using constrained linear regression

    # the equality constraints
    Aequ <- matrix(rep(1,length(controls)),nrow=1, ncol=length(controls))
    # if the values in controls.s and target.s are too large it can happen that we run into
    # overflow errors. For this reason we scale both the vector and the matrix by the Frobenius
    # norm of the matrix
    sc <- CVXR::norm(controls.s,"2")

  if (simplex) { lb <- 0 } else { lb <- NULL }
  # mirror the behavior of the lsqlincon, which explicitly requires
  # quadprog to be loaded which confuses R CMD check
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("quadprog needed for this function to work. Please install it.",
         call. = FALSE)
  }
  # solve with the pracma package, which runs FORTRAN under the hood so it is fast
  weights.opt <- pracma::lsqlincon(controls.s/sc,target.s/sc, A=NULL,b=NULL,
                  Aeq = Aequ, # LHS of equality constraint: matrix of 1s, need coefficients to lie sum up to 1
                  beq = 1, lb, 1 # RHS of equality constraint, unit simplex if lb = 0, otherwise just sum up to 1
                  )


  return(weights.opt)
}
