utils::globalVariables(c("y_col", "id_col", "time_col", "t_col", "group", "x", "y"))

#' @title Distributional Synthetic Controls
#'
#' @description This function implements the distributional synthetic controls (DiSCo) method from \insertCite{gunsilius2023distributional;textual}{DiSCos}.
#' as well as the alternative mixture of distributions approach.
#'
#' @details This function is called for every time period in the DiSCo function. It implements the DiSCo method for a single time period, as well as the mixture of distributions approach.
#' The corresponding results for each time period can be accessed in the `results.periods` list of the output of the DiSCo function. The DiSCo function returns the average weight for each unit across all periods,
#' calculated as a uniform mean, as well as the counterfactual target distribution produced as the weighted average of the control distributions for each period, using these averaged weights.
#'
#' @param df Data frame or data table containing the distributional data for the target and control units. The data table should contain the following columns:
#' \itemize{
#' \item \code{y_col } A numeric vector containing the outcome variable for each unit. Units can be individuals, states, etc., but they should be nested within a larger unit (e.g. individuals or counties within a state)
#' \item \code{id_col } A numeric vector containing the aggregate IDs of the units. This could be, for example, the state if the units are counties or individuals
#' \item \code{time_col } A vector containing the time period of the observation for each unit. This should be a monotonically increasing integer.
#' }
#' @param id_col.target Variable indicating the name of the target unit, as specified in the id_col column of the data table.
#' This variable can be any type, as long as it is the same type as the id_col column of the data table.
#' @param t0 Integer indicating period of treatment.
#' @param M Integer indicating the number of control quantiles to use in the DiSCo method. Default is 1000.
#' @param G Integer indicating the number of grid points for the grid on which the estimated functions are evaluated. Default is 1000.
#' @param num.cores Integer, number of cores to use for parallel computation. Default is 1. If the `permutation` or `CI` arguments are set to TRUE, this can be slow and it is recommended to set this to 4 or more, if possible.
#' If you get an error in "all cores" or similar, try setting num.cores=1 to see the precise error value.
#' @param permutation Logical, indicating whether to use the permutation method for computing the optimal weights. Default is FALSE.
#' @param q_min Numeric, minimum quantile to use. Set this together with `q_max` to restrict the range of quantiles used to construct the synthetic control.
#' Default is 0 (all quantiles). Currently NOT implemented for the `mixture` approach.
#' @param q_max Numeric, maximum quantile to use. Set this together with `q_min` to restrict the range of quantiles used to construct the synthetic control.
#' Default is 1 (all quantiles). Currently NOT implemented for the `mixture` approach.
#' @param CI Logical, indicating whether to compute confidence intervals for the counterfactual quantiles. Default is FALSE.
#' The confidence intervals are computed using the bootstrap procedure described in \insertCite{vandijcke2024rto;textual}{DiSCos}.
#' @param boots Integer, number of bootstrap samples to use for computing confidence intervals. Default is 500.
#' @param replace Logical, indicating whether to sample with replacement when computing the bootstrap samples. Default is TRUE.
#' @param uniform Logical, indicating whether to construct uniform bootstrap confidence intervals. Default is FALSE
#' If FALSE, the confidence intervals are pointwise.
#' @param cl Numeric, confidence level for the (two-sided) confidence intervals.
#' @param graph Logical, indicating whether to plot the permutation graph as in Figure 3 of the paper. Default is FALSE.
#' @param qmethod Character, indicating the method to use for computing the quantiles of the target distribution. The default is NULL, which uses the \code{\link[stats]{quantile}} function from the stats package.
#' Other options are "\code{\link[evmix]{qkden}}" (based on smoothed kernel density function) and "\code{\link[extremeStat:distLquantile]{extreme}}" (based on parametric extreme value distributions).
#' Both are substantially slower than the default method but may be useful for fat-tailed distributions with few data points at the upper quantiles. Alternatively, one could use the q_max option to restrict the range of quantiles used.
#' @param qtype Integer, indicating the type of quantile to compute when using  \code{\link[stats]{quantile}} in the `qmethod` argument.
#' The default 7. See the documentation for the \code{\link[stats]{quantile}} function for more information.
#' @param seed Integer, seed for the random number generator. This needs to be set explicitly in the function call, since it will invoke \code{\link[base]{RNGkind}} which will set the seed for each core
#' when using parallel processes. Default is NULL, which does not set a seed.
#' @param simplex Logical, indicating whether to use to constrain the optimal weights to the unit simplex. Default is FALSE, which only constrains the weights to sum up to 1 but allows them to be negative.
#' @param mixture Logical, indicating whether to use the mixture of distributions approach instead.
#' See Section 4.3. in \insertCite{gunsilius2023distributional;textual}{DiSCos}. This approach minimizes the distance between the CDFs
#' instead of the quantile functions, and is preferred for categorical variables. When working with such variables, one should
#' also provide a list of support points in the `grid.cat` parameter. When that is provided, this parameter is automatically set to TRUE. Default is FALSE.
#' @param grid.cat List, containing the discrete support points for a discrete grid to be used with the mixture of distributions approach.
#' This is useful for constructing synthetic distributions for categorical variables. Default is NULL, which uses a continuous grid based on the other parameters.
#' @return A list containing the following elements:
#' \itemize{
#' \item \code{results.periods} A list containing, for each time period, the elements described in the return argument of \code{\link{DiSCo_iter}}, as well as the following additional elements:
#' \itemize{
#'  \item \code{DiSco}
#'  \itemize{
#'  \item \code{quantile } The counterfactual quantiles for the target unit.
#'  \item \code{weights } The optimal weights for the target unit.
#'  \item \code{cdf } The counterfactual CDF for the target unit.
#'  }
#'  }
#' \item \code{weights} A numeric vector containing the synthetic control weights for the control units, averaged over time.
#' When `mixture` is TRUE, these are the weights for the mixture of distributions, otherwise they are the weights for the quantile-based approach.
#' \item \code{CI} A list containing the confidence intervals for the counterfactual quantiles and CDFs, if `CI` is TRUE.
#' Each element contains two named subelements called `upper`, `lower`, `se` which
#' are the upper and lower confidence bands and the standard error of the estimate, respectively.
#' They are G x T matrices where G is the specified number of grid points and T is the number of time periods.
#' The elements are:
#' \itemize{
#' \item \code{cdf} The bootstrapped CDF
#' \item \code{quantile} The bootstrapped quantile
#' \item \code{quantile_diff} The bootstrapped quantile difference
#' \item \code{cdf_diff} The bootstrapped CDF difference
#' \item \code{bootmat} A list containing the raw bootstrapped samples for the counterfactual quantiles and CDFs, if `CI` is TRUE.
#' These are not meant to be accessed directly, but are used by `DiSCoTEA` to compute aggregated standard errors. Advanced users
#' may wish to access these directly for further analysis. The element names should be self-explanatory.
#' #' \item \code{control_ids} A list containing the control unit IDs used for each time period, which can be used to identify the weights
#' associated with each control as the returned weights have the same order as the control IDs.
#' \item \code{perm } A \code{\link{permut}} object containing the results of the permutation method, if `permutation` is TRUE.
#' Call `summary` on this object to print the overall results of the permutation test.
#' #' \item \code{evgrid} A numeric vector containing the grid points on which the quantiles were evaluated.
#' \item \code{params} A list containing the parameters used in the function call.
#' }
#' }
#' @importFrom Rdpack reprompt
#' @importFrom stats sd quantile rnorm
#' @import data.table ggplot2
#' @references
#'  \insertAllCited()
#' @export
#'
DiSCo <- function(df, id_col.target, t0, M = 1000, G = 1000, num.cores = 1, permutation = FALSE, q_min = 0, q_max = 1,
                  CI = FALSE, boots = 500, replace=TRUE, uniform=FALSE, cl = 0.95, graph = FALSE,
                  qmethod=NULL, qtype=7, seed=NULL, simplex=FALSE, mixture=FALSE, grid.cat=NULL) {

  #---------------------------------------------------------------------------
  ### process inputs
  #---------------------------------------------------------------------------
  # TODO: need to fix discontinuous quantile method for mixture, I think
  # make sure we have a data table
  df <- data.table::as.data.table(df)

  # check the inputs
  checks(df, id_col.target, t0, M, G, num.cores, permutation, q_min, q_max,
         CI, boots, cl, graph, # TODO: update this
         qmethod, seed)

  if (!is.null(grid.cat)) { # TODO: allow restricted grid for mixture
    mixture <- TRUE
    grid.cat <- grid.cat[order(grid.cat)] # order the discrete grid
  }
  df_pres <- data.table::copy(df)

  # if restricted quantile range, subset the data
  # I am passing q_min = 0, q_min = 1 to the other functions below, which is a legacy feature but it's useful to have
  if ((q_min >0) | (q_max < 1)) {
    df[, quantile := data.table::frank(y_col, ties.method = "average") / .N, by = c("id_col", "time_col")]
    df <- df[(quantile >= q_min) & (quantile <= q_max)]
  }



  if (!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG") # to make sure our seeds are the same across parallel processes
    set.seed(seed)
  }


  # create a column for the normalized time period
  t_min <- min(df$time_col)
  df[, t_col := time_col - t_min + 1]
  T0 <- unique(df[time_col == t0]$t_col)  - 1
  T_max <- max(df$t_col)



  # create a list to store the results for each period
  results.periods <- list()

  evgrid = seq(from=0,to=1,length.out=G+1)

  #---------------------------------------------------------------------------
  ###  run the main function in parallel for each period
  #---------------------------------------------------------------------------
  periods <- sort(unique(df$t_col)) # we call the iter function on all periods, but won't calculate weights for the post-treatment periods
                                    # TODO: currently we are actually calculating them, but it's not a huge cost
  controls.id <- unique(df[id_col != id_col.target]$id_col) # list of control ids
  results.periods <- mclapply.hack(periods, DiSCo_iter, df, evgrid, id_col.target = id_col.target, M = M,
                                   G = G, T0 = T0, mc.cores = num.cores, qmethod=qmethod, qtype=qtype, q_min=0, q_max=1,
                                   controls.id=controls.id, simplex=simplex, grid.cat, mixture)

  # turn results.periods into a named list where the name is the period
  names(results.periods) <- as.character(periods)


  #---------------------------------------------------------------------------
  ### calculate average weights across pre-periods
  #---------------------------------------------------------------------------
  Weights_DiSCo_avg <- 0
  Weights_mixture_avg <- 0
  for (yy in 1:T0){
      Weights_DiSCo_avg <- Weights_DiSCo_avg + results.periods[[yy]]$DiSCo$weights
      Weights_mixture_avg <- Weights_mixture_avg + results.periods[[yy]]$mixture$weights
  }
  Weights_DiSCo_avg <- (1/T0) * Weights_DiSCo_avg
  Weights_mixture_avg <- (1/T0) * Weights_mixture_avg


  #---------------------------------------------------------------------------
  ### calculate counterfactual quantile functions and cdfs across pre-periods
  #---------------------------------------------------------------------------
  if (!mixture) {
    bc <- lapply(seq(1:T_max), function(x) DiSCo_bc(controls.q=results.periods[[x]]$controls$quantiles, Weights_DiSCo_avg, evgrid))
  } else {
   cdf <- lapply(seq(1:T_max), function(x) results.periods[[x]]$controls$cdf[,-1] %*% as.vector(Weights_mixture_avg) )
  }

  # calculating the counterfactual target quantiles and CDF
  for (x in seq(1:T_max)) {
    if (!mixture) { # if DiSCo
      bc_x <- bc[[x]]
      results.periods[[x]]$DiSCo$quantile <- bc_x
      grid_ord <- results.periods[[x]]$target$grid
      cdff <- stats::ecdf(bc_x)
      DiSCo_cdf <- cdff(grid_ord)
      results.periods[[x]]$DiSCo$cdf <- DiSCo_cdf
    } else { # if mixture
      cdf_x <- cdf[[x]]
      results.periods[[x]]$DiSCo$cdf <- cdf_x
      grid_ord <- results.periods[[x]]$target$grid
      bc_x <- sapply(evgrid, function(y) grid_ord[which(cdf_x >= (y-(1e-5)))[1]]) # tolerance accounts for inaccuracy (esp != 1)
      results.periods[[x]]$DiSCo$quantile <- bc_x
    }
  }


  if (mixture) { # choose weights for resampling and permutation test
    weights <- Weights_mixture_avg
  } else {
    weights <- Weights_DiSCo_avg
  }

  #---------------------------------------------------------------------------
  ### confidence intervals
  #---------------------------------------------------------------------------


  if (!is.null(grid.cat)) grid <- grid.cat
  if (CI) {
    # extract time-specific data
    controls <- lapply(results.periods, function(x) x$controls$data)
    target <- lapply(results.periods, function(x) x$target$data)
    grid <- lapply(results.periods, function(x) x$target$grid)
    q_disco <- lapply(results.periods, function(x) x$DiSCo$quantile)
    cdf_disco <- lapply(results.periods, function(x) x$DiSCo$cdf)
    q_obs <- lapply(results.periods, function(x) x$target$quantiles)
    cdf_obs <- lapply(results.periods, function(x) x$target$cdf)

    ## bootstrap the estimator
    CI_bootmat <-  mclapply.hack(1:boots, DiSCo_CI, controls=controls, target=target,
                              T_max=T_max, T0=T0, grid=grid, mc.cores=num.cores,
                              evgrid=evgrid, mixture=mixture, M=M, simplex=simplex,
                              qmethod=qmethod, qtype=qtype, replace=replace)
    CI_out <- parseBoots(CI_bootmat, cl, q_disco, cdf_disco, q_obs, cdf_obs, uniform)

  } else {
    CI_out <- NULL
  }


  #---------------------------------------------------------------------------
  ### permutation test
  #---------------------------------------------------------------------------
  if (permutation) {

    # run the permutation test
    perm_obj <- DiSCo_per(results.periods=results.periods, evgrid=evgrid, T0=T0,
                          weights=weights, num.cores=num.cores,
                          graph=graph, qmethod=qmethod, qtype=qtype, M=M, q_min=q_min, q_max=q_max,
                          mixture=mixture, simplex=simplex)


  } else {
    perm_obj <- NULL
  }



  # rename periods for user convenience
  names(results.periods) <- t_min  + seq(1:T_max) - 1

  if (mixture) {
    weights <- Weights_mixture_avg
  } else {
    weights <- Weights_DiSCo_avg
  }

  #---------------------------------------------------------------------------
  return(list(results.periods=results.periods, weights=weights, CI=CI_out,
              control_ids=controls.id, perm=perm_obj, evgrid=evgrid,
              params=list(df=df_pres, id_col.target=id_col.target,
              t0=t0, M=M, G=G, CI=CI, cl=cl, qmethod=qmethod,
              boot=boots, q_min=q_min, q_max=q_max))
         )

}
