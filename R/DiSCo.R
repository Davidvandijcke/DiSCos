
#' @title DiSCo
#'
#' @description This function implements the distributional synthetic controls (DiSCo) method from \insertCite{gunsilius2023distributional;textual}{DiSCo}.
#' as well as the alternative mixture of distributions approach.
#'
#' @details This function is called for every time period in the DiSCo function. It implements the DiSCo method for a single time period, as well as the mixture of distributions approach.
#' The corresponding results for each time period can be accessed in the `results.periods` list of the output of the DiSCo function. The DiSCo function returns the average weight for each unit across all periods,
#' calculated as a uniform mean, as well as the counterfactual target distribution produced as the weighted average of the control distributions for each period, using these averaged weights.
#'
#' @param df Data frame or data table containing the distributional data for the target and control units. The data table should contain the following columns:
#' \itemize{
#' \item{\code{y_col} }{A numeric vector containing the outcome variable for each unit. Units can be individuals, states, etc., but they should be nested within a larger unit (e.g. individuals or counties within a state)}
#' \item{\code{id_col} }{A numeric vector containing the aggregate IDs of the units. This could be, for example, the state if the units are counties or individuals}
#' \item{\code{time_col} }{A vector containing the time period of the observation for each unit. This should be an increasing integer, and the first period should be equal to 1.}
#' }
#' @param id_col.target Variable indicating the name of the target unit, as specified in the id_col column of the data table.
#' This variable can be any type, as long as it is the same type as the id_col column of the data table.
#' @param t0 Integer indicating period of treatment.
#' @param M Integer indicating the number of control quantiles to use in the DiSCo method. Default is 1000.
#' @param G Integer indicating the number of grid points for the grid on which the estimated functions are evaluated. Default is 1000.
#' @param num.cores Integer, number of cores to use for parallel computation. Default is 1 (sequential computation). If the `permutation` or `CI` arguments are set to TRUE, this can be very slow!
#' @param permutation Logical, indicating whether to use the permutation method for computing the optimal weights. Default is FALSE.
#' @param CI Logical, indicating whether to compute confidence intervals for the counterfactual quantiles. Default is FALSE.
#' @param CI_placebo Logical, indicating whether to compute confidence intervals for the pre-treatment periods. Default is TRUE.
#' If you have a lot of pre-treatment periods, setting this to FALSE can speed up the computation.
#' @param boots Integer, number of bootstrap samples to use for computing confidence intervals. Default is 500.
#' @param cl Numeric, confidence level for the (two-sided) confidence intervals.
#' @param graph Logical, indicating whether to plot the permutation graph as in Figure 3 of the paper. Default is FALSE.
#' @param qmethod Character, indicating the method to use for computing the quantiles of the target distribution. The default is NULL, which uses the quantile function from the \code{\}
#'
#' @return A list containing, for each time period, the elements described in the return argument of \code{\link{DiSCo_iter}}, as well as the following additional elements:
#' \itemize{
#'  \item{\code{DiSco}}{
#'  \itemize{
#'  \item{\code{CI} }{A list with the confidence intervals and standard errors for the counterfactual quantiles, if `CI` is TRUE and for the periods specified in `CI_periods`.
#'  See the output of \code{\link{DiSCo_CI}} for details.}
#'  \item{\code{quantile} }{The counterfactual quantiles for the target unit.}
#'  \item{\code{weights} }{The optimal weights for the target unit.}
#'  \item{\code{cdf} }{The counterfactual CDF for the target unit.}
#'  }
#'  }
#' \item{\code{perm} }{A \code{\link{permut}} object containing the results of the permutation method, if `permutation` is TRUE.
#' Call `summary` on this object to print the overall results of the permutation test.}
#' }
#' @importFrom Rdpack reprompt
#' @references
#'  \insertAllCited()
#' @export


DiSCo <- function(df, id_col.target, t0, M = 1000, G = 1000, num.cores = 1, permutation = FALSE,
                  CI = FALSE, CI_placebo=FALSE, boots = 500, cl = 0.95, CI_periods = NULL, graph = FALSE,
                  qmethod=NULL) {

  # make sure we have a data table
  df <- as.data.table(df)

  # check the inputs
  checks(df, id_col.target, t0, M, G, num.cores, permutation)

  # create a column for the normalized time period
  t_min <- min(df$time_col)
  df[, t_col := time_col - t_min + 1]
  T0 <- unique(df[time_col == t0]$t_col)  - 1
  T_max <- max(df$t_col)


  # create a list to store the results for each period
  results.periods <- list()

  evgrid = seq(from=0,to=1,length.out=M+1)

  # run the main function in parallel for each period
  periods <- sort(unique(df$t_col)) # we call the iter function on all periods, but won't calculate weights for the post-treatment periods
  results.periods <- mclapply.hack(periods, DiSCo_iter, df, evgrid, id_col.target = id_col.target, M = M,
                                   G = G, T0 = T0, mc.cores = num.cores, qmethod=qmethod)
  # turn results.periods into a named list where the name is the period
  names(results.periods) <- as.character(periods)


  ###############################
  #obtaining the weights as a uniform mean over all time periods
  Weights_DiSCo_avg <- results.periods$`1`$DiSCo$weights
  Weights_mixture_avg <- results.periods$`1`$mixture$weights
  for (yy in 2:T0){
    Weights_DiSCo_avg <- Weights_DiSCo_avg + results.periods[[yy]]$DiSCo$weights
    Weights_mixture_avg <- Weights_mixture_avg + results.periods[[yy]]$mixture$weights
  }
  Weights_DiSCo_avg <- (1/T0) * Weights_DiSCo_avg
  Weights_mixture_avg <- (1/T0) * Weights_mixture_avg

  bc <- lapply(seq(1:T_max), function(x) DiSCo_bc(controls=results.periods[[x]]$controls$data, controls.q=results.periods[[x]]$controls$quantiles, Weights_DiSCo_avg, evgrid))
  # calculating the counterfactual target quantiles and CDF
  for (x in seq(1:T_max)) {
    bc_x <- bc[[x]]
    results.periods[[x]]$DiSCo$quantile <- bc_x
    grid_ord <- results.periods[[x]]$target$grid
    cdff <- stats::ecdf(bc_x)
    DiSCo_cdf <- cdff(grid_ord)
    results.periods[[x]]$DiSCo$cdf <- DiSCo_cdf
  }

  # calculate confidence intervals for selected time periods
  if (CI & CI_placebo) { # if wants CI for all periods
    CI_periods <- seq(1, T_max)
  } else if (CI & !CI_placebo) {
    CI_periods <- seq(T0+1, T_max)
  }
  for (x in CI_periods) {
    cat(paste0("Computing confidence intervals for period: ", x + t_min - 1, "\n"))
    controls <- results.periods[[x]]$controls$data
    bc_x <- bc[[x]]

    if (CI) {
      CI_temp <- DiSCo_CI(controls=controls, bc=bc_x, weights=Weights_DiSCo_avg,
                          mc.cores=num.cores, cl=cl, num.redraws=boots, evgrid = evgrid, qmethod=qmethod)
    } else {
      CI_temp <- NULL
    }
    results.periods[[x]]$DiSCo$CI <- CI_temp
  }

  # permutation tests
  if (permutation) {
    controls_per <- lapply(seq(1:T_max), function(x) results.periods[[x]]$controls$data)
    target_per <- lapply(seq(1:T_max), function(x) results.periods[[x]]$target$data)
    controls.q <- lapply(seq(1:T_max), function(x) results.periods[[x]]$controls$quantiles)
    target.q <- lapply(seq(1:T_max), function(x) results.periods[[x]]$target$quantiles)
    perm <- DiSCo_per(c_df=controls_per, t_df=target_per, controls.q=controls.q, target.q=target.q, T0=T0, weights=Weights_DiSCo_avg, num_cores=num.cores, evgrid=evgrid,
                      graph=graph, qmethod=qmethod)
    distp <- perm$control.dist
    distt <- perm$target.dist
    perm_values <- DiSCo_per_rank(distt, distp)
    J_1 <- length(distp) + 1
    p_overall <- mean(unlist(perm_values$ranks)[((T0+1):T_max)]) / (J_1)

    perm_obj <- permut(distp, distt, perm_values$ranks, perm_values$p_t, p_overall, J_1)
  } else {
    perm_obj <- NULL
  }



  # rename periods for user convenience
  names(results.periods) <- t_min  + seq(1:T_max) - 1


  return(list(results.periods=results.periods, Weights_DiSCo_avg=Weights_DiSCo_avg,
              Weights_mixture_avg=Weights_mixture_avg, perm=perm_obj, params=list(df=df, id_col.target=id_col.target,
              t0=t0, M=M, G=G, CI=CI, CI_periods=CI_periods, cl=cl, CI_placebo=CI_placebo, qmethod=qmethod)))

}
