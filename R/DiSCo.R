
#' @title DiSCo
#'
#' @description This function implements the distributional synthetic controls (DiSCo) method from \insertref{gunsilius2020distributional}{DiSCo},
#' as well as the alternative mixture of distributions approach.
#'
#' @details This function is called for every time period in the DiSCo function. It implements the DiSCo method for a single time period, as well as the mixture of distributions approach.
#' The corresponding results for each year can be accessed in the `years` list of the output of the DiSCo function. The DiSCo function returns the average weight for each unit across all years,
#' calculated as a uniform mean, as well as the counterfactual target distribution produced as the weighted average of the control distributions for each year, using these averaged weights.
#'
#' @param dt Data frame or data table containing the distributional data for the target and control units. The data table should contain the following columns:
#' \itemize{
#' \item{y_col}{A numeric vector containing the outcome variable for each unit. Units can be individuals, states, etc., but they should be nested within a larger unit (e.g. individuals or counties within a state)}
#' \item{id_col}{A numeric vector containing the aggregate IDs of the units. This could be, for example, the state if the units are counties or individuals}
#' \item{time_col}{A vector containing the time period of the observation for each unit. This should be an increasing integer, and the first period should be equal to 1.}
#' }
#' @param id_col.target Variable indicating the name of the target unit, as specified in the id_col column of the data table.
#' This variable can be any type, as long as it is the same type as the id_col column of the data table.
#' @param t0 Integer indicating period of treatment.
#' @param M Integer indicating the number of control units to use in the DiSCo method. Default is 1000.
#' @param G Integer indicating the number of grid points for the grid on which the estimated functions are evaluated. Default is 1000.
#' @param num.cores Integer, number of cores to use for parallel computation. Default is 1 (sequential computation). If the `permutation` argument is set to TRUE, this can be very slow!
#' @param permutation Logical, indicating whether to use the permutation method for computing the optimal weights. Default is FALSE.
#'
#'
DiSCo <- function(df, id_col.target, t0, M = 1000, G = 1000, num.cores = 1, CI=FALSE, boots=500, cl=0.95, permutation = FALSE,
                  CI_periods = NULL, graph=FALSE) {

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
  start_time <- Sys.time()
  set.seed(1860)

  periods <- sort(unique(df$t_col)) # we call the iter function on all periods, but won't calculate weights for the post-treatment periods

  results.periods <- mclapply.hack(periods, DiSCo_iter, df, evgrid, id_col.target = id_col.target, M = M, G = G, T0 = T0, mc.cores = num.cores)

  end <- Sys.time()
  print(end - start_time)


  # turn results.years into a named list where the name is the year
  names(results.periods) <- as.character(periods)


  #####
  #obtaining the weights as a uniform mean over all time periods
  Weights_DiSCo_avg <- results.periods$`1`$DiSCo_weights
  Weights_mixture_avg <- results.periods$`1`$mixture$weights
  for (yy in 2:T0){
    Weights_DiSCo_avg <- Weights_DiSCo_avg + results.periods[[yy]]$DiSCo_weights
    Weights_mixture_avg <- Weights_mixture_avg + results.periods[[yy]]$mixture$weights
  }
  Weights_DiSCo_avg <- (1/T0) * Weights_DiSCo_avg
  Weights_mixture_avg <- (1/T0) * Weights_mixture_avg

  # calculating the counterfactual target distribution
  bc <- lapply(seq(1:T_max), function(x) DiSCo_bc(results.periods[[x]]$controls$data, results.periods[[x]]$controls.q, Weights_DiSCo_avg, evgrid))
  names(bc) <- t_min  + seq(1:T_max) - 1

  # calculate confidence intervals for selected time periods
  if (is.null(CI_periods) & CI) { # if wants CI for all periods
    CI_periods <- c(1, T0, T_max)
  } else if (!is.null(CI_periods)) { # if wants CI for specific period
    CI_periods <- sort(CI_periods) - t_min + 1
  }
  for (x in CI_periods) {
    controls <- results.periods[[x]]$controls$data
    bc_x <- bc[[x]]

    if (CI) {
      CI <- DiSCo_CI(controls=controls, bc=bc_x, weights=Weights_DiSCo_avg, mc.cores=num.cores, cl=cl, num.redraws=boots, evgrid = evgrid)
    } else {
      CI <- NULL
    }
    results.periods[[x]]$CI <- CI
  }

  # permutation tests
  if (permutation) {
    controls_per <- lapply(seq(1:T_max), function(x) results.periods[[x]]$controls$data)
    target_per <- lapply(seq(1:T_max), function(x) results.periods[[x]]$target$data)
    controls.q <- lapply(seq(1:T_max), function(x) results.periods[[x]]$controls.q)
    perm <- DiSCo_per(c_df=controls_per, t_df=target_per, controls.q=controls.q, T0=T0, weights=Weights_DiSCo_avg, num_cores=num.cores, evgrid=evgrid,
                      graph=graph)
  } else {
    perm <- NULL
  }

  distp <- perm$control.dist
  distt <- perm$target.dist
  perm_values <- DiSCo_per_rank(distt, distp)
  p_overall <- mean(unlist(perm_values$ranks)[((T0+1):T_max)]) / (length(distp) + 1)

  return(list(results.periods=results.periods, Weights_DiSCo_avg=Weights_DiSCo_avg,
              Weights_mixture_avg=Weights_mixture_avg, counterfactual=bc, perm=perm, p_overall=p_overall))

}
