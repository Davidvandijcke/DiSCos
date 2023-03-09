
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
#' \item{t_col}{A vector containing the time period of the observation for each unit. This should be an increasing integer, and the first period should be equal to 1.}
#' }
#' @param id_col.target Variable indicating the name of the target unit, as specified in the id_col column of the data table.
#' This variable can be any type, as long as it is the same type as the id_col column of the data table.
#' @param T0 Integer indicating last period before treatment as counted from 1 (e.g, if treatment year 2003 was the 6th year in the sample, this parameter should be 5).
#' @param M Integer indicating the number of control units to use in the DiSCo method. Default is 1000.
#' @param G Integer indicating the number of grid points for the grid on which the estimated functions are evaluated. Default is 1000.
#' @param num.cores Integer, number of cores to use for parallel computation. Default is 1 (sequential computation). If the `permutation` argument is set to TRUE, this can be very slow!
#' @param permutation Logical, indicating whether to use the permutation method for computing the optimal weights. Default is FALSE.
#'
#'
DiSCo <- function(df, id_col.target, T0, M = 1000, G = 1000, num.cores = 1, permutation = FALSE) {

  # make sure we have a data table
  df <- as.data.table(dt)

  # check the inputs
  checks(df, id_col.target, T0, M, G, num.cores, permutation)

  # create a list to store the results for each period
  results.periods <- list()

  # run the main function in parallel for each period (mclapply doesn't require a cluster so we will run this in parallel if num.cores > 1, but it won't work on windows)
  start_time <- Sys.time()
  periods <- sort(unique(df$t_col)) # we call the iter function on all periods, but won't calculate weights for the post-treatment periods
  results.years <- lapply(periods, DiSCo_iter, df = df, id_col.target = id_col.target, M = M, G = G, T0 = T0) #, mc.cores = num.cores)
  end <- Sys.time()
  print(end - start_time)

  # turn results.years into a named list where the name is the year
  names(results.years) <- sort(unique(as.character(df$t_col)))


  #####
  #obtaining the weights as a uniform mean over all time periods
  Weights_DiSCo_avg <- results.years[[1]][[1]][[1]]
  Weights_mixture_avg <- results.over.years[[1]][[2]][[1]]
  for (yy in 2:7){
    Weights_DiSCo_avg <- Weights_DiSCo_avg + results.over.years[[yy]][[1]][[1]]
    Weights_mixture_avg <- Weights_mixture_avg + results.over.years[[yy]][[2]][[1]]
  }
  Weights_DiSCo_avg <- (1/length(1:7)) * Weights_DiSCo_avg
  Weights_mixture_avg <- (1/length(1:7)) * Weights_mixture_avg

  year.to.plot <- 6
  # generating the counterfactuals
  DiSCo_res2 <- DiSCo_bc(results.over.years[[year.to.plot]][[4]][[1]],
                    Weights_DiSCo_avg,seq(from=0,to=1,length.out=1001))

  # getting the CDF from the quantile function
  DiSCo_res2.cdfF <- ecdf(DiSCo_res2[[2]])
  DiSCo_res2.cdf <- DiSCo_res2.cdfF(results.over.years[[year.to.plot]][[3]][[3]])

  parallel::mclapply(1998:2014, DiSCo_iter)

  #####
  # Solving for the optimal weights in the DiSCo method and the alternative method using mixtures of
  # distributions by year

  # create named list for storing results, one element for each year
  results.by.year <- list()

  # looping over the years from 1998 - 2014




}
