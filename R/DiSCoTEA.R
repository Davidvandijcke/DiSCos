
#' @title DiSCoTEA
#' @description Distributional Synthetic Controls Treatment Effect Aggregator (DiSCoTEA)
#'
#' @details This function takes in the output of the DiSCo_per function and computes aggregate treatment effect using a user-specified aggregation statistic. The default is the average treatment effect (ATE).
#'
#' @param DiSCo_per_output Output of the DiSCo function
#' @param agg String indicating the aggregation statistic to be used. Options are "ATT" (average treatment effect on the treated),
DiSCoTEA <- function(disco, agg="ATT", graph=TRUE){

  # get treatment time window
  df <- disco$params$df
  t_max <- max(df$time_col)
  t_min <- min(df$time_col)
  t0 <- disco$params$t0
  T0 <- unique(df[time_col == t0]$t_col)  - 1
  T_max <- max(df$t_col)
  CI <- disco$params$CI


  if (agg == "ATT") { # average treatment on the treated
    qtiles_centered <- lapply(1:T_max,
                              function(x) disco$results.periods[[x]]$DiSCo$quantile - disco$results.periods[[x]]$target$quantiles)

    avg_treat_time <- lapply(qtiles_centered, mean)

    if (time) { # return at each time point
      treats <- unlist(avg_treat_time)
    } else{
      treats <- list(mean(unlist(avg_treat_time)[T0:T_max]))
    }

    if (CI) { # get confidence intervals if they were calculated before
      qtiles_centered_boot <- lapply(1:T_max, function(x) disco$results.periods[[x]]$DiSCo$CI$bootmat - disco$results.periods[[x]]$target$quantiles)
      avg_treat_time_boot <- sapply(qtiles_centered_boot, function(x) apply(x, 2, mean))

      if (time) { # return at each time point
        sds <- apply(avg_treat_time_boot, 2, sd)

      } else {
        avg_treat_boot <- apply(avg_treat_time_boot[,T0:T_max], 1, mean)
        sds <- list(sd(avg_treat_boot))
      }
    }
    if (graph) {
      if (CI) {
        x <- t_min:t_max
        plot(x=x, y=treats, ylim=c(min(treats - 1.96*sds), max(treats + 1.96*sds)), type="l", xlab="Time", ylab="ATT", main="ATT over time")
        lines(x=x, y=treats + 1.96*sds, col="grey")
        lines(x=x, y=treats - 1.96*sds, col="grey")
        abline(v=t0, col="red", lty="longdash")
        # q: how do I plot a long-dashed vertical line at t0?
        # a: abline(v=t0, col="red", linetype="dashed")



      } else {
        plot(treats, type="l", xlab="Time", ylab="ATT", main="ATT over time")
      }
    }

  } else if (agg == "cdf"){

    avg_treat_time <- lapply(qtiles_centered, mean)

    if (CI) {
      avg_treat_time_boot <- lapply(qtiles_centered_boot, function(x) apply(x, 2, mean))
      sds <- lapply(avg_treat_time_boot, sd)
    }
  }

}

#' @title DiSCoT
#' @description D
DiSCoT <- function() {
  out <- list(distp=distp, distt=distt, rank=rank, p_values=p_values, p_overall=p_overall, J_1=J_1)
  class(out) <- "permut"
  return(out)
}
