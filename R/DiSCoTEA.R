
#' @title DiSCoTEA
#' @description Distributional Synthetic Controls Treatment Effect Aggregator (DiSCoTEA)
#'
#' @details This function takes in the output of the DiSCo_per function and computes aggregate treatment effect using a user-specified aggregation statistic. The default is the average treatment effect (ATE).
#'
#' @param disco Output of the DiSCo function
#' @param agg String indicating the aggregation statistic to be used. Options are "ATT" (average treatment effect on the treated),
#' @param graph Boolean indicating whether to plot graphs. Default is TRUE.
#' @export
DiSCoTEA <- function(disco, agg="ATT", graph=TRUE, time=TRUE){

  # get treatment time window
  df <- disco$params$df
  t_max <- max(df$time_col)
  t_min <- min(df$time_col)
  t0 <- disco$params$t0
  T0 <- unique(df[time_col == t0]$t_col)  - 1
  T_max <- max(df$t_col)
  CI <- disco$params$CI
  cl <- disco$params$cl

  if ((CI) & (length(disco$params$CI_periods) < length(1:T_max))) {
    stop("To aggregate treatment effect uncertainty, we need confidence intervals for all periods.")
  }

  #---------------------------------------------------------------------------
  ### average treatment on the treated
  #---------------------------------------------------------------------------
  if (agg == "ATT") {
    qtiles_centered <- lapply(1:T_max,
                              function(x) disco$results.periods[[x]]$DiSCo$quantile - disco$results.periods[[x]]$target$quantiles)

    avg_treat_time <- lapply(qtiles_centered, mean)

    if (time) { # return at each time point
      treats <- unlist(avg_treat_time)
    } else{
      treats <- as.vector(mean(unlist(avg_treat_time)[T0:T_max]))
    }

    if (CI) { # get confidence intervals if they were calculated before
      qtiles_centered_boot <- lapply(1:T_max, function(x) disco$results.periods[[x]]$DiSCo$CI$bootmat - disco$results.periods[[x]]$target$quantiles)
      avg_treat_time_boot <- sapply(qtiles_centered_boot, function(x) apply(x, 2, mean))

      if (time) { # return at each time point
        sds <- apply(avg_treat_time_boot, 2, sd)
        ci_lower <- apply(avg_treat_time_boot,2,stats::quantile, probs=(1-cl)/2)
        ci_upper <- apply(avg_treat_time_boot,2,stats::quantile, probs=cl+(1-cl)/2)
      } else {
        avg_treat_boot <- apply(avg_treat_time_boot[,T0:T_max], 1, mean)
        sds <- as.vector(sd(avg_treat_boot))
        ci_lower <- stats::quantile(avg_treat_boot, probs=(1-cl)/2)
        ci_upper <- stats::quantile(avg_treat_boot, probs=cl+(1-cl)/2)
      }
    } else {
      sds <- list(NA)
      ci_lower <- NA
      ci_upper <- NA
    }
    if (graph) {
      minx <- minx <- min(c(treats, ci_lower), na.rm = TRUE)
      maxx <- max(c(treats, ci_upper), na.rm = TRUE)
      if (time) {
        xlab <- "Time"
        x <- t_min:t_max
        plot(x=x, y=treats, ylim=c(minx - 1/3*abs(minx), maxx + 1/3*abs(maxx)), type="l", xlab=xlab, ylab="ATT")
        if (CI) {
          lines(x=x, y=ci_lower, col="grey")
          lines(x=x, y=ci_upper, col="grey")
        }
        abline(v=t0, col="red", lty="longdash")
      } else {
        xlab <- ""
        x <- 1
        # plot point without x axis ticks
        plot(x=x, y=treats, ylim=c(minx - 1/3*abs(minx), maxx + 1/3*abs(maxx)), type="p", xlab=xlab, ylab="ATT", xaxt="n")
        # plot whiskers around point
        segments(x0=x, y0=ci_lower, x1=x, y1=ci_upper, col="grey")
      }
    }

  #---------------------------------------------------------------------------
  ### cdfs
  #---------------------------------------------------------------------------
  } else if (agg == "cdf"){

    cdf_centered <- list()
    for (i in 1:T_max) {
      cdff <- stats::ecdf(qtiles_centered[[i]])
      cdf_centered[[i]] <- cdff(disco$results.periods[[i]]$target$grid)
    }
    if (CI) {
      cdf_boot <- list()
      for (i in 1:T_max) {
        # apply ecdf to each column of qtiles_centered_boot[[i]]
        cdf_boot[[i]] <- apply(qtiles_centered_boot[[i]], 2, function(x) stats::ecdf(x)(disco$results.periods[[i]]$target$grid))
      }
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
