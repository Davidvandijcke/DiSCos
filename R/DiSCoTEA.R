
#' @title DiSCoTEA
#' @description Distributional Synthetic Controls Treatment Effect Aggregator (DiSCoTEA)
#'
#' @details This function takes in the output of the DiSCo_per function and computes aggregate treatment effect using a user-specified aggregation statistic. The default is the average treatment effect (ATE).
#'
#' @param disco Output of the DiSCo function
#' @param agg String indicating the aggregation statistic to be used. Options are "ATT" (average treatment effect on the treated),
#' @param graph Boolean indicating whether to plot graphs. Default is TRUE.
#' @param n_per_window Integer indicating the number of periods to include in each plot window, if graph=TRUE. Default is NULL, which means that the entire time window is used.
#' This can lead to margin errors on the plot window if the number of time periods is large, in which case it is recommended to specify a smaller number (e.g. <=10).
#' @export
DiSCoTEA <- function(disco, agg="ATT", graph=TRUE, time=TRUE, n_per_window=NULL) {

  # reconstruct some parameters
  df <- disco$params$df
  t_max <- max(df$time_col)
  t_min <- min(df$time_col)
  t0 <- disco$params$t0
  T0 <- unique(df[time_col == t0]$t_col)  - 1
  T_max <- max(df$t_col)
  CI <- disco$params$CI
  cl <- disco$params$cl
  placebo <- disco$params$CI_placebo
  evgrid = seq(from=0,to=1,length.out=disco$params$M+1)
  qmethod <- disco$params$qmethod


  #
  if (is.null(n_per_window)) n_per_window <- T_max


  if (!placebo) {
    t_start <- t0
    T_start <- T0+1
  } else{
    t_start <- t_min
    T_start <- 1
  }

  if ((CI) & (length(disco$params$CI_periods) < length(1:T_max))) {
    stop("To aggregate treatment effect uncertainty, we need confidence intervals for all periods.")
  }

  # calculate quantile treatment effects
  qtiles_centered <- lapply(T_start:T_max,
                            function(x) disco$results.periods[[x]]$DiSCo$quantile - disco$results.periods[[x]]$target$quantiles)
  if (CI) { # calculate CI quantile treatment effects
    qtiles_centered_boot <- lapply(T_start:T_max, function(x) disco$results.periods[[x]]$DiSCo$CI$bootmat - disco$results.periods[[x]]$target$quantiles)
  }

  qtiles_centered <- lapply(T_start:T_max,
                   function(x) disco$results.periods[[x]]$target$quantiles)
  qtiles_centered_boot <- lapply(T_start:T_max,
                        function(x) disco$results.periods[[x]]$DiSCo$CI$bootmat)

  #---------------------------------------------------------------------------
  ### average treatment on the treated
  #---------------------------------------------------------------------------
  if (agg == "ATT") {

    avg_treat_time <- lapply(qtiles_centered, mean)

    if (time) { # return at each time point
      treats <- unlist(avg_treat_time)
    } else{
      treats <- as.vector(mean(unlist(avg_treat_time)[(T0+1):T_max]))
    }

    if (CI) { # get confidence intervals if they were calculated before
      avg_treat_time_boot <- sapply(qtiles_centered_boot, function(x) apply(x, 2, mean))

      if (time) { # return at each time point
        sds <- apply(avg_treat_time_boot, 2, sd)
        ci_lower <- apply(avg_treat_time_boot,2,myQuant, q=(1-cl)/2, qmethod)
        ci_upper <- apply(avg_treat_time_boot,2,myQuant, q=cl+(1-cl)/2, qmethod)
      } else {
        avg_treat_boot <- apply(avg_treat_time_boot[,T0:T_max], 1, mean)
        sds <- as.vector(sd(avg_treat_boot))
        ci_lower <- myQuant(avg_treat_boot, q=(1-cl)/2, qmethod=qmethod)
        ci_upper <- myQuant(avg_treat_boot, q=cl+(1-cl)/2, qmethod=qmethod)
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
        x <- t_start:t_max
        plot(x=x, y=treats, ylim=c(minx - 1/3*abs(minx), maxx + 1/3*abs(maxx)), type="l", xlab=xlab, ylab="ATT", xaxt="n")
        axis(1, at=x, las=2)
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
  } else if (agg == "cdfTreat"){

    treats <- list()
    grid <- seq(floor(quantile(unlist(qtiles_centered), 0.01)), ceil(quantile(unlist(qtiles_centered), 0.99)), length.out = disco$params$G)

    for (i in 1:length(qtiles_centered)) {
      cdff <- stats::ecdf(qtiles_centered[[i]])
      treats[[i]] <- cdff(grid_cdf)
    }
    if (CI) {
      cdf_boot <- list()
      sds <- list()
      ci_lower <- list()
      ci_upper <- list()

      for (i in 1:length(qtiles_centered)) {
        # apply ecdf to each column of qtiles_centered_boot[[i]]
        cdf_boot[[i]] <- apply(qtiles_centered_boot[[i]], 2, function(x) stats::ecdf(x)(grid_cdf))
        sds[[i]] <- apply(cdf_boot[[i]], 1, sd)
        ci_lower[[i]] <- apply(cdf_boot[[i]],1,myQuant, q=(1-cl)/2, qmethod=qmethod)
        ci_upper[[i]] <- apply(cdf_boot[[i]],1,myQuant, q=cl+(1-cl)/2, qmethod=qmethod)
      }
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA
    }
    if (graph) {
      plotDistOverTime(treats, grid, t_start, t_max, n_per_window, CI, ci_lower, ci_upper)
    }
  #---------------------------------------------------------------------------
  ### quantiles
  #---------------------------------------------------------------------------
  } else if (agg == "quantileTreat") {

    treats <-  qtiles_centered
    grid <- evgrid

    if (CI) {
      sds <- lapply(qtiles_centered_boot, function(x) apply(x, 1, sd))
      ci_lower <- lapply(qtiles_centered_boot, function(x) apply(x,1,myQuant, q=(1-cl)/2, qmethod=NULL))
      ci_upper <- lapply(qtiles_centered_boot, function(x) apply(x,1,myQuant, q=cl+(1-cl)/2, qmethod=NULL))
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA
    }

    if (graph) {
      ymin <- floor(quantile(unlist(qtiles_centered), 0.01))
      ymax <- ceil(quantile(unlist(qtiles_centered), 0.99))
      ymin <- min(unlist(qtiles_centered))
      ymax <- max(unlist(qtiles_centered))
      plotDistOverTime(treats, grid, t_start, t_max, n_per_window, CI, ci_lower, ci_upper, ylim=c(ymin, ymax), xlab="Quantile", ylab="Treatment Effect")
    }

  }

}

#' @title Plot distribution of treatment effects over time
#' @description Plot distribution of treatment effects over time
#' @param cdf_centered list of centered distributional statistics
#' @param grid_cdf grid
#' @param t_start start time
#' @param t_max maximum time
#' @param n_per_window number of windows per plot
#' @param CI logical indicating whether to plot confidence intervals
#' @param ci_lower lower confidence interval
#' @param ci_upper upper confidence interval
#' @return plot of distribution of treatment effects over time
plotDistOverTime <- function(cdf_centered, grid_cdf, t_start, t_max, n_per_window, CI, ci_lower, ci_upper, ylim=c(0,1), cdf=TRUE, xlab="Treatment Effect", ylab="CDF") {
  time_list <- t_start:t_max
  for (i in 1:length(cdf_centered)) {
    if (i %% n_per_window == 1) {
      # Adjust window dimensions as needed
      dev.new(width = 5, height = 5)
      par(mfrow=c(length(i:min(i+n_per_window-1, length(cdf_centered))),1))
      par(mar=c(2, 2, 2, 2))
      par(oma = c(4, 4, 1, 1))
    }
    plot(x=grid_cdf, y=cdf_centered[[i]], type="l", ylim=ylim, col="black", xlab="", ylab="")
    if (CI) {
      lines(x=grid_cdf, y=ci_lower[[i]], col="grey")
      lines(x=grid_cdf, y=ci_upper[[i]], col="grey")
    }
    abline(h=0, col="grey")
    if (cdf) {
      abline(h=1, col="grey")
    }
    t <- time_list[i]
    if (t >=t0) { title(paste("t =", t, sep=" "), col.main = "red") } else { title(paste("t =", t, sep=" "), col.main="blue") }
    mtext(xlab, side=1, padj=2, outer=TRUE)
    mtext(ylab, side=2, padj=-2, outer=TRUE)
  }
}



#' @title DiSCoT
#' @description D
DiSCoT <- function() {
  out <- list(distp=distp, distt=distt, rank=rank, p_values=p_values, p_overall=p_overall, J_1=J_1)
  class(out) <- "permut"
  return(out)
}
