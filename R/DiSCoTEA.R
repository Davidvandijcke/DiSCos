
#' @title DiSCoTEA
#' @description Distributional Synthetic Controls Treatment Effect Aggregator (DiSCoTEA)
#'
#' @details This function takes in the output of the DiSCo_per function and computes aggregate treatment effect using a user-specified aggregation statistic. The default is the average treatment effect (ATE).
#'
#' @param disco Output of the DiSCo function
#' @param agg String indicating the aggregation statistic to be used. Options are
#' \itemize{
#'  \item{\code{quantileDiff} }{Difference in quantiles between the target and the weighted average of the controls.}
#'  \item{\code{quantile} }{.}
#' @param graph Boolean indicating whether to plot graphs. Default is TRUE.
#' @param n_per_window Integer indicating the number of periods to include in each plot window, if graph=TRUE. Default is NULL, which means that the entire time window is used.
#' @param savePlots Boolean indicating whether to save the plots to the current working directory. The plot names will be [agg]_[start_year]_[end_year].pdf. The default is FALSE.
#' @param xlim Optional vector of length 2 indicating the x-axis limits of the plot. This and the `ylim` option can be useful to zoom in on the relevant parts of the distribution for fat-tailed distributions.
#' @param ylim Optional vector of length 2 indicating the y-axis limits of the plot.
#' This can lead to margin errors on the plot window if the number of time periods is large, in which case it is recommended to specify a smaller number (e.g. <=10).
#' @export
DiSCoTEA <- function(disco, agg="quantileDiff", graph=TRUE, time=TRUE, n_per_window=NULL, savePlots=FALSE,
                     xlim=NULL, ylim=NULL, samples=c(0.25, 0.5, 0.75)) {

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
  evgrid = seq(from=0,to=1,length.out=disco$params$G+1)
  qmethod <- disco$params$qmethod
  q_min <- disco$params$q_min
  q_max <- disco$params$q_max


  #
  if (is.null(n_per_window)) n_per_window <- T_max


  if (!placebo) {
    t_start <- t0
    T_start <- T0+1
  } else{
    t_start <- t_min
    T_start <- 1
  }


  ##  calculate quantile treatment effects
  qtiles_centered <- lapply(T_start:T_max,
                            function(x) disco$results.periods[[x]]$DiSCo$quantile - disco$results.periods[[x]]$target$quantiles)
  if (CI) { # calculate CI quantile treatment effects
    qtiles_centered_boot <- lapply(T_start:T_max, function(x) disco$results.periods[[x]]$DiSCo$CI$bootmat - disco$results.periods[[x]]$target$quantiles)

    qtiles_boot <- lapply(T_start:T_max,
                          function(x) disco$results.periods[[x]]$DiSCo$CI$bootmat)
  }
  target_qtiles <- lapply(T_start:T_max,
                          function(x) disco$results.periods[[x]]$target$quantiles)

  ## calculate quantiles
  qtiles <- lapply(T_start:T_max,
                   function(x) disco$results.periods[[x]]$DiSCo$quantile)
  target_qtiles <- lapply(T_start:T_max,
                          function(x) disco$results.periods[[x]]$target$quantiles)

  #---------------------------------------------------------------------------
  ### difference of cdfs
  #---------------------------------------------------------------------------
  if (agg == "cdfDiff"){

    treats <- list()
    treats_boot <- list()
    grid <- disco$results.periods[[1]]$target$grid ## TODO: in future probably wanna address this at root by always taking same grid in DiSco, and don't calculate CDFs

    for (i in 1:length(disco$results.periods)) {
      c_cdf <- stats::ecdf(disco$results.periods[[i]]$DiSCo$quantile)(grid)
      t_cdf <- stats::ecdf(disco$results.periods[[i]]$target$quantile)(grid)
      treats[[i]] <- c_cdf - t_cdf
      boot_cdf <- apply(disco$results.periods[[i]]$DiSCo$CI$bootmat, 2, function(x) stats::ecdf(x)(grid))
      treats_boot[[i]] <- sweep(boot_cdf, 1, t_cdf, "-")
    }

    if (CI){
      sds <- list()
      ci_lower <- list()
      ci_upper <- list()
      for (i in 1:length(disco$results.periods)) {
        sds[[i]] <- apply(treats_boot[[i]], 1, sd)
        ci_lower[[i]] <- apply(treats_boot[[i]],1,stats::quantile, probs=(1-cl)/2)
        ci_upper[[i]] <- apply(treats_boot[[i]],1,stats::quantile, probs=cl+(1-cl)/2)
      }
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA

    }
    if (graph) {
      if (is.null(ylim)) {
        ymin <- quantile(unlist(treats), 0.1)
        ymax <- quantile(unlist(treats), 0.99)
        ylim <- c(ymin, ymax)
      }
      if (is.null(xlim)) {
        xmin <- quantile(unlist(grid), 0.01)
        xmax <- quantile(unlist(grid), 0.99)
        xlim <- c(xmin, xmax)
      }
      plotDistOverTime(treats, grid, t_start, t_max, n_per_window, CI, ci_lower, ci_upper, savePlots=savePlots,
                       plotName=agg, ylim=ylim, xlim=xlim, xlab="Y", ylab="CDF Change")
    }
    #---------------------------------------------------------------------------
    ### cdfs
    #---------------------------------------------------------------------------

  } else if (agg == "cdf"){

    treats <- list()
    target_cdf <- list()
    grid <- seq(min(unlist(qtiles)), max(unlist(qtiles)), length.out = disco$params$G)
    # TODO: obsLine

    for (i in 1:length(disco$results.periods)) {
      treats[[i]] <- stats::ecdf(qtiles[[i]])(grid)
      target_cdf[[i]] <- stats::ecdf(disco$results.periods[[i]]$target$quantiles)(grid)
    }
    if (CI) {
      treats_boot <- list()
      sds <- list()
      ci_lower <- list()
      ci_upper <- list()

      for (i in 1:length(qtiles_centered)) {
        # apply ecdf to each column of qtiles_centered_boot[[i]]
        treats_boot[[i]] <- apply(qtiles_boot[[i]], 2, function(x) stats::ecdf(x)(grid))
        sds[[i]] <- apply(treats_boot[[i]], 1, sd)
        ci_lower[[i]] <- apply(treats_boot[[i]],1,stats::quantile, probs=(1-cl)/2)
        ci_upper[[i]] <- apply(treats_boot[[i]],1,stats::quantile, probs=cl+(1-cl)/2)
      }
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA
    }
    if (is.null(xlim)) xlim <- c(min(grid), max(grid))
    if (graph) {
      plotDistOverTime(treats, grid, t_start, t_max, n_per_window, CI, ci_lower, ci_upper, savePlots=savePlots, plotName=agg,
                       obsLine = target_cdf, xlab="Y", ylab="CDF", lty=1, lty_obs=2, xlim=xlim)
    }
    #---------------------------------------------------------------------------
    ### quantiles of treatment effects
    #---------------------------------------------------------------------------
  } else if (agg == "quantileDiff") {

    treats <-  qtiles_centered
    treats_boot <- qtiles_centered_boot
    grid <- evgrid

    if (CI) {
      sds <- lapply(treats_boot, function(x) apply(x, 1, sd))
      ci_lower <- lapply(treats_boot, function(x) apply(x,1,stats::quantile, probs=(1-cl)/2))
      ci_upper <- lapply(treats_boot, function(x) apply(x,1,stats::quantile, probs=cl+(1-cl)/2))
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA
    }

    if (graph) {
      if (is.null(ylim)) {
        ymin <- quantile(unlist(qtiles_centered), 0.01)
        ymax <- quantile(unlist(qtiles_centered), 0.99)
        ylim <- c(ymin, ymax)
      }
      if (is.null(xlim)) {
        xmin <- min(unlist(evgrid))
        xmax <- max(unlist(evgrid))
        xlim <- c(xmin, xmax)
      }
      plotDistOverTime(treats, grid, t_start, t_max, n_per_window, CI, ci_lower, ci_upper, ylim=ylim,
                       xlab="Quantile", ylab="Treatment Effect", cdf=FALSE, savePlots=savePlots, plotName=agg)
    }
    #---------------------------------------------------------------------------
    ### counterfactual and observed quantiles
    #---------------------------------------------------------------------------
  } else if (agg == "quantile") {

    treats <-  qtiles
    treats_boot <- qtiles_boot
    grid <- evgrid

    if (CI) {
      sds <- lapply(treats_boot, function(x) apply(x, 1, sd))
      ci_lower <- lapply(treats_boot, function(x) apply(x,1,stats::quantile, probs=(1-cl)/2))
      ci_upper <- lapply(treats_boot, function(x) apply(x,1,stats::quantile, probs=cl+(1-cl)/2))
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA
    }

    if (graph) {
      if (is.null(ylim)) {
        ymin <- quantile(unlist(qtiles), 0.01)
        ymax <- quantile(unlist(qtiles), 0.99)
        ylim <- c(ymin, ymax)
      }
      if (is.null(xlim)) {
        xmin <- quantile(unlist(grid), 0.01)
        xmax <- quantile(unlist(grid), 0.99)
        xlim <- c(xmin, xmax)
      }
      plotDistOverTime(treats, grid, t_start, t_max, n_per_window, CI, ci_lower, ci_upper, ylim=ylim, xlab="Quantile",
                       ylab="Treatment Effect", cdf=FALSE, obsLine = target_qtiles, savePlots=savePlots, plotName=agg)
    }

  }
  nams <- names(disco$results.periods)
  if (length(treats) == length(nams)) { # if we have time effects
    names(treats) <- nams
    names(sds) <- nams
    names(ci_lower) <- nams
    names(ci_upper) <- nams
  }

  if (agg %in% c("cdfDiff", "quantileDiff")){

  #---------------------------------------------------------------------------
  ## calculate the aggregated treatment effect at the sample points
  #---------------------------------------------------------------------------
    if (min(samples) >0) samples <- c(0, samples)
    if (max(samples) < 1) samples <- c(samples, 1)

    if (q_min != 0 | q_max != 1) {
      samples <- c(0, 1)
    }

    # get treatment periods
    t_list <- as.numeric(nams)
    # get treatment time periods
    I <- which(t_list >= t0)

    # for treatment time periods, loop over a sample of the distribution and bind to dataframe that we'll print
    grid_q <- samples * (length(grid)-1) + 1
    for (i in I) {
      for (j in 1:(length(grid_q)-1)){
        f <- grid_q[j]
        t <- grid_q[j+1]
        treats_agg <- mean(treats[[i]][f:t])
        treats_boot_agg <- apply(treats_boot[[i]][f:t,], 2, mean)
        sd_agg <- sd(treats_boot_agg)
        ci_lower_agg <- stats::quantile(treats_boot_agg, probs=(1-cl)/2)
        ci_upper_agg <- stats::quantile(treats_boot_agg, probs=cl+(1-cl)/2)

        out_temp <- cbind.data.frame(t=t_list[[i]], grid_lower=grid[f], grid_upper=grid[t], treats = treats_agg, ses = sd_agg,
                                     ci_lower = ci_lower_agg, ci_upper = ci_upper_agg)
        if ((i == I[1]) & (j==1)) {
          out <- out_temp
        } else {
          out <- rbind.data.frame(out, out_temp)
        }
      }
    }
    # header
    if (agg == "cdfDiff") {
      ooi <- "CDF \u0394" # object of interest
    } else if (agg == "quantileDiff") {
      ooi <- "Quantile \u0394"
    }
    # confidence band text
    cband_text1 <- paste0("[", 100*cl,"% ")

    # format the dataframe
    out <- round(out, 4)
    sig <- (out$ci_lower > 0) | (out$ci_upper < 0)
    sig_text <- ifelse(sig, "*", "")
    out <- cbind.data.frame(out, sig_text)
    ooi <- gsub(" \n", "", ooi)
    colnames(out) <- c("Time", "X_from", "X_to",  ooi, "Std. Error", cband_text1, "Conf. Band]", "")

    if (q_min != 0 | q_max != 1) { # if we have a subset of the distribution, we only calculate the effect there
      out$X_from <- q_min
      out$X_to <- q_max
    }


    ## redo permutation for samples
    if (!is.null(disco$perm)) { # if they asked for permutation test in main function

    }

  } else {
    out <- NULL # if they didn't ask for treatment effects
  }


  call <- match.call()
  return(DiSCoT(agg=agg, treats=treats, grid=grid, ses=sds, ci_lower=ci_lower, ci_upper=ci_upper,
                t0=t0, call=call, cl=cl, N=nrow(disco$params$df), J=uniqueN(disco$params$df$id_col)-1, agg_df=out,
                perm=disco$perm))

}


#' @title DiSCoT
#' @description S3 object holding aggregated treatment effects
#' @param agg aggregation method
#' @param treats list of treatment effects
#' @param ses list of standard errors
#' @param ci_lower list of lower confidence intervals
#' @param ci_upper list of upper confidence intervals
#' @param t0 start time
#' @param call call
#' @param cl confidence level
#' @param N number of observations
#' @param J number of treated units
#' @param grid grid
#' @param agg_df dataframe of aggregated treatment effects and their confidence intervals
#' @param perm list of permutation results
#' @export
DiSCoT <- function(agg, treats, ses, grid, ci_lower, ci_upper, t0, call, cl, N, J, agg_df, perm) {
  out <- list(agg=agg, treats=treats, ses=ses, ci_lower=ci_lower, ci_upper=ci_upper, t0=t0, call=call, cl=cl,
              grid=grid, N=N, J=J, agg_df=agg_df, perm=perm)
  class(out) <- "DiSCoT"
  out
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
plotDistOverTime <- function(cdf_centered, grid_cdf, t_start, t_max, n_per_window, CI, ci_lower, ci_upper, ylim=c(0,1), xlim=NULL,
                             cdf=TRUE, xlab="Treatment Effect", ylab="CDF",
                             obsLine = NULL, savePlots=FALSE, plotName=NULL, lty=1, lty_obs=1) {
  time_list <- t_start:t_max
  for (i in 1:length(cdf_centered)) {
    if (i %% n_per_window == 1) {
      # Adjust window dimensions as needed
      if (i > 1 & savePlots) {
        dev.off()
      }
      if (savePlots) {
        pdf(paste0(plotName, time_list[i], "_", time_list[i] + n_per_window-1, ".pdf"), width = 5, height = 5)
      } else {
        dev.new(width = 5, height = 5)
      }
      par(mfrow=c(length(i:min(i+n_per_window-1, length(cdf_centered))),1))
      par(mar=c(2, 2, 2, 2))
      par(oma = c(4, 4, 1, 1))
    }
    if (is.null(xlim)) {
      xlim <- c(min(grid_cdf), max(grid_cdf))
    }
    plot(x=grid_cdf, y=cdf_centered[[i]], type="l", ylim=ylim, xlim=xlim, col="black", xlab="", ylab="", lty=lty)
    if (CI) {
      lines(x=grid_cdf, y=ci_lower[[i]], col="grey", lty=lty)
      lines(x=grid_cdf, y=ci_upper[[i]], col="grey", lty=lty)
    }
    if (!is.null(obsLine)) {
      lines(x=grid_cdf, y=obsLine[[i]], col="blue", lty=lty_obs)
      if (i %% n_per_window == 1) {
        legend("topleft", legend = c("Counterfactual", "Observed"), col = c("black", "blue"), cex = 0.6, lty=1, bty="n")
      }
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





#' @title summary.DiSCoT
#' @description Summary of DiSCoT object
#' @param object DiSCoT object
#' @return summary of DiSCoT object
#' @export
summary.DiSCoT <- function(object) {

  # get results dataframe
  out <- object$agg_df
  # get treatment periods
  t_list <- as.numeric(names(object$treats))
  # get treatment time periods
  I <- which(t_list >= object$t0)

  #------------------------------------------------------------
  # print header
  #------------------------------------------------------------
  # call
  cat("\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")

  # citation
  citation()
  cat("\n")

  if (object$agg %in% c("quantile", "cdf")) {
    cat("No treatment effects to summarize, set graph=TRUE in function call or specify a treatment effect option in `agg`.")
    return(invisible(NULL))
  }


  #------------------------------------------------------------
  # print treatment effects
  #------------------------------------------------------------

  # header
  if (object$agg == "cdfDiff") {
      ooi <- "CDF \u0394 \n" # object of interest
  } else if (object$agg == "quantileDiff") {
    ooi <- "Quantile \u0394 \n"
  }
  cat(paste0("Aggregated Distribution Differences, ", ooi))




  print(out, row.names=FALSE)

  cat("---\n")
  cat("Signif. codes: `*' Confidence band for distribution differences does not cover 0")
  cat("\n\n")

  cat(paste(summary(object$perm), collapse="\n"))

  cat(paste0("Number of pre-treatment periods: ", length(t_list) - length(I)))
  cat("\n")

  cat(paste0("Number of post-treatment periods: ", length(I)))
  cat("\n")


  cat(paste0("N=", formatC(object$N, big.mark=",")))
  cat("\n")
}

