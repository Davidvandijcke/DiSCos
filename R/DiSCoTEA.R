
#' @title Aggregate treatment effects from DiSCo function.
#' @description Function to aggregate treatment effects from the output of the
#' DiSCo function, plot the distribution of the aggregation statistic over time,
#' and report summary tables.
#'
#' @details This function takes in the output of the DiSCo_per function and computes aggregate treatment effect using a user-specified aggregation statistic.
#' The default is the differences between the counterfactual and the observed quantile functions (`quantileDiff`). If `graph` is set to TRUE,
#' the function will plot the distribution of the aggregation statistic over time. The S3 class returned by the function
#' has a `summary` property that will print a selection of aggregated effects (specified by the `samples` parameter) for the chosen `agg` method, by post-treatment year (see examples below).
#' This `summary` call will only print effects if the `agg` parameter requested a distribution difference (`quantileDiff` or `cdfDiff`). The other aggregations are meant to be inspected visually.
#' If the `permutation` parameter was set to TRUE in the original `DiSCo` call, the summary table will include the results of the permutation test.
#' If the original `DiSCo` call was restricted to a range of quantiles smaller than `[0,1]` (i.e. `q_min` > 0 or `q_max` < 1), the `samples` parameter is ignored
#' and only the aggregated differences for the quantile range specified in the original call are returned.
#'
#' @param disco Output of the DiSCo function.
#' @param agg String indicating the aggregation statistic to be used. Options include
#' \itemize{
#'  \item \code{quantileDiff}  Difference in quantiles between the target and the weighted average of the controls.
#'  \item \code{quantile} Plots both the observed and the counterfactual quantile functions. No summary statistics will be produced.
#'  \item \code{cdfDiff}  Difference in CDFs between the target and the weighted average of the controls.
#'  \item \code{cdf} Plots both the observed and the counterfactual CDFs. No summary statistics will be produced.
#' }
#' @param graph Boolean indicating whether to plot graphs (default is TRUE).
#' @param t_plot Optional vector of time periods (`t_col` values in the original dataframe) to be plotted (default is NULL, which plots all time periods).
#' @param savePlots Boolean indicating whether to save the plots to the current working directory (default is FALSE). The plot names will be `[agg]_[start_year]_[end_year].pdf`.
#' @param xlim Optional vector of length 2 indicating the x-axis limits of the plot. Useful for zooming in on relevant parts of the distribution for fat-tailed distributions.
#' @param ylim Optional vector of length 2 indicating the y-axis limits of the plot.
#' @param samples Numeric vector indicating the range of quantiles of the aggregation statistic (`agg`) to be summarized in the `summary` property of the S3 class returned by the function (default is c(0.25, 0.5, 0.75)).
#' For example, if `samples` = c(0.25, 0.5, 0.75), the summary table will include the average effect for the 0-25th, 25-50th, 50-75th and 75-100th quantiles of the distribution of the aggregation statistic over time.
#' @return A \code{\link[DiSCos]{DiSCoT}} object, which is an S3 class that stores a list of treatment effects, their standard errors,
#' the corresponding confidence intervals (if specified), and a dataframe with treatment effects aggregated
#' according to the `agg` input. The S3 class also has a `summary` property that will print a selection of aggregated effects (specified by the `samples` parameter)
#' for the chosen `agg` method, by post-treatment year, as well as the permutation test results, if specified.
#'
#'
#' @importFrom stats sd
#' @export
DiSCoTEA <- function(disco, agg="quantileDiff", graph=TRUE, t_plot=NULL, savePlots=FALSE,
                     xlim=NULL, ylim=NULL, samples=c(0.25, 0.5, 0.75)) {

  # reconstruct some parameters
  df <- disco$params$df
  t_max <- max(df$time_col)
  t_min <- min(df$time_col)
  t0 <- disco$params$t0
  T0 <- unique(df[time_col == t0]$time_col)  - t_min
  T_max <- max(df$time_col) - t_min + 1
  CI <- disco$params$CI
  cl <- disco$params$cl
  evgrid = seq(from=0,to=1,length.out=disco$params$G+1)
  qmethod <- disco$params$qmethod
  q_min <- disco$params$q_min
  q_max <- disco$params$q_max




  t_start <- t_min
  T_start <- 1


  if (is.null(t_plot)) t_plot <- t_start:t_max



  ##  calculate quantile treatment effects
  qtiles_centered <- lapply(T_start:T_max,
                            function(x) disco$results.periods[[x]]$target$quantiles -  disco$results.periods[[x]]$DiSCo$quantile )
  if (CI) { # calculate CI quantile treatment effects
    qtiles_centered_boot <- lapply(T_start:T_max, function(x) disco$results.periods[[x]]$target$quantiles - disco$results.periods[[x]]$DiSCo$CI$bootmat)

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

    # get treatment effects
    for (i in 1:length(disco$results.periods)) { # TODO: left off here
      grid <- disco$results.periods[[i]]$target$grid
      c_cdf <- stats::ecdf(disco$results.periods[[i]]$DiSCo$quantile)(grid)
      t_cdf <- stats::ecdf(disco$results.periods[[i]]$target$quantile)(grid)
      treats[[i]] <- t_cdf - c_cdf
      if (CI) treats_boot[[i]] <- disco$CI$bootmat$cdf_diff[,i,]  # need to grab from the bootmat
    }

    if (CI){
      agg_nam <- "cdf_diff"
      sds <- lapply(1:T_max, function(t) disco$CI[[agg_nam]]$sd[,t])
      ci_lower <- lapply(1:T_max, function(t) disco$CI[[agg_nam]]$lower[,t])
      ci_upper <- lapply(1:T_max, function(t) disco$CI[[agg_nam]]$upper[,t])
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
      p <- plotDistOverTime(treats, grid, t_start, t_max, CI, ci_lower, ci_upper, savePlots=savePlots,
                       plotName=agg, ylim=ylim, xlim=xlim, xlab="Y", ylab="CDF Change", t_plot=t_plot)
    }
    #---------------------------------------------------------------------------
    ### cdfs
    #---------------------------------------------------------------------------

  } else if (agg == "cdf"){

    treats <- list()
    treats_boot <- list()
    target_cdf <- list()

    # get treatment effects
    for (i in 1:length(disco$results.periods)) {
      grid <- disco$results.periods[[i]]$target$grid
      treats[[i]]  <- stats::ecdf(disco$results.periods[[i]]$DiSCo$quantile)(grid)
      target_cdf[[i]] <- stats::ecdf(disco$results.periods[[i]]$target$quantile)(grid)
      if (CI) treats_boot[[i]] <- disco$CI$bootmat$cdf[,i,] # need to grab from the bootmat
    }

    if (CI){
      sds <- lapply(1:T_max, function(t) disco$CI[[agg]]$sd[,t])
      ci_lower <- lapply(1:T_max, function(t) disco$CI[[agg]]$lower[,t])
      ci_upper <- lapply(1:T_max, function(t) disco$CI[[agg]]$upper[,t])
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA

    }


    if (is.null(xlim)) xlim <- c(min(grid), max(grid))
    if (graph) {
      p <- plotDistOverTime(treats, grid, t_start, t_max, CI, ci_lower, ci_upper, savePlots=savePlots, plotName=agg,
                       obsLine = target_cdf, xlab="Y", ylab="CDF", lty=1, lty_obs=2, xlim=xlim, t_plot=t_plot)
    }
    #---------------------------------------------------------------------------
    ### quantiles of treatment effects
    #---------------------------------------------------------------------------
  } else if (agg == "quantileDiff") { # TODO: adopt for discrete mixture

    treats <- list()

    treats_boot <- list()
    grid <- evgrid


    # get treatment effects
    for (i in 1:length(disco$results.periods)) { # TODO: left off here
      c_qtile <- disco$results.periods[[i]]$DiSCo$quantile
      t_qtile <- disco$results.periods[[i]]$target$quantile
      treats[[i]] <- t_qtile - c_qtile
      if (CI) treats_boot[[i]] <- disco$CI$bootmat$quantile_diff[,i,] # need to grab from the bootmat
    }

    if (CI){
      agg_nam <- "quantile_diff"
      sds <- lapply(1:T_max, function(t) disco$CI[[agg_nam]]$sd[,t])
      ci_lower <- lapply(1:T_max, function(t) disco$CI[[agg_nam]]$lower[,t])
      ci_upper <- lapply(1:T_max, function(t) disco$CI[[agg_nam]]$upper[,t])
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA

    }

    if (graph) {
      if (is.null(ylim)) {
        ymin <- quantile(unlist(treats), 0.01)
        ymax <- quantile(unlist(treats), 0.99)
        ylim <- c(ymin, ymax)
      }
      if (is.null(xlim)) {
        xmin <- min(unlist(evgrid))
        xmax <- max(unlist(evgrid))
        xlim <- c(xmin, xmax)
      }
      p <- plotDistOverTime(treats, grid, t_start, t_max, CI, ci_lower, ci_upper, ylim=ylim,
                       xlab="Quantile", ylab="Treatment Effect", cdf=FALSE, savePlots=savePlots,
                       plotName=agg, t_plot=t_plot)
    }
    #---------------------------------------------------------------------------
    ### counterfactual and observed quantiles
    #---------------------------------------------------------------------------
  } else if (agg == "quantile") {

    treats <- list()
    treats_boot <- list()
    target_qtiles <- list()
    grid <- evgrid


    # get treatment effects
    for (i in 1:length(disco$results.periods)) { # TODO: left off here
      c_qtile <- disco$results.periods[[i]]$DiSCo$quantile
      t_qtile <- disco$results.periods[[i]]$target$quantile
      treats[[i]] <- c_qtile
      target_qtiles[[i]] <- t_qtile
      if (CI) treats_boot[[i]] <- disco$CI$bootmat$quantile[,i,]  # need to grab from the bootmat
    }

    if (CI) {
      sds <- lapply(1:T_max, function(t) disco$CI[[agg]]$sd[,t])
      ci_lower <- lapply(1:T_max, function(t) disco$CI[[agg]]$lower[,t])
      ci_upper <- lapply(1:T_max, function(t) disco$CI[[agg]]$upper[,t])
    } else {
      sds <- NA
      ci_lower <- NA
      ci_upper <- NA
    }

    if (graph) {
      if (is.null(ylim)) {
        ymin <- quantile(unlist(treats), 0.01)
        ymax <- quantile(unlist(treats), 0.99)
        ylim <- c(ymin, ymax)
      }
      if (is.null(xlim)) {
        xmin <- quantile(unlist(grid), 0.01)
        xmax <- quantile(unlist(grid), 0.99)
        xlim <- c(xmin, xmax)
      }
      p <- plotDistOverTime(treats, grid, t_start, t_max, CI, ci_lower, ci_upper, ylim=ylim, xlab="Quantile",
                       ylab="Treatment Effect", cdf=FALSE, obsLine = target_qtiles, savePlots=savePlots, plotName=agg, t_plot=t_plot)
    }

  }
  nams <- names(disco$results.periods)
  if (length(treats) == length(nams)) { # if we have time effects
    names(treats) <- nams
    if (CI) {
      names(sds) <- nams
      names(ci_lower) <- nams
      names(ci_upper) <- nams
    }
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

        if (CI) {
          treats_boot_agg <- apply(treats_boot[[i]][f:t,], 2, mean)
          sd_agg <- sd(treats_boot_agg)
          ci_lower_agg <- stats::quantile(treats_boot_agg, probs=(1-cl)/2)
          ci_upper_agg <- stats::quantile(treats_boot_agg, probs=cl+(1-cl)/2)
        } else {
          sd_agg <- NA
          ci_lower_agg <- NA
          ci_upper_agg <- NA
        }

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
    rownames(out) <- NULL

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
  if (!graph) { p <- NULL }
  return(DiSCoT(agg=agg, treats=treats, grid=grid, ses=sds, ci_lower=ci_lower, ci_upper=ci_upper,
                t0=t0, call=call, cl=cl, N=nrow(disco$params$df), J=data.table::uniqueN(disco$params$df$id_col)-1, agg_df=out,
                perm=disco$perm, plot=p))

}


#' @title Store aggregated treatment effects
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
#' @param perm list of per mutation results
#' @param plot a ggplot object containing the plot for the aggregated treatment effects using the `agg` parameter
#' @return S3 object of class `DiSCoT` with associated `summary` and `print` methods
#' @export
DiSCoT <- function(agg, treats, ses, grid, ci_lower, ci_upper, t0, call, cl, N, J, agg_df, perm, plot) {
  out <- list(agg=agg, treats=treats, ses=ses, ci_lower=ci_lower, ci_upper=ci_upper, t0=t0, call=call, cl=cl,
              grid=grid, N=N, J=J, agg_df=agg_df, perm=perm, plot=plot)
  class(out) <- "DiSCoT"
  out
}

#' @title Plot distribution of treatment effects over time
#' @description Plot distribution of treatment effects over time
#' @param cdf_centered list of centered distributional statistics
#' @param grid_cdf grid
#' @param t_start start time
#' @param t_max maximum time
#' @param CI logical indicating whether to plot confidence intervals
#' @param ci_lower lower confidence interval
#' @param ci_upper upper confidence interval
#' @param ylim y limits
#' @param xlim x limits
#' @param cdf logical indicating whether to plot CDF or quantile difference
#' @param xlab x label
#' @param ylab y label
#' @param obsLine optional additional line to plot. Default is NULL which means no line is plotted.
#' @param savePlots logical indicating whether to save plots
#' @param plotName name of plot to save
#' @param lty line type for the main line passed as cdf_centered
#' @param lty_obs line type for the optional additional line passed as obsLine
#' @param t_plot optional vector of times to plot. Default is NULL which means all times are plotted.
#' @keywords internal
#' @return plot of distribution of treatment effects over time
plotDistOverTime <- function(cdf_centered, grid_cdf, t_start, t_max, CI, ci_lower, ci_upper, ylim=c(0,1), xlim=NULL,
                             cdf=TRUE, xlab="Distribution Difference", ylab="CDF",
                             obsLine = NULL, savePlots=FALSE, plotName=NULL, lty=1, lty_obs=1, t_plot = NULL) {

  # Create a data frame for all data points
  df <- data.frame(time = rep(t_start:t_max, each = length(grid_cdf)),
                   x = rep(grid_cdf, times = length(cdf_centered)),
                   y = unlist(cdf_centered),
                   ci_lower = if (CI) unlist(ci_lower) else NA,
                   ci_upper = if (CI) unlist(ci_upper) else NA,
                   obsLine = if (!is.null(obsLine)) unlist(obsLine) else NA)

  if (!is.null(t_plot)) df <- df[df$time %in% t_plot,]

  # df <- df %>% # filter for points where it changes only
  #   mutate(Change = c(FALSE, diff(y) != 0)) %>%
  #   filter(Change) # TODO: implement to avoid jagged plots

  # Create a ggplot
  p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = y), colour = "black", linetype = lty)

  if (!is.null(obsLine)) p <- p + geom_line(aes(y = obsLine), colour = "blue", linetype = lty_obs, show.legend = !is.null(obsLine))

  p <- p +
    geom_hline(yintercept = 0, colour = "grey") +
    labs(title = "Distribution of Treatment Effects Over Time", x = xlab, y = ylab) +
    theme_minimal() +
    coord_cartesian(xlim = xlim, ylim = ylim)
    # # turn grid off
    # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #       panel.background = element_blank(), axis.line = element_line(colour = "black"))


  if (CI) p <- p + geom_line(aes(y = ci_lower), colour = "grey", linetype = lty, show.legend = CI) +
    geom_line(aes(y = ci_upper), colour = "grey", linetype = lty, show.legend = CI)

  if (cdf) {
    p <- p + geom_hline(yintercept = 1, colour = "grey")
  }

  n_row <- data.table::uniqueN(df$time)
  p <- p + facet_wrap(~time, ncol = 1, nrow = n_row)


  # Save or return the plot
  if (savePlots && !is.null(plotName)) {
    ggsave(paste0(plotName, ".pdf"), plot = p, width = 5 * n_row, height = 5)
  } else {
    print(p)
  }
  return(p)
}





#' @title summary.DiSCoT
#' @description Summary of DiSCoT object
#' @param object DiSCoT object
#' @param ... Additional arguments
#' @return summary of DiSCoT object
#' @export
summary.DiSCoT <- function(object, ...) {

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

  if (!is.null(object$perm))  cat(paste(summary(object$perm), collapse="\n")) else cat("No permutation test performed. \n")

  cat(paste0("Number of pre-treatment periods: ", length(t_list) - length(I)))
  cat("\n")

  cat(paste0("Number of post-treatment periods: ", length(I)))
  cat("\n")


  cat(paste0("N=", formatC(object$N, big.mark=",")))
  cat("\n")
}

