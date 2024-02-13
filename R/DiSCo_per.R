
#' @title DiSCo_per
#'
#' @description Function to implement permutation test for Distributional Synthetic Controls
#' @details This program iterates through all units and computes the optimal weights on the other units
#' for replicating the unit of iteration's outcome variable, assuming that it is the treated unit.
#' See Algorithm 1 in \insertCite{gunsilius2023distributional;textual}{DiSCos} for more details. The only modification is that we take the ratio of post- and pre-treatment
#' root mean squared Wasserstein distances to calculate the p-value, rather than the level in each period, following @abadie2010synthetic.
#'
#' @param results.periods List of period-specific results from DiSCo
#' @param T0 Integer indicating first year of treatment as counted from 1 (e.g, if treatment year 2002 was the 5th year in the sample, this parameter should be 5).
#' @param ww Optional vector of weights indicating the relative importance of each time period. If not specified, each time period is weighted equally.
#' @param peridx Optional integer indicating number of permutations. If not specified, by default equal to the number of units in the sample.
#' @param redo_weights Boolean indicating whether to recompute the weights for the "true" treated unit. Set to FALSE by default.
#' @param weights Optional vector of weights to use for the "true" treated unit. `redo_weights` has to be set to FALSE for these weights to be used.
#' @inheritParams DiSCo
#' @return List of matrices containing synthetic time path of the outcome variable
#' for the target unit together with the time paths of the control units
#' @references
#' \insertAllCited{}
#' @keywords internal
#' # TODO add option to select the post-treatment time periods
DiSCo_per <- function(results.periods, T0, ww=0, peridx=0, evgrid=seq(from=0, to=1, length.out=101),
                 graph=TRUE, num.cores = 1, redo_weights=FALSE, weights=NULL, qmethod=NULL, q_min=0, q_max=1, M=1000, simplex=FALSE){

  # slightly hacky way to fix reporting of q_min and q_max
  q_min_report <- q_min
  q_max_report <- q_max

  q_min <- 0
  q_max <- 1

  ## prep
  # grab the raw control and target data
  c_df <- lapply(seq(1:length(results.periods)), function(x) results.periods[[x]]$controls$data)
  t_df <- lapply(seq(1:length(results.periods)), function(x) results.periods[[x]]$target$data)


  # grab the quantiles for the control and target data, for the desired range
  controls.q <- lapply(seq(1:length(results.periods)), function(x) results.periods[[x]]$controls$quantiles) # here it's a matrix hence ,
  target.q <- lapply(seq(1:length(results.periods)), function(x) results.periods[[x]]$target$quantiles) # here it's a vector

  #----------------------------------------#
  # target
  #----------------------------------------#


  if (redo_weights) {
    #calculate lambda_t for t<=T0
    lambda_t=list()


    lambda_t <- mclapply.hack(seq_len(T0), function(t) {
      DiSCo_weights_reg(c_df[[t]], as.vector(t_df[[t]]), M, qmethod=qmethod, q_min=q_min, q_max=q_max, simplex=simplex)
    }, mc.cores = num.cores)

    #calculate the average optimal lambda
    if (length(ww)==1){
      w_t=rep(1/T0, T0)
      lambda.opt=matrix(unlist(lambda_t),ncol=T0)%*%w_t
    } else{
      lambda.opt=matrix(unlist(lambda_t),ncol=T0)%*%ww
    }
  } else if (is.null(weights)){
    stop("Please provide either weights or set redo_weights to TRUE")
  } else {
    lambda.opt=weights
  }

  #calculate the barycenters for each period
  bc_t=list()


  bc_t <- mclapply.hack(1:length(c_df), function(x) {
    DiSCo_bc(controls.q[[x]], lambda.opt, evgrid)
  }, mc.cores = num.cores)




  distt=c()
  for (t in 1:length(c_df)){
    distt[t]=mean((bc_t[[t]]-target.q[[t]])**2)
  }

  #----------------------------------------#
  # permutation
  #----------------------------------------#


  #default permute all controls
  if (peridx==0){
      peridx=1:length(c_df[[1]])
  }


  distp <- mclapply.hack(seq_len(length(peridx)), function(idx) {
    DiSCo_per_iter(c_df=c_df, c_df.q=controls.q, t_df=t_df, T0=T0, ww=ww, peridx=peridx, evgrid=evgrid, idx=idx, qmethod=qmethod, M=M,
                   q_min=q_min, q_max=q_max, simplex=simplex)
  }, mc.cores = num.cores)


  # Convert the list to a nested list (the parallelization messes up the output of foreach when choosing .combine = list)
  distp <- matrix(unlist(distp), ncol = length(c_df), byrow = TRUE)
  distp <- split(distp, seq_len(nrow(distp)))


  #default plot all squared Wasserstein distances
  if (graph==TRUE){
    # Prepare the data for ggplot
    df_distt <- data.frame(x = 1:length(distt), y = distt)
    df_distp <- do.call(rbind, lapply(1:length(distp), function(i) data.frame(x = 1:length(c_df), y = distp[[i]], group = i)))

    # Create the base plot
    p <- ggplot() +
      geom_line(data = df_distp, aes(x = x, y = y, group = group), color = "grey", size = 0.5) +
      geom_line(data = df_distt, aes(x = x, y = y), size = 1) +
      geom_vline(xintercept = T0, linetype = "dashed") +
      labs(x = "Time periods", y = "Squared Wasserstein distance") +
      scale_x_continuous(breaks = seq(1, length(distt), 1)) + # single tick for each time period
      ylim(0, max(c(distt, unlist(distp))) + 0.1*max(c(distt, unlist(distp)))) +
      theme_minimal()

    # Print the plot
    print(p)


  } else {
    p <- NULL
  }


  # calculate the ranks and p-values
  p_val <- DiSCo_per_rank(distt, distp, T0)

  # create the permutation object for easy summarization
  J_1 <- length(distp)
  perm_obj <- permut(distp, distt, p_val, J_1, q_min=q_min_report, q_max=q_max_report, plot=p)


  return(perm_obj)

}



#' @title DiSCo_per_rank
#' @description This function ranks the squared Wasserstein distances and returns the p-values for each time period
#' @param distt List of squared Wasserstein distances between the target unit and the control units
#' @param distp List of squared Wasserstein distances between the control units
#' @return List of p-values for each time period
#' @keywords internal
## rank the squared Wasserstein distances and get the rank of the target unit
# combine distt and distp
DiSCo_per_rank <- function(distt, distp, T0) {
  distall <- distp
  distall$target <- distt
  distall <- matrix(unlist(distall), nrow=length(distall), ncol=length(distall[[1]]), byrow = TRUE)
  J_1 <- nrow(distall)
  rnks <- list()
  p_values <- list()

  R <- apply(distall, 1, function(x) sqrt(mean(x[(T0+1):length(x)])) / sqrt(mean(x[1:T0]))  )
  p_val <- (rank(-R)[length(R)]) / (J_1) # minus cause we want to count the number of tests it is smaller than


  return(p_val)
}


#' @title permut
#' @description Object to hold results of permutation test
#'
#' @param distp List of squared Wasserstein distances between the control units
#' @param distt List of squared Wasserstein distances between the target unit and the control units
#' @param p_overall Overall p-value
#' @param J_1 Number of control units
#' @param q_min Minimum quantile
#' @param q_max Maximum quantile
#' @param plot ggplot object containing plot of squared Wasserstein distances over time for all permutations.
#' @return A list of class permut, with the same elements as the input arguments.
#'
#' @export
permut <- function(distp, distt,p_overall, J_1, q_min, q_max, plot) {
  out <- list(distp=distp, distt=distt, p_overall=p_overall, J_1=J_1, q_min=q_min, q_max=q_max, plot=plot)
  class(out) <- "permut"
  return(out)
}


#' @title print.permut
#'
#' @description Print permutation test results
#'
#' @param x Object of class permut
#' @param ... Additional arguments
#' @return Prints permutation test results
#' @export
#' @keywords internal
print.permut <- function(x, ...) {
  cat(paste0("Permutation test for quantile range: [", x$q_min, ", ", x$q_max, "] \n"))
  cat("P-value: ", format(x$p_overall, digits=3), "\n")
  cat("Number of control units: ", x$J_1, "\n")
}

#' @title summary.permut
#' @description Summarize permutation test results
#' @param object Object of class permut
#' @param ... Additional arguments
#' @return Prints permutation test results
#' @export
summary.permut <- function(object, ...) {
  print.permut(object, digits=3)
}

