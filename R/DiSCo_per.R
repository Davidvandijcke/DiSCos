
#' @title DiSCo_per
#'
#' @description Function to implement permutation test for Distributional Synthetic Controls
#' @details This program iterates through all units and computes the optimal weights on the other units
#' for replicating the unit of iteration's outcome variable, assuming that it is the treated unit.
#' See Algorithm 1 in \insertref{gunsilius2020distributional}{DiSCo} for more details.
#'
#' @param c_df List with matrices of control distributions
#' @param t_df Matrix containing the target distribution
#' @param T0 Integer indicating first year of treatment as counted from 1 (e.g, if treatment year 2002 was the 5th year in the sample, this parameter should be 5).
#' @param ww Optional vector of weights indicating the relative importance of each time period. If not specified, each time period is weighted equally.
#' @param peridx Optional integer indicating number of permutations. If not specified, by default equal to the number of units in the sample.
#' @param evgrid Optional vector containing an evenly spaced grid on [0,1] on which the quantile function for the control units will be evaulated.
#' By default, a grid of 100 points is used.
#' @param graph Boolean indicating whether to plot graphs
#' @param y_name Y axis label of the graph
#' @param x_name X axis label of the graph
#' @param num_cores Integer, number of cores to use for parallel computation. Set to 1 by default (sequential computation), this can be very slow!
#' @return List of matrices containing synthetic time path of the outcome variable
#' for the target unit together with the time paths of the control units
DiSCo_per <- function(c_df, t_df, controls.q, T0, ww=0, peridx=0, evgrid=seq(from=0, to=1, length.out=101),
                 graph=TRUE, num_cores = 1, redo_weights=FALSE, weights=NULL){

  #----------------------------------------#
  # target
  #----------------------------------------#

  c_df=controls_per
  t_df=target_per
  controls.q=controls.q
  T0=T0
  weights=Weights_DiSCo_avg
  num_cores=num.cores

  if (redo_weights) {
    #calculate lambda_t for t<=T0
    lambda_t=list()


    lambda_t <- parallel::mclapply.hack(seq_len(T0), function(t) {
      DiSCo_weights_reg(c_df[[t]], as.vector(t_df[[t]]), 1000)
    }, mc.cores = num_cores)

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
    DiSCo_bc(c_df[[x]], controls.q[[x]], lambda.opt, evgrid)
  }, mc.cores = num_cores)

  #computing the target quantile function
  target_q=list()

  for (t in 1:length(t_df)){
    target_q[[t]] <- mapply(myquant, evgrid, MoreArgs = list(X=t_df[[t]]))
  }

  #squared Wasserstein distance between the target and the corresponding barycenter
  distt=c()
  for (t in 1:length(c_df)){
    distt[t]=mean((bc_t[[t]]-target_q[[t]])**2)
  }

  #----------------------------------------#
  # permutation
  #----------------------------------------#


  #default permute all controls
  if (peridx==0){
      peridx=1:length(c_df[[1]])
  }


  cat("Starting permutation test...")
  distp <- mclapply.hack(seq_len(length(peridx)), function(idx) {
    DiSCo_per_iter(c_df=c_df, c_df.q=controls.q, t_df=t_df, T0=T0, ww=ww, peridx=peridx, evgrid=evgrid, idx=idx)
  }, mc.cores = num_cores)
  cat('Permutation finished!')


  # Convert the list to a nested list (the parallelization messes up the output of foreach when choosing .combine = list)
  distp <- matrix(unlist(distp), ncol = length(c_df), byrow = TRUE)
  distp <- split(distp, seq_len(nrow(distp)))


  #default plot all squared Wasserstein distances
  if (graph==TRUE){
    plot(distt, xlab='',ylab='', type='l', lwd=2)
    for (i in 1:length(distp)){
      lines(1:length(c_df), distp[[i]], col='grey', lwd=1)
    }
    abline(v=T0, lty = 2)
    legend("topleft",legend = c("Target", "Control"),
           col=c("black", "grey"),
           lty= c(1,1), lwd = c(2,2), cex = 1.5)
    title(ylab="Squared Wasserstein distance", line=2.5, cex.lab=1.5)
    title(xlab="Time periods", line=3, cex.lab=1.5)
  }


  return(list(target.dist=distt, control.dist=distp))

}



DiSCo_per_rank <- function(distt, distp) {
  #' @param distt List of squared Wasserstein distances between the target unit and the control units
  #' @param distp List of squared Wasserstein distances between the control units
  #' @return List of p-values for each time period
  #' @export

  ## rank the squared Wasserstein distances and get the rank of the target unit
  # combine distt and distp
  distall <- distp
  distall$target <- distt
  distall <- matrix(unlist(distall), nrow=length(distall), ncol=length(distall[[1]]), byrow = TRUE)
  J_1 <- nrow(distall)
  rnks <- list()
  p_values <- list()
  for (t in 1:ncol(distall)){
    # record rank of target unit
    rnk <- rank(-distall[,t])[J_1]
    p_values[[t]] <- rnk / J_1
    rnks[[t]] <- rnk
  }
  p_values <- unlist(p_values)

  return(list(p_t = p_values, ranks=rnks))
}


