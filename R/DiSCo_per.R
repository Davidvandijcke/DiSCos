
#' @title DiSCo_per
#'
#' @description Function to implement permutation test for Distributional Synthetic Controls
#' @details This program iterates through all units and computes the optimal weights on the other units
#' for replicating the unit of iteration's outcome variable, assuming that it is the treated unit.
#' See Algorithm 1 in \insertCite{gunsilius2023distributional;textual}{DiSCo} for more details. The only modification is that we take the ratio of post- and pre-treatment
#' root mean squared Wasserstein distances to calculate the p-value, rather than the level in each period, following @abadie2010synthetic.
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
#' @references
#' \insertAllCited{}
DiSCo_per <- function(c_df, t_df, controls.q, target.q, T0, ww=0, peridx=0, evgrid=seq(from=0, to=1, length.out=101),
                 graph=TRUE, num_cores = 1, redo_weights=FALSE, weights=NULL, qmethod=NULL, per_q_min=0, per_q_max=1, M){

  #----------------------------------------#
  # target
  #----------------------------------------#


  if (redo_weights) {
    #calculate lambda_t for t<=T0
    lambda_t=list()


    lambda_t <- parallel::mclapply.hack(seq_len(T0), function(t) {
      DiSCo_weights_reg(c_df[[t]], as.vector(t_df[[t]]), M, qmethod=qmethod)
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



  #squared Wasserstein distance between the target and the corresponding barycenter
  # find the closest quantile range that includes [per_q_min, per_q_max]
  per_q_floor <- floor(per_q_min * length(evgrid))
  per_q_ceil <- ceiling(per_q_max * length(evgrid))
  distt=c()
  for (t in 1:length(c_df)){
    distt[t]=mean((bc_t[[t]][per_q_floor:per_q_ceil]-target.q[[t]][per_q_floor:per_q_ceil])**2)
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
    DiSCo_per_iter(c_df=c_df, c_df.q=controls.q, t_df=t_df, T0=T0, ww=ww, peridx=peridx, evgrid=evgrid, idx=idx, qmethod=qmethod, M=M,
                   per_q_min=per_q_min, per_q_max=per_q_max)
  }, mc.cores = num_cores)
  cat('Permutation finished!')


  # Convert the list to a nested list (the parallelization messes up the output of foreach when choosing .combine = list)
  distp <- matrix(unlist(distp), ncol = length(c_df), byrow = TRUE)
  distp <- split(distp, seq_len(nrow(distp)))


  #default plot all squared Wasserstein distances
  if (graph==TRUE){
    plot(distt, xlab='',ylab='', type='l', lwd=2, ylim=c(0, max(c(distt, unlist(distp))) + 0.1*max(c(distt, unlist(distp)))))
    for (i in 1:length(distp)){
      lines(1:length(c_df), distp[[i]], col='grey', lwd=1)
    }
    abline(v=T0, lty = 2)
    title(ylab="Squared Wasserstein distance", line=2.5, cex.lab=1.5)
    title(xlab="Time periods", line=3, cex.lab=1.5)


  }


  return(list(target.dist=distt, control.dist=distp))

}


#' @title DiSCo_per_rank
#' @description This function ranks the squared Wasserstein distances and returns the p-values for each time period
#' @param distt List of squared Wasserstein distances between the target unit and the control units
#' @param distp List of squared Wasserstein distances between the control units
#' @return List of p-values for each time period
#' @export
## rank the squared Wasserstein distances and get the rank of the target unit
# combine distt and distp
DiSCo_per_rank <- function(distt, distp, T0) {
  distall <- distp
  distall$target <- distt
  distall <- matrix(unlist(distall), nrow=length(distall), ncol=length(distall[[1]]), byrow = TRUE)
  J_1 <- nrow(distall)
  rnks <- list()
  p_values <- list()

  R <- apply(distall, 1, function(x) sqrt(mean(x[1:T0])) / sqrt(mean(x[(T0+1):length(x)])) )
  p_val <- rank(-R)[1] / (J_1) # minus cause we want the largest to be the smallest rank


  return(p_val)
}


#' @title permut
#' @description Object to hold results of permutation test
#'
#' @param distp List of squared Wasserstein distances between the control units
#' @param distt List of squared Wasserstein distances between the target unit and the control units
#' @param rank Rank of the target unit
#' @param p_values List of p-values for each time period
#' @param p_overall Overall p-value
#' @param J_1 Number of control units
#' @return A list of class permut, with the following elements
#' \item{distp}{List of squared Wasserstein distances between the control units}
#' \item{distt}{List of squared Wasserstein distances between the target unit and the control units}
#' \item{p_overall}{Overall p-value for post-treatment periods}
#' \item{J_1}{Number of control units}
#'
#' @export
#' @examples
#' permut(distp, distt, rank, p_values, p_overall, J_1)
#' @export
permut <- function(distp, distt,p_overall, J_1) {
  out <- list(distp=distp, distt=distt, p_overall=p_overall, J_1=J_1)
  class(out) <- "permut"
  return(out)
}


#' @ title print.permut
#'
#' @description Print permutation test results
#'
#' @param x Object of class permut
#' @param digits Number of digits to print
#' @return Prints permutation test results
#' @export
#' @examples
#' print(x, digits=3)
#' @export
print.permut <- function(x, digits = 3) {
  cat("Permutation test results:\n")
  cat("P-value: ", format(x$p_overall, digits=digits), "\n")
  cat("Number of control units: ", x$J_1, "\n")
}

#' @ title summary.permut
#' @description Summarize permutation test results
#' @param object Object of class permut
#' @param digits Number of digits to print
#' @return Prints permutation test results
#' @export
#' @examples
#' summary(x, digits=3)
summary.permut <- function(x, digits = 3) {
  print.permut(x, digits=digits)
}

