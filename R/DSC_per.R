
#' Function to implement permutation test for Distributional Synthetic Controls

#' This program iterates through all units and computes the optimal weights on the other units
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
DSC_per <- function(c_df, t_df, T0, ww=0, peridx=0, evgrid=seq(from=0, to=1, length.out=101),
                 graph=TRUE, y_name='y', x_name='x', num_cores = 1){

  #Creating a function for the empirical quantile function
  myquant <- function(X,q){
    # sort if unsorted
    if (is.unsorted(X)) X <- sort(X)
    # compute empirical CDF
    X.cdf <- 1:length(X) / length(X)
    # obtain the corresponding empirical quantile
    return(X[which(X.cdf >= q)[1]])
  }

  #----------------------------------------#
  # target
  #----------------------------------------#

  #calculate lambda_t for t<=T0
  lambda_t=list()


  lambda_t <- parallel::mclapply(seq_len(T0), function(t) {
    DSC_weights_reg(c_df[[t]], as.vector(t_df[[t]]), 1000)
  }, mc.cores = num_cores)

  #calculate the average optimal lambda
  if (length(ww)==1){
    w_t=rep(1/T0, T0)
    lambda.opt=matrix(unlist(lambda_t),ncol=T0)%*%w_t
  } else{
    lambda.opt=matrix(unlist(lambda_t),ncol=T0)%*%ww
  }


  #calculate the barycenters for each period
  bc_t=list()


  bc_t <- parallel::mclapply(c_df, function(x) {
    DSC_bc(x, lambda.opt, evgrid)
  }, mc.cores = num_cores)

  #computing the target quantile function
  target_q=list()

  for (t in 1:length(t_df)){
    target_q[[t]] <- mapply(myquant, evgrid, MoreArgs = list(X=t_df[[t]]))
  }

  #squared Wasserstein distance between the target and the corresponding barycenter
  distt=c()
  for (t in 1:length(c_df)){
    distt[t]=mean((bc_t[[t]][[2]]-target_q[[t]])**2)
  }

  #----------------------------------------#
  # permutation
  #----------------------------------------#

  #list for squared Wasserstein distance
  distp=list()

  cat('Permutation starts')

  #default permute all controls
  if (length(peridx)==1){
    if (peridx==0){
      peridx=1:length(c_df[[1]])
    }
  }



  if (num.cores == 1) { # sequential computation

    distp <- list()
    
    pb <- txtProgressBar(min = 1, max = length(peridx), style = 3, title = "Permutation progress:")

    # START for loop
    for (idx in 1:length(peridx)){
      # one iterations of the permutation test
      distp[[idx]] <- DSC_per_iter(c_df=c_df, t_df=t_df, T0=T0, ww=ww, peridx=peridx, evgrid=evgrid, idx=idx) # the arguments are the same as the function arguments
      
      # update progress bar
      setTxtProgressBar(pb, idx)
    }
    # END for loop

    # close progress bar
    close(pb)

  } else if (num.cores > 1) { # parallel computation: no progress bar since dysfunctional in parallelization

    # register clusters
    cl <- parallel::makeCluster(num_cores)
    doSNOW::registerDoSNOW(cl) # compatible with both windows and Unix platforms

    # START for loop
    distp <- foreach(idx = 1:length(peridx), .combine = c, .export=c("DSC_weights_reg", "DSC_bc", "myquant")) %dopar% {

      # one iterations of the permutation test
      return(DSC_per_iter(c_df=c_df, t_df=t_df, T0=T0, ww=ww, peridx=peridx, evgrid=evgrid, idx=idx)) # the arguments are the same as the function arguments

    }
    # END for loop

    # stop cluster
    snow::stopCluster(cl)
  }
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
    legend("topleft",legend = c("Target", "Control"),
           col=c("black", "grey"),
           lty= c(1,1), lwd = c(2,2), cex = 1.5)
    title(ylab=y_name, line=2, cex.lab=1.5)
    title(xlab=x_name, line=2, cex.lab=1.5)
  }


  return(list(target.dist=distt, control.dist=distp))

}



DSC_per_iter <- function(c_df, t_df, T0, ww, peridx, evgrid, idx){
    # One iteration of the permutation test

    #create new control and target
    pert=list()
    perc=list()
    for (i in 1:length(c_df)){
      perc[[i]]=list()
    }

    for (i in 1:length(perc)){
      perc[[i]][[1]]=t_df[[i]]
    }

    keepcon=peridx[-idx]

    for (i in 1:length(perc)){
      for (j in 1:length(keepcon)){
        perc[[i]][[j+1]]=c_df[[i]][[keepcon[j]]]
      }
    }

    for (i in 1:length(c_df)){
      pert[[i]]=c_df[[i]][[idx]]
    }


    #calculate lambda_t for t<=T0
    lambda_tp=list()

    for (t in 1:T0){
      lambda_tp[[t]] <- DSC_weights_reg(perc[[t]],as.vector(pert[[t]]), 1000)
    }


    #calculate the average optimal lambda
    if (length(ww)==1){
      w_t=rep(1/T0, T0)
      lambda.opt=matrix(unlist(lambda_tp),ncol=T0)%*%w_t
    }else{
      lambda.opt=matrix(unlist(lambda_tp),ncol=T0)%*%ww
    }


    #calculate the barycenters for each period
    bc_t=list()

    for (t in 1:length(perc)){
      bc_t[[t]]=DSC_bc(perc[[t]],lambda.opt,evgrid)
    }


    # computing the target quantile function
    target_q=list()

    for (t in 1:length(pert)){
      target_q[[t]] <- mapply(myquant, evgrid, MoreArgs = list(X=pert[[t]]))
    }


    #squared Wasserstein distance between the target and the corresponding barycenter
    dist=c()
    for (t in 1:length(perc)){
      dist[t]=mean((bc_t[[t]][[2]]-target_q[[t]])**2)
    }
    #setTxtProgressBar(pb, i)

    return(dist)
}