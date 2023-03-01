# Function for computing barycenters in the DSC method at every time period

DSC_bc <- function(controls,weights,evgrid = seq(from=0, to=1, length.out=101)){
  # this program takes in samples from the distributions of the controls and puts out
  # 1. the quantile function of the barycenter associated with the "weights" evaluated at the vector
  #    "evgrid"
  # 2. the quantile functions of the control groups evaluated at the vector "evgrid".
  
  # the control distributions can be given in list form, where each list element contains a
  # vector of observations for the given control unit, in matrix form;
  # in matrix- each column corresponds to one unit and each row is one observation. 
  # The list-form is useful, because the number of draws for each control group can be different. 
  if (!is.list(controls) && !is.matrix(controls)){
    stop ("Controls need to be given in either list-form or matrix form.")
  } 
  
  # if the controls are given in matrix form we turn them into a list
  if (is.matrix(controls)) {
    controls.h <- list()
    for (ii in 1:ncol(controls)){
      controls.h[[ii]] <- as.vector(controls[,ii])
    }
    controls <- controls.h
  }
  
  
  ## Creating a function for the empirical quantile function
  myquant <- function(X,q){
    # sort if unsorted
    if (is.unsorted(X)) X <- sort(X)
    # compute empirical CDF
    X.cdf <- 1:length(X) / length(X)
    # obtain the corresponding empirical quantile
    return(X[which(X.cdf >= q)[1]])
  }
  
  # evaluating the quantile functions on the grid "evgrid":
  controls.q <- matrix(0,nrow = length(evgrid), ncol=length(controls))
  for (jj in 1:length(controls)){
    controls.q[,jj] <- mapply(myquant, evgrid, MoreArgs = list(X=controls[[jj]]))
  }
  
  # Obtaining the Wasserstein barycenter as the average of the quantile functions 
  # weighted by "weights" and evaluate it on the grid "evgrid"
  thebc <- controls.q%*%weights
  
  # returning the quantile functions and the barycenter
  return(list(controls.q,thebc))
}