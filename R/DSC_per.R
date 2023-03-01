# function to implement the permutation test

DSC_per=function(c_df, t_df, T0, ww=0, peridx=0, evgrid=seq(from=0, to=1, length.out=101),
                 graph=TRUE, y_name='y', x_name='x'){
  
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
  
  for (t in 1:T0){
    lambda_t[[t]] <- DSC_weights_reg(c_df[[t]],as.vector(t_df[[t]]), 1000)
  }
  
  #calculate the average optimal lambda
  if (length(ww)==1){
    w_t=rep(1/T0, T0)
    lambda.opt=matrix(unlist(lambda_t),ncol=T0)%*%w_t
  }else{
    lambda.opt=matrix(unlist(lambda_t),ncol=T0)%*%ww
  }
  
  
  #calculate the barycenters for each period
  bc_t=list()
  
  for (t in 1:length(c_df)){
    bc_t[[t]] <- DSC_bc(c_df[[t]],lambda.opt,evgrid)
  }
  
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
  
  #initiate progress bar
  cat('Permutation starts')
  pb <- txtProgressBar(min=0, max=100, style=3) #initiate progress bar
  
  
  #default permute all controls
  if (length(peridx)==1){
    if (peridx==0){
      peridx=1:length(c_df[[1]])
    }
  }
  
  
  for (idx in 1:length(peridx)){
    
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
    
    distp[[idx]]=dist
    
    
    Sys.sleep(0.025)
    setTxtProgressBar(pb, round(idx / length(peridx) * 100)) # progress bar
    cat(if (idx == length(peridx)) '\n')
  }
  
  
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