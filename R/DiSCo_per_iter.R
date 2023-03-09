
DiSCo_per_iter <- function(c_df, t_df, T0, ww, peridx, evgrid, idx){
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
      lambda_tp[[t]] <- DiSCo_weights_reg(perc[[t]],as.vector(pert[[t]]), 1000)
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
      bc_t[[t]]=DiSCo_bc(perc[[t]],lambda.opt,evgrid)
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