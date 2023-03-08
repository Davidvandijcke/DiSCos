
fuck = "state_fips"
y_col = "adj0contpov"
target_id = fips.target

df$id_col = df$state_fips
df$t_col = df$year
df$y_col = df$adj0contpov





getMixture <- function(controls1, target, grid.min, grid.max, grid.rand) {

  ###### The mixture of distributions approach
  # we again only focus on the first half of the data
  # defining the grid on which we define the cumulative distribution functions
  # obtaining the minimal and maximal values among all supports
  # creating a list of controls with only the full data


  # Estimating the empirical CDFs
  CDF.control <- lapply(controls1,ecdf)
  CDF.target <- ecdf(target)


  # Evaluating the CDF on the random grid
  CDF.matrix <- matrix(0,nrow=length(grid.rand), ncol = (length(controls1)+1))
  CDF.matrix[,1] <- CDF.target(grid.rand)
  for (ii in 1:length(controls1)){
    CDF.matrix[,(ii+1)] <- CDF.control[[ii]](grid.rand)
  }

  # Solving the convex problem with CVXR
  # the variable we are after
  theweights <- Variable(length(controls1))
  # the objective function
  objective <- cvxr_norm((CDF.matrix[,2:ncol(CDF.matrix)] %*% theweights - CDF.matrix[,1]))

  # the constraints for the unit simplex
  constraints <- list(theweights>=0, sum_entries(theweights) == 1)
  # the optimization problem
  problem <- Problem(Minimize(objective),constraints)
  # solving the optimization problem
  results <- solve(problem, solver = "SCS")

  # returning the optimal weights and the value function which provides the
  # squared Wasserstein distance between the target and the corresponding barycenter
  theweights.opt <- results$getValue(theweights)
  thedistance.opt <- results$value*1/M*(grid.max - grid.min)


  themean <- CDF.matrix[,2:ncol(CDF.matrix)]%*%theweights.opt


  themean.order <- themean[order(grid.rand, decreasing=FALSE)]
  target.order <- CDF.matrix[order(grid.rand, decreasing=FALSE),1]


  return(list("weights.opt" = theweights.opt, "distance.opt" = thedistance.opt, 
              "mean" = themean.order, "target.order" = target.order, "cdf" = CDF.matrix))

}

getGrid <- function(target, controls1) {
  grid.min <- min(c(min(target),unlist(lapply(controls1,min))))
  grid.max <- max(c(max(target),unlist(lapply(controls1,max))))

  # if grid.max-grid.min<=1 we round to the next 0.1.
  grid.min <- floor(grid.min*10)/10
  grid.max <- ceiling(grid.max*10)/10

  # sampling uniformly on the grid
  grid.rand <- runif(M,grid.min-0.25,grid.max+0.25)

  # ordered grid
  grid.ord <- grid.rand[order(grid.rand)]

  return(list(grid.min, grid.max, grid.rand, grid.ord))
}

DiSCo_time <- function(yy, ...) {

    # obtaining the target state

    # target outcome
    target <- df[id_col == target][[y_col]]

    
    # generate list where each element contains a list of all micro-level outcomes for a control unit
    controls <- list()
    j <- 1
    controls.id <- unique(df[id_col != target_id, id_col])
    for (id in controls.id) {
      controls[j] <- list(df[id_col == id & t_col == yy, y_col])
      j <- j + 1
    }

    # obtaining the optimal weights for the DSC method
    DSC_res_weights <- DSC_weights_reg(controls,as.vector(target), 1000)

    DSC_res2 <- DSC_bc(controls,DSC_res_weights,seq(from=0,to=1,length.out=1001))

    # sample grid
    grid <- list(grid.min = NA, grid.max = NA, grid.rand = NA, grid.ord = NA)
    grid[c("grid.min", "grid.max", "grid.rand", "grid.ord")] <- getGrid(target, controls)

    # getting the CDF from the quantile function
    DSC_res2.cdfF <- ecdf(DSC_res2[[2]])
    DSC_res2.cdf <- DSC_res2.cdfF(grid$grid.ord)

    # assign DSC results to named list
    #DSC_res[(c)
    
 
    # obtaining the optimal weights for the mixture of distributions method
    mixture <- getMixture(controls, target, grid$grid.min, grid$grid.max, grid$grid.rand)

    y_char <- as.character(yy)
    results <- list()
    results[["DSC"]] <- 
      list("weights" = DSC_res_weights, "quantile" = DSC_res2, "cdf" =  DSC_res2.cdf) # DSC estimator
    results[["mixture"]] <- list("weights" = mixture$weights.opt, "distance" = mixture$distance.opt, "mean" = mixture$mean) # mixture of distributions estimator
    results[["target"]] <- list("quantile" = target.s, "cdf" = mixture$target.order, "grid" =  grid.ord, "data" = as.vector(target))
    results[["controls"]] <- list("data" = controls, "cdf" = mixture$CDF.matrix, "id" = controls.id)

}

DiSCo <- function() {



  #####
  # Solving for the optimal weights in the DSC method and the alternative method using mixtures of
  # distributions by year

  # create named list for storing results, one element for each year
  results.by.year <- list()

  # looping over the years from 1998 - 2014



}
