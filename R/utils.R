#' Compute the empirical quantile function
#'
#' @param X A vector containing the data
#' @param q A vector containing the quantiles
#' @return A vector containing the empirical quantile function
#' @keywords internal
myQuant <- function(X,q, qmethod=NULL,...){
  # sort if unsorted
  # if (is.unsorted(X)) X <- sort(X)

  # check if method a valid one
  if (!is.null(qmethod)){
    if (!qmethod %in% c("qkden", "extreme")) stop("Invalid quantile method")
  }
  if (is.null(qmethod)) { # use old-fashioned quantiles
    # obtain the corresponding empirical quantile
    # return(mapply(myquant, q, MoreArgs = list(X=X)))
    return(stats::quantile(X, probs=q))
  } else if (qmethod=="qkden") {
    temp <- evmix::qkden(p=q, kerncentres=X,...)
    temp[1] <- min(c(temp[2]), min(X))
    temp[length(temp)] <- max(c(temp[length(temp)-1]), max(X))
    return(temp)
  } else if (qmethod=="extreme"){
    temp <- extremeStat::distLquantile(X, probs=q, quiet=TRUE)['weighted1',1:length(q)] # take distribution with best fit
    temp[1] <- min(c(temp[2]), min(X))
    temp[length(temp)] <- max(c(temp[length(temp)-1]), max(X))
    return(temp)
  }
}


#' @title getGrid
#' @description Set up a grid for the estimation of the quantile functions and CDFs
#'
#' @param target A vector containing the data for the target unit
#' @param controls A list containing the data for the control units
#' @param G The number of grid points
#' @return A list containing the following elements:
#' \itemize{
#' \item \code{grid.min} The minimum value of the grid
#' \item \code{grid.max} The maximum value of the grid
#' \item \code{grid.rand} A vector containing the grid points
#' \item \code{grid.ord} A vector containing the grid points, ordered
#' }
#' @keywords internal
getGrid <- function(target, controls, G) {
  grid.min <- min(c(min(target),unlist(lapply(controls,min))))
  grid.max <- max(c(max(target),unlist(lapply(controls,max))))

  # if grid.max-grid.min<=1 we round to the next 0.1.
  grid.min <- floor(grid.min*10)/10
  grid.max <- ceiling(grid.max*10)/10

  # sampling uniformly on the grid
  grid.rand <- stats::runif(G,grid.min-0.25,grid.max+0.25)

  # ordered grid
  grid.ord <- grid.rand[order(grid.rand)]

  return(list(grid.min, grid.max, grid.rand, grid.ord))
}

#' @title checks
#' Carry out checks on the inputs
#'
#' @inheritParams DiSCo
#' @param permutation logical, whether to use permutation or not
#' @return NULL
#' @keywords internal
checks <- function(df, id_col.target, t0, M, G, num.cores, permutation, q_min, q_max,
                   CI, CI_placebo, boots, cl, graph,
                   qmethod, seed) {
      # checks on the input data
  if (!id_col.target %in% df$id_col) {
    stop("There is no row in the column `id_col` with the name specified in id_col.target!")
  }
  if (!"time_col" %in% names(df)) {
    stop("time_col is not a column in the data table")
  }
  if (!"y_col" %in% names(df)) {
    stop("y_col is not a column in the data table")
  }
  if (!"id_col" %in% names(df)) {
    stop("id_col is not a column in the data table")
  }

  if ((t0 > max(df$time_col))) {
    stop("First treatment period t0 has to be within time_col range")
  }

  if (t0 <= min(df$time_col)) {
    stop("There are no pre-treatment periods based on provided t0")
  }

  # checks on the input data types
  if (!is.numeric(df$id_col)) {
    stop("id_col must be numeric")
  }
  if (!is.integer(df$time_col)) {
    stop("t_col must be integer")
  }
  if (!is.numeric(df$y_col)) {
    stop("y_col must be numeric")
  }
  if (!is.numeric(t0)) {
    stop("t0 must be numeric")
  }
  if (!is.numeric(M)) {
    stop("M must be numeric")
  }
  if (!is.integer(G)) {
        stop("G must be integer")
  }
  if (!is.integer(num.cores)) {
    stop("num.cores must be integer")
  }
  if (!is.logical(permutation)) {
    stop("permutation must be logical")
  }

  # checks on the input data values
  if ((t0 < min(df$time_col)) | (t0 > max(df$time_col))) {
    stop("T0 must be between minimum and the maximum value of year_col")
  }
  if (M < (uniqueN(df$id_col) -1) ) {
    stop("M must be greater than or equal to the number of control units, and ideally large")
  }
  if (G < 2) {
    stop("G must be greater than or equal to 2, and ideally larger than 100")
  }
  if (num.cores < 1) {
    stop("num.cores must be greater than or equal to 1")
  }

  # check that the number of cores is not greater than the number of available cores
  if (num.cores > parallel::detectCores()) {
    stop("num.cores cannot be greater than the number of available cores")
  }

  if (!is.logical(CI)) {
    stop("CI must be logical")
  }

  if (!is.integer(boots)) {
    stop("boot must be integer")
  }

  if ((CI) & (boots < 2)) {
    stop("IF CI=TRUE, boot must be greater than or equal to 2 (and ideally greater than 100)")
  }

  # check if cl is decimal number between 0 and 1
  if (!is.numeric(cl)) {
    stop("cl must be numeric")
  }
  if ((cl < 0) | (cl > 1)) {
    stop("cl must be between 0 and 1")
  }

  if (!is.logical(graph)) {
    stop("graph must be logical")
  }

  if (!is.numeric(q_min)) {
    stop("q_min must be numeric")
  }

  if (!is.numeric(q_max)) {
    stop("q_max must be numeric")
  }

  if (q_min < 0) {
    stop("q_min must be greater than or equal to 0")
  }

  if (q_max > 1) {
    stop("q_max must be less than or equal to 1")
  }

  if (q_min > q_max) {
    stop("q_min must be less than or equal to q_max")
  }

  if (!is.logical(CI_placebo)) {
    stop("CI_placebo must be logical")
  }

  if (!is.null(qmethod)) {
    if (!qmethod %in% c("qkden", "extreme")) {
      stop("qmethod must be either NULL, 'qkden' or 'extreme'")
    }
  }

  if (!is.null(seed)) {
    if (!is.integer(seed)) {
      stop("seed must be NULL or integer")
    }
  }



}



#' Check if a vector is integer
#'
#' @param x a vector
#' @return TRUE if x is integer, FALSE otherwise
#' @keywords internal
is.integer <- function(x) {
  is.numeric(x) && all(x == as.integer(x))
}




#' @title mclapply.hack
#'
#' @description This function mimics forking (done with mclapply in Mac or Linux) for the
#' Windows environment.  Designed to be used just like mclapply.  Credit goes to
#' Nathan VanHoudnos.
#' @param verbose Should users be warned this is hack-y? Defaults to FALSE.
#' @param mc.cores Number of cores to use. Defaults to 1.
#' @seealso mclapply
#' @keywords internal
#' @importFrom parallel makeCluster detectCores parLapply stopCluster clusterExport
#' @importFrom utils sessionInfo
mclapply.hack <- function(..., verbose=FALSE, mc.cores=1) {

  if (mc.cores == 1) {
    return(lapply(...))
  }

  if( Sys.info()[['sysname']] == 'Windows') {

    ## Create a cluster
    if( is.null(mc.cores) ) {
      size.of.list <- length(list(...)[[1]])
      mc.cores <- min(size.of.list, detectCores())
    }
    ## N.B. setting outfile to blank redirects output to
    ##      the master console, as is the default with
    ##      mclapply() on Linux / Mac
    cl <- makeCluster( mc.cores, outfile="" )

    ## Find out the names of the loaded packages
    loaded.package.names <- c(
      ## Base packages
      utils::sessionInfo()$basePkgs,
      ## Additional packages
      names( sessionInfo()$otherPkgs ))

    tryCatch( {

      ## Copy over all of the objects within scope to
      ## all clusters.
      this.env <- environment()
      while( identical( this.env, globalenv() ) == FALSE ) {
        clusterExport(cl,
                      ls(all.names=TRUE, envir=this.env),
                      envir=this.env)
        this.env <- parent.env(environment())
      }
      clusterExport(cl,
                    ls(all.names=TRUE, envir=globalenv()),
                    envir=globalenv())

      ## Load the libraries on all the clusters
      ## N.B. length(cl) returns the number of clusters
      parLapply( cl, 1:length(cl), function(xx){
        lapply(loaded.package.names, function(yy) {
          require(yy , character.only=TRUE)})
      })

      ## Run the lapply in parallel
      return( parLapply( cl, ...) )
    }, finally = {
      ## Stop the cluster
      stopCluster(cl)
    })

    ## Warn the user if they are using Windows
    if(verbose == TRUE){
      message(paste(
        "\n",
        "   *** Microsoft Windows detected ***\n",
        "   \n",
        "   For technical reasons, the MS Windows version of mclapply()\n",
        "   is implemented as a serial function instead of a parallel\n",
        "   function.",
        "   \n\n",
        "   As a quick hack, we replace this serial version of mclapply()\n",
        "   with a wrapper to parLapply() for this R session. Please see\n\n",
        "     http://www.stat.cmu.edu/~nmv/2014/07/14/implementing-mclapply-on-windows \n\n",
        "   for details.\n\n"))
    }
  } else{
    ## If not on Windows, just call mclapply()
    parallel::mclapply(..., mc.cores=mc.cores)
  }
}


#' @title citation
#'
#' @description print the citation for the relevant paper
#'
#' @keywords internal
citation <- function() {
  cat('Reference: Gunsilius, Florian F. "Distributional synthetic controls." Econometrica 91, no. 3 (2023): 1105-1117. \n')
}


#' @title ex_gmm
#'
#' @description Example data for \code{DiSCo} command.
#' Returns simulated target and control that are mixtures of Gaussian distributions.
#'
#' @param Ts an integer indicating the number of time periods
#' @param num.con an integer indicating the number of control units
#' @param numdraws an integer indicating the number of draws
#' @return
#' \item{\code{target}}{a vector.}
#' \item{\code{control}}{a matrix.}
#' @importFrom MASS mvrnorm
#' @importFrom stats runif
#' @export
ex_gmm <- function(Ts=2, num.con=30, numdraws=1000){

  for (t in 1:Ts) {
    # Mixture of 3 Gaussians
    con <- matrix(0, nrow=num.con, ncol=numdraws)
    for (ii in 1:num.con){
      # generating uniformly distributed weights in the unit simplex
      a1 <- matrix(runif(3), nrow=1)
      a1 <- sweep(a1, 1, rowSums(a1), FUN="/")
      components <- sample(1:3,prob=a1,size=numdraws,replace=TRUE)
      mus <- runif(3,-10,10)
      sigmas <- runif(3,0.5,6)
      sigmas <- (sigmas + t(sigmas))/2
      con[ii,] <- rnorm(numdraws)*sigmas[components]+mus[components]
    }


    # generating the target distribution as a mixture of 4 Normals
    a1 <- matrix(runif(4), nrow=1)
    a1 <- sweep(a1, 1, rowSums(a1), FUN="/")
    components <- sample(1:4,prob=a1,size=numdraws,replace=TRUE)
    mus <- runif(4,-10,10)
    sigmas <- runif(4,0.5,6)
    sigmas <- (sigmas + t(sigmas))/2
    treat <- rnorm(numdraws)*sigmas[components]+mus[components]


    con <- t(con)
    target <- treat

    # stick in dataframe
    df <- data.table(target=target, control=con)
    # convert to long format
    df <- melt(df, id.vars = NULL, measure.vars = c("target", paste0("control.V", 1:num.con)))
    setnames(df, c("variable", "value"), c("id_col", "y_col"))
    df[, time_col := t]

    if (t==1) df_out <- df
    else df_out <- rbind(df_out, df)
  }

  df_out[, id_col := as.numeric(id_col)]

  return(df_out)
}




