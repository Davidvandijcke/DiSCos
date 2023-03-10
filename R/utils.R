#' Compute the empirical quantile function
#' 
#' @param X A vector containing the data
#' @param q A vector containing the quantiles
#' @return A vector containing the empirical quantile function
#' @export
#' @examples
#' set.seed(123)
#' X <- rnorm(100)
#' q <- 0.1
#' myquant(X,q)
myquant <- function(X,q){
  # sort if unsorted
  if (is.unsorted(X)) X <- sort(X)
  # compute empirical CDF
  X.cdf <- 1:length(X) / length(X)
  # obtain the corresponding empirical quantile
  return(X[which(X.cdf >= q)[1]])
}


#' Set up a grid for the estimation of the quantile functions and CDFs
#' 
#' @inheritParams DiSco_iter
#' @return A list containing the following elements:
#' \itemize{
#' \item{grid.min}{The minimum value of the grid}
#' \item{grid.max}{The maximum value of the grid}
#' \item{grid.rand}{A vector containing the grid points}
#' \item{grid.ord}{A vector containing the grid points, ordered}
#' }
#' @export
#' @examples
#' set.seed(123)
#' target <- rnorm(100)
#' controls <- list(rnorm(100),rnorm(100),rnorm(100))
#' grid <- list(grid.min = NA, grid.max = NA, grid.rand = NA, grid.ord = NA)
#' grid[c("grid.min", "grid.max", "grid.rand", "grid.ord")] <- getGrid(target, controls, G)
#' plot(grid$grid.ord, type="l")
getGrid <- function(target, controls, G) {
  grid.min <- min(c(min(target),unlist(lapply(controls,min))))
  grid.max <- max(c(max(target),unlist(lapply(controls,max))))

  # if grid.max-grid.min<=1 we round to the next 0.1.
  grid.min <- floor(grid.min*10)/10
  grid.max <- ceiling(grid.max*10)/10

  # sampling uniformly on the grid
  grid.rand <- runif(G,grid.min-0.25,grid.max+0.25)

  # ordered grid
  grid.ord <- grid.rand[order(grid.rand)]

  return(list(grid.min, grid.max, grid.rand, grid.ord))
}

#' Carry out checks on the inputs
#'
#' @inheritParams DiSCo
#' @param permutation logical, whether to use permutation or not
#' @return NULL
#' @export
#' @examples
#' df <- data.frame(id_col = c(1,1,1,2,2,2,3,3,3), t_col = c(1,2,3,1,2,3,1,2,3), y_col = c(1,1,0,1,1,0,1,1,0))
#' id_col.target <- 1
#' T0 <- 1  
#' M <- 1000  
#' G <- 1000
#' num.cores <- 1
#' permutation <- FALSE
#' checks(df, id_col.target, T0, M, G, num.cores, permutation)
checks <- function(df, id_col.target, T0, M, G, num.cores, permutation) {
      # checks on the input data
  if (!id_col.target %in% df$id_col) {
    stop("There is no row in the column `id_col` with the name specified in id_col.target!")
  }
  if (!"t_col" %in% names(df)) {
    stop("t_col is not a column in the data table")
  }
  if (!"y_col" %in% names(df)) {
    stop("y_col is not a column in the data table")
  }
  if (!"id_col" %in% names(df)) {
    stop("id_col is not a column in the data table")
  }
  
  # checks on the input data types
  if (!is.numeric(df$id_col)) {
    stop("id_col must be numeric")
  }
  if (!is.integer(df$t_col)) {
    stop("t_col must be integer")
  }
  if (!is.numeric(df$y_col)) {
    stop("y_col must be numeric")
  }
  if (!is.numeric(T0)) {
    stop("T0 must be integer")
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
  if (T0 < 1 | T0 > max(df$t_col)) {
    stop("T0 must be between 1 and the maximum value of year_col")
  }
  if (M < 1) {
    stop("M must be greater than or equal to 1")
  }
  if (num.cores < 1) {
    stop("num.cores must be greater than or equal to 1")
  }
  
  # check that the number of cores is not greater than the number of available cores
  if (num.cores > parallel::detectCores()) {
    stop("num.cores cannot be greater than the number of available cores")
  }

  # check that the year_col is a sequence of integers starting at 1
  if (!all(sort(unique(df$t_col)) == 1:max(df$t_col))) {
    stop("t_col must be a sequence of integers starting at 1")
  }

}



#' Check if a vector is integer
#' 
#' @param x a vector
#' @return TRUE if x is integer, FALSE otherwise
#' @export
is.integer <- function(x) {
  is.numeric(x) && all(x == as.integer(x))
}
