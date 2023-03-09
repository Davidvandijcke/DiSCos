
getGrid <- function(target, controls1, G) {
  grid.min <- min(c(min(target),unlist(lapply(controls1,min))))
  grid.max <- max(c(max(target),unlist(lapply(controls1,max))))

  # if grid.max-grid.min<=1 we round to the next 0.1.
  grid.min <- floor(grid.min*10)/10
  grid.max <- ceiling(grid.max*10)/10

  # sampling uniformly on the grid
  grid.rand <- runif(G,grid.min-0.25,grid.max+0.25)

  # ordered grid
  grid.ord <- grid.rand[order(grid.rand)]

  return(list(grid.min, grid.max, grid.rand, grid.ord))
}

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
    stop("year_col must be integer")
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




is.integer <- function(x) {
  is.numeric(x) && all(x == as.integer(x))
}