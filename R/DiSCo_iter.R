
#' @title DiSCo_iter
#'
#' @description This function implements the DiSCo method for a single time period, as well as the mixture of distributions approach.
#'
#' @details This function is called for every time period in the DiSCo function. It implements the DiSCo method for a single time period, as well as the mixture of distributions approach.
#' The corresponding results for each year can be accessed in the `years` list of the output of the DiSCo function. The DiSCo function returns the average weight for each unit across all years,
#' calculated as a uniform mean, as well as the counterfactual target distribution produced as the weighted average of the control distributions for each year, using these averaged weights.
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
#' @return A nested list with the following elements:
#' \itemize{
#' \item{DiSCo}{A list containing the results of the main (DiSCo) method, with the following elements:
#' \itemize{
#' \item{weights}{The optimal weights for the DiSCo method}
#' \item{quantile}{A list containing the quantile functions of the controls and the corresponding barycenter, with the following elements:}
#' \itemize{
#' \item{controls}{A vector containing the quantile functions of the controls}
#' \item{barycenter}{A vector containing the quantile function of the barycenter}
#' }
#' \item{cdf}{The empirical CDFs of the barycenters}
#' }
#' \item{mixture}{A list containing the results of the mixture of distributions approach, with the following elements:
#' \itemize{
#' \item{weights}{The optimal weights for the mixture of distributions approach}
#' \item{distance}{The value of the objective function for the mixture of distributions approach}
#' \item{mean}{The weighted mixture of the controls' CDFs, i.e. the "mixture CDF"}
#' }
#' \item{target}{A list containing the data for the target unit, with the following elements:}
#' \itemize{
#' \item{quantile}{A vector containing the quantile functions of the target}
#' \item{cdf}{A vector containing the empirical CDFs of the target}
#' \item{grid}{A vector containing the grid on which the quantile and CDF functions were evaluated}
#' \item{data}{A vector containing the supplied data for the target unit}
#' }
#' \item{controls}{A list containing the data for the control units, with the following elements:}
#' \itemize{
#' \item{data}{A vector containing the supplied data for the control units}
#' \item{cdf}{A vector containing the empirical CDFs of the controls}
#' \item{id}{A vector containing the IDs of the control units, in the same ordering as the weights returned in the DiSCo and mixture of distributions lists}
#' }
#' }
DiSCo_iter <- function(yy, df, id_col.target, M, G, T0, ...) {

    # obtaining the target state

    # target outcome
    target <- df[id_col == id_col.target, y_col]


    # generate list where each element contains a list of all micro-level outcomes for a control unit
    controls <- list()
    j <- 1

    controls.id <- unique(df[id_col != id_col.target, id_col])
    for (id in controls.id) {
      controls[j] <- list(df[id_col == id & t_col == yy, y_col])
      j <- j + 1
    }

    # evaluating the quantile functions on the grid "evgrid":
    controls.q <- matrix(0,nrow = length(evgrid), ncol=length(controls))
    for (jj in 1:length(controls)){
      controls.q[,jj] <- mapply(myquant, evgrid, MoreArgs = list(X=controls[[jj]]))
    }

    if (yy <= T0) {
      # obtaining the optimal weights for the DiSCo method
      DiSCo_res_weights <- DiSCo_weights_reg(controls, as.vector(target), M)

      DiSCo_res2 <- DiSCo_bc(controls, controls.q, DiSCo_res_weights,seq(from=0,to=1,length.out=M+1))
    }

    # sample grid
    grid <- list(grid.min = NA, grid.max = NA, grid.rand = NA, grid.ord = NA)
    grid[c("grid.min", "grid.max", "grid.rand", "grid.ord")] <- getGrid(target, controls, G)

    # getting the CDF from the quantile function
    DiSCo_res2.cdfF <- ecdf(controls.q)
    DiSCo_res2.cdf <- DiSCo_res2.cdfF(grid$grid.ord)



    # obtaining the optimal weights for the mixture of distributions method
    mixture <- getMixture(controls, target, grid$grid.min, grid$grid.max, grid$grid.rand)

    y_char <- as.character(yy)
    results <- list()
    results[["DiSCo"]] <-
      list("weights" = DiSCo_res_weights, "quantile.barycenter" = DiSCo_res2$barycenter, "cdf" = DiSCo_res2.cdf) # DiSCo estimator
    results[["mixture"]] <- list("weights" = mixture$weights.opt, "distance" = mixture$distance.opt, "mean" = mixture$mean) # mixture of distributions estimator
    results[["target"]] <- list("quantile" = target.s, "cdf" = mixture$target.order, "grid" =  grid.ord, "data" = as.vector(target))
    results[["controls"]] <- list("data" = controls, "cdf" = mixture$CDF.matrix, "id" = controls.id)

}
