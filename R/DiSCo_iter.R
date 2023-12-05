#' @title DiSCo_iter
#'
#' @description This function implements the DiSCo method for a single time period, as well as the mixture of distributions approach.
#'
#' @details This function is part of the DiSCo method, called for each time period. It calculates the optimal weights for the DiSCo method and the mixture of distributions approach for a single time period. The function processes data for both the target and control units, computes the quantile functions, and evaluates these on a specified grid. The function is designed to be used within the broader context of the DiSCo function, which aggregates results across multiple time periods.
#'
#' @param yy Integer indicating the current year being processed.
#' @param df Data frame containing the data for both target and control units.
#' @param evgrid Vector containing an evenly spaced grid on [0,1] for evaluating quantile functions.
#' @param id_col.target String specifying the column name for the target unit's identifier.
#' @param M Integer specifying the number of iterations for optimization.
#' @param G Integer indicating the grid size for evaluation.
#' @param T0 Integer indicating the first year of treatment as counted from 1.
#' @param ... Additional arguments passed to the function.
#'
#' @return A list with the following elements:
#' \itemize{
#'   \item{\code{DiSCo_weights} }{Weights calculated using the DiSCo method.}
#'   \item{\code{mixture} }{
#'     \itemize{
#'       \item{\code{weights} }{Optimal weights for the mixture approach.}
#'       \item{\code{distance} }{Value of the objective function for the mixture approach.}
#'       \item{\code{mean} }{Weighted mixture of the controls' CDFs.}
#'     }
#'   }
#'   \item{\code{target} }{
#'     \itemize{
#'       \item{\code{cdf} }{Empirical CDF of the target.}
#'       \item{\code{grid} }{Grid on which the quantile and CDF functions were evaluated.}
#'       \item{\code{data} }{Original data for the target unit.}
#'       \item{\code{quantiles} }{Quantiles for the target unit, evaluated on the specified grid.}
#'     }
#'   }
#'   \item{\code{controls} }{
#'     \itemize{
#'       \item{\code{data} }{Original data for the control units.}
#'       \item{\code{cdf} }{Empirical CDFs of the control units.}
#'       \item{\code{id} }{IDs of the control units.}
#'       \item{\code{quantiles} }{Quantiles for the control units, evaluated on the specified grid.}
#'     }
#'   }
#'   \item{\code{controls.q} }{Quantiles for the control units, evaluated on the specified grid.}
#' }
#' @export
DiSCo_iter <- function(yy, df, evgrid, id_col.target, M, G, T0, qmethod=NULL, ...) {

    # target
    target <- df[(id_col == id_col.target) & (t_col == yy)]$y_col

    # generate list where each element contains a list of all micro-level outcomes for a control unit
    controls <- list()
    j <- 1
    controls.id <- unique(df[id_col != id_col.target]$id_col)
    for (id in controls.id) {
      controls[[j]] <- df[(id_col == id) & (t_col == yy)]$y_col
      j <- j + 1
    }

    # check whether problem undetermined
    if (length(controls[[1]]) < length(controls)) {
      stop("Problem undetermined: number of data points is smaller than number of weights")
    }

    # evaluating the quantile functions on the grid "evgrid":
    controls.q <- matrix(0,nrow = length(evgrid), ncol=length(controls))
    for (jj in 1:length(controls)){
      # controls.q[,jj] <- mapply(myquant, evgrid, MoreArgs = list(X=controls[[jj]]))
      controls.q[,jj] <- myQuant(controls[[jj]], evgrid, qmethod)
    }

    # sample grid
    grid <- list(grid.min = NA, grid.max = NA, grid.rand = NA, grid.ord = NA)
    grid[c("grid.min", "grid.max", "grid.rand", "grid.ord")] <- getGrid(target, controls, G) # TODO: this can be done just once

    # obtaining the optimal weights for the DiSCo method
    DiSCo_res_weights <- DiSCo_weights_reg(controls, as.vector(target), M, qmethod=qmethod)

    # obtaining the optimal weights for the mixture of distributions method
    mixture <- DiSCo_mixture(controls, target, grid$grid.min, grid$grid.max, grid$grid.rand, M)

    #computing the target quantile function
    target.q <- myQuant(target, evgrid, qmethod)


    results <- list()
    results[["DiSCo"]] <- list("weights" = DiSCo_res_weights) # DiSCo estimator
    results[["mixture"]] <- list("weights" = mixture$weights.opt, "distance" = mixture$distance.opt, "mean" = mixture$mean) # mixture of distributions estimator
    results[["target"]] <- list("cdf" = mixture$target.order, "grid" = grid$grid.ord, "data" = as.vector(target), "quantiles" = target.q)
    results[["controls"]] <- list("cdf" = mixture$CDF.matrix, "data" = controls, "id" = controls.id, "quantiles" = controls.q)
    return(results)
}


