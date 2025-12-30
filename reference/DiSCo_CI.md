# DiSCo_CI

Function for computing the confidence intervals in the DiSCo method
using the bootstrap approach described in

## Usage

``` r
DiSCo_CI(
  redraw,
  controls,
  target,
  T_max,
  T0,
  grid,
  mc.cores = 1,
  evgrid = seq(from = 0, to = 1, length.out = 1001),
  qmethod = NULL,
  qtype = 7,
  M = 1000,
  mixture = FALSE,
  simplex = FALSE,
  replace = TRUE
)
```

## Arguments

- redraw:

  Integer indicating the current bootstrap redraw

- controls:

  A list containing the raw data for the control group

- target:

  A list containing the raw data for the target group

- T_max:

  Index of last time period

- T0:

  Index of the last pre-treatment period

- grid:

  Grid to recompute the CDF on if `mixture` option is chosen

- mc.cores:

  Number of cores to use for parallelization

- qmethod:

  Character, indicating the method to use for computing the quantiles of
  the target distribution. The default is NULL, which uses the
  [`quantile`](https://rdrr.io/r/stats/quantile.html) function from the
  stats package. Other options are
  "[`qkden`](https://rdrr.io/pkg/evmix/man/kden.html)" (based on
  smoothed kernel density function) and
  "[`extreme`](https://rdrr.io/pkg/extremeStat/man/distLquantile.html)"
  (based on parametric extreme value distributions). Both are
  substantially slower than the default method but may be useful for
  fat-tailed distributions with few data points at the upper quantiles.
  Alternatively, one could use the q_max option to restrict the range of
  quantiles used.

- qtype:

  Integer, indicating the type of quantile to compute when using
  [`quantile`](https://rdrr.io/r/stats/quantile.html) in the `qmethod`
  argument. The default 7. See the documentation for the
  [`quantile`](https://rdrr.io/r/stats/quantile.html) function for more
  information.

- M:

  Integer indicating the number of control quantiles to use in the DiSCo
  method. Default is 1000.

- mixture:

  Logical, indicating whether to use the mixture of distributions
  approach instead. See Section 4.3. in Gunsilius (2023) . This approach
  minimizes the distance between the CDFs instead of the quantile
  functions, and is preferred for categorical variables. When working
  with such variables, one should also provide a list of support points
  in the `grid.cat` parameter. When that is provided, this parameter is
  automatically set to TRUE. Default is FALSE.

- simplex:

  Logical, indicating whether to use to constrain the optimal weights to
  the unit simplex. Default is FALSE, which only constrains the weights
  to sum up to 1 but allows them to be negative.

- replace:

  Logical, indicating whether to sample with replacement when computing
  the bootstrap samples. Default is TRUE.

## Value

A list with the following components

- `weights` The bootstrapped weights

- `disco_boot` A list containing the bootstrapped counterfactuals, with
  the following elements, each of which contains named elements called
  `upper` and `lower` which are G x T matrices where G is the specified
  number of grid points and T is the number of time periods
