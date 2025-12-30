# Distributional Synthetic Controls

This function implements the distributional synthetic controls (DiSCo)
method from Gunsilius (2023) . as well as the alternative mixture of
distributions approach.

## Usage

``` r
DiSCo(
  df,
  id_col.target,
  t0,
  M = 1000,
  G = 1000,
  num.cores = 1,
  permutation = FALSE,
  q_min = 0,
  q_max = 1,
  CI = FALSE,
  boots = 500,
  replace = TRUE,
  uniform = FALSE,
  cl = 0.95,
  graph = FALSE,
  qmethod = NULL,
  qtype = 7,
  seed = NULL,
  simplex = FALSE,
  mixture = FALSE,
  grid.cat = NULL
)
```

## Arguments

- df:

  Data frame or data table containing the distributional data for the
  target and control units. The data table should contain the following
  columns:

  - `y_col ` A numeric vector containing the outcome variable for each
    unit. Units can be individuals, states, etc., but they should be
    nested within a larger unit (e.g. individuals or counties within a
    state)

  - `id_col ` A numeric vector containing the aggregate IDs of the
    units. This could be, for example, the state if the units are
    counties or individuals

  - `time_col ` A vector containing the time period of the observation
    for each unit. This should be a monotonically increasing integer.

- id_col.target:

  Variable indicating the name of the target unit, as specified in the
  id_col column of the data table. This variable can be any type, as
  long as it is the same type as the id_col column of the data table.

- t0:

  Integer indicating period of treatment.

- M:

  Integer indicating the number of control quantiles to use in the DiSCo
  method. Default is 1000.

- G:

  Integer indicating the number of grid points for the grid on which the
  estimated functions are evaluated. Default is 1000.

- num.cores:

  Integer, number of cores to use for parallel computation. Default
  is 1. If the `permutation` or `CI` arguments are set to TRUE, this can
  be slow and it is recommended to set this to 4 or more, if possible.
  If you get an error in "all cores" or similar, try setting num.cores=1
  to see the precise error value.

- permutation:

  Logical, indicating whether to use the permutation method for
  computing the optimal weights. Default is FALSE.

- q_min:

  Numeric, minimum quantile to use. Set this together with `q_max` to
  restrict the range of quantiles used to construct the synthetic
  control. Default is 0 (all quantiles). Currently NOT implemented for
  the `mixture` approach.

- q_max:

  Numeric, maximum quantile to use. Set this together with `q_min` to
  restrict the range of quantiles used to construct the synthetic
  control. Default is 1 (all quantiles). Currently NOT implemented for
  the `mixture` approach.

- CI:

  Logical, indicating whether to compute confidence intervals for the
  counterfactual quantiles. Default is FALSE. The confidence intervals
  are computed using the bootstrap procedure described in Van Dijcke et
  al. (2024) .

- boots:

  Integer, number of bootstrap samples to use for computing confidence
  intervals. Default is 500.

- replace:

  Logical, indicating whether to sample with replacement when computing
  the bootstrap samples. Default is TRUE.

- uniform:

  Logical, indicating whether to construct uniform bootstrap confidence
  intervals. Default is FALSE If FALSE, the confidence intervals are
  pointwise.

- cl:

  Numeric, confidence level for the (two-sided) confidence intervals.

- graph:

  Logical, indicating whether to plot the permutation graph as in Figure
  3 of the paper. Default is FALSE.

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

- seed:

  Integer, seed for the random number generator. This needs to be set
  explicitly in the function call, since it will invoke
  [`RNGkind`](https://rdrr.io/r/base/Random.html) which will set the
  seed for each core when using parallel processes. Default is NULL,
  which does not set a seed.

- simplex:

  Logical, indicating whether to use to constrain the optimal weights to
  the unit simplex. Default is FALSE, which only constrains the weights
  to sum up to 1 but allows them to be negative.

- mixture:

  Logical, indicating whether to use the mixture of distributions
  approach instead. See Section 4.3. in Gunsilius (2023) . This approach
  minimizes the distance between the CDFs instead of the quantile
  functions, and is preferred for categorical variables. When working
  with such variables, one should also provide a list of support points
  in the `grid.cat` parameter. When that is provided, this parameter is
  automatically set to TRUE. Default is FALSE.

- grid.cat:

  List, containing the discrete support points for a discrete grid to be
  used with the mixture of distributions approach. This is useful for
  constructing synthetic distributions for categorical variables.
  Default is NULL, which uses a continuous grid based on the other
  parameters.

## Value

A list containing the following elements:

- `results.periods` A list containing, for each time period, the
  elements described in the return argument of
  [`DiSCo_iter`](http://www.davidvandijcke.com/DiSCo/reference/DiSCo_iter.md),
  as well as the following additional elements:

  - `DiSco`

    - `quantile ` The counterfactual quantiles for the target unit.

    - `weights ` The optimal weights for the target unit.

    - `cdf ` The counterfactual CDF for the target unit.

- `weights` A numeric vector containing the synthetic control weights
  for the control units, averaged over time. When `mixture` is TRUE,
  these are the weights for the mixture of distributions, otherwise they
  are the weights for the quantile-based approach.

- `CI` A list containing the confidence intervals for the counterfactual
  quantiles and CDFs, if `CI` is TRUE. Each element contains two named
  subelements called `upper`, `lower`, `se` which are the upper and
  lower confidence bands and the standard error of the estimate,
  respectively. They are G x T matrices where G is the specified number
  of grid points and T is the number of time periods. The elements are:

  - `cdf` The bootstrapped CDF

  - `quantile` The bootstrapped quantile

  - `quantile_diff` The bootstrapped quantile difference

  - `cdf_diff` The bootstrapped CDF difference

  - `bootmat` A list containing the raw bootstrapped samples for the
    counterfactual quantiles and CDFs, if `CI` is TRUE. These are not
    meant to be accessed directly, but are used by `DiSCoTEA` to compute
    aggregated standard errors. Advanced users may wish to access these
    directly for further analysis. The element names should be
    self-explanatory. \#'

  - `control_ids` A list containing the control unit IDs used for each
    time period, which can be used to identify the weights associated
    with each control as the returned weights have the same order as the
    control IDs.

  - `perm ` A
    [`permut`](http://www.davidvandijcke.com/DiSCo/reference/permut.md)
    object containing the results of the permutation method, if
    `permutation` is TRUE. Call `summary` on this object to print the
    overall results of the permutation test. \#'

  - `evgrid` A numeric vector containing the grid points on which the
    quantiles were evaluated.

  - `params` A list containing the parameters used in the function call.

## Details

This function is called for every time period in the DiSCo function. It
implements the DiSCo method for a single time period, as well as the
mixture of distributions approach. The corresponding results for each
time period can be accessed in the `results.periods` list of the output
of the DiSCo function. The DiSCo function returns the average weight for
each unit across all periods, calculated as a uniform mean, as well as
the counterfactual target distribution produced as the weighted average
of the control distributions for each period, using these averaged
weights.

## References

Gunsilius FF (2023). “Distributional synthetic controls.”
*Econometrica*, **91**(3), 1105–1117.  
  
Van Dijcke D, Gunsilius F, Wright AL (2024). “Return to Office and the
Tenure Distribution.” Working Paper 2024-56, University of Chicago,
Becker Friedman Institute for Economics. ()
