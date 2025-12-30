# DiSCo_CI_iter

Function for computing the confidence intervals in the DiSCo method in a
single period

## Usage

``` r
DiSCo_CI_iter(
  t,
  controls_t,
  target_t,
  grid,
  T0,
  M = 1000,
  evgrid = seq(from = 0, to = 1, length.out = 1001),
  qmethod = NULL,
  qtype = 7,
  mixture = FALSE,
  simplex = FALSE,
  replace = TRUE
)
```

## Arguments

- t:

  Time period

- controls_t:

  List of control unit data for given period

- target_t:

  List of target unit data for given period

- grid:

  Grid to recompute the CDF on if `mixture` option is chosen

- T0:

  Index of the last pre-treatment period

- M:

  Integer indicating the number of control quantiles to use in the DiSCo
  method. Default is 1000.

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

The resampled counterfactual barycenter of the target unit
