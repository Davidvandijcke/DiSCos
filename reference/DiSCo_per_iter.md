# DiSCo_per_iter

This function performs one iteration of the permutation test

## Usage

``` r
DiSCo_per_iter(
  c_df,
  c_df.q,
  t_df,
  T0,
  peridx,
  evgrid,
  idx,
  grid_df,
  M = 1000,
  ww = 0,
  qmethod = NULL,
  qtype = 7,
  q_min = 0,
  q_max = 1,
  simplex = FALSE,
  mixture = FALSE
)
```

## Arguments

- c_df:

  List of control units

- c_df.q:

  List of quantiles of control units

- t_df:

  List of target unit

- idx:

  Index of permuted target unit

- grid_df:

  Grids to evaluate CDFs on, only needed when `mixture=TRUE`

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

## Value

List of squared Wasserstein distances between the target unit and the
control units
