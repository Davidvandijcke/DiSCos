# DiSCo_per

Function to implement permutation test for Distributional Synthetic
Controls

## Usage

``` r
DiSCo_per(
  results.periods,
  T0,
  ww = 0,
  peridx = 0,
  evgrid = seq(from = 0, to = 1, length.out = 101),
  graph = TRUE,
  num.cores = 1,
  weights = NULL,
  qmethod = NULL,
  qtype = qtype,
  q_min = 0,
  q_max = 1,
  M = 1000,
  simplex = FALSE,
  mixture = FALSE
)
```

## Arguments

- results.periods:

  List of period-specific results from DiSCo

- T0:

  Integer indicating first year of treatment as counted from 1 (e.g, if
  treatment year 2002 was the 5th year in the sample, this parameter
  should be 5).

- ww:

  Optional vector of weights indicating the relative importance of each
  time period. If not specified, each time period is weighted equally.

- peridx:

  Optional integer indicating number of permutations. If not specified,
  by default equal to the number of units in the sample.

- graph:

  Logical, indicating whether to plot the permutation graph as in Figure
  3 of the paper. Default is FALSE.

- num.cores:

  Integer, number of cores to use for parallel computation. Default
  is 1. If the `permutation` or `CI` arguments are set to TRUE, this can
  be slow and it is recommended to set this to 4 or more, if possible.
  If you get an error in "all cores" or similar, try setting num.cores=1
  to see the precise error value.

- weights:

  Optional vector of weights to use for the "true" treated unit.
  `redo_weights` has to be set to FALSE for these weights to be used.

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

- M:

  Integer indicating the number of control quantiles to use in the DiSCo
  method. Default is 1000.

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

List of matrices containing synthetic time path of the outcome variable
for the target unit together with the time paths of the control units

## Details

This program iterates through all units and computes the optimal weights
on the other units for replicating the unit of iteration's outcome
variable, assuming that it is the treated unit. See Algorithm 1 in
Gunsilius (2023) for more details. The only modification is that we take
the ratio of post- and pre-treatment root mean squared Wasserstein
distances to calculate the p-value, rather than the level in each
period, following @abadie2010synthetic.

## References

Gunsilius FF (2023). “Distributional synthetic controls.”
*Econometrica*, **91**(3), 1105–1117.
