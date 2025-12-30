# checks Carry out checks on the inputs

checks Carry out checks on the inputs

## Usage

``` r
checks(
  df,
  id_col.target,
  t0,
  M,
  G,
  num.cores,
  permutation,
  q_min,
  q_max,
  CI,
  boots,
  cl,
  graph,
  qmethod,
  seed
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

  logical, whether to use permutation or not

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

- seed:

  Integer, seed for the random number generator. This needs to be set
  explicitly in the function call, since it will invoke
  [`RNGkind`](https://rdrr.io/r/base/Random.html) which will set the
  seed for each core when using parallel processes. Default is NULL,
  which does not set a seed.
