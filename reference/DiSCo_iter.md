# Estimate DiSCo in a single period

This function implements the DiSCo method for a single time period, as
well as the mixture of distributions approach. Its return values contain
valuable period-specific estimation outputs.

## Usage

``` r
DiSCo_iter(
  yy,
  df,
  evgrid,
  id_col.target,
  M,
  G,
  T0,
  qmethod = NULL,
  qtype = 7,
  q_min = 0,
  q_max = 1,
  simplex = FALSE,
  controls.id,
  grid.cat,
  mixture
)
```

## Arguments

- yy:

  Integer indicating the current year being processed.

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

- evgrid:

  A vector of grid points on which to evaluate the quantile functions.

- id_col.target:

  Variable indicating the name of the target unit, as specified in the
  id_col column of the data table. This variable can be any type, as
  long as it is the same type as the id_col column of the data table.

- M:

  Integer indicating the number of control quantiles to use in the DiSCo
  method. Default is 1000.

- G:

  Integer indicating the number of grid points for the grid on which the
  estimated functions are evaluated. Default is 1000.

- T0:

  Integer indicating the last pre-treatment period starting from 1.

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

- controls.id:

  List of strings specifying the column names for the control units'
  identifiers.

- grid.cat:

  List, containing the discrete support points for a discrete grid to be
  used with the mixture of distributions approach. This is useful for
  constructing synthetic distributions for categorical variables.
  Default is NULL, which uses a continuous grid based on the other
  parameters.

- mixture:

  Logical, indicating whether to use the mixture of distributions
  approach instead. See Section 4.3. in Gunsilius (2023) . This approach
  minimizes the distance between the CDFs instead of the quantile
  functions, and is preferred for categorical variables. When working
  with such variables, one should also provide a list of support points
  in the `grid.cat` parameter. When that is provided, this parameter is
  automatically set to TRUE. Default is FALSE.

## Value

A list with the following elements:

- `DiSCo_weights ` Weights calculated using the DiSCo method.

- `mixture `

  - `weights ` Optimal weights for the mixture approach.

  - `distance ` Value of the objective function for the mixture
    approach.

  - `mean ` Weighted mixture of the controls' CDFs.

- `target `

  - `cdf ` Empirical CDF of the target. Only computed when
    `mixture=TRUE`.

  - `grid ` Grid on which the quantile and CDF functions were evaluated.

  - `data ` Original data for the target unit.

  - `quantiles ` Quantiles for the target unit, evaluated on the
    specified grid.

- `controls `

  - `data ` Original data for the control units.

  - `cdf ` Empirical CDFs of the control units. Only computed when
    `mixture=TRUE`.

  - `quantiles ` Quantiles for the control units, evaluated on the
    specified grid. .

- `controls.q ` Quantiles for the control units, evaluated on the
  specified grid.

## Details

This function is part of the DiSCo method, called for each time period.
It calculates the optimal weights for the DiSCo method and the mixture
of distributions approach for a single time period. The function
processes data f or both the target and control units, computes the
quantile functions, and evaluates these on a specified grid. The
function is designed to be used within the broader context of the DiSCo
function, which aggregates results across multiple time periods.
