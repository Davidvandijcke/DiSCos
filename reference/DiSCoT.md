# Store aggregated treatment effects

S3 object holding aggregated treatment effects

## Usage

``` r
DiSCoT(
  agg,
  treats,
  ses,
  grid,
  ci_lower,
  ci_upper,
  t0,
  call,
  cl,
  N,
  J,
  agg_df,
  perm,
  plot
)
```

## Arguments

- agg:

  aggregation method

- treats:

  list of treatment effects

- ses:

  list of standard errors

- grid:

  grid

- ci_lower:

  list of lower confidence intervals

- ci_upper:

  list of upper confidence intervals

- t0:

  start time

- call:

  call

- cl:

  confidence level

- N:

  number of observations

- J:

  number of treated units

- agg_df:

  dataframe of aggregated treatment effects and their confidence
  intervals

- perm:

  list of per mutation results

- plot:

  a ggplot object containing the plot for the aggregated treatment
  effects using the `agg` parameter

## Value

S3 object of class `DiSCoT` with associated `summary` and `print`
methods
