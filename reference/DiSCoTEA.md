# Aggregate treatment effects from DiSCo function.

Function to aggregate treatment effects from the output of the DiSCo
function, plot the distribution of the aggregation statistic over time,
and report summary tables.

## Usage

``` r
DiSCoTEA(
  disco,
  agg = "quantileDiff",
  graph = TRUE,
  t_plot = NULL,
  savePlots = FALSE,
  xlim = NULL,
  ylim = NULL,
  samples = c(0.25, 0.5, 0.75)
)
```

## Arguments

- disco:

  Output of the DiSCo function.

- agg:

  String indicating the aggregation statistic to be used. Options
  include

  - `quantileDiff` Difference in quantiles between the target and the
    weighted average of the controls.

  - `quantile` Plots both the observed and the counterfactual quantile
    functions. No summary statistics will be produced.

  - `cdfDiff` Difference in CDFs between the target and the weighted
    average of the controls.

  - `cdf` Plots both the observed and the counterfactual CDFs. No
    summary statistics will be produced.

- graph:

  Boolean indicating whether to plot graphs (default is TRUE).

- t_plot:

  Optional vector of time periods (`t_col` values in the original
  dataframe) to be plotted (default is NULL, which plots all time
  periods).

- savePlots:

  Boolean indicating whether to save the plots to the current working
  directory (default is FALSE). The plot names will be
  `[agg]_[start_year]_[end_year].pdf`.

- xlim:

  Optional vector of length 2 indicating the x-axis limits of the plot.
  Useful for zooming in on relevant parts of the distribution for
  fat-tailed distributions.

- ylim:

  Optional vector of length 2 indicating the y-axis limits of the plot.

- samples:

  Numeric vector indicating the range of quantiles of the aggregation
  statistic (`agg`) to be summarized in the `summary` property of the S3
  class returned by the function (default is c(0.25, 0.5, 0.75)). For
  example, if `samples` = c(0.25, 0.5, 0.75), the summary table will
  include the average effect for the 0-25th, 25-50th, 50-75th and
  75-100th quantiles of the distribution of the aggregation statistic
  over time.

## Value

A [`DiSCoT`](http://www.davidvandijcke.com/DiSCo/reference/DiSCoT.md)
object, which is an S3 class that stores a list of treatment effects,
their standard errors, the corresponding confidence intervals (if
specified), and a dataframe with treatment effects aggregated according
to the `agg` input. The S3 class also has a `summary` property that will
print a selection of aggregated effects (specified by the `samples`
parameter) for the chosen `agg` method, by post-treatment year, as well
as the permutation test results, if specified.

## Details

This function takes in the output of the DiSCo_per function and computes
aggregate treatment effect using a user-specified aggregation statistic.
The default is the differences between the counterfactual and the
observed quantile functions (`quantileDiff`). If `graph` is set to TRUE,
the function will plot the distribution of the aggregation statistic
over time. The S3 class returned by the function has a `summary`
property that will print a selection of aggregated effects (specified by
the `samples` parameter) for the chosen `agg` method, by post-treatment
year (see examples below). This `summary` call will only print effects
if the `agg` parameter requested a distribution difference
(`quantileDiff` or `cdfDiff`). The other aggregations are meant to be
inspected visually. If the `permutation` parameter was set to TRUE in
the original `DiSCo` call, the summary table will include the results of
the permutation test. If the original `DiSCo` call was restricted to a
range of quantiles smaller than `[0,1]` (i.e. `q_min` \> 0 or `q_max` \<
1), the `samples` parameter is ignored and only the aggregated
differences for the quantile range specified in the original call are
returned.
