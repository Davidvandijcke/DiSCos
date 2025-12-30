# Plot distribution of treatment effects over time

Plot distribution of treatment effects over time

## Usage

``` r
plotDistOverTime(
  cdf_centered,
  grid_cdf,
  t_start,
  t_max,
  CI,
  ci_lower,
  ci_upper,
  ylim = c(0, 1),
  xlim = NULL,
  cdf = TRUE,
  xlab = "Distribution Difference",
  ylab = "CDF",
  obsLine = NULL,
  savePlots = FALSE,
  plotName = NULL,
  lty = 1,
  lty_obs = 1,
  t_plot = NULL
)
```

## Arguments

- cdf_centered:

  list of centered distributional statistics

- grid_cdf:

  grid

- t_start:

  start time

- t_max:

  maximum time

- CI:

  logical indicating whether to plot confidence intervals

- ci_lower:

  lower confidence interval

- ci_upper:

  upper confidence interval

- ylim:

  y limits

- xlim:

  x limits

- cdf:

  logical indicating whether to plot CDF or quantile difference

- xlab:

  x label

- ylab:

  y label

- obsLine:

  optional additional line to plot. Default is NULL which means no line
  is plotted.

- savePlots:

  logical indicating whether to save plots

- plotName:

  name of plot to save

- lty:

  line type for the main line passed as cdf_centered

- lty_obs:

  line type for the optional additional line passed as obsLine

- t_plot:

  optional vector of times to plot. Default is NULL which means all
  times are plotted.

## Value

plot of distribution of treatment effects over time
