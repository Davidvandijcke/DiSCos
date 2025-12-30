# bootCounterfactuals

Function for computing the bootstrapped counterfactuals in the DiSCo
method

## Usage

``` r
bootCounterfactuals(result_t, t, mixture, weights, evgrid, grid)
```

## Arguments

- result_t:

  A list containing the results of the DiSCo_CI_iter function

- t:

  The current time period

- mixture:

  Logical, indicating whether to use the mixture of distributions
  approach instead. See Section 4.3. in Gunsilius (2023) . This approach
  minimizes the distance between the CDFs instead of the quantile
  functions, and is preferred for categorical variables. When working
  with such variables, one should also provide a list of support points
  in the `grid.cat` parameter. When that is provided, this parameter is
  automatically set to TRUE. Default is FALSE.

- grid:

  Grid to recompute the CDF on if `mixture` option is chosen

## Value

A list containing the bootstrapped counterfactuals
