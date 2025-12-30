# parseBoots

Function for parsing the bootstrapped counterfactuals in the DiSCo
method

## Usage

``` r
parseBoots(CI_temp, cl, q_disco, cdf_disco, q_obs, cdf_obs, uniform = TRUE)
```

## Arguments

- CI_temp:

  A list containing the bootstrapped counterfactuals

- cl:

  The confidence level

- q_disco:

  The estimated quantiles around which to center

- cdf_disco:

  The estimated cdfs around which to center

- q_obs:

  The observed quantiles

- cdf_obs:

  The observed cdfs

- uniform:

  Whether to use uniform or pointwise confidence intervals

## Value

A list containing the confidence intervals for the quantiles and cdfs
