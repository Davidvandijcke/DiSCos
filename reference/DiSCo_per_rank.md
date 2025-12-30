# DiSCo_per_rank

This function ranks the squared Wasserstein distances and returns the
p-values for each time period

## Usage

``` r
DiSCo_per_rank(distt, distp, T0)
```

## Arguments

- distt:

  List of squared Wasserstein distances between the target unit and the
  control units

- distp:

  List of squared Wasserstein distances between the control units

## Value

List of p-values for each time period
