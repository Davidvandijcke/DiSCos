# permut

Object to hold results of permutation test

## Usage

``` r
permut(distp, distt, p_overall, J_1, q_min, q_max, plot)
```

## Arguments

- distp:

  List of squared Wasserstein distances between the control units

- distt:

  List of squared Wasserstein distances between the target unit and the
  control units

- p_overall:

  Overall p-value

- J_1:

  Number of control units

- q_min:

  Minimum quantile

- q_max:

  Maximum quantile

- plot:

  ggplot object containing plot of squared Wasserstein distances over
  time for all permutations.

## Value

A list of class permut, with the same elements as the input arguments.
