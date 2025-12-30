# DiSCo_mixture_solve

The solver for the alternative mixture of distributions approach in the
paper

## Usage

``` r
DiSCo_mixture_solve(
  c_len,
  CDF.matrix,
  grid.min,
  grid.max,
  grid.rand,
  M,
  simplex
)
```

## Arguments

- c_len:

  The number of controls

- CDF.matrix:

  The matrix of CDFs

- grid.min:

  Minimal value of the grid on which the CDFs are evaluated.

- grid.max:

  Maximal value of the grid on which the CDFs are evaluated.

- grid.rand:

  Random grid on which the CDFs are evaluated.

- M:

  Integer indicating the number of control quantiles to use in the DiSCo
  method. Default is 1000.

- simplex:

  Logical, indicating whether to use to constrain the optimal weights to
  the unit simplex. Default is FALSE, which only constrains the weights
  to sum up to 1 but allows them to be negative.

## Value

A list containing the following elements:

- `distance.opt ` The optimal value of the Wasserstein distance.

- `mean ` The optimal value of the Wasserstein barycenter.

- `target.order ` The target unit, ordered.

- `weights.opt ` The optimal weights.
