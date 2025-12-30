# DiSCo_mixture

The alternative mixture of distributions approach in the paper

## Usage

``` r
DiSCo_mixture(controls1, target, grid.min, grid.max, grid.rand, M, simplex)
```

## Arguments

- controls1:

  A list of controls

- target:

  The target unit

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

- `cdf ` A matrix containing the CDFs of the target and control units
  evaluated on the grid.

- `distance.opt ` The optimal value of the Wasserstein distance.

- `mean ` The optimal value of the Wasserstein barycenter.

- `target.order ` The target unit, ordered.

- `weights.opt ` The optimal weights.
