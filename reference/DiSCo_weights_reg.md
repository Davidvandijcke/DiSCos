# DiSCo_weights_reg

Function for obtaining the weights in the DiSCo method at every time
period

## Usage

``` r
DiSCo_weights_reg(
  controls,
  target,
  M = 500,
  qmethod = NULL,
  qtype = 7,
  simplex = FALSE,
  q_min = 0,
  q_max = 1
)
```

## Arguments

- controls:

  List with matrices of control distributions

- target:

  Matrix containing the target distribution

- M:

  Integer indicating the number of control quantiles to use in the DiSCo
  method. Default is 1000.

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

- simplex:

  Logical, indicating whether to use to constrain the optimal weights to
  the unit simplex. Default is FALSE, which only constrains the weights
  to sum up to 1 but allows them to be negative.

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

## Value

Vector of optimal synthetic control weights

## Details

Estimate the optimal weights for the distributional synthetic controls
method. solving the convex minimization problem in Eq. (2) in Gunsilius
(2023) .. using a regression of the simulated target quantile on the
simulated control quantiles, as in Eq. (3), \\\underset{\vec{\lambda}
\in \Delta^J}{\operatorname{argmin}}\left\\\mathbb{Y}\_t
\vec{\lambda}\_t-\vec{Y}\_{1 t}\right\\\_2^2\\. For the constrained
optimization we rely on the package pracma the control distributions can
be given in list form, where each list element contains a vector of
observations for the given control unit, in matrix form; in matrix- each
column corresponds to one unit and each row is one observation. The
list-form is useful, because the number of draws for each control group
can be different. The target must be given as a vector.

## References

Gunsilius FF (2023). “Distributional synthetic controls.”
*Econometrica*, **91**(3), 1105–1117.
