# Function for computing barycenters in the DiSCo method at every time period

Compute barycenters in the DiSCo method at every time period, as in
Definition 1, Step 4 in Gunsilius (2023) .

## Usage

``` r
DiSCo_bc(controls.q, weights, evgrid = seq(from = 0, to = 1, length.out = 101))
```

## Arguments

- controls.q:

  List with matrices of control quantile functions

- weights:

  Vector of optimal synthetic control weights, computed using the
  DiSCo_weights_reg function.

## Value

The quantile function of the barycenter associated with the "weights"
evaluated at the vector "evgrid"

## References

Gunsilius FF (2023). “Distributional synthetic controls.”
*Econometrica*, **91**(3), 1105–1117.
