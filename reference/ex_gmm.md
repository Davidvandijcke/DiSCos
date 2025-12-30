# ex_gmm

Example data for `DiSCo` command. Returns simulated target and control
that are mixtures of Gaussian distributions.

## Usage

``` r
ex_gmm(Ts = 2, num.con = 30, numdraws = 1000)
```

## Arguments

- Ts:

  an integer indicating the number of time periods

- num.con:

  an integer indicating the number of control units

- numdraws:

  an integer indicating the number of draws

## Value

- `target`:

  a vector.

- `control`:

  a matrix.
