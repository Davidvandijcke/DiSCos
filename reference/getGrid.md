# getGrid

Set up a grid for the estimation of the quantile functions and CDFs

## Usage

``` r
getGrid(target, controls, G)
```

## Arguments

- target:

  A vector containing the data for the target unit

- controls:

  A list containing the data for the control units

- G:

  The number of grid points

## Value

A list containing the following elements:

- `grid.min` The minimum value of the grid

- `grid.max` The maximum value of the grid

- `grid.rand` A vector containing the grid points

- `grid.ord` A vector containing the grid points, ordered
