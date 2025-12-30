# mclapply.hack

This function mimics forking (done with mclapply in Mac or Linux) for
the Windows environment. Designed to be used just like mclapply. Credit
goes to Nathan VanHoudnos.

## Usage

``` r
mclapply.hack(..., verbose = FALSE, mc.cores = 1)
```

## Arguments

- verbose:

  Should users be warned this is hack-y? Defaults to FALSE.

- mc.cores:

  Number of cores to use. Defaults to 1.

## See also

mclapply
