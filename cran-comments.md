# DiSCos 0.1.4

This is a performance release. Bootstrap confidence intervals are roughly 3x
faster, and the permutation test and point estimation are also faster: each
sample is now sorted once and quantiles/CDFs are evaluated on the sorted data
through internal kernels that replicate `stats::quantile(type = 7)`,
`stats::ecdf()`, and the CDF inversion exactly, and the confidence-interval
assembly is vectorized. The random-number draw sequence is unchanged, so
results are numerically identical to 0.1.3 for the same seed. There are no
API changes (one new optional argument with a backward-compatible default on
an internal helper).

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.5.1
* win-builder (R-devel)

## Notes

One test for `q_min`/`q_max` parameters is skipped on macOS CI environments
due to a low-level system interaction issue that only occurs in CI (not
locally). The functionality works correctly on all platforms.
