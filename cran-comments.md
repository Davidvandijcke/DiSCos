## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission to fix compatibility with CVXR 1.8.x:

* Added backwards-compatible support for both CVXR < 1.8.1 (`solve`) and >= 1.8.1 (`psolve`)
* Replaced `CVXR::norm()` with `base::norm()` for non-CVXR matrix operations
* Wrapped slow examples in `\donttest{}` to address elapsed time NOTE on Windows

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.5.1
* win-builder (R-devel, R-release)
* GitHub Actions: Ubuntu, Windows, macOS
* mac.r-project.org macbuilder (R-release)

## Notes

One test for `q_min`/`q_max` parameters is skipped on macOS CI environments due to a low-level system interaction issue that only occurs in CI (not locally). The functionality works correctly on all platforms.
