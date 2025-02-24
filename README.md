
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DiSCos <img src="man/figures/logo.png" align="right" alt="" width="155" />

[![](https://cranlogs.r-pkg.org/badges/grand-total/DiSCos?color=blue)](https://cran.r-project.org/package=DiSCos)
[![](https://cranlogs.r-pkg.org/badges/last-month/DiSCos?color=blue)](https://cran.r-project.org/package=DiSCos)
[![](https://www.r-pkg.org/badges/version/DiSCos?color=blue)](https://cran.r-project.org/package=DiSCos)
[![](https://img.shields.io/badge/devel%20version-0.0.0.9000-blue.svg)](https://github.com/Davidvandijcke/DiSCos)
[![](https://img.shields.io/github/last-commit/Davidvandijcke/DiSCos.svg)](https://github.com/Davidvandijcke/DiSCos/commits/main)

<!-- README.mdƒ is generated from README.Rmd. Please edit that file -->

The **DiSCos** package contains tools for computing counterfactual
quantile functions in a Distributional Synthetic Controls (DiSco)
setting, following the method proposed in Gunsilius (2023).


## Getting Started

Have a look at the vignette replicating the [empirical
application](https://www.davidvandijcke.com/DiSCos/articles/Dube2019.html)
in the paper to get started.

## Installation

To install the latest stable version, run
```r
install.packages("DiSCos")

```

You can install latest development version from GitHub with:

``` r
devtools::install_github("Davidvandijcke/DiSCos")
```

If you find any bugs or have any questions, please email dvdijcke@umich.edu.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-gunsilius2023distributional" class="csl-entry">

If you are using this package in your research, please cite: 

Gunsilius, Florian F. 2023. “Distributional Synthetic Controls.”
*Econometrica* 91 (3): 1105–17.

Van Dijcke, David, Florian Gunsilius, and Austin Wright. "Return to Office and the Tenure Distribution." arXiv preprint arXiv:2405.04352 (2024).

</div>

</div>


## FAQ

`Q`: Why does the code sometimes run slower than expected?
`A`: The approach in DiSCo is non-parametric and requires integrating over entire distributions, which naturally increases computational complexity. The solution uses iterative methods (like in DiSCo_mixture) and can be costly for large datasets or large values of M and G.

`Q`: Is there a faster backend?
`A`: The optimization relies on libraries (e.g., “quadprog” or “CVXR”) that ultimately call C++ or FORTRAN under the hood. Even though these libraries are optimized, fully non-parametric distance-matching remains computationally heavy.

`Q`: How can I reduce the runtime while developing or debugging?
`A`:

Use fewer quantile points by setting M = 100 in DiSCo.
Set a smaller grid, e.g. G = 100, to reduce the evaluation points.
Start with only two time periods (Ts = 2) so there are fewer distributions to match.
Disable time-consuming features, such as confidence intervals (CI=FALSE) and permutation tests (permutation=FALSE).
Increase parallel execution with num.cores in DiSCo; for example, use a cluster or multiple cores to split bootstrap/permutation tasks (see mclapply.hack in R/utils.R).

`Q`: Are there other suggestions for handling large real-world data?
`A`:
• Verify correctness with a minimal working example, then gradually increase M, G, and the number of periods.
• Use discrete or categorical grid options (grid.cat) if suitable for your data to reduce the dimension of the problem.
• Ensure memory and CPU availability align with the size of the task.

`Q`: Can I get a progress bar?
`A`: When the package was implemented, there were no functional solutions for progress bars that work with parallel computation in R. This might have changed, or maybe you want to change it! If you are indeed that hero 🦸, feel free to contribute to this GitHub repo, or even to write your own parallel progress bar package! 
