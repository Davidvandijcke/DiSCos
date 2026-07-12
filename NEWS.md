# DiSCos 0.1.4

* Performance release: the bootstrap confidence intervals are about 3x faster
  and the permutation test and point estimation meaningfully faster, from
  sorting each sample once and evaluating quantiles/CDFs on the sorted data,
  and from vectorizing the confidence-interval assembly. Results are
  numerically identical to 0.1.3 for the same seed.

# DiSCos 0.1.3

* Compatibility with CVXR 1.8.x (`solve` -> `psolve`, `getValue` -> `value`)

# DiSCos 0.1.2

* Compatibility with ggplot 4.0.0

# DiSCos 0.1.1

* Corrected ORCIDs of authors

# DiSCos 0.1.0

Major updates:
* Implemented the fully-fledged bootstrap approach for confidence interval estimation. See Van Dijcke, Gunsilius, and Wright (2024)
* Added mixture of distributions method CDFs, confidence intervals, and permutation test 

Minor updates:
* Fixed a small issue where plot=FALSE option didn't work in DiSCoTEA
* Added option to supply ordinal vector of discrete points for mixture of distributions


# DiSCos 0.0.1

* First version of package, basic DiSco function for estimating counterfactual quantile function, resampling-based CIs, permutation test; DiSCo_TEA function for various aggregation schemes. 


