# DiSCos 0.2.0

* Added `perm_q_range` argument to `DiSCo()`: restrict the permutation test statistic to a quantile
  sub-range (e.g. the upper tail) without changing the synthetic-control fit, so one can test for an
  effect in a specific region of the distribution.
* Added `perm_seed` argument to `DiSCo()`: seed each permutation iteration deterministically by its
  index, making the permutation p-value reproducible and independent of `num.cores`.

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


