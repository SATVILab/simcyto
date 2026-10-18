# simcyto 0.99.0

* Initial submission to Bioconductor.
* Simulation engine for generating synthetic single-cell cytometry data across
  multi-sample, multi-condition experimental designs.
* Supports Gaussian, t-distribution, and t-plus-Gaussian mixture distributions.
* Supports identity, gamma, Gaussian, and skew-normal marker transformations.
* Returns data as `flowCore::flowFrame` objects with per-cell ground-truth
  cluster labels.
* `simCytExperiment()` gains `covVaryBy`, choosing the strata (sample,
  condition, cluster) a random covariance matrix varies by, and `covMatList`,
  fixing each cluster's covariance matrix. The default `covVaryBy` keeps the
  previous per-tube draws.
* `simCytExperiment()` gains `covSd` and `covCor`, building covariance matrices
  from marker SDs and inter-marker correlations, each fixed or drawn from a
  range (for example correlations between 0 and 0.5), overall or per cluster.
* `simCytExperiment()` gains `meanPerturbation`: layers of mean shifts, each
  with an SD and the strata (sample, condition, cluster) it varies by.
* `simCytExperiment()` now errors when a stimulated condition's cluster
  probabilities (`probVecUns` plus its response vector) leave `[0, 1]` or do
  not sum to 1, instead of silently renormalising. Its documentation now
  explains how background and response combine across stimulated conditions.
* New `simCytCovMat()` builds a covariance matrix from marker SDs and
  correlations.
