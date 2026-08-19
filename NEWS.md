# simcyto 0.99.0

* Initial submission to Bioconductor.
* Simulation engine for generating synthetic single-cell cytometry data across
  multi-sample, multi-condition experimental designs.
* Supports Gaussian, t-distribution, and t-plus-Gaussian mixture distributions.
* Supports identity, gamma, Gaussian, and skew-normal marker transformations.
* Returns data as `flowCore::flowFrame` objects with per-cell ground-truth
  cluster labels.
