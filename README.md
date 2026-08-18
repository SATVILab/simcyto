# simcyto

`simcyto` is an R package for simulating cytometry data.

## Installation

You can install the development version of `simcyto` from GitHub with:

```r
# install.packages("pak")
pak::pak("SATVILab/simcyto")
```

## Phenotype Simulation Machinery & Function Audit

The phenotype simulation machinery in `simcyto` derives from the simulation components used in the FAUST benchmarking framework (`functionsForBenchmarking-Pheno.R`). Functions in that script have been audited and classified as follows:

1. **Reusable Simulation Machinery (Migrated to `simcyto`)**:
   - **Scenario Construction**: `simCytBuildScenario()`, `simCytScenarioUnivariate()`, `simCytScenarioBivariate()` construct marker expression grids (`meanExprMat`), cluster labels (`clusterLabelVec`), baseline probabilities (`probVecUns`), and response vectors (`probResponseVecByStimCondition`).
   - **Cluster Data & Mixture Distributions**: `simCytClusterData()` and `simCytCluster()` support `"gaussianOnly"`, `"tOnly"`, and `"tPlusGauss"` (heavy-tailed/skewed) mixture components with positive-definite covariance matrix generation (`.posDef()`).
   - **Perturbations**: Multi-level hierarchy incorporating sample-level, condition-level, and cluster-level perturbations (`samplePerturbationSd`, `conditionPerturbationSd`, `clusterPerturbationSd`).
   - **Ratio Correction**: Upper-population ratio adjustment (`.simCytRatioAdjustUpper()`, `.simCytUsesUpperRatioCorrection()`) preserving distance-to-spread ratios post-transformation for gamma and skew transformations.
   - **Cell Allocation**: `probExact = TRUE` (exact deterministic integer assignment) or `probExact = FALSE` (multinomial sampling).
   - **Experiment Hierarchy**: `simCytExperiment()` -> `simCytSample()` -> `simCytCondition()` -> `simCytCluster()`.

2. **FAUST-Specific Interoperability**:
   - `flowFrame` and cluster label generation (compatible with FAUST phenotype label annotations) returned via `flowFrameList` and `labelsList` in `simCytExperiment()`.

3. **Excluded Non-Simulation Machinery**:
   - Downstream statistical modeling, mixed-model significance testing, and algorithm-evaluation/benchmarking code from `functionsForBenchmarking-Pheno.R` were intentionally excluded to keep `simcyto` focused strictly on cytometry data simulation.

## Simulator Architecture & Roadmap

`simcyto` provides a unified simulation hierarchy where scenario builders produce inputs consumed directly by `simCytExperiment()`:

```r
library(simcyto)

# 1. Build a scenario (e.g. 2-marker / 4-cluster)
sc <- simCytScenarioBivariate()

# 2. Run an experiment simulation
sim_res <- simCytExperiment(
  nSample = 2,
  nMarker = sc$nMarker,
  nCondition = 2,
  nCluster = sc$nCluster,
  nCellByCondition = 100,
  transformationFunc = simCytTransformIdentity(),
  meanExprMat = sc$meanExprMat,
  clusterLabelVec = sc$clusterLabelVec,
  probVecUns = sc$probVecUns,
  probResponseVecByStimCondition = sc$probResponseVecByStimCondition
)
```

Future work will consolidate high-level scenario building directly within `simCytExperiment()` options while preserving explicit parity and backwards compatibility.

## Repository Structure

- `R/`: R package source code.
- `tests/testthat/`: Unit tests for the package.
- `DESCRIPTION`: Package metadata.
- `NAMESPACE`: Exported package functions and dependencies.
