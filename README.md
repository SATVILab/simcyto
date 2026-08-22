# simcyto

`simcyto` is an R package for simulating cytometry data.

## Installation

You can install the development version of `simcyto` from GitHub with:

```r
# install.packages("pak")
pak::pak("SATVILab/simcyto")
```

## Acknowledgment and provenance

`simcyto` acknowledges the upstream provenance line represented by the
`SATVILab/StimGate` project. The reusable simulation machinery in this package is
derived from the same FAUST-inspired stimulation/phenotype simulation workflow
used in `StimGate`; the public upstream source identified in the
`RGLab/faust_manuscript_analyses` repository is
`faustRuns/simulation/functionsForBenchmarking.R` at ref
`f86f1ab19a41fd2690bf180a7dcf483f9552950c`. The downstream project-local
`functionsForBenchmarking-Pheno.R` script is a corresponding adaptation of the
same logic. This provenance is recorded in `inst/COPYRIGHTS` and should be
treated as the package's authoritative redistribution record.

The positive-definite covariance helper used in the cluster simulation core is a
nested third-party provenance item: the construction is the standard R-help /
R programming community algorithm credited to Ravi Varadhan, and it is retained
here as a small numerical utility that is not itself a separate software
package.

## Legacy StimGate regression fixtures

The package keeps a compact, pinned set of regression fixtures under
`tests/testthat/fixtures/` for cross-repository parity checks. These fixtures are
intentionally small enough for routine test runs but cover migration-critical
paths from the legacy FAUST-derived workflow, including exact allocations and
mixed-distribution simulation. Current StimGate transformation semantics are
protected separately by direct formula-parity tests in `tests/testthat/`. The canonical fixture
provenance is recorded as `SATVILab/StimGate` outputs derived from the legacy
FAUST workflow in `RGLab/faust_manuscript_analyses` at
`f86f1ab19a41fd2690bf180a7dcf483f9552950c`, preserving the upstream implementation
reference while keeping the regression expectations self-contained in `simcyto`.

## Phenotype Simulation Machinery & Function Audit

The phenotype simulation machinery in `simcyto` derives from the simulation components used in the FAUST benchmarking framework (`functionsForBenchmarking-Pheno.R`). Functions in that script have been audited and classified as follows:

1. **Reusable Simulation Machinery (Migrated to `simcyto`)**:
   - **Scenario Construction**: `simCytBuildScenario()`, `simCytScenarioUnivariate()`, `simCytScenarioBivariate()` construct marker expression grids (`meanExprMat`), cluster labels (`clusterLabelVec`), baseline probabilities (`probVecUns`), and response vectors (`probResponseVecByStimCondition`).
   - **Cluster Data & Mixture Distributions**: `simCytClusterData()` and
     `simCytCluster()` support `"gaussianOnly"`, `"tOnly"`, and
     `"tPlusGauss"` (heavy-tailed/skewed) mixture components with
     positive-definite covariance matrix generation (`.posDef()`).
   - **Post-simulation Transformations**:
     `simCytTransformGaussian()` reproduces StimGate's no-op Gaussian setting,
     `simCytTransformGamma()` reproduces `gamma(1 + abs(x / 4))`,
     `simCytTransformGammaFixed()` applies that gamma transform then rescales
     to the original mean and SD, and `simCytTransformSkew()` reproduces the
     current StimGate `sinh(...)`/gamma-divisor skew helper.
    - **Perturbations & Condition Mismatches**: Multi-level hierarchy incorporating sample-level, condition-level, and cluster-level perturbations (`samplePerturbationSd`, `conditionPerturbationSd`, `clusterPerturbationSd`), deterministic post-transformation stimulated condition location shifts (`stimMeanShift`, optionally targeted via `stimMeanShiftClusters`), and within-component SD scaling (`stimSdMultiplier`, optionally targeted via `stimSdMultiplierClusters`).
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

# 2. Run an experiment simulation using the scenario contract
sim_res <- simCytExperiment(
  scenario = sc,
  nSample = 2,
  nCondition = 2,
  nCellByCondition = 100,
  transformationFunc = simCytTransformIdentity(),
  mixtureType = "gaussianOnly",
  probExact = FALSE
)

# 3. Inspect expression matrices and per-cell truth labels
names(sim_res$flowFrameList)
head(flowCore::exprs(sim_res$flowFrameList[[1]]))
table(sim_res$labelsList[[1]])
```

High-level simulation now supports a compact scenario contract directly in `simCytExperiment()`: users may pass a scenario list created by a scenario builder and provide only the experiment-level settings that are not redundant with the phenotype grid or probability tables.

`mixtureType` and `transformationFunc` control different parts of the model:

- `mixtureType` chooses the raw simulation distribution (`"gaussianOnly"`,
  `"tOnly"`, or `"tPlusGauss"`).
- `transformationFunc` applies a post-simulation transformation to those raw
  values. For StimGate migration, the public `"gaussian"`, `"gamma"`,
  `"gammaFixed"`, and `"skew"` constructors intentionally match the current
  StimGate helper semantics rather than adding random noise.

## Repository Structure

- `R/`: R package source code.
- `tests/testthat/`: Unit tests for the package.
- `DESCRIPTION`: Package metadata.
- `NAMESPACE`: Exported package functions and dependencies.
