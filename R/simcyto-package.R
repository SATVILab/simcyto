#' simcyto: Cytometry Data Simulation Engine
#'
#' @description
#' The `simcyto` package provides tools for simulating single-cell cytometry data across
#' multi-sample, multi-condition experimental designs with defined phenotype clusters,
#' mixture distributions, perturbations, and transformations.
#'
#' @details
#' The simulation engine acknowledges the upstream provenance line represented by
#' `SATVILab/StimGate`: the reusable simulation machinery in `simcyto` is derived
#' from the FAUST-inspired phenotype simulation workflow used in that project, and
#' specifically from the `functionsForBenchmarking-Pheno.R` source in
#' `RGLab/faust_manuscript_analyses`. Reusable components were audited and
#' migrated into `simcyto`:
#' \itemize{
#'   \item **Scenario Builders**: \code{\link{simCytBuildScenario}}, \code{\link{simCytScenarioUnivariate}},
#'         and \code{\link{simCytScenarioBivariate}} for defining cluster mean matrices and response vectors.
#'   \item **Simulation Core**: \code{\link{simCytExperiment}}, \code{simCytSample},
#'         and \code{simCytCondition} for multi-sample and multi-condition dataset simulation.
#'   \item **Mixture Distributions**: \code{simCytClusterData} supporting \code{"gaussianOnly"},
#'         \code{"tOnly"}, and \code{"tPlusGauss"}.
#'   \item **Transformations & Ratio Corrections**: \code{\link{simCytGetTransformation}} and
#'         upper-population ratio correction for gamma and skew transformations.
#' }
#' Downstream statistical evaluation and mixed-model benchmarking helpers were intentionally
#' excluded to keep the package focused purely on cytometry simulation.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
