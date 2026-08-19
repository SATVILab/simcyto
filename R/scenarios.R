#' Build a simulation scenario for `nMarker` markers
#'
#' Construct the phenotype grid and probability objects consumed by
#' [simCytExperiment()]. The scenario contract is currently tied to binary marker
#' combinations, so the number of populations is fixed at `nCluster = 2^nMarker`.
#'
#' @param nMarker Integer scalar. Number of markers. Must be `>= 1`.
#' @param lowMean Numeric scalar. Mean expression assigned to marker-negative
#'   populations.
#' @param highMean Numeric scalar. Mean expression assigned to marker-positive
#'   populations.
#' @param probUns Optional numeric vector of length `2^nMarker`. Baseline
#'   unstimulated probabilities; must sum to 1. If `NULL`, a uniform distribution
#'   is used.
#' @param probResponse Optional list (or numeric vector) defining stimulated
#'   probability shifts. Each response vector must have length `2^nMarker` and is
#'   added to `probUns` inside [simCytExperiment()]. You are responsible for
#'   ensuring each resulting stimulated probability vector is valid (non-negative
#'   and summing to 1).
#'
#' @return A list with elements `nMarker`, `nCluster`, `meanExprMat`,
#'   `clusterLabelVec`, `probVecUns`, and `probResponseVecByStimCondition`.
#'
#' @examples
#' scenario <- simCytBuildScenario(
#'   nMarker = 2,
#'   lowMean = 0,
#'   highMean = 5,
#'   probUns = c(0.7, 0.1, 0.1, 0.1),
#'   probResponse = c(-0.1, 0, 0, 0.1)
#' )
#' scenario$nCluster
#' scenario$clusterLabelVec
#'
#' @export
simCytBuildScenario <- function(
  nMarker,
  lowMean = 0,
  highMean = 5,
  probUns = NULL,
  probResponse = NULL
) {
  stopifnot(is.numeric(nMarker), length(nMarker) == 1L, nMarker >= 1L)
  nMarker <- as.integer(nMarker)
  nCluster <- 2L^nMarker

  # Binary grid for marker combinations (0 = -, 1 = +)
  grid <- as.matrix(expand.grid(replicate(nMarker, c(0L, 1L), simplify = FALSE)))
  colnames(grid) <- paste0("F", seq_len(nMarker))

  meanExprMat <- matrix(
    ifelse(grid == 1L, highMean, lowMean),
    nrow = nCluster,
    ncol = nMarker
  )

  clusterLabelVec <- apply(grid, 1, function(row) {
    paste0(
      vapply(
        seq_along(row),
        function(i) paste0("F", i, if (row[i] == 1L) "+" else "-"),
        character(1)
      ),
      collapse = ""
    )
  })

  if (is.null(probUns)) {
    probVecUns <- rep(1 / nCluster, nCluster)
  } else {
    stopifnot(is.numeric(probUns), length(probUns) == nCluster)
    stopifnot(abs(sum(probUns) - 1) < 1e-6)
    probVecUns <- probUns
  }

  probResponseVecByStimCondition <- if (is.null(probResponse)) {
    NULL
  } else if (is.numeric(probResponse)) {
    stopifnot(length(probResponse) == nCluster)
    list(probResponse)
  } else if (is.list(probResponse)) {
    stopifnot(all(vapply(probResponse, function(x) is.numeric(x) && length(x) == nCluster, logical(1))))
    probResponse
  } else {
    stopifnot(FALSE)
  }

  list(
    nMarker = nMarker,
    nCluster = nCluster,
    meanExprMat = meanExprMat,
    clusterLabelVec = clusterLabelVec,
    probVecUns = probVecUns,
    probResponseVecByStimCondition = probResponseVecByStimCondition
  )
}

#' Construct a univariate scenario (`nMarker = 1`)
#'
#' Convenience wrapper around [simCytBuildScenario()] for two populations:
#' `F1-` and `F1+`.
#'
#' @inheritParams simCytBuildScenario
#'
#' @return A scenario list from [simCytBuildScenario()] with `nMarker = 1`.
#'
#' @examples
#' scenario <- simCytScenarioUnivariate()
#' scenario$meanExprMat
#' scenario$probResponseVecByStimCondition
#'
#' @export
simCytScenarioUnivariate <- function(
  lowMean = 0,
  highMean = 5,
  probUns = c(0.8, 0.2),
  probResponse = c(-0.1, 0.1)
) {
  simCytBuildScenario(
    nMarker = 1L,
    lowMean = lowMean,
    highMean = highMean,
    probUns = probUns,
    probResponse = probResponse
  )
}

#' Construct a bivariate scenario (`nMarker = 2`)
#'
#' Convenience wrapper around [simCytBuildScenario()] for four populations:
#' `F1-F2-`, `F1+F2-`, `F1-F2+`, and `F1+F2+`.
#'
#' @inheritParams simCytBuildScenario
#'
#' @return A scenario list from [simCytBuildScenario()] with `nMarker = 2`.
#'
#' @examples
#' scenario <- simCytScenarioBivariate()
#' scenario$clusterLabelVec
#' scenario$probVecUns
#'
#' @export
simCytScenarioBivariate <- function(
  lowMean = 0,
  highMean = 5,
  probUns = c(0.7, 0.1, 0.1, 0.1),
  probResponse = c(-0.1, 0.0, 0.0, 0.1)
) {
  simCytBuildScenario(
    nMarker = 2L,
    lowMean = lowMean,
    highMean = highMean,
    probUns = probUns,
    probResponse = probResponse
  )
}
