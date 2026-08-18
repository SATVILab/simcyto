#' Build a simulation scenario for n markers
#'
#' Parameterised builder for constructing cluster mean matrices, cluster labels,
#' unstimulated probabilities, and response probability vectors for n markers.
#'
#' @param nMarker Integer. Number of markers. Must be >= 1.
#' @param lowMean Numeric. Baseline mean for negative marker expressions. Default is 0.
#' @param highMean Numeric. Baseline mean for positive marker expressions. Default is 5.
#' @param probUns Numeric vector. Baseline cluster probabilities for the unstimulated condition
#'   (length 2^nMarker, summing to 1). If NULL, uniform probabilities are assigned.
#' @param probResponse List or numeric vector. Response probabilities for stimulated conditions.
#'   If a numeric vector of length 2^nMarker, it is wrapped in a list for 1 stimulated condition.
#'   Default is NULL.
#'
#' @return A list with elements:
#'   - `nMarker`: Number of markers.
#'   - `nCluster`: Number of clusters (2^nMarker).
#'   - `meanExprMat`: Matrix of dimension `nCluster x nMarker`.
#'   - `clusterLabelVec`: Character vector of cluster labels of length `nCluster`.
#'   - `probVecUns`: Numeric vector of unstimulated cluster probabilities.
#'   - `probResponseVecByStimCondition`: List of response probability vectors or NULL.
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

#' Construct a univariate simulation scenario (1 marker)
#'
#' Helper for 1-marker simulation scenarios.
#'
#' @param lowMean Numeric. Mean for F1- cluster. Default is 0.
#' @param highMean Numeric. Mean for F1+ cluster. Default is 5.
#' @param probUns Numeric vector of length 2. Unstimulated cluster probabilities. Default is `c(0.8, 0.2)`.
#' @param probResponse List or numeric vector of length 2. Response vector for stimulated conditions. Default is `c(-0.1, 0.1)`.
#'
#' @return A scenario list created by `simCytBuildScenario`.
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

#' Construct a bivariate simulation scenario (2 markers)
#'
#' Helper for 2-marker simulation scenarios.
#'
#' @param lowMean Numeric. Mean for negative marker expressions. Default is 0.
#' @param highMean Numeric. Mean for positive marker expressions. Default is 5.
#' @param probUns Numeric vector of length 4. Unstimulated cluster probabilities. Default is `c(0.7, 0.1, 0.1, 0.1)`.
#' @param probResponse List or numeric vector of length 4. Response vector for stimulated conditions. Default is `c(-0.1, 0.0, 0.0, 0.1)`.
#'
#' @return A scenario list created by `simCytBuildScenario`.
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
