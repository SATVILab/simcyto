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

#' Build a baseline probability vector for a binary marker-combination scenario
#'
#' Construct a probability vector of length `2^nMarker` from per-population
#' weights. This is a convenience helper for scenarios with many markers where
#' manually writing a vector of length `2^nMarker` is error-prone.
#'
#' The function accepts an optional named numeric vector `weights` that maps
#' cluster labels (as returned by [simCytBuildScenario()]) to unnormalised
#' weights. Any cluster not named in `weights` receives weight `defaultWeight`.
#' The resulting vector is normalised to sum to 1.
#'
#' @param nMarker Integer scalar. Number of binary markers. Determines the
#'   vector length as `2^nMarker`.
#' @param weights Optional named numeric vector mapping cluster label strings to
#'   unnormalised weights. Cluster labels use the convention `"F1-F2+F3-..."`.
#'   If `NULL`, all populations receive `defaultWeight`.
#' @param defaultWeight Numeric scalar. Weight assigned to populations not
#'   present in `weights`. Default is `1`.
#'
#' @return A normalised numeric probability vector of length `2^nMarker`.
#'
#' @examples
#' # Uniform baseline for 3 markers (8 populations)
#' p <- simCytBuildProbVec(nMarker = 3)
#' length(p) == 8
#' sum(p)
#'
#' # Put most probability on the all-negative population
#' p2 <- simCytBuildProbVec(
#'   nMarker = 3,
#'   weights = c("F1-F2-F3-" = 10),
#'   defaultWeight = 1
#' )
#' p2["F1-F2-F3-"] > 0.5
#'
#' @export
simCytBuildProbVec <- function(nMarker, weights = NULL, defaultWeight = 1) {
  stopifnot(is.numeric(nMarker), length(nMarker) == 1L, nMarker >= 1L)
  nMarker <- as.integer(nMarker)
  nCluster <- 2L^nMarker

  sc <- simCytBuildScenario(nMarker = nMarker)
  labels <- sc$clusterLabelVec

  w <- rep(defaultWeight, nCluster)
  names(w) <- labels

  if (!is.null(weights)) {
    stopifnot(is.numeric(weights), !is.null(names(weights)))
    unknown <- setdiff(names(weights), labels)
    if (length(unknown) > 0L) {
      stop(
        "Unknown cluster labels in `weights`: ",
        paste(unknown, collapse = ", ")
      )
    }
    w[names(weights)] <- weights
  }

  stopifnot(all(w >= 0))
  total <- sum(w)
  stopifnot(total > 0)
  w / total
}

#' Build a response-shift vector for a binary marker-combination scenario
#'
#' Construct a response-shift vector of length `2^nMarker` for use as an
#' element of `probResponseVecByStimCondition` in [simCytExperiment()]. This
#' is a convenience helper for scenarios where only a few populations change
#' in probability under stimulation.
#'
#' Each named entry in `shifts` specifies the additive change in probability
#' for the corresponding cluster. Unspecified clusters receive `0`. The shifts
#' must sum to zero so that the resulting stimulated probability vector remains
#' a valid probability distribution (the caller is responsible for ensuring that
#' adding the shift vector to `probVecUns` yields a non-negative vector summing
#' to 1).
#'
#' @param nMarker Integer scalar. Number of binary markers.
#' @param shifts Named numeric vector mapping cluster label strings to additive
#'   probability shifts. Cluster labels use the convention `"F1-F2+F3-..."`.
#'   Must sum to zero.
#'
#' @return A numeric vector of length `2^nMarker` whose named entries give the
#'   additive probability shifts for each cluster.
#'
#' @examples
#' # Shift 10 % probability from the all-negative to the all-positive population
#' nMarker <- 3
#' sc <- simCytBuildScenario(nMarker = nMarker)
#' neg <- sc$clusterLabelVec[1] # "F1-F2-F3-"
#' pos <- sc$clusterLabelVec[8] # "F1+F2+F3+"
#' delta <- simCytBuildResponseVec(
#'   nMarker = nMarker,
#'   shifts = setNames(c(-0.1, 0.1), c(neg, pos))
#' )
#' sum(delta) # 0
#'
#' @export
simCytBuildResponseVec <- function(nMarker, shifts) {
  stopifnot(is.numeric(nMarker), length(nMarker) == 1L, nMarker >= 1L)
  stopifnot(is.numeric(shifts), !is.null(names(shifts)))
  nMarker <- as.integer(nMarker)
  nCluster <- 2L^nMarker

  sc <- simCytBuildScenario(nMarker = nMarker)
  labels <- sc$clusterLabelVec

  unknown <- setdiff(names(shifts), labels)
  if (length(unknown) > 0L) {
    stop(
      "Unknown cluster labels in `shifts`: ",
      paste(unknown, collapse = ", ")
    )
  }

  delta <- rep(0, nCluster)
  names(delta) <- labels
  delta[names(shifts)] <- shifts

  if (abs(sum(delta)) > 1e-9) {
    stop("`shifts` must sum to zero so that stimulated probabilities remain valid.")
  }

  delta
}
