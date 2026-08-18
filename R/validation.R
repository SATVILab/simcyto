#' @title Validate inputs for simCytExperiment
#' @keywords internal
validateExperimentInputs <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nCellByCondition,
  transformationFunc
) {
  stopifnot(nSample > 0L)
  stopifnot(nMarker > 0L)
  stopifnot(nCondition > 1L)
  stopifnot(nCluster > 0L)
  stopifnot(nCluster == 2^nMarker)
  stopifnot(is.function(transformationFunc))
  stopifnot(is.numeric(nCellByCondition) || is.integer(nCellByCondition))
  stopifnot(length(nCellByCondition) %in% c(1L, nCondition))
  stopifnot(all(nCellByCondition > 0))
}

#' @title Validate inputs for simCytSample
#' @keywords internal
validateSampleInputs <- function(
  nCondition,
  nMarker,
  nCluster,
  nCellByCondition,
  transformationFunc,
  probResponseVecByStimCondition,
  probVecUns,
  probExact,
  conditionPerturbationSd,
  clusterPerturbationSd,
  meanExprMat,
  clusterLabelVec,
  covEvMin,
  covEvMax
) {
  stopifnot(is.logical(probExact))
  stopifnot(is.integer(nCondition))
  stopifnot(is.integer(nMarker))
  stopifnot(nCondition > 1L)
  stopifnot(is.integer(nCluster))
  stopifnot(nCluster > 0L)
  stopifnot(nCluster == 2^nMarker)
  stopifnot(is.function(transformationFunc))
  stopifnot(is.numeric(nCellByCondition) || is.integer(nCellByCondition))
  stopifnot(length(nCellByCondition) %in% c(1L, nCondition))
  stopifnot(all(nCellByCondition > 0))

  if (!is.null(probResponseVecByStimCondition)) {
    stopifnot(is.list(probResponseVecByStimCondition))
    stopifnot(all(sapply(probResponseVecByStimCondition, is.numeric)))
    stopifnot(length(probResponseVecByStimCondition) == (nCondition - 1L))
    stopifnot(all(sapply(probResponseVecByStimCondition, length) == nCluster))
  }

  stopifnot(is.numeric(probVecUns))
  stopifnot(length(probVecUns) == nCluster)
  stopifnot(all(probVecUns >= 0))
  stopifnot(all(probVecUns <= 1))
  stopifnot(abs(sum(probVecUns) - 1) < 1e-6)

  stopifnot(is.numeric(conditionPerturbationSd))
  stopifnot(length(conditionPerturbationSd) == 1L)
  stopifnot(conditionPerturbationSd >= 0)

  stopifnot(is.numeric(clusterPerturbationSd))
  stopifnot(length(clusterPerturbationSd) == 1L)
  stopifnot(clusterPerturbationSd >= 0)

  stopifnot(is.matrix(meanExprMat))
  stopifnot(nrow(meanExprMat) == nCluster)
  stopifnot(ncol(meanExprMat) == nMarker)
  stopifnot(!any(is.na(meanExprMat)))
  stopifnot(is.numeric(meanExprMat))

  stopifnot(is.character(clusterLabelVec))
  stopifnot(length(clusterLabelVec) == nCluster)
}
