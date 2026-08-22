#' @title Validate inputs for simCytExperiment
#' @keywords internal
validateExperimentInputs <- function(
  nSample,
  nMarker,
  nCondition,
  nCluster,
  nCellByCondition,
  transformationFunc,
  stimMeanShift = 0,
  stimSdMultiplier = 1,
  stimMeanShiftClusters = NULL,
  clusterLabelVec = NULL
) {
  stopifnot(is.numeric(nSample), length(nSample) == 1L, nSample > 0, nSample == as.integer(nSample))
  stopifnot(is.numeric(nMarker), length(nMarker) == 1L, nMarker > 0, nMarker == as.integer(nMarker))
  stopifnot(is.numeric(nCondition), length(nCondition) == 1L, nCondition > 1, nCondition == as.integer(nCondition))
  stopifnot(is.numeric(nCluster), length(nCluster) == 1L, nCluster > 0, nCluster == as.integer(nCluster))
  stopifnot(nCluster == 2^nMarker)
  stopifnot(is.function(transformationFunc))
  stopifnot(is.numeric(nCellByCondition) || is.integer(nCellByCondition))
  stopifnot(length(nCellByCondition) %in% c(1L, as.integer(nCondition)))
  stopifnot(all(nCellByCondition > 0))
  stopifnot(is.numeric(stimMeanShift), length(stimMeanShift) == 1L, is.finite(stimMeanShift))
  stopifnot(is.numeric(stimSdMultiplier), length(stimSdMultiplier) == 1L, is.finite(stimSdMultiplier), stimSdMultiplier > 0)

  if (!is.null(stimMeanShiftClusters)) {
    stopifnot(is.character(stimMeanShiftClusters) || is.numeric(stimMeanShiftClusters))
    stopifnot(length(stimMeanShiftClusters) > 0L, !any(is.na(stimMeanShiftClusters)))
    if (is.character(stimMeanShiftClusters)) {
      if (!is.null(clusterLabelVec) && !all(is.na(clusterLabelVec))) {
        stopifnot(all(stimMeanShiftClusters %in% clusterLabelVec))
      }
    } else {
      stopifnot(all(is.finite(stimMeanShiftClusters)))
      stopifnot(all(stimMeanShiftClusters == as.integer(stimMeanShiftClusters)))
      stopifnot(all(as.integer(stimMeanShiftClusters) >= 1L))
      stopifnot(all(as.integer(stimMeanShiftClusters) <= as.integer(nCluster)))
    }
  }
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
  covEvMax,
  stimMeanShift = 0,
  stimSdMultiplier = 1,
  stimMeanShiftClusters = NULL
) {
  stopifnot(is.logical(probExact), length(probExact) == 1L)
  stopifnot(is.numeric(nCondition), length(nCondition) == 1L, nCondition > 1, nCondition == as.integer(nCondition))
  stopifnot(is.numeric(nMarker), length(nMarker) == 1L, nMarker > 0, nMarker == as.integer(nMarker))
  stopifnot(is.numeric(nCluster), length(nCluster) == 1L, nCluster > 0, nCluster == as.integer(nCluster))
  stopifnot(nCluster == 2^nMarker)
  stopifnot(is.function(transformationFunc))
  stopifnot(is.numeric(nCellByCondition) || is.integer(nCellByCondition))
  stopifnot(length(nCellByCondition) %in% c(1L, as.integer(nCondition)))
  stopifnot(all(nCellByCondition > 0))

  if (!is.null(probResponseVecByStimCondition)) {
    stopifnot(is.list(probResponseVecByStimCondition))
    stopifnot(all(vapply(probResponseVecByStimCondition, is.numeric, logical(1))))
    stopifnot(length(probResponseVecByStimCondition) == (as.integer(nCondition) - 1L))
    stopifnot(all(vapply(probResponseVecByStimCondition, length, integer(1)) == as.integer(nCluster)))
  }

  stopifnot(is.numeric(probVecUns))
  stopifnot(length(probVecUns) == as.integer(nCluster))
  stopifnot(all(probVecUns >= 0))
  stopifnot(all(probVecUns <= 1))
  stopifnot(abs(sum(probVecUns) - 1) < 1e-6)

  stopifnot(is.numeric(conditionPerturbationSd))
  stopifnot(length(conditionPerturbationSd) == 1L)
  stopifnot(conditionPerturbationSd >= 0)

  stopifnot(is.numeric(clusterPerturbationSd))
  stopifnot(length(clusterPerturbationSd) == 1L)
  stopifnot(clusterPerturbationSd >= 0)

  stopifnot(is.numeric(covEvMin), length(covEvMin) == 1L, covEvMin > 0)
  stopifnot(is.numeric(covEvMax), length(covEvMax) == 1L, covEvMax >= covEvMin)

  stopifnot(is.numeric(stimMeanShift), length(stimMeanShift) == 1L, is.finite(stimMeanShift))
  stopifnot(is.numeric(stimSdMultiplier), length(stimSdMultiplier) == 1L, is.finite(stimSdMultiplier), stimSdMultiplier > 0)

  if (!is.null(stimMeanShiftClusters)) {
    stopifnot(is.character(stimMeanShiftClusters) || is.numeric(stimMeanShiftClusters))
    stopifnot(length(stimMeanShiftClusters) > 0L, !any(is.na(stimMeanShiftClusters)))
    if (is.character(stimMeanShiftClusters)) {
      stopifnot(all(stimMeanShiftClusters %in% clusterLabelVec))
    } else {
      stopifnot(all(is.finite(stimMeanShiftClusters)))
      stopifnot(all(stimMeanShiftClusters == as.integer(stimMeanShiftClusters)))
      stopifnot(all(as.integer(stimMeanShiftClusters) >= 1L))
      stopifnot(all(as.integer(stimMeanShiftClusters) <= as.integer(nCluster)))
    }
  }

  stopifnot(is.matrix(meanExprMat))
  stopifnot(nrow(meanExprMat) == as.integer(nCluster))
  stopifnot(ncol(meanExprMat) == as.integer(nMarker))
  stopifnot(!any(is.na(meanExprMat)))
  stopifnot(is.numeric(meanExprMat))

  stopifnot(is.character(clusterLabelVec))
  stopifnot(length(clusterLabelVec) == as.integer(nCluster))
}
