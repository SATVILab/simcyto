simCytCluster <- function(
  nMarker,
  nCell,
  meanExprVec,
  perturbationSd = 0,
  mixtureType,
  clusterNumber,
  covEvMin = 1,
  covEvMax = 2
) {
  conditionPerturbationVec <- if (perturbationSd == 0L) {
    meanExprVec
  } else {
    meanExprVec + stats::rnorm(nMarker, mean = 0, sd = perturbationSd)
  }
  currentSigma <- .posDef(nMarker, covEvMin, covEvMax)
  simCytClusterData(
    mixtureType = mixtureType,
    clusterNumber = clusterNumber,
    nCell = nCell,
    muVec = conditionPerturbationVec,
    sigmaMat = currentSigma
  )
}

#' @title Generate simulated data from mixture component
#'
#' @description Sample from multivariate normal or t distribution.
#'
#' @param mixtureType Character. "gaussianOnly", "tOnly", or "tPlusGauss".
#' @param clusterNumber Integer. Used for alternating distributions in "tPlusGauss".
#' @param nCell Integer. Number of samples.
#' @param muVec Numeric vector. Mean vector.
#' @param sigmaMat Numeric matrix. Covariance matrix.
#'
#' @return Numeric matrix of sampled data (nCell x length(muVec)).
#'
#' @keywords internal
simCytClusterData <- function(
  mixtureType,
  clusterNumber,
  nCell,
  muVec,
  sigmaMat
) {
  if (mixtureType == "tPlusGauss") {
    if ((clusterNumber %% 2) == 0) {
      MASS::mvrnorm(
        nCell,
        mu = muVec,
        Sigma = sigmaMat
      )
    } else {
      mvtnorm::rmvt(
        nCell,
        delta = muVec,
        sigma = sigmaMat,
        df = 2
      )
    }
  } else if (mixtureType == "tOnly") {
    mvtnorm::rmvt(
      nCell,
      delta = muVec,
      sigma = sigmaMat,
      df = 2
    )
  } else {
    MASS::mvrnorm(
      nCell,
      mu = muVec,
      Sigma = sigmaMat
    )
  }
}
