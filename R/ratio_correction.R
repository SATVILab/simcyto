#' Detect whether a transformation should receive upper-population ratio correction
#'
#' This is currently applied to the raw gamma transformation and to the skew
#' transformation. The lower component is never rescaled by this correction;
#' only upper components are moved and rescaled so that their observed-scale
#' separation from the pre-perturbation lower reference matches the corresponding
#' raw-scale separation after sample, condition, and cluster perturbations.
#'
#' @keywords internal
.simCytUsesUpperRatioCorrection <- function(transformationFunc) {
  transName <- attr(transformationFunc, "sim_transformation")
  isTRUE(transName %in% c("gamma", "skew")) ||
    (exists("calc_gamma", mode = "function") &&
      isTRUE(identical(transformationFunc, calc_gamma))) ||
    (exists("calc_skew", mode = "function") &&
      isTRUE(identical(transformationFunc, calc_skew)))
}

#' Backwards-compatible alias for older notebooks/scripts
#'
#' @keywords internal
.simCytUsesGammaRatioCorrection <- function(transformationFunc) {
  .simCytUsesUpperRatioCorrection(transformationFunc)
}

#' Lower reference mean after transformation
#'
#' For gamma this is the standard gamma transform of the raw lower reference.
#' For raw skew, calc_skew() is stochastic and depends on the full vector, so
#' the transformed lower reference should be estimated from the realised lower
#' population for that condition and marker.
#'
#' @keywords internal
.simCytLowerMeanTransReference <- function(
  transformationFunc,
  lowerMeanRawReference,
  yLower = NULL
) {
  transName <- attr(transformationFunc, "sim_transformation")
  isSkew <- isTRUE(identical(transName, "skew")) ||
    (exists("calc_skew", mode = "function") &&
      isTRUE(identical(transformationFunc, calc_skew)))

  if (isSkew) {
    yLower <- as.numeric(yLower)
    yLower <- yLower[is.finite(yLower)]
    if (length(yLower) > 0L) {
      return(mean(yLower))
    }
    return(NA_real_)
  }

  as.numeric(transformationFunc(lowerMeanRawReference))[1]
}

#' Solve the correction factor with combined lower and upper variance
#'
#' We adjust only the upper population. If c is the mean-distance multiplier,
#' the upper SD multiplier is 1 / c. c is chosen so that
#'
#'   c * delta / sqrt(sd_lower^2 + (sd_upper / c)^2) == target_ratio.
#'
#' @keywords internal
.simCytCombinedRatioFactor <- function(
  delta,
  sdLower,
  sdUpper,
  targetRatio,
  eps = sqrt(.Machine$double.eps)
) {
  delta <- as.numeric(delta)[1]
  sdLower <- as.numeric(sdLower)[1]
  sdUpper <- as.numeric(sdUpper)[1]
  targetRatio <- as.numeric(targetRatio)[1]

  if (
    !is.finite(delta) ||
      delta <= eps ||
      !is.finite(sdLower) ||
      sdLower < 0 ||
      !is.finite(sdUpper) ||
      sdUpper <= eps ||
      !is.finite(targetRatio) ||
      targetRatio <= eps
  ) {
    return(NA_real_)
  }

  a <- delta^2
  b <- -(targetRatio^2) * sdLower^2
  cc <- -(targetRatio^2) * sdUpper^2

  disc <- b^2 - 4 * a * cc
  if (!is.finite(disc) || disc < 0) {
    return(NA_real_)
  }

  z <- (-b + sqrt(disc)) / (2 * a)
  if (!is.finite(z) || z <= eps) {
    return(NA_real_)
  }

  sqrt(z)
}

#' Apply the upper-population ratio correction
#'
#' The target ratio is measured on the pre-transformation scale as the distance from the
#' upper mean to the lower reference mean, divided by the combined spread of the
#' lower and upper populations. The achieved ratio is measured after the
#' transformation in the same way. Only the upper population is adjusted: its mean
#' distance is multiplied by c and its SD is multiplied by 1 / c.
#'
#' @keywords internal
.simCytRatioAdjustUpper <- function(
  yUpper,
  xUpperRaw,
  xLowerRaw,
  yLower,
  lowerMeanRawReference,
  lowerMeanTransReference,
  eps = sqrt(.Machine$double.eps)
) {
  yUpper <- as.numeric(yUpper)
  xUpperRaw <- as.numeric(xUpperRaw)
  xLowerRaw <- as.numeric(xLowerRaw)
  yLower <- as.numeric(yLower)

  okUpper <- is.finite(yUpper) & is.finite(xUpperRaw)
  okLower <- is.finite(xLowerRaw) & is.finite(yLower)

  if (sum(okUpper) < 2L || sum(okLower) < 2L) {
    return(yUpper)
  }

  yUpperMean <- mean(yUpper[okUpper])
  yUpperSd <- stats::sd(yUpper[okUpper])
  xUpperMean <- mean(xUpperRaw[okUpper])
  xUpperSd <- stats::sd(xUpperRaw[okUpper])
  xLowerSd <- stats::sd(xLowerRaw[okLower])
  yLowerSd <- stats::sd(yLower[okLower])

  if (
    !is.finite(yUpperSd) ||
      yUpperSd <= eps ||
      !is.finite(xUpperSd) ||
      xUpperSd <= eps ||
      !is.finite(xLowerSd) ||
      xLowerSd < 0 ||
      !is.finite(yLowerSd) ||
      yLowerSd < 0 ||
      !is.finite(lowerMeanRawReference) ||
      !is.finite(lowerMeanTransReference)
  ) {
    return(yUpper)
  }

  targetRatio <- (xUpperMean - lowerMeanRawReference) /
    sqrt(xLowerSd^2 + xUpperSd^2)

  deltaTrans <- yUpperMean - lowerMeanTransReference

  cFactor <- .simCytCombinedRatioFactor(
    delta = deltaTrans,
    sdLower = yLowerSd,
    sdUpper = yUpperSd,
    targetRatio = targetRatio,
    eps = eps
  )

  if (!is.finite(cFactor) || cFactor <= eps) {
    return(yUpper)
  }

  yNewMean <- lowerMeanTransReference + cFactor * deltaTrans
  yNewSd <- yUpperSd / cFactor

  out <- yUpper
  out[okUpper] <- (yUpper[okUpper] - yUpperMean) / yUpperSd * yNewSd + yNewMean
  out
}

#' Apply a transformation to a simulated cluster before ratio correction
#'
#' Ratio correction is done at the condition level, because it needs the lower
#' component spread as well as the upper component spread.
#'
#' @keywords internal
.simCytTransformMatrix <- function(simData, transformationFunc) {
  simDataMat <- as.matrix(simData)
  if (is.null(dim(simDataMat))) {
    simDataMat <- matrix(as.numeric(simData), ncol = 1L)
  }
  if (isTRUE(.simCytUsesUpperRatioCorrection(transformationFunc))) {
    minSim <- min(simDataMat, na.rm = TRUE)
    if (is.finite(minSim) && minSim < 0) {
      simDataMat <- simDataMat - minSim
    }
  }
  out <- apply(simDataMat, 2, transformationFunc)
  out <- as.matrix(out)
  if (ncol(out) != ncol(simDataMat)) {
    out <- matrix(out, ncol = ncol(simDataMat))
  }
  out
}
