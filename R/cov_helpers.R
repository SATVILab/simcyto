#' Generate a random positive definite matrix.
#'
#' This helper implements the standard random positive-definite matrix
#' construction described in the R-help / R programming community literature and
#' credited upstream to Ravi Varadhan. It is retained as a small numerical
#' utility inside `simcyto` and is documented in `inst/COPYRIGHTS` as part of
#' the package's provenance trail.
#'
#' @keywords internal
.posDef <- function(n, covEvMin = 1, covEvMax = 2) {
  ev <- stats::runif(n, min = covEvMin, max = covEvMax)
  if (n == 1) {
    return(matrix(ev, 1, 1))
  }
  z <- matrix(ncol = n, stats::rnorm(n^2))
  decomp <- qr(z)
  q <- qr.Q(decomp)
  r <- qr.R(decomp)
  d <- diag(r)
  ph <- d / abs(d)
  o <- q %*% diag(ph, nrow = n)
  z <- t(o) %*% diag(ev, nrow = n) %*% o
  z
}

#' Build a covariance matrix from marker SDs and correlations
#'
#' Convenience constructor for the matrices passed to the `covMatList`
#' argument of [simCytExperiment()].
#'
#' @param sd Numeric vector of positive marker standard deviations. A single
#'   value is recycled to `nMarker` markers.
#' @param rho Numeric. Either a single correlation used for every pair of
#'   markers, or an `nMarker x nMarker` correlation matrix. Default is 0.
#' @param nMarker Integer. Number of markers. Defaults to `length(sd)`, or to
#'   `nrow(rho)` when `rho` is a matrix.
#'
#' @return A symmetric positive-definite `nMarker x nMarker` matrix.
#'
#' @examples
#' # Same per-marker variance as covEvMin = covEvMax = 1.5, now correlated.
#' simCytCovMat(sd = sqrt(1.5), rho = 0.4, nMarker = 2)
#'
#' # Different spreads for the two markers.
#' simCytCovMat(sd = c(1, 2), rho = -0.2)
#'
#' @export
simCytCovMat <- function(sd, rho = 0, nMarker = NULL) {
  stopifnot(is.numeric(sd), length(sd) >= 1L, all(is.finite(sd)), all(sd > 0))
  stopifnot(is.numeric(rho), all(is.finite(rho)))
  if (is.null(nMarker)) {
    nMarker <- if (is.matrix(rho)) nrow(rho) else length(sd)
  }
  nMarker <- as.integer(nMarker)
  stopifnot(length(nMarker) == 1L, nMarker >= 1L)
  if (length(sd) == 1L) {
    sd <- rep(sd, nMarker)
  }
  if (length(sd) != nMarker) {
    stop("`sd` must have length 1 or `nMarker`.", call. = FALSE)
  }
  corMat <- if (is.matrix(rho)) {
    rho
  } else {
    if (length(rho) != 1L) {
      stop(
        "`rho` must be a single correlation or a correlation matrix.",
        call. = FALSE
      )
    }
    m <- matrix(rho, nMarker, nMarker)
    diag(m) <- 1
    m
  }
  isCorMat <- identical(dim(corMat), c(nMarker, nMarker)) &&
    isSymmetric(unname(corMat)) &&
    all(abs(diag(corMat) - 1) <= sqrt(.Machine$double.eps)) &&
    all(abs(corMat) <= 1)
  if (!isCorMat) {
    stop(
      "`rho` must give a symmetric correlation matrix with unit diagonal.",
      call. = FALSE
    )
  }
  out <- diag(sd, nMarker) %*% corMat %*% diag(sd, nMarker)
  .simCytCheckCovMat(out, nMarker, "simCytCovMat()")
  out
}

# Strata a randomly drawn covariance matrix can vary by.
.simCytCovFactors <- c("sample", "condition", "cluster")

.simCytCheckCovMat <- function(m, nMarker, what) {
  isValid <- is.numeric(m) && identical(dim(m), c(nMarker, nMarker)) &&
    all(is.finite(m))
  if (!isValid) {
    stop(
      what, " must be a finite ", nMarker, " x ", nMarker, " numeric matrix.",
      call. = FALSE
    )
  }
  if (!isSymmetric(unname(m))) {
    stop(what, " must be symmetric.", call. = FALSE)
  }
  if (inherits(try(chol(m), silent = TRUE), "try-error")) {
    stop(what, " must be positive definite.", call. = FALSE)
  }
  invisible(m)
}

#' Normalise `covMatList` to one covariance matrix per cluster
#' @keywords internal
.simCytResolveCovMatList <- function(covMatList, nMarker, nCluster) {
  if (is.numeric(covMatList) && !is.list(covMatList)) {
    covMatList <- list(covMatList)
  }
  if (!is.list(covMatList) || !length(covMatList) %in% c(1L, nCluster)) {
    stop(
      "`covMatList` must be a matrix or a list of 1 or `nCluster` (",
      nCluster, ") matrices.",
      call. = FALSE
    )
  }
  covMatList <- lapply(seq_along(covMatList), function(i) {
    m <- covMatList[[i]]
    if (is.numeric(m) && nMarker == 1L && length(m) == 1L) {
      m <- matrix(m, 1L, 1L)
    }
    .simCytCheckCovMat(m, nMarker, paste0("`covMatList[[", i, "]]`"))
    unname(m)
  })
  if (length(covMatList) == 1L) {
    covMatList <- rep(covMatList, nCluster)
  }
  covMatList
}

#' Validate `covVaryBy`
#' @keywords internal
.simCytCheckCovVaryBy <- function(covVaryBy) {
  isValid <- is.character(covVaryBy) && !anyNA(covVaryBy) &&
    !anyDuplicated(covVaryBy) && all(covVaryBy %in% .simCytCovFactors)
  if (!isValid) {
    stop(
      "`covVaryBy` must be a subset (possibly empty) of ",
      paste0('"', .simCytCovFactors, '"', collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(covVaryBy)
}

#' Normalise a `covSd` or `covCor` specification to one entry per cluster
#'
#' Each entry is a single value (fixed) or a length-2 range `c(lo, hi)` from
#' which values are drawn uniformly.
#'
#' @return A list of `nCluster` numeric vectors of length 1 or 2.
#' @keywords internal
.simCytResolveRangeSpec <- function(spec, name, nCluster, lower, upper,
                                    lowerOpen) {
  if (is.numeric(spec) && !is.list(spec)) {
    spec <- list(spec)
  }
  if (!is.list(spec) || !length(spec) %in% c(1L, nCluster)) {
    stop(
      "`", name, "` must be a number, a range `c(lo, hi)` or a list of 1 or ",
      "`nCluster` (", nCluster, ") of these.",
      call. = FALSE
    )
  }
  spec <- lapply(seq_along(spec), function(i) {
    x <- spec[[i]]
    aboveLower <- if (lowerOpen) all(x > lower) else all(x >= lower)
    isValid <- is.numeric(x) && length(x) %in% c(1L, 2L) &&
      all(is.finite(x)) && aboveLower && all(x < upper) &&
      (length(x) == 1L || x[[1L]] <= x[[2L]])
    if (!isValid) {
      stop(
        "`", name, "` entry ", i, " must be a value or an increasing range ",
        "`c(lo, hi)` within (", lower, ", ", upper, ").",
        call. = FALSE
      )
    }
    as.numeric(x)
  })
  if (length(spec) == 1L) {
    spec <- rep(spec, nCluster)
  }
  spec
}

#' Draw one covariance matrix from marker SD and correlation specifications
#'
#' Each marker's SD is `sd` or drawn uniformly from the range `sd`; each pair
#' of markers' correlation is `cor` or drawn uniformly from the range `cor`.
#' Correlation matrices that are not positive definite are redrawn.
#'
#' @keywords internal
.simCytDrawSdCor <- function(nMarker, sd, cor, maxTry = 1000L) {
  sdVec <- if (length(sd) == 1L) {
    rep(sd, nMarker)
  } else {
    stats::runif(nMarker, sd[[1L]], sd[[2L]])
  }
  nPair <- nMarker * (nMarker - 1L) / 2L
  for (try in seq_len(maxTry)) {
    corMat <- diag(nMarker)
    if (nPair > 0L) {
      corVec <- if (length(cor) == 1L) {
        rep(cor, nPair)
      } else {
        stats::runif(nPair, cor[[1L]], cor[[2L]])
      }
      corMat[lower.tri(corMat)] <- corVec
      corMat[upper.tri(corMat)] <- t(corMat)[upper.tri(corMat)]
    }
    if (!inherits(try(chol(corMat), silent = TRUE), "try-error")) {
      return(diag(sdVec, nMarker) %*% corMat %*% diag(sdVec, nMarker))
    }
    if (length(cor) == 1L) {
      break
    }
  }
  stop(
    "`covCor` (", paste(cor, collapse = " to "), ") does not give a ",
    "positive-definite correlation matrix for ", nMarker, " markers.",
    call. = FALSE
  )
}

#' Draw a value for every sample, condition and cluster from shared strata
#'
#' One value (for example a covariance matrix or a mean shift) is drawn with
#' `drawFun(cluster)` for each combination of the factors in `covVaryBy` and
#' reused for every stratum that agrees on those factors. When `covVaryBy`
#' omits `"cluster"`, `drawFun(1)` is used for all clusters.
#'
#' @return A list indexed `[[sample]][[condition]][[cluster]]`.
#' @keywords internal
.simCytDrawByStrata <- function(nSample, nCondition, nCluster, covVaryBy,
                                drawFun) {
  nLevel <- c(sample = nSample, condition = nCondition, cluster = nCluster)
  keyGrid <- if (length(covVaryBy) == 0L) {
    data.frame(row.names = 1L)
  } else {
    do.call(expand.grid, lapply(nLevel[covVaryBy], seq_len))
  }
  keys <- if (length(covVaryBy) == 0L) {
    ""
  } else {
    apply(as.matrix(keyGrid), 1L, paste, collapse = "_")
  }
  drawn <- lapply(seq_along(keys), function(i) {
    drawFun(if ("cluster" %in% covVaryBy) keyGrid$cluster[[i]] else 1L)
  })
  keyOf <- function(idx) paste(idx[covVaryBy], collapse = "_")
  lapply(seq_len(nSample), function(s) {
    lapply(seq_len(nCondition), function(cnd) {
      lapply(seq_len(nCluster), function(cl) {
        key <- keyOf(c(sample = s, condition = cnd, cluster = cl))
        drawn[[match(key, keys)]]
      })
    })
  })
}

#' Normalise `meanPerturbation` to a list of layers
#'
#' @return A list of layers, each a list with numeric `sd` and character
#'   `varyBy`.
#' @keywords internal
.simCytResolveMeanPerturbation <- function(meanPerturbation) {
  isLayer <- function(x) is.list(x) && "sd" %in% names(x)
  if (isLayer(meanPerturbation)) {
    meanPerturbation <- list(meanPerturbation)
  }
  if (!is.list(meanPerturbation) || length(meanPerturbation) == 0L) {
    stop(
      "`meanPerturbation` must be a layer `list(sd = , varyBy = )` or a ",
      "list of such layers.",
      call. = FALSE
    )
  }
  lapply(seq_along(meanPerturbation), function(i) {
    layer <- meanPerturbation[[i]]
    isValid <- isLayer(layer) &&
      all(names(layer) %in% c("sd", "varyBy")) &&
      is.numeric(layer$sd) && length(layer$sd) == 1L &&
      is.finite(layer$sd) && layer$sd >= 0 && "varyBy" %in% names(layer)
    if (!isValid) {
      stop(
        "`meanPerturbation[[", i, "]]` must be `list(sd = , varyBy = )` ",
        "with a single non-negative `sd`.",
        call. = FALSE
      )
    }
    varyBy <- layer$varyBy
    if (is.null(varyBy)) {
      varyBy <- character(0)
    }
    tryCatch(
      .simCytCheckCovVaryBy(varyBy),
      error = function(e) {
        stop(
          "`meanPerturbation[[", i, "]]$varyBy` must be a subset (possibly ",
          "empty) of ", paste0('"', .simCytCovFactors, '"', collapse = ", "),
          ".",
          call. = FALSE
        )
      }
    )
    list(sd = layer$sd, varyBy = varyBy)
  })
}

#' Draw the summed mean shifts of all `meanPerturbation` layers
#'
#' Each layer draws one shift per marker, `N(0, sd^2)`, for each combination of
#' its `varyBy` factors.
#'
#' @return A list indexed `[[sample]][[condition]]` of `nCluster x nMarker`
#'   matrices.
#' @keywords internal
.simCytDrawMeanShiftTable <- function(layers, nSample, nCondition, nCluster,
                                      nMarker) {
  out <- lapply(seq_len(nSample), function(s) {
    lapply(seq_len(nCondition), function(cnd) {
      matrix(0, nrow = nCluster, ncol = nMarker)
    })
  })
  for (layer in layers) {
    if (layer$sd == 0) {
      next
    }
    tbl <- .simCytDrawByStrata(
      nSample, nCondition, nCluster, layer$varyBy,
      function(cl) stats::rnorm(nMarker, mean = 0, sd = layer$sd)
    )
    for (s in seq_len(nSample)) {
      for (cnd in seq_len(nCondition)) {
        out[[s]][[cnd]] <- out[[s]][[cnd]] +
          do.call(rbind, tbl[[s]][[cnd]])
      }
    }
  }
  out
}
