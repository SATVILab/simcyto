#' @title Simulate cytometric data for a single stimulation condition
#'
#' @description Simulate flow cytometric data for a single stimulation condition.
#'
#' @param nMarker Integer. Number of markers/dimensions.
#' @param nCell Integer. Total number of cells to simulate.
#' @param transformationFunc Function. Post-simulation transformation to apply
#'   to simulated data.
#' @param mixtureType Character. Type of mixture distribution.
#' @param meanExprMat Numeric matrix. Cluster mean vectors.
#' @param clusterLabelVec Character vector. Cluster labels.
#' @param probVec Numeric vector. Cluster probabilities.
#' @param probExact Logical. If TRUE, use exact probabilities for cluster assignment; if FALSE, sample from a multinomial distribution.
#' Default is `FALSE`.
#' @param clusterPerturbationSd Numeric. Cluster-level perturbation SD.
#' @param covEvMin Numeric. Minimum eigenvalue for cluster covariance matrices. Default is 1.
#' @param covEvMax Numeric. Maximum eigenvalue for cluster covariance matrices. Default is 2.
#' @param meanExprMatReference Numeric matrix. Reference cluster mean matrix.
#'
#' @return A list with `conditionMatrix` and `conditionLabels`.
#'
#' @keywords internal
simCytCondition <- function(
  nMarker,
  nCell,
  transformationFunc,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVec,
  probExact = FALSE,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2,
  meanExprMatReference = NULL
) {
  numClusters <- nrow(meanExprMat)

  if (is.null(meanExprMatReference)) {
    meanExprMatReference <- meanExprMat
  }
  stopifnot(is.matrix(meanExprMatReference))
  stopifnot(all(dim(meanExprMatReference) == dim(meanExprMat)))

  ratioCorrection <- .simCytUsesUpperRatioCorrection(transformationFunc)
  lowerMeanExprVecReference <- apply(meanExprMatReference, 2, min)

  nCellVec <- if (probExact) {
    initAllocVec <- round(nCell * probVec)
    nAlloc <- sum(initAllocVec)
    if (nAlloc != nCell) {
      diffAlloc <- nCell - nAlloc
      if (diffAlloc > 0) {
        probOrder <- order(probVec, decreasing = TRUE)
        for (i in seq_len(diffAlloc)) {
          idx <- probOrder[((i - 1L) %% length(probOrder)) + 1L]
          initAllocVec[idx] <- initAllocVec[idx] + 1L
        }
      } else {
        for (i in seq_len(-diffAlloc)) {
          candidates <- which(initAllocVec > 0L)
          bestCandidate <- candidates[order(probVec[candidates], decreasing = FALSE)[1L]]
          initAllocVec[bestCandidate] <- initAllocVec[bestCandidate] - 1L
        }
      }
    }
    initAllocVec
  } else {
    as.vector(t(stats::rmultinom(1, nCell, probVec)))
  }
  nCellVecCum <- cumsum(nCellVec)

  nCellVecObserved <- nCellVec[nCellVec > 0L]
  clusterLabelVecObserved <- clusterLabelVec[nCellVec > 0L]
  cellLabelVec <- unlist(lapply(seq_along(nCellVecObserved), function(i) {
    rep(clusterLabelVecObserved[i], nCellVecObserved[i])
  }))

  rawDataList <- vector("list", numClusters)
  transformedDataList <- vector("list", numClusters)
  outDataIndClusterList <- vector("list", numClusters)

  for (clusterNumber in seq_len(numClusters)) {
    nCellCluster <- nCellVec[[clusterNumber]]
    if (nCellCluster == 0L) {
      next
    }
    outDataIndClusterLower <- if (clusterNumber == 1L) {
      1L
    } else {
      nCellVecCum[clusterNumber - 1] + 1
    }
    outDataIndClusterUpper <- nCellVecCum[[clusterNumber]]
    outDataIndClusterVec <- seq.int(
      outDataIndClusterLower,
      outDataIndClusterUpper
    )
    meanExprVec <- as.numeric(meanExprMat[clusterNumber, , drop = TRUE])
    simData <- simCytCluster(
      nMarker = nMarker,
      nCell = nCellCluster,
      meanExprVec = meanExprVec,
      perturbationSd = clusterPerturbationSd,
      mixtureType = mixtureType,
      clusterNumber = clusterNumber,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )
    simData <- as.matrix(simData)
    if (ncol(simData) != nMarker) {
      simData <- matrix(simData, ncol = nMarker)
    }

    rawDataList[[clusterNumber]] <- simData
    transformedDataList[[clusterNumber]] <- .simCytTransformMatrix(
      simData = simData,
      transformationFunc = transformationFunc
    )
    outDataIndClusterList[[clusterNumber]] <- outDataIndClusterVec
  }

  if (isTRUE(ratioCorrection)) {
    for (markerInd in seq_len(nMarker)) {
      lowerMeanRawReference <- lowerMeanExprVecReference[[markerInd]]

      lowerClusterInd <- which(
        meanExprMatReference[, markerInd] <=
          lowerMeanRawReference + sqrt(.Machine$double.eps)
      )
      lowerClusterInd <- lowerClusterInd[
        vapply(
          lowerClusterInd,
          function(i) {
            !is.null(rawDataList[[i]]) && nrow(rawDataList[[i]]) > 0L
          },
          logical(1)
        )
      ]

      if (length(lowerClusterInd) == 0L) {
        next
      }

      xLowerRaw <- unlist(lapply(
        lowerClusterInd,
        function(i) rawDataList[[i]][, markerInd, drop = TRUE]
      ))
      yLower <- unlist(lapply(
        lowerClusterInd,
        function(i) transformedDataList[[i]][, markerInd, drop = TRUE]
      ))

      lowerMeanTransReference <- .simCytLowerMeanTransReference(
        transformationFunc = transformationFunc,
        lowerMeanRawReference = lowerMeanRawReference,
        yLower = yLower
      )
      if (!is.finite(lowerMeanTransReference)) {
        next
      }

      for (clusterNumber in seq_len(numClusters)) {
        if (is.null(rawDataList[[clusterNumber]])) {
          next
        }
        isUpperPopulation <- meanExprMatReference[clusterNumber, markerInd] >
          lowerMeanRawReference + sqrt(.Machine$double.eps)
        if (!isTRUE(isUpperPopulation)) {
          next
        }

        transformedDataList[[clusterNumber]][, markerInd] <-
          .simCytRatioAdjustUpper(
            yUpper = transformedDataList[[clusterNumber]][, markerInd],
            xUpperRaw = rawDataList[[clusterNumber]][, markerInd],
            xLowerRaw = xLowerRaw,
            yLower = yLower,
            lowerMeanRawReference = lowerMeanRawReference,
            lowerMeanTransReference = lowerMeanTransReference
          )
      }
    }
  }

  outData <- if (nMarker == 1L) {
    rep(NA_real_, nCell)
  } else {
    matrix(NA_real_, nrow = nCell, ncol = nMarker)
  }

  for (clusterNumber in seq_len(numClusters)) {
    if (is.null(transformedDataList[[clusterNumber]])) {
      next
    }
    outDataIndClusterVec <- outDataIndClusterList[[clusterNumber]]
    if (nMarker == 1L) {
      outData[outDataIndClusterVec] <- transformedDataList[[clusterNumber]][,
        1L
      ]
    } else {
      outData[outDataIndClusterVec, ] <- transformedDataList[[clusterNumber]]
    }
  }

if (length(outData) > 0L && any(!is.finite(outData))) {
  nonFiniteCount <- sum(!is.finite(outData))
  warning(
    sprintf(
      "Replaced %d non-finite values with 0 before exporting the condition matrix.",
      nonFiniteCount
    ),
    call. = FALSE
  )
  outData[!is.finite(outData)] <- 0
}

reorderVec <- sample.int(nCell)
  if (nMarker == 1L) {
    outData <- outData[reorderVec]
    outData <- matrix(outData, ncol = 1)
  } else {
    outData <- outData[reorderVec, , drop = FALSE]
  }
  colnames(outData) <- paste0("F", seq_len(nMarker))
  cellLabelVec <- cellLabelVec[reorderVec]
  list(
    conditionMatrix = outData,
    conditionLabels = cellLabelVec
  )
}
