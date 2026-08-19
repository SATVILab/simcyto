#' @title Simulate all cytokine combinations for a set of stimulation conditions for a single biological sample
#'
#' @description Simulate a set of stimulation conditions (e.g., stimulated and unstimulated) for
#' a single biological sample, where the first condition is always the unstimulated condition.
#'
#' @param nMarker Integer. Number of markers/dimensions.
#' @param nCondition Integer. Number of stimulation conditions (must be >= 2).
#' @param nCluster Integer. Number of clusters (must be a power of 2 between 2 and 1024).
#' @param nCellByCondition Integer or numeric vector. Number of cells per condition. If a single
#'   value, it is recycled for all conditions.
#' @param transformationFunc Function. Transformation to apply to simulated data (e.g., identity or
#'   gamma transformation).
#' @param mixtureType Character. Type of mixture distribution: "gaussianOnly", "tOnly", or
#'   "tPlusGauss".
#' @param meanExprMat Numeric matrix. Cluster mean vectors (nCluster x nMarker).
#' @param clusterLabelVec Character vector. Labels for each cluster (length nCluster).
#' @param probVecUns Numeric vector. Probability distribution for unstimulated condition
#'   (length nCluster, sums to 1).
#' @param probExact Logical. If TRUE, use exact probabilities for cluster assignment; if FALSE, sample from a multinomial distribution.
#' Default is `FALSE`.
#' @param probResponseVecByStimCondition List. Probability response vectors for each stimulated
#'   condition (length nCondition - 1, each of length nCluster).
#' @param conditionPerturbationSd Numeric. Standard deviation of condition-level perturbations
#'   to cluster means.
#' @param clusterPerturbationSd Numeric. Standard deviation of cluster-level perturbations
#'   within each condition.
#' @param covEvMin Numeric. Minimum eigenvalue for cluster covariance matrices. Default is 1.
#' @param covEvMax Numeric. Maximum eigenvalue for cluster covariance matrices. Default is 2.
#' @param meanExprMatReference Numeric matrix. Reference cluster mean matrix.
#'
#' @return A list with `flowFrameList` and `conditionLabelsList`.
#'
#' @keywords internal
simCytSample <- function(
  nMarker,
  nCondition,
  nCluster,
  nCellByCondition,
  transformationFunc,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVecUns,
  probExact,
  probResponseVecByStimCondition = NULL,
  conditionPerturbationSd = 0,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2,
  meanExprMatReference = NULL
) {
  # Validate inputs using helper
  validateSampleInputs(
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
  )

  # Begin Simulation Logic
  nCellByCondition <- if (length(nCellByCondition) == 1L) {
    rep(nCellByCondition, nCondition)
  } else {
    nCellByCondition
  }

  probVecByCondition <- if (is.null(probResponseVecByStimCondition)) {
    lapply(seq_len(nCondition), function(i) probVecUns)
  } else {
    c(list(probVecUns), lapply(probResponseVecByStimCondition, function(probResponseVec) {
      probVecUns + probResponseVec
    }))
  }

  conditionLabelVec <- c("unstim", paste0("stim", seq_len(nCondition - 1L)))
  flowList <- lapply(seq_len(nCondition), function(i) NULL) |>
    stats::setNames(conditionLabelVec)
  labelsList <- lapply(seq_len(nCondition), function(i) NULL) |>
    stats::setNames(conditionLabelVec)

  lapply(seq_len(nCondition), function(i) {
    meanExprMat <- if (conditionPerturbationSd == 0L) {
      meanExprMat
    } else {
      meanExprMat +
        matrix(
          rep(
            stats::rnorm(nMarker, mean = 0, sd = conditionPerturbationSd),
            each = nCluster
          ),
          byrow = FALSE,
          nrow = nCluster,
          ncol = nMarker
        )
    }

    outListCondition <- simCytCondition(
      nMarker = nMarker,
      nCell = nCellByCondition[[i]],
      transformationFunc = transformationFunc,
      mixtureType = mixtureType,
      meanExprMat = meanExprMat,
      meanExprMatReference = meanExprMatReference,
      clusterLabelVec = clusterLabelVec,
      probVec = probVecByCondition[[i]],
      probExact = probExact,
      clusterPerturbationSd = clusterPerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    # Create annotated data frame
    paramMeta <- data.frame(
      name = paste0("F", seq_len(nMarker)),
      desc = paste0("MarkerF", seq_len(nMarker)),
      range = apply(outListCondition$conditionMatrix, 2, max),
      minRange = apply(outListCondition$conditionMatrix, 2, min),
      maxRange = apply(outListCondition$conditionMatrix, 2, max)
    )
    rownames(paramMeta) <- paste0("$P", seq_len(nMarker))

    paramAnnotated <- Biobase::AnnotatedDataFrame(
      data = paramMeta,
      varMetadata = data.frame(
        labelDescription = c(
          "Name of instrument channel",
          "Actual marker description",
          "Range of values",
          "Minimum binary value",
          "Maximum binary value"
        ),
        row.names = c("name", "desc", "range", "minRange", "maxRange")
      )
    )

    exprMat <- outListCondition$conditionMatrix
    colnames(exprMat) <- paste0("F", seq_len(nMarker))
    ff <- flowCore::flowFrame(
      exprs = exprMat,
      parameters = paramAnnotated
    )
    pDataFF <- Biobase::pData(flowCore::parameters(ff))
    rownames(pDataFF) <- paste0("$P", seq_len(nMarker))
    pDataFF$name <- paste0("F", seq_len(nMarker))
    pDataFF$desc <- paste0("MarkerF", seq_len(nMarker))
    pDataFF$range <- stats::setNames(apply(exprMat, 2, max), paste0("F", seq_len(nMarker)))
    pDataFF$minRange <- stats::setNames(apply(exprMat, 2, min), paste0("F", seq_len(nMarker)))
    pDataFF$maxRange <- stats::setNames(apply(exprMat, 2, max), paste0("F", seq_len(nMarker)))
    Biobase::pData(flowCore::parameters(ff)) <- pDataFF
    colnames(flowCore::exprs(ff)) <- paste0("F", seq_len(nMarker))
    flowList[[i]] <<- ff
    labelsList[[i]] <<- outListCondition$conditionLabels
    NULL
  })

  list(
    "flowFrameList" = flowList,
    "conditionLabelsList" = labelsList
  )
}
