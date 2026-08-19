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

  probVecByCondition <- vector("list", nCondition)
  probVecByCondition[[1L]] <- probVecUns

  if (!is.null(probResponseVecByStimCondition)) {
    for (stimInd in seq_len(nCondition - 1L)) {
      probVecByCondition[[stimInd + 1L]] <- probVecUns + probResponseVecByStimCondition[[stimInd]]
    }
  } else {
    for (stimInd in seq_len(nCondition - 1L)) {
      probVecByCondition[[stimInd + 1L]] <- probVecUns
    }
  }

  conditionLabelVec <- c("unstim", paste0("stim", seq_len(nCondition - 1L)))
  flowList <- stats::setNames(lapply(seq_len(nCondition), function(i) NULL), conditionLabelVec)
  labelsList <- stats::setNames(lapply(seq_len(nCondition), function(i) NULL), conditionLabelVec)

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
    paramAnnotated@data$name <- as.character(paramAnnotated@data$name)
    names(paramAnnotated@data$name) <- NULL
    paramAnnotated@data$desc <- as.character(paramAnnotated@data$desc)
    names(paramAnnotated@data$desc) <- NULL

    exprMat <- outListCondition$conditionMatrix
    colnames(exprMat) <- paste0("F", seq_len(nMarker))
ff <- flowCore::flowFrame(
  exprs = exprMat,
  parameters = paramAnnotated
)
flowCore::exprs(ff) <- exprMat
pDataFF <- Biobase::pData(flowCore::parameters(ff))
rownames(pDataFF) <- paste0("$P", seq_len(nMarker))
pDataFF$name <- paste0("F", seq_len(nMarker))
pDataFF$desc <- paste0("MarkerF", seq_len(nMarker))
colMaxima <- apply(exprMat, 2, max)
colMinima <- apply(exprMat, 2, min)
pDataFF$range <- colMaxima
pDataFF$minRange <- colMinima
pDataFF$maxRange <- colMaxima
Biobase::pData(flowCore::parameters(ff)) <- pDataFF
    flowList[[i]] <<- ff
    labelsList[[i]] <<- outListCondition$conditionLabels
    NULL
  })

  list(
    "flowFrameList" = flowList,
    "conditionLabelsList" = labelsList
  )
}
