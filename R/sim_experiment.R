#' @title Simulate a set of stimulation conditions for multiple biological samples
#'
#' @description Simulate stimulation conditions (e.g., stimulated and unstimulated) for
#' multiple biological samples, where each sample has one unstimulated condition
#' and one or more stimulated conditions. The outputs are returned as lists of `flowFrame`
#' objects representing each sample-condition combination, and matching lists of cellular
#' cluster labels.
#'
#' @param nSample Integer. Number of biological samples to simulate.
#' @param nMarker Integer. Number of markers/dimensions.
#' @param nCondition Integer. Number of conditions per sample. The first
#'   condition is unstimulated and the rest are stimulated.
#' @param nCluster Integer. Number of phenotype clusters. The current simulation
#'   contract explicitly requires `nCluster == 2^nMarker`; this is not a latent
#'   generalisation of the model.
#' @param nCellByCondition Numeric or integer vector. Number of cells per
#'   condition. If length is 1, the value is recycled across all conditions.
#' @param transformationFunc Function. Transformation applied marker-wise to
#'   simulated expression values.
#' @param mixtureType Character. Mixture distribution used for simulation.
#'   Supported values are "gaussianOnly", "tOnly", and "tPlusGauss".
#' @param meanExprMat Numeric matrix. Baseline cluster means with dimensions
#'   `nCluster x nMarker`.
#' @param clusterLabelVec Character vector. Cluster labels of length `nCluster`.
#' @param probVecUns Numeric vector. Baseline cluster probabilities for the
#'   unstimulated condition. Must have length `nCluster` and sum to 1.
#' @param probExact Logical. If TRUE, use exact probabilities for cluster assignment; if FALSE, sample from a multinomial distribution. Default is `FALSE`.
#' @param probResponseVecByStimCondition NULL or list. If provided, must be a list
#'   of length `nCondition - 1`, where each element is a numeric vector of
#'   length `nCluster`. Each vector is added to `probVecUns` to construct the
#'   stimulated-condition cluster probabilities.
#' @param samplePerturbationSd Numeric. Standard deviation of sample-level
#'   perturbations added to cluster means. Default is 0.
#' @param conditionPerturbationSd Numeric. Standard deviation of condition-level
#'   perturbations added to cluster means within each sample. Default is 0.
#' @param clusterPerturbationSd Numeric. Standard deviation of cluster-level
#'   perturbations applied during cell-level simulation. Default is 0.
#' @param covEvMin Numeric. Minimum eigenvalue for cluster covariance matrices. Default is 1.
#' @param covEvMax Numeric. Maximum eigenvalue for cluster covariance matrices. Default is 2.
#' @param scenario Optional scenario list created by a scenario builder such as
#'   `simCytScenarioBivariate()` or `simCytScenarioUnivariate()`. When supplied,
#'   any missing experiment values are filled from the scenario object while
#'   preserving explicit arguments as the authoritative overrides. The scenario
#'   contract is intentionally limited to the phenotype grid and probability
#'   fields consumed by the experiment layer.
#'
#' @return A list with two elements:
#'   - `flowFrameList`: A named list of `flowCore::flowFrame` objects.
#'   - `labelsList`: A named list of character vectors of per-cell cluster labels.
#'   Names of list elements are formatted as `sample001_unstim`, `sample001_stim1`, etc.
#'
#' @export
simCytExperiment <- function(
  nSample = NULL,
  nMarker = NULL,
  nCondition = NULL,
  nCluster = NULL,
  nCellByCondition = NULL,
  transformationFunc = NULL,
  mixtureType = "gaussianOnly",
  meanExprMat = NA,
  clusterLabelVec = NA,
  probVecUns = NULL,
  probExact = FALSE,
  probResponseVecByStimCondition = NULL,
  samplePerturbationSd = 0,
  conditionPerturbationSd = 0,
  clusterPerturbationSd = 0,
  covEvMin = 1,
  covEvMax = 2,
  scenario = NULL
) {
  if (!is.null(scenario)) {
    if (!is.list(scenario)) {
      stop("`scenario` must be a list created by a scenario builder.")
    }
    scenarioRequired <- c(
      "nMarker",
      "nCluster",
      "meanExprMat",
      "clusterLabelVec",
      "probVecUns"
    )
    missingScenarioFields <- setdiff(scenarioRequired, names(scenario))
    if (length(missingScenarioFields) > 0L) {
      stop(
        "`scenario` is missing required fields: ",
        paste(missingScenarioFields, collapse = ", ")
      )
    }

    if (is.null(nMarker)) {
      nMarker <- scenario$nMarker
    }
    if (is.null(nCluster)) {
      nCluster <- scenario$nCluster
    }
    if (is.null(probVecUns)) {
      probVecUns <- scenario$probVecUns
    }
    if (is.null(probResponseVecByStimCondition)) {
      probResponseVecByStimCondition <- scenario$probResponseVecByStimCondition
    }
    if (is.null(meanExprMat) || (is.atomic(meanExprMat) && length(meanExprMat) == 1L && is.na(meanExprMat))) {
      meanExprMat <- scenario$meanExprMat
    }
    if (is.null(clusterLabelVec) || (is.atomic(clusterLabelVec) && length(clusterLabelVec) == 1L && is.na(clusterLabelVec))) {
      clusterLabelVec <- scenario$clusterLabelVec
    }
  }

  if (is.null(nSample) || is.null(nMarker) || is.null(nCondition) || is.null(nCluster) ||
      is.null(nCellByCondition) || is.null(transformationFunc) || is.null(probVecUns)) {
    stop(
      "Missing required experiment inputs. Supply either explicit arguments or a valid `scenario` object."
    )
  }

  # Coerce inputs
  stopifnot(is.numeric(nSample))
  nSample <- as.integer(nSample)
  nMarker <- as.integer(nMarker)
  nCondition <- as.integer(nCondition)
  nCluster <- as.integer(nCluster)

  # Validate inputs using helper
  validateExperimentInputs(
    nSample,
    nMarker,
    nCondition,
    nCluster,
    nCellByCondition,
    transformationFunc
  )

  meanExprMatReference <- meanExprMat

  # Begin Simulation Logic
  nSampleXCondition <- nSample * nCondition
  sampleConditionLabelVec <- lapply(seq_len(nSample), function(currentSample) {
    if (currentSample < 10) {
      sampleName <- paste0("sample00", currentSample)
    } else if (currentSample < 100) {
      sampleName <- paste0("sample0", currentSample)
    } else {
      sampleName <- paste0("sample", currentSample)
    }
    paste0(sampleName, "_", c("unstim", paste0("stim", seq_len(nCondition - 1L))))
  })
  sampleConditionLabelVec <- unlist(sampleConditionLabelVec)

  flowFrameList <- stats::setNames(lapply(seq_len(nSampleXCondition), function(i) NULL), sampleConditionLabelVec)
  labelsList <- stats::setNames(lapply(seq_len(nSampleXCondition), function(i) NULL), sampleConditionLabelVec)

  lapply(seq_len(nSample), function(sampleInd) {
    idxLower <- (sampleInd - 1) * nCondition + 1
    meanExprMatCurrent <- if (samplePerturbationSd == 0L) {
      meanExprMat
    } else {
      meanExprMat +
        matrix(
          rep(
            stats::rnorm(nMarker, mean = 0, sd = samplePerturbationSd),
            each = nCluster
          ),
          byrow = FALSE,
          nrow = nCluster,
          ncol = nMarker
        )
    }

    outListSample <- simCytSample(
      nMarker = nMarker,
      nCondition = nCondition,
      nCluster = nCluster,
      nCellByCondition = nCellByCondition,
      transformationFunc = transformationFunc,
      mixtureType = mixtureType,
      meanExprMat = meanExprMatCurrent,
      meanExprMatReference = meanExprMatReference,
      clusterLabelVec = clusterLabelVec,
      probVecUns = probVecUns,
      probExact = probExact,
      probResponseVecByStimCondition = probResponseVecByStimCondition,
      conditionPerturbationSd = conditionPerturbationSd,
      clusterPerturbationSd = clusterPerturbationSd,
      covEvMin = covEvMin,
      covEvMax = covEvMax
    )

    for (condInd in seq_len(nCondition)) {
      flowFrameList[[idxLower + condInd - 1]] <<- outListSample$flowFrameList[[
        condInd
      ]]
      labelsList[[
        idxLower + condInd - 1
      ]] <<- outListSample$conditionLabelsList[[condInd]]
    }
    NULL
  })

  list(
    flowFrameList = flowFrameList,
    labelsList = labelsList
  )
}
