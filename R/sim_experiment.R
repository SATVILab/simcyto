#' @title Simulate a set of stimulation conditions for multiple biological samples
#'
#' @description
#' Simulate unstimulated/stimulated cytometry experiments across biological
#' samples. The simulator follows the hierarchy
#' `simCytExperiment()` -> `simCytSample()` -> `simCytCondition()` ->
#' `simCytCluster()`.
#'
#' The current phenotype contract is explicit: populations are all binary marker
#' combinations, so `nCluster` must equal `2^nMarker`.
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
#'   If `probExact = TRUE`, each condition allocates exactly this number of cells
#'   by rounded cluster counts. If `probExact = FALSE`, allocations are sampled
#'   from a multinomial draw.
#' @param transformationFunc Function. Post-simulation transformation applied
#'   marker-wise to simulated expression values. This is separate from
#'   `mixtureType`, which controls the raw simulation distribution.
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
#'   stimulated-condition cluster probabilities (response shifts), so each
#'   vector must sum to zero and keep every probability between 0 and 1; an
#'   error is raised otherwise. See Details for how background and response
#'   combine.
#' @param samplePerturbationSd Numeric. Standard deviation of sample-level
#'   perturbations added to cluster means. Default is 0.
#' @param conditionPerturbationSd Numeric. Standard deviation of condition-level
#'   perturbations added to cluster means within each sample. Default is 0.
#' @param clusterPerturbationSd Numeric. Standard deviation of cluster-level
#'   perturbations applied during cell-level simulation. Default is 0. It is
#'   drawn afresh for every cluster of every condition of every sample, so the
#'   conditions of a sample differ in their cluster means; use
#'   `meanPerturbation` with `varyBy = c("sample", "cluster")` for
#'   cluster-level variation between samples only.
#' @param meanPerturbation Optional general mean perturbations (raw scale): a
#'   layer `list(sd = , varyBy = )` or a list of layers. Each layer adds a shift
#'   drawn from `N(0, sd^2)`, independently for each marker, once for each
#'   combination of its `varyBy` factors (any subset of `"sample"`,
#'   `"condition"` and `"cluster"`) and shared by every stratum that agrees on
#'   them. For example, `varyBy = "sample"` matches `samplePerturbationSd`,
#'   `c("sample", "condition")` matches `conditionPerturbationSd`,
#'   `c("sample", "condition", "cluster")` matches `clusterPerturbationSd`,
#'   `c("sample", "cluster")` moves clusters between samples without differences
#'   between a sample's conditions, and `c("condition", "cluster")` gives a fixed
#'   cluster-specific condition effect shared by every sample. Layers add to
#'   each other and to the `*PerturbationSd` arguments. Default is `NULL`.
#' @param covEvMin Numeric. Minimum eigenvalue used when drawing positive-definite
#'   covariance matrices for cluster simulation. Default is 1. Ignored when
#'   `covMatList` or `covSd` is supplied.
#' @param covEvMax Numeric. Maximum eigenvalue used when drawing positive-definite
#'   covariance matrices for cluster simulation. Default is 2. Ignored when
#'   `covMatList` or `covSd` is supplied. Setting `covEvMin == covEvMax` gives the
#'   spherical matrix `covEvMin * I`, whatever `covVaryBy` is.
#' @param covMatList Optional fixed covariance matrices (on the raw,
#'   pre-transformation scale) instead of random draws: a single
#'   `nMarker x nMarker` matrix used for every cluster, or a list of `nCluster`
#'   matrices in the order of `clusterLabelVec`. Each cluster then has the same
#'   matrix in every sample and condition. Build matrices from marker SDs and
#'   correlations with [simCytCovMat()]. For `"tOnly"` and `"tPlusGauss"`
#'   mixtures this is the scale matrix of the t distribution. Default is `NULL`.
#' @param covVaryBy Character vector. For randomly drawn covariance matrices,
#'   the strata the matrix varies by: any subset of `"sample"`, `"condition"`
#'   and `"cluster"`. One matrix is drawn for each combination of the listed
#'   factors and shared by every stratum that agrees on them. For example,
#'   `"cluster"` gives each cluster its own matrix, fixed across samples and
#'   conditions; `c("sample", "cluster")` redraws per sample but keeps a
#'   sample's conditions identical; `"condition"` gives each condition its own
#'   matrix, the same in every sample and cluster; `character(0)` uses one
#'   matrix throughout. The default, all three factors, draws a new matrix for
#'   every cluster of every condition of every sample, which makes the
#'   conditions of a sample differ in covariance unless
#'   `covEvMin == covEvMax`. When `covSd` is supplied the default is instead
#'   `character(0)`, one matrix shared by every cluster, sample and condition,
#'   or `"cluster"` when `covSd` or `covCor` differ by cluster. Cannot be
#'   combined with `covMatList`.
#' @param covSd Optional marker standard deviations (raw scale) for covariance
#'   matrices built from SDs and correlations instead of eigenvalues: a single
#'   value used for every marker, or a range `c(lo, hi)` from which each
#'   marker's SD is drawn uniformly. A list of `nCluster` such entries sets each
#'   cluster separately. Default is `NULL` (eigenvalue draws).
#' @param covCor Inter-marker correlations used with `covSd`: a single value for
#'   every pair of markers, or a range `c(lo, hi)` from which each pair's
#'   correlation is drawn uniformly (redrawn until the correlation matrix is
#'   positive definite). A list of `nCluster` entries sets each cluster
#'   separately. Default is `NULL`, meaning 0 (uncorrelated markers). Requires
#'   `covSd`. Draws are made per `covVaryBy` stratum: by default
#'   `covCor = c(0, 0.5)` draws one correlation structure shared by every
#'   cluster, sample and condition, while `covVaryBy = "cluster"` draws one per
#'   cluster. Per-cluster lists require `"cluster"` in `covVaryBy`.
#' @param stimMeanShift Numeric scalar. Deterministic signed location shift
#'   added to all markers and components in stimulated conditions on the final
#'   transformed expression scale. Default is 0.
#' @param stimSdMultiplier Numeric scalar. Multiplier applied to stimulated
#'   conditions post-transformation standard deviation around component means.
#'   Default is 1.
#' @param stimMeanShiftClusters Optional character or integer vector. Phenotype
#'   cluster labels (from `clusterLabelVec` or scenario) or 1-based cluster
#'   indices to which `stimMeanShift` should be applied. Default is `NULL`
#'   (applies to all clusters in stimulated conditions).
#' @param stimSdMultiplierClusters Optional character or integer vector. Phenotype
#'   cluster labels (from `clusterLabelVec` or scenario) or 1-based cluster
#'   indices to which `stimSdMultiplier` should be applied. Default is `NULL`
#'   (applies to all clusters in stimulated conditions).
#' @param scenario Optional scenario list created by a scenario builder such as
#'   `simCytScenarioBivariate()` or `simCytScenarioUnivariate()`. When supplied,
#'   any missing experiment values are filled from the scenario object while
#'   preserving explicit arguments as the authoritative overrides. The scenario
#'   contract is intentionally limited to the phenotype grid and probability
#'   fields consumed by the experiment layer.
#'
#' @details
#' **Background and response.** Each sample has one unstimulated condition,
#' with cluster probabilities `probVecUns`. Every stimulated condition `k` has
#' probabilities `probVecUns + probResponseVecByStimCondition[[k]]`. All
#' stimulated conditions of a sample therefore share the unstimulated
#' background, with their own response added on top: the net (stimulated minus
#' unstimulated) frequency of a cluster is exactly its response entry. Cell
#' labels name the cluster only, so background and responding cells in a
#' stimulated condition cannot be told apart. Each condition's cluster counts
#' are drawn independently (or rounded when `probExact = TRUE`).
#'
#' @return A list with two elements:
#'   - `flowFrameList`: A named list of `flowCore::flowFrame` objects.
#'   - `labelsList`: A named list of character vectors of per-cell cluster labels.
#'   Names of list elements are formatted as `sample001_unstim`, `sample001_stim1`, etc.
#'
#'   Expression values can be accessed with `flowCore::exprs()` and per-cell truth
#'   labels from `labelsList[[sample_condition_name]]`.
#'
#' @examples
#' scenario <- simCytScenarioBivariate()
#' set.seed(123)
#' sim_res <- simCytExperiment(
#'   scenario = scenario,
#'   nSample = 2,
#'   nCondition = 2,
#'   nCellByCondition = c(100, 120),
#'   transformationFunc = simCytTransformIdentity(),
#'   mixtureType = "gaussianOnly",
#'   probExact = FALSE,
#'   stimMeanShift = 0.5,
#'   stimSdMultiplier = 1.2,
#'   stimMeanShiftClusters = "F1-F2-",
#'   stimSdMultiplierClusters = "F1-F2-"
#' )
#'
#' names(sim_res$flowFrameList)
#' expr_mat <- flowCore::exprs(sim_res$flowFrameList[[1]])
#' head(expr_mat)
#' table(sim_res$labelsList[[1]])
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
  meanPerturbation = NULL,
  covEvMin = 1,
  covEvMax = 2,
  covMatList = NULL,
  covVaryBy = c("sample", "condition", "cluster"),
  covSd = NULL,
  covCor = NULL,
  stimMeanShift = 0,
  stimSdMultiplier = 1,
  stimMeanShiftClusters = NULL,
  stimSdMultiplierClusters = NULL,
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
    transformationFunc,
    stimMeanShift,
    stimSdMultiplier,
    stimMeanShiftClusters,
    stimSdMultiplierClusters,
    clusterLabelVec
  )

  meanExprMatReference <- meanExprMat

  # Covariance matrices indexed [[sample]][[condition]][[cluster]]; NULL keeps
  # the original per-cluster draw inside simCytCluster().
  sdCorMode <- !is.null(covSd) || !is.null(covCor)
  if (sdCorMode) {
    if (is.null(covSd)) {
      stop("`covCor` requires `covSd`.", call. = FALSE)
    }
    if (!is.null(covMatList)) {
      stop("`covSd` and `covCor` cannot be combined with `covMatList`.",
        call. = FALSE
      )
    }
    covVaryByDefault <- missing(covVaryBy)
  }
  .simCytCheckCovVaryBy(covVaryBy)
  stopifnot(is.numeric(covEvMin), length(covEvMin) == 1L, covEvMin > 0)
  stopifnot(is.numeric(covEvMax), length(covEvMax) == 1L, covEvMax >= covEvMin)
  covTable <- if (!is.null(covMatList)) {
    if (!missing(covVaryBy)) {
      stop(
        "`covVaryBy` applies only to randomly drawn covariance matrices; ",
        "it cannot be combined with `covMatList`.",
        call. = FALSE
      )
    }
    covMatList <- .simCytResolveCovMatList(covMatList, nMarker, nCluster)
    rep(list(rep(list(covMatList), nCondition)), nSample)
  } else if (sdCorMode) {
    sdSpec <- .simCytResolveRangeSpec(
      covSd, "covSd", nCluster,
      lower = 0, upper = Inf, lowerOpen = TRUE
    )
    corSpec <- .simCytResolveRangeSpec(
      if (is.null(covCor)) 0 else covCor, "covCor", nCluster,
      lower = -1, upper = 1, lowerOpen = TRUE
    )
    perCluster <- length(unique(sdSpec)) > 1L || length(unique(corSpec)) > 1L
    # By default one matrix is shared by every stratum; per-cluster
    # specifications need one draw per cluster.
    if (covVaryByDefault) {
      covVaryBy <- if (perCluster) "cluster" else character(0)
    }
    if (perCluster && !"cluster" %in% covVaryBy) {
      stop(
        "Per-cluster `covSd` or `covCor` requires \"cluster\" in `covVaryBy`.",
        call. = FALSE
      )
    }
    .simCytDrawByStrata(
      nSample, nCondition, nCluster, covVaryBy,
      function(cl) .simCytDrawSdCor(nMarker, sdSpec[[cl]], corSpec[[cl]])
    )
  } else if (setequal(covVaryBy, .simCytCovFactors)) {
    NULL
  } else {
    .simCytDrawByStrata(
      nSample, nCondition, nCluster, covVaryBy,
      function(cl) .posDef(nMarker, covEvMin, covEvMax)
    )
  }

  # Summed meanPerturbation shifts indexed [[sample]][[condition]]. Drawn only
  # when supplied, so the random-number stream is otherwise unchanged.
  meanShiftTable <- if (is.null(meanPerturbation)) {
    NULL
  } else {
    .simCytDrawMeanShiftTable(
      .simCytResolveMeanPerturbation(meanPerturbation),
      nSample, nCondition, nCluster, nMarker
    )
  }

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
      meanShiftByCondition = meanShiftTable[[sampleInd]],
      covEvMin = covEvMin,
      covEvMax = covEvMax,
      covMatByCondition = covTable[[sampleInd]],
      stimMeanShift = stimMeanShift,
      stimSdMultiplier = stimSdMultiplier,
      stimMeanShiftClusters = stimMeanShiftClusters,
      stimSdMultiplierClusters = stimSdMultiplierClusters
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
