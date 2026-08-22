test_that("simCytExperiment creates expected structure and flowFrames", {
  set.seed(123)
  nSample <- 2
  nMarker <- 2
  nCondition <- 2
  nCluster <- 4
  nCell <- 50
  transFunc <- identity

  meanMat <- matrix(
    c(0, 0, 0, 5, 5, 0, 5, 5),
    nrow = nCluster,
    ncol = nMarker,
    byrow = TRUE
  )
  labels <- c("M1-M2-", "M1-M2+", "M1+M2-", "M1+M2+")
  probUns <- c(0.7, 0.1, 0.1, 0.1)
  probResp <- list(c(-0.2, 0.1, 0.1, 0.0))

  res <- simCytExperiment(
    nSample = nSample,
    nMarker = nMarker,
    nCondition = nCondition,
    nCluster = nCluster,
    nCellByCondition = nCell,
    transformationFunc = transFunc,
    mixtureType = "gaussianOnly",
    meanExprMat = meanMat,
    clusterLabelVec = labels,
    probVecUns = probUns,
    probExact = TRUE,
    probResponseVecByStimCondition = probResp
  )

  expect_named(res, c("flowFrameList", "labelsList"))
  expect_length(res$flowFrameList, nSample * nCondition)
  expect_length(res$labelsList, nSample * nCondition)

  expected_names <- c(
    "sample001_unstim", "sample001_stim1",
    "sample002_unstim", "sample002_stim1"
  )
  expect_named(res$flowFrameList, expected_names)
  expect_named(res$labelsList, expected_names)

  for (nm in expected_names) {
    ff <- res$flowFrameList[[nm]]
    expect_s4_class(ff, "flowFrame")
    expect_equal(nrow(flowCore::exprs(ff)), nCell)
    expect_equal(ncol(flowCore::exprs(ff)), nMarker)
    expect_equal(colnames(flowCore::exprs(ff)), c("F1", "F2"))

    lbls <- res$labelsList[[nm]]
    expect_length(lbls, nCell)
    expect_true(all(lbls %in% labels))
  }
})

test_that("simCytExperiment supports perturbations and ratio correction", {
  set.seed(456)
  nSample <- 1
  nMarker <- 1
  nCondition <- 2
  nCluster <- 2
  nCell <- 100

  # Gamma trans function attribute for ratio correction
  gammaTrans <- function(x) sqrt(abs(x))
  attr(gammaTrans, "sim_transformation") <- "gamma"

  meanMat <- matrix(c(1, 10), nrow = 2, ncol = 1)
  labels <- c("M1-", "M1+")
  probUns <- c(0.5, 0.5)

  res <- simCytExperiment(
    nSample = nSample,
    nMarker = nMarker,
    nCondition = nCondition,
    nCluster = nCluster,
    nCellByCondition = nCell,
    transformationFunc = gammaTrans,
    mixtureType = "tPlusGauss",
    meanExprMat = meanMat,
    clusterLabelVec = labels,
    probVecUns = probUns,
    probExact = FALSE,
    samplePerturbationSd = 0.1,
    conditionPerturbationSd = 0.1,
    clusterPerturbationSd = 0.1
  )

  expect_length(res$flowFrameList, 2)
  ff_unstim <- res$flowFrameList[["sample001_unstim"]]
  expect_s4_class(ff_unstim, "flowFrame")
  expect_equal(nrow(flowCore::exprs(ff_unstim)), nCell)
})

test_that("fixed-seed regression test guarantees exact cell allocation parity and reproducibility", {
  sc <- simCytScenarioUnivariate(probUns = c(0.8, 0.2), probResponse = c(-0.2, 0.2))

  set.seed(777)
  res1 <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 100,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probExact = TRUE,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition
  )

  lbls_unstim <- res1$labelsList[["sample001_unstim"]]
  lbls_stim1 <- res1$labelsList[["sample001_stim1"]]

  expect_equal(sum(lbls_unstim == "F1-"), 80)
  expect_equal(sum(lbls_unstim == "F1+"), 20)
  expect_equal(sum(lbls_stim1 == "F1-"), 60)
  expect_equal(sum(lbls_stim1 == "F1+"), 40)

  set.seed(777)
  res2 <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 100,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probExact = TRUE,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition
  )

  expect_equal(
    flowCore::exprs(res1$flowFrameList[["sample001_unstim"]]),
    flowCore::exprs(res2$flowFrameList[["sample001_unstim"]])
  )
  expect_equal(
    res1$labelsList[["sample001_unstim"]],
    res2$labelsList[["sample001_unstim"]]
  )
})

test_that("fixed-seed regression test validates mixture types and ratio corrections reproducibility", {
  gamma_trans <- simCytTransformGamma()
  mean_mat <- matrix(c(1, 2), nrow = 2, ncol = 1)
  labels <- c("F1-", "F1+")
  prob_uns <- c(0.5, 0.5)
  prob_resp <- list(c(-0.1, 0.1))

  set.seed(999)
  res_gamma1 <- simCytExperiment(
    nSample = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    nCellByCondition = 80,
    transformationFunc = gamma_trans,
    mixtureType = "tOnly",
    meanExprMat = mean_mat,
    clusterLabelVec = labels,
    probVecUns = prob_uns,
    probResponseVecByStimCondition = prob_resp,
    samplePerturbationSd = 0.05,
    conditionPerturbationSd = 0.05,
    clusterPerturbationSd = 0.05
  )

  set.seed(999)
  res_gamma2 <- simCytExperiment(
    nSample = 1,
    nMarker = 1,
    nCondition = 2,
    nCluster = 2,
    nCellByCondition = 80,
    transformationFunc = gamma_trans,
    mixtureType = "tOnly",
    meanExprMat = mean_mat,
    clusterLabelVec = labels,
    probVecUns = prob_uns,
    probResponseVecByStimCondition = prob_resp,
    samplePerturbationSd = 0.05,
    conditionPerturbationSd = 0.05,
    clusterPerturbationSd = 0.05
  )

  expect_equal(
    flowCore::exprs(res_gamma1$flowFrameList[["sample001_stim1"]]),
    flowCore::exprs(res_gamma2$flowFrameList[["sample001_stim1"]])
  )
})

test_that("probExact handles cases where rounded allocations exceed total cells without negative counts", {
  sc <- simCytScenarioBivariate(
    probUns = c(0.05, 0.31, 0.32, 0.32),
    probResponse = c(0, 0, 0, 0)
  )

  res <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 5,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE
  )

  lbls_unstim <- res$labelsList[["sample001_unstim"]]
  lbls_stim1 <- res$labelsList[["sample001_stim1"]]

  expect_length(lbls_unstim, 5L)
  expect_length(lbls_stim1, 5L)
  expect_equal(nrow(flowCore::exprs(res$flowFrameList[["sample001_unstim"]])), 5L)
  expect_equal(nrow(flowCore::exprs(res$flowFrameList[["sample001_stim1"]])), 5L)
})

test_that("simCytExperiment and simCytSample support single-cell simulations", {
  sc1 <- simCytScenarioUnivariate()
  res1 <- simCytExperiment(
    scenario = sc1,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 1,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE
  )
  expect_equal(nrow(flowCore::exprs(res1$flowFrameList[[1]])), 1L)
  expect_equal(ncol(flowCore::exprs(res1$flowFrameList[[1]])), 1L)

  sc2 <- simCytScenarioBivariate()
  res2 <- simCytExperiment(
    scenario = sc2,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 1,
    transformationFunc = simCytTransformIdentity(),
    probExact = FALSE
  )
  expect_equal(nrow(flowCore::exprs(res2$flowFrameList[[1]])), 1L)
  expect_equal(ncol(flowCore::exprs(res2$flowFrameList[[1]])), 2L)
  expect_equal(colnames(flowCore::exprs(res2$flowFrameList[[1]])), c("F1", "F2"))
})

test_that("simCytSample accepts numeric parameters without strict integer type requirements", {
  sc <- simCytScenarioBivariate()
  res <- simcyto:::simCytSample(
    nMarker = 2,
    nCondition = 2,
    nCluster = 4,
    nCellByCondition = 10,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probExact = FALSE
  )
  expect_named(res, c("flowFrameList", "conditionLabelsList"))
  expect_length(res$flowFrameList, 2)
})

test_that("stimMeanShift = 0 preserves exact backward-compatibility", {
  sc <- simCytScenarioBivariate()

  set.seed(42)
  res_default <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE
  )

  set.seed(42)
  res_zero_shift <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE,
    stimMeanShift = 0
  )

  expect_equal(
    flowCore::exprs(res_default$flowFrameList[["sample001_unstim"]]),
    flowCore::exprs(res_zero_shift$flowFrameList[["sample001_unstim"]])
  )
  expect_equal(
    flowCore::exprs(res_default$flowFrameList[["sample001_stim1"]]),
    flowCore::exprs(res_zero_shift$flowFrameList[["sample001_stim1"]])
  )
  expect_equal(res_default$labelsList, res_zero_shift$labelsList)
})

test_that("stimMeanShift applies signed additive shift to stimulated conditions while preserving SD and separation", {
  sc <- simCytScenarioBivariate()
  shifts <- c(1.5, -2.0)

  for (shift in shifts) {
    set.seed(123)
    res_base <- simCytExperiment(
      scenario = sc,
      nSample = 1,
      nCondition = 2,
      nCellByCondition = 100,
      transformationFunc = simCytTransformIdentity(),
      probExact = TRUE,
      stimMeanShift = 0
    )

    set.seed(123)
    res_shifted <- simCytExperiment(
      scenario = sc,
      nSample = 1,
      nCondition = 2,
      nCellByCondition = 100,
      transformationFunc = simCytTransformIdentity(),
      probExact = TRUE,
      stimMeanShift = shift
    )

    # 1. Unstimulated expression is completely unchanged
    expect_equal(
      flowCore::exprs(res_shifted$flowFrameList[["sample001_unstim"]]),
      flowCore::exprs(res_base$flowFrameList[["sample001_unstim"]])
    )

    # 2. Stimulated expression is exactly shifted by `shift`
    expr_base_stim <- flowCore::exprs(res_base$flowFrameList[["sample001_stim1"]])
    expr_shifted_stim <- flowCore::exprs(res_shifted$flowFrameList[["sample001_stim1"]])
    expect_equal(expr_shifted_stim, expr_base_stim + shift, tolerance = 1e-10)

    # 3. Cluster labels and ordering are identical
    lbls_base <- res_base$labelsList[["sample001_stim1"]]
    lbls_shifted <- res_shifted$labelsList[["sample001_stim1"]]
    expect_equal(lbls_shifted, lbls_base)

    # 4. Within-component SD is invariant under additive shift
    for (cluster_label in unique(lbls_base)) {
      idx <- which(lbls_base == cluster_label)
      sd_base <- apply(expr_base_stim[idx, , drop = FALSE], 2, stats::sd)
      sd_shifted <- apply(expr_shifted_stim[idx, , drop = FALSE], 2, stats::sd)
      expect_equal(sd_shifted, sd_base, tolerance = 1e-10)
    }

    # 5. Component separation within stimulated condition is invariant
    unique_clusters <- unique(lbls_base)
    if (length(unique_clusters) >= 2) {
      c1_idx <- which(lbls_base == unique_clusters[1])
      c2_idx <- which(lbls_base == unique_clusters[2])
      mean_diff_base <- colMeans(expr_base_stim[c1_idx, , drop = FALSE]) -
        colMeans(expr_base_stim[c2_idx, , drop = FALSE])
      mean_diff_shifted <- colMeans(expr_shifted_stim[c1_idx, , drop = FALSE]) -
        colMeans(expr_shifted_stim[c2_idx, , drop = FALSE])
      expect_equal(mean_diff_shifted, mean_diff_base, tolerance = 1e-10)
    }
  }
})

test_that("stimMeanShift works after non-linear transformations (gamma and skew)", {
  sc <- simCytScenarioBivariate()
  shift <- 0.75

  for (trans in list(simCytTransformGamma(), simCytTransformSkew())) {
    set.seed(456)
    res_base <- simCytExperiment(
      scenario = sc,
      nSample = 1,
      nCondition = 2,
      nCellByCondition = 80,
      transformationFunc = trans,
      probExact = TRUE,
      stimMeanShift = 0
    )

    set.seed(456)
    res_shifted <- simCytExperiment(
      scenario = sc,
      nSample = 1,
      nCondition = 2,
      nCellByCondition = 80,
      transformationFunc = trans,
      probExact = TRUE,
      stimMeanShift = shift
    )

    # Unstimulated is identical
    expect_equal(
      flowCore::exprs(res_shifted$flowFrameList[["sample001_unstim"]]),
      flowCore::exprs(res_base$flowFrameList[["sample001_unstim"]])
    )

    # Stimulated has exact post-transformation shift added
    expr_base_stim <- flowCore::exprs(res_base$flowFrameList[["sample001_stim1"]])
    expr_shifted_stim <- flowCore::exprs(res_shifted$flowFrameList[["sample001_stim1"]])
    expect_equal(expr_shifted_stim, expr_base_stim + shift, tolerance = 1e-10)
  }
})

test_that("stimMeanShift validates inputs", {
  sc <- simCytScenarioBivariate()
  expect_error(
    simCytExperiment(scenario = sc, nSample = 1, nCondition = 2, nCellByCondition = 10, transformationFunc = identity, stimMeanShift = "invalid"),
    "is.numeric\\(stimMeanShift\\)"
  )
  expect_error(
    simCytExperiment(scenario = sc, nSample = 1, nCondition = 2, nCellByCondition = 10, transformationFunc = identity, stimMeanShift = c(1, 2)),
    "length\\(stimMeanShift\\) == 1L"
  )
  expect_error(
    simCytExperiment(scenario = sc, nSample = 1, nCondition = 2, nCellByCondition = 10, transformationFunc = identity, stimMeanShift = Inf),
    "is.finite\\(stimMeanShift\\)"
  )
})

test_that("stimSdMultiplier = 1 preserves exact backward-compatibility", {
  sc <- simCytScenarioBivariate()

  set.seed(42)
  res_default <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE
  )

  set.seed(42)
  res_mult1 <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE,
    stimSdMultiplier = 1
  )

  expect_equal(
    flowCore::exprs(res_default$flowFrameList[["sample001_unstim"]]),
    flowCore::exprs(res_mult1$flowFrameList[["sample001_unstim"]])
  )
  expect_equal(
    flowCore::exprs(res_default$flowFrameList[["sample001_stim1"]]),
    flowCore::exprs(res_mult1$flowFrameList[["sample001_stim1"]])
  )
  expect_equal(res_default$labelsList, res_mult1$labelsList)
})

test_that("stimSdMultiplier scales within-component SD while preserving component means and separation", {
  sc <- simCytScenarioUnivariate()
  mults <- c(1.5, 0.6)

  for (mult in mults) {
    set.seed(789)
    res_base <- simCytExperiment(
      scenario = sc,
      nSample = 1,
      nCondition = 2,
      nCellByCondition = 100,
      transformationFunc = simCytTransformIdentity(),
      probExact = TRUE,
      stimSdMultiplier = 1
    )

    set.seed(789)
    res_scaled <- simCytExperiment(
      scenario = sc,
      nSample = 1,
      nCondition = 2,
      nCellByCondition = 100,
      transformationFunc = simCytTransformIdentity(),
      probExact = TRUE,
      stimSdMultiplier = mult
    )

    # 1. Unstimulated expression is completely unchanged
    expect_equal(
      flowCore::exprs(res_scaled$flowFrameList[["sample001_unstim"]]),
      flowCore::exprs(res_base$flowFrameList[["sample001_unstim"]])
    )

    # 2. Stimulated component SDs are scaled by mult to numerical tolerance
    expr_base_stim <- flowCore::exprs(res_base$flowFrameList[["sample001_stim1"]])
    expr_scaled_stim <- flowCore::exprs(res_scaled$flowFrameList[["sample001_stim1"]])
    lbls_base <- res_base$labelsList[["sample001_stim1"]]
    lbls_scaled <- res_scaled$labelsList[["sample001_stim1"]]
    expect_equal(lbls_scaled, lbls_base)

    for (cl in unique(lbls_base)) {
      idx <- which(lbls_base == cl)
      sd_base <- stats::sd(expr_base_stim[idx, 1])
      sd_scaled <- stats::sd(expr_scaled_stim[idx, 1])
      expect_equal(sd_scaled, sd_base * mult, tolerance = 1e-10)

      # 3. Component means are unchanged
      mean_base <- mean(expr_base_stim[idx, 1])
      mean_scaled <- mean(expr_scaled_stim[idx, 1])
      expect_equal(mean_scaled, mean_base, tolerance = 1e-10)
    }

    # 4. Component separation within stimulated condition is unchanged
    unique_clusters <- unique(lbls_base)
    if (length(unique_clusters) >= 2) {
      c1_idx <- which(lbls_base == unique_clusters[1])
      c2_idx <- which(lbls_base == unique_clusters[2])
      sep_base <- mean(expr_base_stim[c1_idx, 1]) - mean(expr_base_stim[c2_idx, 1])
      sep_scaled <- mean(expr_scaled_stim[c1_idx, 1]) - mean(expr_scaled_stim[c2_idx, 1])
      expect_equal(sep_scaled, sep_base, tolerance = 1e-10)
    }
  }
})

test_that("stimSdMultiplier scales independently around each component's own mean in bivariate scenarios", {
  sc <- simCytScenarioBivariate()
  mult <- 1.4

  set.seed(999)
  res_base <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 120,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE,
    stimSdMultiplier = 1
  )

  set.seed(999)
  res_scaled <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 120,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE,
    stimSdMultiplier = mult
  )

  expr_base_stim <- flowCore::exprs(res_base$flowFrameList[["sample001_stim1"]])
  expr_scaled_stim <- flowCore::exprs(res_scaled$flowFrameList[["sample001_stim1"]])
  lbls <- res_base$labelsList[["sample001_stim1"]]

  for (cl in unique(lbls)) {
    idx <- which(lbls == cl)
    for (m in seq_len(ncol(expr_base_stim))) {
      sd_b <- stats::sd(expr_base_stim[idx, m])
      sd_s <- stats::sd(expr_scaled_stim[idx, m])
      expect_equal(sd_s, sd_b * mult, tolerance = 1e-10)

      mean_b <- mean(expr_base_stim[idx, m])
      mean_s <- mean(expr_scaled_stim[idx, m])
      expect_equal(mean_s, mean_b, tolerance = 1e-10)
    }
  }
})

test_that("stimSdMultiplier works after non-linear transformations (gamma and skew)", {
  sc <- simCytScenarioBivariate()
  mult <- 1.3

  for (trans in list(simCytTransformGamma(), simCytTransformSkew())) {
    set.seed(333)
    res_base <- simCytExperiment(
      scenario = sc,
      nSample = 1,
      nCondition = 2,
      nCellByCondition = 80,
      transformationFunc = trans,
      probExact = TRUE,
      stimSdMultiplier = 1
    )

    set.seed(333)
    res_scaled <- simCytExperiment(
      scenario = sc,
      nSample = 1,
      nCondition = 2,
      nCellByCondition = 80,
      transformationFunc = trans,
      probExact = TRUE,
      stimSdMultiplier = mult
    )

    expect_equal(
      flowCore::exprs(res_scaled$flowFrameList[["sample001_unstim"]]),
      flowCore::exprs(res_base$flowFrameList[["sample001_unstim"]])
    )

    expr_base_stim <- flowCore::exprs(res_base$flowFrameList[["sample001_stim1"]])
    expr_scaled_stim <- flowCore::exprs(res_scaled$flowFrameList[["sample001_stim1"]])
    lbls <- res_base$labelsList[["sample001_stim1"]]

    for (cl in unique(lbls)) {
      idx <- which(lbls == cl)
      for (m in seq_len(ncol(expr_base_stim))) {
        sd_b <- stats::sd(expr_base_stim[idx, m])
        sd_s <- stats::sd(expr_scaled_stim[idx, m])
        expect_equal(sd_s, sd_b * mult, tolerance = 1e-10)

        mean_b <- mean(expr_base_stim[idx, m])
        mean_s <- mean(expr_scaled_stim[idx, m])
        expect_equal(mean_s, mean_b, tolerance = 1e-10)
      }
    }
  }
})

test_that("stimSdMultiplier and stimMeanShift combine cleanly", {
  sc <- simCytScenarioBivariate()
  mult <- 1.5
  shift <- -1.2

  set.seed(555)
  res_base <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 100,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE,
    stimSdMultiplier = 1,
    stimMeanShift = 0
  )

  set.seed(555)
  res_combined <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 100,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE,
    stimSdMultiplier = mult,
    stimMeanShift = shift
  )

  expr_base_stim <- flowCore::exprs(res_base$flowFrameList[["sample001_stim1"]])
  expr_comb_stim <- flowCore::exprs(res_combined$flowFrameList[["sample001_stim1"]])
  lbls <- res_base$labelsList[["sample001_stim1"]]

  for (cl in unique(lbls)) {
    idx <- which(lbls == cl)
    for (m in seq_len(ncol(expr_base_stim))) {
      sd_b <- stats::sd(expr_base_stim[idx, m])
      sd_c <- stats::sd(expr_comb_stim[idx, m])
      expect_equal(sd_c, sd_b * mult, tolerance = 1e-10)

      mean_b <- mean(expr_base_stim[idx, m])
      mean_c <- mean(expr_comb_stim[idx, m])
      expect_equal(mean_c, mean_b + shift, tolerance = 1e-10)
    }
  }
})

test_that("stimSdMultiplier validates inputs", {
  sc <- simCytScenarioBivariate()
  expect_error(
    simCytExperiment(scenario = sc, nSample = 1, nCondition = 2, nCellByCondition = 10, transformationFunc = identity, stimSdMultiplier = "invalid"),
    "is.numeric\\(stimSdMultiplier\\)"
  )
  expect_error(
    simCytExperiment(scenario = sc, nSample = 1, nCondition = 2, nCellByCondition = 10, transformationFunc = identity, stimSdMultiplier = -0.5),
    "stimSdMultiplier > 0"
  )
  expect_error(
    simCytExperiment(scenario = sc, nSample = 1, nCondition = 2, nCellByCondition = 10, transformationFunc = identity, stimSdMultiplier = 0),
    "stimSdMultiplier > 0"
  )
  expect_error(
    simCytExperiment(scenario = sc, nSample = 1, nCondition = 2, nCellByCondition = 10, transformationFunc = identity, stimSdMultiplier = c(1, 2)),
    "length\\(stimSdMultiplier\\) == 1L"
  )
  expect_error(
    simCytExperiment(scenario = sc, nSample = 1, nCondition = 2, nCellByCondition = 10, transformationFunc = identity, stimSdMultiplier = Inf),
    "is.finite\\(stimSdMultiplier\\)"
  )
})
