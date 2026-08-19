test_that("scenario builders produce deterministic outputs with fixed seeds", {
  set.seed(2026)
  sc1 <- simCytScenarioUnivariate(lowMean = 0, highMean = 5, probUns = c(0.8, 0.2), probResponse = c(-0.1, 0.1))

  set.seed(2026)
  sc1_dup <- simCytScenarioUnivariate(lowMean = 0, highMean = 5, probUns = c(0.8, 0.2), probResponse = c(-0.1, 0.1))

  expect_equal(sc1, sc1_dup)
  expect_equal(sc1$nMarker, 1L)
  expect_equal(sc1$nCluster, 2L)
  expect_equal(sc1$clusterLabelVec, c("F1-", "F1+"))
  expect_equal(sc1$meanExprMat, matrix(c(0, 5), ncol = 1))

  set.seed(2026)
  sc2 <- simCytScenarioBivariate()
  expect_equal(sc2$nMarker, 2L)
  expect_equal(sc2$nCluster, 4L)
  expect_equal(sc2$clusterLabelVec, c("F1-F2-", "F1+F2-", "F1-F2+", "F1+F2+"))

  set.seed(2026)
  sc3 <- simCytBuildScenario(nMarker = 3L)
  expect_equal(sc3$nMarker, 3L)
  expect_equal(sc3$nCluster, 8L)
  expect_equal(dim(sc3$meanExprMat), c(8, 3))
})

test_that("sample and condition naming follows stimgate standards", {
  sc <- simCytScenarioUnivariate()
  set.seed(100)
  res <- simCytExperiment(
    nSample = 2,
    nMarker = sc$nMarker,
    nCondition = 3,
    nCluster = sc$nCluster,
    nCellByCondition = 20,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = list(c(-0.1, 0.1), c(-0.05, 0.05))
  )

  expected_names <- c(
    "sample001_unstim", "sample001_stim1", "sample001_stim2",
    "sample002_unstim", "sample002_stim1", "sample002_stim2"
  )

  expect_named(res$flowFrameList, expected_names)
  expect_named(res$labelsList, expected_names)
})

test_that("flowFrame metadata structure and channel parameters match stimgate expectations", {
  sc <- simCytScenarioBivariate()
  set.seed(101)
  res <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition
  )

  ff <- res$flowFrameList[["sample001_unstim"]]
  expect_s4_class(ff, "flowFrame")

  exprs_mat <- flowCore::exprs(ff)
  expect_true(is.matrix(exprs_mat))
  expect_equal(dim(exprs_mat), c(50, 2))
  expect_equal(unname(colnames(exprs_mat)), c("F1", "F2"))

  pdata <- Biobase::pData(flowCore::parameters(ff))
  expect_equal(rownames(pdata), c("$P1", "$P2"))
  expect_equal(as.character(pdata$name), c("F1", "F2"))
  expect_equal(as.character(pdata$desc), c("MarkerF1", "MarkerF2"))
  expect_equal(unname(pdata$range), c(8, 9))
  expect_equal(unname(pdata$minRange), c(-2.67430615, -2.96089482), tolerance = 1e-6)
  expect_equal(unname(pdata$maxRange), c(7, 8))
})

test_that("cluster cell allocations are deterministic for probExact TRUE and FALSE", {
  sc <- simCytScenarioUnivariate(probUns = c(0.75, 0.25), probResponse = c(-0.15, 0.15))

  # probExact = TRUE
  set.seed(200)
  res_exact <- simCytExperiment(
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

  lbls_unstim <- res_exact$labelsList[["sample001_unstim"]]
  lbls_stim1 <- res_exact$labelsList[["sample001_stim1"]]

  expect_equal(sum(lbls_unstim == "F1-"), 75)
  expect_equal(sum(lbls_unstim == "F1+"), 25)
  expect_equal(sum(lbls_stim1 == "F1-"), 60)
  expect_equal(sum(lbls_stim1 == "F1+"), 40)

  # probExact = FALSE (multinomial) reproducibility
  set.seed(300)
  res_multi1 <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 100,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probExact = FALSE,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition
  )

  set.seed(300)
  res_multi2 <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 100,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probExact = FALSE,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition
  )

  expect_equal(res_multi1$labelsList, res_multi2$labelsList)
  expect_equal(
    flowCore::exprs(res_multi1$flowFrameList[["sample001_unstim"]]),
    flowCore::exprs(res_multi2$flowFrameList[["sample001_unstim"]])
  )
})

test_that("mixture modes (gaussianOnly, tOnly, tPlusGauss) produce numeric parity across seeds", {
  sc <- simCytScenarioBivariate()

  for (mix in c("gaussianOnly", "tOnly", "tPlusGauss")) {
    set.seed(500)
    resA <- simCytExperiment(
      nSample = 1,
      nMarker = sc$nMarker,
      nCondition = 2,
      nCluster = sc$nCluster,
      nCellByCondition = 40,
      transformationFunc = simCytTransformIdentity(),
      mixtureType = mix,
      meanExprMat = sc$meanExprMat,
      clusterLabelVec = sc$clusterLabelVec,
      probVecUns = sc$probVecUns,
      probResponseVecByStimCondition = sc$probResponseVecByStimCondition
    )

    set.seed(500)
    resB <- simCytExperiment(
      nSample = 1,
      nMarker = sc$nMarker,
      nCondition = 2,
      nCluster = sc$nCluster,
      nCellByCondition = 40,
      transformationFunc = simCytTransformIdentity(),
      mixtureType = mix,
      meanExprMat = sc$meanExprMat,
      clusterLabelVec = sc$clusterLabelVec,
      probVecUns = sc$probVecUns,
      probResponseVecByStimCondition = sc$probResponseVecByStimCondition
    )

    expect_equal(
      flowCore::exprs(resA$flowFrameList[["sample001_unstim"]]),
      flowCore::exprs(resB$flowFrameList[["sample001_unstim"]]),
      tolerance = 1e-12
    )
    expect_equal(
      flowCore::exprs(resA$flowFrameList[["sample001_stim1"]]),
      flowCore::exprs(resB$flowFrameList[["sample001_stim1"]]),
      tolerance = 1e-12
    )
  }
})

test_that("transformations and ratio correction exhibit reproducible behaviour", {
  sc <- simCytScenarioBivariate()

  trans_gamma <- simCytTransformGamma(shape = 2, scale = 1.5)
  trans_skew <- simCytTransformSkew(shape = 2, location = 0.5, scale = 1)
  trans_gauss <- simCytTransformGaussian(sd = 0.5)

  expect_true(simcyto:::.simCytUsesUpperRatioCorrection(trans_gamma))
  expect_true(simcyto:::.simCytUsesUpperRatioCorrection(trans_skew))
  expect_false(simcyto:::.simCytUsesUpperRatioCorrection(trans_gauss))

  set.seed(600)
  res_skew1 <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 60,
    transformationFunc = trans_skew,
    mixtureType = "gaussianOnly",
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition
  )

  set.seed(600)
  res_skew2 <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 60,
    transformationFunc = trans_skew,
    mixtureType = "gaussianOnly",
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition
  )

  expect_equal(
    flowCore::exprs(res_skew1$flowFrameList[["sample001_stim1"]]),
    flowCore::exprs(res_skew2$flowFrameList[["sample001_stim1"]]),
    tolerance = 1e-12
  )
})

test_that("multi-level perturbation settings produce reproducible numerical shifts", {
  sc <- simCytScenarioBivariate()

  set.seed(700)
  res_pert1 <- simCytExperiment(
    nSample = 2,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition,
    samplePerturbationSd = 0.2,
    conditionPerturbationSd = 0.1,
    clusterPerturbationSd = 0.1
  )

  set.seed(700)
  res_pert2 <- simCytExperiment(
    nSample = 2,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition,
    samplePerturbationSd = 0.2,
    conditionPerturbationSd = 0.1,
    clusterPerturbationSd = 0.1
  )

  expect_equal(
    flowCore::exprs(res_pert1$flowFrameList[["sample002_stim1"]]),
    flowCore::exprs(res_pert2$flowFrameList[["sample002_stim1"]]),
    tolerance = 1e-12
  )

  # Check that perturbation actually changed the baseline simulation values compared to sd = 0
  set.seed(700)
  res_nopert <- simCytExperiment(
    nSample = 2,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(),
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition,
    samplePerturbationSd = 0,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0
  )

  expect_false(identical(
    flowCore::exprs(res_pert1$flowFrameList[["sample001_unstim"]]),
    flowCore::exprs(res_nopert$flowFrameList[["sample001_unstim"]])
  ))
})
