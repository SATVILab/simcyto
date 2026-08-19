test_that("legacy StimGate exact-allocation fixture matches canonical counts and metadata", {
  fixture <- readRDS(test_path("fixtures", "stimgate_exact_allocation_fixture.rds"))
  sc <- simCytScenarioUnivariate(probUns = c(0.8, 0.2), probResponse = c(-0.2, 0.2))

  expect_equal(fixture$source_repository, "SATVILab/StimGate")
  expect_match(fixture$source_ref, "faust_manuscript_analyses@f86f1ab19a41fd2690bf180a7dcf483f9552950c")

  set.seed(fixture$seed)
  res <- simCytExperiment(
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

  expect_equal(res$labelsList[["sample001_unstim"]], fixture$labels$unstim)
  expect_equal(res$labelsList[["sample001_stim1"]], fixture$labels$stim1)
  expect_equal(sort(as.integer(table(res$labelsList[["sample001_unstim"]]))), sort(fixture$counts$unstim))
  expect_equal(sort(as.integer(table(res$labelsList[["sample001_stim1"]]))), sort(fixture$counts$stim1))
  expect_equal(
    flowCore::exprs(res$flowFrameList[["sample001_unstim"]]),
    fixture$expr$unstim,
    tolerance = 1e-12
  )
  expect_equal(
    flowCore::exprs(res$flowFrameList[["sample001_stim1"]]),
    fixture$expr$stim1,
    tolerance = 1e-12
  )
})

test_that("legacy StimGate mixed t/Gaussian + gamma fixture matches canonical expectations", {
  fixture <- readRDS(test_path("fixtures", "stimgate_mixed_tplusgauss_fixture.rds"))
  sc <- simCytScenarioBivariate()
  gamma_trans <- simCytTransformGamma(shape = 2, scale = 1)

  expect_equal(fixture$scenario, "mixed-2marker-tplusgauss-gamma")
  expect_match(fixture$source_ref, "faust_manuscript_analyses@f86f1ab19a41fd2690bf180a7dcf483f9552950c")

  set.seed(fixture$seed)
  res <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 80,
    transformationFunc = gamma_trans,
    mixtureType = "tPlusGauss",
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition,
    samplePerturbationSd = 0.05,
    conditionPerturbationSd = 0.05,
    clusterPerturbationSd = 0.05
  )

  expect_equal(res$labelsList[["sample001_unstim"]], fixture$labels$unstim)
  expect_equal(res$labelsList[["sample001_stim1"]], fixture$labels$stim1)
  expect_equal(sort(as.integer(table(res$labelsList[["sample001_unstim"]]))), sort(fixture$counts$unstim))
  expect_equal(sort(as.integer(table(res$labelsList[["sample001_stim1"]]))), sort(fixture$counts$stim1))
  expect_equal(
    flowCore::exprs(res$flowFrameList[["sample001_unstim"]]),
    fixture$expr$unstim,
    tolerance = 5e-10
  )
  expect_equal(
    flowCore::exprs(res$flowFrameList[["sample001_stim1"]]),
    fixture$expr$stim1,
    tolerance = 5e-10
  )
})
