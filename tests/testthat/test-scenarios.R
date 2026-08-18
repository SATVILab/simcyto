test_that("simCytBuildScenario constructs valid scenario structures for 1D, 2D, and 3D", {
  sc1 <- simCytScenarioUnivariate()
  expect_equal(sc1$nMarker, 1L)
  expect_equal(sc1$nCluster, 2L)
  expect_equal(dim(sc1$meanExprMat), c(2, 1))
  expect_equal(sc1$clusterLabelVec, c("F1-", "F1+"))
  expect_equal(length(sc1$probVecUns), 2)
  expect_length(sc1$probResponseVecByStimCondition, 1)

  sc2 <- simCytScenarioBivariate()
  expect_equal(sc2$nMarker, 2L)
  expect_equal(sc2$nCluster, 4L)
  expect_equal(dim(sc2$meanExprMat), c(4, 2))
  expect_equal(sc2$clusterLabelVec, c("F1-F2-", "F1+F2-", "F1-F2+", "F1+F2+"))
  expect_equal(length(sc2$probVecUns), 4)

  sc3 <- simCytBuildScenario(nMarker = 3L)
  expect_equal(sc3$nMarker, 3L)
  expect_equal(sc3$nCluster, 8L)
  expect_equal(dim(sc3$meanExprMat), c(8, 3))
  expect_equal(length(sc3$probVecUns), 8)
  expect_null(sc3$probResponseVecByStimCondition)
})

test_that("scenarios can be passed directly into simCytExperiment", {
  set.seed(123)
  sc <- simCytScenarioBivariate()
  f_id <- simCytTransformIdentity()

  res <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 50,
    transformationFunc = f_id,
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition
  )

  expect_named(res, c("flowFrameList", "labelsList"))
  expect_length(res$flowFrameList, 2)
  expect_true(all(res$labelsList[["sample001_unstim"]] %in% sc$clusterLabelVec))
})

test_that("fixed-seed regression test for scenario builders ensures deterministic outputs", {
  set.seed(100)
  sc1 <- simCytScenarioUnivariate(lowMean = 1, highMean = 8, probUns = c(0.9, 0.1), probResponse = c(-0.05, 0.05))
  expect_equal(sc1$meanExprMat, matrix(c(1, 8), ncol = 1))
  expect_equal(sc1$probVecUns, c(0.9, 0.1))
  expect_equal(sc1$probResponseVecByStimCondition, list(c(-0.05, 0.05)))

  set.seed(200)
  sc2 <- simCytBuildScenario(nMarker = 2, lowMean = -1, highMean = 6)
  expect_equal(dim(sc2$meanExprMat), c(4, 2))
  expect_equal(sc2$probVecUns, rep(0.25, 4))
  expect_null(sc2$probResponseVecByStimCondition)
})
