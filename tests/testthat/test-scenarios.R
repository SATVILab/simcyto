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

test_that("simCytBuildScenario supports six markers / 64 binary populations", {
  sc6 <- simCytBuildScenario(nMarker = 6L, lowMean = 0, highMean = 5)
  expect_equal(sc6$nMarker, 6L)
  expect_equal(sc6$nCluster, 64L)
  expect_equal(dim(sc6$meanExprMat), c(64L, 6L))
  expect_length(sc6$clusterLabelVec, 64L)
  expect_equal(sc6$clusterLabelVec[1], "F1-F2-F3-F4-F5-F6-")
  expect_equal(sc6$clusterLabelVec[64], "F1+F2+F3+F4+F5+F6+")
  expect_equal(sum(sc6$probVecUns), 1)
  expect_null(sc6$probResponseVecByStimCondition)
})

test_that("simCytBuildScenario six-marker scenario runs end-to-end in simCytExperiment", {
  sc6 <- simCytBuildScenario(nMarker = 6L, lowMean = 0, highMean = 5)

  neg_label <- sc6$clusterLabelVec[1]
  pos_label <- sc6$clusterLabelVec[64]
  probUns <- simCytBuildProbVec(
    nMarker = 6L,
    weights = setNames(c(50), neg_label),
    defaultWeight = 1
  )
  response <- simCytBuildResponseVec(
    nMarker = 6L,
    shifts = setNames(c(-0.1, 0.1), c(neg_label, pos_label))
  )

  sc6_custom <- simCytBuildScenario(
    nMarker = 6L,
    lowMean = 0,
    highMean = 5,
    probUns = probUns,
    probResponse = response
  )

  set.seed(42)
  res <- simCytExperiment(
    scenario = sc6_custom,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 64,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE
  )

  expect_named(res, c("flowFrameList", "labelsList"))
  expect_length(res$flowFrameList, 2)
  expect_equal(ncol(flowCore::exprs(res$flowFrameList[[1]])), 6L)
  expect_true(all(res$labelsList[[1]] %in% sc6_custom$clusterLabelVec))
  expect_true(all(res$labelsList[[2]] %in% sc6_custom$clusterLabelVec))
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

test_that("simCytExperiment accepts a scenario object as a direct contract", {
  sc <- simCytScenarioBivariate()
  f_id <- simCytTransformIdentity()

  set.seed(123)
  res <- simCytExperiment(
    scenario = sc,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = 50,
    transformationFunc = f_id
  )

  set.seed(123)
  expected <- simCytExperiment(
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
  expect_equal(res$labelsList[["sample001_unstim"]],
    expected$labelsList[["sample001_unstim"]])
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

test_that("simCytBuildProbVec constructs valid probability vectors", {
  p <- simCytBuildProbVec(nMarker = 3L)
  expect_length(p, 8L)
  expect_equal(sum(p), 1)
  expect_true(all(p > 0))

  neg <- "F1-F2-F3-"
  p2 <- simCytBuildProbVec(nMarker = 3L, weights = c("F1-F2-F3-" = 100), defaultWeight = 1)
  expect_equal(sum(p2), 1)
  expect_true(p2[neg] > 0.9)

  expect_error(
    simCytBuildProbVec(nMarker = 2L, weights = c("INVALID" = 1)),
    "Unknown cluster labels"
  )
})

test_that("simCytBuildResponseVec constructs valid shift vectors", {
  sc <- simCytBuildScenario(nMarker = 3L)
  neg <- sc$clusterLabelVec[1]
  pos <- sc$clusterLabelVec[8]

  delta <- simCytBuildResponseVec(
    nMarker = 3L,
    shifts = setNames(c(-0.1, 0.1), c(neg, pos))
  )
  expect_length(delta, 8L)
  expect_equal(sum(delta), 0, tolerance = 1e-9)
  expect_equal(delta[neg], c("F1-F2-F3-" = -0.1))
  expect_equal(delta[pos], c("F1+F2+F3+" = 0.1))
  expect_true(all(delta[setdiff(names(delta), c(neg, pos))] == 0))

  expect_error(
    simCytBuildResponseVec(nMarker = 2L, shifts = c(INVALID = 0.1)),
    "Unknown cluster labels"
  )
  expect_error(
    simCytBuildResponseVec(nMarker = 2L, shifts = c("F1-F2-" = 0.1)),
    "must sum to zero"
  )
})

test_that("simCytBuildProbVec and simCytBuildResponseVec integrate with simCytBuildScenario for 6-marker case", {
  nM <- 6L
  sc <- simCytBuildScenario(nMarker = nM)
  neg <- sc$clusterLabelVec[1]
  pos <- sc$clusterLabelVec[64]

  probUns <- simCytBuildProbVec(nMarker = nM, weights = setNames(50, neg), defaultWeight = 1)
  response <- simCytBuildResponseVec(nMarker = nM, shifts = setNames(c(-0.1, 0.1), c(neg, pos)))

  expect_equal(sum(probUns), 1)
  expect_equal(sum(response), 0, tolerance = 1e-9)
  expect_equal(length(probUns), 64L)
  expect_equal(length(response), 64L)

  # Resulting stimulated probabilities must be non-negative and sum to 1
  stimProbs <- probUns + response
  expect_true(all(stimProbs >= 0))
  expect_equal(sum(stimProbs), 1, tolerance = 1e-9)
})

