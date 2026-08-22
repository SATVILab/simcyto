test_that(".posDef generates valid positive-definite matrix", {
  set.seed(123)
  m1 <- simcyto:::.posDef(1)
  expect_equal(dim(m1), c(1, 1))
  expect_true(m1[1, 1] > 0)

  m3 <- simcyto:::.posDef(3, covEvMin = 1, covEvMax = 2)
  expect_equal(dim(m3), c(3, 3))
  eigs <- eigen(m3, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs >= 0.99))
})

test_that("simCytClusterData supports gaussianOnly, tOnly, and tPlusGauss", {
  set.seed(42)
  mu <- c(1, 2)
  sigma <- diag(2)

  d_gauss <- simCytClusterData("gaussianOnly", clusterNumber = 1, nCell = 10, muVec = mu, sigmaMat = sigma)
  expect_equal(dim(d_gauss), c(10, 2))

  d_t <- simCytClusterData("tOnly", clusterNumber = 1, nCell = 10, muVec = mu, sigmaMat = sigma)
  expect_equal(dim(d_t), c(10, 2))

  d_tg1 <- simCytClusterData("tPlusGauss", clusterNumber = 1, nCell = 10, muVec = mu, sigmaMat = sigma)
  expect_equal(dim(d_tg1), c(10, 2))

  d_tg2 <- simCytClusterData("tPlusGauss", clusterNumber = 2, nCell = 10, muVec = mu, sigmaMat = sigma)
  expect_equal(dim(d_tg2), c(10, 2))
})

test_that("validateExperimentInputs and validateSampleInputs check parameters accurately", {
  identityTrans <- identity

  expect_silent(validateExperimentInputs(
    nSample = 2,
    nMarker = 2,
    nCondition = 2,
    nCluster = 4,
    nCellByCondition = 100,
    transformationFunc = identityTrans
  ))

  expect_error(validateExperimentInputs(
    nSample = 0,
    nMarker = 2,
    nCondition = 2,
    nCluster = 4,
    nCellByCondition = 100,
    transformationFunc = identityTrans
  ))

  expect_error(validateExperimentInputs(
    nSample = 2,
    nMarker = 2,
    nCondition = 2,
    nCluster = 3, # Invalid cluster count for 2 markers
    nCellByCondition = 100,
    transformationFunc = identityTrans
  ))

  meanMat <- matrix(c(0, 0, 0, 5, 5, 0, 5, 5), nrow = 4, ncol = 2, byrow = TRUE)
  labels <- c("c1", "c2", "c3", "c4")
  probUns <- c(0.25, 0.25, 0.25, 0.25)

  expect_silent(validateSampleInputs(
    nCondition = 2L,
    nMarker = 2L,
    nCluster = 4L,
    nCellByCondition = 100L,
    transformationFunc = identityTrans,
    probResponseVecByStimCondition = list(c(0, 0, 0, 0)),
    probVecUns = probUns,
    probExact = FALSE,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    meanExprMat = meanMat,
    clusterLabelVec = labels,
    covEvMin = 1,
    covEvMax = 2
  ))

  # Also accept numeric scalar values (e.g., from interactive calls)
  expect_silent(validateSampleInputs(
    nCondition = 2,
    nMarker = 2,
    nCluster = 4,
    nCellByCondition = 100,
    transformationFunc = identityTrans,
    probResponseVecByStimCondition = list(c(0, 0, 0, 0)),
    probVecUns = probUns,
    probExact = FALSE,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    meanExprMat = meanMat,
    clusterLabelVec = labels,
    covEvMin = 1,
    covEvMax = 2
  ))

  expect_error(validateSampleInputs(
    nCondition = 2L,
    nMarker = 2L,
    nCluster = 4L,
    nCellByCondition = 100L,
    transformationFunc = identityTrans,
    probResponseVecByStimCondition = list(c(0, 0, 0)), # wrong length
    probVecUns = probUns,
    probExact = FALSE,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    meanExprMat = meanMat,
    clusterLabelVec = labels,
    covEvMin = 1,
    covEvMax = 2
  ))

  expect_error(validateSampleInputs(
    nCondition = 2L,
    nMarker = 2L,
    nCluster = 4L,
    nCellByCondition = 100L,
    transformationFunc = identityTrans,
    probResponseVecByStimCondition = list(c(0, 0, 0, 0)),
    probVecUns = probUns,
    probExact = FALSE,
    conditionPerturbationSd = 0,
    clusterPerturbationSd = 0,
    meanExprMat = meanMat,
    clusterLabelVec = labels,
    covEvMin = 2, # covEvMin > covEvMax
    covEvMax = 1
  ))
})
