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
    expect_equal(unname(colnames(flowCore::exprs(ff))), c("F1", "F2"))

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
  gammaTrans <- function(x) x^0.5
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

  expect_equal(flowCore::exprs(res1$flowFrameList[["sample001_unstim"]]),
               flowCore::exprs(res2$flowFrameList[["sample001_unstim"]]))
  expect_equal(res1$labelsList[["sample001_unstim"]],
               res2$labelsList[["sample001_unstim"]])
})

test_that("fixed-seed regression test validates mixture types and ratio corrections reproducibility", {
  sc <- simCytScenarioBivariate()
  gamma_trans <- simCytTransformGamma(shape = 2, scale = 1)

  set.seed(999)
  res_gamma1 <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 80,
    transformationFunc = gamma_trans,
    mixtureType = "tOnly",
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition,
    samplePerturbationSd = 0.05,
    conditionPerturbationSd = 0.05,
    clusterPerturbationSd = 0.05
  )

  set.seed(999)
  res_gamma2 <- simCytExperiment(
    nSample = 1,
    nMarker = sc$nMarker,
    nCondition = 2,
    nCluster = sc$nCluster,
    nCellByCondition = 80,
    transformationFunc = gamma_trans,
    mixtureType = "tOnly",
    meanExprMat = sc$meanExprMat,
    clusterLabelVec = sc$clusterLabelVec,
    probVecUns = sc$probVecUns,
    probResponseVecByStimCondition = sc$probResponseVecByStimCondition,
    samplePerturbationSd = 0.05,
    conditionPerturbationSd = 0.05,
    clusterPerturbationSd = 0.05
  )

  expect_equal(flowCore::exprs(res_gamma1$flowFrameList[["sample001_stim1"]]),
               flowCore::exprs(res_gamma2$flowFrameList[["sample001_stim1"]]))
})
