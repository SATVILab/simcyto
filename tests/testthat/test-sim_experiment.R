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
