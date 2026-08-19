test_that("public high-level API examples run end-to-end", {
  scenario <- simCytBuildScenario(
    nMarker = 2,
    lowMean = 0,
    highMean = 5,
    probUns = c(0.7, 0.1, 0.1, 0.1),
    probResponse = list(c(-0.1, 0, 0, 0.1))
  )

  set.seed(2026)
  sim_res <- simCytExperiment(
    scenario = scenario,
    nSample = 1,
    nCondition = 2,
    nCellByCondition = c(30, 40),
    transformationFunc = simCytTransformIdentity(),
    mixtureType = "gaussianOnly",
    probExact = TRUE
  )

  expect_named(sim_res, c("flowFrameList", "labelsList"))
  expect_equal(names(sim_res$flowFrameList), names(sim_res$labelsList))
  expect_equal(dim(flowCore::exprs(sim_res$flowFrameList[[1]])), c(30, 2))
  expect_true(all(sim_res$labelsList[[1]] %in% scenario$clusterLabelVec))
})
