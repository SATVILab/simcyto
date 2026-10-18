mpTestScenario <- function() {
  simCytBuildScenario(
    nMarker = 2,
    lowMean = 0,
    highMean = 5,
    probUns = c(0.4, 0.2, 0.2, 0.2),
    probResponse = list(c(-0.1, 0, 0, 0.1), c(-0.1, 0.1, 0, 0))
  )
}

# Run simCytExperiment() and record each cluster's mean shift from the
# scenario means, in simulation order (sample, then condition, then cluster).
recordMeanShifts <- function(...) {
  rec <- new.env()
  rec$mu <- list()
  testthat::local_mocked_bindings(
    simCytClusterData = function(mixtureType, clusterNumber, nCell, muVec,
                                 sigmaMat) {
      rec$mu[[length(rec$mu) + 1L]] <- muVec
      matrix(rep(muVec, each = nCell), nrow = nCell)
    }
  )
  scenario <- mpTestScenario()
  simCytExperiment(
    scenario = scenario,
    nSample = 2,
    nCondition = 3,
    nCellByCondition = 40,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE,
    covEvMin = 1,
    covEvMax = 1,
    ...
  )
  idx <- expand.grid(cluster = 1:4, condition = 1:3, sample = 1:2)[, 3:1]
  shift <- lapply(seq_len(nrow(idx)), function(i) {
    rec$mu[[i]] - scenario$meanExprMat[idx$cluster[[i]], ]
  })
  list(idx = idx, shift = shift)
}

test_that("meanPerturbation shares shifts exactly across omitted strata", {
  subsets <- list(
    character(0), "sample", "condition", "cluster",
    c("sample", "condition"), c("sample", "cluster"),
    c("condition", "cluster"), c("sample", "condition", "cluster")
  )
  for (varyBy in subsets) {
    set.seed(1)
    res <- recordMeanShifts(
      meanPerturbation = list(sd = 1, varyBy = varyBy)
    )
    stratum <- if (length(varyBy) == 0L) {
      rep("", nrow(res$idx))
    } else {
      do.call(paste, res$idx[varyBy])
    }
    for (i in seq_len(nrow(res$idx))) {
      for (j in seq_len(i - 1L)) {
        expect_identical(
          isTRUE(all.equal(res$shift[[i]], res$shift[[j]])),
          stratum[[i]] == stratum[[j]],
          info = paste0("varyBy = ", paste(varyBy, collapse = "+"))
        )
      }
    }
  }
})

test_that("meanPerturbation layers add up", {
  set.seed(2)
  res <- recordMeanShifts(
    meanPerturbation = list(
      list(sd = 1, varyBy = "sample"),
      list(sd = 0.5, varyBy = c("condition", "cluster"))
    )
  )
  # Shift = sample part + condition-by-cluster part, so differences between
  # samples are the same for every condition and cluster.
  diffs <- lapply(seq_len(12), function(i) res$shift[[i + 12]] - res$shift[[i]])
  for (d in diffs) {
    expect_equal(d, diffs[[1]])
  }
})

test_that("meanPerturbation leaves the default random-number stream alone", {
  run <- function(...) {
    set.seed(9)
    res <- simCytExperiment(
      scenario = mpTestScenario(), nSample = 2, nCondition = 3,
      nCellByCondition = 30, transformationFunc = simCytTransformIdentity(),
      probExact = TRUE, ...
    )
    lapply(res$flowFrameList, flowCore::exprs)
  }
  expect_identical(run(), run(meanPerturbation = NULL))
  expect_identical(
    run(),
    run(meanPerturbation = list(sd = 0, varyBy = "sample"))
  )
})

test_that("meanPerturbation is validated", {
  run <- function(...) {
    simCytExperiment(
      scenario = mpTestScenario(), nSample = 1, nCondition = 3,
      nCellByCondition = 20, transformationFunc = simCytTransformIdentity(),
      probExact = TRUE, ...
    )
  }
  expect_error(run(meanPerturbation = list(sd = -1, varyBy = "sample")), "sd")
  expect_error(run(meanPerturbation = list(sd = 1)), "sd")
  expect_error(
    run(meanPerturbation = list(sd = 1, varyBy = "tube")),
    "varyBy"
  )
  expect_error(run(meanPerturbation = list(1, 2)), "meanPerturbation")
})

test_that("stimulated probabilities must stay valid", {
  run <- function(probResponse, probExact = TRUE) {
    simCytExperiment(
      nSample = 1, nMarker = 1, nCondition = 2, nCluster = 2,
      nCellByCondition = 50, transformationFunc = simCytTransformIdentity(),
      meanExprMat = matrix(c(0, 5), ncol = 1), clusterLabelVec = c("-", "+"),
      probVecUns = c(0.9, 0.1), probExact = probExact,
      probResponseVecByStimCondition = list(probResponse)
    )
  }
  expect_error(run(c(-0.1, 0.2)), "sum to 1")
  expect_error(run(c(0.2, -0.2)), "between 0 and 1")
  expect_error(run(c(-0.1, 0.2), probExact = FALSE), "sum to 1")
  expect_silent(run(c(-0.1, 0.1)))
  expect_silent(run(c(-0.9, 0.9)))
  expect_silent(run(c(-0.1 + 1e-12, 0.1)))
})
