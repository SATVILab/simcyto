covTestScenario <- function() {
  simCytBuildScenario(
    nMarker = 2,
    lowMean = 0,
    highMean = 5,
    probUns = c(0.4, 0.2, 0.2, 0.2),
    probResponse = list(c(-0.1, 0, 0, 0.1), c(-0.1, 0.1, 0, 0))
  )
}

# Run simCytExperiment() and record the covariance matrix handed to each
# cluster, in simulation order (sample, then condition, then cluster).
recordCovMats <- function(...) {
  rec <- new.env()
  rec$sigma <- list()
  testthat::local_mocked_bindings(
    simCytClusterData = function(mixtureType, clusterNumber, nCell, muVec,
                                 sigmaMat) {
      rec$sigma[[length(rec$sigma) + 1L]] <- sigmaMat
      matrix(rep(muVec, each = nCell), nrow = nCell)
    }
  )
  simCytExperiment(
    scenario = covTestScenario(),
    nSample = 2,
    nCondition = 3,
    nCellByCondition = 40,
    transformationFunc = simCytTransformIdentity(),
    probExact = TRUE,
    ...
  )
  idx <- expand.grid(cluster = 1:4, condition = 1:3, sample = 1:2)[, 3:1]
  list(idx = idx, sigma = rec$sigma)
}

test_that("simCytCovMat builds covariance matrices from SDs and correlations", {
  expect_equal(
    simCytCovMat(sd = sqrt(1.5), rho = 0.4, nMarker = 2),
    matrix(c(1.5, 0.6, 0.6, 1.5), 2)
  )
  expect_equal(
    simCytCovMat(sd = c(1, 2), rho = -0.5),
    matrix(c(1, -1, -1, 4), 2)
  )
  expect_equal(simCytCovMat(sd = 2), matrix(4, 1, 1))
  expect_equal(simCytCovMat(sd = sqrt(1.5), nMarker = 3), diag(1.5, 3))

  corMat <- matrix(c(1, 0.2, 0, 0.2, 1, -0.3, 0, -0.3, 1), 3)
  expect_equal(simCytCovMat(sd = 1, rho = corMat), corMat)

  expect_error(simCytCovMat(sd = -1))
  expect_error(simCytCovMat(sd = c(1, 2, 3), nMarker = 2), "length 1 or")
  expect_error(
    simCytCovMat(sd = 1, rho = -0.9, nMarker = 3),
    "positive definite"
  )
  expect_error(
    simCytCovMat(sd = 1, rho = matrix(c(1, 0.2, 0.3, 1), 2)),
    "correlation matrix"
  )
})

test_that("covVaryBy shares a drawn matrix across the strata it omits", {
  subsets <- list(
    character(0), "sample", "condition", "cluster",
    c("sample", "condition"), c("sample", "cluster"), c("condition", "cluster"),
    c("cluster", "sample", "condition")
  )
  for (covVaryBy in subsets) {
    set.seed(1)
    res <- recordCovMats(covEvMin = 0.5, covEvMax = 3, covVaryBy = covVaryBy)
    expect_length(res$sigma, nrow(res$idx))
    stratum <- if (length(covVaryBy) == 0L) {
      rep("", nrow(res$idx))
    } else {
      do.call(paste, res$idx[covVaryBy])
    }
    for (i in seq_len(nrow(res$idx))) {
      for (j in seq_len(i - 1L)) {
        same <- isTRUE(all.equal(res$sigma[[i]], res$sigma[[j]]))
        expect_identical(
          same, stratum[[i]] == stratum[[j]],
          info = paste0(
            "covVaryBy = ", paste(covVaryBy, collapse = "+"),
            ", rows ", i, " and ", j
          )
        )
      }
    }
  }
})

test_that("the default covVaryBy keeps the original per-tube draws", {
  set.seed(7)
  resDefault <- simCytExperiment(
    scenario = covTestScenario(), nSample = 2, nCondition = 3,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(), probExact = TRUE
  )
  set.seed(7)
  resExplicit <- simCytExperiment(
    scenario = covTestScenario(), nSample = 2, nCondition = 3,
    nCellByCondition = 50,
    transformationFunc = simCytTransformIdentity(), probExact = TRUE,
    covVaryBy = c("cluster", "condition", "sample")
  )
  expect_identical(
    lapply(resDefault$flowFrameList, flowCore::exprs),
    lapply(resExplicit$flowFrameList, flowCore::exprs)
  )

  set.seed(1)
  res <- recordCovMats(covEvMin = 0.5, covEvMax = 3)
  expect_length(unique(res$sigma), nrow(res$idx))
})

test_that("covMatList fixes each cluster's matrix in every sample/condition", {
  covMatList <- list(
    simCytCovMat(sqrt(1.5), rho = 0.5, nMarker = 2),
    simCytCovMat(c(1, 2), rho = 0),
    simCytCovMat(c(2, 1), rho = -0.3),
    simCytCovMat(sqrt(3), rho = 0.1, nMarker = 2)
  )
  res <- recordCovMats(covMatList = covMatList)
  for (i in seq_len(nrow(res$idx))) {
    expect_equal(res$sigma[[i]], covMatList[[res$idx$cluster[[i]]]])
  }

  single <- simCytCovMat(sqrt(1.5), rho = 0.4, nMarker = 2)
  res <- recordCovMats(covMatList = single)
  for (s in res$sigma) {
    expect_equal(s, single)
  }
})

test_that("covMatList sets the within-cluster covariance of the data", {
  sigmaNeg <- simCytCovMat(sqrt(1.5), rho = 0.6, nMarker = 2)
  sigmaPos <- simCytCovMat(c(1, 2), rho = -0.4)
  set.seed(11)
  res <- simCytExperiment(
    scenario = covTestScenario(),
    nSample = 1,
    nCondition = 3,
    nCellByCondition = 2e4,
    transformationFunc = simCytTransformGaussian(),
    probExact = TRUE,
    covMatList = list(sigmaNeg, sigmaPos, sigmaPos, sigmaPos)
  )
  for (nm in names(res$flowFrameList)) {
    x <- flowCore::exprs(res$flowFrameList[[nm]])
    lab <- res$labelsList[[nm]]
    neg <- x[lab == covTestScenario()$clusterLabelVec[[1]], , drop = FALSE]
    pos <- x[lab == covTestScenario()$clusterLabelVec[[2]], , drop = FALSE]
    expect_equal(unname(stats::cov(neg)), sigmaNeg, tolerance = 0.1)
    expect_equal(unname(stats::cov(pos)), sigmaPos, tolerance = 0.1)
  }
})

test_that("covariance arguments are validated", {
  run <- function(...) {
    simCytExperiment(
      scenario = covTestScenario(), nSample = 1, nCondition = 2,
      nCellByCondition = 20,
      transformationFunc = simCytTransformIdentity(), probExact = TRUE, ...
    )
  }
  expect_error(run(covVaryBy = "tube"), "covVaryBy")
  expect_error(run(covVaryBy = c("sample", "sample")), "covVaryBy")
  expect_error(run(covVaryBy = NA_character_), "covVaryBy")
  expect_error(
    run(covMatList = diag(2), covVaryBy = "cluster"),
    "cannot be combined"
  )
  expect_error(run(covMatList = diag(3)), "2 x 2")
  expect_error(run(covMatList = list(diag(2), diag(2))), "nCluster")
  expect_error(run(covMatList = matrix(c(1, 0.5, 0, 1), 2)), "symmetric")
  expect_error(run(covMatList = matrix(c(1, 2, 2, 1), 2)), "positive definite")
  expect_error(run(covEvMin = 2, covEvMax = 1, covVaryBy = "cluster"))
})

test_that("covMatList accepts scalar variances for univariate scenarios", {
  scenario <- simCytScenarioUnivariate()
  set.seed(3)
  res <- simCytExperiment(
    scenario = scenario, nSample = 1, nCondition = 2, nCellByCondition = 2e4,
    transformationFunc = simCytTransformGaussian(), probExact = TRUE,
    covMatList = list(1, 4)
  )
  x <- flowCore::exprs(res$flowFrameList[[1]])[, 1]
  lab <- res$labelsList[[1]]
  labs <- scenario$clusterLabelVec
  expect_equal(stats::var(x[lab == labs[[1]]]), 1, tolerance = 0.05)
  expect_equal(stats::var(x[lab == labs[[2]]]), 4, tolerance = 0.05)
})

covToSdCor <- function(m) {
  sds <- sqrt(diag(m))
  list(sd = sds, cor = m[1, 2] / (sds[[1]] * sds[[2]]))
}

test_that("fixed covSd and covCor give the same matrix everywhere", {
  res <- recordCovMats(covSd = sqrt(1.5), covCor = 0.3)
  target <- simCytCovMat(sqrt(1.5), rho = 0.3, nMarker = 2)
  for (s in res$sigma) {
    expect_equal(s, target)
  }
  res <- recordCovMats(covSd = 2)
  for (s in res$sigma) {
    expect_equal(s, diag(4, 2))
  }
})

test_that("covSd and covCor ranges give one shared matrix by default", {
  set.seed(5)
  res <- recordCovMats(covSd = c(0.5, 2), covCor = c(0, 0.5))
  parts <- covToSdCor(res$sigma[[1]])
  expect_true(all(parts$sd >= 0.5 & parts$sd <= 2))
  expect_true(parts$cor >= 0 && parts$cor <= 0.5)
  for (s in res$sigma) {
    expect_equal(s, res$sigma[[1]])
  }
})

test_that("covSd and covCor ranges can be drawn per cluster", {
  set.seed(5)
  res <- recordCovMats(
    covSd = c(0.5, 2), covCor = c(0, 0.5), covVaryBy = "cluster"
  )
  parts <- lapply(res$sigma, covToSdCor)
  sds <- unlist(lapply(parts, `[[`, "sd"))
  cors <- vapply(parts, `[[`, numeric(1), "cor")
  expect_true(all(sds >= 0.5 & sds <= 2))
  expect_true(all(cors >= 0 & cors <= 0.5))
  for (i in seq_len(nrow(res$idx))) {
    for (j in seq_len(i - 1L)) {
      expect_identical(
        isTRUE(all.equal(res$sigma[[i]], res$sigma[[j]])),
        res$idx$cluster[[i]] == res$idx$cluster[[j]]
      )
    }
  }
})

test_that("covSd and covCor ranges follow an explicit covVaryBy", {
  set.seed(6)
  res <- recordCovMats(
    covSd = c(0.5, 2), covCor = c(-0.3, 0.6),
    covVaryBy = c("sample", "cluster")
  )
  stratum <- paste(res$idx$sample, res$idx$cluster)
  for (i in seq_len(nrow(res$idx))) {
    for (j in seq_len(i - 1L)) {
      expect_identical(
        isTRUE(all.equal(res$sigma[[i]], res$sigma[[j]])),
        stratum[[i]] == stratum[[j]]
      )
    }
  }
})

test_that("covSd and covCor accept one entry per cluster", {
  set.seed(7)
  res <- recordCovMats(
    covSd = sqrt(1.5),
    covCor = list(c(0.2, 0.5), 0, 0, 0)
  )
  for (i in seq_len(nrow(res$idx))) {
    cr <- covToSdCor(res$sigma[[i]])$cor
    if (res$idx$cluster[[i]] == 1L) {
      expect_true(cr >= 0.2 && cr <= 0.5)
    } else {
      expect_equal(cr, 0)
    }
    expect_equal(diag(res$sigma[[i]]), c(1.5, 1.5))
  }
})

test_that("covCor ranges with more markers give positive-definite matrices", {
  set.seed(8)
  for (i in 1:200) {
    m <- simcyto:::.simCytDrawSdCor(4, c(1, 2), c(-0.4, 0.6))
    expect_false(inherits(try(chol(m), silent = TRUE), "try-error"))
    cors <- stats::cov2cor(m)[lower.tri(m)]
    expect_true(all(cors >= -0.4 & cors <= 0.6))
  }
  expect_error(
    simcyto:::.simCytDrawSdCor(3, 1, -0.9),
    "positive-definite"
  )
})

test_that("covSd and covCor are validated", {
  run <- function(...) {
    simCytExperiment(
      scenario = covTestScenario(), nSample = 1, nCondition = 3,
      nCellByCondition = 20, transformationFunc = simCytTransformIdentity(),
      probExact = TRUE, ...
    )
  }
  expect_error(run(covCor = 0.3), "requires `covSd`")
  expect_error(run(covSd = 1, covMatList = diag(2)), "cannot be combined")
  expect_error(run(covSd = -1), "covSd")
  expect_error(run(covSd = c(2, 1)), "covSd")
  expect_error(run(covSd = 1, covCor = 1.2), "covCor")
  expect_error(run(covSd = 1, covCor = list(0, 0)), "nCluster")
  expect_error(
    run(covSd = 1, covCor = list(0.3, 0, 0, 0), covVaryBy = "sample"),
    "requires \"cluster\""
  )
})
