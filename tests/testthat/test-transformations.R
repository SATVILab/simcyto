stimgate_sim_misc_ref <- paste(
  "SATVILab/stimgate",
  "scripts/r/sim-misc.R",
  "@930dd907823ccd23d5acfdf3db7527db5a19241b",
  sep = ""
)

stimgate_tag_transform <- function(fun, name) {
  attr(fun, "sim_transformation") <- name
  fun
}

stimgate_calc_gaussian <- function(x) {
  x
}

stimgate_calc_gamma <- function(x) {
  gamma(1 + abs(x / 4))
}

stimgate_calc_gamma_fixed_mean_and_spread <- function(x) {
  cluster_mean <- mean(x)
  cluster_sd <- stats::sd(x)
  out <- stimgate_calc_gamma(x)
  (out - mean(out)) / stats::sd(out) * cluster_sd + cluster_mean
}

stimgate_calc_skew <- function(x, epsilon = 0.5, delta = 1) {
  cluster_mean <- mean(x)
  weight <- 1 / (1 + exp(1.5 * (cluster_mean - 3.5)))
  epsilon <- epsilon * weight
  out <- sinh(epsilon + delta * asinh(x))
  gamma_divisor <- stats::rgamma(length(x), shape = 5, rate = 5)
  out / sqrt(gamma_divisor)
}

compare_exprs <- function(actual, expected, tolerance = 0) {
  actual_mat <- as.matrix(flowCore::exprs(actual))
  expected_mat <- as.matrix(flowCore::exprs(expected))
  attr(actual_mat, "ranges") <- NULL
  attr(expected_mat, "ranges") <- NULL
  expect_equal(actual_mat, expected_mat, tolerance = tolerance)
}

test_that("simCytGetTransformation creates the expected StimGate-style transformations", {
  expect_match(stimgate_sim_misc_ref, "930dd907823ccd23d5acfdf3db7527db5a19241b")

  f_id <- simCytGetTransformation("identity")
  expect_equal(attr(f_id, "sim_transformation"), "identity")
  expect_equal(f_id(c(1, 2, 3)), c(1, 2, 3))

  f_gauss <- simCytGetTransformation("gaussian")
  expect_equal(attr(f_gauss, "sim_transformation"), "gaussian")
  expect_equal(f_gauss(c(1, 2, 3)), c(1, 2, 3))

  f_gamma <- simCytGetTransformation("gamma")
  expect_equal(attr(f_gamma, "sim_transformation"), "gamma")

  f_gamma_fixed <- simCytGetTransformation("gammaFixed")
  expect_equal(
    attr(f_gamma_fixed, "sim_transformation"),
    "gamma_fixed_mean_and_spread"
  )

  f_skew <- simCytGetTransformation("skew")
  expect_equal(attr(f_skew, "sim_transformation"), "skew")
})

test_that("StimGate gaussian and gamma transformations match the pinned formulas exactly", {
  x <- c(-8, -4, -1, 0, 2, 4, 8)

  expect_equal(simCytTransformGaussian()(x), stimgate_calc_gaussian(x))
  expect_equal(simCytTransformIdentity()(x), stimgate_calc_gaussian(x))
  expect_equal(simCytTransformGamma()(x), stimgate_calc_gamma(x))
})

test_that("StimGate fixed-mean/spread gamma transformation preserves mean and spread", {
  x <- c(-2, -0.5, 1, 3, 5, 9)

  actual <- simCytTransformGammaFixed()(x)
  expected <- stimgate_calc_gamma_fixed_mean_and_spread(x)

  expect_equal(actual, expected)
  expect_equal(mean(actual), mean(x))
  expect_equal(stats::sd(actual), stats::sd(x))
})

test_that("StimGate skew transformation matches the pinned formula under a fixed seed", {
  x <- c(0.1, 0.5, 1, 2, 4, 8)

  set.seed(42)
  actual <- simCytTransformSkew()(x)

  set.seed(42)
  expected <- stimgate_calc_skew(x)

  expect_equal(actual, expected)
})

test_that("Ratio correction logic keys off the StimGate transformation tags", {
  f_gamma <- simCytTransformGamma()
  expect_true(simcyto:::.simCytUsesUpperRatioCorrection(f_gamma))

  f_skew <- simCytTransformSkew()
  expect_true(simcyto:::.simCytUsesUpperRatioCorrection(f_skew))

  f_gauss <- simCytTransformGaussian()
  expect_false(simcyto:::.simCytUsesUpperRatioCorrection(f_gauss))

  f_gamma_fixed <- simCytTransformGammaFixed()
  expect_false(simcyto:::.simCytUsesUpperRatioCorrection(f_gamma_fixed))
})

test_that("simCytExperiment matches pinned StimGate gamma and skew transformations end to end", {
  mean_mat <- matrix(c(1, 8), nrow = 2, ncol = 1)
  labels <- c("F1-", "F1+")
  prob_uns <- c(0.5, 0.5)
  prob_resp <- list(c(-0.15, 0.15))

  run_experiment <- function(transformation_func) {
    simCytExperiment(
      nSample = 1,
      nMarker = 1,
      nCondition = 2,
      nCluster = 2,
      nCellByCondition = 80,
      transformationFunc = transformation_func,
      mixtureType = "tPlusGauss",
      meanExprMat = mean_mat,
      clusterLabelVec = labels,
      probVecUns = prob_uns,
      probExact = FALSE,
      probResponseVecByStimCondition = prob_resp,
      samplePerturbationSd = 0.05,
      conditionPerturbationSd = 0.05,
      clusterPerturbationSd = 0.05
    )
  }

  gamma_ref <- stimgate_tag_transform(stimgate_calc_gamma, "gamma")
  skew_ref <- stimgate_tag_transform(stimgate_calc_skew, "skew")

  set.seed(20260819)
  gamma_actual <- run_experiment(simCytTransformGamma())
  set.seed(20260819)
  gamma_expected <- run_experiment(gamma_ref)

  expect_equal(gamma_actual$labelsList, gamma_expected$labelsList)
  compare_exprs(
    gamma_actual$flowFrameList[["sample001_stim1"]],
    gamma_expected$flowFrameList[["sample001_stim1"]]
  )

  set.seed(20260820)
  skew_actual <- run_experiment(simCytTransformSkew())
  set.seed(20260820)
  skew_expected <- run_experiment(skew_ref)

  expect_equal(skew_actual$labelsList, skew_expected$labelsList)
  compare_exprs(
    skew_actual$flowFrameList[["sample001_stim1"]],
    skew_expected$flowFrameList[["sample001_stim1"]]
  )
})

test_that("simCytTransformGammaFixed handles single element and zero-variance inputs safely", {
  f_gamma_fixed <- simCytTransformGammaFixed()
  expect_equal(f_gamma_fixed(5), 5)
  expect_equal(f_gamma_fixed(c(3, 3, 3)), c(3, 3, 3))
  expect_equal(f_gamma_fixed(numeric(0)), numeric(0))
})
