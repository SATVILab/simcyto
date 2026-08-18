test_that("simCytGetTransformation creates expected transformation functions", {
  f_id <- simCytGetTransformation("identity")
  expect_equal(attr(f_id, "sim_transformation"), "identity")
  expect_equal(f_id(c(1, 2, 3)), c(1, 2, 3))

  f_gauss <- simCytGetTransformation("gaussian", sd = 0)
  expect_equal(attr(f_gauss, "sim_transformation"), "gaussian")
  expect_equal(f_gauss(c(1, 2, 3)), c(1, 2, 3))

  f_gamma <- simCytGetTransformation("gamma", shape = 1, scale = 1)
  expect_equal(attr(f_gamma, "sim_transformation"), "gamma")

  f_gamma_fixed <- simCytGetTransformation("gammaFixed", mean = 2, sd = 1)
  expect_equal(attr(f_gamma_fixed, "sim_transformation"), "gamma")

  f_skew <- simCytGetTransformation("skew", shape = 2, location = 0, scale = 1)
  expect_equal(attr(f_skew, "sim_transformation"), "skew")
})

test_that("Transformations produce reproducible output when seeded", {
  set.seed(42)
  f1 <- simCytTransformGaussian(sd = 0.5)
  x1 <- f1(c(10, 20, 30))

  set.seed(42)
  f2 <- simCytTransformGaussian(sd = 0.5)
  x2 <- f2(c(10, 20, 30))

  expect_equal(x1, x2)
})

test_that("Ratio correction logic detects tagged transformations", {
  f_gamma <- simCytTransformGamma(shape = 2, scale = 1)
  expect_true(simcyto:::.simCytUsesUpperRatioCorrection(f_gamma))

  f_skew <- simCytTransformSkew(shape = 1, location = 0, scale = 1)
  expect_true(simcyto:::.simCytUsesUpperRatioCorrection(f_skew))

  f_id <- simCytTransformIdentity()
  expect_false(simcyto:::.simCytUsesUpperRatioCorrection(f_id))
})
