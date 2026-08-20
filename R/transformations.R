#' Get a simulation transformation function
#'
#' Factory for marker-wise post-simulation transformations used by
#' [simCytExperiment()].
#'
#' `mixtureType` controls how raw cluster values are simulated. The
#' transformation returned here is applied afterwards, marker-wise, to the
#' simulated values.
#'
#' @param type Character scalar. One of `"identity"`, `"gaussian"`, `"gamma"`,
#'   `"gammaFixed"`, or `"skew"`.
#' @param ... Additional arguments passed to the selected constructor.
#'
#' @return A function that takes a numeric vector `x` and returns transformed
#'   values. The returned function carries a `"sim_transformation"` attribute used
#'   by ratio-correction logic in gamma/skew pipelines.
#'
#' @examples
#' f_id <- simCytGetTransformation("identity")
#' f_id(c(1, 2, 3))
#'
#' f_gauss <- simCytGetTransformation("gaussian")
#' f_gauss(c(1, 2, 3))
#'
#' @export
simCytGetTransformation <- function(type = c("identity", "gaussian", "gamma", "gammaFixed", "skew"), ...) {
  type <- match.arg(type)
  switch(
    type,
    identity = simCytTransformIdentity(...),
    gaussian = simCytTransformGaussian(...),
    gamma = simCytTransformGamma(...),
    gammaFixed = simCytTransformGammaFixed(...),
    skew = simCytTransformSkew(...)
  )
}

#' Identity transformation
#'
#' @return Function that returns `x` unchanged.
#'
#' @examples
#' f <- simCytTransformIdentity()
#' f(c(-1, 0, 1))
#'
#' @export
simCytTransformIdentity <- function() {
  f <- function(x) x
  attr(f, "sim_transformation") <- "identity"
  f
}

#' Gaussian/identity transformation
#'
#' Reproduces the current StimGate `"gaussian"` post-simulation setting, which
#' applies no additional transformation.
#'
#' @return Function that returns `x` unchanged.
#'
#' @examples
#' f <- simCytTransformGaussian()
#' f(c(2, 4, 6))
#'
#' @export
simCytTransformGaussian <- function() {
  f <- function(x) x
  attr(f, "sim_transformation") <- "gaussian"
  f
}

#' StimGate gamma transformation
#'
#' Reproduces the current StimGate `calc_gamma()` helper:
#' `gamma(1 + abs(x / 4))`.
#'
#' @return Function that transforms `x` with the StimGate gamma formula.
#'
#' @examples
#' f <- simCytTransformGamma()
#' f(c(-4, 0, 4))
#'
#' @export
simCytTransformGamma <- function() {
  f <- function(x) {
    gamma(1 + abs(x / 4))
  }
  attr(f, "sim_transformation") <- "gamma"
  f
}

#' StimGate gamma transformation with fixed mean and spread
#'
#' Reproduces the current StimGate `calc_gamma_fixed_mean_and_spread()` helper:
#' apply the StimGate gamma transformation, then linearly rescale the result so
#' the output has the same mean and standard deviation as the input vector.
#'
#' @return Function that applies the fixed-mean/spread StimGate gamma
#'   transformation.
#'
#' @examples
#' f <- simCytTransformGammaFixed()
#' f(c(1, 2, 3))
#'
#' @export
simCytTransformGammaFixed <- function() {
  gammaTransform <- simCytTransformGamma()
  f <- function(x) {
    clusterMean <- mean(x)
    clusterSd <- stats::sd(x)
    out <- gammaTransform(x)
    (out - mean(out)) / stats::sd(out) * clusterSd + clusterMean
  }
  attr(f, "sim_transformation") <- "gamma_fixed_mean_and_spread"
  f
}

#' StimGate skew transformation
#'
#' Reproduces the current StimGate `calc_skew()` helper, including the
#' mean-dependent `epsilon` weighting and stochastic gamma-divisor term.
#'
#' @param epsilon Numeric scalar. Base epsilon parameter. Default is `0.5`.
#' @param delta Numeric scalar. Delta parameter. Default is `1`.
#'
#' @return Function that applies the StimGate skew transformation to `x`.
#'
#' @examples
#' f <- simCytTransformSkew()
#' set.seed(1)
#' f(c(1, 1, 1))
#'
#' @export
simCytTransformSkew <- function(epsilon = 0.5, delta = 1) {
  f <- function(x) {
    clusterMean <- mean(x)
    weight <- 1 / (1 + exp(1.5 * (clusterMean - 3.5)))
    epsilonWeighted <- epsilon * weight
    out <- sinh(epsilonWeighted + delta * asinh(x))
    gammaDivisor <- stats::rgamma(length(x), shape = 5, rate = 5)
    out / sqrt(gammaDivisor)
  }
  attr(f, "sim_transformation") <- "skew"
  f
}
