#' Get a simulation transformation function
#'
#' Factory for marker-wise transformations used by [simCytExperiment()].
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
#' f_gauss <- simCytGetTransformation("gaussian", sd = 0.1)
#' set.seed(1)
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

#' Gaussian noise transformation
#'
#' @param sd Numeric scalar. Standard deviation of zero-mean Gaussian noise.
#'
#' @return Function that adds Gaussian noise to `x`.
#'
#' @examples
#' f <- simCytTransformGaussian(sd = 0)
#' f(c(2, 4, 6))
#'
#' @export
simCytTransformGaussian <- function(sd = 1) {
  f <- function(x) {
    x + stats::rnorm(length(x), mean = 0, sd = sd)
  }
  attr(f, "sim_transformation") <- "gaussian"
  f
}

#' Gamma noise transformation
#'
#' @param shape Numeric scalar. Gamma shape parameter.
#' @param scale Numeric scalar. Gamma scale parameter.
#'
#' @return Function that adds gamma-distributed noise to `x`.
#'
#' @examples
#' f <- simCytTransformGamma(shape = 2, scale = 0.5)
#' set.seed(1)
#' f(c(0, 0, 0))
#'
#' @export
simCytTransformGamma <- function(shape = 1, scale = 1) {
  f <- function(x) {
    x + stats::rgamma(length(x), shape = shape, scale = scale)
  }
  attr(f, "sim_transformation") <- "gamma"
  f
}

#' Gamma noise transformation with fixed mean and standard deviation
#'
#' Reparameterises Gamma distribution using mean and standard deviation.
#'
#' @param mean Numeric scalar. Mean of gamma noise (`> 0`).
#' @param sd Numeric scalar. Standard deviation of gamma noise (`> 0`).
#'
#' @return Function that adds gamma-distributed noise to `x`.
#'
#' @examples
#' f <- simCytTransformGammaFixed(mean = 2, sd = 1)
#' set.seed(1)
#' f(c(0, 0, 0))
#'
#' @export
simCytTransformGammaFixed <- function(mean = 1, sd = 1) {
  stopifnot(mean > 0, sd > 0)
  shape <- (mean / sd)^2
  scale <- (sd^2) / mean
  simCytTransformGamma(shape = shape, scale = scale)
}

#' Skew-normal noise transformation
#'
#' Adds skew-normal distributed noise to `x`.
#'
#' @param shape Numeric scalar. Shape (skewness) parameter.
#' @param location Numeric scalar. Location parameter.
#' @param scale Numeric scalar. Scale parameter (`> 0`).
#'
#' @return Function that adds skew-normal noise to `x`.
#'
#' @examples
#' f <- simCytTransformSkew(shape = 2, location = 0, scale = 1)
#' set.seed(1)
#' f(c(1, 1, 1))
#'
#' @export
simCytTransformSkew <- function(shape = 0, location = 0, scale = 1) {
  stopifnot(scale > 0)
  f <- function(x) {
    n <- length(x)
    u0 <- stats::rnorm(n)
    u1 <- stats::rnorm(n)
    delta <- shape / sqrt(1 + shape^2)
    z <- delta * abs(u0) + sqrt(1 - delta^2) * u1
    skewNoise <- location + scale * z
    x + skewNoise
  }
  attr(f, "sim_transformation") <- "skew"
  f
}
