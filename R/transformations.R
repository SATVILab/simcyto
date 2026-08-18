#' Get a simulation transformation function
#'
#' Factory function to create tagged transformation functions used in simcyto
#' simulations.
#'
#' @param type Character. One of "identity", "gaussian", "gamma", "gammaFixed", or "skew".
#' @param ... Additional arguments passed to the specific transformation constructor.
#'
#' @return A function that takes a numeric vector `x` and returns transformed values,
#'   tagged with attribute `"sim_transformation"`.
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

#' Identity transformation for cytometry simulation
#'
#' @return A function returning `x` unchanged, tagged with `sim_transformation = "identity"`.
#' @export
simCytTransformIdentity <- function() {
  f <- function(x) x
  attr(f, "sim_transformation") <- "identity"
  f
}

#' Gaussian noise transformation
#'
#' @param sd Numeric. Standard deviation of Gaussian noise. Default is 1.
#'
#' @return A function that adds zero-mean Gaussian noise to `x`, tagged with `sim_transformation = "gaussian"`.
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
#' @param shape Numeric. Shape parameter for Gamma distribution. Default is 1.
#' @param scale Numeric. Scale parameter for Gamma distribution. Default is 1.
#'
#' @return A function that adds Gamma-distributed noise to `x`, tagged with `sim_transformation = "gamma"`.
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
#' @param mean Numeric. Mean of Gamma noise. Default is 1.
#' @param sd Numeric. Standard deviation of Gamma noise. Default is 1.
#'
#' @return A function that adds Gamma-distributed noise to `x`, tagged with `sim_transformation = "gamma"`.
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
#' @param shape Numeric. Shape (skewness) parameter alpha. Default is 0.
#' @param location Numeric. Location parameter. Default is 0.
#' @param scale Numeric. Scale parameter (>0). Default is 1.
#'
#' @return A function that adds skew-normal noise to `x`, tagged with `sim_transformation = "skew"`.
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
