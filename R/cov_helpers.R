#' Generate a random positive definite matrix
#'
#' @keywords internal
.posDef <- function(n, covEvMin = 1, covEvMax = 2) {
  ev <- stats::runif(n, min = covEvMin, max = covEvMax)
  if (n == 1) {
    return(matrix(ev, 1, 1))
  }
  z <- matrix(ncol = n, rnorm(n^2))
  decomp <- qr(z)
  q <- qr.Q(decomp)
  r <- qr.R(decomp)
  d <- diag(r)
  ph <- d / abs(d)
  o <- q %*% diag(ph, nrow = n)
  z <- t(o) %*% diag(ev, nrow = n) %*% o
  z
}
