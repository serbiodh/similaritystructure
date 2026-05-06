#' Generate skew-normal random values with fixed mean and variance
#'
#' @param n Number of observations.
#' @param alpha Shape parameter of the skew-normal distribution.
#' @param mu_target Target mean.
#' @param var_target Target variance.
#'
#' @return Numeric vector of simulated values.
#'
#' @keywords internal

rskew_fixed <- function(n, alpha, mu_target = 0, var_target = 1){

  if (!requireNamespace("sn", quietly = TRUE)) {
    stop("Package 'sn' is required for rskew_fixed(). Please install it.")
  }

  xi <- 0
  omega <- 1

  # generate standard skew-normal
  x <- sn::rsn(n, xi = xi, omega = omega, alpha = alpha)

  # theoretical moments
  delta <- alpha / sqrt(1 + alpha^2)
  mu_x <- xi + omega * delta * sqrt(2/pi)
  sigma_x <- omega * sqrt(1 - 2*delta^2/pi)

  # recenter y rescale
  y <- (x - mu_x)/sigma_x

  # adjust target mean and variance
  y <- sqrt(var_target)*y + mu_target
  y <- as.vector(y)
  return(y)

}
