# Estimate a binomial GAM using the original Monte Carlo counts.
# This helper is shared by the main analysis and the cached-result updater.
smooth_binomial <- function(N_init, N_fin, propsig, num_repet) {

  if (!is.numeric(num_repet) || length(num_repet) != 1L ||
      !is.finite(num_repet) || num_repet < 1 || num_repet != floor(num_repet)) {
    stop("num_repet must be a positive integer.")
  }
  if (!all(c("N", "num_nosig") %in% names(propsig))) {
    stop("The empirical curve must contain N and num_nosig.")
  }

  propsig <- as.data.frame(propsig)
  propsig <- propsig[order(propsig$N), , drop = FALSE]
  g <- propsig$num_nosig
  if (length(g) < 3L || any(!is.finite(g)) ||
      any(g < 0 | g > num_repet | g != floor(g)) || g[1L] <= 0 ||
      any(!is.finite(propsig$N)) || anyDuplicated(propsig$N)) {
    stop("Invalid empirical counts or subsample sizes for binomial smoothing.")
  }

  vect_N_smooth <- seq.int(N_init, N_fin)
  propsig$num_sig <- num_repet - g

  # Model the non-rejection probability, with no monotonicity constraint.
  # The basis dimension is the same as in the original Gaussian GAM.
  smooth_formula <- stats::as.formula(
    sprintf("cbind(num_nosig, num_sig) ~ s(N, bs = 'tp', k = %d)",
            min(10L, nrow(propsig)))
  )

  smoothdata <- mgcv::gam(
    smooth_formula, data = propsig,
    family = stats::binomial(link = "logit"), method = "REML"
  )
  if (!isTRUE(smoothdata$converged)) {
    stop("Binomial GAM did not converge; inspect the empirical curve.")
  }

  # Probabilities are in [0,1]. Rescaling is only for comparison with the
  # empirical and fitted curves in the plot.
  beta_smooth <- stats::predict(
    smoothdata, newdata = data.frame(N = vect_N_smooth), type = "response"
  )
  fun_smooth <- beta_smooth / (g[1L] / num_repet)

  # Normalize over the finite support; the constant scale cancels.
  total <- sum(beta_smooth)
  if (!is.finite(total) || total <= 0 || any(beta_smooth < 0)) {
    stop("The binomial GAM curve cannot be normalized.")
  }
  Sn <- beta_smooth / total

  list(vect_N_smooth = vect_N_smooth,
       fun_smooth = fun_smooth, Sn = Sn,
       Ns_smooth = sum(vect_N_smooth * Sn))
}
