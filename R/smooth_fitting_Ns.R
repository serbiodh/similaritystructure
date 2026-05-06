#' Smooth and fit a P(N|>alpha) curve
#'
#' Takes an empirical scaled P(N|>alpha) curve and estimates the corresponding
#' similarity structure using two approaches: a smooth GAM-based approximation
#' and a parametric fitted approximation.
#'
#' The smooth approach estimates a smoothed curve over a dense sequence of
#' subsample sizes and normalizes it to obtain an expected similarity size. The
#' fitting approach fits a parametric curve, normalizes it, and computes the
#' corresponding expected similarity size.
#'
#' @param N_init Numeric. Initial subsample size.
#' @param N_fin Numeric. Final subsample size.
#' @param propsig Data frame or `data.table`. Empirical P(N|>alpha) curve, usually
#'   returned by `pnalpha_scaled()`. It must contain at least the columns `N` and
#'   `scaled`.
#' @param plotting Logical. Should the empirical, smoothed and fitted curves be
#'   plotted?
#' @param title_plot Optional character string. Plot title.
#'
#' @return A list with three components:
#' \describe{
#'   \item{smooth}{List containing the smooth approximation, normalized smooth
#'   curve, evaluated subsample sizes and smooth similarity size.}
#'   \item{fitting}{List containing the parametric fit, fitted curve,
#'   goodness-of-fit value, normalized parameters, normalized PN-alpha
#'   distribution and fitted similarity size.}
#'   \item{Ns_diff}{Difference between the fitted and smoothed similarity sizes.}
#' }
#'
#' @keywords internal

smooth_fitting_Ns <- function(N_init, N_fin, propsig,
                              plotting = TRUE, title_plot = NULL) {

  #############
  # Smoothing #
  #############

  vect_N_smooth <- seq(N_init, N_fin, 1)

  smooth_formula <- stats::as.formula("scaled ~ s(N, bs = 'tp')")

  smoothdata <- mgcv::gam(
    smooth_formula,
    data = propsig,
    method = "REML"
  )

  fun_smooth <- stats::predict(
    smoothdata,
    newdata = data.frame(N = vect_N_smooth),
    type = "response"
  )

  fun_smooth[fun_smooth < 0] <- 0

  # Similarity structure (normalization)
  fun_smooth_norm <- fun_smooth/sum(fun_smooth)

  # Similarity size
  Ns_smooth <- sum(vect_N_smooth*fun_smooth_norm)/sum(fun_smooth_norm)


  ###########
  # Fitting #
  ###########

  Output_0 <- fitting_pnalpha(t(propsig$N), t(propsig$scaled), 1)
  # Store the output
  param_0 <- Output_0$Param
  vect_N_fit <- Output_0$vect_N_fit
  fun_fit <- Output_0$func_fit
  Gof_fit <- Output_0$Gof
  Scale <- 1

  # Similarity structure (normalization)
  Output_1 <- fitting_pnalpha_norm(param_0, vect_N_fit)
  param_P_N_alpha <- Output_1$Params
  P_N_alpha <- Output_1$P_N_alpha

  # Similarity size
  Output_2 <- fitting_simil_size(param_P_N_alpha, N_init)
  Ns_fitting <- Output_2$Simil_Size


  ############
  # Plotting #
  ############

  if (plotting == TRUE) {

    graphics::par(mar = c(4, 4, 2, 0.3), bg = "gray")
    graphics::plot(
      propsig$N,
      propsig$scaled,
      type = 'l',
      col = "black",lwd=1.5,
      ylab = "P(N|>alfa) scaled",
      xlab = "Subsample size",
      main = title_plot)
    graphics::lines(vect_N_fit, Scale*fun_fit, type = 'l', col = "darkred",lwd=1.5)
    graphics::lines(vect_N_smooth, fun_smooth, col = "darkgreen")
    graphics::abline(v = 0, lty = "dashed")
    graphics::abline(v = Ns_fitting, col = "darkred")
    graphics::abline(v = Ns_smooth, col = "darkgreen")
    graphics::legend("topright", c("Fitting", "Smooth"), col = c("darkred", "darkgreen"), lty = "solid")

  }


  ###########
  # Results #
  ###########

  return(list(smooth = list(vect_N_smooth = vect_N_smooth,
                            fun_smooth = fun_smooth,
                            fun_smooth_norm = fun_smooth_norm,
                            Ns_smooth = Ns_smooth),
              fitting = list(param_0 = param_0,
                             vect_N_fit = vect_N_fit,
                             fun_fit = fun_fit,
                             Gof_fit = Gof_fit,
                             Scale = Scale,
                             param_P_N_alpha = param_P_N_alpha,
                             P_N_alpha = P_N_alpha,
                             Ns_fitting = Ns_fitting),
              Ns_diff = (Ns_smooth-Ns_fitting)/Ns_fitting*100))

}
