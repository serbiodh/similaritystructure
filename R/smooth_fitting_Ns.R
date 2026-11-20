#' Smooth and fit a scaled beta(n) curve
#'
#' Takes an empirical scaled beta(n) curve and estimates the corresponding
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
#' @param propsig Data frame or `data.table`. Empirical scaled beta(n) curve, usually
#'   returned by `sn_norm()`. It must contain at least the columns `N` and
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
#'   goodness-of-fit value, normalized parameters, similarity structure S(n)
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

  smooth_formula <- stats::as.formula(
    sprintf("scaled ~ s(N, bs = 'tp', k = %d)", min(10L, nrow(propsig)))
  )

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
  if (!is.finite(sum(fun_smooth)) || sum(fun_smooth) <= 0) {
    stop("The smoothed curve cannot be normalized.")
  }
  Sn_smooth <- fun_smooth/sum(fun_smooth)

  # Similarity size
  Ns_smooth <- sum(vect_N_smooth*Sn_smooth)/sum(Sn_smooth)


  ###########
  # Fitting #
  ###########

  fitting <- tryCatch({
    Output_0 <- fitting_sn(t(propsig$N), t(propsig$scaled), 1)
    if (is.null(Output_0$Param)) stop("No valid parametric fit.")

    param_0 <- Output_0$Param
    vect_N_fit <- Output_0$vect_N_fit
    fun_fit <- Output_0$func_fit
    Gof_fit <- Output_0$Gof
    Scale <- 1

    # Similarity structure and similarity size
    Output_1 <- fitting_sn_norm(param_0, vect_N_fit)
    param_Sn <- Output_1$Params
    Sn <- Output_1$Sn
    Ns_fitting <- fitting_simil_size(param_Sn, N_init)$Simil_Size
    if (any(!is.finite(param_Sn)) || any(!is.finite(Sn)) ||
        any(Sn < 0) || !is.finite(Ns_fitting) || Ns_fitting < N_init) {
      stop("Invalid normalized parametric fit.")
    }

    list(param_0 = param_0, vect_N_fit = vect_N_fit,
         fun_fit = fun_fit, Gof_fit = Gof_fit, Scale = Scale,
         param_Sn = param_Sn, Sn = Sn, Ns_fitting = Ns_fitting)
  }, error = function(e) {
    warning("Parametric fitting failed: ", conditionMessage(e),
            ". The smoothed result is retained.", call. = FALSE)
    NULL
  })


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
      ylab = "Scaled beta(n)",
      xlab = "Subsample size",
      main = title_plot)
    if (!is.null(fitting)) {
      graphics::lines(fitting$vect_N_fit, fitting$Scale*fitting$fun_fit,
                      type = 'l', col = "darkred", lwd = 1.5)
      graphics::abline(v = fitting$Ns_fitting, col = "darkred")
    }
    graphics::lines(vect_N_smooth, fun_smooth, col = "darkgreen")
    graphics::abline(v = 0, lty = "dashed")
    graphics::abline(v = Ns_smooth, col = "darkgreen")
    if (is.null(fitting)) {
      graphics::legend("topright", "Smooth", col = "darkgreen", lty = "solid")
    } else {
      graphics::legend("topright", c("Fitting", "Smooth"), col = c("darkred", "darkgreen"), lty = "solid")
    }

  }


  ###########
  # Results #
  ###########

  return(list(smooth = list(vect_N_smooth = vect_N_smooth,
                            fun_smooth = fun_smooth,
                            Sn = Sn_smooth,
                            Ns_smooth = Ns_smooth),
              fitting = fitting,
              Ns_diff = if (is.null(fitting)) NA_real_ else
                (Ns_smooth-fitting$Ns_fitting)/fitting$Ns_fitting*100))

}
