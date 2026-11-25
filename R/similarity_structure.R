#' Estimate similarity structure between two samples
#'
#' @param n1 First sample. It can be a vector or a matrix.
#' @param n2 Second sample. It can be a vector or a matrix.
#' @param N_init Initial subsample size (default 2).
#' @param N_fin Final subsample size.
#' @param num_N Number of approximately equidistant integer subsample sizes (default 60).
#' @param num_repet Number of repetitions for each subsample size (default 300).
#' @param replace Logical. Should sampling be done with replacement?
#' @param test Statistical test to use. Either `"t-test"`, `"ks-test"`, or a custom function.
#' @param test_args Optional list of additional arguments passed to the test.
#' @param alpha Significance level.
#' @param seed Seed used by the parallel backend.
#' @param sav Logical. Should the scaled beta(n) object be saved?
#' @param sav_file Optional file path for saving the empirical beta(n) object.
#' @param plotting Logical. Should the fitted and smoothed curves be plotted?
#' @param title_plot Optional title for the plot.
#'
#' @details The empirical probability of non-rejection is the number of
#'   non-rejections divided by the number of repetitions. The curve used for
#'   parametric fitting is scaled by its first value. The smoothing approach
#'   models the original counts with a binomial GAM (logit link, unrestricted
#'   thin plate spline, REML); predictions are rescaled for plotting. The similarity structure is the
#'   normalized distribution over integer subsample sizes, not this scaled curve.
#'   Smoothing normalizes over the evaluated finite interval; parametric fitting
#'   normalizes over the infinite support starting at `N_init`. The returned
#'   fitted vector covers only the displayed range and need not sum to one.
#'   Compare structures calculated using the same test, significance level and
#'   compatible sample-size support. See [prob_effect()] to compare effects.
#'
#' @return A list with two components:
#' \describe{
#'   \item{scaled_bn}{Data frame with `N`, `num_nosig` and `scaled` (count divided by the first count), with the original number of tests stored in the `num_repet` attribute.}
#'   \item{simil_struct_Ns}{List with `smooth`, `fitting` and `Ns_diff` (signed percentage difference). `smooth` contains `vect_N_smooth`, `fun_smooth`, `Sn` and `Ns_smooth`. `fitting` contains `param_0`, `vect_N_fit`, `fun_fit`, `Gof_fit`, `Scale`, `param_Sn`, `Sn` and `Ns_fitting`, or NULL if fitting fails.}
#' }
#' @examples
#' \dontrun{
#' set.seed(123)
#' A <- rnorm(1000)
#' B <- rnorm(1000, 0.8)
#' res <- similarity_structure(A, B, N_fin = 150, test = "t-test", seed = 123)
#' res$simil_struct_Ns$smooth$Sn
#' res$simil_struct_Ns$smooth$Ns_smooth
#' plot(res$simil_struct_Ns$smooth$vect_N_smooth,
#'      res$simil_struct_Ns$smooth$Sn, type = "l")
#' if (!is.null(res$simil_struct_Ns$fitting)) {
#'   res$simil_struct_Ns$fitting$Sn
#'   res$simil_struct_Ns$fitting$Ns_fitting
#' }
#' 
#' }
#'
#' @export

similarity_structure <- function(n1, n2, N_init = 2, N_fin, num_N = 60, num_repet = 300,
                                 replace = TRUE,
                                 test, test_args = NULL,
                                 alpha = 0.05, seed = TRUE,
                                 sav = FALSE, sav_file = NULL,
                                 plotting = TRUE, title_plot = NULL) {


  ###############
  # Computation #
  ###############

  # Validate the interval and compute approximately equidistant integer sizes
  is_whole <- function(x) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x == floor(x)
  if (!is_whole(N_init) || !is_whole(N_fin) ||
      N_init < 2L || N_fin < N_init) {
    stop("N_init and N_fin must be integers with 2 <= N_init <= N_fin.")
  }
  if (!is_whole(num_N) || num_N < 3L ||
      !is_whole(num_repet) || num_repet < 1L) {
    stop("num_N must be at least 3 and num_repet must be positive (integers).")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a number strictly between 0 and 1.")
  }
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    stop("replace must be TRUE or FALSE.")
  }
  if (!replace && N_fin > min(NROW(n1), NROW(n2))) {
    stop("N_fin exceeds the smaller sample size for sampling without replacement.")
  }
  vector_Nsub <- unique(round(seq(N_init, N_fin, length.out = num_N)))
  if (length(vector_Nsub) < 3L) {
    stop("At least three distinct subsample sizes are required for smoothing.")
  }

  # Number of non significant tests per subsample size
  propsig <- sn_norm(vector_Nsub, n1, n2, num_repet, test = test, test_args = test_args, alpha = alpha,
                            seed = seed, replace = replace)

  # Save data frame with scaled beta(n) in case the next step fails
  if (sav == T) {

    if (is.null(sav_file) == T) {
      saveRDS(propsig, "scaled_bn.RData")

    } else {
      saveRDS(propsig, file = sav_file)

    }
  }

  # Similarity structure and size
  simil_struct <- smooth_fitting_Ns(N_init, N_fin, propsig,
                                    plotting = plotting, title_plot = title_plot,
                                    num_repet = num_repet)

  return(list(scaled_bn = propsig,
              simil_struct_Ns = simil_struct))

}

