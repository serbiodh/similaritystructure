#' Calculate the similar effect probability (pS)
#'
#' @param fitting_A,fitting_B Fitting components from two similarity_structure()
#'   results, obtained with `$simil_struct_Ns$fitting`.
#' @param smooth_A,smooth_B Smoothing components from two similarity_structure()
#'   results, obtained with `$simil_struct_Ns$smooth`.
#' @details Computes the minimum of the two probabilities of a subsample size
#'   being at least as large as the similarity size of the other comparison.
#'   Integer sample sizes are compared to the unrounded similarity sizes.
#'   Use the same test, significance level and compatible sample-size supports.
#'   Smoothing uses finite support; fitting uses the infinite parametric tail.
#'   pS is a descriptive comparison probability, not a hypothesis-test p-value.
#'   A small value indicates dissimilar effects according to this criterion.
#'   A threshold such as 0.05 is a user-defined decision criterion.
#' @return A numeric probability when one method is supplied, or a list
#'   containing `ps_fit` and `ps_smooth` when both methods are supplied.
#' @seealso [similarity_structure()]
#' @examples
#' \dontrun{
#' set.seed(123)
#' A <- rnorm(1000)
#' B <- rnorm(1000, 0.5)
#' C <- rnorm(1000, 0.9)
#' AB <- similarity_structure(A, B, N_fin = 300, test = "t-test", seed = 123)
#' AC <- similarity_structure(A, C, N_fin = 300, test = "t-test", seed = 124)
#' prob_effect(smooth_A = AB$simil_struct_Ns$smooth,
#'             smooth_B = AC$simil_struct_Ns$smooth)
#' if (!is.null(AB$simil_struct_Ns$fitting) &&
#'     !is.null(AC$simil_struct_Ns$fitting)) {
#'   prob_effect(fitting_A = AB$simil_struct_Ns$fitting,
#'               fitting_B = AC$simil_struct_Ns$fitting)
#' }
#' }
#' @export
prob_effect <- function(fitting_A = NULL, fitting_B = NULL,
                        smooth_A = NULL, smooth_B = NULL) {

  if (xor(is.null(fitting_A), is.null(fitting_B)) ||
      xor(is.null(smooth_A), is.null(smooth_B))) {
    stop("Supply both A and B for each requested method.")
  }
  if (is.null(fitting_A) && is.null(smooth_A)) {
    stop("Supply a pair of fitting or smoothing results.")
  }

  # pS from fitting
  if (!is.null(fitting_A)) {
    Param_A <- fitting_A$param_Sn
    Param_B <- fitting_B$param_Sn
    # Accept objects saved with older versions of the package
    if (is.null(Param_A)) Param_A <- fitting_A$param_P_N_alpha
    if (is.null(Param_B)) Param_B <- fitting_B$param_P_N_alpha
    Ns_A <- fitting_A$Ns_fitting
    Ns_B <- fitting_B$Ns_fitting
    m_A <- fitting_A$vect_N_fit[1L]
    m_B <- fitting_B$vect_N_fit[1L]
    if (length(Param_A) < 5L || length(Param_B) < 5L ||
        any(!is.finite(c(Param_A, Param_B))) ||
        Param_A[4L] <= 0 || Param_B[4L] <= 0 ||
        length(Ns_A) != 1L || length(Ns_B) != 1L ||
        length(m_A) != 1L || length(m_B) != 1L ||
        any(!is.finite(c(Ns_A, Ns_B, m_A, m_B))) ||
        Ns_A < m_A || Ns_B < m_B) {
      stop("Invalid fitting components.")
    }
    p_A <- if (Ns_B <= m_A) 1 else ps_fitting(Param_A, Ns_B)$prob
    p_B <- if (Ns_A <= m_B) 1 else ps_fitting(Param_B, Ns_A)$prob
    if (any(!is.finite(c(p_A, p_B))) ||
        any(c(p_A, p_B) < -1e-8 | c(p_A, p_B) > 1 + 1e-8)) {
      stop("Invalid parametric tail probability.")
    }
    Prob_Effect_fit <- max(0, min(1, p_A, p_B))
  }

  # pS from smoothing
  if (!is.null(smooth_A)) {
    fun_A <- smooth_A$Sn
    fun_B <- smooth_B$Sn
    if (is.null(fun_A)) fun_A <- smooth_A$fun_smooth_norm
    if (is.null(fun_B)) fun_B <- smooth_B$fun_smooth_norm
    n_A <- smooth_A$vect_N_smooth
    n_B <- smooth_B$vect_N_smooth
    Ns_A <- smooth_A$Ns_smooth
    Ns_B <- smooth_B$Ns_smooth
    if (!length(fun_A) || !length(fun_B) ||
        length(fun_A) != length(n_A) || length(fun_B) != length(n_B) ||
        length(Ns_A) != 1L || length(Ns_B) != 1L ||
        any(!is.finite(c(fun_A, fun_B, n_A, n_B, Ns_A, Ns_B))) ||
        any(c(fun_A, fun_B) < 0) ||
        abs(sum(fun_A) - 1) > 1e-8 || abs(sum(fun_B) - 1) > 1e-8) {
      stop("Invalid smoothing components.")
    }

    # Index by the actual subsample sizes rather than by vector position.
    p_A <- sum(fun_A[n_A >= Ns_B])
    p_B <- sum(fun_B[n_B >= Ns_A])
    Prob_Effect_smooth <- min(p_A, p_B)
  }

  if (is.null(fitting_A)) return(Prob_Effect_smooth)
  if (is.null(smooth_A)) return(Prob_Effect_fit)
  return(list(ps_fit = Prob_Effect_fit, ps_smooth = Prob_Effect_smooth))
}
