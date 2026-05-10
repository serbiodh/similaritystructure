#' Estimate probability of effect between two similarity structures
#'
#' @param fitting_A Fitting output for group A.
#' @param fitting_B Fitting output for group B.
#' @param smooth_A Smoothing output for group A.
#' @param smooth_B Smoothing output for group B.
#'
#' @return A probability or a list of probabilities comparing A and B.
#'
#' @export

prob_effect <- function(fitting_A = NULL, fitting_B = NULL, smooth_A = NULL, smooth_B = NULL) {

  # ps from fitting
  if (!is.null(fitting_A)) {

    Param_A <- fitting_A$param_P_N_alpha
    Param_B <- fitting_B$param_P_N_alpha
    Ns_A <- fitting_A$Ns_fitting
    Ns_B <- fitting_B$Ns_fitting

    if (Ns_B > Ns_A) {  # The effect A is greater than the effect B, so we obtain the prob of Expect_N_B in the p(n) for Param_A

      Out <- ps_fitting(Param_A, Ns_B)
      Prob_Effect_fit <- Out$prob

    } else {  # The effect B is greater than the effect A, so we obtain the prob of Expect_N_A in the p(n) for Param_B

      Out <- ps_fitting(Param_B, Ns_A)
      Prob_Effect_fit <- Out$prob

    }

  }

  # ps from smooth
  if(!is.null(smooth_A)) {

    fun_A <- smooth_A$fun_smooth_norm
    fun_B <- smooth_B$fun_smooth_norm
    Nss_A <- round(smooth_A$Ns_smooth)
    Nss_B <- round(smooth_B$Ns_smooth)

    Nfin_A <- utils::tail(smooth_A$vect_N_smooth, 1)
    Nfin_B <- utils::tail(smooth_B$vect_N_smooth, 1)

    if (Nss_B > Nss_A) {  # The effect A is greater than the effect B, so we obtain the prob of Expect_N_B in the p(n) for smooth_A

      if (Nss_B >= Nfin_A) {

        Prob_Effect_smooth <- 0

      } else {

        Prob_Effect_smooth <- sum(fun_A[Nss_B:length(fun_A)])
      }

    } else {  # The effect B is greater than the effect A, so we obtain the prob of Expect_N_A in the p(n) for smooth_B

      if (Nss_A >= Nfin_B) {

        Prob_Effect_smooth <- 0

      } else {

        Prob_Effect_smooth <- sum(fun_B[Nss_A:length(fun_B)])
      }

    }

  }

  # Output
  if (!is.null(fitting_A) & is.null(smooth_A)) {

    return(Prob_Effect_fit)

  } else if (is.null(fitting_A) & !is.null(smooth_A)) {

    return(Prob_Effect_smooth)

  } else {

    return(list(ps_fit = Prob_Effect_fit,
                ps_smooth = Prob_Effect_smooth))

  }

}

