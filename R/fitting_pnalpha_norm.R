#' Normalize a fitted P(N|>alpha) curve
#'
#' Converts fitted P(N|>alpha) parameters into normalized parameters so that the
#' fitted curve can be interpreted as a probability distribution over subsample
#' sizes.
#'
#' @param Param_0 Numeric vector. Raw fitted parameters, usually obtained from
#'   `fitting_pnalpha()`.
#' @param vect_N Numeric vector. Subsample sizes over which the normalized
#'   P(N|>alpha) distribution is evaluated.
#'
#' @return A list with two components:
#' \describe{
#'   \item{Params}{Normalized parameter vector.}
#'   \item{P_N_alpha}{Normalized probability values evaluated over `vect_N`.}
#' }
#'
#' @keywords internal


fitting_pnalpha_norm = function(Param_0, vect_N){

  # This script normalizes the experimental p(N | >alpha) to get the final p(N | >alpha).
  # The sum_N( p(>alpha | N) ) goes from m to inf.
  # WARNING! This is done because this final p(N | >alpha) is the probability distribution,
  # i.e., the theoretical law of the distribution of the probability

  # Bayes theorem
  # Parameters
  a = Param_0[1]; b = Param_0[2]; c = Param_0[3];
  d = Param_0[4]; e = Param_0[5]

  # First value of N in the series
  m = vect_N[1]

  # sum( exp(-d*N) ) with N from m to inf:
  t0_inf = (exp(-d*(m - 1)))/(-1+exp(d));
  # sum( N*exp(-d*N) ) with N from m to inf:
  t1_inf = ((m*exp(d) - m + 1)/((1 - exp(d))^2))*exp(-d*(m-1))
  # sum( N^2*exp(-d*N) ) with N from m to inf:
  t2_inf = ( -(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) - ((m - 1)*exp(-d*m))/((-1 + exp(-d))^2) + 2*(-m + m*exp(-d) - exp(-d))*exp(-d*m)/((-1 + exp(-d))^3) )*exp(-d)
  # sum( N^3*exp(-d*N) ) with N from m to inf:
  t3_inf = ( -(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) + ( -((m^2)*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) - 2*(m*(m - 1)*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) + (m*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) + 4*(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^3) + 4*((m-1)*exp(-d*m))/((-1 + exp(-d))^3) -6*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^4) )*exp(-d) - ((m-1)*exp(-d*m))/((-1 + exp(-d))^2) + 2*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^3) )*exp(-d)

  P_alpha = a*t0_inf + b*t1_inf + c*t2_inf + e*t3_inf;  # Eq. [S4] in the ms

  Param = c(a/P_alpha, b/P_alpha, c/P_alpha, d, e/P_alpha)  # Eq. [12]

  P_N_alpha = (Param[1] + Param[2]*vect_N + Param[3]*(vect_N^2) + Param[5]*(vect_N^3))*exp(-Param[4]*vect_N)

  return(list(Params = Param, P_N_alpha = P_N_alpha))
}
