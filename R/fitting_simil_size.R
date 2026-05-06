#' Estimate similarity size from fitted parameters
#'
#' Computes the expected similarity size from the normalized parameters of the
#' fitted P(N|>alpha) distribution.
#'
#' @param Param_Final Numeric vector. Normalized fitted parameters, usually
#'   returned by `fitting_pnalpha_norm()`.
#' @param N_init Numeric. Initial subsample size used as the lower bound of the
#'   similarity-size calculation.
#'
#' @return A list with one component:
#' \describe{
#'   \item{Simil_Size}{Estimated similarity size.}
#' }
#'
#' @keywords internal


fitting_simil_size = function(Param_Final, N_init){

  # This script calculates the expected size <N> by infinite series, according to the equations in the main text and the Supplementary
  # Parameters
  a = Param_Final[1]; b = Param_Final[2]; c = Param_Final[3];
  d = Param_Final[4]; e = Param_Final[5]

  # First value of N in the series
  m = N_init

  # sum( N*exp(-d*N) ) with N from m to inf:
  t1_m_inf = ((m*exp(d) - m + 1)/((1 - exp(d))^2))*exp(-d*(m-1))
  # sum( N^2*exp(-d*N) ) with N from m to inf:
  t2_m_inf = ( -(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) - ((m - 1)*exp(-d*m))/((-1 + exp(-d))^2) + 2*(-m + m*exp(-d) - exp(-d))*exp(-d*m)/((-1 + exp(-d))^3) )*exp(-d)
  # sum( N^3*exp(-d*N) ) with N from m to inf:
  t3_m_inf = ( -(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) + ( -((m^2)*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) - 2*(m*(m - 1)*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) + (m*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) + 4*(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^3) + 4*((m-1)*exp(-d*m))/((-1 + exp(-d))^3) - 6*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^4) )*exp(-d) - ((m-1)*exp(-d*m))/((-1 + exp(-d))^2) + 2*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^3) )*exp(-d)
  # sum( N^4*exp(-d*N) ) with N from 1 to inf:
  t4_m_inf = ( - (m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) +        (  -((m^2)*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2)    - 2*(m*(m - 1)*exp(d)*exp(-d*m))/((-1 + exp(-d))^2)    + (m*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2)    + 4*(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^3)    + 4*((m-1)*exp(-d*m))/((-1 + exp(-d))^3)    - 6*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^4) )*exp(-d)     + ( -2*((m^2)*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) - 4*(m*(m - 1)*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) + 2*(m*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) + 8*(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^3) + ( -((m^3)*(-m + m*exp(-d) - exp(-d))*exp(3*d)*exp(-d*m))/((-1 + exp(-d))^2) - 3*( (m^2)*(m - 1)*exp(2*d)*exp(-d*m) )/((-1 + exp(-d))^2) + (3*(m^2)*(-m + m*exp(-d) - exp(-d))*exp(3*d)*exp(-d*m))/((-1 + exp(-d))^2) + 6*((m^2)*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^3) + 3*(m*(m - 1)*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) - (2*m*(-m + m*exp(-d) - exp(-d))*exp(3*d)*exp(-d*m))/((-1 + exp(-d))^2) + 12*(m*(m - 1)*exp(d)*exp(-d*m))/((-1 + exp(-d))^3) - 6*(m*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^3) - 18*(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^4) - 18*((m - 1)*exp(-d*m))/((-1 + exp(-d))^4) + 24*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^5) )*exp(-d)     + 8*((m - 1)*exp(-d*m))/((-1 + exp(-d))^3) - 12*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^4) )*exp(-d)    - ((m - 1)*exp(-d*m))/((-1 + exp(-d))^2) + 2*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^3)  )*exp(-d)

  Similarity_Size = a*t1_m_inf + b*t2_m_inf + c*t3_m_inf + e*t4_m_inf

  return(list(Simil_Size = Similarity_Size))
}

