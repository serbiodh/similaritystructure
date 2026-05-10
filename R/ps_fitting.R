#' Tail probability from p(n) distribution
#'
#' Computes the probability that the similarity size is greater than or equal to
#' a given threshold `M`, using the normalized p(n) distribution.
#'
#' @param Param_Final Numeric vector. Normalized fitted parameters.
#' @param M Numeric. Lower threshold for the tail probability.
#'
#' @return A list with one component:
#' \describe{
#'   \item{prob}{Estimated tail probability.}
#' }
#'
#' @keywords internal


ps_fitting <- function(Param_Final, M){

  # This script calculates the prob of N >= M given by p(n)

  # This script applies Bayes's theorem numercally, i.e., as
  # p(N | >alpha) = p(>alpha | N)/sum_N( p(>alpha | N) ), the sum of
  # p(>alpha | N) on N is obtained numerically

  # Most of the times the similarity structure cannot be calculated
  # for N small enough, so the sum_N( p(>alpha | N) ) goes from
  # m (>1) to inf.

  # Parameters
  a = Param_Final[1]; b = Param_Final[2]; c = Param_Final[3];
  d = Param_Final[4]; e = Param_Final[5]

  m = M

  # sum( exp(-d*N) ) with N from m to inf:
  t0_inf = (exp(-d*(m - 1)))/(-1+exp(d));
  # sum( N*exp(-d*N) ) with N from m to inf:
  t1_inf = ((m*exp(d) - m + 1)/((1 - exp(d))^2))*exp(-d*(m-1))
  # sum( N^2*exp(-d*N) ) with N from m to inf:
  t2_inf = ( -(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) - ((m - 1)*exp(-d*m))/((-1 + exp(-d))^2) + 2*(-m + m*exp(-d) - exp(-d))*exp(-d*m)/((-1 + exp(-d))^3) )*exp(-d)
  # sum( N^3*exp(-d*N) ) with N from m to inf:
  t3_inf = ( -(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) + ( -((m^2)*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) - 2*(m*(m - 1)*exp(d)*exp(-d*m))/((-1 + exp(-d))^2) + (m*(-m + m*exp(-d) - exp(-d))*exp(2*d)*exp(-d*m))/((-1 + exp(-d))^2) + 4*(m*(-m + m*exp(-d) - exp(-d))*exp(d)*exp(-d*m))/((-1 + exp(-d))^3) + 4*((m-1)*exp(-d*m))/((-1 + exp(-d))^3) -6*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^4) )*exp(-d) - ((m-1)*exp(-d*m))/((-1 + exp(-d))^2) + 2*((-m + m*exp(-d) - exp(-d))*exp(-d*m))/((-1 + exp(-d))^3) )*exp(-d)

  p = a*t0_inf + b*t1_inf + c*t2_inf + e*t3_inf

  return(list(prob = p))
}

