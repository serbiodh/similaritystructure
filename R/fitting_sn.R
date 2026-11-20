#' Select a valid parametric fit for beta(n)
#'
#' Fits and filters candidate parametric beta(n) models. The function calls
#' `fitting_basic()` and then evaluates the fitted parameters to discard
#' unsuitable fits, such as fits with invalid decay parameters or problematic
#' curve shapes.
#'
#' @param vector_N Numeric vector. Subsample sizes used to estimate the empirical
#'   scaled beta(n) curve.
#' @param PROB_N Numeric vector. Empirical scaled beta(n) values corresponding
#'   to `vector_N`.
#' @param Scale Numeric. Scaling factor passed to `fitting_basic()`.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Aprox}{Selected approximation family.}
#'   \item{Param}{Selected parameter vector or matrix.}
#'   \item{Gof}{Goodness-of-fit value for the selected model.}
#'   \item{vect_N_fit}{Sequence of subsample sizes used to evaluate the fitted
#'   curve.}
#'   \item{func_fit}{Values of the fitted beta(n) curve over `vect_N_fit`.}
#' }
#'
#' @keywords internal

fitting_sn = function(vector_N, PROB_N, Scale){

  # This script provides the fitting of a polynomial multiplied by exp(-d*N) to the probabilites beta(n) vs. N's in 'vector_N'

  # To substitute in the fitted function
  vect_N = seq(vector_N[1], 4*vector_N[length(vector_N)], by=1)  # For the fitting the final x-value is four times the largest evaluated size to ensure that the fitted function approaches to zero properly

  # Obtaining p(>alpha | n1* ^ n2) from the fitting (n1* is a fixed value for n1; n2 remains variable)
  # Output_N = Fitting_Prob_Cond_Basic(vector_N, PROB_N)
  Output_N = fitting_basic(vector_N, PROB_N, Scale)

  Aprox_0 = Output_N[[1]]
  Param_0 = Output_N[[2]]
  GOF_0 = Output_N[[3]]
  if (length(GOF_0) == 0L) return(list(Success = 0))

  # Discard those cases with min < 0 and max > 1:
  minim = c()
  maxim = c()
  d_sign = c()
  # local_min = zeros(length(GOF_0), 1)
  local_min = rep(0, length(GOF_0))
  for (i in 1:length(GOF_0)){
    a = Param_0[i, 1]; b = Param_0[i, 2]; c = Param_0[i, 3]; d = Param_0[i, 4]; e = Param_0[i, 5]
    f1 = (a + b*vect_N + c*(vect_N^2) + e*(vect_N^3))*exp(-d*vect_N)
    minim = rbind(minim, min(f1))
    maxim = rbind(maxim, max(f1))
    d_sign = rbind(d_sign, sign(d))

    # local minima
    cc = which(diff(sign(diff(f1)))==2)+1
    if (length(cc) == 1) {  # there is a local minimum
      if (f1[cc] < 0.9) { # the local minimum is not small so we discard this case
        local_min[i] = 1
      }
    }
  }
  ind_pos = d_sign == 1 & minim >= 0 & maxim <= 1.05 & local_min == 0
  if (!any(ind_pos)){ # No fitting fulfilling the three conditions (all FALSE)
    # Relax the conditions - maxim <= 1.1 (numerical)
    ind_pos = d_sign == 1 & minim >= 0 & maxim <= 1.1 & local_min == 0
    if (!any(ind_pos)){ # Relax more the conditions; allow local minima
      ind_pos = d_sign == 1 & minim >= 0 & maxim <= 1.1
    }
  }


  if (!any(ind_pos)) { # No fitting --> stop
    return(list(vector_N = vector_N, PROB_N = PROB_N, Success = 0))
  }
  Aprox = Aprox_0[ind_pos]
  Param = Param_0[ind_pos, ]
  GOF = GOF_0[ind_pos]

  if (!is.null(Aprox)){   # Some fitting(s) was(were) successful; we calculate p(>alpha | n1*) for n1*

    # # Take the prob with the minimum GOF
    ind_min = which.min(GOF)

    GOF_Final = GOF[ind_min]
    Aprox_Final = Aprox[ind_min]
    if (is.matrix(Param)){
      Param_Final = Param[ind_min, ]
    } else {
      Param_Final = Param
    }

  } else{  # No fitting was successful
    return(list(vector_N = vector_N, PROB_N = PROB_N, Success = 0))
  }

  # Fitted function:
  a = Param_Final[1]; b = Param_Final[2]; c = Param_Final[3]; d = Param_Final[4]; e = Param_Final[5]
  f1 = (a + b*vect_N + c*(vect_N^2) + e*(vect_N^3))*exp(-d*vect_N)

  return(list(Aprox = Aprox_Final, Param = Param_Final, Gof = GOF_Final, vect_N_fit = vect_N, func_fit = f1))
}
