#' Fit basic parametric g(n) candidate models
#'
#' Fits several polynomial-exponential candidate models to an empirical g(n)
#' curve using non-linear least squares. The candidate models have the general
#' form of a polynomial in `N` multiplied by an exponential decay term.
#'
#' The function tries different polynomial degrees and initial parameter values,
#' then returns the best valid approximation found among the candidate fits.
#'
#' @param n0 Numeric vector. Subsample sizes.
#' @param Prob Numeric vector. Empirical probabilities or scaled g(n) values
#'   associated with `n0`.
#' @param Scale Numeric. Scaling factor used to build initial parameter values
#'   for the non-linear fits.
#'
#' @return A list with three unnamed elements:
#' \describe{
#'   \item{Aprox}{Integer indicating the selected approximation family.}
#'   \item{Param}{Matrix of fitted parameters.}
#'   \item{GOF}{Goodness-of-fit value, computed from the residual standard
#'   error of the fitted model.}
#' }
#'
#' @keywords internal

fitting_basic = function(n0, Prob, Scale){

  # This fitting does not consider the parameter h for the Aprox = 1 fitting
  # WARNING!! In the function 'nls' (for the fitting), the input should be COLUMNS.
  # Thus, the input variables 'n0' and 'Prob' must be COLUMNS; otherwise, ERROR
  # It should be noted that the ERRORS are not displayed, so be careful!

  # Factor to adjust the initial values of the fitting parameters
  Factor = 1
  if(n0[length(n0)] < 300) {
    Factor = 10
  }

  # Initial parameter values for Aprox = 1
  X01 = rbind(c(Scale*0.9, Scale*Factor*0.1, Factor*0.003),
              c(Scale*0.98, -Scale*Factor*0.0021, Factor*0.0074),
              c(Scale*0.98, Scale*Factor*0.0011, Factor*0.001),
              c(Scale*0.98, Scale*Factor*0.0021, Factor*0.0074),
              c(Scale*0.98, Scale*Factor*0.0001, Factor*0.0054),
              c(Scale*0.98, Scale*Factor*0.001, Factor*0.0004),
              c(Scale*0.98, Scale*Factor*0.0001, Factor*0.0001),
              c(Scale*0.98, Scale*Factor*0.01, Factor*0.0001))

  Aprox = c()
  Param = c()
  GOF = c()
  for (i in c(1:nrow(X01))) {
    skip_to_next <- FALSE
    tryCatch({
      m = minpack.lm::nlsLM(Prob ~ ((a + b*n0)*exp(-d*n0)),
                 start = list(a = X01[i,1], b = X01[i,2], d = X01[i,3]),
                 control = stats::nls.control(maxiter = 300))
    }, error=function(e){skip_to_next <<- TRUE})

    if(skip_to_next) {
      # print("error 1")
      next }
    else { # No error --> save data
      #print("saving")
      # Adjusted parameters:
      Aprox = rbind(Aprox, 1)
      Param0 = unname(m$m$getPars())
      Param = rbind(Param, c(Param0[1], Param0[2], 0, Param0[3], 0, 0, 0, 0))
      GOF = rbind(GOF, stats::sigma(m))
    }
  }

  # Initial parameter values for Aprox = 2
  X02 = rbind(c(Scale*0.95, Scale*Factor*0.00038, Scale*(Factor^2)*6.81e-8, Factor*0.00062),
              c(Scale*1.95, Scale*Factor*0.000038, Scale*(Factor^2)*0.9, Factor*0.0012),
              c(Scale*0.95, Scale*Factor*0.00038, Scale*(Factor^2)*0.1, Factor*0.0002),
              c(-Scale*5, Scale*Factor*4, Scale*(Factor^2)*2e-3, Factor*0.1))

  for (i in c(1:nrow(X02))) {
    skip_to_next <- FALSE
    tryCatch({
      m = minpack.lm::nlsLM(Prob~(a + b*n0 + c*(n0^2))*exp(-d*n0),
                start=list(a = X02[i,1], b = X02[i,2], c = X02[i,3], d = X02[i,4]),
                control=stats::nls.control( maxiter = 300))
    }, error=function(e){skip_to_next <<- TRUE})

    if(skip_to_next) {
      # print("error 2")
      next }
    else { # No error --> save data
      #print("saving")
      # Adjusted parameters:
      Aprox = rbind(Aprox, 2)
      # Param0 = unname(summary(m)$coefficients[,1])
      Param0 = unname(m$m$getPars())
      Param = rbind(Param, c(Param0[1], Param0[2], Param0[3], Param0[4], 0, 0, 0, 0))
      GOF = rbind(GOF, stats::sigma(m))
    }
  }

  # Initial parameter values for Aprox = 3
  X03 = rbind(c(Scale*1, Scale*Factor*0.01, Scale*(Factor^2)*3e-5, Factor*0.07, Scale*(Factor^3)*1.8e-8),
              c(Scale*1.01, Scale*Factor*0.00474, Scale*(Factor^2)*3.33e-7, Factor*0.02, Scale*(Factor^3)*1.8e-7),
              c(Scale*0.98, Scale*Factor*7.96e-3, -Scale*(Factor^2)*2.55e-6, Factor*0.0096, Scale*(Factor^3)*1.04e-8),
              c(Scale*0.98, Scale*Factor*7.96e-3, -Scale*(Factor^2)*2.55e-6, Factor*0.01, Scale*(Factor^3)*1.04e-8),
              c(Scale*0.94, Scale*Factor*0.0014, Scale*(Factor^2)*1e-5, Factor*0.01, Scale*(Factor^3)*2.802e-7),
              c(Scale*0.95, Scale*Factor*0.0003, Scale*(Factor^2)*0.00002, Factor*0.0017, Scale*(Factor^3)*2.802e-7),
              c(Scale*0.95, Scale*Factor*0.0003, Scale*(Factor^2)*0.00002, Factor*0.0027, Scale*(Factor^3)*2.802e-7),
              c(Scale*0.95, Scale*Factor*0.0003, Scale*(Factor^2)*0.00002, Factor*0.0037, Scale*(Factor^3)*2.802e-7),
              c(Scale*1, -Scale*Factor*0.004, Scale*(Factor^2)*0.00001, Factor*0.0045, -Scale*(Factor^3)*2.802e-6),
              c(Scale*1.04, -Scale*Factor*0.00114, Scale*(Factor^2)*0.00001355, Factor*0.01, -Scale*(Factor^3)*2.802e-8),
              c(Scale*1.04, -Scale*Factor*0.00114, Scale*(Factor^2)*0.00001355, Factor*0.0059, -Scale*(Factor^3)*2.802e-8),
              c(Scale*0.96, Scale*Factor*4.78e-4, Scale*(Factor^2)*1.89e-7, Factor*7e-4, -Scale*(Factor^3)*4e-11),
              c(Scale*1.01, Scale*Factor*0.00474, Scale*(Factor^2)*3.33e-7, Factor*0.02, Scale*(Factor^3)*1.8e-7),
              c(Scale*1.09, -Scale*Factor*0.0291, Scale*(Factor^2)*1.82e-4, Factor*0.05, -Scale*(Factor^3)*1e-8))

  for (i in c(1:nrow(X03))) {
    skip_to_next <- FALSE
    tryCatch({
      m = minpack.lm::nlsLM(Prob~(a + b*n0 + c*(n0^2) + e*(n0^3))*exp(-d*n0),
                start=list(a = X03[i,1], b = X03[i,2], c = X03[i,3], d = X03[i,4], e = X03[i,5]),
                control=stats::nls.control( maxiter = 300))
    }, error=function(e){skip_to_next <<- TRUE})

    if(skip_to_next) {
      # print("error 3")
      next }
    else { # No error --> save data
      # print("saving")
      # Adjusted parameters:
      Aprox = rbind(Aprox, 3)
      Param0 = unname(m$m$getPars())
      Param = rbind(Param, c(Param0[1], Param0[2], Param0[3], Param0[4], Param0[5], 0, 0, 0))
      GOF = rbind(GOF, stats::sigma(m))
    }
  }
  return(list(Aprox, Param, GOF))
}
