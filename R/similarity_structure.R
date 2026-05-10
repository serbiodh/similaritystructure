#' Estimate similarity structure between two samples
#'
#' @param n1 First sample. It can be a vector or a matrix.
#' @param n2 Second sample. It can be a vector or a matrix.
#' @param N_init Initial subsample size.
#' @param N_fin Final subsample size.
#' @param num_N Number of equidistant subsample sizes to evaluate.
#' @param num_repet Number of repetitions for each subsample size.
#' @param replace Logical. Should sampling be done with replacement?
#' @param test Statistical test to use. Either `"t-test"`, `"ks-test"`, or a custom function.
#' @param test_args Optional list of additional arguments passed to the test.
#' @param alpha Significance level.
#' @param seed Seed used by the parallel backend.
#' @param sav Logical. Should the scaled g(n) object be saved?
#' @param sav_file Optional file path for saving the empirical PN-alpha object.
#' @param plotting Logical. Should the fitted and smoothed curves be plotted?
#' @param title_plot Optional title for the plot.
#'
#' @return A list with two components:
#' \describe{
#'   \item{scaled_gn}{Scaled g(n) curve.}
#'   \item{simil_struct_Ns}{Estimated similarity structure.}
#' }
#'
#' @export

similarity_structure <- function(n1, n2, N_init, N_fin, num_N, num_repet,
                                 replace = TRUE,
                                 test, test_args = NULL,
                                 alpha = 0.05, seed = TRUE,
                                 sav = FALSE, sav_file = NULL,
                                 plotting = TRUE, title_plot = NULL) {


  ###############
  # Computation #
  ###############

  # Vector of subsample sizes
  vector_Nsub <- unique(round(seq(N_init, N_fin, N_fin / num_N)))

  # Number of non significant tests per subsample size
  propsig <- pnalpha_scaled(vector_Nsub, n1, n2, num_repet, test = test, test_args = test_args, alpha = alpha,
                            seed = seed, replace = replace)

  # Save data frame with scaled g(n) in case the next step fails
  if (sav == T) {

    if (is.null(sav_file) == T) {
      saveRDS(propsig, "scaled_gn.RData")

    } else {
      saveRDS(propsig, file = sav_file)

    }
  }

  # Similarity structure and size
  simil_struct <- smooth_fitting_Ns(N_init, N_fin, propsig,
                                    plotting = plotting, title_plot = title_plot)

  return(list(scaled_gn = propsig,
              simil_struct_Ns = simil_struct))

}

