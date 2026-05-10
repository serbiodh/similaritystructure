#' Scaled g(n)/K curve
#'
#' Computes the empirical proportion of non-significant tests across a sequence
#' of subsample sizes. For each subsample size, the function repeatedly samples
#' from two input samples, applies a statistical test, and counts how many
#' resulting p-values are greater than a given significance level.
#'
#' The resulting number of non-significant tests is scaled by the value obtained
#' for the first subsample size.
#'
#' @param vector_Nsub Numeric vector. Subsample sizes to evaluate.
#' @param n1 First sample. It can be a vector or a matrix.
#' @param n2 Second sample. It can be a vector or a matrix.
#' @param num_repet Integer. Number of repeated subsampling tests for each
#'   subsample size.
#' @param test Statistical test to use. Either `"t-test"`, `"ks-test"`, or a
#'   custom function.
#' @param test_args Optional list. Additional arguments passed to the test
#'   function.
#' @param alpha Numeric. Significance level used to classify a test as
#'   significant or non-significant.
#' @param seed Integer or logical. Seed passed to `future.apply::future_lapply()`
#'   through the `future.seed` argument.
#' @param replace Logical. Should subsampling be performed with replacement?
#'
#' @return A `data.table` with columns:
#' \describe{
#'   \item{N}{Subsample size.}
#'   \item{num_nosig}{Number of p-values greater than `alpha`.}
#'   \item{scaled}{`num_nosig` divided by its first value.}
#' }
#'
#' @keywords internal

pnalpha_scaled <- function(vector_Nsub, n1, n2, num_repet = 60, test = "t-test", test_args = NULL, alpha = 0.05,
                           seed = 12345, replace = TRUE) {

  ######################
  # Parallel computing #
  ######################

  # Number of cores
  workers <- max(1, parallel::detectCores() - 1)

  # Start multi-core work
  future::plan(future::multisession, workers = workers)

  # Show progress bar
  progressr::with_progress({

    progressr::handlers("txtprogressbar")

    p <- progressr::progressor(along = vector_Nsub) # start progress bar

    ################################
    # Estimate P(N|>alfa) directly #
    ################################
    pnalpha <- future.apply::future_lapply(vector_Nsub,
                                           function(x) {
                                             p(sprintf("nsub = %d", x)) # text of progress bar
                                             # print(x)
                                             pvals <- replicate(
                                               num_repet,
                                               subsampling_testing(n1, n2, nsub = x, test = test, test_args = test_args, replace = replace),
                                               simplify = TRUE
                                             )
                                           },
                                           future.seed = seed
    )

  })

  # Stop multi-core working
  future::plan(future::sequential)

  # Turn into data frame
  names(pnalpha) <- vector_Nsub

  pnalpha_df <- data.table::rbindlist(
    lapply(names(pnalpha), function(N) {
      sublist <- pnalpha[[N]]
      data.table::data.table(
        N = as.numeric(N),
        rep = seq_along(sublist),
        pval = sublist
      )
    })
  )

  ######################################################
  # Get number of non-significant tests at alpha level #
  ######################################################
  propsig <- stats::aggregate(
    pval ~ N,
    data = pnalpha_df,
    FUN = function(x) sum(x > alpha)
  )

  names(propsig)[names(propsig) == "pval"] <- "num_nosig"

  # Scale to 1 so initial values of fitting work later on
  scale <- propsig$num_nosig[1]
  # scale <- propsig$num_nosig %>% sum()

  propsig[["scaled"]] <- propsig[["num_nosig"]]/scale

  return(propsig)

}
