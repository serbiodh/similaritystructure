#' Apply a statistical test to two subsamples
#'
#' Draws a subsample of size `nsub` from each of two samples and applies a
#' statistical test to the resulting subsamples. The function supports vectors
#' and matrices. If matrices are provided, rows are sampled.
#'
#' @param n1 First sample. It must be either a vector or a matrix.
#' @param n2 Second sample. It must be of the same type as `n1`.
#' @param nsub Integer. Size of the subsamples to draw from `n1` and `n2`.
#' @param test Statistical test to use. Either `"t-test"`, `"ks-test"`, or a
#'   custom function. Custom functions must return an object containing a
#'   `p.value` element.
#' @param test_args Optional list. Additional arguments passed to the test
#'   function.
#' @param replace Logical. Should subsampling be performed with replacement?
#'
#' @return A numeric value corresponding to the p-value of the selected test.
#'
#' @keywords internal

subsampling_testing <- function(n1, n2, nsub, test = "t-test", test_args = NULL, replace = TRUE) {

  # Subsampling depending on whether n1 and n2 are vectors or matrices
  if (is.vector(n1) && is.vector(n2)){
    n1sub <- sample(n1, nsub, replace = replace)
    n2sub <- sample(n2, nsub, replace = replace)

  } else if (is.matrix(n1) && is.matrix(n2)) {
    # Obtain two random permutations for extracting the random subsamples of sizes n2 and n1
    ind_1 <- sample(nrow(n1), nsub, replace = replace)
    ind_2 <- sample(nrow(n2), nsub, replace = replace)

    # Obtain the random subsamples of size vector_N[i] from Sample_2 and size vector_N[i] from Sample_1, respectively
    n1sub <- n1[ind_1, ]
    n2sub <- n2[ind_2, ]

  } else {
    stop('Mismatching or wrong data dimensions (> 2)')
  }

  # Hypothesis testing
  res <- if (is.character(test)) {
    switch(
      test,
      "t-test"  = do.call(stats::t.test, c(list(n1sub, n2sub), test_args)),
      "ks-test" = suppressWarnings(do.call(stats::ks.test, c(list(n1sub, n2sub), test_args))),
      stop("Unknown test. Use 't-test', 'ks-test', or a function.")
    )

  } else if (is.function(test)) {
    do.call(test, c(list(n1sub, n2sub), test_args))

  } else {
    stop("'test' must be either a character string or a function.")
  }

  if (is.null(res$p.value)) {
    stop("The test function must return an object with a 'p.value' element.")
  }

  res$p.value

}
