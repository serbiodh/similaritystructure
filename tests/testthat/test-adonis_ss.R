test_that("adonis_ss returns a valid PERMANOVA p-value", {
  skip_if_not_installed("vegan")

  set.seed(123)
  observations <- matrix(rnorm(24), ncol = 2)
  D <- as.matrix(dist(observations))
  rownames(D) <- colnames(D) <- paste0("sample", seq_len(nrow(D)))

  # The custom test receives rows from a common square distance matrix.
  first <- D[1:6, , drop = FALSE]
  second <- D[7:12, , drop = FALSE]

  res <- adonis_ss(first, second,
                   subset_choice = "Female", permutations = 19)

  expect_true(is.numeric(res$p.value))
  expect_length(res$p.value, 1)
  expect_true(is.finite(res$p.value))
  expect_gte(res$p.value, 0)
  expect_lte(res$p.value, 1)
})
