test_that("subsampling_testing passes arguments to a custom test", {
  custom_test <- function(a, b, chosen_p) {
    list(p.value = chosen_p)
  }

  res <- similaritystructure:::subsampling_testing(
    1:4, 2:5, nsub = 10,
    replace = TRUE, test = custom_test,
    test_args = list(chosen_p = 0.3)
  )
  expect_equal(res, 0.3)

  expect_error(
    similaritystructure:::subsampling_testing(
      1:4, 2:5, nsub = 4,
      test = function(a, b) list(p.value = 1.5)
    ),
    "finite numeric"
  )
})
