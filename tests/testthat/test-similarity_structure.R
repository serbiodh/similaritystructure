test_that("similarity_structure returns expected structure", {
  set.seed(12345)
  A <- rnorm(1000)
  B <- rnorm(1000, .8)
  res <- similarity_structure(A, B, N_init = 2, N_fin = 150,
                              num_N = 20, num_repet = 40,
                              test = "t-test", seed = 12345,
                              plotting = FALSE)
  expect_named(res, c("scaled_bn", "simil_struct_Ns"))
  expect_named(res$scaled_bn, c("N", "num_nosig", "scaled"))
  expect_equal(range(res$scaled_bn$N), c(2, 150))
  expect_equal(nrow(res$scaled_bn), 20L)
  smooth <- res$simil_struct_Ns$smooth
  expect_equal(sum(smooth$Sn), 1)
  expect_equal(smooth$Ns_smooth, sum(smooth$Sn * smooth$vect_N_smooth))
})

test_that("invalid inputs fail clearly", {
  expect_error(similarity_structure(1:10, 1:10, N_fin = 11,
                                    replace = FALSE, test = "t-test"),
               "exceeds")
  expect_error(similarity_structure(1:10, 1:10, N_init = 1,
                                    N_fin = 10, test = "t-test"), "N_init")
})
