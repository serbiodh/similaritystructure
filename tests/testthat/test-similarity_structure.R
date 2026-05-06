test_that("similarity_structure returns expected structure", {

  seed <- 12345
  set.seed(seed)

  n1 <- rnorm(20000, 0, 1)
  n2 <- n1 + 0.8

  res <- similarity_structure(
    n1 = n1,
    n2 = n2,
    N_init = 2,
    N_fin = round((80/0.8^2)),
    num_N = 60,
    num_repet = 300,
    test = "t-test",
    alpha = 0.05,
    seed <- seed,
    plotting = TRUE
  )

  expect_type(res, "list")
  expect_named(res, c("empiric_PNalpha", "simil_struct_Ns"))
  expect_true("N" %in% names(res$empiric_PNalpha))
  expect_true("scaled" %in% names(res$empiric_PNalpha))
})
