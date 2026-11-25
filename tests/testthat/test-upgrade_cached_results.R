test_that("the updater preserves counts, support and parametric fitting", {
  env <- new.env()
  sys.source(system.file("scripts", "upgrade_cached_results.R",
                         package = "similaritystructure"), envir = env)
  n <- seq(2, 100, by = 2)
  K <- 300L
  g <- round(K * .9 * exp(-(n - 2) / 20))
  empirical <- data.frame(N = n, num_nosig = g, scaled = g / g[1])
  fit <- list(vect_N_fit = 2:400, P_N_alpha = rep(0, 399),
              param_P_N_alpha = c(.1, .01, 0, .1, 0), Ns_fitting = 20)
  old <- list(empiric_PNalpha = empirical,
              simil_struct_Ns = list(smooth = list(vect_N_smooth = 2:101),
                                   fitting = fit, Ns_diff = 0))
  updated <- env$upgrade_ss(old, K)
  expect_equal(updated$scaled_bn$num_nosig, g)
  expect_equal(updated$scaled_bn$N, n)
  expect_equal(updated$simil_struct_Ns$smooth$vect_N_smooth, 2:101)
  expect_equal(updated$simil_struct_Ns$fitting$param_Sn, fit$param_P_N_alpha)
  expect_equal(updated$simil_struct_Ns$fitting$Ns_fitting, fit$Ns_fitting)
  expect_equal(sum(updated$simil_struct_Ns$smooth$Sn), 1)
  expect_error(env$upgrade_ss(old), "Supply the original")
  expect_error(env$upgrade_ss(updated, 301), "disagrees")
})
