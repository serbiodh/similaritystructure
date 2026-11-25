test_that("binomial smoothing matches direct mgcv predictions", {
  n <- unique(round(seq(2, 500, length.out = 60)))
  K <- 300L
  g <- round(K * exp(-(n - 2) / 15))
  empirical <- data.frame(N = n, num_nosig = g, scaled = g / g[1])
  sm <- similaritystructure:::smooth_binomial(2, 500, empirical, K)

  empirical$num_sig <- K - g
  model <- mgcv::gam(cbind(num_nosig, num_sig) ~ s(N, bs = "tp", k = 10),
                     family = stats::binomial(), data = empirical, method = "REML")
  beta <- stats::predict(model, newdata = data.frame(N = 2:500), type = "response")
  expect_equal(sm$Sn, beta / sum(beta), tolerance = 1e-8)
  expect_equal(sum(sm$Sn), 1)
  expect_true(all(sm$Sn >= 0))
  expect_equal(sm$Ns_smooth, sum(sm$vect_N_smooth * sm$Sn))
  expect_lt(sum(sm$Sn[sm$vect_N_smooth >= 100]), .01)
})

test_that("binomial smoothing allows an early maximum", {
  n <- 2:100
  x <- n - 2
  K <- 1000L
  g <- round(K * (.78 * exp(-x / 45) + .20 * (x / 12) * exp(1 - x / 12)))
  empirical <- data.frame(N = n, num_nosig = g, scaled = g / g[1])
  sm <- similaritystructure:::smooth_binomial(2, 100, empirical, K)
  expect_gt(max(sm$Sn), sm$Sn[1])
  expect_gt(sm$vect_N_smooth[which.max(sm$Sn)], 2)
  expect_lt(sm$vect_N_smooth[which.max(sm$Sn)], 25)
  expect_error(similaritystructure:::smooth_binomial(2, 100, empirical, 1),
               "Invalid empirical counts")
})
