test_that("pS uses actual support, includes endpoints and is symmetric", {
  A <- list(vect_N_smooth = 2:7, Sn = c(.3, .25, .2, .12, .08, .05))
  A$Ns_smooth <- sum(A$vect_N_smooth * A$Sn)
  B <- list(vect_N_smooth = 5:7, Sn = c(1, 0, 0), Ns_smooth = 5)
  expect_equal(prob_effect(smooth_A = A, smooth_B = B), .25)
  expect_equal(prob_effect(smooth_A = B, smooth_B = A), .25)

  B$Sn <- c(0, 0, 1)
  B$Ns_smooth <- 7
  expect_equal(prob_effect(smooth_A = A, smooth_B = B), .05)
  B$vect_N_smooth <- 8:10
  B$Sn <- c(1, 0, 0)
  B$Ns_smooth <- 8
  expect_equal(prob_effect(smooth_A = A, smooth_B = B), 0)
})

test_that("pS takes both cross-tails and uses unrounded thresholds", {
  A <- list(vect_N_smooth = 2:10, Sn = c(.8, rep(0, 7), .2))
  A$Ns_smooth <- sum(A$vect_N_smooth * A$Sn)
  B <- list(vect_N_smooth = 2:10, Sn = c(.1, .75, rep(0, 6), .15))
  B$Ns_smooth <- sum(B$vect_N_smooth * B$Sn)
  expected <- min(sum(A$Sn[A$vect_N_smooth >= B$Ns_smooth]),
                  sum(B$Sn[B$vect_N_smooth >= A$Ns_smooth]))
  expect_equal(prob_effect(smooth_A = A, smooth_B = B), expected)
  expect_equal(prob_effect(smooth_A = B, smooth_B = A), expected)
  expect_error(prob_effect(smooth_A = A), "both A and B")
  expect_error(prob_effect(), "Supply a pair")
})

test_that("parametric normalization and mean agree with direct sums", {
  n <- 5:10000
  raw <- c(.8, .01, .0001, .05, .000001)
  obj <- similaritystructure:::fitting_sn_norm(raw, n)
  direct <- (raw[1] + raw[2]*n + raw[3]*n^2 + raw[5]*n^3)*exp(-raw[4]*n)
  direct <- direct/sum(direct)
  expect_equal(obj$Sn, direct, tolerance = 1e-8)
  expect_equal(similaritystructure:::fitting_simil_size(obj$Params, 5)$Simil_Size,
               sum(n*direct), tolerance = 1e-7)
  expect_equal(similaritystructure:::ps_fitting(obj$Params, 24.3)$prob,
               sum(direct[n >= 24.3]), tolerance = 1e-8)
  A <- list(param_Sn = obj$Params, vect_N_fit = n, Ns_fitting = sum(n*direct))
  raw[4] <- .08
  obj_B <- similaritystructure:::fitting_sn_norm(raw, n)
  B <- list(param_Sn = obj_B$Params, vect_N_fit = n,
            Ns_fitting = sum(n*obj_B$Sn))
  expected <- min(sum(direct[n >= B$Ns_fitting]), sum(obj_B$Sn[n >= A$Ns_fitting]))
  expect_equal(prob_effect(fitting_A = A, fitting_B = B), expected, tolerance = 1e-8)
})
