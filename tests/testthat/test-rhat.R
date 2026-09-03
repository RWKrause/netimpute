# Rank-normalized split-R-hat. The golden values below were verified to match
# rstan::Rhat() to 10 decimal places on the post-burn-in draws, so they pin the
# implementation without making rstan a dependency. Note rstan's z_scale uses
# the midpoint plotting position (r - 1/2)/S; the posterior package uses the
# Blom position instead and gives slightly different numbers on short chains.

test_that(".rhat_split: reproduces rstan::Rhat on known inputs", {
  set.seed(99)
  converged <- matrix(rnorm(400), 100, 4)
  expect_equal(.rhat_split(converged), 1.0300376923, tolerance = 1e-9)

  drifting <- matrix(rnorm(200), 50, 4) + outer(seq_len(50) / 5, 1:4)
  expect_equal(.rhat_split(drifting), 3.6857180641, tolerance = 1e-9)

  separated <- matrix(rnorm(200), 50, 4) + rep(c(0, 3, 6, 9), each = 50)
  expect_equal(.rhat_split(separated), 2.6923760509, tolerance = 1e-9)
})

test_that(".rhat_split: converged chains sit near 1, bad chains well above 1.05", {
  set.seed(7)
  expect_lt(.rhat_split(matrix(rnorm(800), 200, 4)), 1.05)
  # a chain that is still climbing at the last iteration
  expect_gt(.rhat_split(matrix(rnorm(400), 100, 4) + seq_len(100) / 10), 1.05)
})

test_that(".rhat_split: degenerate traces give NA rather than erroring", {
  expect_true(is.na(.rhat_split(matrix(3, 40, 4))))        # constant: W = 0
  expect_true(is.na(.rhat_split(matrix(rnorm(6), 3, 2))))  # too few iterations
  # maxit = 4 -> 2 retained -> half-chains of length 1, no within-chain variance
  expect_true(is.na(.rhat_split(matrix(rnorm(8), 4, 2))))
  y <- matrix(rnorm(80), 40, 2)
  y[35, 1] <- NA                                            # in the RETAINED half
  expect_true(is.na(.rhat_split(y)))
  y2 <- matrix(rnorm(80), 40, 2)
  y2[35, 1] <- Inf
  expect_true(is.na(.rhat_split(y2)))
})

test_that(".rhat_split: burn-in discards the leading iterations", {
  set.seed(8)
  # a chain that is wild early and settled late must pass once the burn-in
  # removes the wild part
  x <- matrix(rnorm(400), 100, 4)
  x[1:50, ] <- x[1:50, ] + rep(c(0, 20, 40, 60), each = 50)
  expect_lt(.rhat_split(x), 1.05)
  expect_gt(.rhat_split(x, burnin = 0), 1.05)
})

test_that(".rhat_split: a single chain still works via splitting", {
  set.seed(9)
  expect_false(is.na(.rhat_split(matrix(rnorm(100), 100, 1))))
  expect_gt(.rhat_split(matrix(seq_len(100) + rnorm(100, 0, 0.01), 100, 1)), 1.05)
})

test_that(".warn_rhat: warns only on finite values above the threshold", {
  expect_warning(.warn_rhat(c(a = 1.2, b = 1.0), maxit = 20), "R-hat > 1.05")
  expect_warning(.warn_rhat(c(a = 1.2, b = 1.0), maxit = 20), "Increase `maxit`")
  expect_warning(.warn_rhat(c(a = 1.2), maxit = 20), "currently 20")
  expect_silent(.warn_rhat(c(a = 1.01, b = NA_real_), maxit = 20))
  expect_silent(.warn_rhat(c(a = NA_real_), maxit = 20))
  expect_silent(.warn_rhat(numeric(0), maxit = 20))
})

test_that("netmice: stores an R-hat per tracked quantity, NA for constant ones", {
  set.seed(2)
  n <- 30
  fr <- matrix(rbinom(n * n, 1, 0.15), n, n)
  diag(fr) <- 0
  off <- which(row(fr) != col(fr))
  fr[sample(off, 50)] <- NA
  attrs <- data.frame(age = rnorm(n, 35, 8), score = rnorm(n))
  attrs$age[sample(n, 5)] <- NA
  attrs$score[sample(n, 4)] <- NA
  fit <- suppressWarnings(
    netmice(attrs, list(friends = fr), m = 3, maxit = 30, seed = 1,
            printFlag = FALSE))
  expect_true(all(c("age (mean)", "age (variance)",
                    "friends (density)", "friends (imputed ties, mean)")
                  %in% names(fit$rhat)))
  # no isolates in this network -> the diagnostic is flat at 0 -> undefined
  expect_true(is.na(fit$rhat[["friends (n_isolates)"]]))
  expect_true(any(is.finite(fit$rhat)))
})
