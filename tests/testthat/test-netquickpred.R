# netquickpred(): correlation screening, the quickpred-style missingness
# screen, `steps` recursion, collinearity pruning (pairwise and VIF),
# dyad-level screening for network targets, and the netmice() integration
# (predictor_selection = "quickpred", `targets` dropping, `models`
# protection, precomputed objects).

# y <- x1 <- z chain with calibrated correlations: cor(x1, y) ~ 0.55,
# cor(z, x1) ~ 0.55, cor(z, y) ~ 0.30 - so with mincor = 0.45, z qualifies
# only through the recursion into (missing) x1, never directly.
qp_chain_attrs <- function(n = 600, seed = 11) {
  set.seed(seed)
  z <- rnorm(n)
  x1 <- z + rnorm(n, sd = 1.52)
  y <- x1 + rnorm(n, sd = 2.76)
  df <- data.frame(y = y, x1 = x1, z = z, noise = rnorm(n))
  df$y[sample(n, 120)] <- NA
  df$x1[sample(n, 90)] <- NA
  df
}

test_that("direct screening selects correlated predictors and drops noise", {
  df <- qp_chain_attrs()
  qp <- suppressMessages(netquickpred(df, targets = "y", mincor = 0.45,
                                      steps = 1))
  expect_s3_class(qp, "netquickpred")
  e <- qp$predictors$y
  expect_identical(e$kind, "attribute")
  expect_true("x1" %in% e$predictors)
  expect_false("noise" %in% e$predictors)
  expect_false("z" %in% e$predictors)
  # x1 has missing data and is selected -> it is imputed too, with its own
  # predictor entry
  expect_true("x1" %in% names(qp$predictors))
})

test_that("`steps` recursively adds predictors of missing predictors", {
  df <- qp_chain_attrs()
  qp1 <- suppressMessages(netquickpred(df, targets = "y", mincor = 0.45,
                                       steps = 1))
  qp3 <- suppressMessages(netquickpred(df, targets = "y", mincor = 0.45,
                                       steps = 3))
  expect_false("z" %in% qp1$predictors$y$predictors)
  expect_true("z" %in% qp3$predictors$y$predictors)
  expect_identical(unname(qp3$predictors$y$step[["z"]]), 2L)
  # z is fully observed, so the recursion stops there: nothing enters at
  # step 3 through z
  expect_false(any(qp3$predictors$y$step > 2L))
})

test_that("the missingness-indicator screen works like quickpred's", {
  set.seed(12)
  n <- 300
  w <- rnorm(n)
  df <- data.frame(y = rnorm(n), w = w)
  df$y[w > 0.5] <- NA  # w drives missingness, not the values
  qp_on <- suppressMessages(netquickpred(df, mincor = 0.3))
  qp_off <- suppressMessages(netquickpred(df, mincor = 0.3,
                                          use_missingness = FALSE))
  expect_true("w" %in% qp_on$predictors$y$predictors)
  expect_false("w" %in% qp_off$predictors$y$predictors)
})

test_that("pairwise collinearity pruning keeps one of two near-duplicates", {
  set.seed(13)
  n <- 300
  x1 <- rnorm(n)
  x1b <- x1 + rnorm(n, sd = 0.1)
  df <- data.frame(y = x1 + rnorm(n), x1 = x1, x1b = x1b)
  df$y[sample(n, 60)] <- NA
  qp <- suppressMessages(netquickpred(df, mincor = 0.3))
  expect_identical(sum(c("x1", "x1b") %in% qp$predictors$y$predictors), 1L)
  qp_none <- suppressMessages(netquickpred(df, mincor = 0.3,
                                           collin_method = "none"))
  expect_true(all(c("x1", "x1b") %in% qp_none$predictors$y$predictors))
})

test_that("VIF-based pruning also removes the near-duplicate", {
  set.seed(13)
  n <- 300
  x1 <- rnorm(n)
  x1b <- x1 + rnorm(n, sd = 0.1)
  df <- data.frame(y = x1 + rnorm(n), x1 = x1, x1b = x1b)
  df$y[sample(n, 60)] <- NA
  qp <- suppressMessages(netquickpred(df, mincor = 0.3,
                                      collin_method = "vif"))
  expect_identical(sum(c("x1", "x1b") %in% qp$predictors$y$predictors), 1L)
})

test_that("network measures are candidates for attribute targets, with
           netmice's naming, and ego-value-dependent features are withheld", {
  set.seed(14)
  n <- 60
  net <- matrix(rbinom(n * n, 1, 0.12), n, n); diag(net) <- 0
  y <- colSums(net) + rnorm(n)
  df <- data.frame(y = y, junk = rnorm(n))
  df$y[sample(n, 12)] <- NA
  qp <- suppressMessages(netquickpred(df, list(friends = net),
                                      targets = "y", mincor = 0.3))
  e <- qp$predictors$y
  expect_true("indegree" %in% e$net_features)  # bare name: single network
  expect_false(any(c("y_diff_out", "y_diff_in") %in% e$net_features))

  net2 <- matrix(rbinom(n * n, 1, 0.1), n, n); diag(net2) <- 0
  qp2 <- suppressMessages(netquickpred(df, list(friends = net, advice = net2),
                                       targets = "y", mincor = 0.3))
  expect_true("friends_indegree" %in% qp2$predictors$y$net_features)
})

test_that("network targets are screened at the dyad level with endogenous
           terms always kept", {
  set.seed(15)
  n <- 40
  dept <- sample(c("a", "b"), n, TRUE)
  p <- ifelse(outer(dept, dept, "=="), 0.35, 0.05)
  net <- matrix(rbinom(n * n, 1, p), n, n)
  net[sample(length(net), 120)] <- NA
  diag(net) <- 0
  df <- data.frame(dept = dept, junk = rnorm(n))
  qp <- suppressMessages(netquickpred(df, list(friends = net), mincor = 0.15))
  e <- qp$predictors$friends
  expect_identical(e$kind, "network")
  expect_true(all(c("reciprocity", "twopath") %in% e$dyad_terms))
  expect_true("dept_same" %in% e$dyad_terms)
  expect_false(any(grepl("^junk_", e$dyad_terms)))
})

test_that("cross-network cell terms are screened and pruned like other predictors", {
  set.seed(17)
  n <- 40
  y <- matrix(rbinom(n * n, 1, 0.25), n, n); diag(y) <- 0
  # z is a near-duplicate of y (a few flipped cells): collinear with y_tie
  z <- y
  flip <- sample(which(row(z) != col(z)), 25)
  z[flip] <- 1 - z[flip]
  # the target's ties follow y's, with noise - y_tie is the better predictor
  x <- y
  xflip <- sample(which(row(x) != col(x)), 200)
  x[xflip] <- 1 - x[xflip]
  x[sample(which(row(x) != col(x)), 120)] <- NA
  df <- data.frame(a = rnorm(n))
  qp <- suppressMessages(netquickpred(df, list(x = x, y = y, z = z),
                                      targets = "x", mincor = 0.15))
  e <- qp$predictors$x
  # (a) a cross-network term must pass the mincor screen to enter at all,
  # and (b) of two collinear cross-network carriers, only the one that
  # predicts the target better survives the pruning
  expect_true("y_tie" %in% e$dyad_terms)
  expect_false("z_tie" %in% e$dyad_terms)
  # the target's own endogenous terms remain exempt
  expect_true(all(c("reciprocity", "twopath") %in% e$dyad_terms))
})

test_that("cross-network cell terms below mincor are not selected", {
  set.seed(18)
  n <- 40
  x <- matrix(rbinom(n * n, 1, 0.2), n, n); diag(x) <- 0
  x[sample(which(row(x) != col(x)), 120)] <- NA
  noise <- matrix(rbinom(n * n, 1, 0.2), n, n); diag(noise) <- 0
  df <- data.frame(a = rnorm(n))
  qp <- suppressMessages(netquickpred(df, list(x = x, noise = noise),
                                      targets = "x", mincor = 0.3))
  expect_false("noise_tie" %in% qp$predictors$x$dyad_terms)
})

test_that("`targets` keeps selected predictors and drops the rest", {
  set.seed(16)
  n <- 300
  x1 <- rnorm(n)
  df <- data.frame(y = x1 + rnorm(n), x1 = x1, b = rnorm(n))
  df$y[sample(n, 50)] <- NA
  df$x1[sample(n, 40)] <- NA
  df$b[sample(n, 40)] <- NA
  qp <- suppressMessages(netquickpred(df, targets = "y", mincor = 0.3))
  expect_true("x1" %in% qp$keep$attributes)
  expect_identical(qp$drop$attributes, "b")
  expect_true("x1" %in% names(qp$predictors))
  expect_false("b" %in% names(qp$predictors))
  expect_output(print(qp), "netquickpred")
})

test_that("netmice(predictor_selection = 'quickpred') runs end to end and
           drops non-predictors named by `targets`", {
  set.seed(17)
  n <- 40
  net <- matrix(rbinom(n * n, 1, 0.12), n, n)
  net[sample(length(net), 60)] <- NA
  diag(net) <- 0
  x1 <- rnorm(n)
  df <- data.frame(y = x1 + rnorm(n, sd = 0.8), x1 = x1, b = rnorm(n))
  df$y[sample(n, 8)] <- NA
  df$b[sample(n, 6)] <- NA
  # mincor high enough that no feature derived from b clears the screen by
  # small-sample chance (cor sd ~ 1/sqrt(40) makes 0.3 too permissive here)
  fit <- suppressMessages(netmice(df, list(friends = net), m = 2, maxit = 2,
                                  predictor_selection = "quickpred",
                                  targets = c("y", "friends"), mincor = 0.5,
                                  seed = 42, printFlag = FALSE))
  expect_s3_class(fit, "netmids")
  expect_s3_class(fit$predictor_selection, "netquickpred")
  expect_false("b" %in% names(fit$data))
  expect_false("b" %in% names(fit$imp[[1]]))
  expect_true("x1" %in% names(fit$imp[[1]]))
  for (im in 1:2) {
    expect_false(anyNA(fit$imp[[im]]$y))
    expect_false(anyNA(fit$imp_nets[[im]]$friends))
  }
  expect_output(print(fit), "netquickpred")
})

test_that("variables referenced in `models` are protected from dropping", {
  set.seed(18)
  n <- 40
  net <- matrix(rbinom(n * n, 1, 0.12), n, n)
  net[sample(length(net), 40)] <- NA
  diag(net) <- 0
  x1 <- rnorm(n)
  df <- data.frame(y = x1 + rnorm(n, sd = 0.8), x1 = x1, b = rnorm(n))
  df$y[sample(n, 8)] <- NA
  df$b[sample(n, 6)] <- NA
  fit <- suppressMessages(netmice(df, list(friends = net), m = 1, maxit = 1,
                                  predictor_selection = "quickpred",
                                  targets = c("y", "friends"), mincor = 0.3,
                                  models = list("y ~ b"),
                                  seed = 1, printFlag = FALSE))
  expect_true("b" %in% names(fit$data))
  expect_false(anyNA(fit$imp[[1]]$b))
})

test_that("a precomputed netquickpred object is accepted by netmice()", {
  set.seed(19)
  n <- 40
  net <- matrix(rbinom(n * n, 1, 0.12), n, n)
  net[sample(length(net), 40)] <- NA
  diag(net) <- 0
  x1 <- rnorm(n)
  df <- data.frame(y = x1 + rnorm(n, sd = 0.8), x1 = x1)
  df$y[sample(n, 8)] <- NA
  qp <- suppressMessages(netquickpred(df, list(friends = net),
                                      mincor = 0.3))
  fit <- suppressMessages(netmice(df, list(friends = net), m = 1, maxit = 2,
                                  predictor_selection = qp,
                                  seed = 3, printFlag = FALSE))
  expect_s3_class(fit$predictor_selection, "netquickpred")
  expect_false(anyNA(fit$imp[[1]]$y))
  expect_false(anyNA(fit$imp_nets[[1]]$friends))
})

test_that("`targets` with predictor_selection = 'all' drops unneeded
           missing-data variables", {
  set.seed(20)
  n <- 40
  net <- matrix(rbinom(n * n, 1, 0.12), n, n); diag(net) <- 0
  df <- data.frame(y = rnorm(n), b = rnorm(n))
  df$y[sample(n, 8)] <- NA
  df$b[sample(n, 6)] <- NA
  fit <- suppressMessages(netmice(df, list(friends = net), m = 1, maxit = 1,
                                  targets = "y", seed = 5,
                                  printFlag = FALSE))
  expect_false("b" %in% names(fit$data))
  expect_false(anyNA(fit$imp[[1]]$y))
})

test_that("netquickpred: net_directed overrides inference and is validated", {
  set.seed(9)
  n <- 12
  m <- matrix(rbinom(n * n, 1, 0.2), n, n); diag(m) <- 0
  m[sample(which(row(m) != col(m)), 20)] <- NA
  attrs <- data.frame(age = rnorm(n), grp = sample(c("a", "b"), n, TRUE))
  attrs$age[c(2, 5)] <- NA

  expect_s3_class(netquickpred(attrs, list(m = m)), "netquickpred")
  expect_s3_class(netquickpred(attrs, list(m = m), net_directed = c(m = TRUE)),
                  "netquickpred")
  expect_error(netquickpred(attrs, list(m = m), net_directed = c(TRUE)),
               "named logical vector")
  expect_error(netquickpred(attrs, list(m = m), net_directed = c(m = NA)),
               "named logical vector")
  expect_error(netquickpred(attrs, list(m = m), net_directed = c(other = TRUE)),
               "named logical vector")
})
