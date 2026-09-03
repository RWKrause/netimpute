# `PCA = list(n, ratio)` is the single dimensionality safeguard, replacing the
# former `n_components` argument and the internal max(5, n/3) cap.

test_that(".validate_pca: accepts n, ratio, or both; rejects the rest", {
  expect_equal(.validate_pca(list(n = 5)), list(n = 5L, ratio = NULL))
  expect_equal(.validate_pca(list(ratio = 20)), list(n = NULL, ratio = 20))
  expect_equal(.validate_pca(list(n = 5, ratio = 20)), list(n = 5L, ratio = 20))
  # n is floored to an integer
  expect_equal(.validate_pca(list(n = 5.9))$n, 5L)

  expect_error(.validate_pca(list()), "at least one")
  expect_error(.validate_pca(list(n = NULL, ratio = NULL)), "at least one")
  expect_error(.validate_pca(5), "must be a list")
  expect_error(.validate_pca(list(k = 1)), "only contain")
  expect_error(.validate_pca(list(n = -1)), "positive")
  expect_error(.validate_pca(list(n = 0)), "positive")
  expect_error(.validate_pca(list(ratio = "a")), "positive")
  expect_error(.validate_pca(list(n = c(1, 2))), "positive")
  expect_error(.validate_pca(list(n = NA_real_)), "positive")
})

test_that(".pca_denom: binary targets count events, others count rows", {
  # the whole point: 870 rows but only 130 ties
  y <- c(rep(1, 130), rep(0, 740))
  expect_equal(.pca_denom(y, NULL, binary = TRUE), 130L)
  expect_equal(.pca_denom(y, NULL, binary = FALSE), 870L)
  # the rarer class is what counts, whichever it is
  expect_equal(.pca_denom(c(rep(1, 700), rep(0, 20)), NULL, TRUE), 20L)
  # continuous
  expect_equal(.pca_denom(rnorm(30), NULL, FALSE), 30L)
  # only one class observed -> nothing to discriminate
  expect_equal(.pca_denom(rep(1, 50), NULL, TRUE), 0L)
  # NA and unobserved entries never count
  expect_equal(.pca_denom(c(1, 0, 1, NA), c(TRUE, TRUE, FALSE, TRUE), TRUE), 1L)
  expect_equal(.pca_denom(c(NA, NA), NULL, FALSE), 0L)
  # works on character/factor binaries too
  expect_equal(.pca_denom(c("a", "a", "b"), NULL, TRUE), 1L)
})

test_that(".pca_budget: takes the smaller cap and never returns 0", {
  p <- function(...) .validate_pca(list(...))
  expect_equal(.pca_budget(p(n = 5), 1000), 5L)
  expect_equal(.pca_budget(p(ratio = 10), 130), 13L)
  expect_equal(.pca_budget(p(n = 5, ratio = 10), 130), 5L)    # n binds
  expect_equal(.pca_budget(p(n = 50, ratio = 10), 130), 13L)  # ratio binds
  # the floor: a tiny sample still gets one predictor, not zero
  expect_equal(.pca_budget(p(ratio = 10), 3), 1L)
  expect_equal(.pca_budget(p(ratio = 10), 0), 1L)
})

test_that("netmice: PCA is validated and recorded on the object", {
  set.seed(3)
  n <- 25
  fr <- matrix(rbinom(n * n, 1, 0.2), n, n)
  diag(fr) <- 0
  off <- which(row(fr) != col(fr))
  fr[sample(off, 25)] <- NA
  attrs <- data.frame(age = rnorm(n, 35, 8))
  attrs$age[sample(n, 4)] <- NA
  fit <- suppressWarnings(
    netmice(attrs, list(friends = fr), m = 1, maxit = 2, seed = 1,
            printFlag = FALSE, PCA = list(n = 4, ratio = 20)))
  expect_equal(fit$PCA, list(n = 4L, ratio = 20))
  expect_error(
    netmice(attrs, list(friends = fr), m = 1, maxit = 2, PCA = list(bogus = 1)),
    "only contain")
})

test_that("PCA$n caps the components in dyad_regression and net_predictors", {
  # fx_nets() also carries a signed network, which dyad_regression rejects
  all_nets <- fx_nets(n = 20)
  nets <- all_nets[c("friends_bin", "colleagues_bin", "advice_weighted")]
  attrs <- fx_attrs(n = 20)[c("age", "status")]
  res <- dyad_regression(nets, attrs, target = "friends_bin",
                         other_net_predictors = "pca", PCA = list(n = 2),
                         fit = FALSE)
  expect_lte(sum(grepl("^other_net_PC", names(res$data))), 2)

  pp <- net_predictors(list(nets$friends_bin), list(attrs), output = "pca",
                       PCA = list(n = 3))
  expect_lte(sum(grepl("^PC[0-9]+$", names(pp$components))), 3)
})

# The budget applies on top of `predictor_selection`, and at ordinary network
# sizes it usually still binds - so a user who selected predictors to get named
# coefficients may silently get components instead. That is worth saying.

make_budget_fixture <- function(n = 60) {
  set.seed(9)
  mk <- function(p, s) {
    set.seed(s)
    m <- matrix(rbinom(n * n, 1, p), n, n)
    diag(m) <- 0
    off <- which(row(m) != col(m))
    m[sample(off, 60)] <- NA
    m
  }
  attrs <- data.frame(age = rnorm(n, 35, 8), score = rnorm(n),
                      tenure = rnorm(n), perf = rnorm(n),
                      grp = sample(c("a", "b"), n, TRUE))
  attrs$age[sample(n, 8)] <- NA
  attrs$score[sample(n, 7)] <- NA
  list(attrs = attrs, nets = list(friends = mk(.15, 1), advice = mk(.12, 2)))
}

test_that(".pca_selection_notice: reports each affected target, or nothing", {
  expect_message(
    .pca_selection_notice(list(list(target = "age", n_in = 58, budget = 5))),
    "58 selected vs budget 5")
  expect_message(
    .pca_selection_notice(list(list(target = "age", n_in = 58, budget = 5))),
    "Raise `mincor`")
  expect_silent(.pca_selection_notice(list()))
})

test_that("netmice: notes when a selection is still larger than the PCA budget", {
  fx <- make_budget_fixture()
  expect_message(
    suppressWarnings(netmice(fx$attrs, fx$nets, predictor_selection = "quickpred",
                             m = 1, maxit = 1, seed = 1, printFlag = TRUE)),
    "still leaves more predictors than the `PCA` budget")
})

test_that("netmice: no such note when the selection fits inside the budget", {
  fx <- make_budget_fixture()
  # mincor 0.5 leaves almost nothing, so the budget never binds
  qp <- netquickpred(fx$attrs, fx$nets, mincor = 0.5)
  msgs <- capture_messages(
    suppressWarnings(netmice(fx$attrs, fx$nets, predictor_selection = qp,
                             m = 1, maxit = 1, seed = 1, printFlag = TRUE)))
  expect_false(any(grepl("still leaves more predictors", msgs)))
})

test_that("netmice: the note is not raised for predictor_selection = 'all'", {
  # "all" never promised named coefficients, so components are expected there
  fx <- make_budget_fixture()
  msgs <- capture_messages(
    suppressWarnings(netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 1,
                             printFlag = TRUE)))
  expect_false(any(grepl("still leaves more predictors", msgs)))
})

test_that("netmice: the note is raised once, not once per visit", {
  fx <- make_budget_fixture()
  msgs <- capture_messages(
    suppressWarnings(netmice(fx$attrs, fx$nets, predictor_selection = "quickpred",
                             m = 2, maxit = 3, seed = 1, printFlag = TRUE)))
  expect_equal(sum(grepl("still leaves more predictors", msgs)), 1L)
})

test_that("netmice: printFlag = FALSE silences the note", {
  fx <- make_budget_fixture()
  msgs <- capture_messages(
    suppressWarnings(netmice(fx$attrs, fx$nets, predictor_selection = "quickpred",
                             m = 1, maxit = 1, seed = 1, printFlag = FALSE)))
  expect_false(any(grepl("still leaves more predictors", msgs)))
})
