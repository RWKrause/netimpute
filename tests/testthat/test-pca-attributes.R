# `PCA_attributes` budgets the attribute models separately from the tie
# models, and can switch their PCA collapse off entirely. The default
# (`NULL`) must inherit `PCA` and reproduce earlier behaviour exactly.

# Local copy of the joint attribute+network missingness fixture: two
# networks, so features carry the `<net>_` prefix.
pa_fixture <- function(n = 22, seed = 501) {
  set.seed(seed)
  attrs   <- fx_attrs(n = n, seed = seed)[c("age", "status", "dept")]
  friends <- fx_bin_directed(n = n, seed = seed + 1)
  advice  <- fx_weighted(n = n, seed = seed + 2)

  attrs$age[sample(n, 3)]    <- NA
  attrs$status[sample(n, 2)] <- NA
  attrs$dept[sample(n, 3)]   <- NA
  off <- which(row(friends) != col(friends))
  friends[sample(off, round(length(off) * 0.05))] <- NA

  list(attrs = attrs, nets = list(friends = friends, advice = advice))
}

test_that(".resolve_attribute_pca: NULL inherits, 'none' disables, list validates", {
  PCA <- .validate_pca(list(n = NULL, ratio = 10))

  # inherit -- the backward-compatible default
  expect_identical(.resolve_attribute_pca(NULL, PCA), PCA)

  # NULL return is the "no budget" sentinel .clean_predictor_matrix() reads
  expect_null(.resolve_attribute_pca("none", PCA))

  # an attribute-only budget, validated the same way as `PCA`
  expect_equal(.resolve_attribute_pca(list(n = 3), PCA),
               list(n = 3L, ratio = NULL))
  expect_equal(.resolve_attribute_pca(list(n = 4, ratio = 2), PCA),
               list(n = 4L, ratio = 2))

  # only "none" is accepted as a string
  expect_error(.resolve_attribute_pca("off", PCA), "must be NULL")
  expect_error(.resolve_attribute_pca("no", PCA), "must be NULL")
  # list validation is delegated, so `PCA`'s rules apply unchanged
  expect_error(.resolve_attribute_pca(list(k = 1), PCA), "only contain")
  expect_error(.resolve_attribute_pca(list(n = 0), PCA), "positive")
  expect_error(.resolve_attribute_pca(list(), PCA), "at least one")
})

test_that("a NULL budget disables the collapse in .clean_predictor_matrix()", {
  set.seed(11)
  x  <- matrix(rnorm(40 * 12), 40, 12,
               dimnames = list(NULL, paste0("v", 1:12)))
  ry <- rep(c(TRUE, FALSE), each = 20)

  # a real budget collapses to components
  capped <- .clean_predictor_matrix(x, max_cols = 3, ry = ry)
  expect_lte(ncol(capped), 3)
  expect_false(any(colnames(capped) %in% colnames(x)))

  # NULL imposes no cap: every column keeps its own name and coefficient
  free <- .clean_predictor_matrix(x, max_cols = NULL, ry = ry)
  expect_equal(ncol(free), 12)
  expect_equal(colnames(free), colnames(x))
})

test_that("PCA_attributes = NULL reproduces the inherited-PCA imputations exactly", {
  fx <- pa_fixture()

  a <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 42,
               PCA = list(n = 2), printFlag = FALSE)
  b <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 42,
               PCA = list(n = 2), PCA_attributes = NULL, printFlag = FALSE)

  expect_equal(complete_netmice(a, 1)$data, complete_netmice(b, 1)$data)
  expect_equal(complete_netmice(a, 1)$networks, complete_netmice(b, 1)$networks)
})

test_that("PCA_attributes = 'none' keeps attribute predictors uncollapsed", {
  fx <- pa_fixture()

  # a budget of 1 forces a hard collapse; "none" must undo it, so the two
  # runs cannot agree on the attribute imputations
  tight <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 42,
                   PCA = list(n = 1), printFlag = FALSE)
  free  <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 42,
                   PCA = list(n = 1), PCA_attributes = "none",
                   printFlag = FALSE)

  expect_false(isTRUE(all.equal(complete_netmice(tight, 1)$data$age,
                                complete_netmice(free, 1)$data$age)))
  # and it still produces a usable completion
  expect_false(anyNA(complete_netmice(free, 1)$data$age))
  expect_false(anyNA(complete_netmice(free, 1)$data$dept))
})

test_that("PCA_attributes = 'none' still imputes the networks", {
  fx <- pa_fixture()

  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 42,
                 PCA = list(n = 2), PCA_attributes = "none",
                 printFlag = FALSE)
  net <- complete_netmice(fit, 1)$networks$friends
  off <- which(row(net) != col(net))
  expect_false(anyNA(net[off]))
})

test_that("a hand-built netquickpred selection survives a NULL attribute budget", {
  # this is the netmice.R:952 path: the "selection exceeds budget" notice must
  # not be evaluated when there is no budget for the selection to exceed
  fx <- pa_fixture()

  sel <- structure(list(
    predictors = list(age = list(predictors = "status",
                                 net_features = c("friends_indegree",
                                                  "friends_outdegree"),
                                 step = 1L, relevance = NA_real_)),
    targets = c("age", "status", "dept", "friends"), keep = character(0),
    drop = list(attributes = character(0), networks = character(0)),
    mincor = NA_real_, steps = NA_integer_,
    collin_method = "none", collin_threshold = Inf,
    use_missingness = FALSE, measure_set = "core"
  ), class = "netquickpred")

  expect_no_error(
    fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 42,
                   PCA = list(n = 1), PCA_attributes = "none",
                   predictor_selection = sel, printFlag = FALSE)
  )
  expect_false(anyNA(complete_netmice(fit, 1)$data$age))
  # targets without an entry fall back to the full auto predictor set
  expect_false(anyNA(complete_netmice(fit, 1)$data$dept))
})
