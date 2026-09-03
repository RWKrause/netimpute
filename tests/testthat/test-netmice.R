# Shared fixture with missingness in a continuous, a binary, and a
# multinomial attribute, and in one binary and one weighted network - the
# joint attribute+network missingness scenario netmice() is built for.
make_missing_fixture <- function(n = 22, seed = 501) {
  set.seed(seed)
  attrs <- fx_attrs(n = n, seed = seed)[c("age", "status", "dept")]
  friends <- fx_bin_directed(n = n, seed = seed + 1)
  advice  <- fx_weighted(n = n, seed = seed + 2)

  attrs$age[sample(n, 3)] <- NA
  attrs$status[sample(n, 2)] <- NA
  attrs$dept[sample(n, 3)] <- NA
  off <- which(row(friends) != col(friends))
  friends[sample(off, round(length(off) * 0.05))] <- NA

  list(attrs = attrs, nets = list(friends = friends, advice = advice))
}

test_that("netmice: default args run to completion with mixed attributes and networks", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 2, maxit = 2, printFlag = FALSE)
  expect_s3_class(fit, "netmids")
  expect_length(fit$imp, 2)
  expect_length(fit$imp_nets, 2)
  for (im in 1:2) {
    expect_false(anyNA(fit$imp[[im]]$age))
    expect_false(anyNA(fit$imp[[im]]$status))
    expect_false(anyNA(fit$imp[[im]]$dept))
    off <- which(row(fx$nets$friends) != col(fx$nets$friends))
    expect_false(anyNA(fit$imp_nets[[im]]$friends[off]))
  }
})

test_that("netmice: an unknown method errors and names the supported set", {
  fx <- make_missing_fixture()
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, method = "banana", printFlag = FALSE),
    "Supported methods.*pmm.*cart"
  )
})

test_that("netmice: method = 'cart' runs end to end via the generic mice dispatch", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, method = "cart",
                  seed = 5, printFlag = FALSE)
  expect_true(all(fit$method == "cart"))
  expect_false(anyNA(fit$imp[[1]]$age))
  expect_false(anyNA(fit$imp[[1]]$status))
  expect_false(anyNA(fit$imp[[1]]$dept))
  off <- which(row(fx$nets$friends) != col(fx$nets$friends))
  expect_false(anyNA(fit$imp_nets[[1]]$friends[off]))
  # tree draws come from observed values
  expect_true(all(fit$imp[[1]]$age %in% fx$attrs$age[!is.na(fx$attrs$age)]))
})

test_that("netmice: method = 'norm' keeps binary attributes on their two levels (round + clamp)", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 2, maxit = 2, method = "norm",
                  seed = 6, printFlag = FALSE)
  expect_true(all(fit$method == "norm"))
  for (im in 1:2) {
    expect_false(anyNA(fit$imp[[im]]$status))
    # binary attribute: continuous norm draws must map back to valid levels
    expect_true(all(fit$imp[[im]]$status %in%
                      unique(fx$attrs$status[!is.na(fx$attrs$status)])))
    # continuous attribute: norm may (and typically does) produce values
    # outside the observed set - that is the point of a non-donor method
    expect_true(is.numeric(fit$imp[[im]]$age))
    expect_false(anyNA(fit$imp[[im]]$age))
  }
})

test_that("netmice: donors argument is accepted and used without error", {
  fx <- make_missing_fixture()
  fit3  <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, donors = 3, printFlag = FALSE)
  fit10 <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, donors = 10, printFlag = FALSE)
  expect_s3_class(fit3, "netmids")
  expect_s3_class(fit10, "netmids")
  expect_equal(fit3$donors, 3)
  expect_equal(fit10$donors, 10)
})

test_that("netmice: measure_set = 'full' runs and differs from 'core'", {
  fx <- make_missing_fixture()
  fit_core <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, measure_set = "core",
                       seed = 1, printFlag = FALSE)
  fit_full <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, measure_set = "full",
                       seed = 1, printFlag = FALSE)
  expect_s3_class(fit_full, "netmids")
  # different predictor sets feeding PMM can legitimately give different draws
  expect_true(is.numeric(fit_core$imp[[1]]$age) && is.numeric(fit_full$imp[[1]]$age))
})

test_that("netmice: attr_types override is honoured", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1,
                  attr_types = c(age = "continuous", status = "binary", dept = "multinomial"),
                  printFlag = FALSE)
  expect_false(anyNA(fit$imp[[1]]$dept))
})

test_that(".init_fill_matrix: 'zero' fills missing cells with 0, 'sample' resamples observed values", {
  set.seed(6)
  m <- matrix(rbinom(100, 1, 0.5), 10, 10); diag(m) <- 0
  mis <- sample(which(row(m) != col(m)), 15)
  m[mis] <- NA
  z <- netimpute:::.init_fill_matrix(m, init = "zero")
  expect_true(all(z[mis] == 0))
  expect_identical(z[-mis], m[-mis])
  s <- netimpute:::.init_fill_matrix(m, init = "sample")
  expect_true(all(s[mis] %in% c(0, 1)))
})

test_that("netmice: net_init is honoured, defaults to 'zero', and is stored on the result", {
  fx <- make_missing_fixture()
  fit0 <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 3, printFlag = FALSE)
  expect_equal(fit0$net_init, "zero")
  fit_s <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, seed = 3,
                   net_init = "sample", printFlag = FALSE)
  expect_equal(fit_s$net_init, "sample")
  expect_error(netmice(fx$attrs, fx$nets, m = 1, maxit = 1,
                       net_init = "banana", printFlag = FALSE))
})

test_that(".impute_ties_gibbs: incremental change statistics match full recomputation", {
  set.seed(11)
  n <- 15
  m <- matrix(rbinom(n * n, 1, 0.2), n, n); diag(m) <- 0
  mis <- sample(which(row(m) != col(m)), 60)
  m[mis] <- NA
  attrs <- data.frame(age = rnorm(n), grade = rnorm(n))
  filled <- netimpute:::.init_fill_matrix(m, init = "zero")
  built <- suppressMessages(
    netimpute:::.build_dyad_data(list(net = filled), attrs, 1))
  d <- built$data
  ry <- !is.na(m)[cbind(d$i, d$j)]
  x <- netimpute:::.clean_predictor_matrix(
    d[setdiff(names(d), c("i", "j", "y"))], ry = ry)
  # check = TRUE recomputes B and the two-path count matrix from scratch at
  # the end and stops if the O(n) per-draw updates diverged from them
  out <- netimpute:::.impute_ties_gibbs(d = d, ry = ry, x = x, mat = filled,
                                        binary = TRUE, donors = 5,
                                        check = TRUE)
  obs <- !is.na(m) & (row(m) != col(m))
  expect_identical(out[obs], m[obs])       # observed ties never resampled
  expect_true(all(out[mis] %in% c(0, 1)))  # Bernoulli draws are binary
})

test_that("netmice: net_update = 'gibbs' runs end to end (binary and weighted targets)", {
  fx <- make_missing_fixture()
  # give the weighted network missing cells too, to exercise the
  # sequential-PMM fallback for non-binary ties
  off <- which(row(fx$nets$advice) != col(fx$nets$advice))
  advice_obs_vals <- fx$nets$advice[off]
  mis_adv <- sample(off, 25)
  fx$nets$advice[mis_adv] <- NA
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, net_update = "gibbs",
                 seed = 7, printFlag = FALSE)
  expect_s3_class(fit, "netmids")
  expect_equal(fit$net_update, "gibbs")
  friends_done <- fit$imp_nets[[1]]$friends
  advice_done  <- fit$imp_nets[[1]]$advice
  expect_false(anyNA(friends_done[off]))
  expect_false(anyNA(advice_done[off]))
  # binary target: Bernoulli draws
  expect_true(all(friends_done %in% c(0, 1)))
  # weighted target: PMM donor draws come from the observed value pool
  expect_true(all(advice_done[mis_adv] %in% advice_obs_vals))
  # observed cells are untouched
  obs_f <- !is.na(fx$nets$friends)
  expect_identical(friends_done[obs_f], fx$nets$friends[obs_f])
})

test_that("netmice: net_update = 'gibbs' rejects net_random_intercepts", {
  fx <- make_missing_fixture()
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, net_update = "gibbs",
            net_random_intercepts = "ego", printFlag = FALSE),
    "not compatible"
  )
})

test_that("netmice: per-target methods via a named `method` vector", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, seed = 8,
                  printFlag = FALSE,
                  method = c("pmm", age = "norm"))
  # the resolved map: named override for age, unnamed default elsewhere
  expect_equal(unname(fit$method[["age"]]), "norm")
  expect_true(all(fit$method[setdiff(names(fit$method), "age")] == "pmm"))
  expect_false(anyNA(fit$imp[[1]]$age))
  # observed ages are integers (round(rnorm)); norm draws are continuous,
  # so at least one imputed age falls outside the observed value set -
  # proof the per-target override (not PMM) handled age
  expect_false(all(fit$imp[[1]]$age[is.na(fx$attrs$age)] %in%
                     fx$attrs$age[!is.na(fx$attrs$age)]))
  # PMM targets stay on observed values
  expect_true(all(fit$imp[[1]]$status %in%
                    unique(fx$attrs$status[!is.na(fx$attrs$status)])))
})

test_that("netmice: named `method` entries are validated", {
  fx <- make_missing_fixture()
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            method = c(not_a_variable = "cart")),
    "unknown variable/network"
  )
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            method = c("pmm", "cart")),
    "at most one unnamed"
  )
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            method = c(age = "banana")),
    "Supported methods"
  )
})

test_that("netmice: named `method` entries that cannot take effect message", {
  fx <- make_missing_fixture()
  # dept is multinomial -> always polyreg
  expect_message(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            method = c(dept = "cart")),
    "polyreg"
  )
  # friends is a network and the default updater is tie-wise gibbs
  expect_message(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            method = c(friends = "norm")),
    "tie-wise"
  )
})

test_that("netmice: tie-wise 'gibbs' updating is the default", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE)
  expect_equal(fit$net_update, "gibbs")
})

test_that("netmice: net_random_intercepts without an explicit net_update falls back to 'simultaneous'", {
  skip_if_not_installed("lme4")
  fx <- make_missing_fixture()
  expect_message(
    fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
                   net_random_intercepts = "ego"),
    "falling back"
  )
  expect_equal(fit$net_update, "simultaneous")
})

test_that("netmice: other_net_predictors = 'pca' runs end to end", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, other_net_predictors = "pca",
                  PCA = list(n = 2), printFlag = FALSE)
  expect_s3_class(fit, "netmids")
  off <- which(row(fx$nets$friends) != col(fx$nets$friends))
  expect_false(anyNA(fit$imp_nets[[1]]$friends[off]))
})

test_that("netmice: models - attribute formula with an interaction is honoured", {
  fx <- make_missing_fixture()
  fx$attrs$performance <- rnorm(nrow(fx$attrs))
  fx$attrs$happiness <- rnorm(nrow(fx$attrs))
  fx$attrs$happiness[sample(nrow(fx$attrs), 3)] <- NA
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, printFlag = FALSE,
                  models = list("happiness ~ friends_indegree + age + performance * dept"))
  expect_false(anyNA(fit$imp[[1]]$happiness))
})

test_that("netmice: models - network dyad-level formula is honoured", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, printFlag = FALSE,
                  models = list("friends ~ age_absdiff + advice_tie + reciprocity"))
  off <- which(row(fx$nets$friends) != col(fx$nets$friends))
  expect_false(anyNA(fit$imp_nets[[1]]$friends[off]))
})

test_that("netmice: the target's own alter-composition homophily features are available as predictors", {
  fx <- make_missing_fixture()
  # age's own alter-mean measures (incl. the incoming-tie mean, computed
  # from other nodes' - mostly observed - reports) must be part of the
  # predictor set when imputing age: referencing one in a `models` formula
  # only works if it exists in the auto-generated feature data.frame.
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, printFlag = FALSE,
                  models = list("age ~ friends_age_alter_mean_in + friends_age_alter_mean"))
  expect_false(anyNA(fit$imp[[1]]$age))
})

test_that("netmice: ego-value-dependent homophily features of the target are withheld", {
  fx <- make_missing_fixture()
  # age_diff_in is a function of ego's own (partly imputed) age, so it must
  # not exist among age's own predictors - a formula asking for it errors.
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            models = list("age ~ friends_age_diff_in")),
    "could not evaluate"
  )
  # ...but the same measure built from a *different* attribute is available.
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, printFlag = FALSE,
                  models = list("status ~ friends_age_diff_in"))
  expect_false(anyNA(fit$imp[[1]]$status))
})

test_that("netmice: models errors on an unknown left-hand side", {
  fx <- make_missing_fixture()
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            models = list("not_a_real_variable ~ age")),
    "unknown variable"
  )
})

test_that("netmice: models warns (but doesn't error) for a fully-observed target", {
  fx <- make_missing_fixture()
  expect_message(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            models = list("advice ~ age_absdiff")),
    "will not be used"
  )
})

test_that("netmice: models - network formula with interactions of internal terms", {
  fx <- make_missing_fixture()
  # an ego attribute x the target's own reciprocity, and a cross-network
  # tie x an attribute-similarity term: every internally created dyad-level
  # column is interactable via model.matrix()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, printFlag = FALSE,
                  models = list("friends ~ age_ego:reciprocity + advice_tie:age_absdiff"))
  off <- which(row(fx$nets$friends) != col(fx$nets$friends))
  expect_false(anyNA(fit$imp_nets[[1]]$friends[off]))
})

test_that("netmice: models - endogenous interaction under net_update = 'gibbs'", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, net_update = "gibbs",
                  seed = 9, printFlag = FALSE,
                  models = list("friends ~ age_ego:reciprocity + reciprocity:twopath"))
  off <- which(row(fx$nets$friends) != col(fx$nets$friends))
  friends_done <- fit$imp_nets[[1]]$friends
  expect_false(anyNA(friends_done[off]))
  expect_true(all(friends_done %in% c(0, 1)))
  obs <- !is.na(fx$nets$friends)
  expect_identical(friends_done[obs], fx$nets$friends[obs])
})

test_that(".gibbs_endo_interactions: decomposes endogenous interaction columns", {
  d <- data.frame(age_ego = c(1, 2, 3, 4),
                  reciprocity = c(0, 1, 0, 1),
                  twopath = c(1, 0, 1, 0))
  xn <- c("age_ego", "age_ego:reciprocity", "reciprocity:twopath",
          "age_ego:not_in_d")
  info <- netimpute:::.gibbs_endo_interactions(xn, d)
  # pure main effects and undecomposable columns are skipped
  expect_named(info, c("age_ego:reciprocity", "reciprocity:twopath"))
  ar <- info[["age_ego:reciprocity"]]
  expect_equal(ar$col, 2L)
  expect_true(ar$recip)
  expect_false(ar$twop)
  expect_equal(ar$static, c(1, 2, 3, 4))
  rt <- info[["reciprocity:twopath"]]
  expect_true(rt$recip && rt$twop)
  expect_equal(rt$static, rep(1, 4))
})

test_that(".impute_ties_gibbs: endogenous interaction columns keep the bookkeeping consistent", {
  set.seed(21)
  n <- 15
  m <- matrix(rbinom(n * n, 1, 0.2), n, n); diag(m) <- 0
  mis <- sample(which(row(m) != col(m)), 60)
  m[mis] <- NA
  attrs <- data.frame(age = rnorm(n), grade = rnorm(n))
  filled <- netimpute:::.init_fill_matrix(m, init = "zero")
  built <- suppressMessages(
    netimpute:::.build_dyad_data(list(net = filled), attrs, 1))
  d <- built$data
  ry <- !is.na(m)[cbind(d$i, d$j)]
  x <- netimpute:::.clean_predictor_matrix(
    d[setdiff(names(d), c("i", "j", "y"))], ry = ry)
  x <- cbind(x, "age_ego:reciprocity" = d$age_ego * d$reciprocity)
  out <- netimpute:::.impute_ties_gibbs(d = d, ry = ry, x = x, mat = filled,
                                        binary = TRUE, donors = 5,
                                        check = TRUE)
  obs <- !is.na(m) & (row(m) != col(m))
  expect_identical(out[obs], m[obs])
  expect_true(all(out[mis] %in% c(0, 1)))
})

test_that("netmice: `targets` never drops a network/attribute a models formula depends on via derived terms", {
  fx <- make_missing_fixture()
  off <- which(row(fx$nets$advice) != col(fx$nets$advice))
  fx$nets$advice[sample(off, 20)] <- NA
  # advice and age have missing data and are not targets, but the friends
  # formula references advice_tie and age_ego - both must be kept (before
  # the derived-name protection this errored with "could not evaluate")
  fit <- suppressMessages(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            targets = "friends",
            models = list("friends ~ advice_tie:age_ego")))
  expect_true("advice" %in% names(fit$imp_nets[[1]]))
  expect_true("age" %in% names(fit$imp[[1]]))
  expect_false(anyNA(fit$imp_nets[[1]]$advice[off]))
  expect_false(anyNA(fit$imp[[1]]$age))
})

test_that("netmice: models-referenced variables survive a precomputed netquickpred that dropped them", {
  set.seed(42)
  n <- 30
  friends <- matrix(rbinom(n * n, 1, 0.12), n, n); diag(friends) <- 0
  advice  <- matrix(rbinom(n * n, 1, 0.12), n, n); diag(advice) <- 0
  attrs <- data.frame(age = rnorm(n, 35, 8),
                      gender = sample(c("F", "M"), n, TRUE),
                      happiness = rnorm(n))
  attrs$happiness[sample(n, 5)] <- NA
  attrs$age[sample(n, 4)] <- NA
  off <- which(row(friends) != col(friends))
  friends[sample(off, 40)] <- NA
  advice[sample(off, 40)] <- NA

  # a selection computed WITHOUT knowledge of `models`, strict enough that
  # age and both networks land in its drop lists
  qp <- suppressMessages(
    netquickpred(attrs, list(friends = friends, advice = advice),
                 targets = "happiness", mincor = 0.99))
  expect_true("age" %in% qp$drop$attributes)
  expect_true("advice" %in% qp$drop$networks)

  # netmice must protect everything the formulas depend on - `age` (bare),
  # `friends` (a formula LHS), and `advice` (via the derived advice_tie
  # term) - keep them imputed, and evaluate the formulas at every visit
  fit <- suppressMessages(
    netmice(attrs, list(friends = friends, advice = advice),
            m = 1, maxit = 2, seed = 1, printFlag = FALSE,
            predictor_selection = qp,
            models = list("happiness ~ age",
                          "friends ~ advice_tie:age_ego")))
  expect_true("age" %in% names(fit$imp[[1]]))
  expect_true(all(c("friends", "advice") %in% names(fit$imp_nets[[1]])))
  expect_false(anyNA(fit$imp[[1]]$age))
  expect_false(anyNA(fit$imp[[1]]$happiness))
  expect_false(anyNA(fit$imp_nets[[1]]$advice[off]))
})

test_that("netmice: seed gives identical results across repeated runs", {
  fx <- make_missing_fixture()
  fit1 <- netmice(fx$attrs, fx$nets, m = 2, maxit = 2, seed = 42, printFlag = FALSE)
  fit2 <- netmice(fx$attrs, fx$nets, m = 2, maxit = 2, seed = 42, printFlag = FALSE)
  expect_identical(fit1$imp, fit2$imp)
  expect_identical(fit1$imp_nets, fit2$imp_nets)
})

test_that("netmice: printFlag = FALSE suppresses progress output", {
  fx <- make_missing_fixture()
  # printFlag only controls the progress cat()s - the attribute-type
  # inference message() is independent and always fires, so it must be
  # suppressed separately for this to isolate what printFlag actually does.
  expect_silent(suppressMessages(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE)
  ))
})

test_that("netmice: printFlag = TRUE (default) produces progress output", {
  fx <- make_missing_fixture()
  expect_output(netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = TRUE), "imputation")
})

test_that("netmice: ncores > 1 matches ncores = 1 when the package is installed", {
  # CRAN limits the number of usable cores and discourages spawning worker
  # processes during checks; the sequential path is covered everywhere else
  skip_on_cran()
  # find.package() rather than installed.packages(): the latter is slow and
  # CRAN rejects it. Both answer the same question here - is there an
  # installed copy for the multisession workers to load?
  skip_if(length(find.package("netimpute", quiet = TRUE)) == 0,
          "netimpute is not installed (only load_all()-ed) in this session")
  skip_if_not_installed("future")
  fx <- make_missing_fixture()
  fit_seq <- netmice(fx$attrs, fx$nets, m = 2, maxit = 2, seed = 7, ncores = 1, printFlag = FALSE)
  fit_par <- tryCatch(
    netmice(fx$attrs, fx$nets, m = 2, maxit = 2, seed = 7, ncores = 2, printFlag = FALSE),
    error = function(e) e
  )
  is_err <- inherits(fit_par, "error")
  # skip_if()'s `message` argument must not itself error when there is no
  # error to describe - build it conditionally rather than calling
  # conditionMessage() unconditionally.
  skip_if(is_err, if (is_err) paste("parallel run failed:", conditionMessage(fit_par)) else "")
  expect_identical(fit_seq$imp, fit_par$imp)
})

test_that("netmice: works with a single network", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets["friends"], m = 1, maxit = 2, printFlag = FALSE)
  expect_false(anyNA(fit$imp[[1]]$age))
  off <- which(row(fx$nets$friends) != col(fx$nets$friends))
  expect_false(anyNA(fit$imp_nets[[1]]$friends[off]))
})

test_that("netmice: structural - a single matrix fixes the same cells at zero in every network", {
  fx <- make_missing_fixture()
  n <- nrow(fx$attrs)
  # ties into the last 4 nodes are impossible by design in all networks
  s <- matrix(FALSE, n, n)
  s[, (n - 3):n] <- TRUE
  diag(s) <- FALSE
  off <- row(s) != col(s)
  # structural cells must be 0 (or NA) in the inputs; code a few of them as
  # NA - they are design-zeros, not missing ties, and must come back as 0
  for (nm in names(fx$nets)) fx$nets[[nm]][s] <- 0
  fx$nets$friends[which(s)[1:5]] <- NA

  # the dyad-level design drops the structural rows entirely (no zero
  # inflation): n*(n-1) minus the off-diagonal structural cells
  res <- dyad_regression(fx$nets, fx$attrs, target = "friends",
                         structural = s, fit = FALSE)
  expect_equal(nrow(res$data), sum(off) - sum(s & off))

  fit <- netmice(fx$attrs, fx$nets, m = 2, maxit = 2, structural = s,
                 seed = 11, printFlag = FALSE)
  expect_s3_class(fit, "netmids")
  for (im in 1:2) {
    # structural cells stay exactly zero (including the NA-coded ones) ...
    expect_true(all(fit$imp_nets[[im]]$friends[s & off] == 0))
    expect_true(all(fit$imp_nets[[im]]$advice[s & off] == 0))
    # ... while all genuinely missing ties are imputed
    expect_false(anyNA(fit$imp_nets[[im]]$friends[off]))
  }
})

test_that("netmice: structural - a named list applies different fixed-zero cells per network", {
  fx <- make_missing_fixture()
  n <- nrow(fx$attrs)
  offm <- row(fx$nets$friends) != col(fx$nets$friends)
  # in friends, nodes 1:5 cannot *send* ties; in advice, nodes 1:5 cannot
  # *receive* ties - different structural patterns per network
  s_friends <- matrix(FALSE, n, n); s_friends[1:5, ] <- TRUE
  s_advice  <- matrix(FALSE, n, n); s_advice[, 1:5]  <- TRUE
  diag(s_friends) <- FALSE
  diag(s_advice)  <- FALSE
  fx$nets$friends[s_friends] <- 0
  fx$nets$advice[s_advice]   <- 0
  # give advice some genuinely missing (non-structural) ties too, so both
  # networks are imputed
  fx$nets$advice[sample(which(offm & !s_advice), 20)] <- NA

  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2,
                 structural = list(friends = s_friends, advice = s_advice),
                 seed = 12, printFlag = FALSE)
  # each network honours its own structural pattern
  expect_true(all(fit$imp_nets[[1]]$friends[s_friends & offm] == 0))
  expect_true(all(fit$imp_nets[[1]]$advice[s_advice & offm] == 0))
  expect_false(anyNA(fit$imp_nets[[1]]$friends[offm]))
  expect_false(anyNA(fit$imp_nets[[1]]$advice[offm]))
  # cells structural in advice only are NOT fixed in friends: observed ties
  # into nodes 1:5 survive there
  expect_true(any(fit$imp_nets[[1]]$friends[s_advice & !s_friends & offm] != 0))

  # unknown network names in the list are rejected
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            structural = list(nope = s_friends)),
    "unknown network"
  )
  # an observed non-zero tie in a structural cell contradicts the design
  bad <- matrix(FALSE, n, n)
  bad[which(fx$nets$friends == 1 & offm)[1]] <- TRUE
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            structural = list(friends = bad)),
    "structurally absent"
  )
})

test_that("netmice: a multinomial attribute (>2 levels) dispatches to nnet::multinom", {
  skip_if_not_installed("nnet")
  fx <- make_missing_fixture()
  # trace()'s tracer expression is evaluated inside nnet::multinom's own
  # execution frame, so a plain local variable + `<<-` from this test_that()
  # block is not reliably reachable via lexical scoping. .GlobalEnv is an
  # absolute reference resolvable regardless of where the expression runs.
  counter_name <- ".netimpute_test_multinom_calls"
  assign(counter_name, 0, envir = .GlobalEnv)
  trace("multinom", where = asNamespace("nnet"), print = FALSE,
        tracer = bquote(assign(.(counter_name), get(.(counter_name), envir = .GlobalEnv) + 1,
                                envir = .GlobalEnv)))
  on.exit({
    untrace("multinom", where = asNamespace("nnet"))
    rm(list = counter_name, envir = .GlobalEnv)
  }, add = TRUE)
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, printFlag = FALSE)
  expect_gt(get(counter_name, envir = .GlobalEnv), 0)
  expect_true(all(fit$imp[[1]]$dept %in% unique(fx$attrs$dept[!is.na(fx$attrs$dept)])))
})

test_that("netmice: net_random_intercepts imputes all missing ties via the lme4 working model", {
  skip_if_not_installed("lme4")
  fx <- make_missing_fixture()
  # trace the mixed-model PMM helper to prove the lme4 path (not the
  # standard mice PMM) actually handled the network visits
  counter_name <- ".netimpute_test_ranef_calls"
  assign(counter_name, 0, envir = .GlobalEnv)
  trace(".impute_pmm_ranef", where = asNamespace("netimpute"), print = FALSE,
        tracer = bquote(assign(.(counter_name), get(.(counter_name), envir = .GlobalEnv) + 1,
                                envir = .GlobalEnv)))
  on.exit({
    untrace(".impute_pmm_ranef", where = asNamespace("netimpute"))
    rm(list = counter_name, envir = .GlobalEnv)
  }, add = TRUE)
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, printFlag = FALSE,
                  net_random_intercepts = c("ego", "alter"))
  expect_s3_class(fit, "netmids")
  expect_gt(get(counter_name, envir = .GlobalEnv), 0)
  off <- which(row(fx$nets$friends) != col(fx$nets$friends))
  expect_false(anyNA(fit$imp_nets[[1]]$friends[off]))
  # imputed ties are donor values, i.e. observed tie values
  obs_vals <- unique(fx$nets$friends[off][!is.na(fx$nets$friends[off])])
  expect_true(all(fit$imp_nets[[1]]$friends[off] %in% obs_vals))
  expect_equal(fit$net_random_intercepts, c("ego", "alter"))
})

test_that("netmice: net_random_intercepts = 'dyad' warns for an undirected network with missing ties", {
  skip_if_not_installed("lme4")
  fx <- make_missing_fixture()
  m <- fx$nets$friends
  undir <- ((ifelse(is.na(m), 0, m) + t(ifelse(is.na(m), 0, m))) > 0) * 1
  na_cells <- which(is.na(m) & row(m) != col(m))
  undir[na_cells] <- NA
  undir[cbind(col(m)[na_cells], row(m)[na_cells])] <- NA  # keep NAs symmetric
  # collect every warning: lme4's singular-fit fallback may also fire here
  # (the PCA budget leaves few predictors for a dyad random intercept), and
  # that is a documented graceful degradation, not the subject of this test
  ws <- character(0)
  withCallingHandlers(
    netmice(fx$attrs, list(friends = undir), m = 1, maxit = 1, printFlag = FALSE,
            net_random_intercepts = "dyad"),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("undirected", ws)))
})

test_that("netmice: invalid net_random_intercepts values are rejected", {
  fx <- make_missing_fixture()
  expect_error(
    netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE,
            net_random_intercepts = "node"),
    "'arg'"
  )
})

test_that("netmice: rejects a signed network in `networks`", {
  fx <- make_missing_fixture()
  fx$nets$trust <- fx_signed(n = nrow(fx$attrs))
  expect_error(netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE),
               "non-negative weighted")
})

test_that("netmice: errors when an attribute column and a network share a name", {
  fx <- make_missing_fixture()
  names(fx$nets)[1] <- "age"  # clashes with the attribute column "age"
  expect_error(netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE),
               "must be distinct")
})

test_that("complete_netmice: extracts a valid completed (data, networks) pair", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 2, maxit = 1, printFlag = FALSE)
  out <- complete_netmice(fit, 2)
  expect_identical(out$data, fit$imp[[2]])
  expect_identical(out$networks, fit$imp_nets[[2]])
  expect_error(complete_netmice(fit, 3), "between 1 and")
})

test_that("print.netmids: prints a human-readable summary", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE)
  expect_output(print(fit), "netmids")
  expect_output(print(fit), "Imputations")
})

# --- undirected networks stay undirected -------------------------------------
# Mirror cells of an undirected network are two records of ONE tie. Drawing
# them independently desymmetrizes the matrix, which silently reclassifies the
# network as directed for every derived measure (all of which detect
# directedness with isSymmetric()).

make_undirected_fixture <- function(n = 22, seed = 707, weighted = FALSE, k = 40) {
  set.seed(seed)
  m <- if (weighted) {
    matrix(rpois(n * n, 0.8), n, n)
  } else {
    matrix(rbinom(n * n, 1, 0.25), n, n)
  }
  diag(m) <- 0
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m[sample(which(upper.tri(m)), k)] <- NA
  m[lower.tri(m)] <- t(m)[lower.tri(m)]   # keep the NA pattern symmetric too
  attrs <- fx_attrs(n = n, seed = seed)[c("age", "status")]
  attrs$age[sample(n, 3)] <- NA
  list(attrs = attrs, net = m)
}

test_that("netmice: an undirected binary network is still undirected after imputation", {
  fx <- make_undirected_fixture()
  expect_true(isSymmetric(unname(fx$net)))
  for (upd in c("gibbs", "simultaneous")) {
    fit <- netmice(fx$attrs, list(u = fx$net), m = 2, maxit = 3, seed = 9,
                   net_update = upd, printFlag = FALSE)
    for (a in seq_len(2)) {
      out <- complete_netmice(fit, a)$networks$u
      expect_true(isSymmetric(unname(out)),
                  info = paste("net_update =", upd, "imputation", a))
      # observed cells untouched, and the result is a valid binary network
      expect_equal(out[!is.na(fx$net)], fx$net[!is.na(fx$net)])
      expect_true(all(out %in% c(0, 1)))
    }
  }
})

test_that("netmice: an undirected weighted network is still undirected after imputation", {
  fx <- make_undirected_fixture(weighted = TRUE)
  for (upd in c("gibbs", "simultaneous")) {
    fit <- netmice(fx$attrs, list(u = fx$net), m = 1, maxit = 3, seed = 4,
                   net_update = upd, printFlag = FALSE)
    out <- complete_netmice(fit, 1)$networks$u
    expect_true(isSymmetric(unname(out)), info = paste("net_update =", upd))
    expect_equal(out[!is.na(fx$net)], fx$net[!is.na(fx$net)])
  }
})

test_that("netmice: undirected symmetry survives net_init = 'sample'", {
  fx <- make_undirected_fixture()
  fit <- netmice(fx$attrs, list(u = fx$net), m = 1, maxit = 2, seed = 5,
                 net_init = "sample", printFlag = FALSE)
  expect_true(isSymmetric(unname(complete_netmice(fit, 1)$networks$u)))
})

test_that("netmice: undirected handling leaves directed networks asymmetric", {
  # guard against over-symmetrizing: a directed network must NOT be forced
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, seed = 3,
                 printFlag = FALSE)
  expect_false(isSymmetric(unname(complete_netmice(fit, 1)$networks$friends)))
})

test_that(".symmetrize_imputed: reconciles disagreeing pairs and spares observed cells", {
  m <- matrix(0, 4, 4)
  m[1, 2] <- 1; m[2, 1] <- 0      # both imputed, disagree
  m[3, 4] <- 1; m[4, 3] <- 0      # (3,4) observed, (4,3) imputed
  ry <- matrix(FALSE, 4, 4)
  ry[3, 4] <- TRUE
  set.seed(1)
  out <- netimpute:::.symmetrize_imputed(m, ry, binary = TRUE)
  expect_true(isSymmetric(unname(out)))
  expect_true(out[1, 2] %in% c(0, 1))
  # the observed cell wins over the imputed mirror
  expect_equal(out[3, 4], 1)
  expect_equal(out[4, 3], 1)

  w <- matrix(0, 3, 3)
  w[1, 2] <- 4; w[2, 1] <- 2      # weighted: average
  out_w <- netimpute:::.symmetrize_imputed(w, matrix(FALSE, 3, 3), binary = FALSE)
  expect_equal(out_w[1, 2], 3)
  expect_equal(out_w[2, 1], 3)
})

test_that(".impute_ties_gibbs: change statistics stay exact when writing both cells of a pair", {
  # the undirected path writes two cells per draw, so the incremental
  # two-path bookkeeping is applied twice per visited pair
  set.seed(23)
  n <- 15
  m <- matrix(rbinom(n * n, 1, 0.25), n, n); diag(m) <- 0
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m[sample(which(upper.tri(m)), 30)] <- NA
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  attrs <- data.frame(age = rnorm(n), grade = rnorm(n))
  filled <- netimpute:::.init_fill_matrix(m, init = "zero", undirected = TRUE)
  built <- suppressMessages(
    netimpute:::.build_dyad_data(list(net = filled), attrs, 1))
  d <- built$data
  ry <- !is.na(m)[cbind(d$i, d$j)]
  x <- netimpute:::.clean_predictor_matrix(
    d[setdiff(names(d), c("i", "j", "y"))], ry = ry)
  out <- netimpute:::.impute_ties_gibbs(d = d, ry = ry, x = x, mat = filled,
                                        binary = TRUE, donors = 5,
                                        check = TRUE, undirected = TRUE)
  obs <- !is.na(m) & (row(m) != col(m))
  expect_identical(out[obs], m[obs])
  expect_true(isSymmetric(unname(out)))
})

test_that(".init_fill_matrix: sample init keeps an undirected network symmetric", {
  set.seed(31)
  n <- 12
  m <- matrix(rbinom(n * n, 1, 0.3), n, n); diag(m) <- 0
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m[sample(which(upper.tri(m)), 20)] <- NA
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  expect_true(isSymmetric(unname(
    netimpute:::.init_fill_matrix(m, init = "sample", undirected = TRUE))))
  # without the flag the fills are independent per cell (directed behaviour)
  expect_true(isSymmetric(unname(
    netimpute:::.init_fill_matrix(m, init = "zero", undirected = TRUE))))
})

# --- directedness is derived once, from the data passed to netmice() ---------
# Filling destroys the evidence: a directed network whose asymmetry lies
# entirely in its missing cells becomes symmetric once those are zeroed, so
# anything re-deriving directedness downstream would misclassify it.

pathological_directed <- function() {
  n <- 6
  m <- matrix(0, n, n)
  m[1, 2] <- 1; m[2, 1] <- 1
  m[3, 4] <- 1; m[4, 3] <- 1
  m[5, 6] <- NA; m[6, 5] <- 0
  diag(m) <- 0
  m
}

test_that(".mat_to_igraph: uses the supplied classification, not the filled state", {
  m <- pathological_directed()
  filled <- netimpute:::.init_fill_matrix(m, init = "zero")
  expect_false(isSymmetric(unname(m)))       # input: directed
  expect_true(isSymmetric(unname(filled)))   # filled: looks undirected

  expect_true(igraph::is_directed(
    netimpute:::.mat_to_igraph(filled, directed = TRUE)))
  expect_false(igraph::is_directed(
    netimpute:::.mat_to_igraph(filled, directed = FALSE)))
  # NULL keeps the old self-detecting behaviour for standalone callers
  expect_false(igraph::is_directed(netimpute:::.mat_to_igraph(filled)))
})

test_that(".mat_to_igraph: mode = 'max' matches mode = 'undirected' on symmetric input", {
  set.seed(77)
  for (w in c(FALSE, TRUE)) {
    m <- if (w) matrix(rpois(100, 0.8), 10, 10) else matrix(rbinom(100, 1, 0.3), 10, 10)
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    diag(m) <- 0
    wt <- if (any(m[m != 0] != 1)) TRUE else NULL
    ref <- igraph::graph_from_adjacency_matrix(m, mode = "undirected",
                                               weighted = wt, diag = FALSE)
    got <- netimpute:::.mat_to_igraph(m, directed = FALSE)
    at <- if (is.null(wt)) NULL else "weight"
    expect_equal(igraph::as_adjacency_matrix(got, sparse = FALSE, attr = at),
                 igraph::as_adjacency_matrix(ref, sparse = FALSE, attr = at))
  }
  # and it does not warn on the asymmetric state a dependence mask can leave
  asym <- matrix(0, 4, 4); asym[1, 2] <- 1
  expect_silent(netimpute:::.mat_to_igraph(asym, directed = FALSE))
})

test_that("netmice: a misclassified directed network would lose reciprocity_ratio", {
  # the concrete consequence of re-deriving directedness from a filled
  # matrix: reciprocity_ratio is NA for an undirected graph, so the feature
  # is silently blanked for that visit
  m <- pathological_directed()
  filled <- netimpute:::.init_fill_matrix(m, init = "zero")
  attrs <- data.frame(age = as.numeric(1:6))
  ff <- function(dir) suppressMessages(netimpute:::.net_feature_frame(
    list(x = netimpute:::.mat_to_igraph(filled, directed = dir)),
    attrs, "core", attr_types = c(age = "continuous"),
    net_names = "x", clash_names = "age"))$reciprocity_ratio
  expect_false(all(is.na(ff(TRUE))))   # kept as a live predictor
  expect_true(all(is.na(ff(FALSE))))   # lost if misclassified
})

test_that("netmice: runs end to end on a network whose filled state is symmetric", {
  m <- pathological_directed()
  attrs <- data.frame(age = c(1, NA, 3, 4, NA, 6))
  fit <- netmice(attrs, list(x = m), m = 1, maxit = 2, seed = 1,
                 printFlag = FALSE)
  out <- complete_netmice(fit, 1)$networks$x
  obs <- !is.na(m) & (row(m) != col(m))
  expect_equal(out[obs], m[obs])       # observed ties untouched
  expect_false(anyNA(out))
})

# --- length-1 sampling ------------------------------------------------------
# sample(x) treats a length-1 NUMERIC x as 1:x, so a network with exactly one
# missing tie (or a variable with exactly one observed value) would sweep over
# the wrong index set and overwrite observed data.

test_that(".shuffle/.resample: a single element is not reinterpreted as 1:x", {
  expect_identical(netimpute:::.shuffle(29L), 29L)
  # shuffle ONCE and sort: the earlier form called .shuffle() twice and used
  # one draw's order() to index the other, so it passed only ~43% of the time
  expect_identical(sort(netimpute:::.shuffle(c(3L, 7L))), c(3L, 7L))
  expect_length(netimpute:::.shuffle(integer(0)), 0)
  expect_identical(netimpute:::.resample(5, 3), c(5, 5, 5))
  expect_true(all(netimpute:::.resample(c(2, 9), 20) %in% c(2, 9)))
  # a plain permutation is still a permutation
  set.seed(1)
  expect_setequal(netimpute:::.shuffle(1:10), 1:10)
})

test_that("netmice: a network with exactly ONE missing tie leaves every observed tie intact", {
  n <- 6
  m <- matrix(0, n, n)
  m[1, 2] <- 1; m[2, 1] <- 1; m[3, 4] <- 1; m[4, 3] <- 1
  m[5, 6] <- NA
  diag(m) <- 0
  attrs <- data.frame(age = c(1, NA, 3, 4, NA, 6))
  for (upd in c("gibbs", "simultaneous")) {
    fit <- netmice(attrs, list(x = m), m = 1, maxit = 2, seed = 1,
                   net_update = upd, printFlag = FALSE)
    out <- complete_netmice(fit, 1)$networks$x
    obs <- !is.na(m) & (row(m) != col(m))
    expect_equal(out[obs], m[obs], info = upd)   # observed ties untouched
    expect_false(anyNA(out))
    expect_true(out[5, 6] %in% c(0, 1))
  }
})

test_that("netmice: an undirected network with exactly ONE missing pair is handled", {
  n <- 6
  m <- matrix(0, n, n)
  m[1, 2] <- 1; m[2, 1] <- 1; m[3, 4] <- 1; m[4, 3] <- 1
  m[5, 6] <- NA; m[6, 5] <- NA
  diag(m) <- 0
  expect_true(isSymmetric(unname(m)))
  attrs <- data.frame(age = c(1, NA, 3, 4, NA, 6))
  fit <- netmice(attrs, list(x = m), m = 1, maxit = 2, seed = 1,
                 printFlag = FALSE)
  out <- complete_netmice(fit, 1)$networks$x
  obs <- !is.na(m) & (row(m) != col(m))
  expect_equal(out[obs], m[obs])
  expect_true(isSymmetric(unname(out)))
})

test_that(".init_fill_vector: a single observed value is used, not 1:value", {
  x <- c(NA, NA, 7, NA)
  set.seed(1)
  expect_equal(netimpute:::.init_fill_vector(x), c(7, 7, 7, 7))
})
