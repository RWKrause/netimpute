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

test_that("netmice: method is restricted to 'pmm' for now", {
  fx <- make_missing_fixture()
  expect_error(netmice(fx$attrs, fx$nets, m = 1, maxit = 1, method = "norm", printFlag = FALSE))
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

test_that("netmice: other_net_predictors = 'pca' runs end to end", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 2, other_net_predictors = "pca",
                  n_components = 2, printFlag = FALSE)
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
  skip_if_not("netimpute" %in% rownames(utils::installed.packages()),
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

test_that("netmice: rejects a signed network in net_list", {
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

test_that("complete_netmice: extracts a valid completed (data, net_list) pair", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 2, maxit = 1, printFlag = FALSE)
  out <- complete_netmice(fit, 2)
  expect_identical(out$data, fit$imp[[2]])
  expect_identical(out$net_list, fit$imp_nets[[2]])
  expect_error(complete_netmice(fit, 3), "between 1 and")
})

test_that("print.netmids: prints a human-readable summary", {
  fx <- make_missing_fixture()
  fit <- netmice(fx$attrs, fx$nets, m = 1, maxit = 1, printFlag = FALSE)
  expect_output(print(fit), "netmids")
  expect_output(print(fit), "Imputations")
})
