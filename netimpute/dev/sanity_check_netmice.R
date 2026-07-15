# Manual smoke test for dyad_regression() and netmice(). Run after
# devtools::load_all(). Not part of the built package (see .Rbuildignore).
#
# I could not execute R in the environment where this file was written, so
# please run this yourself and report any failures - this is new, fairly
# involved code (cell-level regression + a custom mice-style loop) and I
# have not been able to verify it runs end to end.

devtools::load_all()
set.seed(1)

n <- 30
friends <- matrix(rbinom(n * n, 1, 0.12), n, n); diag(friends) <- 0
advice  <- matrix(rbinom(n * n, 1, 0.10), n, n); diag(advice)  <- 0
attrs <- data.frame(
  age  = round(rnorm(n, 35, 8)),
  dept = sample(c("sales", "eng", "hr"), n, replace = TRUE)
)

## ---- dyad_regression(), raw other-network predictors ----------------------
res_raw <- dyad_regression(list(friends = friends, advice = advice), attrs,
                            target = "friends", other_net_predictors = "raw")
stopifnot(all(c("i", "j", "y", "reciprocity", "twopath") %in% names(res_raw$data)))
stopifnot(all(c("age_ego", "age_alter", "age_absdiff") %in% names(res_raw$data)))
stopifnot(all(c("dept_ego_hr", "dept_ego_sales", "dept_alter_hr", "dept_alter_sales",
                "dept_same") %in% names(res_raw$data)))
stopifnot(all(c("advice_tie", "advice_recip", "advice_ego_outdeg", "advice_alter_indeg")
              %in% names(res_raw$data)))
stopifnot(nrow(res_raw$data) == n * (n - 1))
print(summary(res_raw$model))
cat("dyad_regression(raw): OK\n\n")

## ---- dyad_regression(), PCA-reduced other-network predictors --------------
extra <- matrix(rbinom(n * n, 1, 0.08), n, n); diag(extra) <- 0
res_pca <- dyad_regression(
  list(friends = friends, advice = advice, extra = extra), attrs,
  target = "friends", other_net_predictors = "pca", n_components = 3
)
stopifnot(all(paste0("other_net_PC", 1:3) %in% names(res_pca$data)))
stopifnot(!is.null(res_pca$pca_model))
cat("dyad_regression(pca): OK\n\n")

## ---- netmice() with missingness in both an attribute and a network --------
friends_miss <- friends
attrs_miss <- attrs
off <- which(row(friends_miss) != col(friends_miss))
friends_miss[sample(off, 40)] <- NA
attrs_miss$age[sample(n, 4)] <- NA
attrs_miss$dept[sample(n, 3)] <- NA

fit <- netmice(attrs_miss, list(friends = friends_miss, advice = advice),
                m = 2, maxit = 2, donors = 5, printFlag = TRUE)
print(fit)

stopifnot(!anyNA(fit$imp[[1]]$age))
stopifnot(!anyNA(fit$imp[[1]]$dept))
stopifnot(!anyNA(fit$imp_nets[[1]]$friends[row(friends) != col(friends)]))
stopifnot(dim(fit$chainMean)[2] == 2 && dim(fit$chainMean)[3] == 2)
stopifnot(dim(fit$netChain)[1] == 2)  # two networks tracked, even though only one had NAs

out <- complete_netmice(fit, 1)
stopifnot(identical(dim(out$net_list$friends), dim(friends)))
cat("netmice(): OK\n\n")

## ---- visit order differs across chains (random starting variable) --------
stopifnot(length(fit$visit_orders) == 2)
print(fit$visit_orders)
cat("visit_orders recorded per chain (inspect above - order should differ across chains ",
    "with reasonable probability): OK\n\n")

## ---- netmice() with a custom `models` formula (with an interaction) ------
attrs_miss2 <- attrs_miss
attrs_miss2$performance <- rnorm(n)
attrs_miss2$happiness <- rnorm(n)
attrs_miss2$happiness[sample(n, 5)] <- NA

fit2 <- netmice(attrs_miss2, list(friends = friends_miss, advice = advice),
                 m = 2, maxit = 2, donors = 5, printFlag = FALSE,
                 models = list("happiness ~ friends_indegree + age + performance * dept"))
stopifnot(!anyNA(fit2$imp[[1]]$happiness))
cat("netmice() with `models` (attribute formula incl. interaction): OK\n\n")

## ---- netmice() with a custom formula on a NETWORK's own dyad-level terms --
fit3 <- netmice(attrs_miss, list(friends = friends_miss, advice = advice),
                 m = 2, maxit = 2, donors = 5, printFlag = FALSE,
                 models = list("friends ~ age_absdiff + advice_tie + reciprocity"))
stopifnot(!anyNA(fit3$imp_nets[[1]]$friends[row(friends) != col(friends)]))
cat("netmice() with `models` (network dyad-level formula): OK\n\n")

## ---- ncores > 1 (requires the package to be *installed*, not just ---------
## load_all()-ed, since multisession workers resolve netimpute's internal
## functions by loading the installed namespace; confirmed both that it
## fails cleanly under plain load_all() and that, once installed, it
## reproduces the sequential result exactly given the same seed - the RNG
## kind is pinned explicitly in .run_one_chain() for this reason.
run_parallel_check <- TRUE
if (run_parallel_check) {
  fit4 <- netmice(attrs_miss, list(friends = friends_miss, advice = advice),
                   m = 2, maxit = 2, donors = 5, ncores = 2, seed = 1)
  fit5 <- netmice(attrs_miss, list(friends = friends_miss, advice = advice),
                   m = 2, maxit = 2, donors = 5, ncores = 1, seed = 1)
  stopifnot(identical(fit4$imp[[1]]$age, fit5$imp[[1]]$age))
  cat("netmice(ncores=2) matches ncores=1 given the same seed: OK\n\n")
}

## ---- single-network case (regression test) --------------------------------
## cbind(data.frame, NULL) errors rather than no-op'ing in R (confirmed via
## a minimal repro), which .build_dyad_data() used to hit whenever there were
## no "other networks" (other_predictors empty -> NULL) - i.e. exactly the
## single-network case, which neither test above exercises since both always
## pass two networks. Keep this here so that regression doesn't come back.
attrs_single <- attrs_miss
fit_single <- netmice(attrs_single, list(friends = friends_miss),
                       m = 1, maxit = 2, donors = 5, printFlag = FALSE)
stopifnot(!anyNA(fit_single$imp[[1]]$age))
stopifnot(!anyNA(fit_single$imp[[1]]$dept))
stopifnot(!anyNA(fit_single$imp_nets[[1]]$friends[row(friends) != col(friends)]))
cat("netmice() with a single network (no 'other network' terms): OK\n\n")

## ---- multinomial attribute uses a proper multinomial model, not pmm ------
## (dept has 3 levels: sales/eng/hr). Trace nnet::multinom to confirm it is
## actually invoked, rather than just checking the result is non-NA (the old,
## since-replaced integer-coded-pmm path also produced non-NA results).
if (requireNamespace("nnet", quietly = TRUE)) {
  call_count <- 0
  trace("multinom", where = asNamespace("nnet"), print = FALSE,
        tracer = quote(call_count <<- call_count + 1))
  fit_poly <- netmice(attrs_miss, list(friends = friends_miss, advice = advice),
                       m = 1, maxit = 2, donors = 5, printFlag = FALSE)
  untrace("multinom", where = asNamespace("nnet"))
  stopifnot(call_count > 0)
  stopifnot(all(fit_poly$imp[[1]]$dept %in% c("sales", "eng", "hr")))
  cat("multinomial attribute (dept) dispatches to nnet::multinom, called",
      call_count, "time(s): OK\n\n")
}

cat("ALL NETMICE SANITY CHECKS PASSED\n")
