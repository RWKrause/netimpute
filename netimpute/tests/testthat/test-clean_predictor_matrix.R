# Guards against mice's "system is exactly singular" crash: predictor
# columns that are constant on the OBSERVED rows (where the model is fit)
# must never reach mice's ridge solver.

test_that(".clean_predictor_matrix drops columns constant on observed rows", {
  set.seed(1)
  ry <- rep(c(TRUE, FALSE), each = 50)
  x <- cbind(good = rnorm(100),
             dead_on_observed = c(rep(0, 50), rbinom(50, 1, 0.5)),
             const_on_observed = c(rep(3, 50), rnorm(50)))
  out <- netimpute:::.clean_predictor_matrix(x, ry = ry)
  expect_identical(colnames(out), "good")
  # without ry, the old all-rows behavior is unchanged
  out_all <- netimpute:::.clean_predictor_matrix(x)
  expect_identical(colnames(out_all), colnames(x))
})

test_that(".clean_predictor_matrix PCA branch never returns rank-deficient columns", {
  set.seed(2)
  # 8 columns but true rank 3: components 4+ have (numerically) zero scores
  x <- matrix(rnorm(300), 100, 3) %*% matrix(rnorm(24), 3, 8)
  colnames(x) <- paste0("c", 1:8)
  out <- netimpute:::.clean_predictor_matrix(x, max_cols = 6)
  expect_lte(ncol(out), 3)
  expect_true(all(apply(out, 2, stats::var) > 1e-8))
})

test_that(".clean_predictor_matrix: keep_raw columns bypass the PCA reduction", {
  set.seed(5)
  x <- matrix(rnorm(100 * 10), 100, 10,
              dimnames = list(NULL, c("reciprocity", "twopath", paste0("c", 1:8))))
  out <- netimpute:::.clean_predictor_matrix(
    x, max_cols = 4, keep_raw = c("reciprocity", "twopath"))
  expect_true(all(c("reciprocity", "twopath") %in% colnames(out)))
  # kept columns are passed through untransformed
  expect_identical(out[, "reciprocity"], x[, "reciprocity"])
  expect_identical(out[, "twopath"], x[, "twopath"])
  # the remaining 8 columns are reduced to max_cols components
  expect_equal(sum(grepl("^PC", colnames(out))), 4)
  # keep_raw names absent from x are ignored, and without keep_raw the
  # old everything-into-PCA behavior is unchanged
  out2 <- netimpute:::.clean_predictor_matrix(x, max_cols = 4,
                                              keep_raw = "no_such_column")
  expect_equal(ncol(out2), 4)
})

test_that("netmice survives another network's ties living only on the target's missing dyads", {
  set.seed(3)
  n <- 30
  A <- matrix(rbinom(n * n, 1, .15), n, n); diag(A) <- 0
  off <- which(row(A) != col(A))
  mis <- sample(off, 120)
  A[mis] <- NA
  # B's only ties sit on A's missing cells -> B_tie is all-zero on the rows
  # A's tie model is fit on; this crashed mice's ridge solver before the
  # ry-aware cleaning ("Lapack routine dgesv: system is exactly singular")
  B <- matrix(0, n, n); B[mis[1:40]] <- 1; diag(B) <- 0
  attrs <- data.frame(age = rnorm(n), happy = rnorm(n))
  attrs$happy[sample(n, 4)] <- NA
  fit <- suppressMessages(
    netmice(attrs, list(A = A, B = B), m = 1, maxit = 2,
            seed = 1, printFlag = FALSE))
  expect_s3_class(fit, "netmids")
  expect_false(anyNA(fit$imp[[1]]$happy))
  expect_false(anyNA(fit$imp_nets[[1]]$A))
})

test_that("a failing visit reports which target it failed at", {
  set.seed(4)
  n <- 20
  A <- matrix(rbinom(n * n, 1, .2), n, n); diag(A) <- 0
  A[sample(which(row(A) != col(A)), 20)] <- NA
  attrs <- data.frame(age = rnorm(n), happy = rnorm(n))
  attrs$happy[sample(n, 3)] <- NA
  expect_error(
    suppressMessages(
      netmice(attrs, list(A = A), m = 1, maxit = 1, seed = 1,
              printFlag = FALSE,
              models = list("happy ~ no_such_predictor"))),
    "imputation failed at target 'happy'"
  )
})
