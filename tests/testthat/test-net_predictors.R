nets_a <- fx_nets(n = 20)
nets_b <- fx_nets(n = 16)
attrs_a <- fx_attrs(n = 20, seed = 301)
attrs_b <- fx_attrs(n = 16, seed = 302)

# One net_predictors() call pools measures across a *list* of (network,
# attributes) pairs, one entry per network - here spanning heterogeneous tie
# types (binary/weighted) and heterogeneous sizes, mirroring "several
# networks, some binary, some weighted" from the brief. Signed networks are
# tested separately below (and in test-net_measures_core.R) since they now
# error at the coercion step.
pooled_nets <- list(a_friends = nets_a$friends_bin, a_advice = nets_a$advice_weighted,
                     b_friends = nets_b$friends_bin, b_advice = nets_b$advice_weighted)
pooled_attrs <- list(attrs_a[c("age", "status", "dept")], attrs_a[c("age", "status", "dept")],
                      attrs_b[c("age", "status", "dept")], attrs_b[c("age", "status", "dept")])

test_that("net_predictors: output = 'measures' stacks across heterogeneous networks", {
  out <- net_predictors(pooled_nets, pooled_attrs, measure_set = "core", output = "measures")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 20 + 20 + 16 + 16)
  expect_true(all(c("network_id", "node_id") %in% names(out)))
  expect_setequal(unique(out$network_id), names(pooled_nets))
})

test_that("net_predictors: measure_set = 'full' produces more columns than 'core'", {
  out_core <- net_predictors(pooled_nets, pooled_attrs, measure_set = "core", output = "measures")
  out_full <- net_predictors(pooled_nets, pooled_attrs, measure_set = "full", output = "measures")
  expect_gt(ncol(out_full), ncol(out_core))
})

test_that("net_predictors: attr_types override is passed through", {
  out <- net_predictors(pooled_nets, pooled_attrs, attr_types = c(age = "continuous"),
                         output = "measures")
  expect_true("age_diff_out" %in% names(out))
})

test_that("net_predictors: use_sna is passed through for measure_set = 'full'", {
  out_no_sna <- net_predictors(pooled_nets, pooled_attrs, measure_set = "full",
                                output = "measures", use_sna = FALSE)
  expect_true(all(is.na(out_no_sna$gilschmidt)))
})

test_that("net_predictors: output = 'pca' returns PCA$n PCs and a prcomp model", {
  res <- net_predictors(pooled_nets, pooled_attrs, output = "pca", PCA = list(n = 4))
  expect_type(res, "list")
  expect_true(all(c("predictors", "pca_model", "raw_measures") %in% names(res)))
  expect_s3_class(res$pca_model, "prcomp")
  expect_equal(sum(grepl("^PC", names(res$predictors))), 4)
  expect_equal(nrow(res$predictors), nrow(res$raw_measures))
})

test_that("net_predictors: the component count is capped at the available columns", {
  res <- net_predictors(pooled_nets, pooled_attrs, output = "pca", PCA = list(n = 9999))
  expect_lt(sum(grepl("^PC", names(res$predictors))), 9999)
})

test_that("net_predictors: output = 'both' requires keep_vars and preserves them untransformed", {
  expect_error(net_predictors(pooled_nets, pooled_attrs, output = "both"),
               "keep_vars")
  res <- net_predictors(pooled_nets, pooled_attrs, output = "both",
                         keep_vars = c("total_degree"), measure_set = "full", PCA = list(n = 3))
  expect_true("total_degree" %in% names(res$predictors))
  expect_equal(sum(grepl("^PC", names(res$predictors))), 3)
  expect_equal(res$predictors$total_degree, res$raw_measures$total_degree)
})

test_that("net_predictors: residualize_kept changes the PCA input relative to plain exclusion", {
  res_plain <- net_predictors(pooled_nets, pooled_attrs, output = "both", measure_set = "full",
                               keep_vars = c("total_degree"), PCA = list(n = 3),
                               residualize_kept = FALSE)
  res_resid <- net_predictors(pooled_nets, pooled_attrs, output = "both", measure_set = "full",
                               keep_vars = c("total_degree"), PCA = list(n = 3),
                               residualize_kept = TRUE)
  expect_false(isTRUE(all.equal(res_plain$pca_model$rotation, res_resid$pca_model$rotation)))
})

test_that("net_predictors: impute_na = 'complete_cases' drops incomplete rows relative to 'mean'", {
  res_mean <- net_predictors(pooled_nets, pooled_attrs, output = "pca", impute_na = "mean")
  res_cc   <- net_predictors(pooled_nets, pooled_attrs, output = "pca", impute_na = "complete_cases")
  expect_lte(nrow(res_cc$predictors), nrow(res_mean$predictors))
})

test_that("net_predictors: scale_pca affects the PCA solution", {
  res_scaled   <- net_predictors(pooled_nets, pooled_attrs, output = "pca", scale_pca = TRUE)
  res_unscaled <- net_predictors(pooled_nets, pooled_attrs, output = "pca", scale_pca = FALSE)
  expect_false(isTRUE(all.equal(res_scaled$pca_model$sdev, res_unscaled$pca_model$sdev)))
})

test_that("net_predictors: id_col is passed through to the per-network measure calls", {
  m <- nets_a$friends_bin
  dimnames(m) <- list(attrs_a$id, attrs_a$id)
  shuffled <- attrs_a[sample(nrow(attrs_a)), c("id", "age", "status", "dept")]
  out <- net_predictors(list(a = m), list(shuffled), id_col = "id", output = "measures")
  out_ref <- net_predictors(list(a = m), list(attrs_a[c("age", "status", "dept")]), output = "measures")
  expect_equal(out$outdegree, out_ref$outdegree)
})

test_that("net_predictors: errors on mismatched networks/attr_list length", {
  expect_error(net_predictors(pooled_nets, pooled_attrs[1:2], output = "measures"),
               "same length")
})

test_that("net_predictors: errors when attr_list schemas differ", {
  bad_attrs <- pooled_attrs
  bad_attrs[[1]] <- bad_attrs[[1]][c("age", "status")]  # drop dept
  expect_error(net_predictors(pooled_nets, bad_attrs, output = "measures"),
               "same attribute set")
})

test_that("net_predictors: errors on unknown keep_vars", {
  expect_error(
    net_predictors(pooled_nets, pooled_attrs, output = "both", keep_vars = "not_a_real_column"),
    "not found among computed measures"
  )
})

test_that("net_predictors: rejects a signed network among the list", {
  bad_nets <- pooled_nets
  bad_nets$a_friends <- fx_signed(n = 20)
  expect_error(net_predictors(bad_nets, pooled_attrs, output = "measures"),
               "non-negative weighted")
})
