nets <- fx_nets(n = 20)
attrs <- fx_attrs(n = 20, seed = 401)[c("age", "status", "dept")]

test_that("dyad_regression: default args (gaussian family) fit and build correctly", {
  res <- dyad_regression(list(friends = nets$friends_bin, advice = nets$advice_weighted),
                          attrs, target = "friends")
  expect_type(res, "list")
  expect_true(all(c("data", "pca_model", "target", "model") %in% names(res)))
  expect_equal(res$target, "friends")
  expect_s3_class(res$model, "glm")
  expect_equal(family(res$model)$family, "gaussian")
  n <- nrow(nets$friends_bin)
  expect_equal(nrow(res$data), n * (n - 1))
  expect_true(all(c("i", "j", "y", "reciprocity", "log_twopath") %in% names(res$data)))
  expect_true(all(c("age_ego", "age_alter", "age_absdiff") %in% names(res$data)))
  expect_true(all(c("status_ego_inactive", "status_alter_inactive", "status_same") %in% names(res$data)) ||
                all(c("status_ego_active", "status_alter_active", "status_same") %in% names(res$data)))
  expect_true(all(c("advice_tie", "advice_recip", "advice_ego_outdeg", "advice_alter_indeg")
                  %in% names(res$data)))
})

test_that("dyad_regression: target may be given by index or by name with identical results", {
  res_name <- dyad_regression(list(friends = nets$friends_bin, advice = nets$advice_weighted),
                               attrs, target = "friends")
  res_idx  <- dyad_regression(list(friends = nets$friends_bin, advice = nets$advice_weighted),
                               attrs, target = 1)
  expect_equal(res_name$data$y, res_idx$data$y)
  expect_equal(coef(res_name$model), coef(res_idx$model))
})

test_that("dyad_regression: family argument is honoured (binomial on a binary target)", {
  res <- dyad_regression(list(friends = nets$friends_bin, advice = nets$advice_weighted),
                          attrs, target = "friends", family = "binomial")
  expect_equal(family(res$model)$family, "binomial")
})

test_that("dyad_regression: other_net_predictors = 'pca' replaces raw other-net terms with PCs", {
  extra <- fx_bin_directed(n = 20, seed = 999)
  res <- dyad_regression(list(friends = nets$friends_bin, advice = nets$advice_weighted, extra = extra),
                          attrs, target = "friends", other_net_predictors = "pca", n_components = 2)
  expect_true(all(paste0("other_net_PC", 1:2) %in% names(res$data)))
  expect_false(any(c("advice_tie", "extra_tie") %in% names(res$data)))
  expect_s3_class(res$pca_model, "prcomp")
})

test_that("dyad_regression: n_components controls the number of PCA terms", {
  extra1 <- fx_bin_directed(n = 20, seed = 991)
  extra2 <- fx_weighted(n = 20, seed = 992)
  res <- dyad_regression(list(friends = nets$friends_bin, advice = nets$advice_weighted,
                               e1 = extra1, e2 = extra2),
                          attrs, target = "friends", other_net_predictors = "pca", n_components = 5)
  expect_equal(sum(grepl("^other_net_PC", names(res$data))), 5)
})

test_that("dyad_regression: id_col aligns shuffled attributes", {
  m <- nets$friends_bin
  dimnames(m) <- list(paste0("n", seq_len(nrow(m))), paste0("n", seq_len(nrow(m))))
  attrs_id <- cbind(id = paste0("n", seq_len(nrow(m))), attrs)
  shuffled <- attrs_id[sample(nrow(attrs_id)), ]
  res <- dyad_regression(list(friends = m), shuffled, target = "friends", id_col = "id")
  res_ref <- dyad_regression(list(friends = m), attrs, target = "friends")
  expect_equal(res$data$y, res_ref$data$y)
  expect_equal(res$data$age_ego, res_ref$data$age_ego)
})

test_that("dyad_regression: fit = FALSE builds data without fitting a model", {
  res <- dyad_regression(list(friends = nets$friends_bin), attrs, target = "friends", fit = FALSE)
  expect_null(res$model)
  expect_true(nrow(res$data) > 0)
})

test_that("dyad_regression: works with a single network (no other-network terms)", {
  res <- dyad_regression(list(friends = nets$friends_bin), attrs, target = "friends")
  expect_s3_class(res$model, "glm")
  expect_false(any(grepl("_tie$|_recip$|_ego_outdeg$|_alter_indeg$", names(res$data))))
})

test_that("dyad_regression: rejects a signed network anywhere in net_list", {
  expect_error(
    dyad_regression(list(friends = nets$friends_bin, trust = fx_signed(n = 20)), attrs,
                     target = "friends"),
    "non-negative weighted"
  )
  expect_error(
    dyad_regression(list(trust = fx_signed(n = 20)), attrs, target = "trust"),
    "non-negative weighted"
  )
})

test_that("dyad_regression: errors clearly when target is not found", {
  expect_error(dyad_regression(list(friends = nets$friends_bin), attrs, target = "nope"),
               "not found")
})
