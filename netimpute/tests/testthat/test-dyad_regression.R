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

test_that(".build_dyad_data: dyad-level values match a hand-built reference (ego = i, alter = j, same = shared category)", {
  m <- matrix(0, 3, 3); m[1, 2] <- 1
  age <- c(10, 20, 30)
  grp <- c("A", "B", "A")
  built <- suppressMessages(.build_dyad_data(
    list(friends = m),
    data.frame(age = age, grp = grp),
    target_idx = 1,
    attr_types = c(age = "continuous", grp = "multinomial")
  ))
  d <- built$data
  # off-diagonal dyads in column-major order
  expect_equal(d$i, c(2L, 3L, 1L, 3L, 1L, 2L))
  expect_equal(d$j, c(1L, 1L, 2L, 2L, 3L, 3L))
  # ego is the sender's (i's) value, alter the receiver's (j's)
  expect_equal(d$age_ego,   age[d$i])
  expect_equal(d$age_alter, age[d$j])
  expect_equal(d$age_absdiff, abs(age[d$i] - age[d$j]))
  # _same must be 1 exactly when i and j hold the identical category, else 0:
  # nodes 1 and 3 are both "A", node 2 is "B"
  expect_equal(d$grp_same, as.numeric(grp[d$i] == grp[d$j]))
  expect_equal(d$grp_same, c(0, 1, 0, 0, 1, 0))
  # level dummies index the right node too
  expect_equal(d$grp_ego_B,   as.numeric(grp[d$i] == "B"))
  expect_equal(d$grp_alter_B, as.numeric(grp[d$j] == "B"))
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
  # the identifier itself must not leak into the predictor set (it would
  # otherwise expand into n-1 ego + n-1 alter dummies plus id_same)
  expect_identical(names(res$data), names(res_ref$data))
  expect_false(any(grepl("^id_", names(res$data))))
})

test_that("dyad_regression: random_intercepts fits an lme4 model with the right grouping factors", {
  skip_if_not_installed("lme4")
  n <- nrow(nets$friends_bin)
  res <- dyad_regression(list(friends = nets$friends_bin, advice = nets$advice_weighted),
                          attrs, target = "friends",
                          random_intercepts = c("ego", "alter"))
  expect_s4_class(res$model, "lmerMod")
  expect_setequal(names(lme4::ranef(res$model)), c(".ego", ".alter"))
  expect_equal(nrow(lme4::ranef(res$model)$.ego), n)
  expect_equal(nrow(lme4::ranef(res$model)$.alter), n)
  # grouping columns are returned so predict() works on res$data
  expect_true(all(c(".ego", ".alter", ".dyad") %in% names(res$data)))
  pred <- predict(res$model, newdata = res$data, allow.new.levels = TRUE)
  expect_length(pred, n * (n - 1))
})

test_that("dyad_regression: a dyad random intercept on a directed target has one group per unordered pair", {
  skip_if_not_installed("lme4")
  n <- nrow(nets$friends_bin)
  expect_message(
    res <- dyad_regression(list(friends = nets$friends_bin), attrs,
                            target = "friends", random_intercepts = "dyad"),
    "reciprocity"
  )
  expect_s4_class(res$model, "lmerMod")
  expect_equal(nrow(lme4::ranef(res$model)$.dyad), n * (n - 1) / 2)
  # the fixed reciprocity term must be gone: with a per-dyad intercept it
  # would make the model exactly singular (b_d = y_ij + y_ji, coef -1)
  expect_false("reciprocity" %in% names(lme4::fixef(res$model)))
})

test_that("dyad_regression: random_intercepts + non-gaussian family uses glmer", {
  skip_if_not_installed("lme4")
  res <- suppressWarnings(
    dyad_regression(list(friends = nets$friends_bin), attrs,
                     target = "friends", family = "binomial",
                     random_intercepts = "ego")
  )
  expect_s4_class(res$model, "glmerMod")
})

test_that("dyad_regression: warns when a dyad intercept is requested for an undirected target", {
  skip_if_not_installed("lme4")
  m <- nets$friends_bin
  m_undir <- ((m + t(m)) > 0) * 1
  # a degenerate dyad intercept on duplicated rows may also trigger lme4
  # convergence warnings - collect everything and look for ours
  w <- capture_warnings(suppressMessages(
    dyad_regression(list(friends = m_undir), attrs, target = "friends",
                     random_intercepts = c("ego", "alter", "dyad"))
  ))
  expect_true(any(grepl("undirected", w)))
})

test_that("dyad_regression: invalid random_intercepts values are rejected", {
  expect_error(
    dyad_regression(list(friends = nets$friends_bin), attrs,
                     target = "friends", random_intercepts = "sender"),
    "'arg'"
  )
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
