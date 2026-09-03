# The `isolate` flag is the missing-data indicator for the alter-based
# features: those are NA for a node with no alters, and
# .clean_predictor_matrix() fills them with the column mean, which hands an
# isolate the average neighbourhood of the connected nodes. The flag must
# therefore survive both predictor selection and the PCA cap - whenever it
# varies.

make_iso_fixture <- function(n = 40, seed = 12, n_iso = 3) {
  set.seed(seed)
  fr <- matrix(rbinom(n * n, 1, 0.15), n, n)
  diag(fr) <- 0
  iso <- seq.int(n - n_iso + 1, n)
  fr[iso, ] <- 0; fr[, iso] <- 0
  off <- which(row(fr) != col(fr))
  fr[sample(setdiff(off, which(row(fr) %in% iso | col(fr) %in% iso)), 40)] <- NA
  attrs <- data.frame(age = rnorm(n, 35, 8), score = rnorm(n),
                      grp = sample(c("a", "b"), n, TRUE))
  attrs$age[sample(n, 6)] <- NA
  attrs$score[sample(n, 5)] <- NA
  list(net = fr, attrs = attrs)
}

test_that(".with_isolate: adds the flag only when alter-based features are present", {
  expect_true("isolate" %in% .with_isolate(c("homophily", "indegree")))
  expect_true("isolate" %in% .with_isolate("homophily"))
  # no homophily block -> no alter measures -> nothing to indicate
  expect_false("isolate" %in% .with_isolate(c("indegree", "betweenness")))
  # already present: unioned, not duplicated
  expect_equal(sum(.with_isolate(c("homophily", "isolate")) == "isolate"), 1L)
})

test_that(".isolate_feat_names: follows the .qp_prefix naming convention", {
  # one network, no clash with an attribute name -> bare
  expect_equal(.isolate_feat_names("friends", c("age", "score")), "isolate")
  # one network, but an attribute is literally called "isolate" -> prefixed
  expect_equal(.isolate_feat_names("friends", c("age", "isolate")), "friends_isolate")
  # several networks -> always prefixed
  expect_equal(.isolate_feat_names(c("friends", "advice"), "age"),
               c("friends_isolate", "advice_isolate"))
  expect_equal(.isolate_feat_names(character(0), "age"), character(0))
})

test_that("netquickpred: the isolate flag bypasses the mincor screen", {
  fx <- make_iso_fixture()
  # mincor 0.35 keeps at least one alter feature (so the flag has something to
  # offset) while sitting far above the flag's own relevance of ~0.15: only
  # the screen bypass can keep it
  qp <- netquickpred(fx$attrs, list(friends = fx$net), mincor = 0.35)
  e <- qp$predictors$age
  expect_true("isolate" %in% e$net_features)
  expect_lt(e$relevance[["isolate"]], 0.35)
  others <- setdiff(names(e$relevance), "isolate")
  expect_true(all(e$relevance[others] >= 0.35))
})

test_that("netquickpred: the flag is dropped when no alter feature needs offsetting", {
  fx <- make_iso_fixture()
  # at mincor 0.5 nothing survives the screen, so no alter measure was
  # mean-filled and the indicator has no job to do
  qp <- netquickpred(fx$attrs, list(friends = fx$net), mincor = 0.5)
  e <- qp$predictors$age
  expect_length(e$net_features, 0)
  expect_false("isolate" %in% e$net_features)
})

test_that("netquickpred: a constant isolate flag is NOT forced in", {
  set.seed(5)
  n <- 40
  fr <- matrix(rbinom(n * n, 1, 0.25), n, n)   # dense: no isolates
  diag(fr) <- 0
  off <- which(row(fr) != col(fr))
  fr[sample(off, 40)] <- NA
  attrs <- data.frame(age = rnorm(n, 35, 8), score = rnorm(n))
  attrs$age[sample(n, 6)] <- NA
  expect_equal(sum(rowSums(fr, na.rm = TRUE) + colSums(fr, na.rm = TRUE) == 0), 0)
  qp <- netquickpred(attrs, list(friends = fr))
  expect_false("isolate" %in% qp$predictors$age$net_features)
})

test_that("netquickpred: near-duplicate isolate flags collapse to one", {
  # The case that motivates this: people with no friends also have no advice,
  # help, popularity or liking ties. Five near-copies of one column are rank
  # 1 and would spend four parameters on nothing, so the flags - exempt from
  # being pruned by ordinary predictors - must still prune each other.
  set.seed(5)
  n <- 50
  core_iso <- 39:50
  mk <- function(p, extra, seed) {
    set.seed(seed)
    m <- matrix(rbinom(n * n, 1, p), n, n)
    diag(m) <- 0
    out <- c(core_iso, extra)
    m[out, ] <- 0; m[, out] <- 0
    off <- which(row(m) != col(m))
    m[sample(setdiff(off, which(row(m) %in% out | col(m) %in% out)), 30)] <- NA
    m
  }
  nets <- list(friends = mk(.12, integer(0), 1), advice = mk(.10, 37, 2),
               help = mk(.11, 38, 3), popularity = mk(.09, integer(0), 4),
               liking = mk(.13, 36, 5))
  attrs <- data.frame(age = rnorm(n, 35, 8), score = rnorm(n),
                      grp = sample(c("a", "b"), n, TRUE))
  attrs$age[sample(n, 6)] <- NA
  attrs$score[sample(n, 5)] <- NA
  kept <- grep("isolate$", netquickpred(attrs, nets)$predictors$age$net_features,
               value = TRUE)
  expect_lt(length(kept), length(nets))
  expect_gte(length(kept), 1)
})

test_that("netquickpred: genuinely different isolate sets keep separate flags", {
  set.seed(11)
  n <- 45
  mk <- function(iso, seed) {
    set.seed(seed)
    m <- matrix(rbinom(n * n, 1, 0.12), n, n)
    diag(m) <- 0
    m[iso, ] <- 0; m[, iso] <- 0
    off <- which(row(m) != col(m))
    m[sample(setdiff(off, which(row(m) %in% iso | col(m) %in% iso)), 25)] <- NA
    m
  }
  # disjoint isolate sets: the two flags carry different information
  nets <- list(friends = mk(1:8, 1), advice = mk(30:40, 2))
  attrs <- data.frame(age = rnorm(n, 35, 8), score = rnorm(n))
  attrs$age[sample(n, 5)] <- NA
  kept <- grep("isolate$", netquickpred(attrs, nets)$predictors$age$net_features,
               value = TRUE)
  expect_setequal(kept, c("friends_isolate", "advice_isolate"))
})

test_that(".dedup_isolate_feats: drops correlated flags, keeps distinct ones", {
  df <- data.frame(a_isolate = c(1, 1, 1, 0, 0, 0, 0, 0),
                   b_isolate = c(1, 1, 1, 0, 0, 0, 0, 0),   # identical to a
                   c_isolate = c(0, 0, 0, 0, 1, 1, 1, 0))   # unrelated
  keep <- .dedup_isolate_feats(names(df), df)
  expect_length(keep, 2)
  expect_true("c_isolate" %in% keep)
  expect_equal(sum(c("a_isolate", "b_isolate") %in% keep), 1)
  # on a tie the earlier network's flag survives, so the result is stable
  expect_true("a_isolate" %in% keep)
  # relevance breaks the tie when it differs
  keep2 <- .dedup_isolate_feats(names(df), df,
                                relevance = c(a_isolate = 0.1, b_isolate = 0.9,
                                              c_isolate = 0.5))
  expect_true("b_isolate" %in% keep2)
  # nothing to do
  expect_equal(.dedup_isolate_feats("a_isolate", df), "a_isolate")
  expect_length(.dedup_isolate_feats(character(0), df), 0)
})

test_that(".clean_predictor_matrix: keep_raw shields the flag from PCA absorption", {
  set.seed(21)
  n <- 60
  x <- cbind(isolate = rep(c(0, 1), c(n - 5, 5)),
             matrix(rnorm(n * 12), n, 12,
                    dimnames = list(NULL, paste0("f", 1:12))))
  # max_cols forces the non-protected columns through prcomp
  out <- .clean_predictor_matrix(x, max_cols = 3, keep_raw = "isolate")
  expect_true("isolate" %in% colnames(out))
  # and it is the untouched original column, not a component
  expect_equal(unname(out[, "isolate"]), unname(x[, "isolate"]))
  expect_true(any(grepl("^PC", colnames(out))))

  # without the shield it is absorbed
  out2 <- .clean_predictor_matrix(x, max_cols = 3)
  expect_false("isolate" %in% colnames(out2))
})

test_that("netmice: an isolate's mean-filled alter measure is offset-able (flag reaches the model)", {
  fx <- make_iso_fixture()
  fit <- netmice(fx$attrs, list(friends = fx$net), m = 2, maxit = 2,
                 printFlag = FALSE, seed = 4)
  expect_s3_class(fit, "netmids")
  expect_false(anyNA(complete_netmice(fit, 1)$data))
})
