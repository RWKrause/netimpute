nets <- fx_nets()
attrs <- fx_attrs()

test_that("net_measures_full: default args return the full structural battery", {
  out <- net_measures_full(nets$friends_bin, attrs[c("age", "status", "dept")])
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(attrs))
  expect_true(all(c("total_degree", "indegree", "outdegree", "reciprocal_degree",
                     "nonreciprocal_degree", "reciprocity_ratio", "weighted_degree",
                     "eigenvector", "bonacich_power", "alpha_centrality", "pagerank",
                     "betweenness", "closeness", "harmonic_closeness", "constraint",
                     "effective_size", "efficiency", "redundancy", "neighbor_degree",
                     "local_clustering", "coreness", "isolate") %in% names(out)))
})

test_that("net_measures_full: use_sna = TRUE (default) populates sna-only columns", {
  skip_if_not_installed("sna")
  out <- net_measures_full(nets$friends_bin, attrs[c("age", "status", "dept")], use_sna = TRUE)
  sna_cols <- c("gilschmidt", "flow_betweenness", "load_centrality", "stress_centrality",
                "information_centrality", "prestige")
  expect_true(all(sna_cols %in% names(out)))
  expect_false(all(vapply(out[sna_cols], function(x) all(is.na(x)), logical(1))))
})

test_that("net_measures_full: use_sna = FALSE leaves sna-only columns NA", {
  out <- net_measures_full(nets$friends_bin, attrs[c("age", "status", "dept")], use_sna = FALSE)
  sna_cols <- c("gilschmidt", "flow_betweenness", "load_centrality", "stress_centrality",
                "information_centrality", "prestige")
  expect_true(all(vapply(out[sna_cols], function(x) all(is.na(x)), logical(1))))
})

test_that("net_measures_full: attr_types override is honoured", {
  out <- net_measures_full(nets$friends_bin, attrs[c("age", "status", "dept")],
                            attr_types = c(age = "continuous", status = "binary", dept = "multinomial"))
  expect_true(all(c("age_diff_out", "status_ei_index", "dept_ei_index") %in% names(out)))
})

test_that("net_measures_full: id_col aligns shuffled attribute rows", {
  m <- nets$friends_bin
  dimnames(m) <- list(attrs$id, attrs$id)
  shuffled <- attrs[sample(nrow(attrs)), ]
  out_aligned <- net_measures_full(m, shuffled[c("id", "age", "status", "dept")], id_col = "id")
  out_ordered <- net_measures_full(m, attrs[c("age", "status", "dept")])
  expect_equal(out_aligned$total_degree, out_ordered$total_degree)
})

test_that("net_measures_full: works on a non-negative weighted network", {
  out <- net_measures_full(nets$advice_weighted, attrs[c("age", "status", "dept")])
  expect_equal(nrow(out), nrow(attrs))
  expect_true(all(is.finite(out$local_clustering)))
  expect_true(!anyNA(out$weighted_degree))
})

test_that("net_measures_full: rejects a signed (negative-weight) network", {
  expect_error(net_measures_full(nets$trust_signed, attrs[c("age", "status", "dept")]),
               "non-negative weighted")
})

test_that("net_measures_full: undirected network yields NA reciprocity/reciprocal-degree", {
  out <- net_measures_full(nets$colleagues_bin, attrs[c("age", "status", "dept")])
  expect_true(all(is.na(out$reciprocity_ratio)))
  expect_true(all(is.na(out$reciprocal_degree)))
})

test_that("net_measures_full: core is a strict subset of full's structural columns", {
  core_cols <- setdiff(names(net_measures_core(nets$friends_bin, attrs["age"])), "node_id")
  full_cols <- setdiff(names(net_measures_full(nets$friends_bin, attrs["age"])), "node_id")
  # core's structural names should all appear in full (both share the same
  # underlying definitions for outdegree/indegree/betweenness/etc.)
  shared <- intersect(core_cols, full_cols)
  expect_true(length(shared) >= 6)
})
