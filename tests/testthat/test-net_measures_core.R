nets <- fx_nets()
attrs <- fx_attrs()

test_that("net_measures_core: default args on a binary directed network", {
  out <- net_measures_core(nets$friends_bin, attrs[c("age", "status", "dept")])
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(attrs))
  expect_true(all(c("node_id", "outdegree", "indegree", "reciprocity_ratio",
                     "bonacich_power", "betweenness", "isolate", "constraint",
                     "harmonic_closeness", "local_clustering") %in% names(out)))
  expect_true(all(c("age_diff_out", "age_diff_in", "age_alter_mean", "age_alter_min",
                     "age_alter_max") %in% names(out)))
  expect_true(all(c("status_ei_index", "status_blau_index", "status_alter_mode") %in% names(out)))
  expect_true(all(c("dept_ei_index", "dept_blau_index", "dept_alter_mode") %in% names(out)))
  # node 1 was forced to isolate in the fixture
  expect_true(out$isolate[1])
  expect_equal(out$outdegree[1], 0)
  expect_equal(out$indegree[1], 0)
})

test_that("net_measures_core: attr_types override changes inference and is honoured", {
  expect_message(
    net_measures_core(nets$friends_bin, attrs[c("age", "status", "dept")]),
    "inferred attribute types"
  )
  out <- net_measures_core(nets$friends_bin, attrs[c("age", "status", "dept")],
                            attr_types = c(age = "continuous", status = "binary", dept = "multinomial"))
  expect_true(all(c("age_diff_out", "status_ei_index", "dept_ei_index") %in% names(out)))
  out2 <- net_measures_core(nets$friends_bin, attrs["dept"], attr_types = c(dept = "multinomial"))
  expect_true(all(c("dept_ei_index", "dept_blau_index", "dept_alter_mode") %in% names(out2)))
})

test_that("net_measures_core: id_col aligns shuffled attribute rows to named-vertex network", {
  m <- nets$friends_bin
  dimnames(m) <- list(attrs$id, attrs$id)
  shuffled <- attrs[sample(nrow(attrs)), ]
  out_aligned <- net_measures_core(m, shuffled[c("id", "age", "status", "dept")], id_col = "id")
  out_ordered <- net_measures_core(m, attrs[c("age", "status", "dept")])
  expect_equal(out_aligned$age_alter_mean, out_ordered$age_alter_mean)
  expect_equal(out_aligned$outdegree, out_ordered$outdegree)
})

test_that("net_measures_core: id_col errors clearly on an unnamed matrix", {
  expect_error(
    net_measures_core(nets$friends_bin, attrs[c("id", "age")], id_col = "id"),
    "no vertex names"
  )
})

test_that("net_measures_core: accepts an igraph object as well as a matrix", {
  g <- igraph::graph_from_adjacency_matrix(nets$friends_bin, mode = "directed", diag = FALSE)
  out_g <- net_measures_core(g, attrs[c("age", "status", "dept")])
  out_m <- net_measures_core(nets$friends_bin, attrs[c("age", "status", "dept")])
  expect_equal(out_g$outdegree, out_m$outdegree)
  expect_equal(out_g$constraint, out_m$constraint)
})

test_that("net_measures_core: undirected network yields NA reciprocity_ratio", {
  out <- net_measures_core(nets$colleagues_bin, attrs[c("age", "status", "dept")])
  expect_true(all(is.na(out$reciprocity_ratio)))
})

test_that("net_measures_core: works on a non-negative weighted network", {
  out <- net_measures_core(nets$advice_weighted, attrs[c("age", "status", "dept")])
  expect_equal(nrow(out), nrow(attrs))
  expect_false(anyNA(out$outdegree))
  expect_true(all(is.finite(out$betweenness)))
  expect_true(all(is.finite(out$harmonic_closeness)))
  expect_true(all(is.finite(out$local_clustering)))
})

test_that("net_measures_core: rejects a signed (negative-weight) network", {
  expect_error(net_measures_core(nets$trust_signed, attrs[c("age", "status", "dept")]),
               "non-negative weighted")
})

test_that("net_measures_core: weighted reciprocity_ratio matches a hand-built reference", {
  # 3-node directed weighted network: 1<->2 asymmetric (4 vs 2), node 3 isolated.
  m <- matrix(0, 3, 3)
  m[1, 2] <- 4; m[2, 1] <- 2
  out <- net_measures_core(m, data.frame(x = c(1, 2, 3)))
  expect_equal(out$reciprocity_ratio[1], abs(4 - 2))
  expect_equal(out$reciprocity_ratio[2], abs(4 - 2))
  expect_true(is.na(out$reciprocity_ratio[3]))
})

test_that("net_measures_core: weighted local_clustering matches a hand-built reference", {
  # 4-node undirected weighted network. Reference computed independently via
  # brute-force nested loops over the Opsahl-Panzarasa geometric-mean formula,
  # rather than re-deriving decimals by hand.
  w <- matrix(0, 4, 4)
  w[1, 2] <- w[2, 1] <- 2
  w[1, 3] <- w[3, 1] <- 8
  w[2, 3] <- w[3, 2] <- 5   # closes the 1-2-3 triangle
  w[2, 4] <- w[4, 2] <- 3   # 4 only connects to 2; 1-4 and 3-4 absent

  ref_clustering <- function(w) {
    n <- nrow(w)
    a <- (w != 0) * 1
    vapply(seq_len(n), function(i) {
      nb <- which(a[i, ] == 1)
      if (length(nb) < 2) return(0)
      denom <- 0; numer <- 0
      pairs <- combn(nb, 2)
      for (p in seq_len(ncol(pairs))) {
        j <- pairs[1, p]; h <- pairs[2, p]
        gm <- sqrt(w[i, j] * w[i, h])
        denom <- denom + gm
        if (a[j, h] == 1) numer <- numer + gm
      }
      numer / denom
    }, numeric(1))
  }

  expected <- ref_clustering(w)
  out <- net_measures_core(w, data.frame(x = 1:4))
  expect_equal(out$local_clustering, expected, tolerance = 1e-8)
  # node 1 has exactly one neighbour pair (2,3) and it IS closed -> ratio 1
  expect_equal(out$local_clustering[1], 1)
})

test_that("net_measures_core: errors informatively on mismatched attribute row count", {
  expect_error(net_measures_core(nets$friends_bin, attrs[1:5, c("age", "status")]),
               "rows but the network has")
})
