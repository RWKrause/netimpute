# Direct unit tests for the internal helpers in R/graph-utils.R - these are
# not exported, but testthat files inside a package's own tests/testthat/
# have direct access to the package namespace, so they can be called by
# their bare (dot-prefixed) name here without :::.

test_that(".reject_negative: errors on any negative value, ignores NA, passes otherwise", {
  expect_silent(.reject_negative(c(0, 1, 2, NA)))
  expect_error(.reject_negative(c(0, -1, 2)), "non-negative weighted")
  expect_error(.reject_negative(c(NA, -0.5)), "non-negative weighted")
})

test_that(".is_weighted_mat: distinguishes binary from weighted matrices", {
  expect_false(.is_weighted_mat(matrix(c(0, 1, 1, 0), 2, 2)))
  expect_true(.is_weighted_mat(matrix(c(0, 2, 1, 0), 2, 2)))
  expect_false(.is_weighted_mat(matrix(0, 3, 3)))  # empty network: vacuously binary
})

test_that(".infer_attr_type: classifies continuous, binary, and multinomial correctly", {
  expect_equal(.infer_attr_type(rnorm(20)), "continuous")
  expect_equal(.infer_attr_type(c(1, 2, 1, 2, NA)), "binary")
  expect_equal(.infer_attr_type(c("a", "b", "a")), "binary")          # 2 unique levels
  expect_equal(.infer_attr_type(c("a", "b", "c", "a")), "multinomial")
  expect_equal(.infer_attr_type(c(1, 2, 3, 4, 5)), "continuous")
  expect_equal(.infer_attr_type(rep(NA, 5)), "binary")  # documented current behaviour:
  # zero observed unique values satisfies length(u) <= 2, so an all-NA column
  # is classified as "binary" by default rather than erroring
})

test_that(".resolve_attr_types: honours overrides, infers the rest, validates values", {
  attrs <- data.frame(a = rnorm(5), b = c("x", "y", "x", "y", "x"))
  resolved <- suppressMessages(.resolve_attr_types(attrs, c(a = "continuous")))
  expect_equal(unname(resolved["a"]), "continuous")
  expect_equal(unname(resolved["b"]), "binary")  # inferred

  expect_message(.resolve_attr_types(attrs, NULL), "inferred attribute types")
  expect_no_message(.resolve_attr_types(attrs, c(a = "continuous", b = "multinomial")))

  expect_error(.resolve_attr_types(attrs, c(a = "not_a_real_type")),
               "must be one of")
})

test_that(".align_attributes: matches by id_col and validates row counts", {
  m <- matrix(0, 3, 3)
  dimnames(m) <- list(c("x", "y", "z"), c("x", "y", "z"))
  g <- igraph::graph_from_adjacency_matrix(m, mode = "directed", diag = FALSE)
  attrs <- data.frame(id = c("z", "x", "y"), v = c(3, 1, 2))
  aligned <- .align_attributes(g, attrs, id_col = "id")
  expect_equal(aligned$v, c(1, 2, 3))  # reordered to match vertex order x, y, z

  expect_error(.align_attributes(g, data.frame(v = 1:2)), "rows but the network has")
})

test_that(".dyad_reciprocity: matches a hand-built binary reference", {
  # 1<->2 mutual, 1->3 one-way (no 3->1), 2 and 3 otherwise unconnected.
  m <- matrix(0, 3, 3)
  m[1, 2] <- 1; m[2, 1] <- 1
  m[1, 3] <- 1
  g <- igraph::graph_from_adjacency_matrix(m, mode = "directed", diag = FALSE)
  out <- .dyad_reciprocity(g)

  expect_equal(out$reciprocal_degree,    c(1, 1, 0))
  expect_equal(out$nonreciprocal_degree, c(1, 0, 1))
  expect_equal(out$reciprocity_ratio,    c(0.5, 1, 0))
})

test_that(".dyad_reciprocity: undirected networks return all-NA", {
  g <- igraph::graph_from_adjacency_matrix(matrix(c(0, 1, 1, 0), 2, 2), mode = "undirected")
  out <- .dyad_reciprocity(g)
  expect_true(all(is.na(out$reciprocity_ratio)))
})

test_that(".bonacich_power: an explicit exponent matches calling igraph directly", {
  m <- matrix(c(0, 1, 1, 1, 0, 0, 1, 0, 0), 3, 3)
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected")
  out <- .bonacich_power(g, exponent = 0.1)
  ref <- as.numeric(igraph::power_centrality(g, exponent = 0.1, rescale = FALSE))
  expect_equal(out, ref)
})

test_that(".bonacich_power: symmetric (star) positions get equal power", {
  # A star: node 1 is the hub, nodes 2-5 are equivalent leaves - by symmetry
  # they must all receive identical Bonacich power regardless of exponent.
  g <- igraph::make_star(5, mode = "undirected")
  out <- .bonacich_power(g)
  expect_equal(out[2], out[3])
  expect_equal(out[3], out[4])
  expect_equal(out[4], out[5])
  expect_gt(out[1], out[2])  # the hub should out-rank any leaf
})

test_that(".structural_holes: matches Burt's canonical redundant-vs-open example", {
  # Ego (node 1) has two alters (2, 3).
  # Open case: alters NOT connected -> no redundancy, effective size = degree.
  m_open <- matrix(0, 3, 3)
  m_open[1, 2] <- m_open[2, 1] <- 1
  m_open[1, 3] <- m_open[3, 1] <- 1
  g_open <- igraph::graph_from_adjacency_matrix(m_open, mode = "undirected")
  h_open <- .structural_holes(g_open)
  expect_equal(h_open$redundancy[1], 0)
  expect_equal(h_open$effective_size[1], 2)
  expect_equal(h_open$efficiency[1], 1)

  # Closed case: alters ARE connected -> fully redundant.
  m_closed <- m_open
  m_closed[2, 3] <- m_closed[3, 2] <- 1
  g_closed <- igraph::graph_from_adjacency_matrix(m_closed, mode = "undirected")
  h_closed <- .structural_holes(g_closed)
  expect_equal(h_closed$redundancy[1], 1)
  expect_equal(h_closed$effective_size[1], 1)
  expect_equal(h_closed$efficiency[1], 0.5)
})

test_that(".structural_holes: an isolate has effective size 0 and NA redundancy/efficiency", {
  m <- matrix(0, 3, 3)
  m[2, 3] <- m[3, 2] <- 1  # node 1 isolated
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected")
  h <- .structural_holes(g)
  expect_equal(h$effective_size[1], 0)
  expect_true(is.na(h$redundancy[1]))
  expect_true(is.na(h$efficiency[1]))
})

test_that(".neighbor_degree: matches a hand-built reference (average alter degree)", {
  # Node 1 connects to node 2 (degree 3: ties to 1, 4, 5) and node 3
  # (degree 1: tie only to 1) -> knn[1] = mean(3, 1) = 2.
  m <- matrix(0, 5, 5)
  m[1, 2] <- m[2, 1] <- 1
  m[1, 3] <- m[3, 1] <- 1
  m[2, 4] <- m[4, 2] <- 1
  m[2, 5] <- m[5, 2] <- 1
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected")
  out <- .neighbor_degree(g)
  expect_equal(out[1], mean(c(igraph::degree(g)[2], igraph::degree(g)[3])))
  expect_equal(out[1], 2)
})

test_that(".neighbor_degree: an isolate returns NA", {
  m <- matrix(0, 3, 3)
  m[2, 3] <- m[3, 2] <- 1
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected")
  out <- .neighbor_degree(g)
  expect_true(is.na(out[1]))
})

test_that(".weighted_local_clustering: reduces to the binary coefficient when all weights are 1", {
  m <- matrix(0, 5, 5)
  m[1, 2] <- m[2, 1] <- 1
  m[1, 3] <- m[3, 1] <- 1
  m[2, 3] <- m[3, 2] <- 1  # closes 1-2-3
  m[1, 4] <- m[4, 1] <- 1  # 4 is not connected to 2 or 3 -> open triplets for node 1
  # .weighted_local_clustering() always reads the "weight" edge attribute, so
  # the comparison graph must actually carry one (all 1s here) rather than
  # being a plain unweighted igraph object.
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected", weighted = TRUE)
  weighted_out <- .weighted_local_clustering(g)
  binary_out <- igraph::transitivity(g, type = "local", isolates = "zero")
  expect_equal(weighted_out, binary_out, tolerance = 1e-8)
})
