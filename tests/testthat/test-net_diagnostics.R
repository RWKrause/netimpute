test_that("net_diagnostics: returns all five finite diagnostics for a normal network", {
  m <- fx_bin_directed(n = 20)
  d <- net_diagnostics(m)
  expect_named(d, c("density", "reciprocity", "transitivity", "n_isolates", "avg_inv_geodesic"))
  expect_true(all(vapply(d, function(x) is.finite(x) || is.na(x), logical(1))))
})

test_that("net_diagnostics: isolates never produce NaN/Inf in avg_inv_geodesic", {
  m <- matrix(0, 10, 10)
  m[1, 2] <- m[2, 3] <- m[3, 1] <- 1  # a 3-cycle among nodes 1-3
  # nodes 4:10 are isolates, entirely disconnected from everything, including
  # each other - this is the scenario the "so isolates don't lead to NaN/Inf"
  # requirement is specifically about.
  d <- net_diagnostics(m)
  expect_true(is.finite(d$avg_inv_geodesic))
  expect_equal(d$n_isolates, 7)
})

test_that("net_diagnostics: undirected network has reciprocity exactly 1", {
  m <- fx_bin_undirected(n = 15)
  d <- net_diagnostics(m)
  expect_equal(d$reciprocity, 1)
})

test_that("net_diagnostics: an empty (no-tie) network doesn't error", {
  m <- matrix(0, 8, 8)
  d <- net_diagnostics(m)
  expect_equal(d$density, 0)
  expect_equal(d$n_isolates, 8)
  expect_equal(d$avg_inv_geodesic, 0)
})

test_that("net_diagnostics: works on a weighted matrix (binarized internally)", {
  m <- fx_weighted(n = 15)
  d1 <- net_diagnostics(m)
  d2 <- net_diagnostics((m != 0) * 1)
  expect_equal(d1, d2)
})

test_that("net_diagnostics: `directed` overrides the isSymmetric() inference", {
  # asymmetric matrix: forcing undirected changes the graph that is built
  m <- matrix(0, 5, 5)
  m[1, 2] <- 1; m[2, 3] <- 1; m[3, 1] <- 1
  auto <- net_diagnostics(m)
  forced_u <- net_diagnostics(m, directed = FALSE)
  forced_d <- net_diagnostics(m, directed = TRUE)
  expect_equal(auto, forced_d)                    # inference agrees: directed
  expect_false(isTRUE(all.equal(forced_u, forced_d)))
  expect_equal(forced_d$reciprocity, 0)           # a 3-cycle: no mutual ties
  expect_equal(forced_u$reciprocity, 1)           # undirected: 1 by convention

  # forcing undirected must not emit igraph's asymmetric-matrix warning
  expect_silent(net_diagnostics(m, directed = FALSE))
})

test_that("net_diagnostics: avg_inv_geodesic follows tie direction when directed", {
  # 1 -> 2 -> 3 -> 1 plus two isolates. Directed: each of the 3 cycle nodes
  # reaches the other two at distances 1 and 2 -> sum(1/d) = 3 * 1.5 = 4.5
  # over 20 ordered pairs = 0.225. Ignoring direction would give 0.3.
  m <- matrix(0, 5, 5)
  m[1, 2] <- 1; m[2, 3] <- 1; m[3, 1] <- 1
  expect_equal(net_diagnostics(m, directed = TRUE)$avg_inv_geodesic, 0.225)
  expect_equal(net_diagnostics(m, directed = FALSE)$avg_inv_geodesic, 0.3)
})

test_that(".imp_tie_stats: mean/variance over the imputed cells only", {
  mat <- matrix(0, 4, 4)
  ry <- matrix(TRUE, 4, 4)
  diag(ry) <- FALSE               # diagonal is never an observation
  ry[1, 2] <- ry[2, 3] <- ry[3, 4] <- FALSE   # three imputed cells
  mat[1, 2] <- 1; mat[2, 3] <- 0; mat[3, 4] <- 1
  mat[4, 1] <- 1                  # an OBSERVED tie: must not be counted
  st <- .imp_tie_stats(mat, ry, undirected = FALSE)
  expect_equal(st$mean, 2 / 3)
  expect_equal(st$var, stats::var(c(1, 0, 1)))
})

test_that(".imp_tie_stats: an undirected dyad counts once, not twice", {
  mat <- matrix(0, 4, 4)
  ry <- matrix(TRUE, 4, 4); diag(ry) <- FALSE
  # one symmetric imputed dyad set to 1, and a second set to 0
  ry[1, 2] <- ry[2, 1] <- FALSE; mat[1, 2] <- mat[2, 1] <- 1
  ry[3, 4] <- ry[4, 3] <- FALSE; mat[3, 4] <- mat[4, 3] <- 0
  st <- .imp_tie_stats(mat, ry, undirected = TRUE)
  expect_equal(st$mean, 0.5)
  expect_equal(st$var, stats::var(c(1, 0)))   # n = 2, not n = 4
})

test_that(".imp_tie_stats: degenerate cell counts give NA, not an error", {
  mat <- matrix(0, 3, 3)
  ry <- matrix(TRUE, 3, 3); diag(ry) <- FALSE
  # nothing imputed at all
  expect_equal(.imp_tie_stats(mat, ry, FALSE), list(mean = NA_real_, var = NA_real_))
  # exactly one imputed cell: a mean but no variance
  ry[1, 2] <- FALSE; mat[1, 2] <- 1
  st <- .imp_tie_stats(mat, ry, FALSE)
  expect_equal(st$mean, 1)
  expect_true(is.na(st$var))
})

test_that(".imp_tie_stats: weighted ties use the raw values", {
  mat <- matrix(0, 3, 3)
  ry <- matrix(TRUE, 3, 3); diag(ry) <- FALSE
  ry[1, 2] <- ry[1, 3] <- FALSE
  mat[1, 2] <- 4; mat[1, 3] <- 2
  st <- .imp_tie_stats(mat, ry, FALSE)
  expect_equal(st$mean, 3)
  expect_equal(st$var, 2)
})
