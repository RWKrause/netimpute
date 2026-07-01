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
