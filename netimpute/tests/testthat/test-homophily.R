# Direct unit tests for the internal helpers in R/homophily.R, called by
# their bare (dot-prefixed) name - see the note at the top of
# test-graph-utils.R about why this works for a package's own tests.

test_that(".ei_index / .blau_index / .alter_mode match a hand-built categorical example", {
  # Undirected star: node 1 (group A) connects to 2 (A), 3 (B), 4 (B), 5 (B).
  m <- matrix(0, 5, 5)
  m[1, 2:5] <- 1; m[2:5, 1] <- 1
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected")
  grp <- c("A", "A", "B", "B", "B")

  # Node 1: alters are A, B, B, B -> 1 internal (same), 3 external (different).
  expect_equal(.ei_index(g, grp)[1], (3 - 1) / 4)
  # p_A = 1/4, p_B = 3/4 -> 1 - (0.25^2 + 0.75^2)
  expect_equal(.blau_index(g, grp)[1], 1 - (0.25^2 + 0.75^2))
  expect_equal(.alter_mode(g, grp)[1], "B")

  # Node 2: its only alter is node 1 (A) -> 1 internal, 0 external.
  expect_equal(.ei_index(g, grp)[2], (0 - 1) / 1)
  expect_equal(.blau_index(g, grp)[2], 0)  # a single, homogeneous alter
  expect_equal(.alter_mode(g, grp)[2], "A")
})

test_that(".ei_index / .blau_index / .alter_mode return NA/NA_character_ for an isolate", {
  m <- matrix(0, 3, 3)
  m[2, 3] <- m[3, 2] <- 1
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected")
  grp <- c("A", "A", "B")
  expect_true(is.na(.ei_index(g, grp)[1]))
  expect_true(is.na(.blau_index(g, grp)[1]))
  expect_true(is.na(.alter_mode(g, grp)[1]))
})

test_that(".ei_index: fully homophilous vs. fully heterophilous bounds are -1 and 1", {
  # Triangle, all same category -> every alter is internal -> ei_index = -1.
  m <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected")
  expect_true(all(.ei_index(g, c("A", "A", "A")) == -1))
  # Same triangle, every node a different category -> every alter external -> +1.
  expect_true(all(.ei_index(g, c("A", "B", "C")) == 1))
})

test_that(".alter_continuous_stats matches a hand-built directed reference via the edge list", {
  # 1 <-> 2 (mutual), 1 -> 3 (one-way only). Confirmed separately that
  # igraph::neighbors(mode = "all") counts a mutual tie twice (once per
  # direction) - the reference below is built the same way, from the edge
  # list, rather than assumed.
  m <- matrix(0, 3, 3)
  m[1, 2] <- 1; m[2, 1] <- 1
  m[1, 3] <- 1
  g <- igraph::graph_from_adjacency_matrix(m, mode = "directed", diag = FALSE)
  x <- c(10, 14, 20)

  el <- igraph::as_edgelist(g, names = FALSE)
  ref <- lapply(seq_len(3), function(i) {
    out_alters <- el[el[, 1] == i, 2]
    in_alters  <- el[el[, 2] == i, 1]
    all_alters <- c(out_alters, in_alters)  # each incident directed edge once
    c(
      diff_out       = if (length(out_alters)) mean(abs(x[i] - x[out_alters])) else NA_real_,
      diff_in        = if (length(in_alters))  mean(abs(x[i] - x[in_alters]))  else NA_real_,
      alter_mean     = if (length(all_alters)) mean(x[all_alters]) else NA_real_,
      alter_mean_in  = if (length(in_alters))  mean(x[in_alters])  else NA_real_,
      alter_mean_out = if (length(out_alters)) mean(x[out_alters]) else NA_real_,
      alter_min      = if (length(all_alters)) min(x[all_alters])  else NA_real_,
      alter_max      = if (length(all_alters)) max(x[all_alters])  else NA_real_
    )
  })
  ref <- as.data.frame(do.call(rbind, ref))

  out <- .alter_continuous_stats(g, x)
  expect_equal(out, ref, tolerance = 1e-8, ignore_attr = TRUE)

  # And confirm node 1's values by direct hand computation, independent of
  # the edge-list reference above:
  # out-neighbours of 1: {2, 3} -> mean(|10-14|, |10-20|) = mean(4, 10) = 7;
  #                                alter mean = mean(14, 20) = 17
  # in-neighbours of 1: {2}     -> |10-14| = 4; alter mean = 14
  # all incident neighbours of 1 (edge-counted, node 2 appears twice since
  # it is both an out- and an in-neighbour): {2, 3, 2} -> x = 14, 20, 14 ->
  # mean 16, min 14, max 20
  expect_equal(out$diff_out[1], 7)
  expect_equal(out$diff_in[1], 4)
  expect_equal(out$alter_mean[1], 16)
  expect_equal(out$alter_mean_in[1], 14)
  expect_equal(out$alter_mean_out[1], 17)
  expect_equal(out$alter_min[1], 14)
  expect_equal(out$alter_max[1], 20)
})

test_that(".alter_continuous_stats: an isolate gives NA throughout", {
  m <- matrix(0, 3, 3)
  m[2, 3] <- 1
  g <- igraph::graph_from_adjacency_matrix(m, mode = "directed", diag = FALSE)
  out <- .alter_continuous_stats(g, c(1, 2, 3))
  expect_true(all(is.na(out[1, ])))
})

test_that(".homophily_block: dispatches continuous vs categorical columns correctly", {
  m <- matrix(0, 4, 4)
  m[1, 2] <- m[2, 1] <- 1
  m[1, 3] <- m[3, 1] <- 1
  g <- igraph::graph_from_adjacency_matrix(m, mode = "undirected")
  attrs <- data.frame(age = c(20, 30, 40, 50), grp = c("A", "A", "B", "B"))

  out <- suppressMessages(.homophily_block(g, attrs, attr_types = NULL))
  expect_true(all(c("age_diff_out", "age_diff_in", "age_alter_mean",
                     "age_alter_mean_in", "age_alter_mean_out",
                     "age_alter_min", "age_alter_max") %in% names(out)))
  expect_true(all(c("grp_ei_index", "grp_blau_index", "grp_alter_mode") %in% names(out)))
  expect_false(any(grepl("^age_ei_index|^grp_diff_out", names(out))))

  # attr_types can force a numeric-looking column to be treated as categorical
  out2 <- .homophily_block(g, attrs["age"], attr_types = c(age = "multinomial"))
  expect_true(all(c("age_ei_index", "age_blau_index", "age_alter_mode") %in% names(out2)))
})
