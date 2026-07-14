nets <- fx_nets()
attrs <- fx_attrs()

test_that("net_measures: 'core' reproduces net_measures_core exactly", {
  expect_equal(
    suppressMessages(net_measures(nets$friends_bin, attrs[c("age", "status", "dept")],
                                  measure_set = "core")),
    suppressMessages(net_measures_core(nets$friends_bin, attrs[c("age", "status", "dept")]))
  )
})

test_that("net_measures: 'full' reproduces net_measures_full exactly", {
  expect_equal(
    suppressMessages(net_measures(nets$friends_bin, attrs[c("age", "status", "dept")],
                                  measure_set = "full")),
    suppressMessages(net_measures_full(nets$friends_bin, attrs[c("age", "status", "dept")]))
  )
})

test_that("net_measures: 'core' plus individual measures adds them to all core columns", {
  core <- suppressMessages(
    net_measures(nets$friends_bin, attrs["age"], measure_set = "core"))
  out <- suppressMessages(
    net_measures(nets$friends_bin, attrs["age"],
                 measure_set = c("core", "total_degree", "pagerank")))
  expect_true(all(names(core) %in% names(out)))
  expect_true(all(c("total_degree", "pagerank") %in% names(out)))
  # the added columns match their 'full' values
  full <- suppressMessages(
    net_measures(nets$friends_bin, attrs["age"], measure_set = "full"))
  expect_equal(out$total_degree, full$total_degree)
  expect_equal(out$pagerank, full$pagerank)
  # and the shared core columns are unchanged by the addition
  expect_equal(out[names(core)], core)
})

test_that("net_measures: an explicit structural list returns exactly those columns", {
  out <- net_measures(nets$friends_bin, measure_set = c("indegree", "betweenness"))
  expect_identical(names(out), c("node_id", "indegree", "betweenness"))
})

test_that("net_measures: attributes are optional without 'homophily', required with it", {
  expect_silent(net_measures(nets$friends_bin, measure_set = "indegree"))
  expect_error(
    net_measures(nets$friends_bin, measure_set = "core"),
    "`attributes` must be supplied"
  )
})

test_that("net_measures: 'homophily' alone returns only the attribute block", {
  out <- suppressMessages(
    net_measures(nets$friends_bin, attrs["age"], measure_set = "homophily"))
  expect_identical(names(out)[1], "node_id")
  expect_true(all(grepl("^node_id$|^age_", names(out))))
  expect_true("age_alter_mean" %in% names(out))
})

test_that("net_measures: duplicate and overlapping requests are collapsed", {
  out <- suppressMessages(
    net_measures(nets$friends_bin, attrs["age"],
                 measure_set = c("core", "indegree", "core")))
  expect_false(anyDuplicated(names(out)) > 0)
})

test_that("net_measures: unknown measure names error informatively", {
  expect_error(
    net_measures(nets$friends_bin, attrs["age"],
                 measure_set = c("core", "not_a_measure")),
    "Unknown entries in `measure_set`: not_a_measure"
  )
  expect_error(net_measures(nets$friends_bin, attrs["age"], measure_set = character(0)),
               "non-empty character")
})

test_that("net_measures: sna-only selection respects use_sna = FALSE", {
  out <- net_measures(nets$friends_bin, measure_set = c("indegree", "prestige"),
                      use_sna = FALSE)
  expect_true(all(is.na(out$prestige)))
  expect_false(anyNA(out$indegree))
})

test_that("net_measures: selection works on a weighted network", {
  out <- suppressMessages(
    net_measures(nets$advice_weighted, attrs["age"],
                 measure_set = c("weighted_degree", "local_clustering", "homophily")))
  expect_true(all(c("weighted_degree", "local_clustering", "age_alter_mean") %in% names(out)))
  expect_true(all(is.finite(out$local_clustering)))
})

test_that(".resolve_measure_set: expansion is idempotent and order-preserving", {
  r1 <- netimpute:::.resolve_measure_set(c("core", "total_degree"))
  expect_identical(netimpute:::.resolve_measure_set(r1), r1)
  # core order preserved, addition appended
  expect_identical(r1[1:2], c("outdegree", "indegree"))
  expect_identical(r1[length(r1)], "total_degree")
  expect_true("homophily" %in% r1)
})
