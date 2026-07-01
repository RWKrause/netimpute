make_plot_fixture <- function(n = 20, seed = 601, m = 2, maxit = 3) {
  set.seed(seed)
  friends <- fx_bin_directed(n = n, seed = seed + 1)
  advice  <- fx_weighted(n = n, seed = seed + 2)
  attrs <- fx_attrs(n = n, seed = seed)[c("age", "dept")]
  attrs$age[sample(n, 3)] <- NA
  attrs$dept[sample(n, 2)] <- NA
  off <- which(row(friends) != col(friends))
  friends[sample(off, 10)] <- NA
  netmice(attrs, list(friends = friends, advice = advice), m = m, maxit = maxit,
          printFlag = FALSE, seed = seed)
}

test_that("netmice()'s default maxit is 20", {
  expect_equal(eval(formals(netmice)$maxit), 20)
})

test_that("plot.netmids: default call (both pages) draws without error", {
  fit <- make_plot_fixture()
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_no_error(plot(fit))
  expect_true(file.exists(path))
})

test_that("plot.netmids: vars/nets = character(0) skip a page each", {
  fit <- make_plot_fixture()
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_no_error(plot(fit, nets = character(0)))
  expect_no_error(plot(fit, vars = character(0)))
})

test_that("plot.netmids: nothing selected messages instead of erroring", {
  fit <- make_plot_fixture()
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_message(plot(fit, vars = character(0), nets = character(0)), "nothing to plot")
})

test_that("plot.netmids: errors clearly on an unknown variable or network name", {
  fit <- make_plot_fixture()
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_error(plot(fit, vars = "not_a_variable"), "Unknown")
  expect_error(plot(fit, nets = "not_a_network"), "Unknown")
})

test_that("plot.netmids: a single-chain (m = 1) fit does not error (dimension-dropping edge case)", {
  fit <- make_plot_fixture(m = 1)
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_no_error(plot(fit))
})
