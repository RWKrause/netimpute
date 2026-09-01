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

test_that("plot.netmids: restores the caller's par() settings", {
  fit <- make_plot_fixture()
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  before <- graphics::par(no.readonly = TRUE)
  invisible(plot(fit))                    # draws BOTH pages
  after <- graphics::par(no.readonly = TRUE)
  for (p in c("mfrow", "mar", "mgp", "tcl", "cex.main", "cex.lab", "cex.axis")) {
    expect_equal(after[[p]], before[[p]], info = p)
  }
})

test_that("plot.netmids: maxit = 1 plots chains as chains, not as iterations", {
  fit <- make_plot_fixture(m = 3, maxit = 1)
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  # `chainMean[v, , ]` drops to a length-m vector here; matplot() must still
  # see a 1 x m matrix (one point per chain), never an m-point single line.
  seen <- NULL
  local_mocked_bindings(
    matplot = function(x, ...) { seen <<- c(seen, list(dim(x))); invisible(NULL) },
    legend = function(...) invisible(NULL),   # nothing was really drawn to key
    .package = "graphics"
  )
  invisible(plot(fit, nets = character(0)))
  expect_true(length(seen) > 0)
  expect_true(all(vapply(seen, function(d) identical(d, c(1L, 3L)), logical(1))))
})

test_that("plot.netmids: `...` overrides the method's own matplot defaults", {
  fit <- make_plot_fixture()
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_no_error(plot(fit, col = "red"))
  expect_no_error(plot(fit, main = "custom"))
  expect_no_error(plot(fit, lty = 2, type = "b"))
})

test_that("netmice(): netImpMean/netImpVar cover the networks with missing ties", {
  fit <- make_plot_fixture(m = 2, maxit = 3)
  # `advice` in the fixture is fully observed; only `friends` has missing ties
  expect_equal(fit$net_missing, "friends")
  for (a in list(fit$netImpMean, fit$netImpVar)) {
    expect_equal(dim(a), c(1L, 3L, 2L))
    expect_equal(dimnames(a)[[1]], "friends")
    expect_false(anyNA(a))
  }
  # a binary network: the mean is the density among the imputed dyads
  expect_true(all(fit$netImpMean >= 0 & fit$netImpMean <= 1))
})

test_that("netmice(): no missing ties gives zero-row imputed arrays that still plot", {
  set.seed(9)
  n <- 20
  friends <- fx_bin_directed(n = n, seed = 91)   # complete, no NAs
  attrs <- fx_attrs(n = n, seed = 9)["age"]
  attrs$age[sample(n, 3)] <- NA
  fit <- netmice(attrs, list(friends = friends), m = 2, maxit = 2,
                 printFlag = FALSE, seed = 9)
  expect_length(fit$net_missing, 0)
  expect_equal(dim(fit$netImpMean), c(0L, 2L, 2L))
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_no_error(plot(fit))     # imputed-tie page simply has nothing to draw
})

test_that("plot.netmids: the imputed-tie page is drawn, and skipped with nets = character(0)", {
  fit <- make_plot_fixture(m = 2, maxit = 3)
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  titles <- NULL
  grab <- function(x, ...) {
    d <- list(...)
    titles <<- c(titles, if (is.null(d$main)) NA_character_ else d$main)
    # still advance the device: a mock that draws nothing leaves par(mfrow)
    # mid-page, and restoring it warns "calling par(new=TRUE) with no plot"
    graphics::plot.new()
    invisible(NULL)
  }
  local_mocked_bindings(matplot = grab, legend = function(...) invisible(NULL),
                        .package = "graphics")
  invisible(plot(fit))
  expect_true(any(grepl("friends - imputed ties, mean", titles, fixed = TRUE)))
  expect_true(any(grepl("friends - imputed ties, variance", titles, fixed = TRUE)))
  # `advice` has no missing ties, so it gets network panels but no imputed page
  expect_true(any(grepl("advice - density", titles, fixed = TRUE)))
  expect_false(any(grepl("advice - imputed", titles, fixed = TRUE)))

  titles <- NULL
  invisible(plot(fit, nets = character(0)))
  expect_false(any(grepl("imputed", titles, fixed = TRUE)))
})

test_that("plot.netmids: an attribute with exactly one missing value does not error", {
  # its chainVar row is all NA (var() of one value), which matplot() rejects
  # with "need finite 'ylim' values" - the panel must be drawn empty instead
  set.seed(31)
  n <- 20
  friends <- fx_bin_directed(n = n, seed = 311)
  off <- which(row(friends) != col(friends))
  friends[sample(off, 10)] <- NA
  attrs <- fx_attrs(n = n, seed = 31)["age"]
  attrs$age[3] <- NA                       # EXACTLY one
  fit <- netmice(attrs, list(friends = friends), m = 2, maxit = 2,
                 printFlag = FALSE, seed = 31)
  expect_true(all(is.na(fit$chainVar["age", , ])))
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 800, height = 600)
  on.exit({ grDevices::dev.off(); unlink(path) }, add = TRUE)
  expect_no_error(plot(fit))
})

test_that("imputed-cell traces are more sensitive than the whole-network ones", {
  # the reason the imputed page exists: netChain mixes in the observed ties,
  # so it moves far less than the cells the sampler actually draws
  set.seed(4)
  n <- 30
  friends <- matrix(rbinom(n * n, 1, 0.15), n, n)
  diag(friends) <- 0
  off <- which(row(friends) != col(friends))
  friends[sample(off, 60)] <- NA            # 60 of 870 dyads = ~7%
  attrs <- data.frame(age = rnorm(n, 35, 8))
  attrs$age[sample(n, 5)] <- NA
  fit <- netmice(attrs, list(friends = friends), m = 3, maxit = 10,
                 printFlag = FALSE, seed = 2)
  whole <- diff(range(fit$netChain["friends", "density", , ]))
  imputed <- diff(range(fit$netImpMean["friends", , ]))
  expect_gt(imputed, 5 * whole)
})
