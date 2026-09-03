# Edgelist input. The conversion runs before the sampler and must produce
# exactly the matrices the user would have built by hand.

el_data <- function() {
  data.frame(pid = c("a", "b", "c", "d"), age = c(20, 30, 40, 50),
             stringsAsFactors = FALSE)
}
el_simple <- function() {
  data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
             stringsAsFactors = FALSE)
}
# collect warnings without letting them abort, so both "warns" and "stays
# silent" can be asserted on the same call shape
warns_of <- function(expr) {
  ws <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  ws
}

# ---- column resolution -----------------------------------------------------

test_that(".el_resolve_names: announces the assumed columns and honours overrides", {
  el <- data.frame(s = "a", r = "b", w = 1)
  expect_message(
    got <- .el_resolve_names(el, NULL, character(0), "long", TRUE),
    "sender = 's', receiver = 'r', value = 'w'")
  expect_equal(got, list(sender = "s", receiver = "r", value = "w"))

  # two columns only -> unweighted
  expect_message(
    got2 <- .el_resolve_names(data.frame(s = "a", r = "b"), NULL,
                              character(0), "long", TRUE),
    "unweighted")
  expect_null(got2$value)

  # explicit names win, and are not announced
  expect_silent(
    got3 <- .el_resolve_names(el, c(sender = "r", receiver = "s"),
                              character(0), "long", TRUE))
  expect_equal(got3$sender, "r")
})

test_that(".el_resolve_names: rejects ambiguous or malformed specifications", {
  el4 <- data.frame(a = 1, b = 2, c = 3, d = 4)
  expect_error(.el_resolve_names(el4, NULL, character(0), "long", FALSE),
               "Cannot tell which column holds the tie value")
  # ...unless the extras are split columns
  expect_silent(.el_resolve_names(el4, NULL, c("c", "d"), "long", FALSE))

  expect_error(.el_resolve_names(el4, c("a", "b"), character(0), "long", FALSE),
               "NAMED vector")
  expect_error(.el_resolve_names(el4, c(sender = "a"), character(0), "long", FALSE),
               "must contain both")
  expect_error(.el_resolve_names(el4, c(sender = "a", receiver = "zz"),
                                 character(0), "long", FALSE),
               "not in")
  expect_error(.el_resolve_names(el4, c(sender = "a", receiver = "b", nope = "c"),
                                 character(0), "long", FALSE),
               "only name")
  # wide takes its values from the split columns, so a value column is ambiguous
  expect_error(.el_resolve_names(data.frame(s = "a", r = "b", v = 1, w = 2),
                                 c(sender = "s", receiver = "r", value = "v"),
                                 "w", "wide", FALSE),
               "ambiguous")
})

# ---- splitting -------------------------------------------------------------

test_that(".el_network_keys: long splits on the interaction of the split columns", {
  el <- data.frame(f = letters[1:4], t = letters[2:5],
                   wave = c(1, 1, 2, 2), question = c("x", "y", "x", "x"))
  k <- .el_network_keys(el, "long", c("wave", "question"))
  expect_setequal(k$names, c("wave_1_question_x", "wave_1_question_y",
                             "wave_2_question_x"))
  expect_equal(k$rows[["wave_2_question_x"]], c(3L, 4L))
  expect_equal(k$components[["wave_1_question_y"]], c("wave_1", "question_y"))

  # no split at all -> a single network
  k0 <- .el_network_keys(el, "long", character(0))
  expect_equal(k0$names, "net1")
})

test_that(".el_network_keys: wide makes one network per split column", {
  el <- data.frame(f = "a", t = "b", friends = 1, advice = 0)
  k <- .el_network_keys(el, "wide", c("friends", "advice"))
  expect_equal(k$names, c("friends", "advice"))
  expect_error(.el_network_keys(el, "wide", character(0)), "needs `edgelist_split`")
  expect_error(.el_network_keys(el, "wide", c("friends", "f")), "must be numeric")
  expect_error(.el_network_keys(el, "long", "nope"), "not in the edgelist")
})

# ---- matrix construction ---------------------------------------------------

test_that("edgelist -> matrices reproduces a hand-built adjacency matrix", {
  dat <- el_data()
  r <- .networks_from_edgelist(el_simple(), list(), dat, id = "pid",
                               printFlag = FALSE)
  want <- matrix(0, 4, 4, dimnames = list(dat$pid, dat$pid))
  want["a", "b"] <- 1; want["b", "c"] <- 1; want["c", "a"] <- 1
  expect_equal(r$mats$net1, want)
  expect_equal(nrow(r$data), 4)
})

test_that("edgelist -> matrices carries tie values when a value column exists", {
  dat <- el_data()
  el <- data.frame(from = c("a", "b"), to = c("b", "c"), w = c(2.5, 4))
  r <- .networks_from_edgelist(el, list(), dat, id = "pid", printFlag = FALSE)
  expect_equal(r$mats$net1["a", "b"], 2.5)
  expect_equal(r$mats$net1["b", "c"], 4)
  expect_equal(r$mats$net1["a", "c"], 0)
})

test_that("wide format: NA in a network's column marks that cell unknown", {
  dat <- el_data()
  el <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
                   friends = c(1, 0, NA), advice = c(NA, 1, 1))
  r <- .networks_from_edgelist(el, list(edgelist_format = "wide",
                                        edgelist_split = c("friends", "advice")),
                               dat, id = "pid", printFlag = FALSE)
  expect_true(is.na(r$mats$friends["c", "a"]))
  expect_equal(r$mats$friends["a", "b"], 1)
  expect_equal(r$mats$friends["b", "c"], 0)   # an observed non-tie, not missing
  expect_true(is.na(r$mats$advice["a", "b"]))
  expect_equal(r$mats$advice["c", "a"], 1)
})

test_that("self-loops are dropped and duplicate pairs are an error", {
  dat <- el_data()
  expect_message(
    r <- .networks_from_edgelist(data.frame(from = c("a", "a"), to = c("a", "b")),
                                 list(), dat, id = "pid", printFlag = TRUE),
    "self-loop")
  expect_equal(r$mats$net1["a", "a"], 0)
  expect_equal(r$mats$net1["a", "b"], 1)

  expect_error(
    .networks_from_edgelist(data.frame(from = c("a", "a"), to = c("b", "b")),
                            list(), dat, id = "pid", printFlag = FALSE),
    "duplicate sender-receiver pair")
})

test_that("an undirected declaration mirrors each listed edge", {
  dat <- el_data()
  r <- .networks_from_edgelist(data.frame(from = "a", to = "b"),
                               list(directed = "undirected"), dat,
                               id = "pid", printFlag = FALSE)
  expect_equal(r$mats$net1["a", "b"], 1)
  expect_equal(r$mats$net1["b", "a"], 1)
  expect_true(isSymmetric(unname(r$mats$net1)))
  expect_false(r$directed[["net1"]])

  # every accepted spelling
  for (w in c("undirected", "graph")) {
    expect_false(.el_resolve_directed(w, list(names = "n", components = list(n = "n")))[["n"]])
  }
  for (w in c("directed", "digraph")) {
    expect_true(.el_resolve_directed(w, list(names = "n", components = list(n = "n")))[["n"]])
  }
  expect_error(.el_resolve_directed("sideways",
                                    list(names = "n", components = list(n = "n"))),
               "must be one of")
})

# ---- nodelist --------------------------------------------------------------

test_that("nodelist: per-network rosters mark off-roster nodes structural", {
  dat <- el_data()
  el <- data.frame(from = c("a", "b", "a"), to = c("b", "c", "c"),
                   wave = c(1, 1, 2))
  r <- .networks_from_edgelist(
    el, list(edgelist_split = "wave",
             nodelist = list(wave_1 = c("a", "b", "c", "d"),
                             wave_2 = c("a", "c"))),
    dat, id = "pid", printFlag = FALSE)
  # b and d are off wave_2's roster. The diagonal is excluded throughout: a
  # self-loop is not a tie, so it is never "structurally absent" either.
  expect_true(all(r$structural$wave_2["b", -2]))
  expect_true(all(r$structural$wave_2[-4, "d"]))
  expect_false(r$structural$wave_2["a", "c"])
  expect_false(any(diag(r$structural$wave_2)))
  # wave_1's roster is everyone, so it gets no structural cells
  expect_null(r$structural$wave_1)
})

test_that("nodelist: several split entries are intersected", {
  dat <- el_data()
  el <- data.frame(from = "a", to = "b", wave = 1, question = "x")
  r <- .networks_from_edgelist(
    el, list(edgelist_split = c("wave", "question"),
             nodelist = list(wave_1 = c("a", "b", "c"),
                             question_x = c("a", "b", "d"))),
    dat, id = "pid", printFlag = FALSE)
  # intersection is {a, b}; c and d are off-roster
  s <- r$structural$wave_1_question_x
  expect_false(s["a", "b"])
  expect_true(all(s["c", -3]))
  expect_true(all(s["d", -4]))
  expect_false(any(diag(s)))
})

test_that("nodelist: a tie reaching off the roster is an informative error", {
  dat <- el_data()
  expect_error(
    .networks_from_edgelist(el_simple(), list(nodelist = c("a", "b")), dat,
                            id = "pid", printFlag = FALSE),
    "not on its node roster")
})

test_that("nodelist: accepts a data column, a vector, and reports gaps", {
  dat <- el_data()
  r1 <- .networks_from_edgelist(el_simple(), list(nodelist = "pid"), dat,
                                id = "pid", printFlag = FALSE)
  expect_null(r1$structural)
  expect_error(
    .networks_from_edgelist(cbind(el_simple(), wave = 1),
                            list(edgelist_split = "wave",
                                 nodelist = list(other = "a")),
                            dat, id = "pid", printFlag = FALSE),
    "does not cover network")
  expect_error(
    .networks_from_edgelist(el_simple(), list(nodelist = list(c("a", "b"))),
                            dat, id = "pid", printFlag = FALSE),
    "must be named")
})

# ---- missing ---------------------------------------------------------------

test_that("missing: a bare vector means outgoing ties", {
  dat <- el_data()
  r <- .networks_from_edgelist(el_simple(), list(missing = "d"), dat,
                               id = "pid", printFlag = FALSE)
  expect_true(all(is.na(r$mats$net1["d", -4])))
  expect_false(anyNA(r$mats$net1[, "d"]))   # incoming ties untouched
})

test_that("missing: out/in are applied to rows and columns respectively", {
  dat <- el_data()
  # c both receives and is marked missing$in, so the direction-matched warning
  # fires by design; it has its own test below
  r <- suppressWarnings(.networks_from_edgelist(
    el_simple(), list(missing = list(out = "d", `in` = "c")),
    dat, id = "pid", printFlag = FALSE))
  expect_true(all(is.na(r$mats$net1["d", -4])))
  # column c is NA except the listed b -> c tie
  expect_equal(r$mats$net1["b", "c"], 1)
  expect_true(is.na(r$mats$net1["a", "c"]))
})

test_that("missing: listed edges stay observed inside a missing row", {
  dat <- el_data()
  r <- suppressWarnings(.networks_from_edgelist(
    el_simple(), list(missing = "a"), dat, id = "pid", printFlag = FALSE))
  expect_equal(r$mats$net1["a", "b"], 1)    # listed: observed
  expect_true(is.na(r$mats$net1["a", "c"])) # everything else: unknown
})

test_that("missing: the conflict warning is strictly direction-matched", {
  dat <- el_data()
  one <- data.frame(from = "a", to = "b", stringsAsFactors = FALSE)
  # a SENDS and is marked missing$out -> warn
  expect_match(
    warns_of(.networks_from_edgelist(one, list(missing = list(out = "a")), dat,
                                     id = "pid", printFlag = FALSE)),
    "marked missing in that same direction", all = FALSE)
  # b only RECEIVES but is marked missing$out -> silent. This is the ordinary
  # actor-non-response case: someone who did not report is still nominated by
  # others, so warning here would fire on essentially every real dataset.
  expect_length(
    warns_of(.networks_from_edgelist(one, list(missing = list(out = "b")), dat,
                                     id = "pid", printFlag = FALSE)), 0)
  # b RECEIVES and is marked missing$in -> warn
  expect_match(
    warns_of(.networks_from_edgelist(one, list(missing = list(`in` = "b")), dat,
                                     id = "pid", printFlag = FALSE)),
    "marked missing in that same direction", all = FALSE)
  # a only SENDS but is marked missing$in -> silent
  expect_length(
    warns_of(.networks_from_edgelist(one, list(missing = list(`in` = "a")), dat,
                                     id = "pid", printFlag = FALSE)), 0)
})

# ---- new nodes -------------------------------------------------------------

test_that("nodes absent from `data` are added, categorised, and warned about", {
  dat <- el_data()
  el <- data.frame(from = c("a", "e"), to = c("e", "b"), stringsAsFactors = FALSE)
  ws <- warns_of(
    r <- .networks_from_edgelist(el,
                                 list(nodelist = c(letters[1:6]), missing = "g"),
                                 dat, id = "pid", printFlag = FALSE))
  expect_length(ws, 1)
  expect_match(ws, "added with connections \\(1\\)")
  expect_match(ws, "added as isolates \\(1\\)")
  expect_match(ws, "added entirely missing \\(1\\)")
  expect_setequal(r$data$pid, c("a", "b", "c", "d", "e", "f", "g"))
  # new rows carry no attribute information
  expect_true(all(is.na(r$data$age[r$data$pid %in% c("e", "f", "g")])))
  expect_equal(dim(r$mats$net1), c(7L, 7L))
})

# ---- id --------------------------------------------------------------------

test_that("an edgelist without `id` is refused", {
  expect_error(
    .networks_from_edgelist(el_simple(), list(), el_data(), id = NULL),
    "`id` must name the column")
  expect_error(
    .networks_from_edgelist(el_simple(), list(), el_data(), id = "nope"),
    "not found in `data`")
  expect_error(
    .networks_from_edgelist(el_simple(), list(),
                            data.frame(pid = c("a", "a"), age = 1:2), id = "pid"),
    "duplicate values")
})

test_that("netmice: `id` reorders matrices to the rows of `data` and drops the column", {
  set.seed(6)
  n <- 12
  ids <- paste0("p", seq_len(n))
  m <- matrix(rbinom(n * n, 1, 0.3), n, n, dimnames = list(ids, ids))
  diag(m) <- 0
  m[1, 3] <- NA
  attrs <- data.frame(pid = ids, age = rnorm(n), stringsAsFactors = FALSE)
  attrs$age[2] <- NA

  # shuffle the matrix's node order: `id` must put it back
  perm <- sample(n)
  shuffled <- m[perm, perm]
  fit <- suppressWarnings(netmice(attrs, list(friends = shuffled), id = "pid",
                                  m = 1, maxit = 2, seed = 1, printFlag = FALSE))
  expect_equal(unname(fit$networks$friends), unname(m))
  expect_false("pid" %in% names(fit$data))   # id is not an attribute

  # a matrix without dimnames cannot be matched by name
  bare <- unname(m)
  expect_error(netmice(attrs, list(friends = bare), id = "pid", m = 1, maxit = 1),
               "no row/column names")
  # nodes that do not line up
  wrong <- m
  dimnames(wrong) <- list(c("zz", ids[-1]), c("zz", ids[-1]))
  expect_error(netmice(attrs, list(friends = wrong), id = "pid", m = 1, maxit = 1),
               "same nodes as")
})

test_that("netmice: positional alignment is announced when `id` is NULL", {
  set.seed(7)
  n <- 10
  m <- matrix(rbinom(n * n, 1, 0.3), n, n)
  diag(m) <- 0
  m[1, 2] <- NA
  attrs <- data.frame(age = rnorm(n))
  attrs$age[3] <- NA
  expect_message(
    suppressWarnings(netmice(attrs, list(friends = m), m = 1, maxit = 2,
                             seed = 1, printFlag = TRUE)),
    "assumed to be in the same order")
})

# ---- end to end ------------------------------------------------------------

test_that("netmice: runs end to end from an edgelist and matches the matrix path", {
  set.seed(21)
  n <- 20
  ids <- paste0("p", seq_len(n))
  attrs <- data.frame(pid = ids, age = rnorm(n, 35, 8), stringsAsFactors = FALSE)
  attrs$age[sample(n, 3)] <- NA

  m <- matrix(rbinom(n * n, 1, 0.15), n, n, dimnames = list(ids, ids))
  diag(m) <- 0
  idx <- which(m == 1, arr.ind = TRUE)
  el <- data.frame(from = ids[idx[, 1]], to = ids[idx[, 2]],
                   stringsAsFactors = FALSE)
  miss <- c("p1", "p2")
  for (p in miss) {
    k <- match(p, ids)
    off <- m[k, ] == 0
    off[k] <- FALSE          # the diagonal is never a tie, missing or not
    m[k, off] <- NA
  }

  fit <- suppressWarnings(
    netmice(attrs, networks = el, id = "pid",
            edgelist_options = list(missing = miss),
            m = 1, maxit = 2, seed = 1, printFlag = FALSE))
  expect_equal(unname(fit$networks$net1), unname(m))
  done <- complete_netmice(fit, 1)
  expect_false(anyNA(done$networks$net1))
  expect_false(anyNA(done$data))
  expect_true("networks" %in% names(done))
})

test_that("netmice: `structural` accepts an edgelist and fixes those cells at 0", {
  set.seed(31)
  n <- 12
  ids <- paste0("p", seq_len(n))
  attrs <- data.frame(pid = ids, age = rnorm(n), stringsAsFactors = FALSE)
  attrs$age[2] <- NA
  m <- matrix(rbinom(n * n, 1, 0.2), n, n, dimnames = list(ids, ids))
  diag(m) <- 0
  m["p1", "p5"] <- 0
  m["p3", "p7"] <- NA
  el <- data.frame(from = ids[which(m == 1, arr.ind = TRUE)[, 1]],
                   to = ids[which(m == 1, arr.ind = TRUE)[, 2]],
                   stringsAsFactors = FALSE)
  struct_el <- data.frame(from = c("p1", "p3"), to = c("p5", "p7"),
                          stringsAsFactors = FALSE)
  fit <- suppressWarnings(
    netmice(attrs, networks = el, id = "pid", structural = struct_el,
            m = 1, maxit = 2, seed = 1, printFlag = FALSE))
  out <- complete_netmice(fit, 1)$networks$net1
  expect_equal(unname(out["p1", "p5"]), 0)
  expect_equal(unname(out["p3", "p7"]), 0)   # was NA, fixed rather than imputed
})

test_that("netmice: edgelist_options is ignored (with a note) for matrix input", {
  set.seed(8)
  n <- 10
  m <- matrix(rbinom(n * n, 1, 0.3), n, n)
  diag(m) <- 0
  m[1, 2] <- NA
  attrs <- data.frame(age = rnorm(n))
  attrs$age[3] <- NA
  expect_message(
    suppressWarnings(netmice(attrs, list(friends = m),
                             edgelist_options = list(edgelist_split = "wave"),
                             m = 1, maxit = 2, seed = 1, printFlag = FALSE)),
    "is ignored")
})

test_that(".validate_edgelist_options: rejects unknown entries and bad formats", {
  expect_error(.validate_edgelist_options(list(nodlist = 1)), "may only contain")
  expect_error(.validate_edgelist_options(list(edgelist_format = "sideways")))
  expect_error(.validate_edgelist_options(list(edgelist_split = 1)),
               "character vector")
  expect_error(.validate_edgelist_options("nope"), "must be a list")
  d <- .validate_edgelist_options(NULL)
  expect_equal(d$edgelist_format, "long")
  expect_equal(d$edgelist_split, character(0))
})
