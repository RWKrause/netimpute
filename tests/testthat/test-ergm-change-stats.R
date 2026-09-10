# Change statistics for the endogenous tie terms.
#
# The reference implementations below are written straight from the ERGM
# statistic definitions and share no code with .ergm_change_stats(); change
# statistics are obtained by finite difference, i.e. by literally toggling the
# cell and re-evaluating. That is the only check that can catch a wrong
# formula, as opposed to a formula applied inconsistently - the sweep's
# eta_base arithmetic is exact whatever the columns contain, so a bad
# statistic produces no arithmetic symptom anywhere.

esp_of  <- function(B) { diag(B) <- 0; B %*% B }
u_gwesp <- function(B, a) { w <- 1 - exp(-a); exp(a) * sum(B * (1 - w^esp_of(B))) }
u_gwod  <- function(B, a) { w <- 1 - exp(-a); exp(a) * sum(1 - w^rowSums(B)) }
u_gwid  <- function(B, a) { w <- 1 - exp(-a); exp(a) * sum(1 - w^colSums(B)) }

delta_ref <- function(u, B, i, j, a) {
  B0 <- B; B0[i, j] <- 0
  B1 <- B; B1[i, j] <- 1
  u(B1, a) - u(B0, a)
}

test_that(".ergm_change_stats matches a finite-difference reference everywhere", {
  for (dens in c(0, 0.15, 0.5, 0.9, 1)) {
    for (a in c(0, 0.25, 0.69, 1, 3)) {
      set.seed(round(1000 * dens + 10 * a))
      n <- 8
      B <- matrix(rbinom(n * n, 1, dens), n, n); diag(B) <- 0
      S <- netimpute:::.ergm_change_stats(
        B, list(esp = a, outdegree = a, indegree = a))
      for (i in seq_len(n)) for (j in seq_len(n)) {
        if (i == j) next
        # sweeping every cell is what exercises the "- b" correction: it only
        # bites where the cell is already a tie
        expect_equal(S$gwesp[i, j],     delta_ref(u_gwesp, B, i, j, a))
        expect_equal(S$gwodegree[i, j], delta_ref(u_gwod,  B, i, j, a))
        expect_equal(S$gwidegree[i, j], delta_ref(u_gwid,  B, i, j, a))
      }
    }
  }
})

test_that(".ergm_change_stats: the undirected identities hold", {
  set.seed(4); n <- 9
  B <- matrix(rbinom(n * n, 1, 0.3), n, n)
  B[lower.tri(B)] <- t(B)[lower.tri(B)]; diag(B) <- 0
  S <- netimpute:::.ergm_change_stats(
    B, list(esp = 0.69, outdegree = 0.69, indegree = 0.69))
  expect_equal(S$gwesp, t(S$gwesp))
  expect_equal(S$gwodegree, t(S$gwidegree))
  # adding an undirected tie raises BOTH endpoints' degrees, so the undirected
  # term is the sum - not gwodegree renamed, which would halve it
  expect_equal(S$gwdegree, S$gwodegree + S$gwidegree)

  # the directed OTP expression IS the undirected ESP change statistic
  u_und <- function(B, a) {
    w <- 1 - exp(-a)
    exp(a) * sum((B * (1 - w^esp_of(B)))[upper.tri(B)])
  }
  for (i in seq_len(n)) for (j in seq_len(n)) if (i < j) {
    B0 <- B; B0[i, j] <- 0; B0[j, i] <- 0
    B1 <- B; B1[i, j] <- 1; B1[j, i] <- 1
    expect_equal(S$gwesp[i, j], u_und(B1, 0.69) - u_und(B0, 0.69))
  }
})

test_that(".ergm_change_stats: decay 0 reduces gwesp's direct term to twopath", {
  set.seed(7)
  B <- matrix(rbinom(64, 1, 0.3), 8, 8); diag(B) <- 0
  S <- netimpute:::.ergm_change_stats(B, list(esp = 0, outdegree = 1, indegree = 1))
  # at decay 0 the direct term is 1(ESP > 0), the statistic twopath encoded
  direct <- exp(0) * (1 - 0^(B %*% B))
  expect_equal(as.numeric(direct), as.numeric((B %*% B) > 0))
  expect_true(all(is.finite(S$gwesp)))
})

test_that(".binarize_target: weights and NAs both count as presence/absence", {
  m <- matrix(0, 3, 3)
  m[2, 1] <- 2.5    # weight
  m[3, 1] <- NA     # unobserved
  m[2, 3] <- -3     # negative weight
  b <- netimpute:::.binarize_target(m)
  expect_true(all(b %in% c(0, 1)))
  expect_equal(b[2, 1], 1)   # any non-zero weight is a tie
  expect_equal(b[3, 1], 0)   # NA counts as absent, and does not propagate
  expect_equal(b[2, 3], 1)   # a negative weight is still a tie
  expect_true(all(diag(b) == 0))
  expect_false(anyNA(b))
})

test_that(".resolve_endo_terms: the term set depends on target type", {
  want <- c("reciprocity", "gwesp", "gwodegree", "gwidegree")
  # binary directed: everything as asked
  expect_setequal(
    netimpute:::.resolve_endo_terms(want, binary = TRUE, undirected = FALSE),
    want)
  # binary undirected: no reciprocity (it equals y), one summed degree term
  expect_setequal(
    netimpute:::.resolve_endo_terms(want, binary = TRUE, undirected = TRUE),
    c("gwesp", "gwdegree"))
  # weighted: the gw terms are binary-ERGM statistics, so twopath stands in
  expect_setequal(
    netimpute:::.resolve_endo_terms(want, binary = FALSE, undirected = FALSE),
    c("reciprocity", "twopath"))
  expect_setequal(
    netimpute:::.resolve_endo_terms(want, binary = FALSE, undirected = TRUE),
    "twopath")
  # an explicit legacy request is honoured as-is for a binary directed target
  expect_setequal(
    netimpute:::.resolve_endo_terms(c("reciprocity", "twopath"),
                                    binary = TRUE, undirected = FALSE),
    c("reciprocity", "twopath"))
})

test_that(".validate_endo_terms / .validate_gw_decay reject bad input", {
  expect_error(netimpute:::.validate_endo_terms("nope"), "may only contain")
  expect_error(netimpute:::.validate_endo_terms(character(0)), "at least one")
  expect_error(netimpute:::.validate_endo_terms("gwdegree"), "not set directly")
  expect_equal(netimpute:::.validate_gw_decay(NULL)$esp, 0.69)
  expect_equal(netimpute:::.validate_gw_decay(1.2),
               list(esp = 1.2, outdegree = 1.2, indegree = 1.2))
  expect_equal(netimpute:::.validate_gw_decay(list(esp = 0.5))$outdegree, 0.69)
  expect_error(netimpute:::.validate_gw_decay(list(esp = -1)), "non-negative")
  expect_error(netimpute:::.validate_gw_decay(list(nope = 1)), "may only contain")
  expect_warning(netimpute:::.validate_gw_decay(list(esp = 40)), "0.05")
})
