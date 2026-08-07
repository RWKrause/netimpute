# Adaptive structural constraints between networks (net_dependence):
# necessary ("if not A then not B") and forbidden ("if A then not B").

make_pair <- function(n = 30, seed = 1) {
  # A random binary; B only where A has a tie (consistent with necessary);
  # C disjoint from A (consistent with forbidden)
  set.seed(seed)
  A <- matrix(rbinom(n * n, 1, .25), n, n); diag(A) <- 0
  B <- A * matrix(rbinom(n * n, 1, .5), n, n); diag(B) <- 0
  C <- (1 - A) * matrix(rbinom(n * n, 1, .1), n, n); diag(C) <- 0
  list(A = A, B = B, C = C, n = n)
}

test_that(".validate_net_dependence rejects malformed specifications", {
  nets <- c("A", "B")
  expect_error(netimpute:::.validate_net_dependence(list(nec = list(c("A", "B"))), nets),
               "'necessary'")
  expect_error(netimpute:::.validate_net_dependence(
    list(necessary = list(c("A", "X"))), nets), "unknown network")
  expect_error(netimpute:::.validate_net_dependence(
    list(forbidden = list(c("A", "A"))), nets), "to itself")
  expect_error(netimpute:::.validate_net_dependence(
    list(necessary = list("A")), nets), "exactly 2")
  # a single unwrapped pair is normalized
  r <- netimpute:::.validate_net_dependence(list(necessary = c("A", "B")), nets)
  expect_identical(r$necessary, list(c("A", "B")))
  expect_null(netimpute:::.validate_net_dependence(NULL, nets))
})

test_that(".deduce_from_dependence fills determined cells and errors on contradictions", {
  p <- make_pair()
  A <- p$A; B <- p$B; C <- p$C
  # punch NAs into cells whose value is determined by the rules
  a0 <- which(A == 0 & row(A) != col(A))[1:5]   # A observed 0 -> B must be 0
  b1 <- which(B == 1)[1:5]                       # B observed 1 -> A must be 1
  a1 <- which(A == 1)[1:5]                       # A observed 1 -> C must be 0
  B[a0] <- NA; A[b1] <- NA; C[a1] <- NA
  rules <- netimpute:::.validate_net_dependence(
    list(necessary = list(c("A", "B")), forbidden = list(c("A", "C"))),
    c("A", "B", "C"))
  out <- netimpute:::.deduce_from_dependence(list(A = A, B = B, C = C), rules)
  expect_true(all(out$B[a0] == 0))
  expect_true(all(out$A[b1] == 1))
  expect_true(all(out$C[a1] == 0))

  # observed contradiction: B tie where A observed 0
  A2 <- p$A; B2 <- p$B
  cell <- which(A2 == 0 & row(A2) != col(A2))[1]
  B2[cell] <- 1
  expect_error(
    netimpute:::.deduce_from_dependence(
      list(A = A2, B = B2),
      netimpute:::.validate_net_dependence(list(necessary = list(c("A", "B"))),
                                           c("A", "B"))),
    "contradict the rule")

  # forbidden contradiction: both observed 1
  C2 <- p$C
  cell2 <- which(p$A == 1)[1]
  C2[cell2] <- 1
  expect_error(
    netimpute:::.deduce_from_dependence(
      list(A = p$A, C = C2),
      netimpute:::.validate_net_dependence(list(forbidden = list(c("A", "C"))),
                                           c("A", "C"))),
    "contradict the rule")
})

test_that("netmice respects necessary and forbidden rules in every completion", {
  p <- make_pair(seed = 7)
  A <- p$A; B <- p$B; C <- p$C; n <- p$n
  off <- which(row(A) != col(A))
  set.seed(8)
  A[sample(off, 60)] <- NA
  B[sample(off, 60)] <- NA
  C[sample(off, 60)] <- NA
  attrs <- data.frame(age = rnorm(n), happy = rnorm(n))
  attrs$happy[sample(n, 4)] <- NA

  fit <- suppressMessages(netmice(
    attrs, list(A = A, B = B, C = C), m = 2, maxit = 3,
    net_dependence = list(necessary = list(c("A", "B")),
                          forbidden = list(c("A", "C"))),
    seed = 11, printFlag = FALSE))

  expect_identical(fit$net_dependence$necessary, list(c("A", "B")))
  for (im in 1:2) {
    com <- complete_netmice(fit, im)$net_list
    expect_false(anyNA(com$A)); expect_false(anyNA(com$B)); expect_false(anyNA(com$C))
    # necessary: B ties only where A has a tie
    expect_true(all(com$B[com$A == 0] == 0))
    # forbidden: A and C never share a tie
    expect_true(all(com$A * com$C == 0))
    # originally observed cells are preserved
    expect_identical(com$A[!is.na(A)], A[!is.na(A)])
    expect_identical(com$B[!is.na(B)], B[!is.na(B)])
    expect_identical(com$C[!is.na(C)], C[!is.na(C)])
  }
})

test_that("net_dependence deductions make deduced cells observed (never re-imputed)", {
  p <- make_pair(seed = 12)
  A <- p$A; B <- p$B
  b1 <- which(B == 1)[1:8]
  A[b1] <- NA                       # deducible: must be 1 (B observed 1)
  # plus genuinely free missing cells so the chain actually runs
  free <- setdiff(which(row(A) != col(A) & B == 0), b1)
  set.seed(13)
  A[sample(free, 15)] <- NA
  fit <- suppressMessages(netmice(
    data.frame(age = rnorm(p$n)),
    list(A = A, B = B), m = 2, maxit = 2,
    net_dependence = list(necessary = list(c("A", "B"))),
    seed = 3, printFlag = FALSE))
  # the deduced cells are already filled in the stored (pre-chain) net_list
  expect_true(all(fit$net_list$A[b1] == 1))
  for (im in 1:2) {
    expect_true(all(complete_netmice(fit, im)$net_list$A[b1] == 1))
  }
})
