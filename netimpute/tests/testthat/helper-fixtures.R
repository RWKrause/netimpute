# Shared fixtures for the netimpute test suite: node attributes covering all
# three supported scales (continuous, binary, multinomial), and networks
# covering all three supported tie types (binary, weighted, signed), both
# directed and undirected.

fx_attrs <- function(n = 24, seed = 101) {
  set.seed(seed)
  data.frame(
    id     = paste0("n", seq_len(n)),
    age    = round(rnorm(n, 40, 10)),
    status = sample(c("active", "inactive"), n, replace = TRUE),
    dept   = sample(c("sales", "eng", "hr", "marketing"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# Binary directed network (0/1), guaranteed non-trivial density and at least
# one isolate, so isolate-handling is always exercised.
fx_bin_directed <- function(n = 24, seed = 201) {
  set.seed(seed)
  m <- matrix(rbinom(n * n, 1, 0.15), n, n)
  diag(m) <- 0
  m[1, ] <- 0
  m[, 1] <- 0
  m
}

# Binary undirected network (symmetric 0/1).
fx_bin_undirected <- function(n = 24, seed = 202) {
  set.seed(seed)
  m <- matrix(rbinom(n * n, 1, 0.15), n, n)
  diag(m) <- 0
  m <- ((m + t(m)) > 0) * 1
  m
}

# Weighted directed network: non-negative integer tie strengths.
fx_weighted <- function(n = 24, seed = 203) {
  set.seed(seed)
  m <- matrix(rpois(n * n, 0.7), n, n)
  diag(m) <- 0
  m
}

# Signed undirected network: ties in {-1, 0, 1} (e.g. distrust/no tie/trust).
fx_signed <- function(n = 24, seed = 204) {
  set.seed(seed)
  m <- matrix(sample(c(-1, 0, 0, 0, 1), n * n, replace = TRUE), n, n)
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  diag(m) <- 0
  m
}

fx_nets <- function(n = 24) {
  list(
    friends_bin    = fx_bin_directed(n),
    colleagues_bin = fx_bin_undirected(n),
    advice_weighted = fx_weighted(n),
    trust_signed   = fx_signed(n)
  )
}
