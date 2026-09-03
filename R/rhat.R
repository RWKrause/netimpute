# Rank-normalized split-R-hat (Vehtari, Gelman, Simpson, Carpenter & Buerkner
# 2021, Bayesian Analysis 16(2), 667-718), implemented here rather than taken
# from rstan/posterior so the package keeps its small dependency footprint.
# Numerically identical to rstan::Rhat on the post-burn-in draws.

#' Rank-normalization: ranks -> standard normal scores
#'
#' Uses the midpoint plotting position (r - 1/2)/S, which is what
#' rstan:::z_scale() applies. (posterior uses the Blom position
#' (r - 3/8)/(S + 1/4) instead; the two differ by a few 1e-2 in R-hat on
#' short chains, so it matters which reference the numbers are meant to
#' reproduce.)
#' @noRd
.rhat_z_scale <- function(x) {
  s <- length(x)
  r <- rank(x, ties.method = "average")
  z <- stats::qnorm((r - 1 / 2) / s)
  dim(z) <- dim(x)
  z
}

#' Plain (already-split, already-normalized) Gelman-Rubin statistic
#' @noRd
.rhat_core <- function(z) {
  n <- nrow(z)
  m <- ncol(z)
  if (n < 2L || m < 2L) return(NA_real_)
  chain_means <- colMeans(z)
  chain_vars <- apply(z, 2, stats::var)
  w <- mean(chain_vars)
  # A constant trace has W = 0 and R-hat = 0/0. This is not hypothetical:
  # n_isolates is flat at 0 in many networks. NA, never a warning.
  if (!is.finite(w) || w <= 0) return(NA_real_)
  b <- n * stats::var(chain_means)
  var_hat <- ((n - 1) / n) * w + b / n
  sqrt(var_hat / w)
}

#' Rank-normalized split-R-hat for one tracked quantity
#'
#' Discards the first `burnin` fraction of iterations, splits each remaining
#' chain in half (so a within-chain trend shows up as between-chain
#' disagreement), and returns the larger of the bulk statistic and the
#' folded-tail statistic, as in Vehtari et al. (2021).
#'
#' @param x An `iterations x chains` matrix of tracked values.
#' @param burnin Fraction of the leading iterations to discard (default 0.5).
#' @return A single R-hat, or `NA_real_` when the trace is constant, has
#'   non-finite values, or is too short to split.
#' @noRd
.rhat_split <- function(x, burnin = 0.5) {
  x <- as.matrix(x)
  if (nrow(x) < 4L || ncol(x) < 1L) return(NA_real_)
  keep <- seq.int(floor(nrow(x) * burnin) + 1L, nrow(x))
  x <- x[keep, , drop = FALSE]
  half <- floor(nrow(x) / 2)
  # each half-chain needs at least 2 draws for a within-chain variance
  if (half < 2L) return(NA_real_)
  if (!all(is.finite(x))) return(NA_real_)
  # split: the two halves of every chain become separate chains, so one chain
  # drifting across iterations inflates B even with a single chain supplied.
  # An odd number of retained iterations drops the middle one, as rstan does.
  split_chains <- function(mm) {
    cbind(mm[seq_len(half), , drop = FALSE],
          mm[seq.int(nrow(mm) - half + 1L, nrow(mm)), , drop = FALSE])
  }
  bulk <- .rhat_core(.rhat_z_scale(split_chains(x)))
  # folded: |x - median(x)| targets the tails/scale, catching chains that
  # agree on the centre while disagreeing on spread. Folded on the UNSPLIT
  # draws, so the median matches the one rstan::Rhat uses.
  folded <- .rhat_core(.rhat_z_scale(split_chains(abs(x - stats::median(x)))))
  if (is.na(bulk) && is.na(folded)) return(NA_real_)
  max(bulk, folded, na.rm = TRUE)
}

#' R-hat for every quantity a netmids object tracks
#'
#' @param x A `netmids` object.
#' @return Named numeric vector of R-hats (`NA` where undefined).
#' @noRd
.netmids_rhat <- function(x) {
  out <- numeric(0)
  add <- function(vals, nms) {
    if (length(vals)) out <<- c(out, stats::setNames(vals, nms))
  }
  # attributes: mean and variance of the imputed values
  for (v in x$var_missing) {
    add(.rhat_split(x$chainMean[v, , ]), paste0(v, " (mean)"))
    add(.rhat_split(x$chainVar[v, , ]),  paste0(v, " (variance)"))
  }
  # whole-network diagnostics
  nets <- dimnames(x$netChain)[[1]]
  diags <- dimnames(x$netChain)[[2]]
  for (nm in nets) {
    for (dn in diags) {
      add(.rhat_split(x$netChain[nm, dn, , ]), paste0(nm, " (", dn, ")"))
    }
  }
  # imputed ties
  for (nm in x$net_missing) {
    add(.rhat_split(x$netImpMean[nm, , ]), paste0(nm, " (imputed ties, mean)"))
    add(.rhat_split(x$netImpVar[nm, , ]),  paste0(nm, " (imputed ties, variance)"))
  }
  out
}

#' Warn when any R-hat exceeds the threshold
#' @noRd
.warn_rhat <- function(rhat, maxit, threshold = 1.05) {
  bad <- rhat[is.finite(rhat) & rhat > threshold]
  if (!length(bad)) return(invisible(FALSE))
  bad <- sort(bad, decreasing = TRUE)
  shown <- utils::head(bad, 5L)
  warning(
    "netimpute: ", length(bad), " tracked quantity/quantities have R-hat > ",
    threshold, ", so the chains have not converged: ",
    paste0(names(shown), " = ", format(round(shown, 3), nsmall = 3),
           collapse = ", "),
    if (length(bad) > length(shown))
      paste0(", and ", length(bad) - length(shown), " more"),
    ". Increase `maxit` (currently ", maxit, ") and refit.",
    call. = FALSE
  )
  invisible(TRUE)
}
