# How many predictors a model is allowed to carry. Replaces the former
# `n_components` argument and the internal .safe_max_cols() cap with one rule.

#' Validate and normalize a `PCA` specification
#'
#' @param PCA A list with `n` (a fixed maximum number of components) and/or
#'   `ratio` (rows-or-events per predictor). At least one must be set.
#' @return A list with `n` (integer or `NULL`) and `ratio` (numeric or `NULL`).
#' @noRd
.validate_pca <- function(PCA) {
  if (!is.list(PCA)) {
    stop("`PCA` must be a list, e.g. PCA = list(n = 5) for at most five ",
         "components, PCA = list(ratio = 20) for a 20:1 budget, or both.",
         call. = FALSE)
  }
  bad <- setdiff(names(PCA), c("n", "ratio"))
  if (length(bad)) {
    stop("`PCA` may only contain 'n' and 'ratio'; got: ", toString(bad), ".",
         call. = FALSE)
  }
  check_num <- function(v, nm) {
    if (is.null(v)) return(NULL)
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v <= 0) {
      stop("`PCA$", nm, "` must be a single positive number, or NULL.",
           call. = FALSE)
    }
    v
  }
  n <- check_num(PCA$n, "n")
  ratio <- check_num(PCA$ratio, "ratio")
  if (is.null(n) && is.null(ratio)) {
    stop("`PCA` must set at least one of `n` (a fixed number of components) ",
         "and `ratio` (rows or events per predictor).", call. = FALSE)
  }
  list(n = if (is.null(n)) NULL else as.integer(floor(n)), ratio = ratio)
}

#' Effective sample size backing one model's predictors
#'
#' For a binary target this is the *minority class count* among the observed
#' values - the events-per-variable rule (Peduzzi et al. 1996; Harrell 2015).
#' A logistic model is limited by its rarer outcome, not by its row count: a
#' 30-node network at density 0.15 has ~870 dyad rows but only ~130 ties, and
#' budgeting against 870 would allow enough predictors to separate the data
#' perfectly. Continuous and multinomial targets use the observed row count.
#'
#' @param y Target values.
#' @param ry Logical mask of observed entries, or `NULL` if `y` is already
#'   restricted to the observed ones.
#' @param binary Whether the target is binary.
#' @noRd
.pca_denom <- function(y, ry, binary) {
  obs <- if (is.null(ry)) y else y[ry]
  obs <- obs[!is.na(obs)]
  if (!length(obs)) return(0L)
  if (isTRUE(binary)) {
    tb <- table(obs)
    # only one class observed: nothing to discriminate, no budget to spend
    if (length(tb) < 2L) return(0L)
    return(as.integer(min(tb)))
  }
  length(obs)
}

#' Predictor budget for one model
#'
#' @param PCA A validated `PCA` list.
#' @param denom Effective sample size from \code{.pca_denom()}.
#' @return A positive integer: the most predictors this model may carry. Never
#'   0 - an intercept-only fallback is handled downstream, and a budget of 0
#'   would silently empty every design matrix on tiny samples.
#' @noRd
.pca_budget <- function(PCA, denom) {
  caps <- integer(0)
  if (!is.null(PCA$n)) caps <- c(caps, PCA$n)
  if (!is.null(PCA$ratio)) caps <- c(caps, as.integer(floor(denom / PCA$ratio)))
  max(1L, if (length(caps)) min(caps) else 1L)
}
