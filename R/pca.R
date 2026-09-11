# How many predictors a model is allowed to carry. Replaces the former
# `n_components` argument and the internal .safe_max_cols() cap with one rule.

#' Validate and normalize a `PCA` specification
#'
#' @param PCA A list with `n` (a fixed maximum number of components) and/or
#'   `ratio` (rows-or-events per predictor). At least one must be set. `n = 0`
#'   is allowed and means "no components at all": only the protected columns
#'   survive (see `.pca_budget()`).
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
  check_num <- function(v, nm, allow_zero = FALSE) {
    if (is.null(v)) return(NULL)
    ok <- is.numeric(v) && length(v) == 1 && !is.na(v) &&
      (if (allow_zero) v >= 0 else v > 0)
    if (!ok) {
      stop("`PCA$", nm, "` must be a single ",
           if (allow_zero) "non-negative" else "positive", " number, or NULL.",
           if (allow_zero) " (`n = 0` means no components: only the protected columns are kept.)",
           call. = FALSE)
    }
    v
  }
  # n = 0 is a deliberate request for no collapsed predictors at all; a ratio
  # of 0 has no sensible meaning (infinitely many predictors per event)
  n <- check_num(PCA$n, "n", allow_zero = TRUE)
  ratio <- check_num(PCA$ratio, "ratio")
  if (is.null(n) && is.null(ratio)) {
    stop("`PCA` must set at least one of `n` (a fixed number of components) ",
         "and `ratio` (rows or events per predictor).", call. = FALSE)
  }
  # only an explicit 0 means "no components": a fractional n below 1 floors
  # to 1, as it always did, rather than silently becoming a zero request
  n_int <- if (is.null(n)) NULL else if (n == 0) 0L else
    max(1L, as.integer(floor(n)))
  list(n = n_int, ratio = ratio)
}

#' Resolve the attribute-model PCA budget
#'
#' `PCA` is otherwise a single call-level setting shared by the attribute and
#' the tie models. Their predictor sets are not comparable, though: an
#' attribute model is budgeted against observed *rows* (order 10^2 here) while
#' a tie model is budgeted against observed *events* (order 10^3), so one
#' `ratio` cannot be right for both. `PCA_attributes` lets the attribute side
#' be set independently, including switching the collapse off entirely -
#' useful when the attribute predictors are a short, deliberately chosen set
#' whose individual coefficients are the point.
#'
#' @param PCA_attributes `NULL` to inherit `PCA` (the default, and the only
#'   value that reproduces pre-0.9 behaviour); `"none"` to leave attribute
#'   predictors uncollapsed; or a `list(n=, ratio=)` used for attribute
#'   models only.
#' @param PCA The already-validated call-level `PCA` list.
#' @return A validated PCA list, or `NULL` meaning "impose no budget".
#'   `.clean_predictor_matrix()` treats `max_cols = NULL` as no cap, so `NULL`
#'   needs no special handling downstream.
#' @noRd
.resolve_attribute_pca <- function(PCA_attributes, PCA) {
  if (is.null(PCA_attributes)) return(PCA)
  if (is.character(PCA_attributes)) {
    if (!identical(PCA_attributes, "none")) {
      stop("`PCA_attributes` must be NULL (inherit `PCA`), \"none\" (no ",
           "collapse for attribute models), or a list like ",
           "PCA = list(n = 5).", call. = FALSE)
    }
    return(NULL)
  }
  .validate_pca(PCA_attributes)
}

#' Resolve the tie-model PCA budget
#'
#' The mirror of \code{\link{.resolve_attribute_pca}} for the other half of the
#' split it describes. Until 1.1.0 only the attribute side could be overridden,
#' which had it backwards: an attribute model is budgeted against observed rows
#' (order 10^2) while a tie model is budgeted against observed *events* (order
#' 10^3 but as few as ~10^2 on a sparse network), and it is the tie models that
#' run closest to their budget - a network target's protected dyad block grows
#' by two columns per additional network and cannot be collapsed.
#'
#' @param PCA_networks `NULL` to inherit `PCA` (the default, and the only value
#'   that reproduces pre-1.1.0 behaviour); `"none"` to impose no budget on tie
#'   models; or a `list(n=, ratio=)` used for tie models only.
#' @param PCA The already-validated call-level `PCA` list.
#' @return A validated PCA list, or `NULL` meaning "impose no budget".
#' @noRd
.resolve_network_pca <- function(PCA_networks, PCA) {
  if (is.null(PCA_networks)) return(PCA)
  if (is.character(PCA_networks)) {
    if (!identical(PCA_networks, "none")) {
      stop("`PCA_networks` must be NULL (inherit `PCA`), \"none\" (no ",
           "collapse for tie models), or a list like ",
           "PCA = list(n = 5).", call. = FALSE)
    }
    return(NULL)
  }
  .validate_pca(PCA_networks)
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
#' @return A non-negative integer: the most predictors this model may carry.
#'   A budget derived from `ratio` is never 0 - on a tiny sample that would
#'   silently empty every design matrix - so it floors at 1. The one way to get
#'   0 is to ask for it explicitly with `n = 0`, which keeps only the protected
#'   columns (a network target's endogenous and cross-network dyad terms, an
#'   attribute model's isolate flags and `models` terms) and collapses nothing.
#' @noRd
.pca_budget <- function(PCA, denom) {
  if (!is.null(PCA$n) && PCA$n == 0L) return(0L)
  caps <- integer(0)
  if (!is.null(PCA$n)) caps <- c(caps, PCA$n)
  if (!is.null(PCA$ratio)) caps <- c(caps, as.integer(floor(denom / PCA$ratio)))
  max(1L, if (length(caps)) min(caps) else 1L)
}

#' Validate and normalize a `net_ridge` specification
#'
#' The tie working model is a logistic (or linear) regression fitted on the
#' observed dyads and then applied to the missing ones, and it is routinely
#' fitted at a handful of events per predictor: a network target's protected
#' dyad block cannot be collapsed, so on a sparse network it can fill the whole
#' `PCA` budget on its own. An unpenalised fit there produces a linear
#' predictor whose spread is far too large - probabilities pinned near 0 and 1
#' - and because `plogis()` is convex below 0.5, an over-dispersed `eta`
#' inflates the *mean* imputed tie probability wherever ties are rare. A small
#' ridge penalty shrinks that spread, and shrinks the Bayesian coefficient draw
#' with it.
#'
#' @param net_ridge A list with `lambda` (>= 0; 0 disables the penalty) and
#'   `scale`, either `"fixed"` (use `lambda` as given) or `"epv"` (multiply it
#'   by predictors-per-event, so a model with little information to spend is
#'   penalised harder than a well-identified one).
#' @return A list with `lambda` (numeric) and `scale` (character).
#' @noRd
.validate_net_ridge <- function(net_ridge) {
  if (is.null(net_ridge)) return(list(lambda = 0, scale = "fixed"))
  if (is.numeric(net_ridge) && length(net_ridge) == 1L) {
    net_ridge <- list(lambda = net_ridge)
  }
  if (!is.list(net_ridge)) {
    stop("`net_ridge` must be a list, e.g. net_ridge = list(lambda = 0.01), ",
         "or a single number giving the penalty.", call. = FALSE)
  }
  bad <- setdiff(names(net_ridge), c("lambda", "scale"))
  if (length(bad)) {
    stop("`net_ridge` may only contain 'lambda' and 'scale'; got: ",
         toString(bad), ".", call. = FALSE)
  }
  lambda <- if (is.null(net_ridge$lambda)) 0.01 else net_ridge$lambda
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) ||
      lambda < 0) {
    stop("`net_ridge$lambda` must be a single non-negative number ",
         "(0 disables the penalty).", call. = FALSE)
  }
  scale <- if (is.null(net_ridge$scale)) "fixed" else net_ridge$scale
  if (!identical(scale, "fixed") && !identical(scale, "epv")) {
    stop("`net_ridge$scale` must be \"fixed\" or \"epv\".", call. = FALSE)
  }
  list(lambda = as.numeric(lambda), scale = scale)
}

#' The penalty actually applied to one tie model
#'
#' `scale = "epv"` turns the user's `lambda` into `lambda * p / events`, so the
#' penalty tracks how thin the data are for *this* target rather than being one
#' number across networks that differ in density by an order of magnitude.
#' @param net_ridge A validated `net_ridge` list.
#' @param p Number of non-intercept predictors in the design.
#' @param events Effective sample size from \code{.pca_denom()}.
#' @noRd
.net_ridge_lambda <- function(net_ridge, p, events) {
  lambda <- net_ridge$lambda
  if (lambda <= 0 || p < 1L) return(0)
  if (identical(net_ridge$scale, "epv")) {
    if (!is.finite(events) || events < 1) return(lambda)
    lambda <- lambda * p / events
  }
  lambda
}

#' Validate `net_endo_terms` and `net_gw_decay`
#'
#' `net_endo_terms` is a wish list, not a promise: which terms a given target
#' actually receives depends on whether it is binary and whether it is
#' undirected (see `.resolve_endo_terms()`). Filtering is silent rather than an
#' error because `networks` is routinely a heterogeneous list - a binary
#' directed nomination network beside a weighted or undirected one - and no
#' single set can be literally correct for all of them. Erroring would make the
#' default unusable on exactly the mixed input the package is built for.
#' @noRd
.validate_endo_terms <- function(net_endo_terms) {
  allowed <- c("reciprocity", "twopath", "gwesp", "gwodegree", "gwidegree")
  if (!is.character(net_endo_terms) || !length(net_endo_terms)) {
    stop("`net_endo_terms` must be a character vector with at least one of: ",
         toString(allowed), ".", call. = FALSE)
  }
  bad <- setdiff(net_endo_terms, allowed)
  if (length(bad)) {
    stop("`net_endo_terms` may only contain ", toString(allowed),
         "; got: ", toString(bad),
         ". (\"gwdegree\" is not set directly - an undirected binary target ",
         "gets it automatically in place of gwodegree/gwidegree.)",
         call. = FALSE)
  }
  unique(net_endo_terms)
}

#' @noRd
.validate_gw_decay <- function(net_gw_decay) {
  default <- list(esp = 0.69, outdegree = 0.69, indegree = 0.69)
  if (is.null(net_gw_decay)) return(default)
  if (is.numeric(net_gw_decay) && length(net_gw_decay) == 1L) {
    net_gw_decay <- list(esp = net_gw_decay, outdegree = net_gw_decay,
                         indegree = net_gw_decay)
  }
  if (!is.list(net_gw_decay)) {
    stop("`net_gw_decay` must be a list like ",
         "list(esp = 0.69, outdegree = 0.69, indegree = 0.69), or a single ",
         "number used for all three.", call. = FALSE)
  }
  bad <- setdiff(names(net_gw_decay), names(default))
  if (length(bad)) {
    stop("`net_gw_decay` may only contain 'esp', 'outdegree' and 'indegree'; ",
         "got: ", toString(bad), ".", call. = FALSE)
  }
  out <- default
  for (nm in names(net_gw_decay)) {
    v <- net_gw_decay[[nm]]
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v < 0) {
      stop("`net_gw_decay$", nm, "` must be a single non-negative number.",
           call. = FALSE)
    }
    out[[nm]] <- as.numeric(v)
  }
  # w = 1 - exp(-decay). As decay grows w -> 1 and every change statistic
  # collapses to a constant column (dropped, or near-constant with an enormous
  # standard error); as it shrinks w -> 0 and the degree terms degenerate into
  # an isolate indicator the package already provides elsewhere.
  extreme <- vapply(out, function(v) v > 0 && (v < 0.05 || v > 5), logical(1))
  if (any(extreme)) {
    warning("netimpute: `net_gw_decay` outside roughly [0.05, 5] for: ",
            toString(names(out)[extreme]),
            ". Very small decays reduce the degree terms to an isolate ",
            "indicator; very large ones make every change statistic nearly ",
            "constant. Typical ERGM practice is 0.25-1.5.",
            call. = FALSE)
  }
  out
}
