# Adaptive (network-dependent) structural constraints between networks:
# "necessary" (a tie in A is required for a tie in B - B's ties can only
# exist on A's support) and "forbidden" (ties in A and B are mutually
# exclusive). Unlike `structural` - a fixed set of zero-by-design cells -
# these constraints are re-evaluated against the CURRENT state of the
# partner network at every visit of the imputation chain, so the set of
# determined cells adapts as the partner is itself re-imputed.
#
# Three mechanisms, all driven by the same rule list:
#   1. Deduction (once, before the chain): missing cells that are logically
#      determined by OBSERVED cells are filled in and thereafter count as
#      observed (necessary: A observed 0 => B = 0, B observed tie => A = 1;
#      forbidden: A observed tie => B = 0 and vice versa). Observed cells
#      that contradict a rule raise an error.
#   2. Dynamic exclusion (every visit of a constrained target): cells of the
#      target currently determined to be 0 by the partner's current
#      (possibly imputed) state are fixed at 0 and excluded from the
#      dyad-level working model, exactly like `structural` cells - they are
#      design-zeros GIVEN the partner, so fitting on them would inflate the
#      zeros of the model.
#   3. Propagation (after every network re-imputation): freshly imputed ties
#      elsewhere may newly determine cells of dependent networks; those
#      (non-observed) cells are re-zeroed to a fixpoint so the current state
#      always satisfies every rule.

#' Validate and normalize a `net_dependence` specification
#'
#' Returns `NULL` (no rules) or a list with elements `necessary` and
#' `forbidden`, each a (possibly empty) list of length-2 character vectors.
#' For `necessary`, the order is `c(parent, child)`: a tie in `parent` is
#' required for a tie in `child`. `forbidden` pairs are unordered.
#' @noRd
.validate_net_dependence <- function(net_dependence, net_names) {
  if (is.null(net_dependence)) return(NULL)
  if (!is.list(net_dependence) || is.null(names(net_dependence)) ||
      length(setdiff(names(net_dependence), c("necessary", "forbidden")))) {
    stop("`net_dependence` must be a named list with elements 'necessary' ",
         "and/or 'forbidden' (each a list of length-2 character vectors of ",
         "network names).", call. = FALSE)
  }
  norm_pairs <- function(x, label) {
    if (is.null(x)) return(list())
    if (is.character(x) && length(x) == 2) x <- list(x)
    if (!is.list(x)) {
      stop("`net_dependence$", label, "` must be a list of length-2 ",
           "character vectors (or a single such vector).", call. = FALSE)
    }
    lapply(x, function(p) {
      if (!is.character(p) || length(p) != 2) {
        stop("Each `net_dependence$", label, "` rule must be a character ",
             "vector of exactly 2 network names.", call. = FALSE)
      }
      unknown <- setdiff(p, net_names)
      if (length(unknown)) {
        stop("`net_dependence` (", label, ") references unknown network(s): ",
             toString(unknown), call. = FALSE)
      }
      if (p[1] == p[2]) {
        stop("`net_dependence` (", label, ") relates network '", p[1],
             "' to itself.", call. = FALSE)
      }
      p
    })
  }
  out <- list(necessary = norm_pairs(net_dependence$necessary, "necessary"),
              forbidden = norm_pairs(net_dependence$forbidden, "forbidden"))
  if (!length(out$necessary) && !length(out$forbidden)) return(NULL)
  out
}

#' Fill missing cells that are logically determined by observed cells
#'
#' Iterates the deduction rules to a fixpoint (a fill can enable further
#' fills through chained rules). Errors when observed cells contradict a
#' rule. Deducing a *required tie* (necessary rule, child observed tie,
#' parent missing) is only possible for binary parents - the tie's existence
#' is implied but a weight is not - so for weighted parents that deduction
#' is skipped with a warning.
#' @noRd
.deduce_from_dependence <- function(mats, rules) {
  if (is.null(rules)) return(mats)
  off <- row(mats[[1]]) != col(mats[[1]])
  is_binary <- vapply(mats, function(m) all(m[!is.na(m)] %in% c(0, 1)),
                      logical(1))

  for (p in rules$necessary) {
    if (!is_binary[p[1]] &&
        any(off & !is.na(mats[[p[2]]]) & mats[[p[2]]] != 0 & is.na(mats[[p[1]]]))) {
      warning("netimpute: `net_dependence` (necessary): network '", p[1],
              "' is weighted, so required ties implied by observed ties in '",
              p[2], "' cannot be deduced (their weights are unknown); those ",
              "cells stay missing and the rule is not guaranteed for them.",
              call. = FALSE)
    }
  }

  repeat {
    changed <- FALSE
    for (p in rules$necessary) {
      A <- mats[[p[1]]]; B <- mats[[p[2]]]
      bad <- off & !is.na(A) & A == 0 & !is.na(B) & B != 0
      if (any(bad)) {
        stop("`net_dependence` (necessary): network '", p[2], "' has ",
             sum(bad), " observed tie(s) in cell(s) where '", p[1],
             "' is observed 0 - the observed data contradict the rule.",
             call. = FALSE)
      }
      fill0 <- off & !is.na(A) & A == 0 & is.na(B)
      if (any(fill0)) {
        B[fill0] <- 0; mats[[p[2]]] <- B; changed <- TRUE
      }
      need1 <- off & !is.na(B) & B != 0 & is.na(A)
      if (any(need1) && is_binary[p[1]]) {
        A[need1] <- 1; mats[[p[1]]] <- A; changed <- TRUE
      }
    }
    for (p in rules$forbidden) {
      A <- mats[[p[1]]]; B <- mats[[p[2]]]
      bad <- off & !is.na(A) & A != 0 & !is.na(B) & B != 0
      if (any(bad)) {
        stop("`net_dependence` (forbidden): networks '", p[1], "' and '",
             p[2], "' have ", sum(bad),
             " cell(s) where both have an observed tie - the observed data ",
             "contradict the rule.", call. = FALSE)
      }
      fillB <- off & !is.na(A) & A != 0 & is.na(B)
      if (any(fillB)) {
        B[fillB] <- 0; mats[[p[2]]] <- B; changed <- TRUE
      }
      fillA <- off & !is.na(B) & B != 0 & is.na(A)
      if (any(fillA)) {
        A[fillA] <- 0; mats[[p[1]]] <- A; changed <- TRUE
      }
    }
    if (!changed) break
  }
  mats
}

#' Cells of `target` currently determined to be 0 by its partners
#'
#' Evaluated against the *current* (complete, possibly imputed) state of
#' the other networks. Returns `NULL` when no cell is determined.
#' @noRd
.dependence_zero_mask <- function(cur_mats, rules, target) {
  if (is.null(rules)) return(NULL)
  mask <- matrix(FALSE, nrow(cur_mats[[target]]), ncol(cur_mats[[target]]))
  for (p in rules$necessary) {
    if (p[2] == target) mask <- mask | (cur_mats[[p[1]]] == 0)
  }
  for (p in rules$forbidden) {
    if (p[1] == target) mask <- mask | (cur_mats[[p[2]]] != 0)
    if (p[2] == target) mask <- mask | (cur_mats[[p[1]]] != 0)
  }
  if (!any(mask)) return(NULL)
  mask
}

#' Re-zero non-observed cells now determined 0, to a fixpoint
#'
#' Observed cells are never touched (deduction already guaranteed they are
#' consistent); only imputed/filled cells are re-zeroed. Returns the updated
#' matrices and the names of the networks that changed (so callers can
#' refresh derived igraph views).
#' @noRd
.enforce_dependence <- function(cur_mats, rules, ry_nets) {
  changed_nets <- character(0)
  if (is.null(rules)) return(list(mats = cur_mats, changed = changed_nets))
  repeat {
    changed <- FALSE
    for (nm in names(cur_mats)) {
      mask <- .dependence_zero_mask(cur_mats, rules, nm)
      if (is.null(mask)) next
      tozero <- mask & !ry_nets[[nm]] & cur_mats[[nm]] != 0
      if (any(tozero)) {
        cur_mats[[nm]][tozero] <- 0
        changed_nets <- union(changed_nets, nm)
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  list(mats = cur_mats, changed = changed_nets)
}
