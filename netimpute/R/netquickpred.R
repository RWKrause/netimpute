# Per-target predictor selection: a network-aware, recursive generalization
# of mice::quickpred(). Instead of feeding every variable (or principal
# components of everything) into each univariate imputation model,
# netquickpred() screens candidates by absolute pairwise-complete Pearson
# correlation with each target - and, like quickpred, with the target's
# missingness indicator - then recursively adds predictors-of-predictors up
# to `steps` levels, and prunes highly collinear survivors per target.
# Attribute targets are screened at the node level (other attributes plus
# network-derived measures); network targets at the dyad level (the
# dyad_regression() design matrix). netmice() consumes the result via
# `predictor_selection = "quickpred"` (or a precomputed object), replacing
# the everything-plus-PCA-safeguard default with lean, target-specific
# models.

#' Apply netmice()'s network-feature naming rule to raw measure names
#'
#' Measure names are prefixed with the network name only when there is more
#' than one network, or when the bare name would collide with a name in
#' `clash_names` (the raw attribute columns) - matching the convention used
#' inside the netmice() attribute visits, so names selected here and names
#' built there always agree.
#' @keywords internal
.qp_prefix <- function(raw, net_name, multi, clash_names) {
  if (multi) return(paste0(net_name, "_", raw))
  ifelse(raw %in% clash_names, paste0(net_name, "_", raw), raw)
}

#' The homophily-block column names net_measures() produces for one attribute
#' (unprefixed), by attribute type - must mirror .homophily_block()
#' @keywords internal
.qp_homophily_names <- function(nm, type) {
  if (type == "continuous") {
    paste0(nm, "_", c("diff_out", "diff_in", "alter_mean", "alter_mean_in",
                      "alter_mean_out", "alter_min", "alter_max"))
  } else {
    paste0(nm, "_", c("ei_index", "blau_index", "alter_mode"))
  }
}

#' Node-level network-feature frame across all networks, with netmice()'s
#' naming convention
#'
#' Computes net_measures() for every network and cbinds the results, with
#' names prefixed per `.qp_prefix()`. `exclude` names (unprefixed) are
#' removed per network before prefixing - used by the netmice() attribute
#' visits to withhold the target's ego-value-dependent homophily measures.
#' `clash_names` fixes the prefix-on-collision rule to a stable attribute
#' set, so feature names do not change when netmice() drops unselected
#' attribute columns.
#' @keywords internal
.net_feature_frame <- function(gs, data, measure_set, attr_types, net_names,
                               exclude = character(0),
                               clash_names = names(data)) {
  lst <- lapply(seq_along(gs), function(k) {
    out <- net_measures(gs[[k]], data, measure_set = measure_set,
                        attr_types = attr_types)
    out <- out[setdiff(names(out), c("node_id", exclude))]
    names(out) <- .qp_prefix(names(out), net_names[k],
                             multi = length(net_names) > 1,
                             clash_names = clash_names)
    out
  })
  do.call(cbind, lst)
}

#' Expand a data.frame into a numeric screening matrix with a parent map
#'
#' Continuous columns pass through as one numeric column; categorical
#' columns become one dummy per non-reference level (`NA` propagates, so
#' pairwise-complete correlations see the true missingness). `parent[j]`
#' names the original column that produced matrix column j, so correlations
#' can be aggregated back to whole variables (a variable's relevance is the
#' max over its columns).
#' @keywords internal
.qp_expand_columns <- function(df, types = NULL) {
  cols <- list()
  parent <- character(0)
  for (nm in names(df)) {
    x <- df[[nm]]
    if (is.logical(x)) x <- as.numeric(x)
    type <- if (!is.null(types)) types[[nm]] else
      if (is.numeric(x)) "continuous" else "multinomial"
    if (type == "continuous") {
      cols[[nm]] <- as.numeric(x)
      parent <- c(parent, nm)
    } else {
      xf <- factor(x)
      lv <- levels(xf)
      use_lv <- if (length(lv) > 1) lv[-1] else lv
      for (l in use_lv) {
        cn <- paste0(nm, "_", l)
        cols[[cn]] <- as.numeric(xf == l)
        parent <- c(parent, nm)
      }
    }
  }
  if (!length(cols)) {
    return(list(mat = matrix(numeric(0), nrow(df), 0), parent = character(0)))
  }
  mat <- do.call(cbind, cols)
  colnames(mat) <- names(cols)
  list(mat = mat, parent = parent)
}

#' Dyad-level term names and their parent variables for one network target -
#' must mirror .build_dyad_data()'s naming
#' @keywords internal
.qp_dyad_parents <- function(term_names, data, types, net_names, target) {
  par <- stats::setNames(rep(NA_character_, length(term_names)), term_names)
  kind <- stats::setNames(rep("endo", length(term_names)), term_names)
  for (nm in names(data)) {
    tn <- if (types[[nm]] == "continuous") {
      paste0(nm, c("_ego", "_alter", "_absdiff"))
    } else {
      lv <- levels(factor(data[[nm]]))[-1]
      c(paste0(nm, "_ego_", lv), paste0(nm, "_alter_", lv), paste0(nm, "_same"))
    }
    hit <- intersect(tn, term_names)
    par[hit] <- nm
    kind[hit] <- "attr"
  }
  for (onm in setdiff(net_names, target)) {
    hit <- intersect(paste0(onm, c("_tie", "_recip", "_ego_outdeg",
                                   "_ego_indeg", "_alter_indeg",
                                   "_alter_outdeg")),
                     term_names)
    par[hit] <- onm
    kind[hit] <- "net"
  }
  list(parent = par, kind = kind)
}

#' Prune highly collinear predictor units for one target
#'
#' A "unit" is one whole predictor (an attribute spanning several dummy
#' columns, one network feature, or one dyad-level term). `exempt` units are
#' never dropped (a network target's endogenous reciprocity/twopath terms
#' and the cross-network `*_tie`/`*_recip` terms, whose named coefficients
#' carry substantive dependence structure). Two methods:
#' \describe{
#'   \item{pairwise}{While any unit pair's max absolute between-column
#'     correlation exceeds `threshold` (default 0.9), drop the member of the
#'     worst pair whose relevance for the target is weaker
#'     (caret::findCorrelation logic).}
#'   \item{vif}{While any unit's variance-inflation factor (max over its
#'     columns; VIF computed as the diagonal of the inverse correlation
#'     matrix of all remaining columns) exceeds `threshold` (default 10),
#'     drop the worst unit. Catches multi-variable collinearity that
#'     pairwise correlations miss.}
#' }
#' @keywords internal
.qp_prune_units <- function(units, unit_cols, X, Cmat, relevance, exempt,
                            method, threshold) {
  if (method == "none" || length(units) < 2) return(units)

  if (method == "pairwise") {
    k <- length(units)
    pc <- matrix(0, k, k, dimnames = list(units, units))
    for (a in seq_len(k - 1)) {
      for (b in seq(a + 1, k)) {
        v <- max(abs(Cmat[unit_cols[[units[a]]], unit_cols[[units[b]]]]))
        pc[a, b] <- pc[b, a] <- v
      }
    }
    keep <- units
    repeat {
      if (length(keep) < 2) break
      sub <- pc[keep, keep, drop = FALSE]
      both_ex <- outer(keep %in% exempt, keep %in% exempt, "&")
      sub[both_ex] <- 0
      mx <- max(sub)
      if (mx <= threshold) break
      idx <- which(sub == mx, arr.ind = TRUE)[1, ]
      a <- keep[idx[1]]
      b <- keep[idx[2]]
      drop_u <- if (a %in% exempt) b else if (b %in% exempt) a else
        if (relevance[[a]] >= relevance[[b]]) b else a
      keep <- setdiff(keep, drop_u)
    }
    return(keep)
  }

  # method == "vif"
  keep <- units
  repeat {
    non_ex <- setdiff(keep, exempt)
    if (length(keep) < 2 || !length(non_ex)) break
    u_of_col <- rep(keep, times = vapply(unit_cols[keep], length, integer(1)))
    Xs <- X[, unlist(unit_cols[keep], use.names = FALSE), drop = FALSE]
    Xs <- apply(Xs, 2, function(cl) {
      if (all(is.na(cl))) return(rep(0, length(cl)))
      cl[is.na(cl)] <- mean(cl, na.rm = TRUE)
      cl
    })
    v0 <- apply(Xs, 2, stats::var)
    ok <- is.finite(v0) & v0 > 0
    Xs <- Xs[, ok, drop = FALSE]
    u_of_col <- u_of_col[ok]
    if (ncol(Xs) < 2) break
    R <- stats::cor(Xs)
    inv <- tryCatch(solve(R), error = function(e) {
      tryCatch(solve(R + diag(1e-8, ncol(R))), error = function(e2) NULL)
    })
    if (is.null(inv)) {
      warning("netimpute: VIF computation failed (singular correlation matrix ",
              "even after ridging); stopping the collinearity prune early.",
              call. = FALSE)
      break
    }
    vif_col <- diag(inv)
    vif_unit <- tapply(vif_col, factor(u_of_col, levels = keep), max)
    vif_unit[is.na(vif_unit)] <- 0
    cand <- vif_unit[setdiff(names(vif_unit), exempt)]
    if (!length(cand) || max(cand) <= threshold) break
    worst <- names(cand)[which.max(cand)]
    keep <- setdiff(keep, worst)
  }
  keep
}

#' Per-target predictor selection for netmice() (network-aware quickpred)
#'
#' A network-aware, recursive generalization of \code{mice::quickpred()}:
#' for every target with missing values, select the predictors worth
#' entering that target's univariate imputation model, so \code{netmice()}
#' can fit lean, target-specific models instead of feeding everything (or
#' principal components of everything) into every model.
#'
#' @details
#' \strong{Screening (step 1).} A candidate qualifies as a direct predictor
#' of a target when its absolute pairwise-complete Pearson correlation with
#' the target's observed values, \emph{or} - like \code{mice::quickpred()} -
#' with the target's missingness indicator (\code{is.na(target)}), reaches
#' `mincor`. The missingness screen keeps variables that predict \emph{why}
#' values are missing, which supports the MAR assumption even when they
#' barely predict the value itself; disable with
#' \code{use_missingness = FALSE}. Categorical variables are screened
#' through their level dummies (a variable qualifies if any of its columns
#' does). For \emph{attribute} targets the candidates are all other
#' attributes plus every node-level network measure
#' (\code{\link{net_measures}} with `measure_set`, computed on the observed
#' ties with missing cells set to 0) - including the target's own
#' alter-composition homophily measures (e.g. the mean of the target among
#' a node's alters), which carry homophily/influence signal; only measures
#' that are direct functions of ego's own value of the target (its
#' \code{*_ei_index}, \code{*_diff_out}, \code{*_diff_in}) are withheld.
#' For \emph{network} targets, screening happens at the dyad level on
#' \code{\link{dyad_regression}}'s design terms (ego/alter/similarity terms
#' per attribute, cross-network tie/reciprocity/degree terms); the target's
#' own endogenous `reciprocity` and `twopath` terms are always included
#' without screening.
#'
#' \strong{Recursion (`steps`).} Missing values in a selected predictor are
#' themselves imputed inside the chained-equations loop, so the quality of
#' the target's model also depends on the predictors of its predictors.
#' With `steps = 3` (default), the direct predictors, their predictors, and
#' the predictors of those are all added to the target's model (marked with
#' their step in the result). The recursion only descends into predictors
#' that themselves have missing data - the predictors of a fully observed
#' predictor carry no information about the target beyond that predictor
#' itself - and, beyond step 1, adds attribute predictors only (for a
#' network target, a deeper-step attribute contributes its standard dyad
#' terms).
#'
#' \strong{Collinearity pruning.} Per target, near-duplicate predictors are
#' pruned so the model does not need a PCA fallback: `"pairwise"` (default)
#' repeatedly drops the less relevant member of any predictor pair whose
#' absolute correlation exceeds `collin_threshold` (default 0.9);
#' `"vif"` repeatedly drops the predictor with the largest
#' variance-inflation factor above `collin_threshold` (default 10),
#' catching multi-variable collinearity; `"none"` disables pruning. A
#' network target's `reciprocity`/`twopath` and cross-network
#' `*_tie`/`*_recip` terms are never pruned.
#'
#' \strong{Targets and dropping.} `targets` names the variables/networks
#' whose imputations you actually need. Every other variable or network
#' \emph{with missing data} that is not selected (at any step) as a
#' predictor for some target is listed in `$drop` - \code{\link{netmice}}
#' removes those from the imputation entirely, saving computation. Selected
#' predictors with missing data are kept and imputed too (using their own
#' direct predictors, restricted to the kept variables); fully observed
#' variables are never dropped.
#'
#' @param data A data.frame of node attributes, `NA` marking missing values.
#'   Every column is treated as an attribute (remove identifier columns
#'   first).
#' @param net_list Optional (preferably named) list of networks aligned to
#'   `data` (adjacency matrices with `NA` for unknown ties, as in
#'   \code{\link{netmice}}). When supplied, node-level network measures are
#'   candidate predictors for attribute targets, and networks with missing
#'   ties are themselves targets/candidates at the dyad level.
#' @param targets Character vector of attribute/network names to treat as
#'   the targets of the imputation, or `NULL` (default) for every variable
#'   and network with missing data.
#' @param mincor Minimum absolute correlation for the screening (default
#'   0.1, as in \code{mice::quickpred()}).
#' @param steps Number of predictor levels to add (default 3): 1 = direct
#'   predictors only, 2 = plus their predictors, etc.
#' @param collin_method Collinearity pruning method: `"pairwise"` (default),
#'   `"vif"`, or `"none"` - see Details.
#' @param collin_threshold Pruning threshold; `NULL` (default) means 0.9 for
#'   `"pairwise"` and 10 for `"vif"`.
#' @param use_missingness Logical (default `TRUE`); also screen candidates
#'   against each target's missingness indicator.
#' @param measure_set Which network measures to consider as candidate
#'   predictors for attribute targets, as in \code{\link{net_measures}}
#'   (default "core").
#' @param attr_types Optional named attribute-type overrides, as in
#'   \code{\link{netmice}}.
#' @param structural Optional structural-zero specification, as in
#'   \code{\link{netmice}}; structurally absent dyads are excluded from the
#'   dyad-level screening.
#'
#' @return An object of class `"netquickpred"`: a list with
#' \describe{
#'   \item{predictors}{Named list with one entry per variable/network that
#'     will be imputed. Attribute entries hold `predictors` (attribute
#'     names), `net_features` (node-level network measure names), `step`
#'     (the step at which each entered), and `relevance` (each predictor's
#'     max absolute screening correlation). Network entries hold
#'     `dyad_terms` (dyad-level term names, always including `reciprocity`
#'     and `twopath`), `step`, `relevance`, and the parent `attrs`/`nets`
#'     behind those terms.}
#'   \item{targets}{The resolved target names.}
#'   \item{keep, drop}{Attribute and network names kept for / dropped from
#'     the imputation (see Details).}
#'   \item{mincor, steps, collin_method, collin_threshold, use_missingness,
#'     measure_set}{The settings used.}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 60
#' friends <- matrix(rbinom(n * n, 1, 0.1), n, n); diag(friends) <- 0
#' attrs <- data.frame(age = rnorm(n, 35, 8))
#' attrs$performance <- attrs$age / 10 + rnorm(n)
#' attrs$happiness <- attrs$performance + rnorm(n)
#' attrs$happiness[sample(n, 12)] <- NA
#' attrs$performance[sample(n, 8)] <- NA
#' qp <- netquickpred(attrs, list(friends = friends), targets = "happiness")
#' qp
#' qp$predictors$happiness
netquickpred <- function(data,
                         net_list = NULL,
                         targets = NULL,
                         mincor = 0.1,
                         steps = 3,
                         collin_method = c("pairwise", "vif", "none"),
                         collin_threshold = NULL,
                         use_missingness = TRUE,
                         measure_set = "core",
                         attr_types = NULL,
                         structural = NULL) {
  collin_method <- match.arg(collin_method)
  if (is.null(collin_threshold)) {
    collin_threshold <- switch(collin_method, pairwise = 0.9, vif = 10,
                               none = Inf)
  }
  if (!is.numeric(mincor) || length(mincor) != 1 || mincor < 0 || mincor > 1) {
    stop("`mincor` must be a single value in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(steps) || length(steps) != 1 || steps < 1) {
    stop("`steps` must be a single integer >= 1.", call. = FALSE)
  }
  measure_set <- .resolve_measure_set(measure_set)

  data <- as.data.frame(data)
  n <- nrow(data)
  var_names <- names(data)
  types <- .resolve_attr_types(data, attr_types)

  has_nets <- !is.null(net_list) && length(net_list) > 0
  net_names <- character(0)
  mats <- list()
  struct_list <- NULL
  if (has_nets) {
    net_names <- names(net_list)
    if (is.null(net_names)) {
      net_names <- paste0("net", seq_along(net_list))
      names(net_list) <- net_names
    }
    mats <- lapply(net_list, .as_matrix_generic)
    for (k in seq_along(mats)) diag(mats[[k]]) <- 0
    if (!all(vapply(mats, nrow, integer(1)) == n)) {
      stop("Every network must have as many nodes as there are rows in `data`.",
           call. = FALSE)
    }
    struct_list <- .validate_structural(structural, net_names, n)
    mats <- .apply_structural_zeros(mats, struct_list)
  }

  vm <- var_names[vapply(data, anyNA, logical(1))]
  nm_miss <- net_names[vapply(mats, function(m) {
    any(is.na(m[row(m) != col(m)]))
  }, logical(1))]

  if (is.null(targets)) {
    tar <- c(vm, nm_miss)
  } else {
    unknown <- setdiff(targets, c(var_names, net_names))
    if (length(unknown)) {
      stop("`targets` references unknown variable/network(s): ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    no_miss <- setdiff(targets, c(vm, nm_miss))
    if (length(no_miss)) {
      message("netimpute: target(s) ", paste(no_miss, collapse = ", "),
              " have no missing values and are skipped.")
    }
    tar <- intersect(targets, c(vm, nm_miss))
    if (!length(tar)) {
      stop("None of the requested `targets` has missing values to impute.",
           call. = FALSE)
    }
  }
  if (!length(c(vm, nm_miss))) {
    stop("Nothing to select predictors for: no attribute or network has ",
         "missing values.", call. = FALSE)
  }

  ## ---- node-level screening matrix (attributes + network features) ----
  attr_exp <- .qp_expand_columns(data, types)
  feat_df <- NULL
  feat_exp <- list(mat = matrix(numeric(0), n, 0), parent = character(0))
  feat_src_attr <- character(0)
  feat_src_net <- character(0)
  multi <- length(net_names) > 1
  if (has_nets) {
    mats_zero <- lapply(mats, function(m) { m[is.na(m)] <- 0; m })
    gs <- lapply(mats_zero, .mat_to_igraph)
    feat_df <- .net_feature_frame(gs, data, measure_set, attr_types = types,
                                  net_names = net_names,
                                  clash_names = var_names)
    feat_exp <- .qp_expand_columns(feat_df)
    # generative source map: which attribute/network each feature column
    # derives from (structural measures have no source attribute)
    feat_src_attr <- stats::setNames(rep(NA_character_, ncol(feat_df)),
                                     names(feat_df))
    feat_src_net <- stats::setNames(rep(NA_character_, ncol(feat_df)),
                                    names(feat_df))
    for (k in seq_along(net_names)) {
      nk <- net_names[k]
      smn <- .qp_prefix(setdiff(measure_set, "homophily"), nk, multi, var_names)
      feat_src_net[intersect(smn, names(feat_df))] <- nk
      for (nm in var_names) {
        hn <- .qp_prefix(.qp_homophily_names(nm, types[[nm]]), nk, multi,
                         var_names)
        hit <- intersect(hn, names(feat_df))
        feat_src_net[hit] <- nk
        feat_src_attr[hit] <- nm
      }
    }
  }

  M <- cbind(attr_exp$mat, feat_exp$mat)
  parent <- c(attr_exp$parent, feat_exp$parent)
  # attribute names and (possibly prefixed) feature names cannot collide:
  # a feature whose bare name matches an attribute is prefixed by .qp_prefix
  unit_col_map <- split(seq_len(ncol(M)), factor(parent, levels = unique(parent)))
  feat_units <- if (is.null(feat_df)) character(0) else names(feat_df)

  C <- suppressWarnings(stats::cor(M, use = "pairwise.complete.obs"))
  C[!is.finite(C)] <- 0
  Cm <- NULL
  if (use_missingness && length(vm)) {
    miss_ind <- vapply(vm, function(v) as.numeric(is.na(data[[v]])), numeric(n))
    Cm <- suppressWarnings(stats::cor(M, miss_ind,
                                      use = "pairwise.complete.obs"))
    Cm[!is.finite(Cm)] <- 0
  }

  # the target's own ego-value-dependent homophily measures are withheld
  # when screening for that target (a variable must not be predicted by a
  # transform of itself); its alter-composition measures remain candidates
  ego_feat_names <- function(v) {
    if (!has_nets) return(character(0))
    unlist(lapply(net_names, function(nk) {
      .qp_prefix(paste0(v, "_", c("ei_index", "diff_out", "diff_in")),
                 nk, multi, var_names)
    }), use.names = FALSE)
  }

  node_cache <- new.env(parent = emptyenv())
  node_rel <- function(v) {
    if (exists(v, envir = node_cache)) return(get(v, envir = node_cache))
    tcols <- unit_col_map[[v]]
    rel_of <- function(cols) {
      r <- if (length(cols) && length(tcols)) max(abs(C[cols, tcols])) else 0
      if (!is.null(Cm) && v %in% colnames(Cm)) {
        r <- max(r, if (length(cols)) max(abs(Cm[cols, v])) else 0)
      }
      if (!is.finite(r)) 0 else r
    }
    a_units <- setdiff(var_names, v)
    f_units <- setdiff(feat_units, ego_feat_names(v))
    out <- list(
      attrs = vapply(a_units, function(a) rel_of(unit_col_map[[a]]),
                     numeric(1)),
      feats = vapply(f_units, function(f) rel_of(unit_col_map[[f]]),
                     numeric(1))
    )
    assign(v, out, envir = node_cache)
    out
  }

  ## ---- dyad-level screening (network targets/predictor networks) ----
  dyad_cache <- new.env(parent = emptyenv())
  dyad_screen <- function(g) {
    if (exists(g, envir = dyad_cache)) return(get(g, envir = dyad_cache))
    k <- match(g, net_names)
    built <- .build_dyad_data(mats, data, target_idx = k, attr_types = types,
                              other_net_predictors = "raw",
                              structural = if (is.null(struct_list)) NULL else
                                struct_list[[g]])
    d <- built$data
    ry <- !is.na(d$y)
    term_cols <- setdiff(names(d), c("i", "j", "y"))
    X <- as.matrix(d[term_cols])
    r1 <- abs(suppressWarnings(stats::cor(d$y, X,
                                          use = "pairwise.complete.obs")))
    r1[!is.finite(r1)] <- 0
    rel <- drop(r1)
    if (use_missingness) {
      r2 <- abs(suppressWarnings(stats::cor(as.numeric(!ry), X,
                                            use = "pairwise.complete.obs")))
      r2[!is.finite(r2)] <- 0
      rel <- pmax(rel, drop(r2))
    }
    names(rel) <- term_cols
    Cd <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
    Cd[!is.finite(Cd)] <- 0
    pmap <- .qp_dyad_parents(term_cols, data, types, net_names, g)
    out <- list(rel = rel, cor = Cd, X = X,
                parent = pmap$parent, kind = pmap$kind)
    assign(g, out, envir = dyad_cache)
    out
  }

  direct_attr_preds <- function(v) {
    r <- node_rel(v)
    names(r$attrs)[r$attrs >= mincor]
  }

  ## ---- per-target selection: screen, recurse, prune ----
  build_attr_entry <- function(t, allowed_attrs = NULL, allowed_feats = NULL,
                               depth = steps) {
    r <- node_rel(t)
    sel_attrs <- names(r$attrs)[r$attrs >= mincor]
    sel_feats <- names(r$feats)[r$feats >= mincor]
    if (!is.null(allowed_attrs)) sel_attrs <- intersect(sel_attrs, allowed_attrs)
    if (!is.null(allowed_feats)) sel_feats <- intersect(sel_feats, allowed_feats)
    step <- stats::setNames(rep(1L, length(sel_attrs) + length(sel_feats)),
                            c(sel_attrs, sel_feats))
    frontier <- sel_attrs
    s <- 2L
    while (s <= depth && length(frontier)) {
      nxt <- character(0)
      # only predictors that are themselves missing pull in their own
      # predictors: the predictors of a fully observed predictor carry no
      # information about the target beyond the predictor itself
      for (p in intersect(frontier, vm)) {
        np <- setdiff(direct_attr_preds(p), c(t, sel_attrs, nxt))
        if (!is.null(allowed_attrs)) np <- intersect(np, allowed_attrs)
        nxt <- c(nxt, np)
      }
      if (length(nxt)) {
        step[nxt] <- s
        sel_attrs <- c(sel_attrs, nxt)
      }
      frontier <- nxt
      s <- s + 1L
    }
    units <- c(sel_attrs, sel_feats)
    rel_all <- c(r$attrs, r$feats)[units]
    kept <- .qp_prune_units(units, unit_col_map, X = M, Cmat = C,
                            relevance = rel_all, exempt = character(0),
                            method = collin_method,
                            threshold = collin_threshold)
    list(kind = "attribute",
         predictors = kept[kept %in% sel_attrs],
         net_features = kept[kept %in% sel_feats],
         step = step[kept],
         relevance = rel_all[kept])
  }

  build_net_entry <- function(t, allowed_attrs = NULL, allowed_nets = NULL,
                              depth = steps) {
    scr <- dyad_screen(t)
    endo <- names(scr$kind)[scr$kind == "endo"]
    cand <- setdiff(names(scr$rel)[scr$rel >= mincor], endo)
    if (!is.null(allowed_attrs)) {
      cand <- cand[!(scr$kind[cand] == "attr" &
                       !(scr$parent[cand] %in% allowed_attrs))]
    }
    if (!is.null(allowed_nets)) {
      cand <- cand[!(scr$kind[cand] == "net" &
                       !(scr$parent[cand] %in% allowed_nets))]
    }
    sel_terms <- c(endo, cand)
    step <- stats::setNames(rep(1L, length(sel_terms)), sel_terms)
    sel_parents <- unique(scr$parent[cand][scr$kind[cand] == "attr"])
    frontier <- sel_parents
    s <- 2L
    while (s <= depth && length(frontier)) {
      nxt <- character(0)
      for (p in intersect(frontier, vm)) {
        np <- setdiff(direct_attr_preds(p), c(sel_parents, nxt))
        if (!is.null(allowed_attrs)) np <- intersect(np, allowed_attrs)
        nxt <- c(nxt, np)
      }
      if (length(nxt)) {
        # a deeper-step attribute contributes its standard dyad terms
        new_terms <- names(scr$parent)[scr$kind == "attr" &
                                         scr$parent %in% nxt]
        new_terms <- setdiff(new_terms, sel_terms)
        if (length(new_terms)) {
          step[new_terms] <- s
          sel_terms <- c(sel_terms, new_terms)
        }
        sel_parents <- c(sel_parents, nxt)
      }
      frontier <- nxt
      s <- s + 1L
    }
    # endogenous terms and cross-network cell terms keep their coefficients
    # unconditionally (see the tie-inflation safeguards in netmice)
    exempt <- c(endo,
                sel_terms[scr$kind[sel_terms] == "net" &
                            grepl("_(tie|recip)$", sel_terms)])
    term_cols <- stats::setNames(as.list(match(sel_terms, names(scr$rel))),
                                 sel_terms)
    kept <- .qp_prune_units(sel_terms, term_cols, X = scr$X, Cmat = scr$cor,
                            relevance = scr$rel[sel_terms], exempt = exempt,
                            method = collin_method,
                            threshold = collin_threshold)
    list(kind = "network",
         dyad_terms = kept,
         step = step[kept],
         relevance = scr$rel[kept],
         attrs = unique(scr$parent[kept][scr$kind[kept] == "attr"]),
         nets = unique(scr$parent[kept][scr$kind[kept] == "net"]))
  }

  preds <- list()
  for (t in tar) {
    preds[[t]] <- if (t %in% var_names) build_attr_entry(t) else
      build_net_entry(t)
  }

  ## ---- keep/drop: which missing variables/networks stay in the call ----
  sel_attrs_all <- unique(unlist(lapply(preds, function(e) {
    if (e$kind == "attribute") {
      c(e$predictors, feat_src_attr[e$net_features])
    } else {
      e$attrs
    }
  }), use.names = FALSE))
  sel_attrs_all <- sel_attrs_all[!is.na(sel_attrs_all)]
  sel_nets_all <- unique(unlist(lapply(preds, function(e) {
    if (e$kind == "attribute") feat_src_net[e$net_features] else e$nets
  }), use.names = FALSE))
  sel_nets_all <- sel_nets_all[!is.na(sel_nets_all)]

  imputed_attrs <- union(intersect(tar, var_names), intersect(sel_attrs_all, vm))
  imputed_nets <- union(intersect(tar, net_names), intersect(sel_nets_all, nm_miss))
  drop_attrs <- setdiff(vm, imputed_attrs)
  drop_nets <- setdiff(nm_miss, imputed_nets)
  kept_attrs <- setdiff(var_names, drop_attrs)
  kept_nets <- setdiff(net_names, drop_nets)

  # kept non-target variables are imputed too (they feed the targets'
  # models): give each its own direct-predictor model, restricted to the
  # kept universe so no model references a dropped column
  allowed_feats <- feat_units[feat_src_net[feat_units] %in% kept_nets &
                                (is.na(feat_src_attr[feat_units]) |
                                   feat_src_attr[feat_units] %in% kept_attrs)]
  for (w in setdiff(imputed_attrs, tar)) {
    preds[[w]] <- build_attr_entry(w, allowed_attrs = setdiff(kept_attrs, w),
                                   allowed_feats = allowed_feats, depth = 1L)
  }
  for (g in setdiff(imputed_nets, tar)) {
    preds[[g]] <- build_net_entry(g, allowed_attrs = kept_attrs,
                                  allowed_nets = kept_nets, depth = 1L)
  }

  structure(
    list(
      predictors = preds,
      targets = tar,
      keep = list(attributes = kept_attrs, networks = kept_nets),
      drop = list(attributes = drop_attrs, networks = drop_nets),
      mincor = mincor,
      steps = steps,
      collin_method = collin_method,
      collin_threshold = collin_threshold,
      use_missingness = use_missingness,
      measure_set = measure_set
    ),
    class = "netquickpred"
  )
}

#' @export
print.netquickpred <- function(x, ...) {
  cat("Class: netquickpred (per-target predictor selection)\n")
  cat("mincor:", x$mincor, " steps:", x$steps,
      " collinearity:", x$collin_method,
      if (is.finite(x$collin_threshold)) paste0("(threshold ", x$collin_threshold, ")"),
      " missingness screen:", x$use_missingness, "\n")
  cat("Targets:", paste(x$targets, collapse = ", "), "\n")
  if (length(x$drop$attributes) || length(x$drop$networks)) {
    cat("Dropped (missing data, not needed):",
        paste(c(x$drop$attributes, x$drop$networks), collapse = ", "), "\n")
  }
  for (nm in names(x$predictors)) {
    e <- x$predictors[[nm]]
    if (e$kind == "attribute") {
      cat("- ", nm, " (attribute): ", length(e$predictors), " attribute(s)",
          if (length(e$net_features)) paste0(" + ", length(e$net_features),
                                             " network feature(s)"),
          "\n", sep = "")
      if (length(e$step)) {
        cat("    ", paste0(names(e$step), " [step ", e$step, "]",
                           collapse = ", "), "\n", sep = "")
      }
    } else {
      cat("- ", nm, " (network): ", length(e$dyad_terms), " dyad term(s)\n",
          sep = "")
      if (length(e$dyad_terms)) {
        cat("    ", paste(e$dyad_terms, collapse = ", "), "\n", sep = "")
      }
    }
  }
  invisible(x)
}
