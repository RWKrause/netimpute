# A mice-style chained-equations loop that jointly imputes a nodal-attribute
# data.frame AND a list of networks. Within one iteration, every attribute
# and every network with missingness is visited *once*, in a single
# interleaved, randomly-ordered sequence (a different random order per
# imputation chain m) - not "all attributes then all networks" - so that
# whichever was updated most recently (whether that was an attribute or a
# network) always feeds into the next target's predictors. Network-derived
# predictors for an attribute (net_measures()) and
# attribute-derived predictors for a network (dyad_regression()'s terms) are
# therefore recomputed from scratch at every visit, never cached. The
# univariate imputation step is dispatched through .impute_univariate(),
# generically over the mice::mice.impute.<method> functions listed in
# .nm_methods - adding another numeric-response mice method is a one-line
# change there.

#' @keywords internal
.init_fill_vector <- function(x) {
  miss <- is.na(x)
  if (!any(miss)) return(x)
  obs <- x[!miss]
  if (!length(obs)) stop("A variable is entirely missing; cannot initialize imputation.", call. = FALSE)
  x[miss] <- sample(obs, sum(miss), replace = TRUE)
  x
}

#' @keywords internal
.init_fill_matrix <- function(mat, structural = NULL, init = c("zero", "sample")) {
  init <- match.arg(init)
  diag(mat) <- 0
  off <- row(mat) != col(mat)
  # structurally absent cells are fixed zeros: they are neither filled nor
  # allowed into the donor pool (their design-zeros would inflate the zeros
  # sampled into the genuinely missing cells)
  if (!is.null(structural)) off <- off & !structural
  miss <- off & is.na(mat)
  if (any(miss)) {
    if (init == "zero") {
      # conservative start: every missing tie begins as "no tie", so the
      # first sweep's endogenous predictors (reciprocity, two-paths, other
      # networks' ties/degrees) are built from observed ties only, instead
      # of from random fills that can seed a self-reinforcing surplus of
      # imputed ties across iterations
      mat[miss] <- 0
    } else {
      obs_vals <- mat[off & !is.na(mat)]
      if (!length(obs_vals)) stop("A network is entirely missing; cannot initialize imputation.", call. = FALSE)
      mat[miss] <- sample(obs_vals, sum(miss), replace = TRUE)
    }
  }
  mat
}

#' @keywords internal
.mat_to_igraph <- function(mat) {
  directed <- !isSymmetric(unname(mat))
  igraph::graph_from_adjacency_matrix(
    mat, mode = if (directed) "directed" else "undirected",
    weighted = if (any(mat[mat != 0] != 1)) TRUE else NULL, diag = FALSE
  )
}

#' Network-level diagnostics tracked across netmice() iterations
#'
#' Density, (network-level) reciprocity, global transitivity, isolate count,
#' and average inverse geodesic distance (global efficiency). The last is
#' computed as the mean of 1/d(i,j) over all ordered pairs i != j; since
#' `1/Inf == 0` in R, disconnected pairs (including isolates) contribute 0
#' rather than `NaN`/`Inf` - only the (excluded) diagonal would otherwise
#' cause a division by zero.
#'
#' @param mat An adjacency matrix (weights, if any, are binarized: `!= 0`).
#' @return A named list: density, reciprocity, transitivity, n_isolates,
#'   avg_inv_geodesic.
#' @export
#'
#' @examples
#' set.seed(1)
#' m <- matrix(rbinom(400, 1, 0.15), 20, 20)
#' diag(m) <- 0
#' net_diagnostics(m)
net_diagnostics <- function(mat) {
  b <- (mat != 0) * 1
  diag(b) <- 0
  directed <- !isSymmetric(unname(b))
  g <- igraph::graph_from_adjacency_matrix(b, mode = if (directed) "directed" else "undirected", diag = FALSE)
  d <- igraph::distances(g, mode = "all")
  off <- row(d) != col(d)
  inv <- 1 / d[off]
  inv[is.infinite(inv)] <- 0
  tr <- igraph::transitivity(g, type = "global")
  list(
    density          = igraph::edge_density(g, loops = FALSE),
    reciprocity      = if (directed) igraph::reciprocity(g) else 1,
    transitivity     = if (is.nan(tr)) NA_real_ else tr,
    n_isolates       = sum(igraph::degree(g, mode = "all") == 0),
    avg_inv_geodesic = mean(inv)
  )
}

#' Univariate imputation methods netmice() accepts, all dispatched to the
#' corresponding `mice::mice.impute.<method>()` function. Restricted to
#' methods whose working model handles a *numeric* response, since netmice()
#' passes numeric targets (continuous attributes, integer-coded binary
#' attributes, tie values) - factor-only methods like "logreg"/"polr" are
#' excluded ("polyreg" is applied automatically to multinomial attributes
#' and is not part of this user-selectable set). To support another mice
#' method, add its name here (and, if it has donor-type tuning arguments
#' beyond `donors`, extend the dispatch in `.impute_univariate()`).
#' @keywords internal
.nm_methods <- c("pmm", "midastouch", "sample", "cart", "rf",
                 "norm", "norm.nob", "norm.boot", "norm.predict", "mean")

#' @keywords internal
.validate_method <- function(method) {
  if (!is.character(method) || length(method) != 1 || is.na(method)) {
    stop("`method` must be a single character string.", call. = FALSE)
  }
  if (!method %in% .nm_methods) {
    stop("Unsupported imputation method: '", method, "'. Supported methods ",
         "(all dispatched to mice::mice.impute.<method>): ",
         paste(.nm_methods, collapse = ", "), ".", call. = FALSE)
  }
  method
}

#' Resolve `method` into a per-target method map, mice-style
#'
#' Accepts a single unnamed string (one method for everything - the historic
#' interface), or a named character vector assigning methods to individual
#' attributes/networks: named entries override, and at most one unnamed
#' entry sets the default for everything without its own entry ("pmm" when
#' no unnamed entry is given). E.g. `c("cart", performance = "norm")` uses
#' `norm` for performance and `cart` elsewhere. Returns a full named vector
#' with one (validated) entry per attribute and network.
#' @keywords internal
.resolve_methods <- function(method, var_names, net_names) {
  if (!is.character(method) || length(method) < 1 || anyNA(method)) {
    stop("`method` must be a character vector (a single method, or a named ",
         "vector of per-target methods).", call. = FALSE)
  }
  nms <- names(method)
  if (is.null(nms)) nms <- rep("", length(method))
  unnamed <- method[!nzchar(nms)]
  named <- method[nzchar(nms)]
  if (length(unnamed) > 1) {
    stop("`method` may contain at most one unnamed entry (the default ",
         "method); name every other entry after an attribute or network.",
         call. = FALSE)
  }
  default <- if (length(unnamed)) unname(unnamed) else "pmm"
  for (m in c(default, unname(named))) .validate_method(m)
  unknown <- setdiff(names(named), c(var_names, net_names))
  if (length(unknown)) {
    stop("`method` names unknown variable/network(s): ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- stats::setNames(rep(default, length(var_names) + length(net_names)),
                         c(var_names, net_names))
  out[names(named)] <- named
  out
}

#' Dispatch one univariate imputation to `mice::mice.impute.<method>()`
#'
#' Generic over every method in `.nm_methods`: the mice univariate methods
#' share the signature `(y, ry, x, ...)`, so dispatch needs no per-method
#' code; `donors` is forwarded only to methods that declare it (currently
#' `pmm`). Method-specific tuning parameters (e.g. `cart`'s `minbucket`,
#' `rf`'s `ntree`) keep mice's defaults. `"rf"` additionally requires the
#' \pkg{ranger} or \pkg{randomForest} package (checked by mice itself).
#' `"polyreg"` is reserved for the internal multinomial-attribute path and
#' is not part of the user-selectable set.
#' @keywords internal
.impute_univariate <- function(method, y, ry, x, donors) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required for netmice() imputation. Install with install.packages('mice').",
         call. = FALSE)
  }
  if (is.null(dim(x)) || ncol(x) == 0) {
    x <- matrix(1, length(y), 1, dimnames = list(NULL, "(Intercept)"))
  }
  if (method == "polyreg") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Package 'nnet' is required for multinomial imputation of attributes with more than ",
           "two categories. Install with install.packages('nnet').", call. = FALSE)
    }
    return(mice::mice.impute.polyreg(y = y, ry = ry, x = x))
  }
  .validate_method(method)
  # some mice internals (e.g. the logging on a singular fit) look these up
  # dynamically in the calling frame
  state <- list(it = 0L, im = 0L, dep = "y", meth = method, log = FALSE)
  loggedEvents <- NULL
  fn <- get(paste0("mice.impute.", method), envir = asNamespace("mice"))
  args <- list(y = y, ry = ry, x = x)
  if ("donors" %in% names(formals(fn))) args$donors <- donors
  do.call(fn, args)
}

#' PMM for network ties with lme4 random intercepts
#'
#' Fits a linear mixed working model on the observed dyads - fixed effects
#' for every predictor column plus random intercepts per ego (sender), alter
#' (receiver), and/or unordered dyad, as requested - then performs classic
#' predictive-mean matching on the conditional fitted values (fixed +
#' predicted random effects): each missing dyad receives the observed `y` of
#' one of its `donors` nearest neighbours in fitted-value space, drawn at
#' random. A linear working model is used regardless of tie type, consistent
#' with the PMM step for binary ties (the fitted values only rank/match
#' donors). Unlike `mice::mice.impute.pmm()` this matching is "type-0" (no
#' Bayesian draw of the coefficients); between-imputation variability comes
#' from the donor draw and the chain stochasticity. If the mixed model fails
#' to fit, falls back to `mice::mice.impute.pmm()` with a warning.
#' @keywords internal
.impute_pmm_ranef <- function(y, ry, x, donors, ego, alter, random_intercepts) {
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for `net_random_intercepts`. ",
         "Install with install.packages('lme4').", call. = FALSE)
  }
  x <- as.matrix(x)
  if (ncol(x)) colnames(x) <- paste0("V", seq_len(ncol(x)))
  df <- as.data.frame(x)
  df$y <- y
  df <- .add_dyad_groups(cbind(df, i = ego, j = alter))
  re_terms <- c(ego = "(1 | .ego)", alter = "(1 | .alter)",
                dyad = "(1 | .dyad)")[random_intercepts]
  fml <- stats::reformulate(c(colnames(x), re_terms), response = "y")

  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      lme4::lmer(fml, data = df[ry, , drop = FALSE], REML = FALSE)
    )),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    warning("netimpute: mixed-model PMM failed to fit (",
            conditionMessage(fit),
            "); falling back to standard PMM for this visit.", call. = FALSE)
    return(mice::mice.impute.pmm(y = y, ry = ry, x = x, donors = donors))
  }

  # allow.new.levels covers nodes/dyads whose every tie is missing: they get
  # the population-level (zero) random effect, i.e. the fixed-effect fit.
  yhat <- stats::predict(fit, newdata = df, allow.new.levels = TRUE)
  yhat_obs <- yhat[ry]
  y_obs <- y[ry]
  vapply(yhat[!ry], function(yh) {
    nn <- order(abs(yhat_obs - yh))[seq_len(min(donors, length(yhat_obs)))]
    y_obs[nn][sample.int(length(nn), 1)]
  }, numeric(1))
}

#' Decompose interaction predictor columns that involve an endogenous
#' statistic, for the sequential (Gibbs) tie updater
#'
#' A `models` formula can interact the target network's endogenous terms
#' (`reciprocity`, `twopath`) with any other dyad-level column, producing
#' `model.matrix()` columns named e.g. `age_ego:reciprocity`. Inside the
#' sequential sweep those columns must track the live matrix state exactly
#' like the pure `reciprocity`/`twopath` main effects, so each is split
#' here into its static part (the product of its non-endogenous factors,
#' recomputed from the dyad data `d`) and flags for which endogenous
#' statistics it carries. Every dyad-level column is numeric, so an
#' interaction column is an elementwise product and the split is exact. A
#' `:`-named column whose non-endogenous factors are not all columns of `d`
#' cannot be decomposed and is left alone (it stays static during the
#' sweep, like any other predictor).
#'
#' Returns a list, one entry per decomposable endogenous interaction
#' column: `col` (its index in the predictor matrix), `static` (per-row
#' static product), `recip`/`twop` (logical flags).
#' @keywords internal
.gibbs_endo_interactions <- function(xnames, d) {
  out <- list()
  if (is.null(xnames)) return(out)
  for (k in seq_along(xnames)) {
    nm <- xnames[k]
    if (!grepl(":", nm, fixed = TRUE)) next
    parts <- strsplit(nm, ":", fixed = TRUE)[[1]]
    endo <- parts %in% c("reciprocity", "twopath")
    if (!any(endo)) next
    static_parts <- parts[!endo]
    if (!all(static_parts %in% names(d))) next
    static <- if (length(static_parts)) {
      Reduce(`*`, lapply(static_parts, function(p) as.numeric(d[[p]])))
    } else {
      rep(1, nrow(d))
    }
    # mirror .clean_predictor_matrix()'s NA handling so the recomputed
    # static product matches the fitted column
    if (anyNA(static)) static[is.na(static)] <- mean(static, na.rm = TRUE)
    out[[nm]] <- list(col = k, static = static,
                      recip = "reciprocity" %in% parts,
                      twop = "twopath" %in% parts)
  }
  out
}

#' Sequential ("Gibbs") single-tie imputation for one network visit
#'
#' The default tie updater (the simultaneous PMM update is the alternative,
#' `net_update = "simultaneous"`). A
#' working model is fit once per visit on the observed dyads (mice's
#' estimation rule: the target's own imputed cells never enter the fit,
#' while the predictor columns - attribute terms, other-network terms, and
#' the endogenous statistics - are evaluated from the full current, partly
#' imputed state) - logistic
#' regression for a binary network, linear regression otherwise - and one
#' coefficient vector is drawn from its asymptotic posterior (mirroring
#' mice's Bayesian beta-draw, so between-imputation variability reflects
#' estimation uncertainty). The missing cells are then visited ONE AT A TIME
#' in random order; each cell's linear predictor is evaluated with the
#' CURRENT values of the endogenous statistics - reciprocity (the transpose
#' cell, read off the matrix, O(1)) and the two-path indicator (from a
#' running two-path count matrix) - which are refreshed after every single
#' draw via change statistics: toggling cell (i, j) only alters two-path
#' counts of pairs running through i or j, i.e. two O(n) vector updates,
#' never a full O(n^3) recomputation. A binary tie is drawn from
#' Bernoulli(plogis(eta)); a weighted tie by PMM - the `donors` observed
#' dyads with the closest fitted values donate one of their observed values
#' at random (donor fitted values are taken at fit time and not re-refreshed
#' during the sweep).
#'
#' Each draw therefore conditions on all previously imputed ties: this is
#' exactly the full-conditional (Gibbs) update of an ERGM whose parameters
#' were estimated by pseudo-likelihood, restricted to the missing cells
#' (observed ties are never resampled). It removes the synchronous-update
#' artifact of the simultaneous scheme, where every missing cell jumps
#' together off the same stale snapshot of the statistics. Note it does not,
#' by itself, rule out ERGM-style degeneracy - that is a property of the
#' statistics; the bounded 0/1 `twopath` indicator is the main safeguard.
#'
#' All static (non-endogenous) predictor columns - attribute terms,
#' other-network terms, any `models` extras, possibly PCA components - are
#' evaluated once at visit start and held fixed during the sweep; only
#' `reciprocity` and `twopath` are updated per draw - including any `models`
#' interaction column that involves them (e.g. `age_ego:reciprocity`), whose
#' static factors are recomputed from `d` and whose endogenous factors are
#' read off the live matrix state per cell, exactly like the main effects.
#'
#' `check = TRUE` recomputes the binarized adjacency and the two-path count
#' matrix from scratch at the end and errors if the incremental bookkeeping
#' diverged - used by the unit tests, skipped in production calls.
#' @keywords internal
.impute_ties_gibbs <- function(d, ry, x, mat, binary, donors, check = FALSE) {
  endo_int <- .gibbs_endo_interactions(colnames(x), d)
  x <- as.matrix(x)
  recip_col <- match("reciprocity", colnames(x))
  twop_col  <- match("twopath", colnames(x))

  # syntactic stand-in names: predictor columns may be PCA components or
  # model.matrix terms with characters that formulas cannot carry verbatim
  vnames <- if (ncol(x)) paste0("V", seq_len(ncol(x))) else character(0)
  df <- stats::setNames(as.data.frame(x), vnames)
  df$y <- d$y
  fml <- if (length(vnames)) stats::reformulate(vnames, response = "y") else
    stats::as.formula("y ~ 1")
  fit <- if (binary) {
    suppressWarnings(stats::glm(fml, data = df[ry, , drop = FALSE],
                                family = stats::binomial()))
  } else {
    stats::lm(fml, data = df[ry, , drop = FALSE])
  }
  b <- stats::coef(fit)
  ok <- !is.na(b)
  Vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (!is.null(Vc)) {
    Rch <- tryCatch(chol(Vc), error = function(e) NULL)
    if (!is.null(Rch)) {
      b[ok] <- b[ok] + drop(crossprod(Rch, stats::rnorm(nrow(Rch))))
    }
  }
  b[!ok] <- 0  # aliased (rank-deficient) terms contribute nothing
  intercept <- b[["(Intercept)"]]
  beta_x <- b[vnames]

  eta_all <- intercept + drop(x %*% beta_x)
  b_recip <- if (!is.na(recip_col)) beta_x[[recip_col]] else 0
  b_twop  <- if (!is.na(twop_col))  beta_x[[twop_col]]  else 0
  # strip the endogenous contributions once; they are re-added per cell
  # from the live matrix state inside the sweep
  eta_base <- eta_all
  if (!is.na(recip_col)) eta_base <- eta_base - b_recip * x[, recip_col]
  if (!is.na(twop_col))  eta_base <- eta_base - b_twop  * x[, twop_col]
  for (ei in endo_int) {
    eta_base <- eta_base - beta_x[[ei$col]] * x[, ei$col]
  }

  y_obs <- d$y[ry]
  yhat_obs <- if (!binary) eta_all[ry] else NULL

  B <- (mat != 0) * 1
  Tp <- B %*% B  # two-path counts; maintained incrementally below

  for (r in sample(which(!ry))) {
    i <- d$i[r]
    j <- d$j[r]
    eta <- eta_base[r] + b_recip * mat[j, i] +
      b_twop * as.numeric(Tp[i, j] > 0)
    for (ei in endo_int) {
      eta <- eta + beta_x[[ei$col]] * ei$static[r] *
        (if (ei$recip) mat[j, i] else 1) *
        (if (ei$twop) as.numeric(Tp[i, j] > 0) else 1)
    }
    v <- if (binary) {
      stats::rbinom(1, 1, stats::plogis(eta))
    } else {
      nn <- order(abs(yhat_obs - eta))[seq_len(min(donors, length(yhat_obs)))]
      y_obs[nn][sample.int(length(nn), 1)]
    }
    if (v != mat[i, j]) {
      mat[i, j] <- v
      db <- as.numeric(v != 0) - B[i, j]
      if (db != 0) {
        # change statistics: B[i,j] enters Tp[a,b] = sum_k B[a,k] B[k,b]
        # only as B[i,j]B[j,b] (row i) and B[a,i]B[i,j] (column j); the
        # updates read row j / column i of B, which do not contain B[i,j],
        # so their order is immaterial
        Tp[i, ] <- Tp[i, ] + db * B[j, ]
        Tp[, j] <- Tp[, j] + db * B[, i]
        B[i, j] <- B[i, j] + db
      }
    }
  }

  if (check) {
    stopifnot(isTRUE(all.equal(B, (mat != 0) * 1, check.attributes = FALSE)),
              isTRUE(all.equal(Tp, B %*% B, check.attributes = FALSE)))
  }
  mat
}

#' Clean an auto-generated predictor matrix for PMM, and - if `max_cols` is
#' given and exceeded - collapse it to `max_cols` principal components.
#'
#' This exists because the auto-generated predictor sets can plausibly
#' approach or exceed the number of *observed* cases available to fit on
#' (e.g. two networks' worth of core measures plus attribute homophily terms
#' can be 30+ columns, which is not safely smaller than a few dozen complete
#' cases). mice's own ridge fallback inside estimice() only nudges a
#' near-singular X'X - it is not wrapped in its own try(), and does not
#' reliably rescue a design matrix that is *exactly* rank-deficient (p >= n),
#' which then crashes with an uncaught "system is exactly singular" error
#' instead of degrading gracefully. Preemptively reducing dimensionality
#' here keeps p comfortably below n by construction, the same fix already
#' offered explicitly via net_predictors(output = "pca") - just applied
#' automatically inside the imputation loop where the user has no direct
#' control over the predictor count.
#'
#' `ry` (the observed-response indicator) matters because the model is fit
#' on the observed rows only: a column can vary over all rows yet be
#' constant on the observed subset (e.g. another network's ties that occur
#' only on dyads where the target is missing). Such a column adds nothing to
#' the fit, and when it is constant *zero* mice's ridge penalty
#' (ridge * diag(X'X)) is zero for it, so solve() crashes with "system is
#' exactly singular" instead of being rescued - hence constant-on-observed
#' columns are dropped here. The PCA branch needs the analogous guard:
#' components beyond the input matrix's true rank have (numerically) zero
#' scores everywhere and crash the same way, so components are capped at the
#' effective rank and re-checked against `ry`.
#'
#' `keep_raw` names columns that must never be absorbed into the PCA: they
#' are exempted from the reduction and appended unchanged (still subject to
#' the NA-fill and constant-column screens). Used for a network target's own
#' endogenous terms (`reciprocity`, `twopath`) and the cross-network
#' dyad terms (`*_tie`, `*_recip`), whose named coefficients carry
#' substantive dependence structure that a composite component would dilute.
#' @keywords internal
.clean_predictor_matrix <- function(x, max_cols = NULL, ry = NULL,
                                    keep_raw = NULL) {
  x <- as.matrix(x)
  # an empty selection (predictor_selection with no survivors) is legal:
  # .impute_univariate() falls back to an intercept-only model
  if (ncol(x) == 0) return(matrix(numeric(0), nrow(x), 0))
  x <- apply(x, 2, function(col) {
    if (all(is.na(col))) return(rep(0, length(col)))
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  })
  keep_varying <- function(m) {
    var_all <- apply(m, 2, stats::var)
    ok <- is.finite(var_all) & var_all > 0
    if (!is.null(ry) && any(ry)) {
      var_obs <- apply(m[ry, , drop = FALSE], 2, stats::var)
      ok <- ok & is.finite(var_obs) & var_obs > 0
    }
    m[, ok, drop = FALSE]
  }
  x <- keep_varying(x)

  if (!is.null(max_cols) && max_cols >= 1 && ncol(x) > max_cols) {
    keep_idx <- if (is.null(keep_raw) || is.null(colnames(x))) {
      rep(FALSE, ncol(x))
    } else {
      colnames(x) %in% keep_raw
    }
    kept <- x[, keep_idx, drop = FALSE]
    rest <- x[, !keep_idx, drop = FALSE]
    if (ncol(rest) > max_cols) {
      pca <- stats::prcomp(rest, center = TRUE, scale. = TRUE)
      n_comp <- min(max_cols, sum(pca$sdev > pca$sdev[1] * 1e-8))
      rest <- keep_varying(pca$x[, seq_len(n_comp), drop = FALSE])
    }
    x <- cbind(kept, rest)
  }
  x
}

#' @keywords internal
.safe_max_cols <- function(n_obs) {
  max(5, floor(n_obs / 3))
}

#' Evaluate one target visit, rethrowing any error with the target's name
#' and position in the chain. printFlag only prints the per-iteration visit
#' order, which does not identify *which* target a mid-sweep failure came
#' from - this does.
#' @keywords internal
.visit_context <- function(label, expr) {
  tryCatch(expr, error = function(e) {
    stop("netimpute: imputation failed at ", label, ": ",
         conditionMessage(e), call. = FALSE)
  })
}

#' Attributes/networks a set of `models` formulas depends on
#'
#' Every variable referenced in a formula is mapped back to the attribute(s)
#' and network(s) it is derived from, so those can be protected from being
#' dropped via `targets` (the formula must stay evaluable at every visit).
#' Bare names are taken as-is; internally created derived names are matched
#' by their naming convention: a leading `<name>_` or an embedded `_<name>_`
#' marks a dependency on that attribute/network (e.g. `advice_tie` depends
#' on the network `advice`; `friends_age_alter_mean_in` depends on the
#' network `friends` AND the attribute `age`). Matching is deliberately
#' permissive - protecting too much only keeps a variable imputed, whereas
#' protecting too little would error mid-chain on an unevaluable formula.
#' @keywords internal
.model_referenced_bases <- function(model_map, var_names, net_names) {
  universe <- c(var_names, net_names)
  vars <- unique(unlist(lapply(model_map, all.vars)))
  bases <- intersect(vars, universe)
  for (v in setdiff(vars, universe)) {
    for (nm in universe) {
      if (startsWith(v, paste0(nm, "_")) ||
          grepl(paste0("_", nm, "_"), v, fixed = TRUE)) {
        bases <- c(bases, nm)
      }
    }
  }
  unique(bases)
}

#' Parse a `models` list into a named-by-LHS lookup of formulas
#' @keywords internal
.parse_models <- function(models) {
  if (is.null(models) || length(models) == 0) return(list())
  out <- list()
  for (mstr in models) {
    f <- stats::as.formula(mstr)
    lhs <- all.vars(f)[1]
    out[[lhs]] <- f
  }
  out
}

#' Build extra (formula-specified) predictor columns, NA-safe and
#' intercept-free, for cbind-ing onto the auto-generated predictor matrix.
#' @keywords internal
.model_extra_terms <- function(formula, data, ry = NULL) {
  rhs <- stats::delete.response(stats::terms(formula))
  mf <- tryCatch(
    stats::model.frame(rhs, data = data, na.action = stats::na.pass),
    error = function(e) {
      stop("netimpute: could not evaluate a `models` formula (",
           conditionMessage(e), "). ",
           "Every predictor named in a formula must exist in `data`, or - for network-derived ",
           "terms in an attribute's formula, or for any term in a network's formula - match the ",
           "auto-generated predictor names (e.g. 'friends_indegree' for an attribute's formula ",
           "when multiple networks are supplied, or 'age_absdiff'/'friends_tie' for a network's ",
           "own dyad-level formula). See ?netmice.", call. = FALSE)
    }
  )
  mm <- stats::model.matrix(rhs, data = mf)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  .clean_predictor_matrix(mm, ry = ry)
}

#' Run one full imputation chain (one of the `m` imputations)
#' @keywords internal
.run_one_chain <- function(im,
                           data,
                           mats0,
                           var_names,
                           net_names,
                           attr_types,
                           var_levels,
                           ry_vars,
                           ry_nets,
                           structural,
                           net_dependence,
                           vm_names,
                           net_missing_names,
                           method,
                           donors,
                           measure_set,
                           other_net_predictors,
                           n_components,
                           net_random_intercepts,
                           model_map,
                           qp_pred,
                           feat_clash_names,
                           net_init,
                           net_update,
                           net_binary,
                           maxit,
                           printFlag,
                           seed) {
  if (!is.na(seed)) set.seed(seed + im, kind = "Mersenne-Twister")

  cur_data <- data
  for (v in vm_names) {
    cur_data[[v]] <- .init_fill_vector(cur_data[[v]])
    }
  cur_mats <- mats0
  for (nm in net_missing_names) {
    cur_mats[[nm]] <- .init_fill_matrix(cur_mats[[nm]], structural[[nm]],
                                        init = net_init)
  }
  # random init fills may violate the network-dependence rules; re-zero the
  # determined cells so every chain starts from a consistent state
  cur_mats <- .enforce_dependence(cur_mats, net_dependence, ry_nets)$mats

  targets <- sample(c(vm_names, net_missing_names))

  # igraph views of the current networks, for attribute-visit predictors.
  # Converted once here and refreshed per network right after that network
  # is re-imputed - not rebuilt for every network at every attribute visit.
  cur_gs <- stats::setNames(lapply(cur_mats, .mat_to_igraph), net_names)

  net_diag_names <- c("density", "reciprocity", "transitivity", "n_isolates", "avg_inv_geodesic")
  chainMean <- matrix(NA_real_, length(vm_names), maxit,
                      dimnames = list(vm_names, seq_len(maxit)))
  chainVar  <- chainMean
  netChain  <- array(NA_real_, dim = c(length(net_names),
                                       length(net_diag_names), maxit),
                      dimnames = list(net_names, net_diag_names,
                                      seq_len(maxit)))

  for (it in seq_len(maxit)) {
    if (printFlag) {
      cat(" imputation", im, "- iter", it, "- order:",
          paste(targets, collapse = " -> "), "\n")
    }
    for (tgt in targets) .visit_context(
      sprintf("target '%s' (iteration %d, imputation %d)", tgt, it, im), {
      if (tgt %in% vm_names) {
        v <- tgt
        ry <- ry_vars[[v]]

        # ALL attributes - including `v` itself - feed the homophily-type
        # network features: alter-composition measures of v (e.g. the mean of
        # v among one's alters, and especially among in-neighbours, whose
        # nominations are typically observed reports by others) are functions
        # of OTHER nodes' values and carry real homophily/influence signal
        # for v. Only the measures that are direct functions of ego's own
        # current value of v (v's E-I index and avg-|diff| terms) are dropped
        # for the target: they would regress v on a transform of itself,
        # anchoring imputations to their own previous draws. The raw `v`
        # column is likewise excluded from the predictors - it is the target.
        other_data <- cur_data[setdiff(var_names, v)]
        ego_feats <- paste0(v, "_", c("ei_index", "diff_out", "diff_in"))
        # feature naming (network-name prefixing) lives in
        # .net_feature_frame(); `feat_clash_names` pins the prefix-on-clash
        # rule to the ORIGINAL attribute set, so names stay identical to the
        # ones netquickpred() selected even after `targets` dropped columns
        net_feats_df <- .net_feature_frame(cur_gs, cur_data, measure_set,
                                           attr_types, net_names,
                                           exclude = ego_feats,
                                           clash_names = feat_clash_names)
        sel <- if (!is.null(qp_pred)) qp_pred[[v]] else NULL
        if (is.null(sel)) {
          # net_feats_df is NULL when every network was dropped via `targets`
          auto_x <- .clean_predictor_matrix(
            cbind(.prep_pca_matrix(other_data),
                  if (!is.null(net_feats_df)) .prep_pca_matrix(net_feats_df)),
            max_cols = .safe_max_cols(sum(ry)),
            ry = ry
          )
        } else {
          # netquickpred() selection: only this target's selected attributes
          # and network features enter the model. The PCA cap stays on as a
          # rank fallback only (selection normally keeps p well below it).
          a_cols <- intersect(sel$predictors, names(other_data))
          f_cols <- intersect(sel$net_features, names(net_feats_df))
          parts <- list()
          if (length(a_cols)) {
            parts <- c(parts, list(.prep_pca_matrix(other_data[a_cols])))
          }
          if (length(f_cols)) {
            parts <- c(parts, list(.prep_pca_matrix(net_feats_df[f_cols])))
          }
          raw_x <- if (length(parts)) do.call(cbind, parts) else
            matrix(numeric(0), nrow(cur_data), 0)
          auto_x <- .clean_predictor_matrix(raw_x,
                                            max_cols = .safe_max_cols(sum(ry)),
                                            ry = ry)
        }

        if (!is.null(model_map[[v]])) {
          combined_df <- cbind(other_data, net_feats_df)
          extra <- .model_extra_terms(model_map[[v]], combined_df, ry = ry)
          x <- cbind(auto_x, extra)
          x <- x[, !duplicated(colnames(x)), drop = FALSE]
        } else {
          x <- auto_x
        }

        lvls <- var_levels[[v]]
        if (is.numeric(cur_data[[v]])) {
          imp_vals <- .impute_univariate(method[[v]],
                                         y = cur_data[[v]],
                                         ry = ry,
                                         x = x,
                                         donors = donors)
        } else if (length(lvls) > 2) {
          imp_vals <- .impute_univariate("polyreg",
                                         y = cur_data[[v]],
                                         ry = ry,
                                         x = x,
                                         donors = donors)
        } else {
          y_codes <- as.integer(factor(cur_data[[v]], levels = lvls))
          imp_codes <- .impute_univariate(method[[v]],
                                          y = y_codes,
                                          ry = ry,
                                          x = x,
                                          donors = donors)
          # donor-based methods return observed codes, but regression-type
          # methods (e.g. "norm") can return values outside {1, 2}: round
          # and clamp to the valid level range before mapping back
          imp_codes <- pmin(pmax(round(as.numeric(imp_codes)), 1),
                            length(lvls))
          imp_vals <- lvls[imp_codes]
        }
        cur_data[[v]][!ry] <- imp_vals

      } else {
        k <- match(tgt, net_names)
        # cells determined 0 by the network-dependence rules given the
        # partners' CURRENT state are treated exactly like structural cells
        # for this visit: fixed at 0, not fit on, never imputed. The mask is
        # re-evaluated here every visit because the partners are themselves
        # re-imputed as the chain progresses.
        dep_mask <- .dependence_zero_mask(cur_mats, net_dependence, tgt)
        struct_eff <- structural[[tgt]]
        if (!is.null(dep_mask)) {
          stale <- dep_mask & !ry_nets[[tgt]] & cur_mats[[tgt]] != 0
          if (any(stale)) cur_mats[[tgt]][stale] <- 0
          struct_eff <- if (is.null(struct_eff)) dep_mask else
            (struct_eff | dep_mask)
        }
        # structural cells of the target never enter the dyad data: they are
        # not fit on (no zero inflation) and never receive imputed values.
        # With a netquickpred() selection, the "other networks" terms are
        # built raw and then subset by name - the selection replaces the PCA
        # reduction as the dimensionality control.
        sel <- if (!is.null(qp_pred)) qp_pred[[tgt]] else NULL
        built <- .build_dyad_data(cur_mats,
                                  cur_data,
                                  target_idx = k,
                                  attr_types = attr_types,
                                  other_net_predictors = if (is.null(sel))
                                    other_net_predictors else "raw",
                                  n_components = n_components,
                                  structural = struct_eff)
        d <- built$data
        ry <- ry_nets[[tgt]][cbind(d$i, d$j)]
        # With a dyad random intercept the fixed `reciprocity` term (the
        # transpose of y) makes the working model exactly singular (see
        # dyad_regression); the dyad effect is the SRM reciprocity term and
        # replaces it.
        drop_cols <- c("i", "j", "y",
                       if ("dyad" %in% net_random_intercepts) "reciprocity")
        # the target's endogenous terms and the cross-network dyad terms
        # keep their own coefficients even when the dimensionality
        # safeguard collapses the rest to principal components
        keep_raw <- c("reciprocity", "twopath",
                      paste0(rep(setdiff(net_names, tgt), each = 2),
                             c("_tie", "_recip")))
        if (is.null(sel)) {
          auto_x <- .clean_predictor_matrix(d[setdiff(names(d), drop_cols)],
                                             max_cols = .safe_max_cols(sum(ry)),
                                             ry = ry,
                                             keep_raw = keep_raw)
        } else {
          want <- setdiff(unique(c("reciprocity", "twopath", sel$dyad_terms)),
                          drop_cols)
          auto_x <- .clean_predictor_matrix(d[intersect(names(d), want)],
                                             max_cols = .safe_max_cols(sum(ry)),
                                             ry = ry,
                                             keep_raw = keep_raw)
        }

        if (!is.null(model_map[[tgt]])) {
          extra <- .model_extra_terms(model_map[[tgt]], d, ry = ry)
          x <- cbind(auto_x, extra)
          x <- x[, !duplicated(colnames(x)), drop = FALSE]
        } else {
          x <- auto_x
        }

        y <- d$y
        if (net_update == "gibbs") {
          cur_mats[[tgt]] <- .impute_ties_gibbs(d = d,
                                                ry = ry,
                                                x = x,
                                                mat = cur_mats[[tgt]],
                                                binary = net_binary[[tgt]],
                                                donors = donors)
        } else {
          if (length(net_random_intercepts)) {
            imp_vals <- .impute_pmm_ranef(y = y,
                                          ry = ry,
                                          x = x,
                                          donors = donors,
                                          ego = d$i,
                                          alter = d$j,
                                          random_intercepts = net_random_intercepts)
          } else {
            imp_vals <- .impute_univariate(method[[tgt]],
                                           y = y,
                                           ry = ry,
                                           x = x,
                                           donors = donors)
          }
          mis_idx <- which(!ry)
          cur_mats[[tgt]][cbind(d$i[mis_idx], d$j[mis_idx])] <- imp_vals
        }
        cur_gs[[tgt]] <- .mat_to_igraph(cur_mats[[tgt]])
        # freshly imputed ties may newly determine cells of dependent
        # networks; propagate so the current state always satisfies the rules
        if (!is.null(net_dependence)) {
          enf <- .enforce_dependence(cur_mats, net_dependence, ry_nets)
          cur_mats <- enf$mats
          for (nm2 in enf$changed) {
            cur_gs[[nm2]] <- .mat_to_igraph(cur_mats[[nm2]])
          }
        }
      }
    })

    ## ---- diagnostics (once per full sweep) ----
    for (v in vm_names) {
      vals <- cur_data[[v]][!ry_vars[[v]]]
      if (!is.numeric(vals)) vals <- as.numeric(
        factor(vals, levels = var_levels[[v]]))
      chainMean[v, it] <- mean(vals, na.rm = TRUE)
      chainVar[v, it]  <- stats::var(vals, na.rm = TRUE)
    }
    for (nm in net_names) {
      dk <- net_diagnostics(cur_mats[[nm]])
      netChain[nm, , it] <- unlist(dk[net_diag_names])
    }
  }

  list(data = cur_data,
       net_list = cur_mats,
       chainMean = chainMean,
       chainVar = chainVar,
       netChain = netChain,
       visit_order = targets)
}

#' Joint multiple imputation of nodal attributes and network ties
#'
#' A mice-style chained-equations algorithm (Van Buuren's sampler
#' architecture, generalized to also cycle through network ties) that
#' imputes a nodal-attribute data.frame and a list of networks together.
#' Within each iteration, every attribute and every network with
#' missingness is visited exactly once, in a single interleaved sequence
#' whose order is randomized independently for each of the `m` imputation
#' chains (see Details). Each attribute is imputed conditional on the other
#' attributes and on network-derived predictors from whichever networks were
#' most recently updated (via \code{\link{net_measures}} with the requested
#' `measure_set`) - including the alter-composition
#' homophily measures of the target attribute itself (e.g. the mean of the
#' attribute among a node's alters, overall and over incoming/outgoing ties
#' separately), which are functions of other nodes' values and typically
#' among the strongest predictors under homophily or social influence. Only
#' the target's homophily measures that are direct functions of ego's own
#' value (its E-I index and average-absolute-difference terms) are withheld
#' when imputing that attribute, since regressing a variable on a transform
#' of itself would anchor imputations to their own previous draws. Each
#' network's missing ties are imputed
#' conditional on whichever attributes and other networks were most recently
#' updated (via \code{\link{dyad_regression}}'s predictor set). The
#' univariate draw dispatches to \code{mice::mice.impute.<method>()}:
#' predictive mean matching (`method = "pmm"`, the default) or any of the
#' other supported numeric-response mice methods (`"midastouch"`,
#' `"sample"`, `"cart"`, `"rf"`, the `"norm"` family, `"mean"`) - see the
#' `method` argument for exactly where the method applies and where a
#' different working model is fixed by design.
#'
#' @details
#' \strong{Interleaving and re-computation.} Because network-derived
#' attribute predictors (e.g. an alter's average attribute value) and
#' attribute-derived network predictors (e.g. ego/alter values) both depend
#' on *the other* data type, they cannot be computed once and reused: they
#' are rebuilt from scratch at every single visit in the sequence, using
#' whatever the current state of attributes/networks happens to be at that
#' point. A variable that depends on both a network and an attribute is
#' therefore effectively recomputed twice per full iteration - once when the
#' network it depends on is (re-)imputed, and once when the attribute it
#' depends on is (re-)imputed - and, because the visit order is randomized,
#' not always in the same order.
#'
#' \strong{Estimation rule (identical to mice).} At every visit, the
#' working model is estimated on the entries where the \emph{target} is
#' observed - the target's own current imputations never enter its fit, at
#' any sweep - while every predictor column is evaluated from the full
#' current state of the data: the current imputations of all \emph{other}
#' attributes and networks, and, for a network target, the endogenous
#' statistics (`reciprocity`, `twopath`) computed from its own partly
#' imputed matrix. Previously imputed information therefore informs every
#' model at every sweep through the predictors, exactly as in
#' \pkg{mice}'s fully conditional specification - netimpute simply
#' generalizes which predictors exist, not how the models are estimated.
#'
#' \strong{Automatic dimensionality safeguard.} The auto-generated predictor
#' set (structural + homophily measures across every network, for attribute
#' imputation; ego/alter/reciprocity/etc. terms, for network imputation) can
#' plausibly approach or exceed the number of observed cases available to
#' impute a given target from, especially with several networks and a modest
#' number of nodes. Rather than risk an exactly rank-deficient design matrix
#' (which crashes mice's internal ridge fallback instead of degrading
#' gracefully), the auto-generated predictors are capped at
#' \code{max(5, floor(n_observed / 3))} columns, reducing to that many
#' principal components when exceeded - the same fix
#' \code{\link{net_predictors}(output = "pca")} offers explicitly, applied
#' automatically here since the user has no direct handle on the predictor
#' count inside the imputation loop. This cap applies only to the
#' auto-generated predictors; any extra terms from `models` are added on top
#' and are the user's own, deliberate choice.
#'
#' \strong{Per-target predictor selection (`predictor_selection`,
#' `targets`).} As an alternative to the everything-plus-PCA default,
#' `predictor_selection = "quickpred"` runs \code{\link{netquickpred}} once
#' per call: each imputed attribute gets only its correlation-selected
#' attributes and network measures, each imputed network only its selected
#' dyad-level terms (its own `reciprocity`/`twopath` always included), so
#' every coefficient keeps its own name instead of being folded into a
#' principal component. The selection is computed once, on the pre-chain
#' data (missing ties treated as 0 for the network measures), and reused at
#' every visit of every chain; the PCA cap remains only as a rank fallback.
#' `targets` additionally drops unneeded missing-data variables/networks
#' from the call - see the argument documentation.
#'
#' \strong{Custom imputation models (`models`).} By default every attribute
#' and network is imputed from the full auto-generated predictor set (other
#' attributes/networks' structural+homophily measures, or dyad-level ego/
#' alter/reciprocity/etc. terms). To protect a specific hypothesized
#' relationship (e.g. `happiness ~ indegree + age + performance * gender`)
#' from being diluted/omitted, supply it via `models`. The left-hand side
#' identifies which attribute or network the formula applies to; the
#' right-hand side is evaluated with `model.matrix()` (so `*`/`:` interaction
#' syntax works) and its columns are appended to - never substituted for -
#' the auto-generated predictors, exactly satisfying "regardless of the
#' formula, all other predictors should still be added". Any attribute or
#' network without an entry in `models` uses the standard, fully
#' auto-generated predictor set. For a network's own formula, right-hand-side
#' terms must be dyad-level predictor names as produced by
#' \code{\link{dyad_regression}} (e.g. `age_absdiff`, `friends_tie`), not raw
#' node attributes. For an attribute's formula referencing a network-derived
#' term, use the bare measure name (e.g. `indegree`) if only one network is
#' supplied, or the network-prefixed name (e.g. `friends_indegree`) if more
#' than one network is supplied (prefixing is also applied, even with a
#' single network, if the bare name would otherwise collide with a raw
#' attribute of the same name) - matching \code{\link{net_measures_core}}'s
#' naming.
#'
#' \strong{Interactions in `models`.} Because the right-hand side goes
#' through `model.matrix()`, `*`/`:` interaction syntax works for network
#' targets exactly as it does for attributes, and any internally created
#' dyad-level term can be interacted with any other - including the target
#' network's own endogenous terms and other networks' terms. For example,
#' `friends ~ age_ego:reciprocity` lets the tendency to reciprocate depend
#' on the sender's age; `friends ~ advice_tie:gender_same` interacts another
#' network's tie with a same-category indicator;
#' `friends ~ reciprocity:twopath` interacts the two endogenous terms with
#' each other. Main effects the auto-generated predictor set already
#' contains are not duplicated. Under `net_update = "gibbs"`, an interaction
#' involving `reciprocity` or `twopath` is - like those main effects -
#' re-evaluated from the current matrix state after every single tie draw,
#' so it stays live throughout the sequential sweep. The section
#' \emph{Dyad-level term names} below lists every name a network's formula
#' can reference; the section \emph{Implemented network measures} lists
#' every measure name an attribute's formula can reference.
#'
#' @section Dyad-level term names (for a network's `models` formula):
#' When a network is the imputation target, the working model's rows are its
#' off-diagonal cells (i = sender, j = receiver) and the following predictor
#' columns are created internally (by \code{\link{dyad_regression}}'s
#' builder). A `models` formula for a network may reference any of them,
#' alone or in interactions:
#' \describe{
#'   \item{\code{reciprocity}}{the target's transpose cell - the tie value
#'     from j back to i.}
#'   \item{\code{twopath}}{0/1 indicator for at least one two-path
#'     i -> k -> j in the target network (at least one shared contact).}
#'   \item{\code{<attr>_ego}, \code{<attr>_alter}, \code{<attr>_absdiff}}{
#'     for every \emph{continuous} attribute: the sender's value, the
#'     receiver's value, and their absolute difference.}
#'   \item{\code{<attr>_ego_<level>}, \code{<attr>_alter_<level>},
#'     \code{<attr>_same}}{for every \emph{binary/multinomial} attribute:
#'     sender and receiver dummies for each non-reference level, and a 0/1
#'     same-category indicator.}
#'   \item{\code{<net>_tie}, \code{<net>_recip}, \code{<net>_ego_outdeg},
#'     \code{<net>_ego_indeg}, \code{<net>_alter_indeg},
#'     \code{<net>_alter_outdeg}}{for every \emph{other} network: its cell
#'     (i, j), its transpose cell (j, i), and the sender's and receiver's
#'     out- and in-degrees in it. With `other_net_predictors = "pca"` the
#'     four degree terms are collapsed into components named
#'     \code{other_net_PC<k>}; the \code{_tie}/\code{_recip} terms always
#'     stay raw.}
#' }
#'
#' @section Implemented network measures (for `measure_set` and an attribute's `models` formula):
#' These are the node-level measures \code{\link{net_measures}} implements;
#' any of them can be requested via `measure_set` (individually or through
#' the named sets) and - once part of `measure_set` - referenced in an
#' attribute's `models` formula (bare, or network-prefixed as described
#' above):
#' \describe{
#'   \item{Named sets}{\code{"core"} = outdegree, indegree,
#'     reciprocity_ratio, bonacich_power, betweenness, isolate, constraint,
#'     harmonic_closeness, local_clustering, homophily;
#'     \code{"full"} = every measure below plus homophily.}
#'   \item{Degree family}{\code{total_degree}, \code{indegree},
#'     \code{outdegree}, \code{reciprocal_degree},
#'     \code{nonreciprocal_degree}, \code{reciprocity_ratio},
#'     \code{weighted_degree}}
#'   \item{Influence}{\code{eigenvector}, \code{bonacich_power},
#'     \code{alpha_centrality}, \code{pagerank}, \code{gilschmidt}
#'     (requires \pkg{sna})}
#'   \item{Brokerage/position}{\code{betweenness}, \code{load_centrality},
#'     \code{stress_centrality}, \code{flow_betweenness} (the latter three
#'     require \pkg{sna}), \code{constraint}, \code{effective_size},
#'     \code{efficiency}, \code{redundancy}}
#'   \item{Distance-based}{\code{closeness}, \code{harmonic_closeness},
#'     \code{information_centrality} (requires \pkg{sna})}
#'   \item{Prestige}{\code{prestige} (requires \pkg{sna})}
#'   \item{Position/closure}{\code{neighbor_degree},
#'     \code{local_clustering}, \code{coreness}, \code{isolate}}
#'   \item{\code{"homophily"} (attribute-based block)}{per
#'     binary/multinomial attribute: \code{<attr>_ei_index},
#'     \code{<attr>_blau_index}, \code{<attr>_alter_mode}; per continuous
#'     attribute: \code{<attr>_diff_out}, \code{<attr>_diff_in},
#'     \code{<attr>_alter_mean}, \code{<attr>_alter_mean_in},
#'     \code{<attr>_alter_mean_out}, \code{<attr>_alter_min},
#'     \code{<attr>_alter_max}. When imputing an attribute, its own
#'     ego-value-dependent measures (\code{<attr>_ei_index},
#'     \code{<attr>_diff_out}, \code{<attr>_diff_in}) are withheld from its
#'     predictor set (see above) and cannot be referenced in its formula.}
#' }
#'
#' \strong{Parallelization (`ncores`).} When `ncores > 1`, the `m` chains run
#' as independent \pkg{future} futures (`future::multisession`), each with
#' its own random visit order and RNG stream (`seed + im`, reproducible
#' regardless of `ncores`). \strong{This requires the package to be
#' installed} (`devtools::install()`/`R CMD INSTALL`), not merely loaded via
#' `devtools::load_all()`, because multisession workers are separate R
#' processes that resolve `netimpute`'s internal functions by loading the
#' installed namespace. Use `ncores = 1` (default) during development.
#'
#' @param data A data.frame/tibble of node attributes, `NA` marking missing
#'   values, rows aligned to the node order of every network in `net_list`.
#' @param net_list A (preferably named) list of networks, all with the same
#'   number of nodes and node order as `data`. Supply networks with unknown
#'   ties as adjacency matrices with `NA` marking those cells (igraph/`network`
#'   objects are treated as fully observed, since neither format can encode
#'   an "unknown" tie).
#' @param m Number of independent imputations (default 5, as in `mice`).
#' @param maxit Number of chained-equations iterations per imputation
#'   (default 20 - higher than `mice`'s own default of 5, since jointly
#'   cycling through both attributes and networks is a slower-mixing,
#'   more complex conditional model than attribute-only imputation).
#' @param method Univariate imputation method(s), dispatched to the
#'   corresponding \code{mice::mice.impute.<method>()} function. Each
#'   method is one of `"pmm"` (default), `"midastouch"`, `"sample"`,
#'   `"cart"`, `"rf"` (requires \pkg{ranger} or \pkg{randomForest}),
#'   `"norm"`, `"norm.nob"`, `"norm.boot"`, `"norm.predict"`, or `"mean"` -
#'   the mice univariate methods whose working model handles a numeric
#'   response. A single unnamed string applies one method globally. For
#'   per-target methods, supply a named character vector, mice-style: named
#'   entries assign a method to individual attributes/networks, and at most
#'   one unnamed entry sets the default for every target without its own
#'   entry (`"pmm"` if omitted). E.g.
#'   `method = c("cart", performance = "norm", friends = "pmm")` imputes
#'   `performance` with `norm`, the `friends` network with PMM, and
#'   everything else with trees. Names must be attribute columns or
#'   networks; entries for targets whose working model is fixed by design
#'   (see below) are ignored with a message. Method-specific
#'   tuning parameters (e.g. `cart`'s `minbucket`, `rf`'s `ntree`) use
#'   mice's defaults. The method applies to numeric attributes, to binary
#'   attributes (integer-coded; a regression-type method's continuous draws
#'   are rounded and clamped back to the two codes - donor- and tree-based
#'   methods return valid codes as-is), and to network ties under
#'   `net_update = "simultaneous"`. It does NOT apply where another working
#'   model is fixed by design: nominal attributes with more than two
#'   categories are *always* imputed with a proper multinomial logistic
#'   regression (`mice::mice.impute.polyreg()`, via `nnet::multinom()`) - a
#'   numeric-response method would have to impose an arbitrary ordering on
#'   the categories; ties under the default `net_update = "gibbs"` are drawn
#'   from the sequential updater's own logistic/linear working model; and
#'   ties under `net_random_intercepts` use the lme4 mixed-model PMM step.
#'   Note that for *binary networks* under the simultaneous scheme, a
#'   regression-type method (the `norm` family, `"mean"`) can impute
#'   non-0/1 tie values; prefer a donor/tree method - or the default
#'   tie-wise updater, which draws proper Bernoulli ties - there. `"mean"`
#'   and `"norm.predict"` are deterministic and yield improper multiple
#'   imputations (no between-imputation variability from the draw itself);
#'   mice's documentation carries the same caveat.
#' @param donors Number of donors for donor-based matching (used by
#'   `method = "pmm"`, by the tie-wise updater's weighted-tie PMM step, and
#'   by the `net_random_intercepts` mixed-model PMM step; default 5).
#' @param measure_set Which network measures supply predictors for attribute
#'   imputation, passed to \code{\link{net_measures}}: any combination of
#'   "core" (default, faster), "full", "homophily", and individual measure
#'   names - named sets and single measures are unioned, e.g.
#'   \code{c("core", "total_degree")}. The \emph{Implemented network
#'   measures} section below lists every available measure name. Note the
#'   "homophily" block (part of
#'   both "core" and "full") supplies the alter-composition predictors that
#'   carry homophily/influence signal; dropping it by selecting structural
#'   measures only will typically weaken attribute imputation.
#' @param attr_types Optional named attribute-type overrides.
#' @param other_net_predictors,n_components Passed to
#'   \code{\link{dyad_regression}}'s dyad-data builder for the "other
#'   networks" terms - "raw" (all terms) or "pca" (dyad-level `*_tie`/
#'   `*_recip` cross-network terms kept raw, node-level degree terms
#'   reduced to the first `n_components` components). Neither mode ever
#'   PCA-transforms the target network's own `reciprocity`/`twopath`
#'   terms or the cross-network cell terms: those keep their own
#'   coefficients even when the automatic dimensionality safeguard (see
#'   Details) collapses the remaining predictors to components.
#' @param net_random_intercepts `NULL` (default) or a character vector - any
#'   of `"ego"`, `"alter"`, `"dyad"`. When set, the working model for
#'   *network-tie* imputation is a linear mixed model fit with
#'   \code{lme4::lmer()} instead of mice's linear PMM regression, adding a
#'   random intercept per sender node (`"ego"`), per receiver node
#'   (`"alter"`), and/or per unordered node pair (`"dyad"`); missing ties are
#'   then imputed by predictive-mean matching on the conditional fitted
#'   values (fixed + predicted random effects). Ego/alter intercepts capture
#'   actor-level activity/popularity heterogeneity beyond the degree-based
#'   fixed effects, as in the social relations model. The `"dyad"` intercept
#'   is the SRM relationship effect and is meaningful only for *directed*
#'   networks, where cells (i,j) and (j,i) are two genuine observations per
#'   pair; for an undirected network those two rows are exact duplicates, so
#'   a dyad intercept degenerates into a residual term - a warning is issued
#'   if any network to impute is undirected. When `"dyad"` is included, the
#'   target's fixed `reciprocity` term is dropped from the working model
#'   (it is the transpose of `y` and exactly singular with a per-dyad
#'   intercept - the dyad effect is the SRM reciprocity term and replaces
#'   it). Note the matching is "type-0"
#'   (no Bayesian coefficient draw, unlike `mice::mice.impute.pmm()`);
#'   between-imputation variability comes from the donor draw and the chain
#'   stochasticity. Attribute imputation is unaffected. Requires \pkg{lme4}.
#'   Random intercepts exist only in the simultaneous update scheme: setting
#'   this argument while leaving `net_update` at its default switches the
#'   tie updating to `"simultaneous"` (with a message); combining it with an
#'   explicit `net_update = "gibbs"` is an error.
#' @param structural `NULL` (default), one n x n logical matrix, or a named
#'   list of n x n logical matrices whose names match networks in `net_list`.
#'   `TRUE` marks a cell whose tie is *structurally absent* - zero by design
#'   (e.g. a nomination that was impossible by the study design) rather than
#'   a genuine "no tie" observation. A single matrix applies the same
#'   structural cells to every network; a named list lets them differ per
#'   network (networks without an entry have none). Structural cells must
#'   hold 0 or `NA` in the supplied networks (an observed non-zero value
#'   there errors); they are fixed at 0 throughout: never imputed, excluded
#'   from the initialization donor pool, and removed from the dyad-level
#'   regression rows whenever that network is the imputation target, so the
#'   design-zeros cannot inflate the zeros of the working model.
#' @param net_dependence `NULL` (default) or a named list with elements
#'   `necessary` and/or `forbidden`, each a list of length-2 character
#'   vectors of network names - *adaptive* structural constraints between
#'   networks:
#'   \describe{
#'     \item{necessary}{`c("A", "B")` means a tie in A is required for a tie
#'       in B ("if not A then not B"): B's ties can only exist where A
#'       currently has a tie. Order matters - the first name is the
#'       prerequisite (parent) network.}
#'     \item{forbidden}{`c("A", "B")` means ties in A and B are mutually
#'       exclusive ("if A then not B, if B then not A"). Order is
#'       irrelevant.}
#'   }
#'   E.g. \code{list(necessary = list(c("interaction", "friendship")),
#'   forbidden = list(c("friend_pos", "friend_neg")))}. Three things happen:
#'   (1) before the chain, missing cells that are logically determined by
#'   *observed* cells are deduced and filled (e.g. an observed tie in B
#'   implies the missing tie in A is 1; an observed tie in A implies the
#'   missing cell in a forbidden partner is 0) - deduced cells count as
#'   observed and are never imputed; observed data that contradict a rule
#'   raise an error. (2) whenever a constrained network is the imputation
#'   target, its cells currently determined to be 0 by the partners' current
#'   (possibly imputed) state are fixed at 0 and excluded from the working
#'   model, exactly like `structural` cells - but re-evaluated at every
#'   visit, since the partners are themselves re-imputed. (3) after every
#'   network re-imputation the rules are propagated (non-observed cells of
#'   dependent networks re-zeroed), so every completed network list
#'   satisfies every rule. Deduction of *required* ties (necessary rule) is
#'   only possible when the parent network is binary; for a weighted parent
#'   those cells stay missing, with a warning.
#' @param models Optional list/character vector of formula strings (or
#'   formula objects), one per attribute/network you want a custom model
#'   for - see Details. E.g.
#'   `list("happiness ~ indegree + age + performance * gender")`. Formulas
#'   may interact any of the internally created terms, for networks
#'   included - e.g. `list("friends ~ age_ego:reciprocity")`; see the
#'   \emph{Interactions in `models`} paragraph and the \emph{Dyad-level
#'   term names} section.
#' @param predictor_selection How each target's predictor set is chosen.
#'   `"all"` (default; the behavior to date): every auto-generated predictor
#'   enters every model, with the automatic PCA safeguard capping the count.
#'   `"quickpred"`: \code{\link{netquickpred}} is run once, before the
#'   chains, and each imputed attribute/network gets its own lean predictor
#'   set - correlation-screened (against the target's values and, like
#'   \code{mice::quickpred()}, its missingness indicator), recursively
#'   expanded over `steps` levels of predictors-of-predictors, and pruned
#'   for collinearity - instead of everything-plus-PCA. Any extra `models`
#'   terms are still appended on top. Alternatively, pass a precomputed
#'   `netquickpred` object (built on the same data/networks) for full
#'   control and inspection; targets/variables without an entry in it fall
#'   back to the full predictor set.
#' @param targets Optional character vector naming the attributes/networks
#'   whose imputations you actually need. Every *other* variable or network
#'   with missing data that is not selected as a predictor (at any step) by
#'   the predictor selection is dropped from the call entirely - it is not
#'   imputed, and it appears in neither the predictors nor the returned
#'   completed data - saving computation. Selected predictors with missing
#'   data are kept and imputed alongside the targets; fully observed
#'   variables are never dropped, and variables/networks a `models` formula
#'   depends on - by bare name or as the base of a derived term such as
#'   `advice_tie` or `friends_age_alter_mean` - or networks bound by
#'   `net_dependence` rules are never dropped either.
#'   With `predictor_selection = "all"`, non-target missing-data variables
#'   are dropped unconditionally (there is no selection to keep them).
#'   Default `NULL`: everything with missing data is imputed.
#' @param mincor,steps,collin_method,collin_threshold Passed to
#'   \code{\link{netquickpred}} when `predictor_selection = "quickpred"`:
#'   the screening correlation threshold (default 0.1), the number of
#'   predictor levels (default 3 - direct predictors, their predictors, and
#'   theirs again), the collinearity pruning method (`"pairwise"`, `"vif"`,
#'   or `"none"`), and its threshold (`NULL` = 0.9 for pairwise, 10 for
#'   VIF). Ignored otherwise.
#' @param net_update How a network's missing cells are updated at each
#'   visit. `"gibbs"` (default): sequential single-tie ("tie-wise")
#'   imputation - a working model is fit once per visit on the observed
#'   dyads (logistic regression for a binary network, linear regression
#'   otherwise), one coefficient vector is drawn from its asymptotic
#'   posterior, and the missing cells are then visited one at a time in
#'   random order, each draw using the *current* values of the endogenous
#'   statistics (`reciprocity`, `twopath`), which are refreshed after every
#'   single draw via change statistics (O(1) for reciprocity, two O(n)
#'   vector updates on the running two-path count matrix - never a full
#'   recomputation). Binary ties are drawn from `Bernoulli(plogis(eta))`;
#'   weighted ties by PMM donor matching on the fitted values. Each tie is
#'   thus drawn conditional on all previously imputed ties - the
#'   full-conditional (Gibbs) update of an ERGM-style model with
#'   pseudo-likelihood-estimated parameters, restricted to the missing
#'   cells - which yields properly dependent within-network imputations and
#'   avoids the synchronous-update artifact of the simultaneous scheme.
#'   Static predictors (attribute terms, other-network terms, `models`
#'   extras) are evaluated once per visit and held fixed during the sweep;
#'   `models` interaction terms involving `reciprocity`/`twopath` are
#'   refreshed per draw like the main effects. `"simultaneous"`: one PMM
#'   working model per visit imputes every missing cell at once, all from
#'   the same snapshot of the endogenous statistics. `net_random_intercepts`
#'   requires the simultaneous scheme: combining it with an *explicit*
#'   `net_update = "gibbs"` errors, while setting `net_random_intercepts`
#'   alone (leaving `net_update` at its default) falls back to
#'   `"simultaneous"` with a message. Attribute imputation is unaffected by
#'   this argument.
#' @param net_init How each chain initializes a network's missing cells
#'   before the first sweep. `"zero"` (default) starts every missing tie at 0
#'   ("no tie"), so the first sweep's endogenous predictors (reciprocity,
#'   two-paths, other networks' tie/degree terms) are built from observed
#'   ties only - a deliberately conservative start that biases the very
#'   first sweep's structural terms downward but cannot seed a
#'   self-reinforcing surplus of random initial ties. `"sample"` restores
#'   the previous behavior of filling missing cells by resampling observed
#'   off-diagonal values (marginal-density start).
#' @param ncores Number of parallel workers for the `m` chains via
#'   \pkg{future} (default 1 = sequential). See Details for the
#'   installed-package requirement when `ncores > 1`.
#' @param seed Optional RNG seed (chain `im` uses `seed + im`, so results are
#'   identical whether run sequentially or in parallel).
#' @param printFlag Logical; print progress (imputation/iteration/visit
#'   order) as in `mice`.
#'
#' @return An object of class `"netmids"` with elements: `data`, `net_list`
#'   (original, NA-preserving inputs, except that structurally absent cells
#'   are fixed at 0 and cells deduced from `net_dependence` are filled),
#'   `m`, `maxit`, `method` (the resolved per-target method map: one named
#'   entry per attribute/network), `donors`, `net_init`, `net_update`,
#'   `structural` (the
#'   per-network list of structural-zero matrices, or `NULL`),
#'   `net_dependence` (the normalized rule list, or `NULL`),
#'   `models`, `predictor_selection` (`"all"`, or the `netquickpred` object
#'   used - inspect it to see each target's selected predictors), `targets`
#'   (as supplied; when variables were dropped via `targets`, `data`,
#'   `net_list`, and the completed datasets contain only the kept
#'   columns/networks), `imp` (list of `m` completed attribute data.frames),
#'   `imp_nets` (list of `m` lists of completed adjacency matrices),
#'   `visit_orders` (the randomized target order used in each chain),
#'   `chainMean`/`chainVar` (variable x iteration x m arrays of the mean/
#'   variance of the *imputed* values only, as in `mice`), and `netChain`
#'   (network x diagnostic x iteration x m array from
#'   \code{\link{net_diagnostics}}). Use \code{\link{complete_netmice}} to
#'   extract a completed dataset.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 30
#' friends <- matrix(rbinom(n * n, 1, 0.15), n, n)
#' diag(friends) <- 0
#' # unknown ties are marked NA (off the diagonal)
#' off <- which(row(friends) != col(friends))
#' friends[sample(off, 40)] <- NA
#'
#' attrs <- data.frame(age = rnorm(n, 35, 8),
#'                     performance = rnorm(n),
#'                     gender = sample(c("F", "M"), n, TRUE))
#' attrs$performance[sample(n, 4)] <- NA
#'
#' # a few seconds: loading the imputation engine dominates the runtime
#' \donttest{
#' fit <- netmice(attrs, list(friends = friends), m = 2, maxit = 3,
#'                seed = 1, printFlag = FALSE,
#'                models = list("performance ~ indegree + age * gender"))
#' fit
#' completed <- complete_netmice(fit, 1)
#' anyNA(completed$data)
#' }
netmice <- function(data,
                    net_list,
                    m = 5, maxit = 20,
                    method = "pmm",
                    donors = 5,
                    measure_set = "core",
                    attr_types = NULL,
                    other_net_predictors = c("raw", "pca"),
                    n_components = 3,
                    net_random_intercepts = NULL,
                    structural = NULL,
                    net_dependence = NULL,
                    models = NULL,
                    predictor_selection = "all",
                    targets = NULL,
                    mincor = 0.1,
                    steps = 3,
                    collin_method = c("pairwise", "vif", "none"),
                    collin_threshold = NULL,
                    net_init = c("zero", "sample"),
                    net_update = c("gibbs", "simultaneous"),
                    ncores = 1L,
                    seed = NA,
                    printFlag = TRUE) {

  net_init <- match.arg(net_init)
  net_update_missing <- missing(net_update)
  net_update <- match.arg(net_update)
  measure_set <- .resolve_measure_set(measure_set)
  other_net_predictors <- match.arg(other_net_predictors)
  if (!is.null(net_random_intercepts)) {
    net_random_intercepts <- match.arg(net_random_intercepts,
                                       c("ego", "alter", "dyad"),
                                       several.ok = TRUE)
  }
  if (net_update == "gibbs" && length(net_random_intercepts)) {
    if (net_update_missing) {
      # the sequential sampler has its own logistic/linear working model, so
      # random intercepts only exist in the simultaneous PMM scheme: an
      # explicit `net_random_intercepts` beats the implicit gibbs default
      message("netimpute: `net_random_intercepts` is set - falling back to ",
              "net_update = 'simultaneous' (the default sequential Gibbs ",
              "updater does not support random intercepts).")
      net_update <- "simultaneous"
    } else {
      stop("net_update = 'gibbs' is not compatible with `net_random_intercepts`: ",
           "the sequential sampler uses its own logistic/linear working model. ",
           "Use one or the other.", call. = FALSE)
    }
  }

  data <- as.data.frame(data)
  n <- nrow(data)
  var_names <- names(data)

  net_names <- names(net_list)
  if (is.null(net_names)) {
    net_names <- paste0("net", seq_along(net_list))
    names(net_list) <- net_names
  }

  name_clash <- intersect(var_names, net_names)
  if (length(name_clash)) {
    stop("Attribute column name(s) and network name(s) must be distinct - both are used as: ",
         paste(name_clash, collapse = ", "),
         ". Rename one or the other.", call. = FALSE)
  }

  # per-target method map (mice-style named vector, or one global method).
  # `method_named` remembers which targets the user explicitly named, for
  # the ignored-entry messages below.
  method_named <- if (is.null(names(method))) character(0) else
    names(method)[nzchar(names(method))]
  method <- .resolve_methods(method, var_names, net_names)

  mats0 <- lapply(net_list, .as_matrix_generic)
  for (k in seq_along(mats0)) diag(mats0[[k]]) <- 0
  if (!all(vapply(mats0, nrow, integer(1)) == n)) {
    stop("Every network must have as many nodes as there are rows in `data`.", call. = FALSE)
  }

  # fix structurally absent cells at 0 up front: they then count as observed
  # (never missing, never imputed) and only need special handling where the
  # dyad-level regression rows are built
  struct_list <- .validate_structural(structural, net_names, n)
  mats0 <- .apply_structural_zeros(mats0, struct_list)

  # adaptive structural constraints between networks: validate the rules,
  # error on observed contradictions, and fill missing cells already
  # determined by observed cells - the filled cells then count as observed
  # (never imputed), exactly like structural zeros
  dep_rules <- .validate_net_dependence(net_dependence, net_names)
  mats0 <- .deduce_from_dependence(mats0, dep_rules)

  attr_types <- .resolve_attr_types(data, attr_types)

  model_map <- .parse_models(models)
  bad_models <- setdiff(names(model_map), c(var_names, net_names))
  if (length(bad_models)) {
    stop("`models` references unknown variable/network(s): ",
         paste(bad_models, collapse = ", "), call. = FALSE)
  }

  var_missing <- vapply(data, function(x) any(is.na(x)), logical(1))
  net_missing <- vapply(mats0, function(mat) any(
    is.na(mat[row(mat) != col(mat)])), logical(1))
  vm_names <- var_names[var_missing]
  net_missing_names <- net_names[net_missing]

  ## ---- per-target predictor selection (netquickpred) and `targets` ----
  # Runs once per netmice() call, on the pre-chain data/network state. The
  # resulting per-target predictor lists replace the everything-plus-PCA
  # default inside the chains; missing-data variables/networks that are
  # neither targets nor selected predictors are dropped from the call.
  if (!is.null(targets)) {
    unknown_t <- setdiff(targets, c(var_names, net_names))
    if (length(unknown_t)) {
      stop("`targets` references unknown variable/network(s): ",
           paste(unknown_t, collapse = ", "), call. = FALSE)
    }
  }
  qp <- NULL
  if (inherits(predictor_selection, "netquickpred")) {
    qp <- predictor_selection
    if (!is.null(targets) && !setequal(targets, qp$targets)) {
      message("netimpute: a precomputed netquickpred object was supplied; ",
              "its own targets (", paste(qp$targets, collapse = ", "),
              ") take precedence over the `targets` argument.")
    }
  } else {
    predictor_selection <- match.arg(predictor_selection,
                                     c("all", "quickpred"))
    if (predictor_selection == "quickpred") {
      qp <- netquickpred(data, mats0, targets = targets,
                         mincor = mincor, steps = steps,
                         collin_method = collin_method,
                         collin_threshold = collin_threshold,
                         measure_set = measure_set,
                         attr_types = attr_types,
                         structural = struct_list)
    }
  }

  # variables/networks that must never be dropped, whatever the selection
  # says: anything a `models` formula depends on - directly by name or as
  # the base of a derived term like `advice_tie` or `friends_age_alter_mean`
  # (the formula must stay evaluable) - and any network bound by a
  # `net_dependence` rule
  protected <- .model_referenced_bases(model_map, var_names, net_names)
  if (!is.null(dep_rules)) {
    protected <- union(protected, unique(unlist(dep_rules)))
  }

  drop_attrs <- character(0)
  drop_nets <- character(0)
  if (!is.null(qp)) {
    drop_attrs <- setdiff(intersect(qp$drop$attributes, vm_names), protected)
    drop_nets <- setdiff(intersect(qp$drop$networks, net_missing_names),
                         protected)
  } else if (!is.null(targets)) {
    # `targets` without selection: everything not targeted (and with missing
    # data) is dropped - with the full predictor set there is no "selected
    # as predictor" escape hatch
    drop_attrs <- setdiff(vm_names, c(targets, protected))
    drop_nets <- setdiff(net_missing_names, c(targets, protected))
  }
  # feature naming inside the chains keeps using the ORIGINAL attribute set
  # for the prefix-on-clash rule, so netquickpred()'s selected names still
  # match after columns are dropped
  feat_clash_names <- var_names
  if (length(drop_attrs) || length(drop_nets)) {
    message("netimpute: dropping from the imputation (missing data, but ",
            "neither target nor selected predictor): ",
            paste(c(drop_attrs, drop_nets), collapse = ", "))
    if (length(drop_attrs)) {
      data <- data[setdiff(var_names, drop_attrs)]
      var_names <- names(data)
      attr_types <- attr_types[var_names]
    }
    if (length(drop_nets)) {
      mats0 <- mats0[setdiff(net_names, drop_nets)]
      net_names <- names(mats0)
      if (!is.null(struct_list)) struct_list <- struct_list[net_names]
    }
    var_missing <- vapply(data, function(x) any(is.na(x)), logical(1))
    net_missing <- vapply(mats0, function(mat) any(
      is.na(mat[row(mat) != col(mat)])), logical(1))
    vm_names <- var_names[var_missing]
    net_missing_names <- net_names[net_missing]
  }

  # binary networks get the Bernoulli draw under net_update = "gibbs";
  # weighted ones fall back to sequential PMM donor matching
  net_binary <- vapply(mats0, function(m) {
    v <- m[!is.na(m)]
    all(v %in% c(0, 1))
  }, logical(1))

  var_levels <- lapply(data, function(x) if (
    !is.numeric(x)) sort(unique(x[!is.na(x)])) else NULL)

  # explicitly named `method` entries that cannot take effect: multinomial
  # attributes always use polyreg, and network targets consult `method`
  # only under the simultaneous scheme without random intercepts
  multi_named <- intersect(method_named, vm_names)
  multi_named <- multi_named[vapply(multi_named, function(v)
    length(var_levels[[v]]) > 2, logical(1))]
  if (length(multi_named)) {
    message("netimpute: `method` entries for ",
            paste(multi_named, collapse = ", "),
            " will not be used - attributes with more than two categories ",
            "are always imputed with polyreg.")
  }
  net_named <- intersect(method_named, net_missing_names)
  if (length(net_named) &&
      (net_update == "gibbs" || length(net_random_intercepts))) {
    message("netimpute: `method` entries for network(s) ",
            paste(net_named, collapse = ", "), " will not be used - ties ",
            "are drawn ",
            if (net_update == "gibbs") {
              "by the tie-wise (gibbs) updater's own working model."
            } else {
              "by the mixed-model PMM step (`net_random_intercepts`)."
            })
  }

  if ("dyad" %in% net_random_intercepts) {
    undirected <- net_missing_names[vapply(
      mats0[net_missing_names],
      function(mat) isSymmetric(unname(mat)), logical(1))]
    if (length(undirected)) {
      warning("net_random_intercepts includes 'dyad' but network(s) ",
              paste(undirected, collapse = ", "),
              " are undirected: each dyad's two rows (i,j)/(j,i) are exact ",
              "duplicates, so a dyad random intercept has effectively one ",
              "unique observation per group and will absorb the residual ",
              "variance. Consider c('ego', 'alter') instead.", call. = FALSE)
    }
  }

  ry_vars <- lapply(data, function(x) !is.na(x))
  ry_nets <- lapply(mats0, function(mat) !is.na(mat) & (row(mat) != col(mat)))

  unused_models <- setdiff(names(model_map), c(vm_names, net_missing_names))
  if (length(unused_models)) {
    message("netimpute: `models` formula(s) for ",
            paste(unused_models, collapse = ", "),
            " will not be used - that variable/network has no missing values to impute.")
  }

  run_chain <- function(im) {
    .run_one_chain(
      im = im,
      data = data,
      mats0 = mats0,
      var_names = var_names,
      net_names = net_names,
      attr_types = attr_types,
      var_levels = var_levels,
      ry_vars = ry_vars,
      ry_nets = ry_nets,
      structural = struct_list,
      net_dependence = dep_rules,
      vm_names = vm_names,
      net_missing_names = net_missing_names,
      method = method,
      donors = donors,
      measure_set = measure_set,
      other_net_predictors = other_net_predictors,
      n_components = n_components,
      net_random_intercepts = net_random_intercepts,
      model_map = model_map,
      qp_pred = if (!is.null(qp)) qp$predictors else NULL,
      feat_clash_names = feat_clash_names,
      net_init = net_init,
      net_update = net_update,
      net_binary = net_binary,
      maxit = maxit,
      printFlag = printFlag,
      seed = seed
    )
  }

  if (ncores > 1) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required for ncores > 1. Install with install.packages('future').",
           call. = FALSE)
    }
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = ncores)
    if (printFlag) {
      cat("Running", m, "chains on", ncores,
          "workers via future::multisession\n")
    }
    futures <- lapply(seq_len(m), function(im) {
      future::future(run_chain(im), seed = TRUE, packages = "netimpute")
    })
    chains <- lapply(futures, future::value)
  } else {
    chains <- lapply(seq_len(m), run_chain)
  }

  imp_data <- lapply(chains, `[[`, "data")
  imp_nets <- lapply(chains, `[[`, "net_list")
  visit_orders <- lapply(chains, `[[`, "visit_order")

  net_diag_names <- c("density", "reciprocity",
                      "transitivity", "n_isolates", "avg_inv_geodesic")
  chainMean <- array(NA_real_, dim = c(length(vm_names), maxit, m),
                      dimnames = list(vm_names, seq_len(maxit), seq_len(m)))
  chainVar <- chainMean
  netChain <- array(NA_real_, dim = c(length(net_names),
                                      length(net_diag_names), maxit, m),
                     dimnames = list(net_names,
                                     net_diag_names,
                                     seq_len(maxit),
                                     seq_len(m)))
  for (im in seq_len(m)) {
    if (length(vm_names)) {
      chainMean[, , im] <- chains[[im]]$chainMean
      chainVar[, , im]  <- chains[[im]]$chainVar
    }
    netChain[, , , im] <- chains[[im]]$netChain
  }

  structure(
    list(
      data = data,
      net_list = mats0,
      m = m,
      maxit = maxit,
      method = method,
      donors = donors,
      net_random_intercepts = net_random_intercepts,
      structural = struct_list,
      net_dependence = dep_rules,
      models = model_map,
      predictor_selection = if (is.null(qp)) "all" else qp,
      targets = targets,
      net_init = net_init,
      net_update = net_update,
      ncores = ncores,
      imp = imp_data,
      imp_nets = imp_nets,
      visit_orders = visit_orders,
      chainMean = chainMean,
      chainVar = chainVar,
      netChain = netChain,
      var_missing = vm_names,
      net_missing = net_missing_names
    ),
    class = "netmids"
  )
}

#' Extract one completed data.frame/network-list from a netmids object
#' @param x A `netmids` object from \code{\link{netmice}}.
#' @param action Which imputation to extract (integer, 1..m).
#' @return A list with `data` (completed data.frame) and `net_list`
#'   (completed list of adjacency matrices).
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 25
#' friends <- matrix(rbinom(n * n, 1, 0.15), n, n)
#' diag(friends) <- 0
#' friends[sample(which(row(friends) != col(friends)), 30)] <- NA
#' attrs <- data.frame(age = rnorm(n, 35, 8),
#'                     gender = sample(c("F", "M"), n, TRUE))
#' attrs$age[sample(n, 4)] <- NA
#'
#' \donttest{
#' fit <- netmice(attrs, list(friends = friends), m = 2, maxit = 2,
#'                seed = 1, printFlag = FALSE)
#' completed <- complete_netmice(fit, 1)
#' str(completed$data)
#' anyNA(completed$net_list$friends)
#' }
complete_netmice <- function(x, action = 1) {
  if (!inherits(x, "netmids")) stop("`x` must be a netmids object from netmice().", call. = FALSE)
  if (action < 1 || action > x$m) stop("`action` must be between 1 and ", x$m, ".", call. = FALSE)
  list(data = x$imp[[action]], net_list = x$imp_nets[[action]])
}

#' Print a summary of a netmids object
#'
#' Reports the imputation settings (number of imputations and iterations,
#' univariate method(s), donors, cores), which attributes and networks had
#' missing values, and any non-default features in use (tie-wise updating,
#' random intercepts, structural zeros, custom `models`, per-target
#' predictor selection).
#'
#' @param x A `netmids` object from \code{\link{netmice}}.
#' @param ... Ignored, for consistency with the generic.
#' @return `x`, invisibly.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 25
#' friends <- matrix(rbinom(n * n, 1, 0.15), n, n)
#' diag(friends) <- 0
#' friends[sample(which(row(friends) != col(friends)), 30)] <- NA
#' attrs <- data.frame(age = rnorm(n, 35, 8),
#'                     gender = sample(c("F", "M"), n, TRUE))
#' attrs$age[sample(n, 4)] <- NA
#'
#' \donttest{
#' fit <- netmice(attrs, list(friends = friends), m = 2, maxit = 2,
#'                seed = 1, printFlag = FALSE)
#' print(fit)
#' }
print.netmids <- function(x, ...) {
  cat("Class: netmids\n")
  meths <- unique(unname(x$method))
  cat("Imputations (m):", x$m, "  Iterations (maxit):", x$maxit,
      "  Method:", if (length(meths) == 1) meths else "per-target",
      "  Donors:", x$donors, "  Cores:", x$ncores, "\n")
  if (length(meths) > 1) {
    mm <- x$method[intersect(names(x$method),
                             c(x$var_missing, x$net_missing))]
    cat("Per-target methods:",
        paste(names(mm), mm, sep = " = ", collapse = ", "), "\n")
  }
  cat("Variables with missingness:", paste(x$var_missing, collapse = ", "), "\n")
  cat("Networks with missingness:", paste(x$net_missing, collapse = ", "), "\n")
  if (identical(x$net_update, "gibbs")) {
    cat("Network-tie updating: sequential Gibbs (single-tie draws with",
        "change statistics)\n")
  }
  if (length(x$net_random_intercepts)) {
    cat("Network-tie random intercepts (lme4):",
        paste(x$net_random_intercepts, collapse = ", "), "\n")
  }
  if (length(x$structural)) {
    has_struct <- !vapply(x$structural, is.null, logical(1))
    cat("Structurally absent ties (fixed zeros) in:",
        paste(names(x$structural)[has_struct], collapse = ", "), "\n")
  }
  if (length(x$models)) {
    cat("Custom models for:", paste(names(x$models), collapse = ", "), "\n")
  }
  if (inherits(x$predictor_selection, "netquickpred")) {
    qp <- x$predictor_selection
    cat("Per-target predictor selection (netquickpred): mincor", qp$mincor,
        "- steps", qp$steps, "- collinearity", qp$collin_method, "\n")
    if (length(qp$drop$attributes) || length(qp$drop$networks)) {
      cat("Dropped (not target/predictor):",
          paste(c(qp$drop$attributes, qp$drop$networks), collapse = ", "),
          "\n")
    }
  }
  invisible(x)
}
