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
# univariate imputation step is dispatched through .impute_univariate() so
# methods beyond "pmm" can be added later without changing the call surface.

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
  state <- list(it = 0L, im = 0L, dep = "y", meth = method, log = FALSE)
  loggedEvents <- NULL
  switch(method,
    pmm = mice::mice.impute.pmm(y = y, ry = ry, x = x, donors = donors),
    stop("Unsupported imputation method: '", method, "'. Currently only 'pmm' is implemented.",
         call. = FALSE)
  )
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

#' Sequential ("Gibbs") single-tie imputation for one network visit
#'
#' The ERGM-flavoured alternative to the simultaneous PMM tie update. A
#' working model is fit once per visit on the observed dyads - logistic
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
#' `reciprocity` and `twopath` are updated per draw.
#'
#' `check = TRUE` recomputes the binarized adjacency and the two-path count
#' matrix from scratch at the end and errors if the incremental bookkeeping
#' diverged - used by the unit tests, skipped in production calls.
#' @keywords internal
.impute_ties_gibbs <- function(d, ry, x, mat, binary, donors, check = FALSE) {
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

  y_obs <- d$y[ry]
  yhat_obs <- if (!binary) eta_all[ry] else NULL

  B <- (mat != 0) * 1
  Tp <- B %*% B  # two-path counts; maintained incrementally below

  for (r in sample(which(!ry))) {
    i <- d$i[r]
    j <- d$j[r]
    eta <- eta_base[r] + b_recip * mat[j, i] +
      b_twop * as.numeric(Tp[i, j] > 0)
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
        # Prefix measure names with the network name only when there is more
        # than one network (or a bare name would collide with a raw attribute
        # of the same name) - with a single, non-colliding network, bare
        # names like "indegree" are unambiguous and nicer to use in `models`.
        net_feats_list <- lapply(seq_along(cur_gs), function(k) {
          out <- net_measures(cur_gs[[k]], cur_data, measure_set = measure_set,
                              attr_types = attr_types)
          out <- out[setdiff(names(out), c("node_id", ego_feats))]
          if (length(net_names) > 1) {
            names(out) <- paste0(net_names[k], "_", names(out))
          } else {
            clash <- names(out) %in% names(cur_data)
            names(out)[clash] <- paste0(net_names[k], "_", names(out)[clash])
          }
          out
        })
        net_feats_df <- do.call(cbind, net_feats_list)
        auto_x <- .clean_predictor_matrix(
          cbind(.prep_pca_matrix(other_data), .prep_pca_matrix(net_feats_df)),
          max_cols = .safe_max_cols(sum(ry)),
          ry = ry
        )

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
          imp_vals <- .impute_univariate(method,
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
          imp_codes <- .impute_univariate(method,
                                          y = y_codes,
                                          ry = ry,
                                          x = x,
                                          donors = donors)
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
        # not fit on (no zero inflation) and never receive imputed values
        built <- .build_dyad_data(cur_mats,
                                  cur_data,
                                  target_idx = k,
                                  attr_types = attr_types,
                                  other_net_predictors = other_net_predictors,
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
        auto_x <- .clean_predictor_matrix(d[setdiff(names(d), drop_cols)],
                                           max_cols = .safe_max_cols(sum(ry)),
                                           ry = ry,
                                           keep_raw = keep_raw)

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
            imp_vals <- .impute_univariate(method,
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
#' updated (via \code{\link{dyad_regression}}'s predictor set). Only
#' predictive mean matching is implemented for now (`method = "pmm"`);
#' `method` is kept as an argument so alternatives can be added later without
#' changing the call signature.
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
#' @param method Univariate imputation method for numeric and binary
#'   attributes (and for network ties); only `"pmm"` is currently
#'   implemented. Nominal attributes with more than two categories are
#'   *always* imputed with a proper multinomial logistic regression
#'   (`mice::mice.impute.polyreg()`, via `nnet::multinom()`) regardless of
#'   `method` - a linear-model-based method like `"pmm"` would otherwise
#'   have to impose an arbitrary numeric ordering on the categories.
#' @param donors Number of PMM donors (default 5).
#' @param measure_set Which network measures supply predictors for attribute
#'   imputation, passed to \code{\link{net_measures}}: any combination of
#'   "core" (default, faster), "full", "homophily", and individual measure
#'   names - named sets and single measures are unioned, e.g.
#'   \code{c("core", "total_degree")}. Note the "homophily" block (part of
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
#'   `list("happiness ~ indegree + age + performance * gender")`.
#' @param net_update How a network's missing cells are updated at each
#'   visit. `"simultaneous"` (default; the behavior to date): one PMM
#'   working model per visit imputes every missing cell at once, all from
#'   the same snapshot of the endogenous statistics. `"gibbs"`: sequential
#'   single-tie imputation - a working model is fit once per visit on the
#'   observed dyads (logistic regression for a binary network, linear
#'   regression otherwise), one coefficient vector is drawn from its
#'   asymptotic posterior, and the missing cells are then visited one at a
#'   time in random order, each draw using the *current* values of the
#'   endogenous statistics (`reciprocity`, `twopath`), which are refreshed
#'   after every single draw via change statistics (O(1) for reciprocity,
#'   two O(n) vector updates on the running two-path count matrix - never a
#'   full recomputation). Binary ties are drawn from
#'   `Bernoulli(plogis(eta))`; weighted ties by PMM donor matching on the
#'   fitted values. Each tie is thus drawn conditional on all previously
#'   imputed ties - the full-conditional (Gibbs) update of an ERGM-style
#'   model with pseudo-likelihood-estimated parameters, restricted to the
#'   missing cells - which yields properly dependent within-network
#'   imputations and avoids the synchronous-update artifact of the
#'   simultaneous scheme. Static predictors (attribute terms, other-network
#'   terms, `models` extras) are evaluated once per visit and held fixed
#'   during the sweep. Not compatible with `net_random_intercepts`;
#'   attribute imputation is unaffected.
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
#'   `m`, `maxit`, `method`, `donors`, `net_init`, `net_update`,
#'   `structural` (the
#'   per-network list of structural-zero matrices, or `NULL`),
#'   `net_dependence` (the normalized rule list, or `NULL`),
#'   `models`, `imp` (list of `m` completed attribute data.frames),
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
#' \dontrun{
#' set.seed(1)
#' n <- 40
#' friends <- matrix(rbinom(n * n, 1, 0.1), n, n); diag(friends) <- 0
#' friends[sample(length(friends), 50)] <- NA
#' attrs <- data.frame(age = rnorm(n, 35, 8),
#'                      performance = rnorm(n),
#'                      gender = sample(c("F", "M"), n, TRUE),
#'                      happiness = rnorm(n))
#' attrs$happiness[sample(n, 5)] <- NA
#' fit <- netmice(attrs, list(friends = friends), m = 3, maxit = 3,
#'                models = list("happiness ~ indegree + age + performance * gender"))
#' out <- complete_netmice(fit, 1)
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
                    net_init = c("zero", "sample"),
                    net_update = c("simultaneous", "gibbs"),
                    ncores = 1L,
                    seed = NA,
                    printFlag = TRUE) {

  method <- match.arg(method, choices = "pmm")
  net_init <- match.arg(net_init)
  net_update <- match.arg(net_update)
  measure_set <- .resolve_measure_set(measure_set)
  other_net_predictors <- match.arg(other_net_predictors)
  if (!is.null(net_random_intercepts)) {
    net_random_intercepts <- match.arg(net_random_intercepts,
                                       c("ego", "alter", "dyad"),
                                       several.ok = TRUE)
  }
  if (net_update == "gibbs" && length(net_random_intercepts)) {
    stop("net_update = 'gibbs' is not compatible with `net_random_intercepts`: ",
         "the sequential sampler uses its own logistic/linear working model. ",
         "Use one or the other.", call. = FALSE)
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

  # binary networks get the Bernoulli draw under net_update = "gibbs";
  # weighted ones fall back to sequential PMM donor matching
  net_binary <- vapply(mats0, function(m) {
    v <- m[!is.na(m)]
    all(v %in% c(0, 1))
  }, logical(1))

  attr_types <- .resolve_attr_types(data, attr_types)
  var_levels <- lapply(data, function(x) if (
    !is.numeric(x)) sort(unique(x[!is.na(x)])) else NULL)

  var_missing <- vapply(data, function(x) any(is.na(x)), logical(1))
  net_missing <- vapply(mats0, function(mat) any(
    is.na(mat[row(mat) != col(mat)])), logical(1))
  vm_names <- var_names[var_missing]
  net_missing_names <- net_names[net_missing]

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

  model_map <- .parse_models(models)
  bad_models <- setdiff(names(model_map), c(var_names, net_names))
  if (length(bad_models)) {
    stop("`models` references unknown variable/network(s): ",
         paste(bad_models, collapse = ", "), call. = FALSE)
  }
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
complete_netmice <- function(x, action = 1) {
  if (!inherits(x, "netmids")) stop("`x` must be a netmids object from netmice().", call. = FALSE)
  if (action < 1 || action > x$m) stop("`action` must be between 1 and ", x$m, ".", call. = FALSE)
  list(data = x$imp[[action]], net_list = x$imp_nets[[action]])
}

#' @export
print.netmids <- function(x, ...) {
  cat("Class: netmids\n")
  cat("Imputations (m):", x$m, "  Iterations (maxit):", x$maxit,
      "  Method:", x$method, "  Donors:", x$donors, "  Cores:", x$ncores, "\n")
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
  invisible(x)
}
