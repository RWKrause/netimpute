# A mice-style chained-equations loop that jointly imputes a nodal-attribute
# data.frame AND a list of networks. Within one iteration, every attribute
# and every network with missingness is visited *once*, in a single
# interleaved, randomly-ordered sequence (a different random order per
# imputation chain m) - not "all attributes then all networks" - so that
# whichever was updated most recently (whether that was an attribute or a
# network) always feeds into the next target's predictors. Network-derived
# predictors for an attribute (net_measures_core()/full()) and
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
.init_fill_matrix <- function(mat) {
  diag(mat) <- 0
  off <- row(mat) != col(mat)
  miss <- off & is.na(mat)
  if (any(miss)) {
    obs_vals <- mat[off & !is.na(mat)]
    if (!length(obs_vals)) stop("A network is entirely missing; cannot initialize imputation.", call. = FALSE)
    mat[miss] <- sample(obs_vals, sum(miss), replace = TRUE)
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

  # `method` here is either the user-facing method for numeric/binary
  # targets ("pmm", for now) or the literal string "polyreg", used
  # internally (never user-selectable via netmice()'s `method` argument) for
  # nominal attributes with more than two categories - see .run_one_chain().
  if (method == "polyreg") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Package 'nnet' is required for multinomial imputation of attributes with more than ",
           "two categories. Install with install.packages('nnet').", call. = FALSE)
    }
    # mice::mice.impute.polyreg() fits a proper multinomial logistic
    # regression (nnet::multinom) and draws from the predicted class
    # probabilities - unlike pmm's linear-fit-then-nearest-donor approach,
    # it doesn't impose a spurious numeric ordering on unordered categories.
    # It works standalone (no updateLog()/state dependency like pmm below),
    # and takes/returns the category labels directly, so no integer coding
    # round-trip is needed.
    return(mice::mice.impute.polyreg(y = y, ry = ry, x = x))
  }

  # mice::mice.impute.pmm() is normally only ever called from inside a full
  # mice() run. Whenever its internal estimice() hits a (near-)singular
  # design matrix - routine here, given the auto-generated predictor sets
  # are often collinear - it calls mice's updateLog(), which walks up the
  # call stack (via parent.frame()) looking for `state`/`loggedEvents`
  # objects that mice()'s own sampler() normally creates. Calling pmm
  # standalone skips that setup entirely, so updateLog() fails with
  # "object 'state' not found". Defining minimal stand-ins right here (this
  # frame is within updateLog's search range) satisfies that lookup without
  # needing to go through mice()'s full machinery.
  state <- list(it = 0L, im = 0L, dep = "y", meth = method, log = FALSE)
  loggedEvents <- NULL
  switch(method,
    pmm = mice::mice.impute.pmm(y = y, ry = ry, x = x, donors = donors),
    stop("Unsupported imputation method: '", method, "'. Currently only 'pmm' is implemented.",
         call. = FALSE)
  )
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
#' instead of degrading gracefully. Pre-emptively reducing dimensionality
#' here keeps p comfortably below n by construction, the same fix already
#' offered explicitly via net_predictors(output = "pca") - just applied
#' automatically inside the imputation loop where the user has no direct
#' control over the predictor count.
#' @keywords internal
.clean_predictor_matrix <- function(x, max_cols = NULL) {
  x <- as.matrix(x)
  x <- apply(x, 2, function(col) {
    if (all(is.na(col))) return(rep(0, length(col)))
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  })
  var0 <- apply(x, 2, stats::var)
  x <- x[, is.finite(var0) & var0 > 0, drop = FALSE]

  if (!is.null(max_cols) && max_cols >= 1 && ncol(x) > max_cols) {
    pca <- stats::prcomp(x, center = TRUE, scale. = TRUE)
    x <- pca$x[, seq_len(max_cols), drop = FALSE]
  }
  x
}

#' @keywords internal
.safe_max_cols <- function(n_obs) {
  max(5, floor(n_obs / 3))
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
.model_extra_terms <- function(formula, data) {
  rhs <- stats::delete.response(stats::terms(formula))
  mf <- tryCatch(
    stats::model.frame(rhs, data = data, na.action = stats::na.pass),
    error = function(e) {
      stop("netimpute: could not evaluate a `models` formula (", conditionMessage(e), "). ",
           "Every predictor named in a formula must exist in `data`, or - for network-derived ",
           "terms in an attribute's formula, or for any term in a network's formula - match the ",
           "auto-generated predictor names (e.g. 'friends_indegree' for an attribute's formula ",
           "when multiple networks are supplied, or 'age_absdiff'/'friends_tie' for a network's ",
           "own dyad-level formula). See ?netmice.", call. = FALSE)
    }
  )
  mm <- stats::model.matrix(rhs, data = mf)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  .clean_predictor_matrix(mm)
}

#' Run one full imputation chain (one of the `m` imputations)
#' @keywords internal
.run_one_chain <- function(im, data, mats0, var_names, net_names, attr_types, var_levels,
                            ry_vars, ry_nets, vm_names, net_missing_names,
                            method, donors, measure_set, other_net_predictors, n_components,
                            model_map, maxit, printFlag, seed) {
  # Pin the RNG algorithm explicitly, not just the seed value: future's
  # seed = TRUE (used for ncores > 1) switches multisession workers to
  # "L'Ecuyer-CMRG" for proper parallel streams, so the *same* seed integer
  # would otherwise produce a *different* random sequence there than in the
  # calling session's default "Mersenne-Twister" - silently breaking the
  # sequential-vs-parallel reproducibility this is meant to guarantee.
  if (!is.na(seed)) set.seed(seed + im, kind = "Mersenne-Twister")

  cur_data <- data
  for (v in vm_names) cur_data[[v]] <- .init_fill_vector(cur_data[[v]])
  cur_mats <- mats0
  for (nm in net_missing_names) cur_mats[[nm]] <- .init_fill_matrix(cur_mats[[nm]])

  # A single interleaved, randomly-ordered visit sequence over BOTH
  # attributes and networks, fixed for this chain (different per chain m).
  targets <- sample(c(vm_names, net_missing_names))

  net_diag_names <- c("density", "reciprocity", "transitivity", "n_isolates", "avg_inv_geodesic")
  chainMean <- matrix(NA_real_, length(vm_names), maxit, dimnames = list(vm_names, seq_len(maxit)))
  chainVar  <- chainMean
  netChain  <- array(NA_real_, dim = c(length(net_names), length(net_diag_names), maxit),
                      dimnames = list(net_names, net_diag_names, seq_len(maxit)))

  for (it in seq_len(maxit)) {
    if (printFlag) cat(" imputation", im, "- iter", it, "- order:", paste(targets, collapse = " -> "), "\n")

    for (tgt in targets) {
      if (tgt %in% vm_names) {
        v <- tgt
        ry <- ry_vars[[v]]
        # Network-derived features are rebuilt from cur_mats *as they stand
        # right now* - i.e. reflecting any network updated earlier in this
        # same sweep - not a snapshot from the start of the iteration.
        cur_gs <- stats::setNames(lapply(cur_mats, .mat_to_igraph), net_names)

        # `v` is excluded from the attribute set used to build homophily-type
        # network features: those (e.g. v's own E-I index / alter mean) are a
        # function of v's current, not-yet-finalized value at each node, so
        # using them as predictors for v would be circular. Other attributes'
        # homophily features are fine and are kept.
        other_data <- cur_data[setdiff(var_names, v)]
        # Prefix measure names with the network name only when there is more
        # than one network (or a bare name would collide with a raw attribute
        # of the same name) - with a single, non-colliding network, bare
        # names like "indegree" are unambiguous and nicer to use in `models`.
        net_feats_list <- lapply(seq_along(cur_gs), function(k) {
          fn <- if (measure_set == "full") net_measures_full else net_measures_core
          out <- fn(cur_gs[[k]], other_data, attr_types = attr_types)
          out <- out[setdiff(names(out), "node_id")]
          if (length(net_names) > 1) {
            names(out) <- paste0(net_names[k], "_", names(out))
          } else {
            clash <- names(out) %in% names(other_data)
            names(out)[clash] <- paste0(net_names[k], "_", names(out)[clash])
          }
          out
        })
        net_feats_df <- do.call(cbind, net_feats_list)
        auto_x <- .clean_predictor_matrix(
          cbind(.prep_pca_matrix(other_data), .prep_pca_matrix(net_feats_df)),
          max_cols = .safe_max_cols(sum(ry))
        )

        if (!is.null(model_map[[v]])) {
          combined_df <- cbind(other_data, net_feats_df)
          extra <- .model_extra_terms(model_map[[v]], combined_df)
          x <- cbind(auto_x, extra)
          x <- x[, !duplicated(colnames(x)), drop = FALSE]
        } else {
          x <- auto_x
        }

        lvls <- var_levels[[v]]
        if (is.numeric(cur_data[[v]])) {
          imp_vals <- .impute_univariate(method, y = cur_data[[v]], ry = ry, x = x, donors = donors)
        } else if (length(lvls) > 2) {
          # Nominal, >2 categories: a proper multinomial model, not a linear-
          # regression-on-integer-codes hack (see .impute_univariate()) -
          # avoids imposing a spurious numeric ordering on the categories.
          imp_vals <- .impute_univariate("polyreg", y = cur_data[[v]], ry = ry, x = x, donors = donors)
        } else {
          # Binary: mice::mice.impute.pmm() requires a numeric y (it hits
          # lm.fit() internally), so we code the two categories to 0/1, run
          # PMM in that space, and map back. Harmless here - a linear
          # probability model on a 2-level outcome, and PMM only ever
          # returns an *observed* donor's value either way.
          y_codes <- as.integer(factor(cur_data[[v]], levels = lvls))
          imp_codes <- .impute_univariate(method, y = y_codes, ry = ry, x = x, donors = donors)
          imp_vals <- lvls[imp_codes]
        }
        cur_data[[v]][!ry] <- imp_vals

      } else {
        k <- match(tgt, net_names)
        built <- .build_dyad_data(cur_mats, cur_data, target_idx = k, attr_types = attr_types,
                                   other_net_predictors = other_net_predictors,
                                   n_components = n_components)
        d <- built$data
        ry <- ry_nets[[tgt]][cbind(d$i, d$j)]
        auto_x <- .clean_predictor_matrix(d[setdiff(names(d), c("i", "j", "y"))],
                                           max_cols = .safe_max_cols(sum(ry)))

        if (!is.null(model_map[[tgt]])) {
          extra <- .model_extra_terms(model_map[[tgt]], d)
          x <- cbind(auto_x, extra)
          x <- x[, !duplicated(colnames(x)), drop = FALSE]
        } else {
          x <- auto_x
        }

        y <- d$y
        imp_vals <- .impute_univariate(method, y = y, ry = ry, x = x, donors = donors)
        mis_idx <- which(!ry)
        cur_mats[[tgt]][cbind(d$i[mis_idx], d$j[mis_idx])] <- imp_vals
      }
    }

    ## ---- diagnostics (once per full sweep) ----
    for (v in vm_names) {
      vals <- cur_data[[v]][!ry_vars[[v]]]
      if (!is.numeric(vals)) vals <- as.numeric(factor(vals, levels = var_levels[[v]]))
      chainMean[v, it] <- mean(vals, na.rm = TRUE)
      chainVar[v, it]  <- stats::var(vals, na.rm = TRUE)
    }
    for (nm in net_names) {
      dk <- net_diagnostics(cur_mats[[nm]])
      netChain[nm, , it] <- unlist(dk[net_diag_names])
    }
  }

  list(data = cur_data, net_list = cur_mats, chainMean = chainMean, chainVar = chainVar,
       netChain = netChain, visit_order = targets)
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
#' most recently updated (via \code{\link{net_measures_core}}/
#' \code{\link{net_measures_full}}); each network's missing ties are imputed
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
#' @param measure_set "core" (default, faster) or "full" - which
#'   network-measure function supplies predictors for attribute imputation.
#' @param attr_types Optional named attribute-type overrides.
#' @param other_net_predictors,n_components Passed to
#'   \code{\link{dyad_regression}}'s dyad-data builder for the "other
#'   networks" terms - "raw" (all terms) or "pca" (first `n_components`
#'   components).
#' @param models Optional list/character vector of formula strings (or
#'   formula objects), one per attribute/network you want a custom model
#'   for - see Details. E.g.
#'   `list("happiness ~ indegree + age + performance * gender")`.
#' @param ncores Number of parallel workers for the `m` chains via
#'   \pkg{future} (default 1 = sequential). See Details for the
#'   installed-package requirement when `ncores > 1`.
#' @param seed Optional RNG seed (chain `im` uses `seed + im`, so results are
#'   identical whether run sequentially or in parallel).
#' @param printFlag Logical; print progress (imputation/iteration/visit
#'   order) as in `mice`.
#'
#' @return An object of class `"netmids"` with elements: `data`, `net_list`
#'   (original, NA-preserving inputs), `m`, `maxit`, `method`, `donors`,
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
netmice <- function(data, net_list,
                     m = 5, maxit = 20, method = "pmm", donors = 5,
                     measure_set = c("core", "full"),
                     attr_types = NULL,
                     other_net_predictors = c("raw", "pca"), n_components = 3,
                     models = NULL,
                     ncores = 1L,
                     seed = NA, printFlag = TRUE) {

  method <- match.arg(method, choices = "pmm")
  measure_set <- match.arg(measure_set)
  other_net_predictors <- match.arg(other_net_predictors)

  data <- as.data.frame(data)
  n <- nrow(data)
  var_names <- names(data)

  net_names <- names(net_list)
  if (is.null(net_names)) net_names <- paste0("net", seq_along(net_list))
  names(net_list) <- net_names

  name_clash <- intersect(var_names, net_names)
  if (length(name_clash)) {
    stop("Attribute column name(s) and network name(s) must be distinct - both are used as: ",
         paste(name_clash, collapse = ", "), ". Rename one or the other.", call. = FALSE)
  }

  mats0 <- lapply(net_list, .as_matrix_generic)
  for (k in seq_along(mats0)) diag(mats0[[k]]) <- 0
  if (!all(vapply(mats0, nrow, integer(1)) == n)) {
    stop("Every network must have as many nodes as there are rows in `data`.", call. = FALSE)
  }

  attr_types <- .resolve_attr_types(data, attr_types)
  var_levels <- lapply(data, function(x) if (!is.numeric(x)) sort(unique(x[!is.na(x)])) else NULL)

  var_missing <- vapply(data, function(x) any(is.na(x)), logical(1))
  net_missing <- vapply(mats0, function(mat) any(is.na(mat[row(mat) != col(mat)])), logical(1))
  vm_names <- var_names[var_missing]
  net_missing_names <- net_names[net_missing]

  ry_vars <- lapply(data, function(x) !is.na(x))
  ry_nets <- lapply(mats0, function(mat) !is.na(mat) & (row(mat) != col(mat)))

  model_map <- .parse_models(models)
  bad_models <- setdiff(names(model_map), c(var_names, net_names))
  if (length(bad_models)) {
    stop("`models` references unknown variable/network(s) on the left-hand side: ",
         paste(bad_models, collapse = ", "), call. = FALSE)
  }
  unused_models <- setdiff(names(model_map), c(vm_names, net_missing_names))
  if (length(unused_models)) {
    message("netimpute: `models` formula(s) for ", paste(unused_models, collapse = ", "),
            " will not be used - that variable/network has no missing values to impute.")
  }

  run_chain <- function(im) {
    .run_one_chain(
      im = im, data = data, mats0 = mats0, var_names = var_names, net_names = net_names,
      attr_types = attr_types, var_levels = var_levels, ry_vars = ry_vars, ry_nets = ry_nets,
      vm_names = vm_names, net_missing_names = net_missing_names,
      method = method, donors = donors, measure_set = measure_set,
      other_net_predictors = other_net_predictors, n_components = n_components,
      model_map = model_map, maxit = maxit, printFlag = printFlag, seed = seed
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
    if (printFlag) cat("Running", m, "chains on", ncores, "workers via future::multisession\n")
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

  net_diag_names <- c("density", "reciprocity", "transitivity", "n_isolates", "avg_inv_geodesic")
  chainMean <- array(NA_real_, dim = c(length(vm_names), maxit, m),
                      dimnames = list(vm_names, seq_len(maxit), seq_len(m)))
  chainVar <- chainMean
  netChain <- array(NA_real_, dim = c(length(net_names), length(net_diag_names), maxit, m),
                     dimnames = list(net_names, net_diag_names, seq_len(maxit), seq_len(m)))
  for (im in seq_len(m)) {
    if (length(vm_names)) {
      chainMean[, , im] <- chains[[im]]$chainMean
      chainVar[, , im]  <- chains[[im]]$chainVar
    }
    netChain[, , , im] <- chains[[im]]$netChain
  }

  structure(
    list(
      data = data, net_list = mats0, m = m, maxit = maxit, method = method, donors = donors,
      models = model_map, ncores = ncores,
      imp = imp_data, imp_nets = imp_nets, visit_orders = visit_orders,
      chainMean = chainMean, chainVar = chainVar, netChain = netChain,
      var_missing = vm_names, net_missing = net_missing_names
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
  if (length(x$models)) {
    cat("Custom models for:", paste(names(x$models), collapse = ", "), "\n")
  }
  invisible(x)
}
