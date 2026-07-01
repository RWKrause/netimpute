# Cell-level (dyadic) regression: vectorize adjacency matrices and regress a
# target network's ties on ego/alter/(diff|same) nodal-attribute terms,
# reciprocity, two-path count, and terms derived from all other networks in
# the list. This is the point estimate ("MR-QAP without the permutation
# test") used both standalone and as the workhorse inside netmice().
#
# IMPORTANT: because a missing tie has no natural representation in an
# igraph/network edge list, supply networks with unknown ties as plain
# (possibly asymmetric) adjacency matrices with NA marking unknown cells.
# igraph/network objects are still accepted but are treated as fully observed.

#' @keywords internal
.as_matrix_generic <- function(net) {
  if (is.matrix(net)) {
    .reject_negative(net)  # .reject_negative() already ignores NA (unknown ties)
    return(net)
  }
  g <- .as_igraph(net)  # already validates non-negativity
  as.matrix(igraph::as_adjacency_matrix(
    g, sparse = FALSE,
    attr = if ("weight" %in% igraph::edge_attr_names(g)) "weight" else NULL
  ))
}

#' Build the dyad-level (cell-level) design matrix for one target network
#'
#' @param mats named list of adjacency matrices (all n x n, same node order)
#' @param attributes data.frame of node attributes aligned to the matrices
#' @param target_idx index of the target network within `mats`
#' @param attr_types optional named attribute-type overrides
#' @param other_net_predictors "raw" (all terms) or "pca" (first `n_components`
#'   principal components of all "other network" terms)
#' @param n_components number of PCA components to retain when
#'   `other_net_predictors = "pca"`
#' @keywords internal
.build_dyad_data <- function(mats, attributes, target_idx, attr_types = NULL,
                              other_net_predictors = c("raw", "pca"),
                              n_components = 3) {
  other_net_predictors <- match.arg(other_net_predictors)
  net_names <- names(mats)
  n <- nrow(mats[[1]])
  if (!all(vapply(mats, nrow, integer(1)) == n) || !all(vapply(mats, ncol, integer(1)) == n)) {
    stop("All networks in the list must be square matrices of the same size.", call. = FALSE)
  }

  target_mat <- mats[[target_idx]]
  other_mats <- mats[-target_idx]

  types <- .resolve_attr_types(attributes, attr_types)

  base <- expand.grid(i = seq_len(n), j = seq_len(n))
  keep <- base$i != base$j
  base <- base[keep, ]

  y <- as.vector(target_mat)[keep]
  recip <- as.vector(t(target_mat))[keep]
  twopath <- as.vector(target_mat %*% target_mat)[keep]

  attr_predictors <- list()
  for (nm in names(attributes)) {
    x <- attributes[[nm]]
    type <- types[[nm]]
    if (type == "continuous") {
      x <- as.numeric(x)
      ego   <- matrix(x, n, n)[keep]
      alter <- matrix(x, n, n, byrow = TRUE)[keep]
      # NOTE: we use |ego - alter|, not the signed difference. The signed
      # difference (ego - alter) is an exact linear combination of the ego
      # and alter columns already in the model, which would make the design
      # matrix singular. The absolute difference is the standard homophily
      # term in dyadic/QAP models precisely because it is *not* redundant
      # with the two level terms.
      attr_predictors[[paste0(nm, "_ego")]]     <- ego
      attr_predictors[[paste0(nm, "_alter")]]   <- alter
      attr_predictors[[paste0(nm, "_absdiff")]] <- abs(ego - alter)
    } else {
      xf <- factor(x)
      xc <- as.character(xf)
      same <- (matrix(xc, n, n)[keep] == matrix(xc, n, n, byrow = TRUE)[keep]) * 1
      # dummy-code ego/alter category membership (nodefactor-style), baseline dropped
      for (lv in levels(xf)[-1]) {
        dv <- as.numeric(xc == lv)
        attr_predictors[[paste0(nm, "_ego_", lv)]]   <- matrix(dv, n, n)[keep]
        attr_predictors[[paste0(nm, "_alter_", lv)]] <- matrix(dv, n, n, byrow = TRUE)[keep]
      }
      attr_predictors[[paste0(nm, "_same")]] <- same
    }
  }

  other_predictors <- list()
  for (k in seq_along(other_mats)) {
    onm <- names(other_mats)[k]
    if (is.null(onm) || onm == "") onm <- paste0("other", k)
    mk <- other_mats[[k]]
    outdeg_k <- rowSums(mk, na.rm = TRUE)
    indeg_k  <- colSums(mk, na.rm = TRUE)
    other_predictors[[paste0(onm, "_tie")]]          <- as.vector(mk)[keep]
    other_predictors[[paste0(onm, "_recip")]]        <- as.vector(t(mk))[keep]
    other_predictors[[paste0(onm, "_ego_outdeg")]]   <- matrix(outdeg_k, n, n)[keep]
    other_predictors[[paste0(onm, "_alter_indeg")]]  <- matrix(indeg_k, n, n, byrow = TRUE)[keep]
  }

  pca_model <- NULL
  if (length(other_predictors)) {
    if (other_net_predictors == "pca") {
      om <- as.matrix(as.data.frame(other_predictors))
      om[is.na(om)] <- 0
      var0 <- apply(om, 2, stats::var)
      om <- om[, is.finite(var0) & var0 > 0, drop = FALSE]
      if (ncol(om) > 0) {
        k_comp <- min(n_components, ncol(om))
        pca_model <- stats::prcomp(om, center = TRUE, scale. = TRUE)
        comps <- as.data.frame(pca_model$x[, seq_len(k_comp), drop = FALSE])
        names(comps) <- paste0("other_net_PC", seq_len(k_comp))
        other_predictors <- as.list(comps)
      } else {
        other_predictors <- list()
      }
    }
  }

  # NOTE: cbind(data.frame, NULL) is NOT a no-op - it errors ("arguments
  # imply differing number of rows: n, 0"), so the single-network case
  # (other_predictors empty -> NULL) must never be passed to cbind directly.
  # Building a list and dropping NULLs first (do.call) sidesteps that.
  parts <- list(
    base,
    y = y,
    reciprocity = recip,
    log_twopath = log1p(twopath),
    as.data.frame(attr_predictors, check.names = FALSE),
    if (length(other_predictors)) as.data.frame(other_predictors, check.names = FALSE) else NULL
  )
  data <- do.call(cbind, Filter(Negate(is.null), parts))
  rownames(data) <- NULL

  list(data = data, pca_model = pca_model, target_name = net_names[target_idx])
}

#' Dyadic (cell-level) regression on a vectorized adjacency matrix
#'
#' Regresses the vectorized off-diagonal cells of a target network on
#' ego/alter/similarity terms for every nodal attribute, reciprocity (the
#' transpose of the target network), the log number of two-paths, and terms
#' derived from every other network supplied (tie value, reciprocity,
#' ego out-degree, alter in-degree). This is the point-estimate analogue of
#' MR-QAP: the same predictor set and model, but fit once by OLS/GLM rather
#' than by permutation, since here the goal is prediction (for imputation)
#' rather than a permutation-based significance test.
#'
#' @param net_list A (preferably named) list of networks - igraph/`network`
#'   objects (treated as fully observed) or adjacency matrices (`NA` marks an
#'   unknown/missing tie). All must have the same number of nodes and node
#'   order.
#' @param attributes A data.frame/tibble of node attributes aligned to the
#'   networks' node order (or matched via `id_col`).
#' @param target Which network in `net_list` is the response - an index or
#'   (if `net_list` is named) a name.
#' @param attr_types Optional named attribute-type overrides, as in
#'   \code{\link{net_measures_core}}.
#' @param family Passed to `glm()`; default `"gaussian"`. As the authors of
#'   PMM note, PMM only uses the fitted values to rank/match donors, so the
#'   exact link/family rarely changes the imputations much - this argument
#'   matters for standalone inferential use of this function (a QAP-style
#'   tie-formation model), not for the PMM step inside \code{\link{netmice}},
#'   which always uses a linear working model regardless of tie type.
#' @param other_net_predictors "raw" (default) includes all four terms per
#'   other network; "pca" replaces them with the first `n_components`
#'   principal components, useful when many other networks are supplied.
#' @param n_components Number of PCA components when `other_net_predictors = "pca"`.
#' @param id_col Optional attribute column used to match node order.
#' @param fit Logical; if `FALSE`, only build and return the design data
#'   (no model fit) - e.g. for feeding to a custom imputation routine.
#'
#' @return A list with `data` (dyad-level data.frame including `i`, `j`
#'   indices and, when applicable, `NA` in `y` for missing ties), `pca_model`
#'   (or `NULL`), `target` (network name), and `model` (the fitted `glm`, or
#'   `NULL` if `fit = FALSE`). The model is fit only on dyads with observed
#'   `y`; predictions for all dyads (including missing ones) can be obtained
#'   with `predict(result$model, newdata = result$data)`.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 30
#' friends <- matrix(rbinom(n * n, 1, 0.1), n, n); diag(friends) <- 0
#' advice  <- matrix(rbinom(n * n, 1, 0.15), n, n); diag(advice) <- 0
#' attrs <- data.frame(age = rnorm(n, 35, 8), dept = sample(letters[1:3], n, TRUE))
#' res <- dyad_regression(list(friends = friends, advice = advice), attrs, target = "friends")
#' summary(res$model)
dyad_regression <- function(net_list, attributes, target,
                             attr_types = NULL,
                             family = "gaussian",
                             other_net_predictors = c("raw", "pca"),
                             n_components = 3,
                             id_col = NULL,
                             fit = TRUE) {
  other_net_predictors <- match.arg(other_net_predictors)

  net_names <- names(net_list)
  if (is.null(net_names)) net_names <- paste0("net", seq_along(net_list))
  names(net_list) <- net_names
  mats <- lapply(net_list, .as_matrix_generic)

  target_idx <- if (is.character(target)) match(target, net_names) else target
  if (is.na(target_idx) || length(target_idx) == 0) {
    stop("`target` not found in `net_list`.", call. = FALSE)
  }

  ref_g <- if (inherits(net_list[[1]], "igraph")) net_list[[1]] else
    igraph::graph_from_adjacency_matrix(mats[[1]] != 0, mode = "directed", diag = FALSE)
  attributes <- .align_attributes(ref_g, attributes, id_col = id_col)

  built <- .build_dyad_data(mats, attributes, target_idx, attr_types,
                             other_net_predictors, n_components)

  result <- list(data = built$data, pca_model = built$pca_model,
                  target = built$target_name, model = NULL)

  if (fit) {
    fit_rows <- !is.na(built$data$y)
    # y ~ . (rather than pasting predictor names into a formula string) so
    # that predictor names arising from factor levels with non-syntactic
    # characters (spaces, hyphens, "&", ...) don't break formula parsing.
    model_data <- built$data[fit_rows, setdiff(names(built$data), c("i", "j")), drop = FALSE]
    result$model <- stats::glm(y ~ ., data = model_data, family = family)
  }

  result
}
