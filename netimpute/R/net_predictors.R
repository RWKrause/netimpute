#' Build node-level predictors (raw measures and/or PCA components) across a
#' list of networks
#'
#' Applies \code{\link{net_measures}} (with the requested `measure_set`)
#' to each network/attribute pair in `net_list`/`attr_list`, stacks the
#' results, and optionally reduces the stacked measures to principal
#' components suitable as auxiliary predictors in a multiple-imputation model
#' (e.g. `mice` with `method = "pmm"`).
#'
#' This is a single entry point with three modes, controlled by `output`:
#' \describe{
#'   \item{"measures"}{Return the stacked raw measures only (no PCA).}
#'   \item{"pca"}{Return principal components computed from the stacked
#'     measures (all of them, or all minus `keep_vars` if supplied).}
#'   \item{"both"}{Return the PCA components *and* the raw columns named in
#'     `keep_vars`, un-transformed. Use this when a specific measure (e.g.
#'     degree) is itself part of a substantive hypothesis and must not be
#'     diluted/absorbed into a composite component - the PCA is fit on the
#'     remaining measures only (optionally residualized on `keep_vars`, see
#'     `residualize_kept`), and `keep_vars` are appended untouched.}
#' }
#'
#' Because the auxiliary predictors built here are meant to feed a subsequent
#' imputation model rather than to be interpreted themselves, missing values
#' in the measure matrix (e.g. `alter_mean` for isolates, or
#' `reciprocity_ratio` in an undirected network) are mean-imputed column-wise
#' before PCA by default (see `impute_na`). This single mean-imputation of
#' auxiliary variables prior to PCA is standard practice and is far less
#' consequential than mean-imputing a substantive analysis variable.
#'
#' @param net_list A list of networks (igraph / network / adjacency matrix),
#'   optionally named - names are carried through as `network_id`.
#' @param attr_list A list of attribute data.frames, same length and order as
#'   `net_list`, sharing the same column names/types ("the same attribute set").
#' @param attr_types Optional named character vector overriding auto-detected
#'   attribute types, applied to every network.
#' @param measure_set Which network measures to compute, passed to
#'   \code{\link{net_measures}}: any combination of "core", "full",
#'   "homophily", and individual measure names (named sets and single
#'   measures are unioned, e.g. \code{c("core", "total_degree")}). Default
#'   "core".
#' @param output One of "measures", "pca", "both" - see Details.
#' @param n_components Number of principal components to keep (default 5).
#' @param keep_vars Character vector of raw measure column names to preserve
#'   untransformed. Required (non-empty) when `output = "both"`.
#' @param residualize_kept Logical (default `FALSE`). If `TRUE`, every
#'   variable entering the PCA is first residualized on `keep_vars` (via
#'   `lm()`), so the components capture variance orthogonal to the
#'   preserved predictors instead of variance that overlaps with them.
#' @param impute_na How to handle missing values in the measure matrix before
#'   PCA: "mean" (default) or "complete_cases" (drop rows with any NA).
#' @param scale_pca Logical (default `TRUE`); center and scale measures
#'   before `prcomp()`. Strongly recommended given the very different scales
#'   of counts, ratios, and continuous differences produced here.
#' @param id_col,use_sna Passed through to the underlying measure function.
#'
#' @return
#' If `output = "measures"`: a data.frame.
#' If `output = "pca"` or `"both"`: a list with
#' \describe{
#'   \item{predictors}{data.frame of `network_id`, `node_id`, PCs (and
#'     `keep_vars` if `output = "both"`) - ready to hand to an imputation model.}
#'   \item{pca_model}{the `prcomp` object, for scree/loadings diagnostics.}
#'   \item{raw_measures}{the full stacked measures data.frame, for reference.}
#' }
#' @export
#'
#' @examples
#' set.seed(1)
#' nets <- list(
#'   school_a = igraph::sample_gnp(30, 0.08, directed = TRUE),
#'   school_b = igraph::sample_gnp(35, 0.06, directed = TRUE)
#' )
#' attrs <- lapply(c(30, 35), function(n) {
#'   data.frame(age = round(rnorm(n, 35, 8)),
#'              gender = sample(c("F", "M"), n, replace = TRUE))
#' })
#' # keep total_degree raw, reduce the remaining measures to components
#' res <- net_predictors(nets, attrs, measure_set = c("core", "total_degree"),
#'                        output = "both", keep_vars = c("total_degree"))
#' head(res$predictors)
net_predictors <- function(net_list,
                           attr_list,
                           attr_types = NULL,
                           measure_set = "core",
                           output = c("measures", "pca", "both"),
                           n_components = 5,
                           keep_vars = NULL,
                           residualize_kept = FALSE,
                           impute_na = c("mean", "complete_cases"),
                           scale_pca = TRUE,
                           id_col = NULL,
                           use_sna = TRUE) {
  measure_set <- .resolve_measure_set(measure_set)
  output <- match.arg(output)
  impute_na <- match.arg(impute_na)

  if (length(net_list) != length(attr_list)) {
    stop("`net_list` and `attr_list` must have the same length.", call. = FALSE)
  }
  if (output == "both" && length(keep_vars) == 0) {
    stop("`keep_vars` must be supplied (non-empty) when output = 'both'.",
         call. = FALSE)
  }

  ref_cols <- names(attr_list[[1]])
  if (!all(vapply(attr_list, function(a) identical(names(a), ref_cols), logical(1)))) {
    stop("All data.frames in `attr_list` must share the same column names ",
         "('the same attribute set').", call. = FALSE)
  }

  net_ids <- names(net_list)
  if (is.null(net_ids)) net_ids <- as.character(seq_along(net_list))

  measures_list <- Map(function(net, attrs, nid) {
    out <- net_measures(net, attrs, measure_set = measure_set,
                        attr_types = attr_types, id_col = id_col,
                        use_sna = use_sna)
    out$network_id <- nid
    out
  }, net_list, attr_list, net_ids)

  all_cols <- Reduce(union, lapply(measures_list, names))
  measures_list <- lapply(measures_list, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (m in missing) df[[m]] <- NA
    df[all_cols]
  })
  measures_df <- do.call(rbind, measures_list)
  rownames(measures_df) <- NULL
  front <- c("network_id", "node_id")
  measures_df <- measures_df[c(front, setdiff(names(measures_df), front))]

  if (output == "measures") return(measures_df)

  if (length(keep_vars)) {
    missing_keep <- setdiff(keep_vars, names(measures_df))
    if (length(missing_keep)) {
      stop("keep_vars not found among computed measures: ",
           toString(missing_keep), call. = FALSE)
    }
  }

  pca_input <- .prep_pca_matrix(measures_df, exclude = c(front, keep_vars))

  # NA-handling always happens BEFORE residualizing on keep_vars: lm()/
  # residuals() drop incomplete-case rows per column by default, so if
  # different pca_input columns had different NA patterns, residualizing
  # each one separately would return residual vectors of different lengths
  # and apply() could no longer reassemble them into a matrix.
  if (impute_na == "mean") {
    pca_input <- apply(pca_input, 2, function(col) {
      if (all(is.na(col))) return(col)
      col[is.na(col)] <- mean(col, na.rm = TRUE)
      col
    })
    keep_rows <- rep(TRUE, nrow(pca_input))
  } else {
    keep_rows <- stats::complete.cases(pca_input)
    if (!all(keep_rows)) {
      message("netimpute: dropping ", sum(!keep_rows),
              " row(s) with missing measures (impute_na = 'complete_cases').")
    }
    pca_input <- pca_input[keep_rows, , drop = FALSE]
  }

  if (residualize_kept && length(keep_vars)) {
    keep_mat <- as.matrix(measures_df[keep_rows, keep_vars, drop = FALSE])
    keep_mat <- apply(keep_mat, 2, function(col) {
      col[is.na(col)] <- mean(col, na.rm = TRUE)
      col
    })
    pca_input <- apply(pca_input, 2, function(col) {
      stats::residuals(stats::lm(col ~ keep_mat))
    })
  }

  var0 <- apply(pca_input, 2, function(col) stats::var(col, na.rm = TRUE))
  bad <- names(var0)[is.na(var0) | var0 == 0]
  if (length(bad)) {
    message(
      "netimpute: dropping zero-variance / all-NA measure(s) before PCA: ",
            toString(bad))
    pca_input <- pca_input[, setdiff(colnames(pca_input), bad), drop = FALSE]
  }

  n_components <- min(n_components, ncol(pca_input))
  pca_fit <- stats::prcomp(pca_input, center = TRUE, scale. = scale_pca)
  comps <- as.data.frame(pca_fit$x[, seq_len(n_components), drop = FALSE])
  names(comps) <- paste0("PC", seq_len(n_components))

  base_df <- measures_df[keep_rows, front, drop = FALSE]
  if (output == "both") {
    predictors <- cbind(base_df,
                        measures_df[keep_rows, keep_vars, drop = FALSE],
                        comps)
  } else {
    predictors <- cbind(base_df, comps)
  }
  rownames(predictors) <- NULL

  list(predictors = predictors, pca_model = pca_fit, raw_measures = measures_df)
}

#' Coerce a stacked measures data.frame into a numeric matrix suitable for PCA:
#' logicals become 0/1, character/factor columns become dummy indicators
#' (first level dropped).
#' @noRd
.prep_pca_matrix <- function(df, exclude = character(0)) {
  df <- df[setdiff(names(df), exclude)]
  df[] <- lapply(df, function(col) if (is.logical(col)) as.numeric(col) else col)

  cat_cols <- names(df)[vapply(df, function(col) is.character(col) || is.factor(col), logical(1))]
  num_df <- df[setdiff(names(df), cat_cols)]

  if (length(cat_cols)) {
    dummy_list <- lapply(cat_cols, function(cn) {
      f <- factor(df[[cn]])
      if (nlevels(f) < 2) return(NULL)
      mf <- stats::model.frame(~f - 1, data = data.frame(f = f),
                               na.action = stats::na.pass)
      mm <- stats::model.matrix(~f - 1, data = mf)
      lvls <- levels(f)[-1]
      mm <- mm[, paste0("f", lvls), drop = FALSE]
      colnames(mm) <- paste0(cn, "_", lvls)
      mm
    })
    dummy_list <- Filter(Negate(is.null), dummy_list)
    if (length(dummy_list)) num_df <- cbind(num_df, do.call(cbind, dummy_list))
  }

  as.matrix(num_df)
}
