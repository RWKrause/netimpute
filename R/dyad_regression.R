# Cell-level (dyadic) regression: vectorize adjacency matrices and regress a
# target network's ties on ego/alter/(diff|same) nodal-attribute terms,
# reciprocity, a two-path indicator, and terms derived from all other
# networks in the list. This is the point estimate ("MR-QAP without the permutation
# test") used both standalone and as the workhorse inside netmice().
#
# IMPORTANT: because a missing tie has no natural representation in an
# igraph/network edge list, supply networks with unknown ties as plain
# (possibly asymmetric) adjacency matrices with NA marking unknown cells.
# igraph/network objects are still accepted but are treated as fully observed.

#' @noRd
.as_matrix_generic <- function(net) {
  if (is.matrix(net)) {
    .reject_negative(net)
    return(net)
  }
  g <- .as_igraph(net)
  as.matrix(igraph::as_adjacency_matrix(
    g, sparse = FALSE,
    attr = if ("weight" %in% igraph::edge_attr_names(g)) "weight" else NULL
  ))
}

#' Validate a `structural` specification against the networks
#'
#' Accepts `NULL`, a single n x n logical matrix (recycled for every
#' network), or a named list of n x n logical matrices whose names match
#' networks in `net_names`. Returns `NULL` or a full-length named list
#' (one entry per network, `NULL` where no structural cells were given).
#' @noRd
.validate_structural <- function(structural, net_names, n) {
  if (is.null(structural)) return(NULL)
  check_one <- function(mat, label) {
    if (is.numeric(mat) && all(mat %in% c(0, 1))) mat <- mat == 1
    if (!is.matrix(mat) || !is.logical(mat)) {
      stop("`structural` (", label, ") must be a logical matrix ",
           "(TRUE = structurally absent cell).", call. = FALSE)
    }
    if (nrow(mat) != n || ncol(mat) != n) {
      stop("`structural` (", label, ") must be ", n, " x ", n,
           " to match the networks.", call. = FALSE)
    }
    if (anyNA(mat)) {
      stop("`structural` (", label, ") must not contain NA: every cell is ",
           "either structurally absent (TRUE) or not (FALSE).", call. = FALSE)
    }
    mat
  }
  if (is.matrix(structural)) {
    m <- check_one(structural, "single matrix")
    return(stats::setNames(rep(list(m), length(net_names)), net_names))
  }
  if (is.list(structural)) {
    if (is.null(names(structural)) || !all(nzchar(names(structural)))) {
      stop("When `structural` is a list, every element must be named after ",
           "a network in `networks`.", call. = FALSE)
    }
    unknown <- setdiff(names(structural), net_names)
    if (length(unknown)) {
      stop("`structural` names unknown network(s): ", toString(unknown),
           ". When `structural` is an edgelist, its `edgelist_split` must ",
           "produce the same network names as `networks` does.",
           call. = FALSE)
    }
    out <- stats::setNames(vector("list", length(net_names)), net_names)
    for (nm in names(structural)) {
      out[[nm]] <- check_one(structural[[nm]], nm)
    }
    return(out)
  }
  stop("`structural` must be NULL, one n x n logical matrix, or a named ",
       "list of n x n logical matrices.", call. = FALSE)
}

#' Fix structurally absent cells at zero in a list of adjacency matrices
#'
#' A structurally absent tie is zero by design, so an observed non-zero
#' value there contradicts the specification and errors; `NA` cells (coded
#' missing, but actually impossible) are silently set to 0 so they are never
#' treated as missing ties to impute.
#' @noRd
.apply_structural_zeros <- function(mats, struct_list) {
  if (is.null(struct_list)) return(mats)
  for (nm in names(mats)) {
    s <- struct_list[[nm]]
    if (is.null(s)) next
    m <- mats[[nm]]
    cells <- s & (row(m) != col(m))
    bad <- cells & !is.na(m) & m != 0
    if (any(bad)) {
      stop("Network '", nm, "' has ", sum(bad), " observed non-zero tie(s) ",
           "in cell(s) marked structurally absent by `structural`. ",
           "Structurally absent cells must be 0 (or NA).", call. = FALSE)
    }
    m[cells] <- 0
    mats[[nm]] <- m
  }
  mats
}

# Every endogenous statistic a tie model can carry. Which of them a given
# target actually gets is decided by .resolve_endo_terms(); the full list is
# needed in several places (protecting them from the PCA collapse, telling a
# user that a `models` formula named one this target does not have).
.ENDO_TERM_NAMES <- c("reciprocity", "twopath", "gwesp", "gwodegree",
                      "gwidegree", "gwdegree")

#' Binarize a target matrix for the endogenous statistics
#'
#' The endogenous terms are change statistics of a *binary* ERGM, so they are
#' computed on the presence/absence pattern rather than on tie weights. Before
#' 1.1.0 the two-path count was formed as `target_mat %*% target_mat` on the
#' raw matrix and only then thresholded, which happened to give the right
#' indicator for a binary target but silently gave a weighted matrix's squared
#' weights - and propagated a single `NA` across a whole row and column.
#' Unobserved cells count as absent: inside `netmice()` the matrix is always
#' filled before this is reached, but `dyad_regression()` is exported and does
#' receive matrices with `NA`.
#' @noRd
.binarize_target <- function(mat) {
  b <- (!is.na(mat) & mat != 0) * 1
  diag(b) <- 0
  b
}

#' ERGM change statistics for every dyad of a binary network
#'
#' The change statistic of cell (i, j) is the amount the network statistic
#' would move if that cell went from 0 to 1, holding every other cell fixed.
#' It is evaluated on the network with (i, j) *removed*, which is what makes
#' the degree terms leave-one-out by construction - the outcome is never a
#' summand of its own predictor - and is the reason the exponents below
#' subtract `b = B[i, j]`.
#'
#' With `w = 1 - exp(-decay)`:
#' \itemize{
#'   \item `gwodegree(i,j) = w^(outdeg_i - b)`
#'   \item `gwidegree(i,j) = w^(indeg_j - b)`
#'   \item `gwesp(i,j) = exp(a) (1 - w^ESP_ij)`
#'     `+ sum_{k: B[i,k] & B[j,k]} w^(ESP_ik - b)`
#'     `+ sum_{k: B[k,j] & B[k,i]} w^(ESP_kj - b)`
#' }
#' where `ESP` is the directed (transitive, "OTP") shared-partner count
#' `B %*% B`. The direct term is the edge's own contribution; the two sums are
#' the existing edges whose shared-partner count rises when `i -> j` appears.
#' `ESP_ij` itself needs no correction: cell (i, j) can enter `sum_k B[i,k]
#' B[k,j]` only through a diagonal entry, and the diagonal is zero.
#'
#' Both degree statistics are bounded in (0, 1]. `gwesp` is bounded by
#' `exp(decay) + 2(n - 2)`, i.e. by a quantity that grows with the network -
#' unlike the 0/1 `twopath` indicator it replaces. That is why the sweep
#' clamps its linear predictor.
#'
#' For an undirected (symmetric) `B` the same OTP expression is the undirected
#' ESP change statistic, so no separate branch is needed. The degree side does
#' differ: adding an undirected tie raises the degree of *both* endpoints, so
#' the undirected `gwdegree` is the sum of the two directed terms, not a
#' rename of either (see `.build_dyad_data()`).
#'
#' @param B binarized adjacency, zero diagonal.
#' @param decay list with `esp`, `outdegree`, `indegree`.
#' @return list of n x n matrices plus the `Tp`, `od`, `id` state the
#'   sequential sweep seeds itself from.
#' @noRd
.ergm_change_stats <- function(B, decay) {
  n <- nrow(B)
  Tp <- B %*% B
  od <- rowSums(B)
  id <- colSums(B)

  w_esp <- 1 - exp(-decay$esp)
  w_od  <- 1 - exp(-decay$outdegree)
  w_id  <- 1 - exp(-decay$indegree)

  Gwod <- w_od^(matrix(od, n, n) - B)
  Gwid <- w_id^(matrix(id, n, n, byrow = TRUE) - B)

  Wtp  <- w_esp^Tp
  # the -1 exponent is read only where B[i,j] = 1, and there the shared
  # partner being counted is reached through the (i,j) tie itself, so the
  # true exponent is >= 1; pmax() only guards entries that are never used.
  # Dividing by w_esp instead would be Inf at decay = 0.
  Wtp1 <- w_esp^pmax(Tp - 1, 0)
  M    <- B * Wtp
  M1   <- B * Wtp1
  tB   <- t(B)
  P0 <- M  %*% tB       # sum_k B[i,k] B[j,k] w^ESP_ik
  P1 <- M1 %*% tB
  Q0 <- tB %*% M        # sum_k B[k,i] B[k,j] w^ESP_kj
  Q1 <- tB %*% M1
  edge <- B == 1
  Gwesp <- exp(decay$esp) * (1 - Wtp) +
    ifelse(edge, P1, P0) + ifelse(edge, Q1, Q0)
  diag(Gwesp) <- 0
  diag(Gwod) <- 0
  diag(Gwid) <- 0

  list(gwesp = Gwesp, gwodegree = Gwod, gwidegree = Gwid,
       gwdegree = Gwod + Gwid, Tp = Tp, od = od, id = id)
}

#' Which endogenous terms a given target actually gets
#'
#' Two restrictions, both closing latent problems rather than adding policy:
#' \itemize{
#'   \item `reciprocity` is dropped for an **undirected** target, where `y_ji`
#'     *is* `y_ij`. The column is a perfect predictor there and the working
#'     model separates on it; it survived only because the fit warning was
#'     suppressed and the resulting `chol()` failure tolerated. There is
#'     precedent: a `"dyad"` random intercept already drops it.
#'   \item The `gw*` terms are dropped for a **weighted** target. They are
#'     change statistics of a binary ERGM, so on a weighted network they would
#'     describe a binarization that is not the quantity being modelled, and
#'     the PMM donor step would match a live linear predictor that can swing
#'     by O(n) against fitted values frozen at fit time.
#' }
#' A weighted or undirected target therefore keeps the bounded `twopath`
#' indicator as its closure term, which is why `twopath` is retained in the
#' package rather than removed outright.
#' @noRd
.resolve_endo_terms <- function(endo_terms, binary, undirected) {
  gw <- c("gwesp", "gwodegree", "gwidegree", "gwdegree")
  out <- endo_terms
  if (undirected) out <- setdiff(out, "reciprocity")
  if (!binary) {
    out <- setdiff(out, gw)
    if (!"twopath" %in% out && any(endo_terms %in% gw)) out <- c(out, "twopath")
  }
  if (undirected && binary) {
    # one tie, two endpoints: the undirected degree term is the sum of the
    # directed pair, not either of them renamed
    if (any(c("gwodegree", "gwidegree") %in% out)) {
      out <- setdiff(out, c("gwodegree", "gwidegree"))
      out <- c(out, "gwdegree")
    }
  } else {
    out <- setdiff(out, "gwdegree")
  }
  unique(out)
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
#' @param structural optional n x n logical matrix marking the *target*
#'   network's structurally absent (zero-by-design) cells; those dyads are
#'   dropped from the design data entirely, like the diagonal
#' @noRd
.build_dyad_data <- function(mats,
                             attributes,
                             target_idx,
                             attr_types = NULL,
                             other_net_predictors = c("raw", "pca"),
                             n_components = 3,
                             structural = NULL,
                             endo_terms = c("reciprocity", "twopath"),
                             gw_decay = list(esp = 0.69, outdegree = 0.69,
                                             indegree = 0.69),
                             undirected = FALSE,
                             binary = TRUE) {
  other_net_predictors <- match.arg(other_net_predictors)
  net_names <- names(mats)
  n <- nrow(mats[[1]])
  if (!all(vapply(mats, nrow, integer(1)) == n) || !all(vapply(mats, ncol, integer(1)) == n)) {
    stop("All networks in the list must be square matrices of the same size.",
         call. = FALSE)
  }

  target_mat <- mats[[target_idx]]
  other_mats <- mats[-target_idx]

  types <- .resolve_attr_types(attributes, attr_types)

  base <- expand.grid(i = seq_len(n), j = seq_len(n))
  keep <- base$i != base$j
  # structurally absent cells are zero by design: they are no more an
  # observation of the tie process than the diagonal is, so they are dropped
  # from the dyad data (not fit on, not predicted) - keeping them would
  # inflate the zeros in the working model. expand.grid varies `i` fastest,
  # matching as.vector()'s column-major order used for y/recip below.
  if (!is.null(structural)) keep <- keep & !as.vector(structural)
  base <- base[keep, ]

  y <- as.vector(target_mat)[keep]

  # The endogenous terms are statistics of the binary tie pattern, so they are
  # built from a binarized, NA-safe copy rather than from the raw (possibly
  # weighted) matrix. Which of them a target gets depends on whether it is
  # binary and whether it is undirected - see .resolve_endo_terms().
  B_target <- .binarize_target(target_mat)
  terms_used <- .resolve_endo_terms(endo_terms, binary = binary,
                                    undirected = undirected)
  gw_names <- c("gwesp", "gwodegree", "gwidegree", "gwdegree")
  endo_cols <- list()
  if ("reciprocity" %in% terms_used) {
    endo_cols$reciprocity <- as.vector(t(target_mat))[keep]
  }
  if ("twopath" %in% terms_used) {
    # bounded 0/1 closure indicator: "at least one shared contact"
    endo_cols$twopath <-
      as.numeric(as.vector(B_target %*% B_target)[keep] > 0)
  }
  if (any(gw_names %in% terms_used)) {
    gw_stats <- .ergm_change_stats(B_target, gw_decay)
    for (gnm in intersect(gw_names, terms_used)) {
      endo_cols[[gnm]] <- as.vector(gw_stats[[gnm]])[keep]
    }
  }

  attr_predictors <- list()
  for (nm in names(attributes)) {
    x <- attributes[[nm]]
    type <- types[[nm]]
    if (type == "continuous") {
      x <- as.numeric(x)
      ego   <- x[base$i]
      alter <- x[base$j]
      attr_predictors[[paste0(nm, "_ego")]]     <- ego
      attr_predictors[[paste0(nm, "_alter")]]   <- alter
      attr_predictors[[paste0(nm, "_absdiff")]] <- abs(ego - alter)
    } else {
      xf <- factor(x)
      # comparing the factor's integer codes: 1 = i and j hold the same
      # category, 0 = different (NA if either value is NA)
      xi <- as.integer(xf)
      same <- (xi[base$i] == xi[base$j]) * 1
      for (lv in levels(xf)[-1]) {
        dv <- as.numeric(xf == lv)
        attr_predictors[[paste0(nm, "_ego_", lv)]]   <- dv[base$i]
        attr_predictors[[paste0(nm, "_alter_", lv)]] <- dv[base$j]
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
    other_predictors[[paste0(onm, "_ego_outdeg")]]   <- outdeg_k[base$i]
    other_predictors[[paste0(onm, "_ego_indeg")]]    <- indeg_k[base$i]
    other_predictors[[paste0(onm, "_alter_indeg")]]  <- indeg_k[base$j]
    other_predictors[[paste0(onm, "_alter_outdeg")]] <- outdeg_k[base$j]
  }

  pca_model <- NULL
  if (length(other_predictors)) {
    if (other_net_predictors == "pca") {
      # The dyad-level cross-network terms (the other network's cell y_ij
      # and its transpose y_ji) are NOT collapsed into components: x_ij
      # being predicted by y_ij is exactly the cross-network association the
      # model should carry as its own, named coefficient. Only the
      # node-level degree terms are reduced to principal components.
      dyadic <- grepl("_(tie|recip)$", names(other_predictors))
      keep_list <- other_predictors[dyadic]
      om_list <- other_predictors[!dyadic]
      other_predictors <- keep_list
      if (length(om_list)) {
        om <- as.matrix(as.data.frame(om_list))
        om[is.na(om)] <- 0
        var0 <- apply(om, 2, stats::var)
        om <- om[, is.finite(var0) & var0 > 0, drop = FALSE]
        if (ncol(om) > 0) {
          k_comp <- min(n_components, ncol(om))
          pca_model <- stats::prcomp(om, center = TRUE, scale. = TRUE)
          comps <- as.data.frame(pca_model$x[, seq_len(k_comp), drop = FALSE])
          names(comps) <- paste0("other_net_PC", seq_len(k_comp))
          other_predictors <- c(other_predictors, as.list(comps))
        }
      }
    }
  }

  parts <- list(
    base,
    y = y,
    if (length(endo_cols)) as.data.frame(endo_cols, check.names = FALSE) else NULL,
    as.data.frame(attr_predictors, check.names = FALSE),
    if (length(other_predictors)) as.data.frame(other_predictors,
                                                check.names = FALSE) else NULL
  )
  data <- do.call(cbind, Filter(Negate(is.null), parts))
  rownames(data) <- NULL

  list(data = data, pca_model = pca_model, target_name = net_names[target_idx])
}

#' Dyadic (cell-level) regression on a vectorized adjacency matrix
#'
#' Regresses the vectorized off-diagonal cells of a target network on
#' ego/alter/similarity terms for every nodal attribute, the target's own
#' endogenous terms (see `endo_terms` - by default `reciprocity` plus the
#' geometrically weighted shared-partner and degree change statistics), and
#' terms
#' derived from every other network supplied (tie value, reciprocity, and
#' the sender's and receiver's out- and in-degrees in that network). This
#' is the point-estimate analogue of
#' MR-QAP: the same predictor set and model, but fit once by OLS/GLM rather
#' than by permutation, since here the goal is prediction (for imputation)
#' rather than a permutation-based significance test.
#'
#' @section Term names:
#' The returned `data` (and hence the fitted model) contains, besides the
#' `i`/`j` indices and the response `y`, exactly these predictor columns -
#' the same names a \code{\link{netmice}} `models` formula for a network can
#' reference:
#' \describe{
#'   \item{\code{reciprocity}}{the target's transpose cell - the tie value
#'     from j back to i.}
#'   \item{\code{gwesp}, \code{gwodegree}, \code{gwidegree}}{for a binary
#'     target: the geometrically weighted shared-partner and degree **change
#'     statistics** - how much each ERGM statistic would move if this cell
#'     went from 0 to 1, evaluated on the network with the cell removed. An
#'     undirected binary target gets \code{gwesp} and a single summed
#'     \code{gwdegree} instead, and no \code{reciprocity}.}
#'   \item{\code{twopath}}{0/1 indicator for at least one two-path
#'     i -> k -> j (at least one shared contact). The default closure term
#'     before 1.1.0, and still the one used for a \emph{weighted} target,
#'     where the binary-ERGM \code{gw*} statistics do not apply.}
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
#' @param networks A (preferably named) list of networks - igraph/`network`
#'   objects (treated as fully observed) or adjacency matrices (`NA` marks an
#'   unknown/missing tie). All must have the same number of nodes and node
#'   order.
#' @param attributes A data.frame/tibble of node attributes aligned to the
#'   networks' node order (or matched via `id_col`).
#' @param target Which network in `networks` is the response - an index or
#'   (if `networks` is named) a name.
#' @param attr_types Optional named attribute-type overrides, as in
#'   \code{\link{net_measures_core}}.
#' @param family Passed to `glm()`; default `"gaussian"`. As the authors of
#'   PMM note, PMM only uses the fitted values to rank/match donors, so the
#'   exact link/family rarely changes the imputations much - this argument
#'   matters for standalone inferential use of this function (a QAP-style
#'   tie-formation model), not for the PMM step inside \code{\link{netmice}},
#'   which always uses a linear working model regardless of tie type.
#' @param other_net_predictors "raw" (default) includes all six terms per
#'   other network; "pca" keeps the dyad-level cross-network terms
#'   (`*_tie`, `*_recip`) as raw predictors and replaces only the node-level
#'   degree terms (`*_ego_outdeg`, `*_ego_indeg`, `*_alter_indeg`,
#'   `*_alter_outdeg`) with the first
#'   principal components allowed by `PCA`, useful when many other networks
#'   are supplied. The cross-network cell terms are never absorbed into
#'   components: the target's `x_ij` being predicted by the other network's
#'   `y_ij` is a substantive association that should keep its own
#'   coefficient.
#' @param PCA Predictor budget when `other_net_predictors = "pca"`: a list
#'   with `n` (a fixed maximum number of components) and/or `ratio` (observed
#'   rows per predictor, or observed *events* - the rarer of ties/non-ties -
#'   when the target is binary). The default caps components at 3 as before
#'   and additionally refuses to exceed a 10:1 budget.
#' @param structural `NULL` (default), one n x n logical matrix, or a named
#'   list of n x n logical matrices (names matching networks in `networks`).
#'   `TRUE` marks a cell whose tie is *structurally absent* - zero by design
#'   (e.g. ties into another school class that respondents could not
#'   nominate) rather than a genuine observation of "no tie". A single matrix
#'   applies the same structural cells to every network; a named list lets
#'   them differ per network (networks without an entry have none).
#'   Structural cells must hold 0 or `NA` in the supplied networks (a
#'   non-zero observed value there errors); they are fixed at 0 and the
#'   *target* network's structural dyads are removed from the returned
#'   `data` - and hence from the regression - like the diagonal, so they
#'   cannot inflate the zeros of the fitted model.
#' @param endo_terms Endogenous (network-derived) terms to include, any of
#'   `"reciprocity"`, `"twopath"`, `"gwesp"`, `"gwodegree"`, `"gwidegree"`.
#'   Default `c("reciprocity", "twopath")`; use
#'   `c("reciprocity", "gwesp", "gwodegree", "gwidegree")` for the
#'   geometrically weighted change statistics, which are implemented and
#'   tested but not the default - see `?netmice`. The set is
#'   filtered by what the target supports - a weighted target gets `twopath`
#'   rather than the `gw*` binary-ERGM statistics, and an undirected one has
#'   no `reciprocity` (there `y_ji` is `y_ij`, which separates the fit) and a
#'   single summed `gwdegree`. See `?netmice` for the full table.
#' @param gw_decay Decay for the `gw*` statistics: a list with `esp`,
#'   `outdegree`, `indegree`, or one number for all three (default `0.69`).
#' @param random_intercepts `NULL` (default) fits an ordinary `glm()`.
#'   Otherwise a character vector - any of `"ego"`, `"alter"`, `"dyad"` - and
#'   the model is fit with \pkg{lme4} instead, adding a random intercept per
#'   sender node (`"ego"`, i.e. `(1 | i)`), per receiver node (`"alter"`,
#'   `(1 | j)`), and/or per unordered node pair (`"dyad"`). Ego and alter
#'   intercepts capture actor-level activity/popularity heterogeneity beyond
#'   the degree terms, as in the social relations model; the dyad intercept
#'   is the SRM-style relationship effect and only makes sense for *directed*
#'   networks, where the ordered cells (i,j) and (j,i) are two genuine
#'   observations of the pair. For an undirected target those two rows are
#'   exact duplicates of one another, so a dyad intercept has effectively one
#'   unique observation per group and degenerates into a residual term - a
#'   warning is issued and ego/alter intercepts are recommended instead.
#'   When `"dyad"` is included, the fixed `reciprocity` term is dropped from
#'   the model: it is the transpose of `y`, so together with a per-dyad
#'   intercept the model could reproduce `y` exactly (exact singularity) -
#'   in SRM terms the dyad effect *is* the reciprocity/relationship term and
#'   replaces it. Gaussian-identity models use `lme4::lmer()`; any other
#'   `family` uses `lme4::glmer()`.
#' @param id_col Optional attribute column used to match node order.
#' @param fit Logical; if `FALSE`, only build and return the design data
#'   (no model fit) - e.g. for feeding to a custom imputation routine.
#'
#' @return A list with `data` (dyad-level data.frame including `i`, `j`
#'   indices and, when applicable, `NA` in `y` for missing ties), `pca_model`
#'   (or `NULL`), `target` (network name), and `model` (the fitted `glm` -
#'   or `merMod` when `random_intercepts` is used - or `NULL` if
#'   `fit = FALSE`). The model is fit only on dyads with observed `y`;
#'   predictions for all dyads (including missing ones) can be obtained with
#'   `predict(result$model, newdata = result$data)` (for a mixed model, the
#'   grouping columns `.ego`/`.alter`/`.dyad` are already present in `data`,
#'   and `allow.new.levels = TRUE` may be needed if some node never appears
#'   in an observed dyad). Dyads marked structurally absent via `structural`
#'   do not appear in `data` at all.
#' @seealso \code{\link{netmice}} for the joint imputation loop this
#'   dyad-level model is the tie engine of; \code{\link{netquickpred}} to
#'   screen its predictor set; \code{\link{net_measures}} for the node-level
#'   measures that complement these dyad-level terms.
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
dyad_regression <- function(networks,
                            attributes,
                            target,
                            attr_types = NULL,
                            family = "gaussian",
                            other_net_predictors = c("raw", "pca"),
                            PCA = list(n = 3, ratio = 10),
                            structural = NULL,
                            endo_terms = c("reciprocity", "twopath"),
                            gw_decay = list(esp = 0.69, outdegree = 0.69,
                                            indegree = 0.69),
                            random_intercepts = NULL,
                            id_col = NULL,
                            fit = TRUE) {
  net_list <- networks   # internal name, kept to avoid churning the body
  other_net_predictors <- match.arg(other_net_predictors)
  PCA <- .validate_pca(PCA)
  if (!is.null(random_intercepts)) {
    random_intercepts <- match.arg(random_intercepts,
                                   c("ego", "alter", "dyad"),
                                   several.ok = TRUE)
  }

  net_names <- names(net_list)
  if (is.null(net_names)) net_names <- paste0("net", seq_along(net_list))
  names(net_list) <- net_names
  mats <- lapply(net_list, .as_matrix_generic)

  target_idx <- if (is.character(target)) match(target, net_names) else target
  if (is.na(target_idx) || length(target_idx) == 0) {
    stop("`target` not found in `networks`.", call. = FALSE)
  }

  struct_list <- .validate_structural(structural, net_names, nrow(mats[[1]]))
  mats <- .apply_structural_zeros(mats, struct_list)

  # reference graph only supplies node count/names for attribute alignment;
  # unknown (NA) ties must not break its construction
  ref_m <- mats[[1]] != 0
  ref_m[is.na(ref_m)] <- FALSE
  ref_g <- if (inherits(net_list[[1]], "igraph")) net_list[[1]] else
    igraph::graph_from_adjacency_matrix(ref_m, mode = "directed",
                                        diag = FALSE)
  attributes <- .align_attributes(ref_g, attributes, id_col = id_col)

  # budget from the target's own observed cells: for a binary network the
  # rarer of ties/non-ties, otherwise the observed row count
  tgt_mat <- mats[[target_idx]]
  tgt_struct <- struct_list[[net_names[target_idx]]]
  obs_mask <- !is.na(tgt_mat) & (row(tgt_mat) != col(tgt_mat))
  if (!is.null(tgt_struct)) obs_mask <- obs_mask & !tgt_struct
  tgt_obs <- tgt_mat[obs_mask]
  n_comp <- .pca_budget(PCA, .pca_denom(
    tgt_obs, NULL, binary = all(tgt_obs[!is.na(tgt_obs)] %in% c(0, 1))))

  tgt_bin <- all(tgt_obs[!is.na(tgt_obs)] %in% c(0, 1))
  built <- .build_dyad_data(mats, attributes, target_idx, attr_types,
                             other_net_predictors, n_comp,
                             structural = struct_list[[net_names[target_idx]]],
                             endo_terms = .validate_endo_terms(endo_terms),
                             gw_decay = .validate_gw_decay(gw_decay),
                             undirected = isSymmetric(unname(
                               mats[[target_idx]])),
                             binary = tgt_bin)

  if (length(random_intercepts)) {
    built$data <- .add_dyad_groups(built$data)
  }

  result <- list(data = built$data, pca_model = built$pca_model,
                  target = built$target_name, model = NULL)

  if (fit) {
    fit_rows <- !is.na(built$data$y)
    group_cols <- c(".ego", ".alter", ".dyad")
    model_data <- built$data[fit_rows,
                             setdiff(names(built$data),
                                     c("i", "j", group_cols)),
                             drop = FALSE]

    if (length(random_intercepts)) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' is required for `random_intercepts`. ",
             "Install with install.packages('lme4').", call. = FALSE)
      }
      if ("dyad" %in% random_intercepts &&
          isSymmetric(unname(mats[[target_idx]]))) {
        warning("The target network is undirected, so each dyad's two rows ",
                "(i,j)/(j,i) are exact duplicates: a dyad random intercept ",
                "has effectively one unique observation per group and will ",
                "absorb the residual variance. Consider random_intercepts = ",
                "c('ego', 'alter') instead.", call. = FALSE)
      }
      md <- cbind(model_data,
                  built$data[fit_rows, group_cols, drop = FALSE])
      fixed <- setdiff(names(model_data), "y")
      if ("dyad" %in% random_intercepts) {
        # The fixed `reciprocity` term is the transpose of y, so together
        # with a per-dyad intercept the model can reproduce y exactly
        # (b_d = y_ij + y_ji with coefficient -1 on reciprocity): the
        # residual variance collapses and lme4 aborts. In SRM terms the
        # dyad effect *is* the reciprocity/relationship term, so it
        # replaces the fixed effect rather than coexisting with it.
        fixed <- setdiff(fixed, "reciprocity")
        message("netimpute: dropping the fixed `reciprocity` term - it is ",
                "absorbed by (and singular with) the dyad random intercept.")
      }
      re_terms <- c(ego = "(1 | .ego)", alter = "(1 | .alter)",
                    dyad = "(1 | .dyad)")[random_intercepts]
      fml <- stats::reformulate(
        c(sprintf("`%s`", fixed), re_terms),
        response = "y"
      )
      fam <- family
      if (is.character(fam)) fam <- get(fam, mode = "function")
      if (is.function(fam)) fam <- fam()
      result$model <- if (fam$family == "gaussian" && fam$link == "identity") {
        lme4::lmer(fml, data = md)
      } else {
        lme4::glmer(fml, data = md, family = fam)
      }
    } else {
      result$model <- stats::glm(y ~ ., data = model_data, family = family)
    }
  }

  result
}

#' Append `.ego`, `.alter`, `.dyad` grouping factors (from the `i`/`j`
#' indices) to a dyad-level data.frame, for lme4 random intercepts.
#' @noRd
.add_dyad_groups <- function(d) {
  d$.ego   <- factor(d$i)
  d$.alter <- factor(d$j)
  d$.dyad  <- factor(paste0(pmin(d$i, d$j), "_", pmax(d$i, d$j)))
  d
}
