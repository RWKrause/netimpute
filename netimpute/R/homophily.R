# Node-level, attribute-scale-appropriate homophily/alter-similarity measures.
# This is the "homophily" block of net_measures(), so every measure set uses
# exactly the same definitions.

#' Integer neighbor lists for all vertices in one igraph call
#'
#' Per-node igraph::neighbors() calls carry heavy per-call argument-checking
#' overhead and were the single largest cost in netmice()'s profile (the
#' homophily block calls them per node x per measure x per attribute).
#' as_adj_list() returns the identical neighbor sets - including counting a
#' mutual tie twice under mode = "all", matching neighbors() - in one call.
#' @keywords internal
.adj_int_list <- function(g, mode) {
  lapply(igraph::as_adj_list(g, mode = mode), as.integer)
}

#' All three neighbor-list variants (out / in / all) for a graph
#'
#' For undirected graphs the out- and in-lists are the all-list itself (no
#' copies), mirroring how neighbors(mode = "out"/"in") behaves there.
#' @keywords internal
.neighbor_lists <- function(g) {
  all_nb <- .adj_int_list(g, "all")
  if (igraph::is_directed(g)) {
    list(out = .adj_int_list(g, "out"), `in` = .adj_int_list(g, "in"), all = all_nb)
  } else {
    list(out = all_nb, `in` = all_nb, all = all_nb)
  }
}

#' E-I index per node: (external - internal ties) / total ties
#' Negative = homophilous (mostly ties to same category), positive = heterophilous.
#' @keywords internal
.ei_index <- function(g, x, nb_all = NULL) {
  if (is.null(nb_all)) nb_all <- .adj_int_list(g, "all")
  vapply(seq_along(nb_all), function(i) {
    nb <- nb_all[[i]]
    if (length(nb) == 0) return(NA_real_)
    same <- x[i] == x[nb]
    internal <- sum(same, na.rm = TRUE)
    external <- sum(!same, na.rm = TRUE)
    (external - internal) / length(nb)
  }, numeric(1))
}

#' Blau's index of heterogeneity of alters' attribute values: 1 - sum(p_k^2)
#' @keywords internal
.blau_index <- function(g, x, nb_all = NULL) {
  if (is.null(nb_all)) nb_all <- .adj_int_list(g, "all")
  # tabulate() on the factor's integer codes replaces a per-node table()
  # call; like table(), it ignores NA values, and dividing by length(nb)
  # (NA alters included) matches the previous definition exactly.
  xf <- factor(x)
  xi <- as.integer(xf)
  k <- nlevels(xf)
  vapply(nb_all, function(nb) {
    if (length(nb) == 0) return(NA_real_)
    p <- tabulate(xi[nb], nbins = k) / length(nb)
    1 - sum(p^2)
  }, numeric(1))
}

#' Modal category among a node's alters
#' @keywords internal
.alter_mode <- function(g, x, nb_all = NULL) {
  if (is.null(nb_all)) nb_all <- .adj_int_list(g, "all")
  # factor levels are sorted like table()'s dimnames, so which.max() breaks
  # ties toward the alphabetically first category exactly as before.
  xf <- factor(x)
  xi <- as.integer(xf)
  lv <- levels(xf)
  vapply(nb_all, function(nb) {
    if (length(nb) == 0) return(NA_character_)
    cnt <- tabulate(xi[nb], nbins = length(lv))
    if (all(cnt == 0L)) return(NA_character_)  # every alter's value is NA
    lv[which.max(cnt)]
  }, character(1))
}

#' Continuous-attribute alter statistics: directional absolute difference plus
#' mean/min/max of alters' values, and the directional alter means
#' (over incoming and outgoing ties separately). The directional means are
#' pure alter-composition measures - unlike diff_in/diff_out they do not use
#' ego's own value, so they remain valid predictors even for the attribute
#' being imputed; the incoming-tie mean is especially valuable there, since
#' incoming ties are typically reported by (observed) others.
#' @keywords internal
.alter_continuous_stats <- function(g, x, nb_lists = NULL) {
  if (is.null(nb_lists)) nb_lists <- .neighbor_lists(g)
  rows <- lapply(seq_along(nb_lists$all), function(i) {
    nb_out <- nb_lists$out[[i]]
    nb_in  <- nb_lists$`in`[[i]]
    nb_all <- nb_lists$all[[i]]
    # min/max (unlike mean) return +/-Inf with a warning when every alter's
    # value is NA - e.g. screening on pre-imputation data - so guard on
    # having at least one non-NA alter value
    any_all <- any(!is.na(x[nb_all]))
    c(
      diff_out       = if (length(nb_out)) mean(abs(x[i] - x[nb_out]),
                                                na.rm = TRUE) else NA_real_,
      diff_in        = if (length(nb_in))  mean(abs(x[i] - x[nb_in]),
                                                na.rm = TRUE) else NA_real_,
      alter_mean     = if (length(nb_all)) mean(x[nb_all],
                                                na.rm = TRUE) else NA_real_,
      alter_mean_in  = if (length(nb_in))  mean(x[nb_in],
                                                na.rm = TRUE) else NA_real_,
      alter_mean_out = if (length(nb_out)) mean(x[nb_out],
                                                na.rm = TRUE) else NA_real_,
      alter_min      = if (any_all) min(x[nb_all], na.rm = TRUE) else NA_real_,
      alter_max      = if (any_all) max(x[nb_all], na.rm = TRUE) else NA_real_
    )
  })
  as.data.frame(do.call(rbind, rows))
}

#' Build the homophily/alter-similarity block for every attribute column,
#' using the measure appropriate to that attribute's declared/inferred scale:
#'  - binary / multinomial -> E-I index, Blau's index, modal alter category
#'  - continuous           -> avg |diff| on outgoing ties, on incoming ties,
#'                             mean/min/max of alters' values, and the alter
#'                             mean over incoming and outgoing ties separately
#'
#' @param g igraph object
#' @param attributes data.frame of node attributes, aligned to V(g)
#' @param attr_types named character vector overriding auto-detected types
#' @keywords internal
.homophily_block <- function(g, attributes, attr_types = NULL) {
  types <- .resolve_attr_types(attributes, attr_types)
  # neighbor lists are attribute-independent: build them once here and reuse
  # across every attribute and measure, rather than per node inside each
  # helper (netmice() calls this block per network at every attribute visit)
  nb <- .neighbor_lists(g)
  blocks <- lapply(names(attributes), function(nm) {
    x <- attributes[[nm]]
    type <- types[[nm]]
    if (type == "continuous") {
      x <- as.numeric(x)
      df <- .alter_continuous_stats(g, x, nb_lists = nb)
      names(df) <- paste0(nm, "_", c("diff_out", "diff_in", "alter_mean",
                                     "alter_mean_in", "alter_mean_out",
                                     "alter_min", "alter_max"))
    } else {
      x <- as.character(x)
      df <- data.frame(
        ei_index   = .ei_index(g, x, nb_all = nb$all),
        blau_index = .blau_index(g, x, nb_all = nb$all),
        alter_mode = .alter_mode(g, x, nb_all = nb$all),
        stringsAsFactors = FALSE
      )
      names(df) <- paste0(nm, "_", c("ei_index", "blau_index", "alter_mode"))
    }
    df
  })
  do.call(cbind, blocks)
}
