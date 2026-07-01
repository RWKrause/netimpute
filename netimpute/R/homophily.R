# Node-level, attribute-scale-appropriate homophily/alter-similarity measures.
# Shared by net_measures_core() and net_measures_full() so both functions use
# exactly the same definitions.

#' E-I index per node: (external - internal ties) / total ties
#' Negative = homophilous (mostly ties to same category), positive = heterophilous.
#' @keywords internal
.ei_index <- function(g, x) {
  n <- igraph::vcount(g)
  vapply(seq_len(n), function(i) {
    nb <- as.integer(igraph::neighbors(g, i, mode = "all"))
    if (length(nb) == 0) return(NA_real_)
    same <- x[i] == x[nb]
    internal <- sum(same, na.rm = TRUE)
    external <- sum(!same, na.rm = TRUE)
    (external - internal) / length(nb)
  }, numeric(1))
}

#' Blau's index of heterogeneity of alters' attribute values: 1 - sum(p_k^2)
#' @keywords internal
.blau_index <- function(g, x) {
  n <- igraph::vcount(g)
  vapply(seq_len(n), function(i) {
    nb <- as.integer(igraph::neighbors(g, i, mode = "all"))
    if (length(nb) == 0) return(NA_real_)
    p <- table(x[nb]) / length(nb)
    1 - sum(p^2)
  }, numeric(1))
}

#' Modal category among a node's alters
#' @keywords internal
.alter_mode <- function(g, x) {
  n <- igraph::vcount(g)
  vapply(seq_len(n), function(i) {
    nb <- as.integer(igraph::neighbors(g, i, mode = "all"))
    if (length(nb) == 0) return(NA_character_)
    tab <- table(x[nb])
    names(tab)[which.max(tab)]
  }, character(1))
}

#' Continuous-attribute alter statistics: directional absolute difference plus
#' mean/min/max of alters' values.
#' @keywords internal
.alter_continuous_stats <- function(g, x) {
  n <- igraph::vcount(g)
  directed <- igraph::is_directed(g)
  rows <- lapply(seq_len(n), function(i) {
    nb_out <- as.integer(igraph::neighbors(g, i, mode = if (directed) "out" else "all"))
    nb_in  <- as.integer(igraph::neighbors(g, i, mode = if (directed) "in"  else "all"))
    nb_all <- as.integer(igraph::neighbors(g, i, mode = "all"))
    c(
      diff_out   = if (length(nb_out)) mean(abs(x[i] - x[nb_out]), na.rm = TRUE) else NA_real_,
      diff_in    = if (length(nb_in))  mean(abs(x[i] - x[nb_in]),  na.rm = TRUE) else NA_real_,
      alter_mean = if (length(nb_all)) mean(x[nb_all], na.rm = TRUE) else NA_real_,
      alter_min  = if (length(nb_all)) min(x[nb_all],  na.rm = TRUE) else NA_real_,
      alter_max  = if (length(nb_all)) max(x[nb_all],  na.rm = TRUE) else NA_real_
    )
  })
  as.data.frame(do.call(rbind, rows))
}

#' Build the homophily/alter-similarity block for every attribute column,
#' using the measure appropriate to that attribute's declared/inferred scale:
#'  - binary / multinomial -> E-I index, Blau's index, modal alter category
#'  - continuous           -> avg |diff| on outgoing ties, on incoming ties,
#'                             and mean/min/max of alters' values
#'
#' @param g igraph object
#' @param attributes data.frame of node attributes, aligned to V(g)
#' @param attr_types named character vector overriding auto-detected types
#' @keywords internal
.homophily_block <- function(g, attributes, attr_types = NULL) {
  types <- .resolve_attr_types(attributes, attr_types)
  blocks <- lapply(names(attributes), function(nm) {
    x <- attributes[[nm]]
    type <- types[[nm]]
    if (type == "continuous") {
      x <- as.numeric(x)
      df <- .alter_continuous_stats(g, x)
      names(df) <- paste0(nm, "_", c("diff_out", "diff_in", "alter_mean", "alter_min", "alter_max"))
    } else {
      x <- as.character(x)
      df <- data.frame(
        ei_index   = .ei_index(g, x),
        blau_index = .blau_index(g, x),
        alter_mode = .alter_mode(g, x),
        stringsAsFactors = FALSE
      )
      names(df) <- paste0(nm, "_", c("ei_index", "blau_index", "alter_mode"))
    }
    df
  })
  do.call(cbind, blocks)
}
