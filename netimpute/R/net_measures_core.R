#' Core node-level network measures
#'
#' Computes a compact, low-collinearity set of node-level measures
#' (see Details) plus attribute-scale-appropriate homophily/alter-similarity
#' measures for every column in `attributes`.
#'
#' Structural measures returned:
#' \itemize{
#'   \item \code{outdegree}, \code{indegree}
#'   \item \code{reciprocity_ratio} (NA for undirected networks; for binary
#'         networks the usual mutual/either dyad-count ratio in \[0, 1\],
#'         higher = more reciprocal; for non-negative weighted networks the
#'         average \code{|w_ij - w_ji|} over alters tied in either direction -
#'         in the original weight units, *lower* = more reciprocal)
#'   \item \code{bonacich_power} (Bonacich power centrality, auto-scaled exponent)
#'   \item \code{betweenness} (normalized)
#'   \item \code{isolate} (logical)
#'   \item \code{constraint} (Burt's constraint)
#'   \item \code{harmonic_closeness}
#'   \item \code{local_clustering} (local transitivity for binary networks;
#'         for non-negative weighted networks, the Opsahl-Panzarasa (2009)
#'         geometric-mean weighted clustering coefficient)
#' }
#'
#' Only binary (0/1) and non-negative weighted networks are currently
#' supported; a network with negative tie values (signed networks) raises an
#' error.
#'
#' For each column of `attributes`, homophily columns are added:
#' \itemize{
#'   \item binary / multinomial: \code{<attr>_ei_index}, \code{<attr>_blau_index},
#'         \code{<attr>_alter_mode}
#'   \item continuous: \code{<attr>_diff_out}, \code{<attr>_diff_in},
#'         \code{<attr>_alter_mean}, \code{<attr>_alter_min}, \code{<attr>_alter_max}
#' }
#'
#' @param net An igraph object, a `network`-class object, or an adjacency matrix.
#' @param attributes A data.frame/tibble of node attributes. Either aligned by
#'   row position to the vertices of `net`, or matched via `id_col`.
#' @param attr_types Optional named character vector overriding auto-detected
#'   attribute types, e.g. \code{c(age = "continuous", dept = "multinomial")}.
#'   Allowed values: "continuous", "binary", "multinomial".
#' @param id_col Optional name of a column in `attributes` holding node
#'   identifiers matching `igraph::V(net)$name`, used to align rows.
#'
#' @return A data.frame with one row per node: a `node_id` column followed by
#'   the structural measures and the per-attribute homophily columns.
#' @export
#'
#' @examples
#' set.seed(1)
#' g <- igraph::sample_gnp(40, p = 0.08, directed = TRUE)
#' attrs <- data.frame(
#'   age = round(rnorm(40, 35, 8)),
#'   gender = sample(c("F", "M"), 40, replace = TRUE),
#'   dept = sample(c("sales", "eng", "hr"), 40, replace = TRUE)
#' )
#' head(net_measures_core(g, attrs))
net_measures_core <- function(net, attributes, attr_types = NULL, id_col = NULL) {
  g <- .as_igraph(net)
  attributes <- .align_attributes(g, attributes, id_col = id_col)
  directed <- igraph::is_directed(g)
  n <- igraph::vcount(g)

  recip <- .dyad_reciprocity(g)
  adj_w <- as.matrix(igraph::as_adjacency_matrix(
    g, sparse = FALSE,
    attr = if ("weight" %in% igraph::edge_attr_names(g)) "weight" else NULL
  ))
  weighted <- .is_weighted_mat(adj_w)

  structural <- data.frame(
    node_id            = if (!is.null(igraph::V(g)$name)) igraph::V(g)$name else seq_len(n),
    outdegree          = igraph::degree(g, mode = "out"),
    indegree           = igraph::degree(g, mode = "in"),
    reciprocity_ratio  = recip$reciprocity_ratio,
    bonacich_power     = .bonacich_power(g),
    betweenness        = igraph::betweenness(g, directed = directed, normalized = TRUE),
    isolate            = igraph::degree(g, mode = "all") == 0,
    constraint         = as.numeric(igraph::constraint(g)),
    harmonic_closeness = igraph::harmonic_centrality(g, mode = "all", normalized = TRUE),
    local_clustering   = if (weighted) .weighted_local_clustering(g) else
      igraph::transitivity(g, type = "local", isolates = "zero"),
    stringsAsFactors = FALSE
  )

  homophily <- .homophily_block(g, attributes, attr_types)

  cbind(structural, homophily)
}
