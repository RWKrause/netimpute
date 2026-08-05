#' Full node-level network measure battery
#'
#' Computes the full battery of node-level structural measures (the
#' superset of \code{\link{net_measures_core}}) plus, for every column of
#' `attributes`, the same attribute-scale-appropriate homophily/alter
#' measures used by \code{net_measures_core}. Equivalent to
#' \code{\link{net_measures}} with `measure_set = "full"`; use
#' \code{net_measures()} directly to select individual measures.
#'
#' Structural measures returned (all self-contained via igraph unless noted):
#' \itemize{
#'   \item Degree family: \code{total_degree}, \code{indegree}, \code{outdegree},
#'         \code{reciprocal_degree}, \code{nonreciprocal_degree},
#'         \code{reciprocity_ratio}, \code{weighted_degree}
#'   \item Influence: \code{eigenvector}, \code{bonacich_power}, \code{alpha_centrality},
#'         \code{pagerank}, \code{gilschmidt} (requires \pkg{sna})
#'   \item Brokerage/position: \code{betweenness}, \code{load_centrality},
#'         \code{stress_centrality}, \code{flow_betweenness} (latter three require \pkg{sna}),
#'         \code{constraint}, \code{effective_size}, \code{efficiency}, \code{redundancy}
#'   \item Distance-based: \code{closeness}, \code{harmonic_closeness},
#'         \code{information_centrality} (requires \pkg{sna})
#'   \item Prestige: \code{prestige} (requires \pkg{sna})
#'   \item Position/closure: \code{neighbor_degree}, \code{local_clustering},
#'         \code{coreness}, \code{isolate}
#' }
#'
#' Measures that require \pkg{sna} are filled with `NA` and reported in a
#' single message if \pkg{sna} is not installed.
#'
#' Only binary (0/1) and non-negative weighted networks are currently
#' supported (a network with negative tie values raises an error). For
#' non-negative weighted networks, \code{reciprocity_ratio} and
#' \code{local_clustering} use weighted redefinitions rather than the binary
#' formulas - see \code{\link{net_measures_core}} for details.
#'
#' @inheritParams net_measures_core
#' @param use_sna Logical; if `TRUE` (default) and \pkg{sna} is installed, also
#'   compute the sna-only measures (gilschmidt, flow/load/stress centrality,
#'   information centrality, prestige).
#'
#' @return A data.frame with one row per node, `node_id` followed by all
#'   structural measures and the per-attribute homophily columns.
#' @export
#'
#' @examples
#' set.seed(1)
#' g <- igraph::sample_gnp(25, p = 0.1, directed = TRUE)
#' attrs <- data.frame(
#'   age = round(rnorm(25, 35, 8)),
#'   gender = sample(c("F", "M"), 25, replace = TRUE)
#' )
#' head(net_measures_full(g, attrs))
net_measures_full <- function(net, attributes, attr_types = NULL, id_col = NULL, use_sna = TRUE) {
  net_measures(net, attributes, measure_set = "full",
               attr_types = attr_types, id_col = id_col, use_sna = use_sna)
}
