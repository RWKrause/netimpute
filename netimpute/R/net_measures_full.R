#' Full node-level network measure battery
#'
#' Computes the full battery of node-level structural measures (the
#' superset of \code{\link{net_measures_core}}) plus, for every column of
#' `attributes`, the same attribute-scale-appropriate homophily/alter
#' measures used by \code{net_measures_core}.
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
#' g <- igraph::sample_gnp(40, p = 0.08, directed = TRUE)
#' attrs <- data.frame(
#'   age = round(rnorm(40, 35, 8)),
#'   gender = sample(c("F", "M"), 40, replace = TRUE)
#' )
#' head(net_measures_full(g, attrs))
net_measures_full <- function(net, attributes, attr_types = NULL, id_col = NULL, use_sna = TRUE) {
  g <- .as_igraph(net)
  attributes <- .align_attributes(g, attributes, id_col = id_col)
  directed <- igraph::is_directed(g)
  n <- igraph::vcount(g)

  recip <- .dyad_reciprocity(g)
  holes <- .structural_holes(g)
  adj_w <- as.matrix(igraph::as_adjacency_matrix(
    g, sparse = FALSE,
    attr = if ("weight" %in% igraph::edge_attr_names(g)) "weight" else NULL
  ))
  weighted <- .is_weighted_mat(adj_w)

  alpha_exp <- {
    adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
    lambda1 <- tryCatch(max(Re(eigen(adj, only.values = TRUE)$values)), error = function(e) NA_real_)
    if (is.finite(lambda1) && lambda1 > 0) 0.75 / lambda1 else 0.01
  }

  structural <- data.frame(
    node_id               = if (!is.null(igraph::V(g)$name)) igraph::V(g)$name else seq_len(n),
    total_degree          = igraph::degree(g, mode = "all"),
    indegree              = igraph::degree(g, mode = "in"),
    outdegree             = igraph::degree(g, mode = "out"),
    reciprocal_degree     = recip$reciprocal_degree,
    nonreciprocal_degree  = recip$nonreciprocal_degree,
    reciprocity_ratio     = recip$reciprocity_ratio,
    weighted_degree       = igraph::strength(g, mode = "all"),
    eigenvector           = igraph::eigen_centrality(g, directed = directed)$vector,
    bonacich_power        = .bonacich_power(g),
    alpha_centrality      = tryCatch(
      as.numeric(igraph::alpha_centrality(g, alpha = alpha_exp)),
      error = function(e) rep(NA_real_, n)
    ),
    pagerank              = igraph::page_rank(g, directed = directed)$vector,
    betweenness           = igraph::betweenness(g, directed = directed, normalized = TRUE),
    closeness             = tryCatch(
      as.numeric(igraph::closeness(g, mode = "all", normalized = TRUE)),
      error = function(e) rep(NA_real_, n)
    ),
    harmonic_closeness    = igraph::harmonic_centrality(g, mode = "all", normalized = TRUE),
    constraint            = as.numeric(igraph::constraint(g)),
    effective_size        = holes$effective_size,
    efficiency            = holes$efficiency,
    redundancy            = holes$redundancy,
    neighbor_degree       = .neighbor_degree(g),
    local_clustering      = if (weighted) .weighted_local_clustering(g) else
      igraph::transitivity(g, type = "local", isolates = "zero"),
    coreness               = igraph::coreness(g, mode = "all"),
    isolate               = igraph::degree(g, mode = "all") == 0,
    stringsAsFactors = FALSE
  )

  sna_ok <- use_sna && requireNamespace("sna", quietly = TRUE)
  if (use_sna && !sna_ok) {
    message("netimpute: package 'sna' not installed - gilschmidt, flow/load/stress ",
            "centrality, information centrality and prestige will be NA. ",
            "Install with install.packages('sna') to compute them.")
  }

  sna_measures <- data.frame(
    gilschmidt            = rep(NA_real_, n),
    flow_betweenness      = rep(NA_real_, n),
    load_centrality       = rep(NA_real_, n),
    stress_centrality     = rep(NA_real_, n),
    information_centrality = rep(NA_real_, n),
    prestige              = rep(NA_real_, n)
  )

  if (sna_ok) {
    m <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
    gmode <- if (directed) "digraph" else "graph"
    try_sna <- function(label, expr) {
      tryCatch(as.numeric(expr), error = function(e) {
        warning("netimpute: sna::", label, " failed (", conditionMessage(e),
                "); returning NA for this measure. This may indicate a signature ",
                "mismatch with your installed sna version - check ?sna::", label, ".",
                call. = FALSE)
        rep(NA_real_, n)
      })
    }
    sna_measures$gilschmidt        <- try_sna("gilschmidt", sna::gilschmidt(m, g = 1))
    sna_measures$flow_betweenness  <- try_sna("flowbet", sna::flowbet(m, gmode = gmode, rescale = TRUE))
    sna_measures$load_centrality   <- try_sna("loadcent", sna::loadcent(m, gmode = gmode, rescale = TRUE))
    sna_measures$stress_centrality <- try_sna("stresscent", sna::stresscent(m, gmode = gmode, rescale = TRUE))
    sna_measures$information_centrality <- try_sna("infocent", sna::infocent(m, gmode = gmode))
    sna_measures$prestige          <- try_sna("prestige", sna::prestige(m, cmode = "indegree.rownorm", gmode = gmode))
  }

  homophily <- .homophily_block(g, attributes, attr_types)

  cbind(structural, sna_measures, homophily)
}
