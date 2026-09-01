# Selectable node-level network measures.
#
# Every structural measure lives as its own small function in `.nm_registry`;
# `net_measures()` evaluates an arbitrary subset of them. Expensive
# intermediates that several measures share (dyad reciprocity, structural
# holes, adjacency matrices) sit in a lazily-evaluated context environment,
# so each is computed at most once per call - and not at all when no
# requested measure touches it.

#' Lazily-evaluated shared context for measure functions
#'
#' `delayedAssign()` promises: each intermediate is computed on first access
#' by a measure function and cached in the environment afterwards.
#' @noRd
.nm_context <- function(g) {
  ctx <- new.env(parent = emptyenv())
  ctx$g <- g
  ctx$n <- igraph::vcount(g)
  ctx$directed <- igraph::is_directed(g)
  delayedAssign("recip", .dyad_reciprocity(ctx$g), assign.env = ctx)
  delayedAssign("holes", .structural_holes(ctx$g), assign.env = ctx)
  delayedAssign("adj", as.matrix(igraph::as_adjacency_matrix(ctx$g, sparse = FALSE)),
                assign.env = ctx)
  delayedAssign("adj_w", as.matrix(igraph::as_adjacency_matrix(
    ctx$g, sparse = FALSE,
    attr = if ("weight" %in% igraph::edge_attr_names(ctx$g)) "weight" else NULL
  )), assign.env = ctx)
  delayedAssign("weighted", .is_weighted_mat(ctx$adj_w), assign.env = ctx)
  delayedAssign("gmode", if (ctx$directed) "digraph" else "graph", assign.env = ctx)
  ctx
}

#' @noRd
.nm_try_sna <- function(label, expr, n) {
  tryCatch(as.numeric(expr), error = function(e) {
    warning("netimpute: sna::", label, " failed (", conditionMessage(e),
            "); returning NA for this measure. This may indicate a signature ",
            "mismatch with your installed sna version - check ?sna::", label, ".",
            call. = FALSE)
    rep(NA_real_, n)
  })
}

#' Registry of node-level structural measures
#'
#' Named list in the canonical ("full"-battery) column order. Each entry maps
#' a measure name to a `function(ctx)` returning one vector of length
#' `ctx$n`. To add a measure: add an entry here, and list it in the
#' \code{\link{net_measures}} documentation (plus \code{.nm_core_set} if it
#' should be part of "core").
#' @noRd
.nm_registry <- list(
  total_degree         = function(ctx) igraph::degree(ctx$g, mode = "all"),
  indegree             = function(ctx) igraph::degree(ctx$g, mode = "in"),
  outdegree            = function(ctx) igraph::degree(ctx$g, mode = "out"),
  reciprocal_degree    = function(ctx) ctx$recip$reciprocal_degree,
  nonreciprocal_degree = function(ctx) ctx$recip$nonreciprocal_degree,
  reciprocity_ratio    = function(ctx) ctx$recip$reciprocity_ratio,
  weighted_degree      = function(ctx) igraph::strength(ctx$g, mode = "all"),
  eigenvector          = function(ctx) igraph::eigen_centrality(
    ctx$g, directed = ctx$directed)$vector,
  bonacich_power       = function(ctx) .bonacich_power(ctx$g),
  alpha_centrality     = function(ctx) {
    lambda1 <- tryCatch(max(Re(eigen(ctx$adj, only.values = TRUE)$values)),
                        error = function(e) NA_real_)
    alpha_exp <- if (is.finite(lambda1) && lambda1 > 0) 0.75 / lambda1 else 0.01
    tryCatch(as.numeric(igraph::alpha_centrality(ctx$g, alpha = alpha_exp)),
             error = function(e) rep(NA_real_, ctx$n))
  },
  pagerank             = function(ctx) igraph::page_rank(
    ctx$g, directed = ctx$directed)$vector,
  betweenness          = function(ctx) igraph::betweenness(
    ctx$g, directed = ctx$directed, normalized = TRUE),
  closeness            = function(ctx) tryCatch(
    as.numeric(igraph::closeness(ctx$g, mode = "all", normalized = TRUE)),
    error = function(e) rep(NA_real_, ctx$n)
  ),
  harmonic_closeness   = function(ctx) igraph::harmonic_centrality(
    ctx$g, mode = "all", normalized = TRUE),
  constraint           = function(ctx) as.numeric(igraph::constraint(ctx$g)),
  effective_size       = function(ctx) ctx$holes$effective_size,
  efficiency           = function(ctx) ctx$holes$efficiency,
  redundancy           = function(ctx) ctx$holes$redundancy,
  neighbor_degree      = function(ctx) .neighbor_degree(ctx$g),
  local_clustering     = function(ctx) if (ctx$weighted) {
    .weighted_local_clustering(ctx$g)
  } else {
    igraph::transitivity(ctx$g, type = "local", isolates = "zero")
  },
  coreness             = function(ctx) igraph::coreness(ctx$g, mode = "all"),
  isolate              = function(ctx) igraph::degree(ctx$g, mode = "all") == 0,
  gilschmidt           = function(ctx) .nm_try_sna(
    "gilschmidt", sna::gilschmidt(ctx$adj, g = 1), ctx$n),
  flow_betweenness     = function(ctx) .nm_try_sna(
    "flowbet", sna::flowbet(ctx$adj, gmode = ctx$gmode, rescale = TRUE), ctx$n),
  load_centrality      = function(ctx) .nm_try_sna(
    "loadcent", sna::loadcent(ctx$adj, gmode = ctx$gmode, rescale = TRUE), ctx$n),
  stress_centrality    = function(ctx) .nm_try_sna(
    "stresscent", sna::stresscent(ctx$adj, gmode = ctx$gmode, rescale = TRUE), ctx$n),
  information_centrality = function(ctx) .nm_try_sna(
    "infocent", sna::infocent(ctx$adj, gmode = ctx$gmode), ctx$n),
  prestige             = function(ctx) .nm_try_sna(
    "prestige", sna::prestige(ctx$adj, cmode = "indegree.rownorm", gmode = ctx$gmode), ctx$n)
)

#' Measures that require the sna package
#' @noRd
.nm_sna_measures <- c("gilschmidt", "flow_betweenness", "load_centrality",
                      "stress_centrality", "information_centrality", "prestige")

#' What "core" and "full" expand to ("homophily" is the per-attribute
#' homophily/alter block; it is a member of both named sets)
#' @noRd
.nm_core_set <- c("outdegree", "indegree", "reciprocity_ratio", "bonacich_power",
                  "betweenness", "isolate", "constraint", "harmonic_closeness",
                  "local_clustering", "homophily")

#' @noRd
.nm_full_set <- c(names(.nm_registry), "homophily")

#' Expand and validate a `measure_set` specification
#'
#' Accepts any mix of `"core"`, `"full"`, `"homophily"`, and individual
#' measure names; named sets are expanded and the union is returned, keeping
#' first-appearance order (so `"core"` alone reproduces the historical core
#' column order). Idempotent: the returned vector is itself a valid
#' `measure_set`.
#' @noRd
.resolve_measure_set <- function(measure_set) {
  if (!is.character(measure_set) || length(measure_set) == 0 || anyNA(measure_set)) {
    stop("`measure_set` must be a non-empty character vector.", call. = FALSE)
  }
  valid <- c("core", "full", "homophily", names(.nm_registry))
  bad <- setdiff(measure_set, valid)
  if (length(bad)) {
    stop("Unknown entries in `measure_set`: ", toString(bad),
         ". Valid entries are 'core', 'full', 'homophily', and the individual ",
         "measures: ", toString(names(.nm_registry)), ".",
         call. = FALSE)
  }
  expanded <- unlist(lapply(measure_set, function(s) {
    switch(s, core = .nm_core_set, full = .nm_full_set, s)
  }), use.names = FALSE)
  unique(expanded)
}

#' Node-level network measures with a selectable measure set
#'
#' Computes any subset of the node-level structural measures implemented in
#' netimpute, plus (optionally) the per-attribute homophily/alter block.
#' \code{\link{net_measures_core}} and \code{\link{net_measures_full}} are
#' thin wrappers around this function with `measure_set = "core"` and
#' `measure_set = "full"`.
#'
#' `measure_set` may combine named sets and individual measures freely; the
#' result is their union. For example
#' `measure_set = c("core", "total_degree", "pagerank")` returns every core
#' measure plus `total_degree` and `pagerank`.
#'
#' Available entries:
#' \describe{
#'   \item{Named sets}{\code{"core"} (compact, low-collinearity set:
#'     outdegree, indegree, reciprocity_ratio, bonacich_power, betweenness,
#'     isolate, constraint, harmonic_closeness, local_clustering, homophily),
#'     \code{"full"} (every measure below plus homophily).}
#'   \item{Degree family}{\code{total_degree}, \code{indegree},
#'     \code{outdegree}, \code{reciprocal_degree}, \code{nonreciprocal_degree},
#'     \code{reciprocity_ratio}, \code{weighted_degree}}
#'   \item{Influence}{\code{eigenvector}, \code{bonacich_power},
#'     \code{alpha_centrality}, \code{pagerank}, \code{gilschmidt}
#'     (requires \pkg{sna})}
#'   \item{Brokerage/position}{\code{betweenness}, \code{load_centrality},
#'     \code{stress_centrality}, \code{flow_betweenness} (the latter three
#'     require \pkg{sna}), \code{constraint}, \code{effective_size},
#'     \code{efficiency}, \code{redundancy}}
#'   \item{Distance-based}{\code{closeness}, \code{harmonic_closeness},
#'     \code{information_centrality} (requires \pkg{sna})}
#'   \item{Prestige}{\code{prestige} (requires \pkg{sna})}
#'   \item{Position/closure}{\code{neighbor_degree}, \code{local_clustering},
#'     \code{coreness}, \code{isolate}}
#'   \item{Attribute-based}{\code{"homophily"} - the per-attribute
#'     homophily/alter-similarity block (E-I index, Blau index, alter
#'     mode/means/min/max, absolute differences; which columns appear depends
#'     on each attribute's type). Included in both \code{"core"} and
#'     \code{"full"}; when an explicit list of structural measures is given
#'     without \code{"homophily"}, no attribute-based columns are returned
#'     and `attributes` may be omitted.}
#' }
#'
#' Requested measures that need \pkg{sna} are filled with `NA` (with a single
#' message) if \pkg{sna} is not installed or `use_sna = FALSE`.
#'
#' Shared intermediate quantities (e.g. the dyad-reciprocity tabulation used
#' by \code{reciprocal_degree}, \code{nonreciprocal_degree} and
#' \code{reciprocity_ratio}) are computed lazily and at most once per call,
#' so the cost of a call scales with the measures actually requested.
#'
#' Only binary (0/1) and non-negative weighted networks are currently
#' supported (a network with negative tie values raises an error). For
#' non-negative weighted networks, \code{reciprocity_ratio} and
#' \code{local_clustering} use weighted redefinitions rather than the binary
#' formulas - see \code{\link{net_measures_core}} for details.
#'
#' @param net An igraph object, a `network`-class object, or an adjacency
#'   matrix.
#' @param attributes A data.frame/tibble of node attributes, used by the
#'   `"homophily"` block. Either aligned by row position to the vertices of
#'   `net`, or matched via `id_col`. May be `NULL` when `"homophily"` is not
#'   among the requested measures.
#' @param measure_set Character vector: any combination of `"core"`,
#'   `"full"`, `"homophily"`, and individual measure names (see Details).
#'   Default `"core"`.
#' @param attr_types Optional named character vector overriding auto-detected
#'   attribute types, e.g. \code{c(age = "continuous", dept = "multinomial")}.
#'   Allowed values: "continuous", "binary", "multinomial".
#' @param id_col Optional name of a column in `attributes` holding node
#'   identifiers matching `igraph::V(net)$name`, used to align rows.
#' @param use_sna Logical; if `TRUE` (default) and \pkg{sna} is installed,
#'   compute any requested sna-only measures (gilschmidt, flow/load/stress
#'   centrality, information centrality, prestige); otherwise those columns
#'   are `NA`.
#'
#' @return A data.frame with one row per node: `node_id`, the requested
#'   structural measures (in a canonical column order), and - if
#'   `"homophily"` was requested - the per-attribute homophily columns.
#' @seealso \code{\link{net_measures_core}} and \code{\link{net_measures_full}}
#'   for the two ready-made measure sets; \code{\link{net_predictors}} to
#'   turn measures from several networks into one predictor frame;
#'   \code{\link{netmice}}, which computes these as attribute predictors.
#' @export
#'
#' @examples
#' set.seed(1)
#' g <- igraph::sample_gnp(40, p = 0.08, directed = TRUE)
#' attrs <- data.frame(
#'   age = round(rnorm(40, 35, 8)),
#'   gender = sample(c("F", "M"), 40, replace = TRUE)
#' )
#' # core battery plus two extra measures
#' head(net_measures(g, attrs, measure_set = c("core", "total_degree", "pagerank")))
#' # a bespoke structural-only selection (no attributes needed)
#' head(net_measures(g, measure_set = c("indegree", "betweenness")))
net_measures <- function(net, attributes = NULL, measure_set = "core",
                         attr_types = NULL, id_col = NULL, use_sna = TRUE) {
  measures <- .resolve_measure_set(measure_set)
  g <- .as_igraph(net)

  want_homophily <- "homophily" %in% measures
  if (want_homophily) {
    if (is.null(attributes)) {
      stop("`attributes` must be supplied when the 'homophily' block is part ",
           "of `measure_set` (it is included in both 'core' and 'full').",
           call. = FALSE)
    }
    attributes <- .align_attributes(g, attributes, id_col = id_col)
  }

  ctx <- .nm_context(g)
  structural_names <- setdiff(measures, "homophily")

  sna_requested <- intersect(structural_names, .nm_sna_measures)
  sna_ok <- use_sna && requireNamespace("sna", quietly = TRUE)
  if (length(sna_requested) && use_sna && !sna_ok) {
    message("netimpute: package 'sna' not installed - ",
            toString(sna_requested),
            " will be NA. Install with install.packages('sna') to compute them.")
  }

  cols <- lapply(structural_names, function(nm) {
    if (nm %in% .nm_sna_measures && !sna_ok) return(rep(NA_real_, ctx$n))
    .nm_registry[[nm]](ctx)
  })
  names(cols) <- structural_names

  node_id <- if (!is.null(igraph::V(g)$name)) igraph::V(g)$name else seq_len(ctx$n)
  out <- do.call(data.frame,
                 c(list(node_id = node_id), cols, list(stringsAsFactors = FALSE)))

  if (want_homophily) {
    out <- cbind(out, .homophily_block(g, attributes, attr_types))
  }
  out
}
