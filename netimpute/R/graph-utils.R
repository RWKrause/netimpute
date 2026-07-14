# Internal helpers: graph coercion, attribute typing, dyad reciprocity,
# Bonacich power, and structural-holes measures. Nothing here is exported.

#' Reject negative tie values
#'
#' netimpute currently supports binary (0/1) and non-negative weighted
#' networks only. Signed networks (negative tie values) are not yet
#' supported - several measures (e.g. path-based centrality, which requires
#' non-negative edge costs; the weighted transitivity/reciprocity formulas
#' below, which assume weights are magnitudes) are not well-defined or not
#' implemented for them.
#' @keywords internal
.reject_negative <- function(vals) {
  if (any(vals < 0, na.rm = TRUE)) {
    stop("netimpute currently only supports binary (0/1) or non-negative weighted ",
         "networks. Found negative tie value(s) - signed network support is not ",
         "implemented yet. Try splitting the network into positive and negative ",
         "ties and imputing these separately. (Signed nature here not guaranteed.)",
         call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
.as_igraph <- function(net) {
  if (inherits(net, "igraph")) {
    if ("weight" %in% igraph::edge_attr_names(net)) .reject_negative(igraph::E(net)$weight)
    return(net)
  }

  if (inherits(net, "network")) {
    if (!requireNamespace("network", quietly = TRUE)) {
      stop("Package 'network' is needed to convert a `network`-class object.", call. = FALSE)
    }
    m <- as.matrix(net, matrix.type = "adjacency")
    .reject_negative(m)
    directed <- network::is.directed(net)
    weighted <- any(m[m != 0] != 1)
    return(igraph::graph_from_adjacency_matrix(
      m, mode = if (directed) "directed" else "undirected",
      weighted = if (weighted) TRUE else NULL, diag = FALSE
    ))
  }

  if (is.matrix(net)) {
    .reject_negative(net)
    directed <- !isSymmetric(unname(net))
    weighted <- any(net[net != 0] != 1)
    return(igraph::graph_from_adjacency_matrix(
      net, mode = if (directed) "directed" else "undirected",
      weighted = if (weighted) TRUE else NULL, diag = FALSE
    ))
  }

  stop("`net` must be an igraph object, a network-class object, or an adjacency matrix.",
       call. = FALSE)
}

#' @keywords internal
.is_weighted_mat <- function(m) {
  any(m[m != 0] != 1)
}

#' @keywords internal
.align_attributes <- function(g, attributes, id_col = NULL) {
  n <- igraph::vcount(g)
  attributes <- as.data.frame(attributes)

  if (!is.null(id_col)) {
    vnames <- igraph::V(g)$name
    if (is.null(vnames)) {
      stop("`id_col` was supplied but the network has no vertex names to match against.",
           call. = FALSE)
    }
    if (!id_col %in% names(attributes)) {
      stop(sprintf("id_col '%s' not found in `attributes`.", id_col), call. = FALSE)
    }
    ord <- match(vnames, as.character(attributes[[id_col]]))
    if (anyNA(ord)) {
      stop("Not every vertex in the network could be matched to a row in `attributes` via id_col.",
           call. = FALSE)
    }
    # Once rows are aligned the identifier has served its purpose; leaving it
    # in would have every downstream consumer treat it as an n-level
    # multinomial attribute (n-1 ego + n-1 alter dummies in the dyad design,
    # meaningless homophily measures in net_measures_*).
    attributes <- attributes[ord, setdiff(names(attributes), id_col),
                             drop = FALSE]
  } else if (nrow(attributes) != n) {
    stop(sprintf(
      "`attributes` has %d rows but the network has %d nodes. Supply `id_col` to match by name, or align rows by position.",
      nrow(attributes), n
    ), call. = FALSE)
  }

  rownames(attributes) <- NULL
  attributes
}

#' Infer whether an attribute is continuous, binary, or multinomial
#' @keywords internal
.infer_attr_type <- function(x) {
  u <- unique(x[!is.na(x)])
  if (length(u) <= 2) return("binary")
  if (is.numeric(x) && !is.factor(x)) return("continuous")
  return("multinomial")
}

#' @keywords internal
.resolve_attr_types <- function(attributes, attr_types) {
  nms <- names(attributes)
  resolved <- stats::setNames(vector("character", length(nms)), nms)
  for (nm in nms) {
    if (!is.null(attr_types) && nm %in% names(attr_types)) {
      resolved[nm] <- attr_types[[nm]]
    } else {
      resolved[nm] <- .infer_attr_type(attributes[[nm]])
    }
  }
  bad <- setdiff(unique(resolved), c("continuous", "binary", "multinomial"))
  if (length(bad)) {
    stop("attr_types must be one of 'continuous', 'binary', 'multinomial'. Got: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  inferred <- setdiff(nms, names(attr_types))
  if (length(inferred)) {
    message("netimpute: inferred attribute types -> ",
            paste(sprintf("%s = %s", inferred, resolved[inferred]), collapse = ", "))
  }
  resolved
}

#' Node-level dyadic reciprocity
#'
#' Binary networks: the standard mutual/either dyad-count ratio (in \[0, 1\],
#' higher = more reciprocal). Weighted networks: the average \code{|w_ij -
#' w_ji|} over alters with a tie in either direction - this is *not* a \[0,
#' 1\] ratio and is in the original tie-weight units, with the direction of
#' interpretation flipped (lower = more reciprocal, 0 = perfectly symmetric
#' weighted ties). \code{reciprocal_degree}/\code{nonreciprocal_degree} are
#' binary-only concepts (discrete dyad counts) and are `NA` for weighted
#' networks.
#' @keywords internal
.dyad_reciprocity <- function(g) {
  n <- igraph::vcount(g)
  if (!igraph::is_directed(g)) {
    return(data.frame(
      reciprocal_degree = rep(NA_real_, n),
      nonreciprocal_degree = rep(NA_real_, n),
      reciprocity_ratio = rep(NA_real_, n)
    ))
  }
  adj <- as.matrix(igraph::as_adjacency_matrix(
    g, sparse = FALSE,
    attr = if ("weight" %in% igraph::edge_attr_names(g)) "weight" else NULL
  ))
  diag(adj) <- 0

  if (.is_weighted_mat(adj)) {
    either <- (adj != 0) | (t(adj) != 0)
    either_n <- rowSums(either)
    absdiff_sum <- rowSums(abs(adj - t(adj)) * either)
    return(data.frame(
      reciprocal_degree = rep(NA_real_, n),
      nonreciprocal_degree = rep(NA_real_, n),
      reciprocity_ratio = ifelse(either_n == 0, NA_real_,
                                 absdiff_sum / either_n)
    ))
  }

  bin <- adj
  bin[bin != 0] <- 1
  mutual_mat <- bin * t(bin)
  either_mat <- ((bin + t(bin)) > 0) * 1
  mutual_n <- rowSums(mutual_mat)
  either_n <- rowSums(either_mat)
  data.frame(
    reciprocal_degree = mutual_n,
    nonreciprocal_degree = either_n - mutual_n,
    reciprocity_ratio = ifelse(either_n == 0, NA_real_, mutual_n / either_n)
  )
}

#' Weighted local clustering coefficient (Opsahl & Panzarasa 2009,
#' geometric-mean formulation)
#'
#' For node i and each pair of its neighbours (j, k): the "two-path" value is
#' the geometric mean \code{sqrt(w_ij * w_ik)}. The denominator is the sum of
#' this value over *every* pair of neighbours (i.e. every two-path centred on
#' i); the numerator is the same sum restricted to pairs where the triangle
#' is actually closed (j and k are themselves tied). The ratio reduces to the
#' ordinary (binary) local clustering coefficient when all weights are 1.
#' Direction is collapsed first (mean of the two directed weights) since
#' clustering is an undirected-triad concept here.
#' @keywords internal
.weighted_local_clustering <- function(g) {
  n <- igraph::vcount(g)
  gu <- if (igraph::is_directed(g)) {
    igraph::as_undirected(g, mode = "collapse", edge.attr.comb = list(weight = "mean"))
  } else {
    g
  }
  w <- as.matrix(igraph::as_adjacency_matrix(gu, sparse = FALSE, attr = "weight"))
  diag(w) <- 0
  a <- (w != 0) * 1

  vapply(seq_len(n), function(i) {
    nb <- which(a[i, ] == 1)
    k <- length(nb)
    if (k < 2) return(0)
    denom <- 0
    numer <- 0
    for (p in seq_len(k - 1)) {
      for (q in (p + 1):k) {
        j <- nb[p]; h <- nb[q]
        gm <- sqrt(w[i, j] * w[i, h])
        denom <- denom + gm
        if (a[j, h] == 1) numer <- numer + gm
      }
    }
    if (denom == 0) 0 else numer / denom
  }, numeric(1))
}

#' Bonacich power centrality with an automatically convergent exponent
#' @keywords internal
.bonacich_power <- function(g, exponent = NULL) {
  n <- igraph::vcount(g)
  adj <- as.matrix(igraph::as_adjacency_matrix(g, sparse = FALSE))
  if (is.null(exponent)) {
    lambda1 <- tryCatch(max(Re(eigen(adj, only.values = TRUE)$values)),
                         error = function(e) NA_real_)
    exponent <- if (is.finite(lambda1) && lambda1 > 0) 0.75 / lambda1 else 0.01
  }
  out <- tryCatch(
    as.numeric(igraph::power_centrality(g, exponent = exponent, rescale = FALSE)),
    error = function(e) {
      warning("Bonacich power centrality failed (", e$message, "); returning NA.", call. = FALSE)
      rep(NA_real_, n)
    }
  )
  out
}

#' Burt's structural holes: effective size, redundancy, efficiency
#' @keywords internal
.structural_holes <- function(g) {
  gu <- igraph::as_undirected(g, mode = "collapse")
  n <- igraph::vcount(g)
  eff_size <- redundancy <- efficiency <- numeric(n)
  for (i in seq_len(n)) {
    nb <- as.integer(igraph::neighbors(gu, i))
    d <- length(nb)
    if (d == 0) {
      eff_size[i] <- 0; redundancy[i] <- NA_real_; efficiency[i] <- NA_real_
      next
    }
    sub <- igraph::induced_subgraph(gu, nb)
    t_i <- igraph::ecount(sub)
    redundancy[i] <- 2 * t_i / d
    eff_size[i] <- d - redundancy[i]
    efficiency[i] <- eff_size[i] / d
  }
  data.frame(effective_size = eff_size, redundancy = redundancy, efficiency = efficiency)
}

#' Average nearest-neighbour (alter) degree
#' @keywords internal
.neighbor_degree <- function(g) {
  deg <- igraph::degree(g, mode = "all")
  n <- igraph::vcount(g)
  vapply(seq_len(n), function(i) {
    nb <- as.integer(igraph::neighbors(g, i, mode = "all"))
    if (length(nb) == 0) return(NA_real_)
    mean(deg[nb])
  }, numeric(1))
}
