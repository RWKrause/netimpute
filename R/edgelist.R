# Edgelist input for netmice(). Everything here runs BEFORE the sampler and
# produces the same named list of aligned adjacency matrices that the matrix
# input path produces, so no downstream code needs to know which was used.

#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Colour text firebrick4 when the console can show it
#'
#' Used to mark the nodes an edgelist adds to `data`: those have no observed
#' attributes at all, so they carry the most imputation variance and are the
#' thing a user most needs to notice in a long warning.
#' @noRd
.firebrick <- function(txt) {
  supported <- interactive() &&
    !identical(Sys.getenv("TERM"), "dumb") &&
    identical(Sys.getenv("NO_COLOR"), "")
  if (!supported) return(txt)
  rgb <- grDevices::col2rgb("firebrick4")
  paste0("\033[38;2;", rgb[1], ";", rgb[2], ";", rgb[3], "m", txt, "\033[39m")
}

#' Resolve the sender/receiver/value columns of an edgelist
#'
#' With `edgelist_names = NULL` the first two non-split columns are taken as
#' sender and receiver and a third (if any) as the tie value, and that
#' assumption is printed. Any further non-split column is an error: guessing
#' which of several unexplained columns holds the value would silently build
#' the wrong networks.
#' @noRd
.el_resolve_names <- function(el, edgelist_names, split_cols, format,
                              printFlag = TRUE, what = "`networks`") {
  cn <- names(el)
  if (!is.null(edgelist_names)) {
    if (is.null(names(edgelist_names))) {
      stop("`edgelist_options$edgelist_names` must be a NAMED vector with ",
           "entries 'sender' and 'receiver' (optionally 'value').",
           call. = FALSE)
    }
    bad <- setdiff(names(edgelist_names), c("sender", "receiver", "value"))
    if (length(bad)) {
      stop("`edgelist_options$edgelist_names` may only name 'sender', ",
           "'receiver' and 'value'; got: ", toString(bad), ".", call. = FALSE)
    }
    if (!all(c("sender", "receiver") %in% names(edgelist_names))) {
      stop("`edgelist_options$edgelist_names` must contain both 'sender' ",
           "and 'receiver'.", call. = FALSE)
    }
    miss <- setdiff(unname(edgelist_names), cn)
    if (length(miss)) {
      stop("`edgelist_options$edgelist_names` refers to column(s) not in ",
           what, ": ", toString(miss), ".", call. = FALSE)
    }
    sender <- unname(edgelist_names[["sender"]])
    receiver <- unname(edgelist_names[["receiver"]])
    value <- if ("value" %in% names(edgelist_names))
      unname(edgelist_names[["value"]]) else NULL
  } else {
    non_split <- setdiff(cn, split_cols)
    if (length(non_split) < 2) {
      stop("An edgelist needs a sender and a receiver column; ", what,
           " has ", length(non_split), " column(s) outside `edgelist_split`.",
           call. = FALSE)
    }
    sender <- non_split[1]
    receiver <- non_split[2]
    rest <- non_split[-(1:2)]
    if (length(rest) > 1) {
      stop("Cannot tell which column holds the tie value: ", what, " has ",
           length(rest), " columns beyond sender ('", sender,
           "') and receiver ('", receiver, "') that are not in ",
           "`edgelist_split`: ", toString(rest), ". Name the columns ",
           "explicitly via `edgelist_options$edgelist_names`, add them to ",
           "`edgelist_split`, or remove them.", call. = FALSE)
    }
    value <- if (length(rest)) rest else NULL
    if (printFlag) {
      message("netimpute: reading ", what, " as sender = '", sender,
              "', receiver = '", receiver, "'",
              if (is.null(value)) " (unweighted: a listed edge is a 1)" else
                paste0(", value = '", value, "'"),
              ". Set `edgelist_options$edgelist_names` to override.")
    }
  }
  if (identical(format, "wide") && !is.null(value)) {
    stop("`edgelist_format = 'wide'` takes the tie values from the ",
         "`edgelist_split` columns themselves, so a separate 'value' column ",
         "is ambiguous. Drop it, or use `edgelist_format = 'long'`.",
         call. = FALSE)
  }
  list(sender = sender, receiver = receiver, value = value)
}

#' Network names and the edgelist rows behind them
#'
#' long: one network per observed combination of the `edgelist_split` columns,
#' named `<col>_<value>` joined by `_`. wide: one network per split column,
#' named after the column, with every row belonging to all of them.
#' @noRd
.el_network_keys <- function(el, format, split_cols) {
  if (identical(format, "wide")) {
    if (!length(split_cols)) {
      stop("`edgelist_format = 'wide'` needs `edgelist_split` to name the ",
           "value columns, one per network.", call. = FALSE)
    }
    miss <- setdiff(split_cols, names(el))
    if (length(miss)) {
      stop("`edgelist_split` names column(s) not in the edgelist: ",
           toString(miss), ".", call. = FALSE)
    }
    non_num <- split_cols[!vapply(el[split_cols], is.numeric, logical(1))]
    if (length(non_num)) {
      stop("With `edgelist_format = 'wide'` every `edgelist_split` column ",
           "must be numeric (the tie value, or NA): ", toString(non_num),
           " is not.", call. = FALSE)
    }
    return(list(
      names = split_cols,
      rows = stats::setNames(rep(list(seq_len(nrow(el))), length(split_cols)),
                             split_cols),
      components = stats::setNames(as.list(split_cols), split_cols)))
  }
  if (!length(split_cols)) {
    return(list(names = "net1",
                rows = list(net1 = seq_len(nrow(el))),
                components = list(net1 = character(0))))
  }
  miss <- setdiff(split_cols, names(el))
  if (length(miss)) {
    stop("`edgelist_split` names column(s) not in the edgelist: ",
         toString(miss), ".", call. = FALSE)
  }
  parts <- lapply(split_cols, function(cl)
    paste0(cl, "_", as.character(el[[cl]])))
  key <- do.call(paste, c(parts, sep = "_"))
  uk <- unique(key)
  comp <- lapply(uk, function(k) {
    i <- match(k, key)
    vapply(parts, function(p) p[i], character(1))
  })
  list(names = uk,
       rows = stats::setNames(lapply(uk, function(k) which(key == k)), uk),
       components = stats::setNames(comp, uk))
}

#' Resolve the `nodelist` option into a per-network node set
#'
#' Accepts a column of `data`, a plain vector, or a named list keyed by
#' network, by a single split value ('wave_1'), or by a full combination
#' ('wave_1_question_1'). When several split values address one network their
#' node sets are intersected, so a node must appear in every one of them.
#' @noRd
.el_resolve_nodelist <- function(nodelist, data, keys, base_ids) {
  nms <- keys$names
  if (is.null(nodelist)) {
    return(stats::setNames(rep(list(base_ids), length(nms)), nms))
  }
  if (is.character(nodelist) && length(nodelist) == 1 &&
      nodelist %in% names(data)) {
    v <- unique(as.character(data[[nodelist]]))
    return(stats::setNames(rep(list(v), length(nms)), nms))
  }
  if (!is.list(nodelist)) {
    v <- unique(as.character(nodelist))
    return(stats::setNames(rep(list(v), length(nms)), nms))
  }
  if (is.null(names(nodelist)) || !all(nzchar(names(nodelist)))) {
    stop("When `edgelist_options$nodelist` is a list, every element must be ",
         "named - after a network, after one split value (e.g. 'wave_1'), or ",
         "after a full combination (e.g. 'wave_1_question_1').", call. = FALSE)
  }
  out <- stats::setNames(vector("list", length(nms)), nms)
  for (nm in nms) {
    if (nm %in% names(nodelist)) {
      out[[nm]] <- unique(as.character(nodelist[[nm]]))
      next
    }
    comp <- keys$components[[nm]]
    hit <- intersect(comp, names(nodelist))
    if (!length(hit)) {
      stop("`edgelist_options$nodelist` does not cover network '", nm,
           "'. Give it an entry named '", nm, "', or entries named after its ",
           "split value(s): ", toString(comp), ". Available names: ",
           toString(names(nodelist)), ".", call. = FALSE)
    }
    sets <- lapply(nodelist[hit], function(v) unique(as.character(v)))
    out[[nm]] <- Reduce(intersect, sets)
    if (!length(out[[nm]])) {
      stop("The `nodelist` entries for network '", nm, "' (", toString(hit),
           ") have no nodes in common, so that network would be empty.",
           call. = FALSE)
    }
  }
  out
}

#' Resolve the `missing` option into per-network out/in node vectors
#'
#' A bare vector means *outgoing* ties are missing - the actor-non-response
#' case, where a person did not report but others could still nominate them.
#' @noRd
.el_resolve_missing <- function(missing, keys) {
  nms <- keys$names
  empty <- stats::setNames(
    rep(list(list(out = character(0), `in` = character(0))), length(nms)), nms)
  if (is.null(missing)) return(empty)

  as_dirs <- function(x, label) {
    if (!is.list(x)) {
      return(list(out = unique(as.character(x)), `in` = character(0)))
    }
    bad <- setdiff(names(x), c("out", "in"))
    if (length(bad)) {
      stop("`edgelist_options$missing` (", label, ") may only contain 'out' ",
           "and 'in'; got: ", toString(bad), ".", call. = FALSE)
    }
    list(out = unique(as.character(x$out %||% character(0))),
         `in` = unique(as.character(x[["in"]] %||% character(0))))
  }
  is_dir_list <- is.list(missing) && !is.null(names(missing)) &&
    length(names(missing)) && all(names(missing) %in% c("out", "in"))
  if (!is.list(missing) || is_dir_list) {
    d <- as_dirs(missing, "all networks")
    return(stats::setNames(rep(list(d), length(nms)), nms))
  }
  if (is.null(names(missing)) || !all(nzchar(names(missing)))) {
    stop("When `edgelist_options$missing` is a per-network list, every ",
         "element must be named after a network, a split value, or a full ",
         "combination.", call. = FALSE)
  }
  out <- empty
  for (nm in nms) {
    hit <- if (nm %in% names(missing)) nm else
      intersect(keys$components[[nm]], names(missing))
    if (!length(hit)) next
    ds <- lapply(hit, function(h) as_dirs(missing[[h]], h))
    out[[nm]] <- list(
      out = unique(unlist(lapply(ds, `[[`, "out"), use.names = FALSE)),
      `in` = unique(unlist(lapply(ds, `[[`, "in"), use.names = FALSE)))
  }
  out
}

#' Resolve the `directed` option into a named logical vector (NA = infer)
#' @noRd
.el_resolve_directed <- function(directed, keys) {
  nms <- keys$names
  if (is.null(directed)) return(NULL)
  word_to_flag <- function(w, label) {
    if (is.logical(w) && length(w) == 1 && !is.na(w)) return(w)
    w <- tolower(as.character(w))
    if (length(w) != 1) {
      stop("`edgelist_options$directed` (", label, ") must be a single word.",
           call. = FALSE)
    }
    if (w %in% c("directed", "digraph")) return(TRUE)
    if (w %in% c("undirected", "graph")) return(FALSE)
    stop("`edgelist_options$directed` (", label, ") must be one of ",
         "'directed', 'undirected', 'digraph', 'graph'; got '", w, "'.",
         call. = FALSE)
  }
  if (length(directed) == 1 && is.null(names(directed))) {
    return(stats::setNames(rep(word_to_flag(directed, "all networks"),
                               length(nms)), nms))
  }
  if (is.null(names(directed))) {
    stop("`edgelist_options$directed` must be a single word applying to all ",
         "networks, or be named per network / split value.", call. = FALSE)
  }
  out <- stats::setNames(rep(NA, length(nms)), nms)
  for (nm in nms) {
    hit <- if (nm %in% names(directed)) nm else
      intersect(keys$components[[nm]], names(directed))
    if (!length(hit)) next
    flags <- vapply(hit, function(h) word_to_flag(directed[[h]], h),
                    logical(1))
    if (length(unique(flags)) > 1) {
      stop("`edgelist_options$directed` gives conflicting directions for ",
           "network '", nm, "' via ", toString(hit), ".", call. = FALSE)
    }
    out[[nm]] <- flags[[1]]
  }
  out
}

#' Validate and fill in an `edgelist_options` list
#' @noRd
.validate_edgelist_options <- function(opts) {
  if (is.null(opts)) opts <- list()
  if (!is.list(opts)) {
    stop("`edgelist_options` must be a list.", call. = FALSE)
  }
  allowed <- c("edgelist_names", "edgelist_format", "edgelist_split",
               "nodelist", "missing", "directed")
  bad <- setdiff(names(opts), allowed)
  if (length(bad)) {
    stop("`edgelist_options` may only contain ", toString(allowed),
         "; got: ", toString(bad), ".", call. = FALSE)
  }
  fmt <- opts$edgelist_format %||% "long"
  fmt <- match.arg(fmt, c("long", "wide"))
  split_cols <- opts$edgelist_split
  if (!is.null(split_cols) && !is.character(split_cols)) {
    stop("`edgelist_options$edgelist_split` must be a character vector of ",
         "column names.", call. = FALSE)
  }
  list(edgelist_names = opts$edgelist_names,
       edgelist_format = fmt,
       edgelist_split = split_cols %||% character(0),
       nodelist = opts$nodelist,
       missing = opts$missing,
       directed = opts$directed)
}

#' Build one adjacency matrix from the edgelist rows of one network
#'
#' @param el The full edgelist.
#' @param rows Row indices belonging to this network.
#' @param cols Resolved sender/receiver/value column names.
#' @param value_col For wide format, the split column holding this network's
#'   values; `NULL` in long format (where `cols$value` is used instead).
#' @param nodes Node ids, in the final row order.
#' @param nm Network name, for messages.
#' @noRd
.el_one_matrix <- function(el, rows, cols, value_col, nodes, roster, nm,
                           printFlag = TRUE) {
  n <- length(nodes)
  mat <- matrix(0, n, n, dimnames = list(nodes, nodes))
  if (!length(rows)) return(mat)

  s <- as.character(el[[cols$sender]])[rows]
  r <- as.character(el[[cols$receiver]])[rows]
  v <- if (!is.null(value_col)) {
    as.numeric(el[[value_col]])[rows]
  } else if (!is.null(cols$value)) {
    as.numeric(el[[cols$value]])[rows]
  } else {
    rep(1, length(rows))
  }

  # wide format: a row is only an observation of THIS network where the
  # network's own value column is filled in
  if (!is.null(value_col)) {
    na_rows <- is.na(v)
  } else {
    na_rows <- rep(FALSE, length(v))
  }

  self <- s == r
  if (any(self, na.rm = TRUE)) {
    if (printFlag) {
      message("netimpute: dropped ", sum(self, na.rm = TRUE),
              " self-loop(s) from network '", nm,
              "' - the diagonal is not a tie.")
    }
    keep <- !self
    s <- s[keep]; r <- r[keep]; v <- v[keep]; na_rows <- na_rows[keep]
  }

  # Endpoints are checked against this network's ROSTER, not merely against
  # the full node set: a tie involving somebody the nodelist says was not in
  # this wave/question means the nodelist is wrong, and silently keeping the
  # tie (or silently dropping it) would hide that.
  outside <- !(s %in% roster) | !(r %in% roster)
  if (any(outside)) {
    bad <- utils::head(unique(paste0(s[outside], " -> ", r[outside])), 5)
    stop("Network '", nm, "' has tie(s) whose endpoints are not on its node ",
         "roster: ", toString(bad),
         if (sum(outside) > 5) paste0(", and ", sum(outside) - 5, " more"),
         ". Either `edgelist_options$nodelist` is missing those nodes, or ",
         "those rows belong to a different split.", call. = FALSE)
  }

  idx <- cbind(match(s, nodes), match(r, nodes))
  dup <- duplicated(idx)
  if (any(dup)) {
    d <- unique(paste0(s[dup], " -> ", r[dup]))
    stop("Network '", nm, "' has duplicate sender-receiver pair(s): ",
         toString(utils::head(d, 5)),
         if (length(d) > 5) paste0(", and ", length(d) - 5, " more"),
         ". Add the distinguishing column to `edgelist_split`, or aggregate ",
         "the duplicates before calling netmice().", call. = FALSE)
  }

  # NA-valued rows in wide format mark the cell unknown rather than absent
  if (any(na_rows)) mat[idx[na_rows, , drop = FALSE]] <- NA
  if (any(!na_rows)) {
    mat[idx[!na_rows, , drop = FALSE]] <- v[!na_rows]
  }
  diag(mat) <- 0
  mat
}

#' Apply the `missing` masks to one network
#'
#' A node whose outgoing ties are missing gets an all-NA row, EXCEPT the cells
#' the edgelist actually lists for it. The direction-matched warning below
#' matters: under actor non-response a non-reporting node still appears as a
#' *receiver* throughout the edgelist because others nominated it, so warning
#' on that would fire on essentially every real dataset.
#' @noRd
.el_apply_missing <- function(mat, miss, listed_out, listed_in, nm) {
  nodes <- rownames(mat)
  conflict_out <- intersect(miss$out, listed_out)
  conflict_in <- intersect(miss[["in"]], listed_in)
  if (length(conflict_out) || length(conflict_in)) {
    parts <- c(
      if (length(conflict_out))
        paste0("outgoing ties listed for ",
               toString(utils::head(conflict_out, 8))),
      if (length(conflict_in))
        paste0("incoming ties listed for ",
               toString(utils::head(conflict_in, 8))))
    warning(
      "netimpute: in network '", nm, "', ", paste(parts, collapse = "; "),
      ", yet those nodes are marked missing in that same direction. An ",
      "edgelist can only say which ties EXIST, so netimpute treats the ",
      "listed edges as observed and every other cell in those rows/columns ",
      "as missing - the zeros are not observed zeros. If some of those ",
      "zeros are genuinely observed, supply the data as ",
      "`edgelist_format = 'wide'` or as matrices instead.", call. = FALSE)
  }
  if (length(miss$out)) {
    i <- match(intersect(miss$out, nodes), nodes)
    i <- i[!is.na(i)]
    for (k in i) {
      keep <- !is.na(mat[k, ]) & mat[k, ] != 0
      mat[k, !keep] <- NA
    }
  }
  if (length(miss[["in"]])) {
    j <- match(intersect(miss[["in"]], nodes), nodes)
    j <- j[!is.na(j)]
    for (k in j) {
      keep <- !is.na(mat[, k]) & mat[, k] != 0
      mat[!keep, k] <- NA
    }
  }
  diag(mat) <- 0
  mat
}

#' Report the nodes an edgelist adds to `data`
#' @noRd
.el_warn_new_nodes <- function(connected, isolates, all_missing) {
  if (!length(c(connected, isolates, all_missing))) return(invisible(FALSE))
  fmt <- function(v) if (!length(v)) "none" else
    .firebrick(toString(utils::head(v, 20)))
  warning(
    "netimpute: ", length(c(connected, isolates, all_missing)),
    " node(s) appear in the edgelist/nodelist/missing arguments but not in ",
    "`data`. They have been added as new rows whose attributes are all ",
    "missing, and will be imputed from the network alone.\n",
    "  added with connections (", length(connected), "): ", fmt(connected), "\n",
    "  added as isolates (", length(isolates), "): ", fmt(isolates), "\n",
    "  added entirely missing (", length(all_missing), "): ",
    fmt(all_missing), "\n",
    "These nodes carry the most uncertainty in the imputations, since none ",
    "of their attributes were observed.", call. = FALSE)
  invisible(TRUE)
}

#' Convert an edgelist into netmice()'s aligned list of adjacency matrices
#'
#' @return A list with `mats` (named list of n x n matrices), `data` (possibly
#'   extended with rows for new nodes), `structural` (named list of logical
#'   matrices marking off-roster cells, or `NULL`), and `directed` (named
#'   logical vector, or `NULL` to infer).
#' @noRd
.networks_from_edgelist <- function(el, edgelist_options, data, id,
                                    printFlag = TRUE, what = "`networks`",
                                    allow_missing = TRUE) {
  if (is.null(id)) {
    stop(what, " is an edgelist, so `id` must name the column of `data` ",
         "holding the node identifiers used in its sender/receiver columns. ",
         "Row order alone cannot resolve node names.", call. = FALSE)
  }
  if (!id %in% names(data)) {
    stop("`id` column '", id, "' not found in `data`.", call. = FALSE)
  }
  el <- as.data.frame(el)
  opts <- .validate_edgelist_options(edgelist_options)
  keys <- .el_network_keys(el, opts$edgelist_format, opts$edgelist_split)
  cols <- .el_resolve_names(el, opts$edgelist_names, opts$edgelist_split,
                            opts$edgelist_format, printFlag, what)

  base_ids <- as.character(data[[id]])
  if (anyDuplicated(base_ids)) {
    stop("`id` column '", id, "' has duplicate values; node identifiers must ",
         "be unique.", call. = FALSE)
  }
  rosters <- .el_resolve_nodelist(opts$nodelist, data, keys, base_ids)
  miss_by_net <- if (allow_missing)
    .el_resolve_missing(opts$missing, keys) else
      .el_resolve_missing(NULL, keys)
  directed <- .el_resolve_directed(opts$directed, keys)

  el_nodes <- unique(c(as.character(el[[cols$sender]]),
                       as.character(el[[cols$receiver]])))
  roster_nodes <- unique(unlist(rosters, use.names = FALSE))
  miss_nodes <- unique(unlist(lapply(miss_by_net, function(d)
    c(d$out, d[["in"]])), use.names = FALSE))

  # nodes the edgelist knows about that `data` does not
  new_nodes <- setdiff(unique(c(el_nodes, roster_nodes, miss_nodes)), base_ids)
  if (length(new_nodes)) {
    connected <- intersect(new_nodes, el_nodes)
    rest <- setdiff(new_nodes, connected)
    all_missing <- intersect(rest, miss_nodes)
    isolates <- setdiff(rest, all_missing)
    .el_warn_new_nodes(connected, isolates, all_missing)
    extra <- data[rep(NA_integer_, length(new_nodes)), , drop = FALSE]
    extra[[id]] <- new_nodes
    rownames(extra) <- NULL
    data <- rbind(data, extra)
    base_ids <- as.character(data[[id]])
  }

  nodes <- base_ids                      # final row order, matching `data`
  mats <- stats::setNames(vector("list", length(keys$names)), keys$names)
  structural <- stats::setNames(vector("list", length(keys$names)), keys$names)
  any_struct <- FALSE

  for (nm in keys$names) {
    rows <- keys$rows[[nm]]
    value_col <- if (identical(opts$edgelist_format, "wide")) nm else NULL
    if (!is.null(value_col)) {
      # a row belongs to this network only where its value column is filled
      # in; NA there means "unknown", which .el_one_matrix records as NA
      rows <- rows
    }
    roster <- intersect(nodes, rosters[[nm]])
    # build on the FULL node set so every network is n x n, then mark the
    # off-roster cells structural
    m <- .el_one_matrix(el, rows, cols, value_col, nodes, roster, nm, printFlag)

    undirected <- !is.null(directed) && isFALSE(directed[[nm]])
    if (undirected) {
      # an undirected edgelist lists each edge once: mirror it, taking the
      # non-zero/non-NA side where only one direction was given
      tm <- t(m)
      filled <- (!is.na(tm) & tm != 0) & (is.na(m) | m == 0)
      m[filled] <- tm[filled]
    }

    listed_out <- unique(as.character(el[[cols$sender]])[rows])
    listed_in <- unique(as.character(el[[cols$receiver]])[rows])
    if (undirected) {
      both <- unique(c(listed_out, listed_in))
      listed_out <- both
      listed_in <- both
    }
    m <- .el_apply_missing(m, miss_by_net[[nm]], listed_out, listed_in, nm)

    off <- setdiff(nodes, roster)
    if (length(off)) {
      sm <- matrix(FALSE, length(nodes), length(nodes),
                   dimnames = list(nodes, nodes))
      sm[off, ] <- TRUE
      sm[, off] <- TRUE
      diag(sm) <- FALSE
      structural[[nm]] <- sm
      any_struct <- TRUE
      m[sm] <- 0
      if (printFlag) {
        message("netimpute: ", length(off), " node(s) are not on network '",
                nm, "'s roster; their ties there are fixed at 0 ",
                "(structurally absent), never imputed.")
      }
    }
    mats[[nm]] <- m
  }

  list(mats = mats,
       data = data,
       structural = if (any_struct) structural else NULL,
       directed = directed)
}
