# Convergence trace plots for netmids objects, in the spirit of
# mice::plot.mids(): one line per imputation chain (m), across iterations,
# for each tracked quantity. Well-mixed chains look like unstructured,
# overlapping noise around a stable level; a trend or chains that fail to
# overlap indicate maxit should be increased.

#' Plot convergence diagnostics for a \code{netmids} object
#'
#' For every attribute in \code{x$var_missing}, plots the mean and the
#' variance of its imputed (not observed) values across iterations, one line
#' per imputation chain - exactly what \code{mice::plot.mids()} shows for an
#' ordinary \code{mids} object. For every network named in \code{nets},
#' plots each of the five \code{\link{net_diagnostics}} (density,
#' reciprocity, transitivity, isolate count, average inverse geodesic
#' distance) across iterations, again one line per chain. Attributes and
#' networks are drawn as separate pages; call \code{par(ask = TRUE)}
#' yourself (or rely on the \code{ask} argument) to step through them
#' interactively.
#'
#' As with any MCMC-style diagnostic, this is a visual aid, not a formal
#' test: chains that overlap with no trend across iterations are consistent
#' with (but do not prove) convergence, while chains that drift or fail to
#' overlap are clear evidence that \code{maxit} should be increased.
#'
#' @param x A \code{netmids} object from \code{\link{netmice}}.
#' @param vars Character vector of attribute names to plot (default: all of
#'   \code{x$var_missing}). Pass \code{character(0)} to skip the attributes
#'   page entirely.
#' @param nets Character vector of network names to plot (default: all
#'   networks in \code{x}, regardless of whether they had missing ties,
#'   since density/reciprocity/etc. are tracked for every network every
#'   iteration). Pass \code{character(0)} to skip the networks page.
#' @param ask Logical; prompt before each page of plots (default: only when
#'   the session is interactive, as for other multi-page base R plots).
#' @param max_rows Maximum number of grid rows per page (default 4). With
#'   many attributes or networks the plots are split across additional
#'   pages instead of being squeezed into ever-smaller panels (which
#'   eventually errors with "figure margins too large").
#' @param ... Passed on to the underlying \code{matplot()} calls.
#'
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 25
#' friends <- matrix(rbinom(n * n, 1, 0.15), n, n)
#' diag(friends) <- 0
#' off <- which(row(friends) != col(friends))
#' friends[sample(off, 30)] <- NA
#' attrs <- data.frame(age = rnorm(n, 35, 8),
#'                     gender = sample(c("F", "M"), n, TRUE))
#' attrs$age[sample(n, 4)] <- NA
#'
#' \donttest{
#' fit <- netmice(attrs, list(friends = friends), m = 2, maxit = 5,
#'                seed = 1, printFlag = FALSE)
#' plot(fit)                       # attribute and network pages
#' plot(fit, nets = character(0))  # attributes only
#' }
plot.netmids <- function(x, vars = NULL, nets = NULL, ask = interactive(), max_rows = 4L, ...) {
  if (!inherits(x, "netmids")) stop("`x` must be a netmids object from netmice().", call. = FALSE)

  vars <- if (is.null(vars)) x$var_missing else vars
  all_nets <- dimnames(x$netChain)[[1]]
  nets <- if (is.null(nets)) all_nets else nets
  bad_nets <- setdiff(nets, all_nets)
  if (length(bad_nets)) {
    stop("Unknown network name(s) in `nets`: ", paste(bad_nets, collapse = ", "), call. = FALSE)
  }
  bad_vars <- setdiff(vars, x$var_missing)
  if (length(bad_vars)) {
    stop("Unknown/fully-observed attribute name(s) in `vars`: ", paste(bad_vars, collapse = ", "),
         call. = FALSE)
  }
  if (length(vars) == 0 && length(nets) == 0) {
    message("netimpute: nothing to plot (no attributes or networks selected).")
    return(invisible(x))
  }

  max_rows <- as.integer(max_rows)
  if (length(max_rows) != 1 || is.na(max_rows) || max_rows < 1) {
    stop("`max_rows` must be a single positive integer.", call. = FALSE)
  }

  old_ask <- grDevices::devAskNewPage(ask)
  on.exit(grDevices::devAskNewPage(old_ask), add = TRUE)

  m <- x$m
  chain_cols <- grDevices::hcl.colors(max(m, 3), palette = "Dark 3")[seq_len(m)]

  # Compact panels: tight margins and scaled-down text, so even a full
  # max_rows-page fits comfortably on an ordinary device.
  compact_par <- function(nrow, ncol) {
    graphics::par(mfrow = c(nrow, ncol), mar = c(2.2, 2.4, 1.6, 0.5),
                  mgp = c(1.3, 0.4, 0), tcl = -0.3,
                  cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.75)
  }
  # Split a vector of row labels into pages of at most max_rows each.
  paginate <- function(items) split(items, ceiling(seq_along(items) / max_rows))

  if (length(vars)) {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    for (page in paginate(vars)) {
      compact_par(length(page), 2)
      for (i in seq_along(page)) {
        v <- page[i]
        graphics::matplot(x$chainMean[v, , ], type = "l", lty = 1, col = chain_cols,
                           xlab = "Iteration", ylab = "mean", main = paste(v, "- mean"), ...)
        if (i == 1) {
          graphics::legend("topright", legend = paste("chain", seq_len(m)), col = chain_cols,
                            lty = 1, cex = 0.6, bty = "n")
        }
        graphics::matplot(x$chainVar[v, , ], type = "l", lty = 1, col = chain_cols,
                           xlab = "Iteration", ylab = "variance", main = paste(v, "- variance"), ...)
      }
    }
  }

  if (length(nets)) {
    diag_names <- dimnames(x$netChain)[[2]]
    old_par2 <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par2), add = TRUE)
    for (page in paginate(nets)) {
      compact_par(length(page), length(diag_names))
      first <- TRUE
      for (nm in page) {
        for (dn in diag_names) {
          graphics::matplot(x$netChain[nm, dn, , ], type = "l", lty = 1, col = chain_cols,
                             xlab = "Iteration", ylab = dn, main = paste(nm, "-", dn), ...)
          if (first) {
            graphics::legend("topright", legend = paste("chain", seq_len(m)), col = chain_cols,
                              lty = 1, cex = 0.6, bty = "n")
            first <- FALSE
          }
        }
      }
    }
  }

  invisible(x)
}
