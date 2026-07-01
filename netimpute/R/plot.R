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
#' @param ... Passed on to the underlying \code{matplot()} calls.
#'
#' @return \code{x}, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- netmice(attrs_miss, list(friends = friends_miss), m = 3, maxit = 20)
#' plot(fit)                       # both pages
#' plot(fit, nets = character(0))  # attributes only
#' }
plot.netmids <- function(x, vars = NULL, nets = NULL, ask = interactive(), ...) {
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

  old_ask <- grDevices::devAskNewPage(ask)
  on.exit(grDevices::devAskNewPage(old_ask), add = TRUE)

  m <- x$m
  chain_cols <- grDevices::hcl.colors(max(m, 3), palette = "Dark 3")[seq_len(m)]

  if (length(vars)) {
    old_par <- graphics::par(mfrow = c(length(vars), 2), mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
    on.exit(graphics::par(old_par), add = TRUE)
    for (i in seq_along(vars)) {
      v <- vars[i]
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

  if (length(nets)) {
    diag_names <- dimnames(x$netChain)[[2]]
    old_par2 <- graphics::par(mfrow = c(length(nets), length(diag_names)),
                               mar = c(3, 3, 2, 1), mgp = c(1.8, 0.6, 0))
    on.exit(graphics::par(old_par2), add = TRUE)
    first <- TRUE
    for (nm in nets) {
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

  invisible(x)
}
