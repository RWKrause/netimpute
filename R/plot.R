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
#' distance) across iterations, again one line per chain. Finally, for every
#' network that had missing ties, plots the mean and variance of its
#' *imputed tie values* - the tie-side analogue of the attribute page, and
#' the more sensitive of the two network diagnostics: the
#' \code{\link{net_diagnostics}} above describe the whole completed network,
#' so with little missingness they are dominated by the observed ties and can
#' look flat while the imputed portion drifts. Attributes, networks and
#' imputed ties are drawn as separate pages; call \code{par(ask = TRUE)}
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
#'   iteration). Pass \code{character(0)} to skip both network pages. The
#'   imputed-tie page additionally covers only those selected networks that
#'   actually had missing ties (\code{x$net_missing}).
#' @param ask Logical; prompt before each page of plots (default: only when
#'   the session is interactive, as for other multi-page base R plots).
#' @param max_rows Maximum number of grid rows per page (default 4). With
#'   many attributes or networks the plots are split across additional
#'   pages instead of being squeezed into ever-smaller panels (which
#'   eventually errors with "figure margins too large").
#' @param ... Passed on to the underlying \code{matplot()} calls, and
#'   overriding the defaults this method sets (\code{type}, \code{lty},
#'   \code{col}, \code{xlab}, \code{ylab}, \code{main}) where they clash.
#'
#' @return \code{x}, invisibly.
#' @seealso \code{\link{netmice}} to create the object;
#'   \code{\link{net_diagnostics}} for the network statistics plotted;
#'   \code{\link{print.netmids}} for a settings summary.
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
#' plot(fit)                       # attribute, network and imputed-tie pages
#' plot(fit, nets = character(0))  # attributes only
#' }
plot.netmids <- function(x, vars = NULL, nets = NULL, ask = interactive(), max_rows = 4L, ...) {
  if (!inherits(x, "netmids")) stop("`x` must be a netmids object from netmice().", call. = FALSE)

  vars <- if (is.null(vars)) x$var_missing else vars
  all_nets <- dimnames(x$netChain)[[1]]
  nets <- if (is.null(nets)) all_nets else nets
  bad_nets <- setdiff(nets, all_nets)
  if (length(bad_nets)) {
    stop("Unknown network name(s) in `nets`: ", toString(bad_nets), call. = FALSE)
  }
  bad_vars <- setdiff(vars, x$var_missing)
  if (length(bad_vars)) {
    stop("Unknown/fully-observed attribute name(s) in `vars`: ", toString(bad_vars),
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
  # Save par ONCE, before anything touches it, and restore it once on exit.
  # (Saving again between pages would restore the compact settings rather
  # than the caller's, leaving the device in a modified state.)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  m <- x$m
  maxit <- dim(x$chainMean)[2L]
  chain_cols <- grDevices::hcl.colors(max(m, 3), palette = "Dark 3")[seq_len(m)]
  dots <- list(...)

  # Compact panels: tight margins and scaled-down text, so even a full
  # max_rows-page fits comfortably on an ordinary device.
  compact_par <- function(nrow, ncol) {
    graphics::par(mfrow = c(nrow, ncol), mar = c(2.2, 2.4, 1.6, 0.5),
                  mgp = c(1.3, 0.4, 0), tcl = -0.3,
                  cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.75)
  }
  # Split a vector of row labels into pages of at most max_rows each.
  paginate <- function(items) split(items, ceiling(seq_along(items) / max_rows))

  # One trace panel. `y` is always coerced to an explicit maxit x m matrix:
  # `[v, , ]` drops to a plain vector when maxit == 1 or m == 1, and matplot()
  # would then read a length-m vector as ONE chain of m iterations - silently
  # plotting chains along the iteration axis. Arguments in `...` override the
  # defaults instead of colliding with them ("formal argument matched by
  # multiple actual arguments").
  trace_panel <- function(y, ylab, main) {
    y <- matrix(as.numeric(y), nrow = maxit, ncol = m)
    # matplot() dies with "need finite 'ylim' values" on an all-NA series.
    # That happens whenever a variance trace has nothing to vary: an
    # attribute with exactly ONE missing value, or a network with exactly one
    # imputed dyad, both give var() = NA at every iteration. Draw the panel
    # empty rather than taking the whole plot down.
    if (!any(is.finite(y))) {
      graphics::plot.new()
      graphics::plot.window(c(0, 1), c(0, 1))
      graphics::title(main = main, xlab = "Iteration", ylab = ylab)
      graphics::box()
      graphics::text(0.5, 0.5, "no finite values
(single imputed cell)",
                     cex = 0.7, col = "grey40")
      return(invisible(NULL))
    }
    args <- list(x = y, type = "l", lty = 1, col = chain_cols,
                 xlab = "Iteration", ylab = ylab, main = main)
    args[names(dots)] <- dots
    do.call(graphics::matplot, args)
  }
  chain_legend <- function() {
    graphics::legend("topright", legend = paste("chain", seq_len(m)),
                     col = chain_cols, lty = 1, cex = 0.6, bty = "n")
  }

  if (length(vars)) {
    for (page in paginate(vars)) {
      compact_par(length(page), 2)
      for (i in seq_along(page)) {
        v <- page[i]
        trace_panel(x$chainMean[v, , ], "mean", paste(v, "- mean"))
        if (i == 1) chain_legend()
        trace_panel(x$chainVar[v, , ], "variance", paste(v, "- variance"))
      }
    }
  }

  if (length(nets)) {
    diag_names <- dimnames(x$netChain)[[2]]
    for (page in paginate(nets)) {
      compact_par(length(page), length(diag_names))
      first <- TRUE
      for (nm in page) {
        for (dn in diag_names) {
          trace_panel(x$netChain[nm, dn, , ], dn, paste(nm, "-", dn))
          if (first) {
            chain_legend()
            first <- FALSE
          }
        }
      }
    }

    # Imputed-cell page: the same mean/variance layout as the attributes
    # page, because these ARE the mice-style imputed-value traces for ties.
    # Only networks that had missing ties appear; `nets = character(0)`
    # already skipped both network pages above.
    imp_nets <- intersect(nets, x$net_missing)
    for (page in paginate(imp_nets)) {
      compact_par(length(page), 2)
      for (i in seq_along(page)) {
        nm <- page[i]
        trace_panel(x$netImpMean[nm, , ], "mean",
                    paste(nm, "- imputed ties, mean"))
        if (i == 1) chain_legend()
        trace_panel(x$netImpVar[nm, , ], "variance",
                    paste(nm, "- imputed ties, variance"))
      }
    }
  }

  invisible(x)
}
