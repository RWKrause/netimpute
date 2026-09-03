# Package-local mutable state. Used only for "tell the user this once per
# session" notices, never for anything the results depend on - a chain must
# behave identically whether or not a notice has already fired.
.netimpute_env <- new.env(parent = emptyenv())

#' Emit a message at most once per session
#'
#' @param key Identifier for the notice; the first call with a given key
#'   prints, every later one is silent.
#' @param ... Passed to `message()`.
#' @return `TRUE` if the message was emitted, `FALSE` if it was suppressed.
#' @noRd
.notify_once <- function(key, ...) {
  if (isTRUE(.netimpute_env[[key]])) return(invisible(FALSE))
  assign(key, TRUE, envir = .netimpute_env)
  message(...)
  invisible(TRUE)
}

#' Reset the once-per-session notices (test support)
#' @noRd
.reset_notices <- function() {
  rm(list = ls(envir = .netimpute_env), envir = .netimpute_env)
  invisible(NULL)
}
