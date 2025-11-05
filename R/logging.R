# This is a facade over littlelogger. It will make it easier if the logger needs
# to change later.

#' Rebind the global logger, e.g., to take into account values taken from the
#' config file.
#'
#' This function should only be called _once_, and near the top of the
#' \code{pipeline} function after the directory in which you will be logging to
#' has been created. IGNORE THIS ADVICE AT YOUR OWN PERIL!
#'
rebind_global_logger <- function(
  console_level,
  file_level = NULL,
  file_path = NULL
) {
  if (is.null(file_level) && is.null(file_path)) {
    file_logger <- NULL
  } else if (!is.null(file_level) && !is.null(file_path)) {
    file_logger <- littlelogger::new_logger(
      log_level = file_level,
      file = file_path
    )
  } else {
    stop(
      "file_level and file path must both be given, or neither should be given"
    )
  }

  if (exists("global_logger", envir = asNamespace("rCRUXMini"))) {
    unlockBinding("global_logger", asNamespace("rCRUXMini"))
  }

  assign(
    "global_logger",
    value = list(
      console = littlelogger::new_logger(log_level = console_level),
      file = file_logger
    ),
    envir = asNamespace("rCRUXMini")
  )

  lockBinding("global_logger", asNamespace("rCRUXMini"))
}

log_unknown <- function(...) {
  if (!is.null(global_logger$console)) {
    global_logger$console$unknown(...)
  }

  if (!is.null(global_logger$file)) {
    global_logger$file$unknown(...)
  }
}

log_fatal <- function(...) {
  if (!is.null(global_logger$console)) {
    global_logger$console$fatal(...)
  }

  if (!is.null(global_logger$file)) {
    global_logger$file$fatal(...)
  }
}

log_error <- function(...) {
  if (!is.null(global_logger$console)) {
    global_logger$console$error(...)
  }

  if (!is.null(global_logger$file)) {
    global_logger$file$error(...)
  }
}

log_warning <- function(...) {
  if (!is.null(global_logger$console)) {
    global_logger$console$warning(...)
  }

  if (!is.null(global_logger$file)) {
    global_logger$file$warning(...)
  }
}

log_info <- function(...) {
  if (!is.null(global_logger$console)) {
    global_logger$console$info(...)
  }

  if (!is.null(global_logger$file)) {
    global_logger$file$info(...)
  }
}

log_debug <- function(...) {
  if (!is.null(global_logger$console)) {
    global_logger$console$debug(...)
  }

  if (!is.null(global_logger$file)) {
    global_logger$file$debug(...)
  }
}

log_trace <- function(...) {
  if (!is.null(global_logger$console)) {
    global_logger$console$trace(...)
  }

  if (!is.null(global_logger$file)) {
    global_logger$file$trace(...)
  }
}
