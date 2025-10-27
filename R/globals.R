# Suppress .data warning from lintr
utils::globalVariables(".data")

global_logger <- list(
  console = littlelogger::new_logger(log_level = "info"),
  file = NULL
)
