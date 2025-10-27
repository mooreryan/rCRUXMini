# Suppress .data warning from lintr
utils::globalVariables(".data")
utils::globalVariables("global_logger")
utils::globalVariables("rebind_global_logger")
utils::globalVariables("log_unknown")
utils::globalVariables("log_fatal")
utils::globalVariables("log_error")
utils::globalVariables("log_warning")
utils::globalVariables("log_info")
utils::globalVariables("log_debug")
utils::globalVariables("log_trace")
utils::globalVariables("abort_rcrux_mini_error")

global_logger <- list(
  console = littlelogger::new_logger(log_level = "info"),
  file = NULL
)
