# Very simple result struct. Doesn't even use s3, but that might make things
# clearer. We'll see.

result__ok <- function(value) {
  list(status = "ok", value = value)
}

result__error <- function(value) {
  list(status = "error", value = value)
}

result__is_ok <- function(result) {
  result$status == "ok"
}

result__is_error <- function(result) {
  result$status == "error"
}

result__unwrap <- function(result, default) {
  if (result__is_ok(result)) {
    result$value
  } else {
    default
  }
}

result__partition <- function(results) {
  list(
    ok = purrr::keep(results, result__is_ok),
    error = purrr::keep(results, result__is_error)
  )
}
