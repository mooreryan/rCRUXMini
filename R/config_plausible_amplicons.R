new_plausible_amplicons_params <- function(params = NULL) {
  params |>
    .validate_plausible_amplicons_params() |>
    .apply_plausible_amplicons_param_defaults() |>
    structure(class = "plausible_amplicons_params")
}

assert_plausible_amplicons_params <- function(object) {
  checkmate::assert_class(object, "plausible_amplicons_params")
}

.validate_plausible_amplicons_params <- function(params = NULL) {
  checkmate::assert_list(
    params,
    # No missing values are allowed, i.e., no list(a=NULL)
    any.missing = FALSE,
    null.ok = TRUE
  )

  if (is.null(params)) {
    return(params)
  }

  checkmate::assert_names(
    names(params),
    # The user doesn't have to include all of these, or any of them
    subset.of = c(
      "minimum_length",
      "maximum_length",
      "maximum_mismatches",
      "ambiguous_run_limit"
    )
  )

  # `null.ok = TRUE` is nice, but I want to ensure the user doesn't do something
  # like `minimum_length: NULL` in the yaml file. In other words, if the key/name is
  # present, then NULL is _not_ okay, but none of the keys/names are required to
  # be present.

  params_names <- names(params)

  if ("minimum_length" %in% params_names) {
    params$minimum_length |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("maximum_length" %in% params_names) {
    params$maximum_length |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("maximum_mismatches" %in% params_names) {
    params$maximum_mismatches |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("ambiguous_run_limit" %in% params_names) {
    params$ambiguous_run_limit |>
      checkmate::assert_count()
  }

  # Cross-parameter validation: ensure minimum_length < maximum_length
  if (all(c("minimum_length", "maximum_length") %in% params_names)) {
    if (params$minimum_length >= params$maximum_length) {
      # TODO: raise the same type of error as the other validations
      stop("minimum_length must be less than maximum_length")
    }
  }

  # If we get here, then we know the params is valid, so return it.
  params
}

.apply_plausible_amplicons_param_defaults <- function(params = NULL) {
  checkmate::assert_list(
    params,
    # No missing values are allowed, i.e., no list(a=NULL)
    any.missing = FALSE,
    null.ok = TRUE
  )

  defaults <- list(
    minimum_length = 150,
    maximum_length = 650,
    maximum_mismatches = 4,
    ambiguous_run_limit = 5
  )

  if (is.null(params)) {
    return(defaults)
  }

  # Merge config with defaults, with user-provided values taking precedence
  modifyList(defaults, params)
}
