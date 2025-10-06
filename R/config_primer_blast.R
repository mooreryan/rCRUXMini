new_primer_blast_params <- function(params = NULL) {
  params |>
    .validate_primer_blast_params() |>
    .apply_primer_blast_param_defaults() |>
    structure(class = "primer_blast_params")
}

assert_primer_blast_params <- function(object) {
  checkmate::assert_class(object, "primer_blast_params")
}

# This assumes that the names are the same that the BLAST program CLI expects.
primer_blast_params_to_cli_args <- function(primer_blast_params) {
  assert_primer_blast_params(primer_blast_params)

  param_names <- names(primer_blast_params)

  seq_along(primer_blast_params) |>
    purrr::map(function(i) {
      # Formats the name the way the blast CLI arg should be
      param_name <- stringr::str_glue("-{name}", name = param_names[[i]])
      # Don't force numeric values to scientific notation. E.g.,
      #
      # > format(10000000)
      # [1] "1e+07"
      # > format(10000000, scientific=FALSE)
      # [1] "10000000"
      param_value <- format(primer_blast_params[[i]], scientific = FALSE)
      c(param_name, param_value)
    }) |>
    purrr::list_c()
}

.validate_primer_blast_params <- function(params = NULL) {
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
      "evalue",
      "num_alignments",
      "num_threads",
      "perc_identity",
      "qcov_hsp_perc",
      "reward",
      "task",
      "word_size"
    )
  )

  # `null.ok = TRUE` is nice, but I want to ensure the user doesn't do something
  # like `evalue: NULL` in the yaml file. In other words, if the key/name is
  # present, then NULL is _not_ okay, but none of the keys/names are required to
  # be present.

  params_names <- names(params)

  if ("evalue" %in% params_names) {
    params$evalue |>
      checkmate::assert_number(finite = TRUE)
  }

  if ("num_alignments" %in% params_names) {
    params$num_alignments |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("num_threads" %in% params_names) {
    params$num_threads |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("perc_identity" %in% params_names) {
    params$perc_identity |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("qcov_hsp_perc" %in% params_names) {
    params$qcov_hsp_perc |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("reward" %in% params_names) {
    params$reward |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("task" %in% params_names) {
    params$task |>
      checkmate::assert_choice(
        # These are all the legal options for blastn command task
        choices = c(
          "blastn",
          "blastn-short",
          "dc-megablast",
          "megablast",
          "rmblastn"
        )
      )
  }

  if ("word_size" %in% params_names) {
    params$word_size |>
      checkmate::assert_count(positive = TRUE)
  }

  # If we get here, then we know the params is valid, so return it.
  params
}

.apply_primer_blast_param_defaults <- function(params = NULL) {
  checkmate::assert_list(
    params,
    # No missing values are allowed, i.e., no list(a=NULL)
    any.missing = FALSE,
    null.ok = TRUE
  )

  defaults <- list(
    evalue = 3e7,
    num_alignments = 10000000,
    num_threads = 1,
    perc_identity = 50,
    qcov_hsp_perc = 90,
    reward = 2,
    task = "blastn-short",
    word_size = 7
  )

  if (is.null(params)) {
    return(defaults)
  }

  # Merge config with defaults, with user-provided values taking precedence
  modifyList(defaults, params)
}
