.blast_config_class_name <- function() {
  "rcrux_blast_config"
}

new_primer_blast_config <- function(params = NULL) {
  params |>
    .validate_blast_params() |>
    .apply_primer_blast_defaults() |>
    structure(class = .blast_config_class_name())
}

new_amplicon_blast_config <- function(params = NULL) {
  params |>
    .validate_blast_params() |>
    .apply_amplicon_blast_defaults() |>
    structure(class = .blast_config_class_name())
}

assert_blast_config <- function(object) {
  checkmate::assert_class(object, .blast_config_class_name())
}

# This assumes that the names are the same that the BLAST program CLI expects.
blast_config_to_cli_args <- function(blast_config) {
  assert_blast_config(blast_config)

  param_names <- names(blast_config)

  seq_along(blast_config) |>
    purrr::map(function(i) {
      # Formats the name the way the blast CLI arg should be
      param_name <- stringr::str_glue("-{name}", name = param_names[[i]])
      # Don't force numeric values to scientific notation. E.g.,
      #
      # > format(10000000)
      # [1] "1e+07"
      # > format(10000000, scientific=FALSE)
      # [1] "10000000"
      param_value <- format(blast_config[[i]], scientific = FALSE)
      c(param_name, param_value)
    }) |>
    purrr::list_c()
}

.validate_blast_params <- function(params = NULL) {
  checkmate::assert_list(
    params,
    # No missing values are allowed, i.e., no list(a=NULL)
    any.missing = FALSE,
    null.ok = TRUE
  )

  # `null.ok = TRUE` is nice, but I want to ensure the user doesn't do something
  # like `evalue: NULL` in the yaml file. In other words, if the key/name is
  # present, then NULL is _not_ okay, but none of the keys/names are required to
  # be present.

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

  params_names <- names(params)

  # NOTE: the parsing is done _again_ because the yaml reader treats certain
  # forms of scientific notation as characters, e.g., 3e7 will parse to "3e7"
  # string.

  if ("evalue" %in% params_names) {
    params$evalue <- .parse_numeric(params$evalue)

    params$evalue |>
      checkmate::assert_number(finite = TRUE)
  }

  if ("num_alignments" %in% params_names) {
    params$num_alignments <- .parse_integer(params$num_alignments)

    params$num_alignments |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("num_threads" %in% params_names) {
    params$num_threads <- .parse_integer(params$num_threads)

    params$num_threads |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("perc_identity" %in% params_names) {
    params$perc_identity <- .parse_integer(params$perc_identity)

    params$perc_identity |>
      checkmate::assert_int(lower = 0, upper = 100)
  }

  if ("qcov_hsp_perc" %in% params_names) {
    params$qcov_hsp_perc <- .parse_integer(params$qcov_hsp_perc)

    params$qcov_hsp_perc |>
      checkmate::assert_int(lower = 0, upper = 100)
  }

  if ("reward" %in% params_names) {
    params$reward <- .parse_integer(params$reward)

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
    params$word_size <- .parse_integer(params$word_size)

    params$word_size |>
      checkmate::assert_count(positive = TRUE)
  }

  # At this point, params are valid.
  params
}

.apply_primer_blast_defaults <- function(params = NULL) {
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

  modifyList(defaults, params)
}

# TODO: update amplicon blast defaults.
# Currently using the same defaults as primer blast as a placeholder.
.apply_amplicon_blast_defaults <- function(params = NULL) {
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

  modifyList(defaults, params)
}
