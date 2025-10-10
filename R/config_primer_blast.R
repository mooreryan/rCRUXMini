#' Create a new primer BLAST configuration object
#'
#' Constructs a validated configuration object for primer BLAST operations.
#' User-provided parameters are validated and merged with sensible defaults.
#'
#' @param params A named list of BLAST parameters, or NULL. Valid parameter
#'   names include:
#'   \describe{
#'     \item{evalue}{Expectation value threshold (numeric). Default: 3e7}
#'     \item{num_alignments}{Maximum number of alignments to keep (positive integer). Default: 10000000}
#'     \item{num_threads}{Number of threads to use (positive integer). Default: 1}
#'     \item{perc_identity}{Minimum percent identity (positive integer). Default: 50}
#'     \item{qcov_hsp_perc}{Minimum query coverage per HSP (positive integer). Default: 90}
#'     \item{reward}{Reward for a nucleotide match (positive integer). Default: 2}
#'     \item{task}{BLAST task type (character). Must be one of: "blastn", "blastn-short",
#'       "dc-megablast", "megablast", or "rmblastn". Default: "blastn-short"}
#'     \item{word_size}{Word size for wordfinder algorithm (positive integer). Default: 7}
#'   }
#'
#' @return An object of class "rcrux_primer_blast_config" containing the
#'   validated and complete set of BLAST parameters.
#'
#' @examples
#' # Use all defaults
#' config <- new_primer_blast_config()
#'
#' # Override specific parameters
#' config <- new_primer_blast_config(params = list(
#'   num_threads = 4,
#'   perc_identity = 80
#' ))
#'
new_primer_blast_config <- function(params = NULL) {
  params |>
    .validate_primer_blast_params() |>
    .apply_primer_blast_param_defaults() |>
    structure(class = .primer_blast_config_class_name())
}

assert_primer_blast_params <- function(object) {
  checkmate::assert_class(object, .primer_blast_config_class_name())
}

.primer_blast_config_class_name <- function() {
  "rcrux_primer_blast_config"
}

# This assumes that the names are the same that the BLAST program CLI expects.
primer_blast_config_to_cli_args <- function(primer_blast_config) {
  assert_primer_blast_params(primer_blast_config)

  param_names <- names(primer_blast_config)

  seq_along(primer_blast_config) |>
    purrr::map(function(i) {
      # Formats the name the way the blast CLI arg should be
      param_name <- stringr::str_glue("-{name}", name = param_names[[i]])
      # Don't force numeric values to scientific notation. E.g.,
      #
      # > format(10000000)
      # [1] "1e+07"
      # > format(10000000, scientific=FALSE)
      # [1] "10000000"
      param_value <- format(primer_blast_config[[i]], scientific = FALSE)
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

  # TODO: validate percentages are between 0-100

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
      checkmate::assert_count(positive = TRUE)
  }

  if ("qcov_hsp_perc" %in% params_names) {
    params$perc_identity <- .parse_integer(params$perc_identity)

    params$qcov_hsp_perc |>
      checkmate::assert_count(positive = TRUE)
  }

  if ("reward" %in% params_names) {
    params$perc_identity <- .parse_integer(params$perc_identity)

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
