#' Create a new general configuration object
#'
#' Constructs a validated configuration object for general pipeline parameters.
#' User-provided parameters are validated and merged with sensible defaults.
#'
#' @param params A named list of general parameters, or NULL. Valid parameter
#'   names include:
#'   \describe{
#'     \item{forward_primers}{Forward primer sequence(s) (character vector, required). Non-empty string(s).}
#'     \item{reverse_primers}{Reverse primer sequence(s) (character vector, required). Non-empty string(s).}
#'     \item{output_directory}{Path to output directory (character, required). Must not already exist.}
#'     \item{taxonomy_db_path}{Path to taxonomy database (character, required). Must exist.}
#'     \item{ncbi_bin_directory}{Path to NCBI bin directory (character, optional). Must exist if provided. Default: NULL}
#'     \item{blast_db_paths}{Path(s) to BLAST database(s) (character vector, required). Must refer to a valid BLAST db.}
#'     \item{query_chunk_count}{Number of query chunks for parallel processing (positive integer, optional). Default: 1}
#'   }
#'
#' @return An object of class "rcrux_general_config" containing the
#'   validated and complete set of general parameters.
#'
#' @examples
#' # With required parameters and some optional ones
#' config <- new_general_config(params = list(
#'   forward_primers = "ACGTACGT",
#'   reverse_primers = "TGCATGCA",
#'   output_directory = "/path/to/output",
#'   taxonomy_db_path = "/path/to/taxonomy.db",
#'   blast_db_paths = "/path/to/blast_db",
#'   query_chunk_count = 4
#' ))
#'
new_general_config <- function(params = NULL) {
  params |>
    .apply_general_param_defaults() |>
    .validate_general_params()
}

.validate_general_params <- function(params = NULL) {
  checkmate::assert_list(
    params,
    # No missing values are allowed, i.e., no list(a=NULL)
    any.missing = FALSE,
    null.ok = FALSE
  )

  if (is.null(params)) {
    return(params)
  }

  required_top_level_params <- c(
    "forward_primers",
    "reverse_primers",
    "output_directory",
    "taxonomy_db_path",
    "blast_db_paths"
  )

  optional_top_level_params <- c(
    "ncbi_bin_directory",
    "query_chunk_count",
    # Users aren't required to specify params for the nested params.
    "primer_blast",
    "plausible_amplicons"
  )

  checkmate::assert_names(
    names(params),
    must.include = required_top_level_params
  )

  checkmate::assert_names(
    names(params),
    subset.of = c(
      required_top_level_params,
      optional_top_level_params
    )
  )

  params_names <- names(params)

  # Validate forward_primers
  if ("forward_primers" %in% params_names) {
    params$forward_primers |>
      checkmate::assert_character(
        min.len = 1,
        any.missing = FALSE,
        min.chars = 1
      )
  }

  # Validate reverse_primers
  if ("reverse_primers" %in% params_names) {
    params$reverse_primers |>
      checkmate::assert_character(
        min.len = 1,
        any.missing = FALSE,
        min.chars = 1
      )
  }

  # Validate output_directory (must NOT exist)
  if ("output_directory" %in% params_names) {
    params$output_directory |>
      checkmate::assert_string(min.chars = 1)

    if (
      file.exists(params$output_directory) ||
        dir.exists(params$output_directory)
    ) {
      abort_rcrux_mini_error(
        stringr::str_glue(
          "Output directory already exists: {path}",
          path = params$output_directory
        )
      )
    }
  }

  # Validate taxonomy_db_path (must exist)
  if ("taxonomy_db_path" %in% params_names) {
    params$taxonomy_db_path |>
      checkmate::assert_string(min.chars = 1)

    if (!file.exists(params$taxonomy_db_path)) {
      abort_rcrux_mini_error(
        stringr::str_glue(
          "Taxonomy database path does not exist: {path}",
          path = params$taxonomy_db_path
        )
      )
    }
  }

  # Validate ncbi_bin_directory (optional, must exist if provided)
  if ("ncbi_bin_directory" %in% params_names) {
    # The default value in this case is for the user to NOT provide an bin
    # directory. But we use `""` to represent that.
    if (
      !is.null(params$ncbi_bin_directory) && params$ncbi_bin_directory != ""
    ) {
      params$ncbi_bin_directory |>
        checkmate::assert_string(min.chars = 1)

      if (!dir.exists(params$ncbi_bin_directory)) {
        abort_rcrux_mini_error(
          stringr::str_glue(
            "NCBI bin directory does not exist: {path}",
            path = params$ncbi_bin_directory
          )
        )
      }
    }
  }

  # Validate blast_db_paths (must exist)
  if ("blast_db_paths" %in% params_names) {
    params$blast_db_paths |>
      checkmate::assert_character(
        min.len = 1,
        any.missing = FALSE,
        min.chars = 1
      )

    # TODO: This won't work because blast db paths aren't actually real paths.
    # They point to the path and name. Check snail blast for code to ensure these
    # exist.
    #
    # Check that all paths exist
    # missing_paths <- params$blast_db_paths[!file.exists(params$blast_db_paths)]
    # if (length(missing_paths) > 0) {
    #   abort_rcrux_mini_error(
    #     stringr::str_glue(
    #       "BLAST database path(s) do not exist: {paths}",
    #       paths = paste(missing_paths, collapse = ", ")
    #     )
    #   )
    # }
  }

  # Validate query_chunk_count
  if ("query_chunk_count" %in% params_names) {
    params$query_chunk_count <- .parse_integer(params$query_chunk_count)

    params$query_chunk_count |>
      checkmate::assert_count(positive = TRUE)
  }

  # If we get here, then we know the params is valid, so return it.
  params
}

.apply_general_param_defaults <- function(params = NULL) {
  checkmate::assert_list(
    params,
    # No missing values are allowed, i.e., no list(a=NULL)
    any.missing = FALSE,
    null.ok = TRUE
  )

  defaults <- list(
    query_chunk_count = 1,
    # NOTE: we use an empty string rather than NULL here, because we don't allow
    # users to enter NULL-like values as valid paramters in the yaml files.
    ncbi_bin_directory = ""
  )

  if (is.null(params)) {
    return(defaults)
  }

  # Merge config with defaults, with user-provided values taking precedence
  modifyList(defaults, params)
}
