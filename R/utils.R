enumerate_ambiguity <- function(dna_sequence) {
  checkmate::assert_string(dna_sequence)

  if (dna_sequence == "") {
    return("")
  }

  iupac <- list(
    "A" = "A",
    "C" = "C",
    "G" = "G",
    "T" = "T",
    "M" = list("A", "C"),
    "R" = list("A", "G"),
    "W" = list("A", "T"),
    "S" = list("C", "G"),
    "Y" = list("C", "T"),
    "K" = list("G", "T"),
    "V" = list("A", "C", "G"),
    "H" = list("A", "C", "T"),
    "D" = list("A", "G", "T"),
    "B" = list("C", "G", "T"),
    "N" = list("A", "C", "G", "T"),
    "I" = list("A", "T", "C")
  )

  # Split dna_sequence into individual characters
  chars <- strsplit(dna_sequence, "")[[1]]

  # Replace each character with its possible values
  possibilities <- lapply(chars, function(char) {
    if (char %in% names(iupac)) {
      iupac[[char]]
    } else {
      # Return as-is if not an IUPAC code
      char
    }
  })

  # Generate all combinations
  combinations <- expand.grid(possibilities, stringsAsFactors = FALSE)

  # Collapse each row into a sequence string
  result <- apply(combinations, 1, paste, collapse = "")

  checkmate::assert_character(result)

  result
}

enumerate_ambiguities <- function(dna_sequences) {
  if (checkmate::test_string(dna_sequences)) {
    return(enumerate_ambiguity(dna_sequences))
  }

  checkmate::assert_character(dna_sequences)

  result <- lapply(dna_sequences, enumerate_ambiguity) |> unlist()

  checkmate::assert_character(result)

  result
}


read_file <- function(file_name) {
  if (file.exists(file_name) && file.info(file_name)$size > 0) {
    result <- readLines(file_name, warn = FALSE) |>
      paste(collapse = "\n")
  } else {
    result <- ""
  }

  checkmate::assert_string(result)

  result
}

# TODO: need to replace this with processx::run
# Create standardized object/class for output rather than list???
# run_command <- function(command, args = character(), timeout = 0) {
#   checkmate::assert_string(command)
#   checkmate::assert(
#     checkmate::check_string(args),
#     checkmate::check_character(args),
#     checkmate::check_null(args),
#     combine = "or"
#   )
#   checkmate::assert_int(timeout)
#   checkmate::assert_true(timeout >= 0)

#   # Create temporary files to capture stdout and stderr
#   stdout_file <- tempfile()
#   stderr_file <- tempfile()

#   on.exit(
#     {
#       unlink(c(stdout_file, stderr_file))
#     },
#     add = TRUE
#   )

#   # Run the command
#   warnings <- NULL
#   exit_code <- withCallingHandlers(
#     {
#       system2(
#         command = command,
#         args = args,
#         stdout = stdout_file,
#         stderr = stderr_file,
#         wait = TRUE,
#         timeout = timeout
#       )
#     },
#     warning = function(w) {
#       warnings <<- w$message
#       invokeRestart("muffleWarning")
#     }
#   )

#   checkmate::assert_int(exit_code)

#   stdout_message <- read_file(stdout_file)
#   stderr_message <- read_file(stderr_file)

#   result <- list(
#     exit_code = exit_code,
#     stdout = stdout_message,
#     stderr = stderr_message,
#     warnings = ifelse(rlang::is_null(warnings), "", warnings),
#     full_command = paste(c(command, args), collapse = " ")
#   )

#   result
# }

to_fasta_string <- function(id, sequences) {
  id <- ifelse(id == "", "seq", id)

  strings <- lapply(
    seq_along(sequences),
    function(i) {
      stringr::str_glue(">{id}_{i}\n{sequences[[i]]}")
    }
  )

  paste(strings, collapse = "\n")
}

write_fasta <- function(
  df,
  to_filename,
  id_column,
  sequence_column
) {
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c(id_column, sequence_column)
  )

  fasta_lines <- df |>
    dplyr::mutate(
      header = paste0(">", .data[[id_column]]),
      entry = paste(.data$header, .data[[sequence_column]], sep = "\n")
    ) |>
    dplyr::pull(.data$entry)

  writeLines(fasta_lines, to_filename)
}

# TODO: I want to wrap all the checkmate errors in an rcrux_mini_assertion_error
# or something like that so it is easier for the user to catch the errors.
abort_rcrux_mini_error <- function(message, ...) {
  rlang::abort(message = message, class = "rcrux_mini_error", ...)
}

#' Find executable or abort with error
#'
#' Like \code{\link{sys_which}} but throws an error if the command is not found
#' instead of returning an empty string.
#'
#' @param command \code{[character(1)]}\cr
#'   The command name to look for.
#' @param ncbi_bin_directory \code{[character(1) | NULL]}\cr
#'   Optional directory to search in. If provided, searches only in this
#'   directory instead of PATH.
#'
#' @return A single string: the path to the command.
#'
#' @throws \code{rcrux_mini_error} if the command is not found.
#'
#'
sys_which_or_abort <- function(command, ncbi_bin_directory) {
  result <- SnailBLAST::sys_which(command, ncbi_bin_directory)

  if (result == "") {
    abort_rcrux_mini_error(
      stringr::str_glue("{command} doesn't exist"),
      command = command,
      ncbi_bin_directory = ncbi_bin_directory
    )
  }

  result
}

# TODO: this run_command is different from the old one, so you will need to check all `run_command` callsites and ensure they are compatible with this one

#' Like processx::run, except that it captures error and interrupt conditions and
#' rethrows them as \code{rcrux_mini_error} conditions.
#'
#' \code{error_on_status = TRUE} is explicitly provided to \code{processx::run},
#' so don't include it in the dots.
#'
#'
run_command <- function(command, args, ...) {
  rlang::try_fetch(
    processx::run(command = command, args = args, ..., error_on_status = TRUE),
    system_command_interrupt = function(condition) {
      abort_rcrux_mini_error(
        message = stringr::str_glue("Command '{command}' was interrupted'"),
        parent = condition
      )
    },
    system_command_timeout_error = function(condition) {
      abort_rcrux_mini_error(
        message = stringr::str_glue("Command '{command}' timed out'"),
        parent = condition
      )
    },
    system_command_status_error = function(condition) {
      # TODO: need to log errors

      abort_rcrux_mini_error(
        message = stringr::str_glue("Command '{command}' failed'"),
        parent = condition
      )
    }
  )
}
