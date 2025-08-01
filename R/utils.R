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

# Create standardized object/class for output rather than list???
run_command <- function(command, args = character(), timeout = 0) {
  checkmate::assert_string(command)
  checkmate::assert(
    checkmate::check_string(args),
    checkmate::check_character(args),
    checkmate::check_null(args),
    combine = "or"
  )
  checkmate::assert_int(timeout)
  checkmate::assert_true(timeout >= 0)

  # Create temporary files to capture stdout and stderr
  stdout_file <- tempfile()
  stderr_file <- tempfile()

  on.exit(
    {
      unlink(c(stdout_file, stderr_file))
    },
    add = TRUE
  )

  # Run the command
  warnings <- NULL
  exit_code <- withCallingHandlers(
    {
      system2(
        command = command,
        args = args,
        stdout = stdout_file,
        stderr = stderr_file,
        wait = TRUE,
        timeout = timeout
      )
    },
    warning = function(w) {
      warnings <<- w$message
      invokeRestart("muffleWarning")
    }
  )

  checkmate::assert_int(exit_code)

  stdout_message <- read_file(stdout_file)
  stderr_message <- read_file(stderr_file)

  result <- list(
    exit_code = exit_code,
    stdout = stdout_message,
    stderr = stderr_message,
    warnings = ifelse(rlang::is_null(warnings), "", warnings),
    full_command = paste(c(command, args), collapse = " ")
  )

  result
}
