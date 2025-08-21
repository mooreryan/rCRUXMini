empty_primer_blast_result <- function() {
  tibble::tibble(
    query_id = character(),
    subject_gi = character(),
    subject_accession_version = character(),
    mismatches = integer(),
    subject_alignment_start = integer(),
    subject_alignment_end = integer(),
    subject_taxonomy_ids = character()
  )
}

primer_blast_column_types <- list(
  "query_id" = readr::col_character(),
  "subject_gi" = readr::col_character(),
  "subject_accession_version" = readr::col_character(),
  "mismatches" = readr::col_integer(),
  "subject_alignment_start" = readr::col_integer(),
  "subject_alignment_end" = readr::col_integer(),
  "subject_taxonomy_ids" = readr::col_character()
)

run_primer_blast <- function(
  blastn_executable,
  query_path,
  db_paths,
  jobs = 1,
  # These arguments are named after their blastn arg names.
  word_size = 7,
  evalue = 3e7,
  num_alignments = 10000000,
  qcov_hsp_perc = 90,
  perc_identity = 50,
  reward = 2,
  num_threads = 1
) {
  checkmate::assert_string(query_path)
  checkmate::assert_file_exists(query_path)
  checkmate::assert_character(db_paths, min.len = 1)
  checkmate::assert_int(jobs, lower = 1)
  checkmate::assert_character(db_paths)

  # Note that using %f will make something like this:
  #   sprintf("%f", -1e10) => "-10000000000.000000"
  fmt <- function(x, fmt) sprintf(fmt, x)

  shared_blastn_arguments <- c(
    "-task",
    "blastn-short",
    "-outfmt",
    '"6 qseqid sgi saccver mismatch sstart send staxids"',
    "-word_size",
    word_size |> checkmate::assert_int(lower = 4) |> fmt("%d"),
    "-evalue",
    evalue |> checkmate::assert_number() |> fmt("%f"),
    "-num_alignments",
    num_alignments |> checkmate::assert_int(lower = 0) |> fmt("%d"),
    "-qcov_hsp_perc",
    qcov_hsp_perc |>
      checkmate::assert_number(lower = 0, upper = 100) |>
      fmt("%f"),
    "-perc_identity",
    perc_identity |>
      checkmate::assert_number(lower = 0, upper = 100) |>
      fmt("%f"),
    "-reward",
    reward |>
      checkmate::assert_int(lower = 0) |>
      fmt("%f"),
    "-num_threads",
    num_threads
  )

  # TODO: any way to check the blast db before running?
  # TODO: this falls back to serial execution on windows
  parallel::mclapply(
    X = db_paths,
    mc.preschedule = FALSE,
    mc.cores = jobs,
    FUN = function(db_path) {
      result <- processx::run(
        command = blastn_executable,
        args = append(
          c(
            "-db",
            db_path,
            "-query",
            query_path
          ),
          shared_blastn_arguments
        ),
        error_on_status = FALSE
      )

      # TODO: also write out some log file about any blast that failed

      if (rlang::is_na(result$status) || result$status != 0) {
        # The job failed. Either need to retry or do something else.
        empty_primer_blast_result()
      } else {
        # Try to read the BLAST data
        rlang::try_fetch(
          expr = {
            blast_hit_data <- result$stdout |>
              readr::read_tsv(
                col_names = names(primer_blast_column_types),
                col_types = primer_blast_column_types
              )
            readr::stop_for_problems(blast_hit_data)
            # NOTE: since we specify col_names and col_types, if there are no
            # rows, this will be the good default value.
            blast_hit_data
          },
          # If there are any problems reading the data, report the error and
          # return
          # the default value.
          error = function(condition) {
            # TODO: need real logging
            warning(condition)
            empty_primer_blast_result()
          }
        )
      }
    }
  ) |>
    # Calling bind_rows on a single data frame returns that data frame, so the
    # `length(db_paths) == 1` case will be fine
    dplyr::bind_rows() |>
    # TODO: pull this function out. it's repeated in the parse function too
    checkmate::assert_data_frame(
      # TODO: need to turn the case of 0 rows into a real error, log it, and
      # shutdown gracefully.
      min.rows = 1,
      col.names = names(primer_blast_column_types),
      # Checkmate needs the types in a different form than readr, so we don't
      # use the data from primer_blast_column_types directly.
      types = c(
        "character",
        "character",
        "character",
        "integer",
        "integer",
        "integer",
        "character"
      )
    )
}

# TODO: you will need to adapt the tests for the original function which parsed
#       a real file to work with this one that takes the data frame directly.
parse_primer_blast_results <- function(
  primer_blast_data,
  maximum_mismatches = 4
) {
  # TODO: does this lower value make sense?
  checkmate::assert_int(maximum_mismatches, lower = 0)
  checkmate::assert_data_frame(
    primer_blast_data,
    min.rows = 1,
    col.names = names(primer_blast_column_types),
    # Checkmate needs the types in a different form than readr, so we don't
    # use the data from primer_blast_column_types directly.
    types = c(
      "character",
      "character",
      "character",
      "integer",
      "integer",
      "integer",
      "character"
    )
  )

  result <- primer_blast_data |>
    # TODO: pretty sure the original code was confused about whether it should
    #       be 4 or 6 max mismatches.
    #
    # We only care about hits with less than 4 mismatches.
    dplyr::filter(.data$mismatch < maximum_mismatches) |>
    # Group by target sequence and target sequence start position
    dplyr::group_by(.data$saccver, .data$sstart) |>
    # Keep hits with the lowest number of mismatches
    dplyr::filter(.data$mismatch == min(.data$mismatch)) |>
    # Remove any remaining hits with duplicate saacver-sstart pairs, ignoring
    # any other differences like send.
    dplyr::slice(1)

  result <- result |>
    # Keep only the accessions that have both a forward and a reverse primer hit
    dplyr::group_by(.data$saccver) |>
    dplyr::filter(
      any(grepl("forward", .data$qseqid)) & any(grepl("reverse", .data$qseqid))
    ) |>
    # Remove duplicate end positions per accession
    dplyr::distinct(.data$saccver, .data$send, .keep_all = TRUE)

  # Ungroup and sort
  result <- result |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$saccver, .data$sstart, .data$mismatch)

  # TODO: checkmate on results
  checkmate::assert_data_frame(
    result,
    min.rows = 1,
    # The result is potentially removing rows, so it should never have more
    # rows than the starting data does.
    max.rows = nrow(primer_blast_data)
  )

  result
}
