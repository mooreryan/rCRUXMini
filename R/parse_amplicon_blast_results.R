#' @param ambiguous_run_limit \code{integer(1)}\cr
#'   Sets the max allowed number of ambiguous bases in a row that a sequence
#'   can have before it is thrown out.
#'
#' @details We keep the longest alignment per query accession. The "longest" is
#' based on the original \code{alignment_length} as reported by BLAST. After
#' this, the alignment is degapped and the alignment length is recalculated. So
#' be aware that the sequences are filtered based on original reported alignment
#' length, but then degapped alignment length is used going forward. This is the
#' behavior that matches the original rCRUX.
parse_amplicon_blast_results <- function(
  amplicon_blast_result,
  ambiguous_run_limit,
  # We need this to "expand" the multi tax ID hits
  blastdbcmd,
  blast_db_paths,
  taxonomy_db_path
) {
  # TODO: checkmate the inputs
  #

  checkmate::assert_names(
    names(amplicon_blast_result),
    must.include = c(
      "query_accession",
      "subject_accession_version",
      "percent_identical_matches",
      "alignment_length",
      "expect_value",
      "subject_sequence_length",
      "subject_alignment_start",
      "subject_alignment_end",
      "subject_aligned_sequence",
      "unique_subject_taxonomy_ids"
    )
  )

  amplicon_blast_result <- amplicon_blast_result |>
    .add_degapped_subject_alignment()

  checkmate::assert_names(
    names(amplicon_blast_result),
    must.include = c(
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    )
  )

  amplicon_blast_result <- amplicon_blast_result |>
    .keep_longest_degapped_alignment_per_accession()

  #RYAN
  checkmate::assert_names(
    names(amplicon_blast_result),
    must.include = c(
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    )
  )

  # We want to report the most info to the user as possible for all the sequences. So we need to filter and add data in the correct order:
  #
  # - split single/multi taxids, expand multi taxids, join expanded to single
  # - add taxonomy info to all these
  # - add diagnostic info to all of these: maybe a problems column? or maybe a function that detects problems?
  #

  # tmp <- .split_by_ambiguous_sequences(
  #   amplicon_blast_result,
  #   ambiguous_run_limit = ambiguous_run_limit
  # )
  # if (nrow(tmp$ok) == 0) {
  #   abort_rcrux_mini_error(
  #     "There were no sequences remaining after filtering for ambiguous bases"
  #   )
  # }

  # amplicon_blast_result <- tmp$ok

  # # TODO: need to write these to a file or something
  # hits_with_too_many_ambiguous_bases <- tmp$ambiguous

  # NOTE: either of these may have 0 rows, which is not an error
  tmp <- .split_by_multiple_unique_subject_taxonomy_ids(amplicon_blast_result)
  hits_with_single_taxonomy_id <- tmp$single
  hits_with_multiple_taxonomy_ids <- tmp$multiple

  #RYAN
  checkmate::assert_names(
    names(hits_with_single_taxonomy_id),
    must.include = c(
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    )
  )

  #RYAN
  checkmate::assert_names(
    names(hits_with_multiple_taxonomy_ids),
    must.include = c(
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    )
  )

  # This is safe to run if the input has 0 rows, as long as it has the correct
  # column names.
  hits_with_multiple_taxonomy_ids_expanded <- expand_multi_tax_ids(
    blastdbcmd = blastdbcmd,
    hits_with_multiple_taxonomy_ids = hits_with_multiple_taxonomy_ids,
    blast_db_paths = blast_db_paths
  )

  #RYAN
  checkmate::assert_names(
    names(hits_with_multiple_taxonomy_ids_expanded),
    must.include = c(
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    )
  )

  if (
    !all.equal(
      colnames(hits_with_single_taxonomy_id),
      colnames(hits_with_multiple_taxonomy_ids_expanded)
    )
  ) {
    # TODO: this should be an assertion error
    abort_rcrux_mini_error(
      "hits_with_single_taxonomy_id and hits_with_multiple_taxonomy_ids_expanded should have identical column names"
    )
  }

  # These are the good hits ready to go
  hits <- dplyr::bind_rows(
    hits_with_single_taxonomy_id,
    hits_with_multiple_taxonomy_ids_expanded
  )

  #RYAN
  checkmate::assert_names(
    names(hits),
    must.include = c(
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    )
  )

  # This returns taxonomy info for the target/subject sequences that were hit.
  taxonomy_info_of_hits <- accession_to_taxonomy3(
    accessions = hits$subject_accession_version,
    taxonomy_db_path = taxonomy_db_path,
    # This will be the name of the accession column in the output. We can join
    # it later.
    accession_column_name = "subject_accession_version"
  )

  result <- dplyr::left_join(
    hits,
    taxonomy_info_of_hits,
    by = c("subject_accession_version")
  )

  # NOTE: right here, the original rCRUX does another fitler on the amplicon
  # length. However, amplicons that didn't have the correct length should never
  # even end up in the amplicon BLAST in the first place. So as long as you
  # don't call this function on random inputs, rather call it in the correct
  # spot in the pipeline and you won't have that problem. Not to mention the
  # amplicon length isn't even available here. In order to do it, we would need
  # to include the query sequence length in the blast output for the amplicon
  # blast.
  #
  # OHHH this is just an exapmle of the weird naming used by the original rCRUX.
  # They aren't filtering on the length of the query (ie the putative amplicon),
  # but rather they are using the same min and max amplicon length to filter the
  # degapped, alignment length of the hit. Okay, so I do need to do that.

  #RYAN
  checkmate::assert_names(
    names(result),
    must.include = c(
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    )
  )

  result
}

.keep_longest_degapped_alignment_per_accession <- function(blast_data) {
  checkmate::assert_data_frame(blast_data, min.rows = 1)
  checkmate::assert_names(
    names(blast_data),
    must.include = c("subject_accession_version", "degapped_alignment_length")
  )
  checkmate::assert_character(
    blast_data$subject_accession_version,
    any.missing = FALSE
  )
  checkmate::assert_integer(
    blast_data$degapped_alignment_length,
    any.missing = FALSE,
    lower = 0
  )

  blast_data |>
    dplyr::group_by(.data$subject_accession_version) |> # nolint: object_usage_linter
    dplyr::slice_max(
      .data$degapped_alignment_length,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()
}

# Adds degapped_subject_aligned_sequence and degapped_alignment_length
#
# Empty data frames are fine for this function
.add_degapped_subject_alignment <- function(blast_data) {
  checkmate::assert_data_frame(blast_data)
  checkmate::assert_names(
    names(blast_data),
    must.include = c("subject_aligned_sequence", "alignment_length")
  )
  checkmate::assert_character(
    blast_data$subject_aligned_sequence,
    any.missing = FALSE
  )

  blast_data |>
    dplyr::mutate(
      degapped_subject_aligned_sequence = stringr::str_remove_all(
        .data$subject_aligned_sequence,
        stringr::fixed("-")
      ),
      degapped_alignment_length = stringr::str_length(
        .data$subject_aligned_sequence
      )
    )
}

#' Split BLAST hits by runs of ambiguous bases in the subject-aligned sequence
#'
#' Splits a BLAST result data frame into two data frames based on whether
#' the subject_aligned_sequence contains a contiguous run of ambiguous bases
#' (the letter N, case-insensitive) longer than the allowed limit.
#'
#' @param blast_data data.frame
#'   BLAST results containing at least a character column named
#'   subject_aligned_sequence (no missing values).
#' @param ambiguous_run_limit integer(1)
#'   Maximum allowed number of ambiguous bases in a row. A value of 0 means
#'   no ambiguous bases are allowed (i.e., any sequence containing at least
#'   one N will be flagged as ambiguous).
#'
#' @return list
#'   A named list with two data frames:
#'   - ok: rows whose subject_aligned_sequence does not contain a run of N
#'     longer than ambiguous_run_limit (i.e., kept).
#'   - ambiguous: rows whose subject_aligned_sequence contains a run of N
#'     of length ambiguous_run_limit + 1 or greater (i.e., removed).
#'
#'   Each data frame preserves the original columns and row order; either
#'   element may be empty (0 rows).
#'
#' @details
#' Detection is performed using a fixed, case-insensitive pattern consisting
#' of the letter N repeated ambiguous_run_limit + 1 times. Only the character
#' N is treated as an ambiguous base here.
#'
#' @examples
#' # tibble::tibble(subject_aligned_sequence = c("ACGT", "ACNGT", "ACGNNN")) |>
#' #   .split_by_ambiguous_sequences(ambiguous_run_limit = 1)
#'
#' @seealso
#' .remove_ambiguous_sequences()
#'
#' @noRd
#'
#'
.split_by_ambiguous_sequences <- function(blast_data, ambiguous_run_limit) {
  checkmate::assert_data_frame(blast_data, min.rows = 1)
  checkmate::assert_names(
    names(blast_data),
    must.include = "subject_aligned_sequence"
  )
  checkmate::assert_character(
    blast_data$subject_aligned_sequence,
    any.missing = FALSE
  )

  checkmate::assert_int(ambiguous_run_limit, lower = 0)

  pattern <- stringr::fixed(
    stringr::str_dup("N", ambiguous_run_limit + 1),
    ignore_case = TRUE
  )

  too_many_ambiguous <- stringr::str_detect(
    string = blast_data$subject_aligned_sequence,
    pattern = pattern
  )

  list(
    ok = blast_data[!too_many_ambiguous, , drop = FALSE],
    ambiguous = blast_data[too_many_ambiguous, , drop = FALSE]
  )
}

#' Split BLAST hits by presence of multiple unique subject taxonomy IDs
#'
#' Splits a BLAST result data frame into two data frames based on whether
#' the unique_subject_taxonomy_ids column contains a semicolon, indicating
#' multiple taxonomy IDs.
#'
#' @param blast_data data.frame
#'   BLAST results containing at least a character column named
#'   unique_subject_taxonomy_ids (no missing values).
#'
#' @return list
#'   A named list with two data frames:
#'   - single: rows with a single taxonomy ID (no semicolon present).
#'   - multiple: rows with multiple taxonomy IDs (semicolon present).
#'
#'   Each data frame preserves the original columns and row order; either
#'   element may be empty (0 rows).
#'
#' @examples
#' # tibble::tibble(unique_subject_taxonomy_ids = c("12345", "12345;67890")) |>
#' #   .split_by_multiple_unique_subject_taxonomy_ids()
#'
#' @noRd
#'
#'
.split_by_multiple_unique_subject_taxonomy_ids <- function(blast_data) {
  checkmate::assert_data_frame(blast_data, min.rows = 1)
  checkmate::assert_names(
    names(blast_data),
    must.include = "unique_subject_taxonomy_ids"
  )
  checkmate::assert_character(
    blast_data$unique_subject_taxonomy_ids,
    any.missing = FALSE
  )

  has_multiple_ids <- stringr::str_detect(
    string = blast_data$unique_subject_taxonomy_ids,
    pattern = stringr::fixed(";")
  )

  list(
    single = blast_data[!has_multiple_ids, , drop = FALSE],
    multiple = blast_data[has_multiple_ids, , drop = FALSE]
  )
}

# TODO: might be good to save original gapped alignment and gapped alignment length

.check_for_long_ambiguous_runs <- function(
  aligned_sequences,
  ambiguous_run_limit
) {
  checkmate::assert_character(aligned_sequences, any.missing = FALSE)
  checkmate::assert_int(ambiguous_run_limit, lower = 0)

  pattern <- stringr::fixed(
    stringr::str_dup("N", ambiguous_run_limit + 1),
    ignore_case = TRUE
  )

  stringr::str_detect(
    string = aligned_sequences,
    pattern = pattern
  )
}

.check_for_missing_taxonomy <- function(taxonomy_ids) {
  checkmate::assert_integer(taxonomy_ids)
  taxonomy_ids |> purrr::map_lgl(.f = rlang::is_na)
}

# TODO: this should be degapped
.check_for_bad_subject_alignment_length <- function(
  subject_alignment_lengths,
  minimum_length,
  maximum_length
) {
  subject_alignment_lengths |>
    purrr::map_lgl(.f = function(subject_alignment_length) {
      subject_alignment_length < minimum_length ||
        subject_alignment_length > maximum_length
    })
}
