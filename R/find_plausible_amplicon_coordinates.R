#' Calculate amplicon product length based on primer positions
#'
#' Note that the product length represents the full PCR product length which
#' includes the primers.
#'
#' TODO: these terms are pretty ambiguous...is it the length of the PCR product, or is it the part with the primers stripped off?
#'
#' PCR product includes the primers
#' target sequence is the good stuff between the primers
#' amplicon is ambiguous
#'
#' @param forward_start Numeric vector of forward primer start positions
#' @param forward_stop Numeric vector of forward primer stop positions
#' @param reverse_start Numeric vector of reverse primer start positions
#' @param reverse_stop Numeric vector of reverse primer stop positions
#'
#' @return Numeric vector of product lengths (NA values represent invalid orientations)
#' @keywords internal
.calculate_product_length <- function(
  forward_start,
  forward_stop,
  reverse_start,
  reverse_stop
) {
  checkmate::assert_numeric(forward_start, any.missing = FALSE)
  checkmate::assert_numeric(forward_stop, any.missing = FALSE)
  checkmate::assert_numeric(reverse_start, any.missing = FALSE)
  checkmate::assert_numeric(reverse_stop, any.missing = FALSE)
  checkmate::assert_true(
    length(forward_start) == length(forward_stop) &&
      length(forward_start) == length(reverse_start) &&
      length(forward_start) == length(reverse_stop)
  )

  # TODO: is this product length is off by 2 (off by one on both sides?)
  dplyr::case_when(
    #  F ---------->                        <------------ R  Primers
    # ====================================================== DNA Target
    #    ---------->************************<------------    Product
    #    ^ Forward start                  Reverse start ^
    (forward_start < reverse_start &
      forward_start < forward_stop &
      reverse_stop < reverse_start) ~
      (reverse_start - forward_start),
    #  R ---------->                        <------------ F  Primers
    # ====================================================== DNA Target
    #    ---------->************************<------------    Product
    #    ^ Reverse start                  Forward start ^
    #
    (forward_start > reverse_start &
      forward_start > forward_stop &
      reverse_stop > reverse_start) ~
      (forward_start - reverse_start),
  )
}

.find_all_product_lengths <- function(forward_hits, reverse_hits) {
  checkmate::assert_data_frame(
    forward_hits,
    any.missing = FALSE,
    types = c(
      "character",
      "character",
      "character",
      "numeric",
      "numeric",
      "numeric"
    )
  )
  checkmate::assert_names(
    colnames(forward_hits),
    must.include = c(
      "subject_accession_version",
      "subject_gi",
      "unique_subject_taxonomy_ids",
      "forward_start",
      "forward_stop",
      "forward_mismatch"
    )
  )

  checkmate::assert_data_frame(
    reverse_hits,
    any.missing = FALSE,
    types = c(
      "character",
      "character",
      "character",
      "numeric",
      "numeric",
      "numeric"
    )
  )
  checkmate::assert_names(
    colnames(reverse_hits),
    must.include = c(
      "subject_accession_version",
      "subject_gi",
      "unique_subject_taxonomy_ids",
      "reverse_start",
      "reverse_stop",
      "reverse_mismatch"
    )
  )

  result <- dplyr::inner_join(
    forward_hits,
    reverse_hits,
    # TODO: why bother including staxids here? Seems redundant.
    by = c(
      "subject_accession_version",
      "subject_gi",
      "unique_subject_taxonomy_ids"
    ),
    relationship = "many-to-many"
  ) |>
    dplyr::mutate(
      product_length = .calculate_product_length(
        forward_start = .data$forward_start,
        forward_stop = .data$forward_stop,
        reverse_start = .data$reverse_start,
        reverse_stop = .data$reverse_stop
      )
    )

  checkmate::assert_names(
    names(result),
    must.include = c(
      "subject_accession_version",
      "subject_gi",
      "unique_subject_taxonomy_ids",
      "forward_start",
      "forward_stop",
      "forward_mismatch",
      "reverse_start",
      "reverse_stop",
      "reverse_mismatch",
      "product_length"
    )
  )

  result
}

# NOTE: this returns the location of the amplicon on the target sequence, NOT
#       the actual sequence itself
#
find_plausible_amplicon_coordinates <- function(
  primer_blast_results,
  minimum_length,
  maximum_length,
  maximum_mismatches = 4
) {
  checkmate::assert_data_frame(
    primer_blast_results,
    min.rows = 1,
    any.missing = FALSE,
    types = c(
      "character",
      "character",
      "character",
      "numeric",
      "numeric",
      "numeric",
      "character"
    )
  )
  checkmate::assert_names(
    colnames(primer_blast_results),
    must.include = c(
      "qseqid",
      "sgi",
      "saccver",
      "mismatch",
      "sstart",
      "send",
      "staxids"
    )
  )
  checkmate::assert_int(minimum_length, lower = 1, upper = maximum_length)
  checkmate::assert_int(maximum_length, lower = minimum_length)
  checkmate::assert_int(maximum_mismatches, lower = 0)

  filtered_primer_blast_results <- primer_blast_results |>
    dplyr::filter(.data$mismatch <= maximum_mismatches)

  # Separate forward and reverse hits
  forward_hits <- filtered_primer_blast_results |>
    dplyr::filter(stringr::str_detect(.data$qseqid, "forward")) |>
    dplyr::select(
      subject_accession_version = "saccver",
      subject_gi = "sgi",
      unique_subject_taxonomy_ids = "staxids",
      forward_start = "sstart",
      forward_stop = "send",
      forward_mismatch = "mismatch"
    )

  reverse_hits <- filtered_primer_blast_results |>
    dplyr::filter(stringr::str_detect(.data$qseqid, "reverse")) |>
    dplyr::select(
      subject_accession_version = "saccver",
      subject_gi = "sgi",
      unique_subject_taxonomy_ids = "staxids",
      reverse_start = "sstart",
      reverse_stop = "send",
      reverse_mismatch = "mismatch"
    )

  # TODO: note that the starts and stops may be a little different that you might expect given the generated sequences depending on what is around the primers. (e.g., see sequence_11)

  # TODO: we need to process these target by target, rather than building the whole join at once to avoid the cartesian explosion ... some targets may have 10,000 hits
  #
  # moorer | biomix14 | big_nt_test__split -> sed -E 's/^ +//' YO | awk '$1 > 10' | sort -nr | head
  # 11930 2440392461
  # 10539 2705591226
  # 8999 2020514959
  # 7410 2734211036
  # 7132 2440392458
  # 6843 2020514973
  # 6569 2211227459
  # 6041 2955233968
  # 5999 2626227570
  # 5957 2955233952
  # Join and calculate amplicon length in one step
  amplicons <- .find_all_product_lengths(forward_hits, reverse_hits) |>
    dplyr::filter(
      # product length will be NA if none of the case_when conditions match
      !rlang::is_na(.data$product_length),
      dplyr::between(
        .data$product_length,
        minimum_length,
        maximum_length
      )
    ) |>
    # What if there is more than one plausible amplicon,
    # shouldn't we select the best?
    dplyr::distinct(.data$subject_accession_version, .keep_all = TRUE)

  checkmate::assert_names(
    names(amplicons),
    must.include = c(
      "subject_accession_version",
      "subject_gi",
      "unique_subject_taxonomy_ids",
      "forward_start",
      "forward_stop",
      "forward_mismatch",
      "reverse_start",
      "reverse_stop",
      "reverse_mismatch",
      "product_length"
    )
  )

  amplicons
}
