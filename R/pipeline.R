# TODO: the names of the output files need fixing
# TODO: I don't think I actually need all the output files anymore
#
# START HERE: take the rest of the params from the config

#' Run the rCRUXMini pipeline
#'
#' @param forward_primers A character vector of forward primers.
#' @param reverse_primers A character vector of reverse primers.
#' @param output_directory_path A string representing the output directory path.
#' @param blast_db_paths A character vector of BLAST database paths.
#' @param taxonomy_db_path A string representing the taxonomy database path.
#' @param query_chunk_count An integer representing the query chunk count.
#' @param ncbi_bin_directory A string representing the NCBI binary directory path.
#' @return A list containing the pipeline results.
pipeline <- function(
  forward_primers,
  reverse_primers,
  output_directory_path,
  blast_db_paths,
  taxonomy_db_path,
  config,
  query_chunk_count = 1,
  ncbi_bin_directory = NULL
) {
  checkmate::assert_character(forward_primers, min.len = 1, min.chars = 1)
  checkmate::assert_character(reverse_primers, min.len = 1, min.chars = 1)
  checkmate::assert_string(output_directory_path, min.chars = 1)
  checkmate::assert_character(blast_db_paths, min.len = 1, min.chars = 1)
  checkmate::assert_string(taxonomy_db_path, min.chars = 1)
  checkmate::assert_count(query_chunk_count, positive = TRUE)
  checkmate::assert_string(ncbi_bin_directory, min.chars = 1, null.ok = TRUE)

  # From here on, we assume that the params was created properly
  # and has the required keys/names.
  assert_config_class(config)

  if (!dir.exists(output_directory_path)) {
    dir.create(
      output_directory_path,
      showWarnings = FALSE,
      mode = "0750",
      recursive = TRUE
    )
  }

  enumerated_forward_primers <- enumerate_ambiguities(forward_primers)
  enumerated_reverse_primers <- enumerate_ambiguities(reverse_primers)

  primers_fasta_path <- file.path(
    output_directory_path,
    "primers.fasta"
  )

  .write_primers(
    forward = enumerated_forward_primers,
    reverse = enumerated_reverse_primers,
    to_file = primers_fasta_path
  )

  primer_blast_data <- .run_blastn(
    blast_executable_directory = ncbi_bin_directory,
    query_paths = primers_fasta_path,
    db_paths = blast_db_paths,
    outfmt_specifiers = "qseqid sgi saccver mismatch sstart send staxids",
    extra_blast_arguments = primer_blast_config_to_cli_args(
      config$primer_blast
    ),
    use_long_names_in_parsed_result = FALSE
  )

  # TODO: checkmate the names of the output

  primer_blast_data |>
    readr::write_tsv(file.path(
      output_directory_path,
      "primer_blast.tsv"
    ))

  # TODO: shutdown gracefully
  if (nrow(primer_blast_data) == 0) {
    abort_rcrux_mini_error("No hits found in the primer blast data")
  }

  primer_blast_results <- parse_primer_blast_results(
    primer_blast_data,
    maximum_mismatches = 4
  )

  plausible_amplicons_coordinates <- find_plausible_amplicon_coordinates(
    primer_blast_results,
    minimum_length = config$plausible_amplicons$minimum_length,
    maximum_length = config$plausible_amplicons$maximum_length,
    maximum_mismatches = config$plausible_amplicons$maximum_mismatches
  )

  # TODO: gracefully handle if there are no plausible amplicons.
  #       This throws the assertion error, which isn't very nice for the user.
  checkmate::assert_data_frame(plausible_amplicons_coordinates)
  if (nrow(plausible_amplicons_coordinates) == 0) {
    abort_rcrux_mini_error("No plausible amplicons found")
  }

  # TODO: this data frame has a completely different naming scheme than all the ones that I made later in time. These are a lot of the old rCRUX names and they no longer match up with the other stuff
  plausible_amplicons_coordinates_with_taxonomy <- accession_to_taxonomy(
    plausible_amplicons_coordinates,
    taxonomy_db_path
  )

  # We shouldn't lose any rows
  checkmate::assert_true(
    nrow(plausible_amplicons_coordinates) ==
      nrow(plausible_amplicons_coordinates_with_taxonomy)
  )

  plausible_amplicons_coordinates_distinct_taxonomic_ranks <- distinct_taxonomic_ranks(
    plausible_amplicons_coordinates_with_taxonomy
  )

  blastdbcmd <- sys_which_or_abort(
    "blastdbcmd",
    ncbi_bin_directory = ncbi_bin_directory
  )

  amplicon_data <- pull_amplicons(
    blastdbcmd = blastdbcmd,
    amplicon_coordinates = plausible_amplicons_coordinates_with_taxonomy,
    chunks = query_chunk_count,
    blast_db_paths = blast_db_paths,
    ncbi_bin_directory = ncbi_bin_directory
  )

  if (nrow(amplicon_data) == 0) {
    abort_rcrux_mini_error(
      "no amplicon data TODO i thought this was already checked!"
    )
  }

  amplicon_data_tsv_path <- file.path(
    output_directory_path,
    "amplicon_data.tsv"
  )
  amplicon_data |> readr::write_tsv(file = amplicon_data_tsv_path)

  checkmate::assert_names(
    names(amplicon_data),
    # TODO: would be nice to always use accession version for everything
    must.include = c("accession", "sequence")
  )

  amplicon_query_fasta_paths <- .write_pieces_to_tempfastas(
    amplicon_data,
    chunks = query_chunk_count,
    id_column = "accession",
    sequence_column = "sequence"
  )

  amplicon_blast_result <- .run_blastn(
    blast_executable_directory = ncbi_bin_directory,
    query_paths = amplicon_query_fasta_paths,
    db_paths = blast_db_paths,
    # These are the outfmt specifiers from the original rCRUX
    #
    # TODO: why these? why not the usual + taxids?
    outfmt_specifiers = "qacc saccver pident length evalue slen sstart send sseq staxids",
    extra_blast_arguments = c("-num_threads", "1"),
    use_long_names_in_parsed_result = TRUE
  )

  # TODO: include the blast ordinal ID in the blast results, it will make it
  # 1000x easier to do the expand multi taxids step ... watch out though there
  # is another todo with a note about whether this would actually make it easier
  # or not.
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

  if (nrow(amplicon_blast_result) == 0) {
    abort_rcrux_mini_error("No hits found in the amplicon blast result")
  }

  # TODO take these from the config
  parsed_amplicon_blast_result <- parse_amplicon_blast_results(
    amplicon_blast_result = amplicon_blast_result,
    blastdbcmd = blastdbcmd,
    blast_db_paths = blast_db_paths,
    taxonomy_db_path = taxonomy_db_path
  ) |>
    checkmate::assert_data_frame()

  parsed_amplicon_blast_result <- parsed_amplicon_blast_result |>
    # TODO: these should probably be set in the parse function above
    dplyr::mutate(
      problem_missing_taxonomy = .check_for_missing_taxonomy(.data$taxonomy_id),
      problem_long_ambiguous_runs = .check_for_long_ambiguous_runs(
        aligned_sequences = .data$degapped_subject_aligned_sequence,
        ambiguous_run_limit = config$plausible_amplicons$ambiguous_run_limit
      ),
      # TODO: this should probably be its own function
      problem_bad_degapped_alignment_length = purrr::map_lgl(
        .x = .data$degapped_alignment_length,
        .f = function(len) {
          len < config$plausible_amplicons$minimum_length ||
            len > config$plausible_amplicons$maximum_length
        }
      ),
      # Same with this one...
      has_any_problem = .data$problem_missing_taxonomy |
        .data$problem_long_ambiguous_runs |
        .data$problem_bad_degapped_alignment_length
    )

  if (nrow(parsed_amplicon_blast_result) == 0) {
    # TODO: informative error message
    abort_rcrux_mini_error(
      "No hits remained after parsing the amplicon blast"
    )
  }

  # TODO: it would be nice to have the amplicon location on the original query as well
  # TODO: need the DB path in this data as well

  # NOTE: At this point you can have duplicate query accessions, duplicate
  # target accessions, etc. This is because we are allowing searches against
  # multiple databases, and those DBs may have overlaps, which means that a
  # subject sequence may be present in multiple DBs. Additionally, there is a
  # natural ability of queries to hit multiple targets, and targets to be hit by
  # multiple queries. Also tricky is that the subject sequences can have hits in
  # multiple locations along the subject sequence. This is okay! However, it
  # will make the counts different from the original rCRUX.
  #
  #

  # This is sort of like the old summary.csv file combined with some of the
  # other diagnostic files that reported sequences/amplicons that had any
  # problems.
  parsed_amplicon_blast_result |>
    readr::write_tsv(file.path(
      output_directory_path,
      "parsed_amplicon_blast_results.tsv"
    ))

  parsed_amplicon_blast_result |>
    write_fasta(
      to_filename = file.path(
        output_directory_path,
        "parsed_amplicon_blast_results.fasta"
      ),
      # Note that the _subject_ sequence represents the amplicon here, not the
      # query
      id_column = "subject_accession_version",
      sequence_column = "degapped_subject_aligned_sequence"
    )

  #  Count distinct taxonomic ranks again on the taxonomy info (rank counts
  # -- this is also done in the primer blasting step) (Need to ask Seanie about
  # this one: Since subjects can be present multiple times, their taxonomy lines
  # appear multiple times, so I probably want to count the distinct taxanomic
  # groups after taking unique subject+taxonomy pairs. The reason the taxonomy
  # needs to be in the identity calculation is that the same subject ID could
  # theoretically be present in multiple different DBs, and those DBs could have
  # a different taxonomy for each of them.)
  parsed_amplicon_blast_result_distinct_taxonomic_ranks <- distinct_taxonomic_ranks2(
    parsed_amplicon_blast_result,
    id_column = "subject_accession_version"
  )

  parsed_amplicon_blast_result_taxonomy <- parsed_amplicon_blast_result |>
    dplyr::select(c(
      "subject_accession_version",
      "superkingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    )) |>
    # See the previous notes for why we do distinct here...(multi DB...)
    dplyr::distinct() |>
    readr::write_tsv(file.path(
      output_directory_path,
      "parsed_amplicon_blast_result_taxonomy.tsv"
    ))

  result <- list(
    primer_blast_results = primer_blast_results,
    plausible_amplicons_coordinates_with_taxonomy = plausible_amplicons_coordinates_with_taxonomy,
    plausible_amplicons_coordinates_distinct_taxonomic_ranks = plausible_amplicons_coordinates_distinct_taxonomic_ranks,
    amplicon_blast_result = amplicon_blast_result,
    parsed_amplicon_blast_result = parsed_amplicon_blast_result,
    parsed_amplicon_blast_result_distinct_taxonomic_ranks = parsed_amplicon_blast_result_distinct_taxonomic_ranks,
    parsed_amplicon_blast_result_taxonomy = parsed_amplicon_blast_result_taxonomy
  )

  saveRDS(result, file.path(output_directory_path, "pipeline_results.rds"))

  readr::write_tsv(
    x = primer_blast_results,
    file = file.path(output_directory_path, "primer_blast_results.tsv")
  )
  readr::write_tsv(
    x = plausible_amplicons_coordinates_with_taxonomy,
    file = file.path(
      output_directory_path,
      "plausible_amplicons_coordinates_with_taxonomy.tsv"
    )
  )
  readr::write_tsv(
    x = plausible_amplicons_coordinates_distinct_taxonomic_ranks,
    file = file.path(
      output_directory_path,
      "plausible_amplicons_coordinates_distinct_taxonomic_ranks.tsv"
    )
  )

  parsed_amplicon_blast_result_distinct_taxonomic_ranks |>
    readr::write_tsv(file.path(
      output_directory_path,
      "parsed_amplicon_blast_result_distinct_taxonomic_ranks.tsv"
    ))

  result
}

# NOTE: You may have duplicate accessions in this output, IF they have different
# coordinates.
#
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
    # TODO: you need the other function to check col names
    # col.names = names(primer_blast_column_types),
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

  # TODO: column names assert

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
      # TODO: pull these forward and reverse names into a global -- it's also
      #       used in to_fasta_string
      any(grepl("forward", .data$qseqid)) & any(grepl("reverse", .data$qseqid))
    ) |>
    # Remove duplicate end positions per accession
    dplyr::distinct(.data$saccver, .data$send, .keep_all = TRUE)

  # Ungroup and sort
  result <- result |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$saccver, .data$sstart, .data$mismatch)

  # TODO: min rows shouldn't be a checkmate thing right? or is it just to ensure
  #       we didn't lose anything in the parsing

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
      accession = "saccver",
      gi = "sgi",
      "staxids",
      forward_start = "sstart",
      forward_stop = "send",
      forward_mismatch = "mismatch"
    )

  reverse_hits <- filtered_primer_blast_results |>
    dplyr::filter(stringr::str_detect(.data$qseqid, "reverse")) |>
    dplyr::select(
      accession = "saccver",
      gi = "sgi",
      "staxids",
      reverse_start = "sstart",
      reverse_stop = "send",
      reverse_mismatch = "mismatch"
    )

  # TODO: note that the starts and stops may be a little different that you might expect given the generated sequences depending on what is around the primers. (e.g., see sequence_11)

  # Join and calculate amplicon length in one step
  amplicons <- dplyr::inner_join(
    forward_hits,
    reverse_hits,
    by = c("accession", "gi", "staxids"),
    relationship = "many-to-many"
  ) |>
    dplyr::mutate(
      # TODO: this product length is off by 2 (off by one on both sides?)
      product_length = dplyr::case_when(
        #  F ---------->                        <------------ R  Primers
        # ====================================================== DNA Target
        #    ---------->************************<------------    Amplicon
        #    ^ Forward start                  Reverse start ^
        (forward_start < reverse_start &
          forward_start < forward_stop &
          reverse_stop < reverse_start) ~
          (reverse_start - forward_start),
        #  R ---------->                        <------------ F  Primers
        # ====================================================== DNA Target
        #    ---------->************************<------------    Amplicon
        #    ^ Reverse start                  Forward start ^
        #
        (forward_start > reverse_start &
          forward_start > forward_stop &
          reverse_stop > reverse_start) ~
          (forward_start - reverse_start),
      )
    ) |>
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
    dplyr::distinct(.data$accession, .keep_all = TRUE)

  # TODO: switch to assert data frame
  checkmate::assert_names(
    names(amplicons),
    must.include = c(
      "accession",
      "gi",
      "staxids",
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

# This function will throw if taxonomy data cannot be found for the given
# accessions.
#
# TODO (important): do I really want this to throw? Probably not. Rather, it
# should report accessions that didn't have taxonomy, and return two things, one
# with good accessions with taxonomy and one with missing accessions.
accession_to_taxonomy <- function(input_data, taxonomy_db_path) {
  checkmate::assert_data_frame(input_data, min.rows = 1)
  checkmate::assert_names(names(input_data), must.include = c("accession"))
  checkmate::assert_character(input_data$accession)
  checkmate::assert_file_exists(taxonomy_db_path)

  # TODO: What does this do for accessions that don't have taxonomy IDs?
  taxonomy_ids <- taxonomizr::accessionToTaxa(
    input_data$accession,
    taxonomy_db_path
  )
  checkmate::assert_numeric(taxonomy_ids, len = length(input_data$accession))

  # TODO: What does this do for taxonomy IDs that do not exist?
  taxonomy_data <- taxonomizr::getTaxonomy(
    taxonomy_ids,
    taxonomy_db_path,
    desiredTaxa = c(
      "species",
      "superkingdom",
      "kingdom",
      "phylum",
      "subphylum",
      "superclass",
      "class",
      "subclass",
      "order",
      "family",
      "subfamily",
      "genus",
      "infraorder",
      "subcohort",
      "superorder",
      "superfamily",
      "tribe",
      "subspecies",
      "subgenus",
      "species group",
      "parvorder",
      "varietas"
    )
  )
  checkmate::assert_matrix(taxonomy_data, nrows = length(input_data$accession))

  result <- input_data |>
    dplyr::mutate(taxonomy_id = taxonomy_ids, data.frame(taxonomy_data)) |>
    dplyr::arrange(
      "superkingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    )

  checkmate::assert_data_frame(result, nrows = length(input_data$accession))

  result
}

accession_to_taxonomy2 <- function(
  accessions,
  taxonomy_db_path,
  accession_column_name
) {
  checkmate::assert_character(accessions)
  checkmate::assert_string(accession_column_name)
  checkmate::assert_file_exists(taxonomy_db_path)

  if (length(accessions) == 0) {
    return(tibble::tibble(
      species = character(0),
      superkingdom = character(0),
      kingdom = character(0),
      phylum = character(0),
      subphylum = character(0),
      superclass = character(0),
      class = character(0),
      subclass = character(0),
      order = character(0),
      family = character(0),
      subfamily = character(0),
      genus = character(0),
      infraorder = character(0),
      subcohort = character(0),
      superorder = character(0),
      superfamily = character(0),
      tribe = character(0),
      subspecies = character(0),
      subgenus = character(0),
      `species group` = character(0),
      parvorder = character(0),
      varietas = character(0)
    ))
  }

  # If any of the accessions do NOT have a match in the taxonomy DB, then those
  # will have `NA` in the result.
  taxonomy_ids <- taxonomizr::accessionToTaxa(accessions, taxonomy_db_path)

  taxonomy_df <- tibble::tibble(
    # This will splice in the string held by the accession_column_name parameter
    # as a column name
    !!accession_column_name := accessions,
    taxonomy_id = taxonomy_ids
  )

  # TODO: I don't actually want this split up.

  tmp <- .split_by_has_taxonomy_id(taxonomy_df)
  taxonomy_df <- tmp$has_taxonomy_id
  missing_taxonomy_df <- tmp$missing_taxonomy_id
  # TODO: log the accessions with missing taxonomy IDs

  # TODO: now I need to split out rows that have NA for taxonomy ID.

  # For taxonomy IDs that don't exist, all elements of that row will be NA.
  # This will only contain columns that are specified in the desiredTaxa arg.
  desired_taxa <- c(
    "species",
    "superkingdom",
    "kingdom",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "order",
    "family",
    "subfamily",
    "genus",
    "infraorder",
    "subcohort",
    "superorder",
    "superfamily",
    "tribe",
    "subspecies",
    "subgenus",
    "species group",
    "parvorder",
    "varietas"
  )
  taxonomy_details <- taxonomizr::getTaxonomy(
    taxonomy_df$taxonomy_id,
    taxonomy_db_path,
    desiredTaxa = desired_taxa
  ) |>
    tibble::as_tibble()

  checkmate::assert_data_frame(
    taxonomy_details,
    nrows = length(taxonomy_df$taxonomy_id)
  )
  checkmate::assert_names(names(taxonomy_details), identical.to = desired_taxa)

  result <- dplyr::bind_cols(taxonomy_df, taxonomy_details) |>
    dplyr::arrange(
      "superkingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    )

  list(has_taxonomy_id = result, missing_taxonomy_id = missing_taxonomy_df)
}

accession_to_taxonomy3 <- function(
  accessions,
  taxonomy_db_path,
  accession_column_name
) {
  checkmate::assert_character(accessions)
  checkmate::assert_string(accession_column_name)
  checkmate::assert_file_exists(taxonomy_db_path)

  if (length(accessions) == 0) {
    # TODO: this def needs to be hit in the tests.
    return(tibble::tibble(
      !!accession_column_name := character(0),
      taxonomy_id = integer(0),
      species = character(0),
      superkingdom = character(0),
      kingdom = character(0),
      phylum = character(0),
      subphylum = character(0),
      superclass = character(0),
      class = character(0),
      subclass = character(0),
      order = character(0),
      family = character(0),
      subfamily = character(0),
      genus = character(0),
      infraorder = character(0),
      subcohort = character(0),
      superorder = character(0),
      superfamily = character(0),
      tribe = character(0),
      subspecies = character(0),
      subgenus = character(0),
      `species group` = character(0),
      parvorder = character(0),
      varietas = character(0)
    ))
  }

  # If any of the accessions do NOT have a match in the taxonomy DB, then those
  # will have `NA` in the result.
  taxonomy_ids <- taxonomizr::accessionToTaxa(accessions, taxonomy_db_path)

  taxonomy_df <- tibble::tibble(
    # This will splice in the string held by the accession_column_name parameter
    # as a column name
    !!accession_column_name := accessions,
    # Some of these may be NA
    taxonomy_id = taxonomy_ids
  )

  # For taxonomy IDs that don't exist, all elements of that row will be NA.
  # This will only contain columns that are specified in the desiredTaxa arg.
  desired_taxa <- c(
    "species",
    "superkingdom",
    "kingdom",
    "phylum",
    "subphylum",
    "superclass",
    "class",
    "subclass",
    "order",
    "family",
    "subfamily",
    "genus",
    "infraorder",
    "subcohort",
    "superorder",
    "superfamily",
    "tribe",
    "subspecies",
    "subgenus",
    "species group",
    "parvorder",
    "varietas"
  )
  taxonomy_details <- taxonomizr::getTaxonomy(
    taxonomy_df$taxonomy_id,
    taxonomy_db_path,
    desiredTaxa = desired_taxa
  ) |>
    tibble::as_tibble()

  checkmate::assert_data_frame(
    taxonomy_details,
    nrows = length(taxonomy_df$taxonomy_id)
  )
  checkmate::assert_names(names(taxonomy_details), identical.to = desired_taxa)

  result <- dplyr::bind_cols(taxonomy_df, taxonomy_details) |>
    dplyr::arrange(
      "superkingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    )

  checkmate::assert_names(
    names(result),
    # TODO: is permutation.of fine here? I think so...
    identical.to = c(
      accession_column_name,
      "taxonomy_id",
      desired_taxa
    )
  )

  result
}

distinct_taxonomic_ranks <- function(df) {
  taxonomy_levels <- c(
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
  checkmate::assert_data_frame(df, min.rows = 1)
  checkmate::assert_names(names(df), must.include = taxonomy_levels)

  # TODO: This is not how the LCA script treats taxonomic ranks
  result <- df |>
    dplyr::summarise(dplyr::across(
      dplyr::all_of(taxonomy_levels),
      .fns = dplyr::n_distinct
    ))

  checkmate::assert_names(names(result), must.include = taxonomy_levels)

  result
}

distinct_taxonomic_ranks2 <- function(df, id_column) {
  taxonomy_levels <- c(
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )

  identity_columns <- c(id_column, taxonomy_levels)

  checkmate::assert_data_frame(df, min.rows = 1)
  checkmate::assert_names(names(df), must.include = identity_columns)

  # TODO: This is not how the LCA script treats taxonomic ranks
  result <- df |>
    dplyr::summarise(dplyr::across(
      dplyr::all_of(identity_columns),
      .fns = dplyr::n_distinct
    ))

  checkmate::assert_names(names(result), must.include = identity_columns)

  result
}


# TODO: not used....
# NOTE: It is *possible* that you have duplicate accession-range pairs, if the
# user provided BLAST DBs that have some overlapping sequences. This means that
# you could in theory have exact duplicates in the resulting amplicon
# files/data. (TODO may need to address this -- not sure if it's really a big
# deal, because when you blast those amplicons, you would get "duplicate" hits
# too...so potentially as long as we dedup the amplicon blast results, we can
# ignore this). The question is, do I care to attach the range to the resulting
# IDs? I think we will just leave them off for now.
write_plausible_amplicon_sequences <- function(
  plausible_amplicons_coordinates_with_taxonomy,
  blastdbcmd_entry_batch_path,
  plausible_amplicon_path,
  ncbi_bin_directory = NULL
) {
  # Empty data frames shouldn't be passed in here, they should be handled before
  # you get to this function
  checkmate::assert_data_frame(
    plausible_amplicons_coordinates_with_taxonomy,
    min.rows = 1
  )
  checkmate::assert_names(
    names(plausible_amplicons_coordinates_with_taxonomy),
    must.include = c(
      "accession",
      "forward_start",
      "forward_stop",
      "reverse_start",
      "reverse_stop"
    )
  )

  # TODO: not sure about correct orientation, e.g., if the reverse primer is
  # actually provided in the reversecomplement orientation then the reverse
  # start will be < reverse stop, but if they provide it in the other direction,
  # reverse start > reverse stop, and we won't know which way the user will have
  # provided it. so as long as the forward is to the left of the reverse its
  # probably fine?

  plausible_amplicons_coordinates_with_taxonomy |>
    # # Take only the columns that we need to construct the entry_batch file
    # dplyr::select("accession", "forward_start", "reverse_stop") |>
    dplyr::mutate(
      far_left_coordinate = min(
        .data$forward_start,
        .data$forward_stop,
        .data$reverse_start,
        .data$reverse_stop
      ),
      far_right_coordinate = max(
        .data$forward_start,
        .data$forward_stop,
        .data$reverse_start,
        .data$reverse_stop
      )
    ) |>
    dplyr::select("accession", "far_left_coordinate", "far_right_coordinate") |>
    # Original code uses forward_stop and reverse_stop, which I think is a bug.
    # Since that would exclude the forward primer and include the reverse. This
    # way is including both.
    tidyr::unite(
      col = "range",
      "far_left_coordinate",
      "far_right_coordinate",
      sep = "-",
      remove = TRUE
    ) |>
    # This should already be all the columns that are left, but just to be
    # sure....
    dplyr::select("accession", "range") |>
    readr::write_delim(
      file = blastdbcmd_entry_batch_path,
      delim = " ",
      col_names = FALSE
    )

  # blastdbcmd -entry all -db nt -outfmt '%a | %s'
  #

  blastdbcmd <- SnailBLAST::sys_which("blastdbcmd", ncbi_bin_directory)
  # TODO: need a better user facing error message for this
  checkmate::assert_string(blastdbcmd, min.chars = 1)
  result <- processx::run(
    command = blastdbcmd,
    args = c(
      "-db",
      blast_db_path,
      "-dbtype",
      "nucl",
      "-entry_batch",
      blastdbcmd_entry_batch_path,
      "-outfmt",
      # This is sequence ID and sequence
      "%a %s"
    ),
    error_on_status = FALSE
  )

  if (rlang::is_na(result$status) || result$status != 0) {
    # TODO: custom error class and better error message
    stop("blastdbcmd to create amplicons failed")
  }

  result$stdout |>
    I() |>
    readr::read_delim(col_names = c("id", "sequence"), delim = " ")
}

write_plausible_amplicon_sequences2 <- function(
  plausible_amplicons_coordinates_with_taxonomy,
  blast_db_path,
  out_dir,
  prefix = "amplicons",
  num_splits = 1
) {
  checkmate::assert_data_frame(plausible_amplicons_coordinates_with_taxonomy)
  checkmate::assert_names(
    names(plausible_amplicons_coordinates_with_taxonomy),
    must.include = c("accession", "forward_start", "reverse_stop")
  )
  checkmate::assert_string(blast_db_path)
  checkmate::assert_directory_exists(out_dir)
  checkmate::assert_count(num_splits, positive = TRUE)

  # Prepare data for entry_batch: accession and start-end range
  df <- plausible_amplicons_coordinates_with_taxonomy |>
    dplyr::select("accession", "forward_start", "reverse_stop") |>
    tidyr::unite(
      col = "range",
      "forward_start",
      "reverse_stop",
      sep = "-",
      remove = TRUE
    ) |>
    dplyr::select("accession", "range")

  # TODO
  n <- nrow(df)
  if (n == 0) {
    stop("No amplicon coordinates provided (0 rows).")
  }

  k <- max(1, min(num_splits, n)) # don't create empty splits
  # Split indices into k roughly equal groups
  grp_id <- cut(seq_len(n), breaks = k, labels = FALSE)
  splits <- split(df, grp_id)

  results <- vector("list", length(splits))

  for (i in seq_along(splits)) {
    batch_path <- file.path(
      out_dir,
      sprintf("%s.entry_batch.%03d.txt", prefix, i)
    )
    query_path <- file.path(out_dir, sprintf("%s.query.%03d.fasta", prefix, i))

    # Write entry_batch with no header, space-delimited: "accession range"
    readr::write_delim(
      splits[[i]],
      file = batch_path,
      delim = " ",
      col_names = FALSE
    )

    # Run blastdbcmd for this batch
    # -outfmt %f ensures FASTA
    exit_code <- system2(
      command = "blastdbcmd",
      args = c(
        "-db",
        blast_db_path,
        "-dbtype",
        "nucl",
        "-entry_batch",
        batch_path,
        "-outfmt",
        "%f",
        "-out",
        query_path
      ),
      stdout = TRUE,
      stderr = TRUE
    )

    # You can enhance error handling by inspecting attr(exit_code, "status")
    # and the captured stdout/stderr
    status <- attr(exit_code, "status")
    if (!is.null(status) && status != 0) {
      stop(sprintf(
        "blastdbcmd failed for split %d (status %s). See stderr above.",
        i,
        status
      ))
    }

    results[[i]] <- list(
      split = i,
      entry_batch_path = batch_path,
      query_path = query_path
    )
  }

  tibble::tibble(
    split = vapply(results, `[[`, integer(1), "split"),
    entry_batch_path = vapply(results, `[[`, character(1), "entry_batch_path"),
    query_path = vapply(results, `[[`, character(1), "query_path")
  )
}


#' Checks if `path` is a file that exists and that it has a non-zero size.
#'
#' If `path` refers to a directory this will fail.
test_file_non_empty <- function(path) {
  checkmate::test_file_exists(path) && file.size(path) > 0
}


.write_primers <- function(forward, reverse, to_file) {
  write(
    x = to_fasta_string("forward", forward),
    file = to_file,
    append = FALSE
  )
  write(
    x = to_fasta_string("reverse", reverse),
    file = to_file,
    append = TRUE
  )
}


# TODO: actual logging for failures and such
# TODO: would be nice to include the blast_db_path in the output here.
#
#' Wrapper for \code{SnailBLAST::crawl} to run \code{blastn} that specifies
#' standard rCRUX error handlers.
#'
.run_blastn <- function(
  blast_executable_directory,
  query_paths,
  db_paths,
  outfmt_specifiers,
  extra_blast_arguments,
  use_long_names_in_parsed_result
) {
  SnailBLAST::crawl(
    "blastn",
    blast_executable_directory = blast_executable_directory,
    query_paths = query_paths,
    db_paths = db_paths,
    # These are the original rCRUX specifiers TODO i think this needs to be adjusted
    outfmt_specifiers = outfmt_specifiers,
    extra_blast_arguments = extra_blast_arguments,
    use_long_names_in_parsed_result = use_long_names_in_parsed_result,
    job_failed_callback = function(
      query_path,
      db_path,
      exit_status,
      command,
      args,
      stderr
    ) {
      message(paste(
        "Job failed for query",
        query_path,
        "and database",
        db_path,
        "with exit status",
        exit_status,
        "and command",
        command,
        "and arguments",
        args,
        "and stderr",
        stderr
      ))
    },
    parse_failed_callback = function(
      query_path,
      db_path,
      error_condition,
      command,
      args,
      stderr
    ) {
      message(paste(
        "Parsing output for query",
        query_path,
        "and database",
        db_path,
        "with error condition",
        error_condition,
        "and command",
        command,
        "and arguments",
        args,
        "and stderr",
        stderr
      ))
    }
  )
}
