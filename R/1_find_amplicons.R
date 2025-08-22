find_amplicons <- function(
  forward_primers,
  reverse_primers,
  output_directory_path,
  blast_db_paths,
  taxonomy_db_path,
  random_seed,
  ncbi_bin_directory = NULL
) {
  # TODO: checkmate arguments
  print(forward_primers)
  print(reverse_primers)
  print(output_directory_path)
  print(blast_db_paths)
  print(taxonomy_db_path)
  print(random_seed)

  set.seed(random_seed)

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

  primer_blast_output_path <- file.path(
    output_directory_path,
    "primer_blast.tsv"
  )

  write(
    x = to_fasta_string("forward", enumerated_forward_primers),
    file = primers_fasta_path,
    append = FALSE
  )
  write(
    x = to_fasta_string("reverse", enumerated_reverse_primers),
    file = primers_fasta_path,
    append = TRUE
  )

  # TODO: Should this can throw, so it should be run in a rlang::try_fetch
  primer_blast_data <- run_primer_blast(
    blastn_executable = sys_which("blastn", directory = ncbi_bin_directory),
    query_path = primers_fasta_path,
    db_paths = blast_db_paths,
    output_path = primer_blast_output_path
  )

  primer_blast_results <- parse_primer_blast_results(
    primer_blast_data,
    maximum_mismatches = 4
  )

  # TODO: take these options from the config
  plausible_amplicons_coordinates <- find_plausible_amplicon_coordinates(
    primer_blast_results,
    minimum_length = 150,
    maximum_length = 650,
    maximum_mismatches = 4
  )

  # TODO: gracefully handle if there are no plausible amplicons.
  #       This throws the assertion error, which isn't very nice for the user.
  checkmate::assert_data_frame(plausible_amplicons_coordinates, min.rows = 1)

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

  result <- list(
    primer_blast_results = primer_blast_results,
    plausible_amplicons_coordinates_with_taxonomy = plausible_amplicons_coordinates_with_taxonomy,
    plausible_amplicons_coordinates_distinct_taxonomic_ranks = plausible_amplicons_coordinates_distinct_taxonomic_ranks
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

  # TODO: need to validate all these values
}

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

# NOTE: this returns the location of the amplicon on the target sequence, NOT
#       the actual sequence itself
#
# TODO: this should be differently named, like find_plausible_amplicon_coordinates or
#       something like that
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

  # Join and calculate amplicon length in one step
  amplicons <- dplyr::inner_join(
    forward_hits,
    reverse_hits,
    by = c("accession", "gi", "staxids"),
    relationship = "many-to-many"
  ) |>
    dplyr::mutate(
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
