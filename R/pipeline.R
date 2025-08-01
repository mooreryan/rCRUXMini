# TODO: original rCRUX has mismatches 4 (as a constant) and 6 (as a parameter)
#       in get_seeds_local
# TODO: grepl to str_detect?
# TODO: I think I will need a run ID to make it easier to hook parts into the
#       multidb pipeline
# TODO: Use drop = FALSE for all data frame slicing
# TODO: Original get_seeds_local has an unfiltered output and a primers chosen

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

# Defaults are from the original rCRUX
run_primer_blast <- function(
  query_path,
  target_path,
  output_path,
  word_size = 7,
  evalue = "3e+07",
  num_alignments = "10000000",
  qcov_hsp_perc = 90,
  perc_identity = 50,
  reward = 2,
  num_threads = 1,
  ncbi_bin = NULL
) {
  if (rlang::is_null(ncbi_bin)) {
    blastn <- "blastn"
  } else {
    blastn <- file.path(ncbi_bin, "blastn")
  }

  args <-
    c(
      "-db",
      target_path,
      "-task",
      "blastn-short",
      "-query",
      query_path,
      "-out",
      output_path,
      "-outfmt",
      # staxids means unique Subject Taxonomy ID(s), separated by a ';' (in
      # numerical order) -- TODO: probably need to account for this, the
      # original code doesn't though.
      '"6 qseqid sgi saccver mismatch sstart send staxids"',
      "-evalue",
      evalue,
      "-num_alignments",
      num_alignments,
      "-qcov_hsp_perc",
      qcov_hsp_perc,
      "-perc_identity",
      perc_identity,
      "-reward",
      reward,
      "-word_size",
      word_size,
      "-num_threads",
      num_threads
    )

  run_command(command = blastn, args = args)
}


parse_primer_blast_results <- function(file_name, maximum_mismatches = 4) {
  file_data <- readr::read_tsv(
    file_name,
    col_names = c(
      "qseqid",
      "sgi",
      "saccver",
      "mismatch",
      "sstart",
      "send",
      "staxids"
    ),
    col_types = c("ccciiic")
  )

  if (nrow(file_data) < 1) {
    stop("NO PRIMER BLAST RESULTS")
  }

  result <- file_data |>
    # TODO: max mismatches should be a parameter
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

  result
}

find_plausible_amplicons <- function(
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

# Search the primers against the given database to identify plausible amplicons.
# These will be used later in the `blast_seeds` step, or will need to be
# collated if the user wants to run the search in parallel across multiple DBs.
# TODO: rename to something like "identify_plausible_amplicons" or something
#       like that.
get_seeds <- function(
  forward_primers,
  reverse_primers,
  output_directory_path,
  blast_db_path,
  taxonomy_db_path,
  random_seed
) {
  print(forward_primers)
  print(reverse_primers)
  print(output_directory_path)
  print(blast_db_path)
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

  # TODO: Should check this exit code
  run_primer_blast(
    query_path = primers_fasta_path,
    target_path = blast_db_path,
    output_path = primer_blast_output_path
  )

  primer_blast_results <- parse_primer_blast_results(
    primer_blast_output_path,
    maximum_mismatches = 4
  )

  # TODO: take these options from the config
  plausible_amplicons <- find_plausible_amplicons(
    primer_blast_results,
    minimum_length = 150,
    maximum_length = 650,
    maximum_mismatches = 4
  )

  # TODO: handle if there are no plausible amplicons

  plausible_amplicons_with_taxonomy <- accession_to_taxonomy(
    plausible_amplicons,
    taxonomy_db_path
  )

  # We shouldn't lose any rows
  checkmate::assert_true(
    nrow(plausible_amplicons) == nrow(plausible_amplicons_with_taxonomy)
  )

  plausible_amplicons_distinct_taxonomic_ranks <- distinct_taxonomic_ranks(
    plausible_amplicons_with_taxonomy
  )

  result <- list(
    primer_blast_results = primer_blast_results,
    plausible_amplicons_with_taxonomy = plausible_amplicons_with_taxonomy,
    plausible_amplicons_distinct_taxonomic_ranks = plausible_amplicons_distinct_taxonomic_ranks
  )

  saveRDS(result, file.path(output_directory_path, "pipeline_results.rds"))

  readr::write_tsv(
    x = primer_blast_results,
    file = file.path(output_directory_path, "primer_blast_results.tsv")
  )
  readr::write_tsv(
    x = plausible_amplicons_with_taxonomy,
    file = file.path(
      output_directory_path,
      "plausible_amplicons_with_taxonomy.tsv"
    )
  )
  readr::write_tsv(
    x = plausible_amplicons_distinct_taxonomic_ranks,
    file = file.path(
      output_directory_path,
      "plausible_amplicons_distinct_taxonomic_ranks.tsv"
    )
  )
}

#################################
# blast_seeds ###################
#################################

# TODO: come back to this after all the blasts are done.
expand_multi_taxids <- function(blastn_result) {
  checkmate::assert_names(
    names(blastn_result),
    must.include = c("accession", "BLAST_db_taxids")
  )

  multi_taxid_hits <- blastn_result |>
    dplyr::filter(stringr::str_detect(
      string = .data$BLAST_db_taxids,
      pattern = stringr::fixed(";")
    ))

  entry_batch_path <- "TODO"
  accession_taxid_path <- "TODO"

  multi_taxid_hits$accession |> readr::write_lines(file = entry_batch_path)

  makeblastdb_result <- run_command(
    command = "makeblastdb",
    args = c(
      "-db",
      blast_db_path,
      "-dbtype",
      "nucl",
      "-entry_batch",
      entry_batch_path,
      "-out",
      accession_taxid_path,
      "-outfmt",
      # accession[space]taxid
      "%a %T"
    )
  )
  # TODO: check exit code of the command
  #

  accession_taxids <- readr::read_delim(
    file = accession_taxid_path,
    delim = " ",
    col_names = c("accession", "BLAST_db_taxids"),
    col_types = list(
      "accession" = readr::col_character(),
      # TODO: elsewhere these might be treated as integers
      "BLAST_db_taxids" = readr::col_character()
    )
  )
  readr::stop_for_problems(accession_taxids)

  single_taxid_hits <- blastn_result |> dplyr::setdiff(multi_taxid_hits)
}

# plausible_amplicons_with_taxonomy = blast_seeds which comes from reading the
# seeds_output_path in the blast_seeds original function.
blast_seeds <- function(
  plausible_amplicons_with_taxonomy,
  blast_db_path,
  taxonomy_db_path
) {
  checkmate::assert_names(
    names(plausible_amplicons_with_taxonomy),
    must.include = c(
      # These columns come from `find_plausible_amplicons`
      "accession",
      "gi",
      "staxids",
      "forward_start",
      "forward_stop",
      "forward_mismatch",
      "reverse_start",
      "reverse_stop",
      "reverse_mismatch",
      "product_length",
      # These columns come from `accession_to_taxonomy`
      "superkingdom",
      "phylum",
      "class",
      "order",
      "family",
      "genus",
      "species"
    )
  )

  # TODO: check blast exe
  # TODO: check accession taxa
  #

  entry_batch_path <- "TODO"

  plausible_amplicons_with_taxonomy |>
    # Take only the columns that we need to construct the entry_batch file
    dplyr::select("accession", "forward_start", "reverse_stop") |>
    # Original code uses forward_stop and reverse_stop, which I think is a bug.
    # Since that would exclude the forward primer and include the reverse. This
    # way is including both.
    dplyr::unite(
      col = "range",
      "forward_start",
      "reverse_stop",
      remove = TRUE
    ) |>
    readr::write_delim(file = entry_batch_path, delim = " ")

  # blastdbcmd example:
  # $ blastdbcmd -entry K02672.1 -db /usr/local/blastdb_v5/nt -range 100-200
  # >K02672.1:100-200 Escherichia coli ribonucleoside diphosphate reductase operon, complete sequence
  # CAGCGTGTAGTCATACGCTCCCATATTGACGTAGCCGCTGCCTGCCAGACTGAACTGCGCATTGTTTTGAAACTGCGTGG
  # GATCAGCTGCGAAATCAACCA

  query_sequence_path <- "TODO"

  makeblastdb_result <- run_command(
    command = "makeblastdb",
    args = c(
      "-db",
      blast_db_path,
      "-dbtype",
      "nucl",
      "-entry_batch",
      entry_batch_path,
      "-out",
      query_sequence_path
    )
  )
  # TODO: check exit code of the command

  blastn_results_path <- "TODO"

  blastn_result <- run_command(
    command = "blastn",
    args = c(
      "-db",
      blast_db_path,
      "-query",
      query_sequence_path,
      "-out",
      blastn_results_path,
      "-outfmt",
      "'6 saccver length pident qacc slen sstart send sseq evalue staxids'",
      "-evalue",
      "1e-6",
      "-num_alignments",
      "TODO",
      "-qcov_hsp_perc",
      "TODO",
      "-perc_identity",
      "TODO",
      "-num_threads ",
      "TODO"
    )
  )
  # TODO: check exit code

  # TODO: these are the column names from the original. Rename them or just use
  #       default blast names.
  blastn_result <- readr::read_tsv(
    blastn_results_path,
    col_names = c(
      "accession",
      "amplicon_length",
      "pident",
      "query_accession",
      "accession_sequence_length",
      "amplicon_start",
      "amplicon_stop",
      "sequence",
      "evalue",
      "BLAST_db_taxids"
    ),
    col_types = c(
      "accession" = readr::col_character(),
      "amplicon_length" = readr::col_integer(),
      "pident" = readr::col_double(),
      "query_accession" = readr::col_character(),
      "accession_sequence_length" = readr::col_integer(),
      "amplicon_start" = readr::col_integer(),
      "amplicon_stop" = readr::col_integer(),
      "sequence" = readr::col_character(),
      "evalue" = readr::col_double(),
      "BLAST_db_taxids" = readr::col_character()
    )
  )
  readr::stop_for_problems(blastn_result)

  # Keep only one sequence per accession (the longest)
  blastn_result <- blastn_result |>
    dplyr::group_by(.data$accession) |> # nolint: object_usage_linter
    dplyr::slice_max(.data$amplicon_length, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()

  # TODO: why do they need this?
  # Remove hyphens from the sequence and recalculate the sequence length
  blastn_result <- blastn_result |>
    dplyr::mutate(
      sequence = stringr::str_remove_all(
        # TODO: does .data work here?
        string = .data$sequence,
        pattern = stringr::fixed("-")
      )
    ) |>
    dplyr::mutate(amplicon_length = stringr::str_length(sequence))

  # TODO: track the sequences that are removed due to too many wildcards in
  #       addition to only removing them
  blastn_result <- blastn_result |>
    dplyr::filter(stringr::str_detect(
      # TODO: does .data work here?
      string = .data$sequence,
      pattern = stringr::fixed(wildcards)
    ))

  # Get any accessions that have multiple taxonomy IDs
  multiple_taxids <- blastn_result |>
    dplyr::filter(stringr::str_detect(
      string = .data$BLAST_db_taxids,
      pattern = stringr::fixed(";")
    ))

  # TODO: handle case where there are none with multiple taxids
}
