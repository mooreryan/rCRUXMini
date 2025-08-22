# plausible_amplicons_with_taxonomy = blast_seeds which comes from reading the
# seeds_output_path in the blast_seeds original function.
search_amplicons <- function(
  plausible_amplicons_with_taxonomy,
  blast_db_path,
  taxonomy_db_path
) {
  # TODO: switch to the assert data frame probably
  checkmate::assert_names(
    names(plausible_amplicons_with_taxonomy),
    must.include = c(
      # These columns come from `find_plausible_amplicon_coordinates`
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

  # This is for constructing the query file
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
