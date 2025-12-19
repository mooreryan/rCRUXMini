# TODO: rename this to expand_multi_tax_ids
#
expand_multi_tax_ids2 <- function(
  blastdbcmd,
  hits_with_multiple_taxonomy_ids,
  blast_db_paths
) {
  checkmate::assert_string(blastdbcmd, min.chars = 1)

  checkmate::assert_character(
    blast_db_paths,
    # At least one
    min.len = 1,
    # Non empty strings
    min.chars = 1,
    all.missing = FALSE
  )

  checkmate::assert_data_frame(hits_with_multiple_taxonomy_ids)
  # TODO: a function below expects more columns I think....
  checkmate::assert_names(
    names(hits_with_multiple_taxonomy_ids),
    must.include = c("subject_accession_version")
  )

  # Get the subject accessions. We will pull these from the BLAST DB with
  # blastdbcmd.
  subject_accession_versions <-
    hits_with_multiple_taxonomy_ids$subject_accession_version

  checkmate::assert_character(subject_accession_versions)

  # If there aren't any, exit early, returning the original hit.
  if (length(subject_accession_versions) == 0) {
    log_debug("subject_accession_versions had zero entries")
    return(NULL)
  }

  # Write the entry batch file
  entry_batch_path <- tempfile()
  tibble::tibble(subject_accession_versions) |>
    readr::write_delim(file = entry_batch_path, delim = " ", col_names = FALSE)

  ############################################################################

  # Run the blastdbcmd jobs in parallel

  # These are dependencies of the run_blastdbcmd itself
  dependencies <- list(result__ok = result__ok, result__error = result__error)
  result <- furrr::future_map(
    .x = blast_db_paths,
    .options = furrr::furrr_options(globals = names(dependencies)),
    .env_globals = list2env(dependencies),
    .f = function(
      blastdb_path,
      run_blastdbcmd,
      blastdbcmd,
      entry_batch_path,
      outfmt
    ) {
      # NOTE: We aren't wrapping this in try_fetch, because it should not fail
      # under normal usage. If it does crash, then it is a bug.
      run_blastdbcmd(
        blastdbcmd = blastdbcmd,
        entry_batch_path = entry_batch_path,
        blast_db_path = blastdb_path,
        outfmt = outfmt
      )
    },
    run_blastdbcmd = run_blastdbcmd,
    blastdbcmd = blastdbcmd,
    entry_batch_path = entry_batch_path,
    outfmt = "%o ~~~ %a ~~~ %h ~~~ %T"
  )

  ############################################################################

  # checkmate list result
  #

  # Each of the jobs can either be ok or error. If any of them were errors, then
  # we quit early.
  partitioned_result <- result__partition(result)
  if (length(partitioned_result$error) > 0) {
    log_error("Oh no! Some of the blastdbcmds failed!")
    # TODO: Not a very user friendly string repr here of the error, but its okay
    # for now.
    cat(str(partitioned_result$error), file = stderr())
    abort_rcrux_mini_error(
      "blastdbcmd in expand_multi_tax_ids2 failed at least some of the time"
    )
  }

  # Column spec for the blastdbcmd output
  sequence_hash_column_specification <- readr::cols(
    # You will need this for grouping sequences later. It's an int but let's
    # treat it as a character.
    blast_ordinal_id = readr::col_character(),
    # This name matches the blast hit name
    subject_accession_version = readr::col_character(),
    sequence_hash_value = readr::col_character(),
    # This name matches the blast hit name
    unique_subject_taxonomy_ids = readr::col_character()
  )

  blastdbcmd_data <- partitioned_result$ok |>
    purrr::map(.f = function(result) {
      ok_data <- result__unwrap(
        result,
        # You shouldn't actually hit this default value, since we partitioned
        # the results into ok/error before, and we are only iterating on the OK
        # values.
        #
        # However, the stdout itself *could* still be empty!
        default = list(
          blast_db_path = "NA",
          stdout = ""
        )
      )

      blastdbcmd_data <- ok_data$stdout |>
        # TODO: try_catch with the reader error?
        #
        # NOTE: some of these stdouts may be the empty string, so need the I()
        # to treat it as data.
        I() |>
        readr::read_delim(
          col_names = names(sequence_hash_column_specification$cols),
          col_types = sequence_hash_column_specification,
          delim = " ~~~ "
        ) |>
        tibble::rowid_to_column("index") |>
        dplyr::mutate(blast_db_path = ok_data$blast_db_path)

      checkmate::assert_data_frame(blastdbcmd_data)
      checkmate::assert_names(
        names(blastdbcmd_data),
        must.include = c(
          "index",
          "blast_ordinal_id",
          "subject_accession_version",
          "sequence_hash_value",
          "unique_subject_taxonomy_ids",
          "blast_db_path"
        )
      )

      blastdbcmd_data <- .add_reference_subject_accession_versions(
        blastdbcmd_data
      )

      blastdbcmd_data
    }) |>
    purrr::list_rbind()

  checkmate::assert_data_frame(blastdbcmd_data)
  checkmate::assert_names(
    names(blastdbcmd_data),
    must.include = c(
      "index",
      "blast_ordinal_id",
      "subject_accession_version",
      "sequence_hash_value",
      "unique_subject_taxonomy_ids",
      "blast_db_path",
      "reference_subject_accession_version"
    )
  )

  result <- .expand_hits(
    blastdbcmd_result = blastdbcmd_data,
    blast_result = hits_with_multiple_taxonomy_ids
  )

  # TODO: this is obscure, but I need this check here too because I can't test
  # the full pipeline without mocking blast itself, because we can't construct a
  # blastdb that gives multi taxIDs...This function is tested in isolation, but
  # it needs the shape of the data to be correct.
  checkmate::assert_names(
    names(result),
    identical.to = c(
      "query_accession",
      "subject_accession_version",
      "percent_identical_matches",
      "alignment_length",
      "expect_value",
      "subject_sequence_length",
      "subject_alignment_start",
      "subject_alignment_end",
      "subject_aligned_sequence",
      "unique_subject_taxonomy_ids",
      # TODO: these degapped stuff must get tacked on at a later time
      "degapped_subject_aligned_sequence",
      "degapped_alignment_length"
    )
  )

  result
}

.add_reference_subject_accession_versions <- function(.data_frame) {
  # If the input has no rows, add the column that would have been added and
  # return.
  if (nrow(.data_frame) == 0) {
    return(dplyr::mutate(
      .data_frame,
      reference_subject_accession_version = character(0)
    ))
  }

  checkmate::assert_names(
    names(.data_frame),
    must.include = c(
      "subject_accession_version",
      "index",
      "blast_ordinal_id",
      "sequence_hash_value"
    )
  )

  .data_frame |>
    # Group by both to be sure its the same sequence cluster
    dplyr::group_by(.data$blast_ordinal_id, .data$sequence_hash_value) |>
    # Within a group/cluster, the lowest index 'should' be the
    # reference/centroid, i.e., the one's on which we actaully wanted to expand.
    #
    # This *should* hold as long as the blastdb cmd doesn't change out it
    # outputs data.
    #
    dplyr::mutate(
      reference_subject_accession_version = dplyr::first(
        .data$subject_accession_version,
        order_by = .data$index
      )
    ) |>
    dplyr::ungroup()
}

#' Ensure you have run \code{.add_reference_subject_accession_versions} on the
#' \code{blastdbcmd_result} before running this function. If you don't, you will
#' hit an assertion error.
#'
#'
#' NOTE: that we are "pretending" the original reference accessions only have a
#' single tax id after this function runs. That way things aren't counted more
#' than they should be, since they are now "expanded". TODO: need to double
#' check if the original rCRUX did this or not.
.expand_hits <- function(blastdbcmd_result, blast_result) {
  # TODO: given the way the multidb works, there could be duplicates
  # parent/representative IDs in the blastdbcmd results, and potentially in the
  # blast_result as well.

  checkmate::assert_data_frame(
    blastdbcmd_result,
    min.rows = 1,
    any.missing = FALSE
  )
  # TODO: shouldn't this also must include this: reference_subject_accession_version? If that original note above about adding in the reference subject accession versions is correct....
  #
  checkmate::assert_names(
    names(blastdbcmd_result),
    must.include = c(
      "index",
      "blast_ordinal_id",
      "subject_accession_version",
      "sequence_hash_value",
      "unique_subject_taxonomy_ids",
      "blast_db_path"
    )
  )

  checkmate::assert_data_frame(blast_result, min.rows = 1, any.missing = FALSE)
  checkmate::assert_names(
    names(blast_result),
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

  # We want the output columns of this function return value to match the original
  output_column_name_order <- names(blast_result)

  result <- dplyr::left_join(
    x = blastdbcmd_result,
    y = blast_result,
    by = dplyr::join_by(
      x$reference_subject_accession_version == y$subject_accession_version
    ),
    # Return every match detected in `y`
    multiple = "all",
    # Check for unmatched keys in `y` that could potentially drop rows. If any
    # keys are present in `x`, that are not present in `y`, that is a bug, and
    # we should signal an error.
    unmatched = "error",
    # Each row in `x` may have many rows in `y`, and vice versa.
    relationship = "many-to-many"
  ) |>
    dplyr::select(
      !c(
        "unique_subject_taxonomy_ids.y",
        "index",
        "blast_ordinal_id",
        "sequence_hash_value",
        "blast_db_path",
        "reference_subject_accession_version"
      )
    ) |>
    dplyr::rename(
      unique_subject_taxonomy_ids = "unique_subject_taxonomy_ids.x"
    ) |>
    # Put the columns back in the original order that BLAST returned
    dplyr::relocate(dplyr::all_of(output_column_name_order))

  checkmate::assert_data_frame(
    result,
    # Every row in the original blastdbcmd result should be present in the
    # output, but also, all rows in the original blast result should be either
    # present, or expanded. So we should never see fewer rows in the result than
    # either of the inputs.
    min.rows = max(nrow(blastdbcmd_result), nrow(blast_result)),
    # Missing values are not allowed AT ALL. If there are missing values, then
    # we have made a mistake in our assumptions, and it is a programmer error
    # that needs to be fixed.
    any.missing = FALSE
  )

  # The exact ordering of names in the blast result should be used in the joined
  # result here
  checkmate::assert_names(
    names(result),
    identical.to = output_column_name_order
  )

  result
}
