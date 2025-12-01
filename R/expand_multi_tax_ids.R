#' If \code{hits_with_multiple_taxonomy_ids} has zero rows, this function will
#' simply return that value back. That way you can join it with the single hits
#' and nothing will break, avoiding more conditional expressions.
expand_multi_tax_ids <- function(
  # TODO: checkmate
  blastdbcmd,
  # TODO: need to double check, but will this break if you pass in hits with
  # only a single taxonomy ID?
  hits_with_multiple_taxonomy_ids,
  blast_db_paths
) {
  checkmate::assert_string(blastdbcmd, min.chars = 1)

  checkmate::assert_data_frame(hits_with_multiple_taxonomy_ids)
  checkmate::assert_names(
    names(hits_with_multiple_taxonomy_ids),
    must.include = c("subject_accession_version")
  )

  checkmate::assert_character(
    blast_db_paths,
    # At least one
    min.len = 1,
    # Non empty strings
    min.chars = 1,
    all.missing = FALSE
  )

  if (nrow(hits_with_multiple_taxonomy_ids) == 0) {
    return(hits_with_multiple_taxonomy_ids)
  }

  # This is a function like the one in the original TODO it needs a better name though
  #
  #
  # after you run the blast db cmd, you need to then match up the IDs and stuff back to the original accession that has the sequence info, and take all that info and attatch it to these IDs and stuff

  entry_batch_path <- tempfile()
  hits_with_multiple_taxonomy_ids |>
    dplyr::select("subject_accession_version") |>
    readr::write_delim(file = entry_batch_path, delim = " ", col_names = FALSE)

  # TODO: almost all of this function can be replaced by pull_amplicons, it just needs renaming

  result <- furrr::future_map(
    .x = blast_db_paths,
    .pull_sequence_hashes_from_blast_db = .pull_sequence_hashes_from_blast_db,
    empty_tibble_from_col_spec = empty_tibble_from_col_spec,
    .f = function(blast_db_path) {
      rlang::try_fetch(
        .pull_sequence_hashes_from_blast_db(
          blastdbcmd = blastdbcmd,
          entry_batch_path = entry_batch_path,
          blast_db_path = blast_db_path
        ),
        # TODO: need to check for the rcrux_error instead
        error = function(condition) {
          # TODO: log this error
          #
          rlang::warn("something bad happend in .pull_sequences_from_blast_db")

          empty_tibble_from_col_spec(.sequence_hash_column_specification) |>
            tibble::rowid_to_column("index")
        }
      )
    }
  )

  result <- result |>
    purrr::list_rbind(names_to = "blast_db_path")

  if (nrow(result) == 0) {
    abort_rcrux_mini_error(
      "Found no hits when searching for multi ID stuff TODO"
    )
  }

  # result <- .add_reference_accessions(result)

  checkmate::assert_names(
    names(result),
    must.include = c(
      "index",
      "blast_ordinal_id",
      "subject_accession_version",
      "sequence_hash_value",
      "blast_db_path",
      "reference_subject_accession_version",
      "unique_subject_taxonomy_ids"
    )
  )

  # result <- result |> .add_reference_accessions()

  # Now that we have the expanded accessions, we need to link some of this data back to the original BLAST hits.  For each row in the blastdbcmd result, we need to convert it to a blast hit row.
  #
  # Blast hit rows have these columns:
  #
  # c(
  #   "query_accession",
  #   "subject_accession_version",
  #   "percent_identical_matches",
  #   "alignment_length",
  #   "expect_value",
  #   "subject_sequence_length",
  #   "subject_alignment_start",
  #   "subject_alignment_end",
  #   "subject_aligned_sequence",
  #   "unique_subject_taxonomy_ids"
  # )
  #
  # So you need to join the blastdbcmd reference_accession maps to blast result subject_accession_version. These should ALL be present. If it isn't there there is a bug somewhere, or I have misunderstood the way the blastdbcmd handles version numbers. (but they should be present in the blastdbcmd output)
  #
  # The unique_subject_taxonomy_ids in the blast output must be replaced by the taxonomy_id from the dbcmd output
  # The subject_accession_version in the blast output must be replaced by the accession in the dbcmd output
  # All the rest of the dbcmd columsn should drop out

  dplyr::left_join(
    x = result,
    y = hits_with_multiple_taxonomy_ids,
    by = dplyr::join_by(
      x$cluster_reference_accession,
      y$subject_accession_version
    )
  )

  # TODO: do I need to fix the output rows?
  # TODO: checkmate assert the column names
}

# TODO: I don't really want this as a top level item.
# TODO: fix the name
.sequence_hash_column_specification <- readr::cols(
  # You will need this for grouping sequences later. It's an int but let's
  # treat it as a character.
  blast_ordinal_id = readr::col_character(),
  # This name matches the blast hit name
  subject_accession_version = readr::col_character(),
  sequence_hash_value = readr::col_character(),
  # This name matches the blast hit name
  unique_subject_taxonomy_ids = readr::col_character()
)


# TODO: the only difference between this and the .pull_sequences_from_blast_db is the outfmt option
.pull_sequence_hashes_from_blast_db <- function(
  blastdbcmd,
  entry_batch_path,
  blast_db_path
) {
  # TODO: log the command that was run
  result <- processx::run(
    command = blastdbcmd,
    args = c(
      "-db",
      blast_db_path,
      "-dbtype",
      "nucl",
      "-entry_batch",
      entry_batch_path,
      "-outfmt",
      # This is blast ordinal ID (OID), accession, sequence hash value, taxid
      "%o ~~~ %a ~~~ %h ~~~ %T"
    ),
    # TODO: log the command
    echo_cmd = FALSE,
    error_on_status = FALSE
  )

  if (rlang::is_na(result$status) || result$status != 0) {
    # blastdbcmd is a bit different with respect to the exit status. Say you are
    # searching for 1000 accessions. If even *one* of these accessions is not
    # found in the DB, then the exit code will be non-zero. The stdout will
    # still be populated with data that was found, the stderr will have the
    # errors. You need to check if there were actually critical errors or not.
    # If something was "skipped" or not found it will be something like
    # `Error: [blastdbcmd] Skipped ACCESSION`.
    no_true_errors <- stringr::str_split_1(result$stderr, pattern = "\r?\n") |>
      purrr::map(.f = function(stderr_line) {
        # If the line is empty or it matches the "skipped" then it's okay
        stderr_line == "" ||
          stringr::str_detect(
            string = stderr_line,
            pattern = stringr::regex(
              "error.*blastdbcmd.*skipped",
              ignore_case = TRUE
            )
          )
      }) |>
      all()

    if (!no_true_errors) {
      # Then at least some of the errors were ones that we can't recover from, so fail.
      abort_rcrux_mini_error(
        "There was some error when running blastdbcmd",
        processx_object = result
      )
    }
  }

  # TODO: I think if we're here the stdout should be okay now, but would be nice
  #       to do some sanity check.
  result$stdout |>
    I() |>
    # TODO: handle parsing problems
    readr::read_delim(
      col_names = names(.sequence_hash_column_specification$cols),
      col_types = .sequence_hash_column_specification,
      delim = " ~~~ "
    ) |>
    tibble::rowid_to_column("index")
}


.add_reference_subject_accession_versions <- function(.data_frame) {
  checkmate::assert_data_frame(.data_frame, min.rows = 1)
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
