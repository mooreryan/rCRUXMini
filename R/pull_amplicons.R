# NOTE: It is *possible* that you have duplicate accession-range pairs, if the
# user provided BLAST DBs that have some overlapping sequences. This means that
# you could in theory have exact duplicates in the resulting amplicon
# files/data. (TODO may need to address this -- not sure if it's really a big
# deal, because when you blast those amplicons, you would get "duplicate" hits
# too...so potentially as long as we dedup the amplicon blast results, we can
# ignore this). The question is, do I care to attach the range to the resulting
# IDs? I think we will just leave them off for now.

# START HERE: need to split the amplicons into multiple files as well so that way you are blasting multiple amplicon files each against all of the DBs

# blastdbcmd <- sys_which_or_abort(
#   command = "blastdbcmd",
#   ncbi_bin_directory = ncbi_bin_directory
# )

# This returns a tibble with id, sequence, and blast_db_path from which the
# amplicon sequence originated.
#
# TODO: actually, I can reuse this function to do the expand multi tax ids step as well
pull_amplicons <- function(
  blastdbcmd,
  amplicon_coordinates,
  chunks,
  blast_db_paths,
  ncbi_bin_directory = NULL
) {
  print("amplicon_coordinates")
  print(amplicon_coordinates)
  entry_batch_data <- .prepare_blastdbcmd_entry_batch_data(amplicon_coordinates)

  print("entry_batch_data")
  print(entry_batch_data)

  entry_batch_paths <- .write_pieces_to_tempfiles(
    entry_batch_data,
    chunks = chunks,
    delim = " "
  )

  # These need to be unique as we will be using them as names.
  checkmate::assert_vector(blast_db_paths, min.len = 1, unique = TRUE)
  # These names will be used to let the user know from which DB an amplicon
  # originated
  names(blast_db_paths) <- blast_db_paths

  # TODO: I actually don't think I need to pull this part into chunks. We would rather run the blastdbcmd

  result <- furrr::future_map2(
    .x = entry_batch_paths,
    .y = blast_db_paths,
    .f = function(entry_batch_path, blast_db_path) {
      rlang::try_fetch(
        .pull_sequences_from_blast_db(
          blastdbcmd = blastdbcmd,
          entry_batch_path = entry_batch_path,
          blast_db_path = blast_db_path
        ),
        # TODO: need to check for the rcrux_error instead
        error = function(condition) {
          rlang::warn("something bad happend in .pull_sequences_from_blast_db")
          print(condition)
          # TODO: log this error

          # TODO: this type is connected to the read_delim call in
          #      .pull_sequences_from_blast_db, so it would be nice to connect
          #      them
          tibble::tibble(
            index = integer(0),
            blast_ordinal_id = integer(0),
            accession = character(0),
            sequence_hash_value = character(0),
            sequence = character(0)
          )
        }
      )
    }
  )

  result <- result |>
    purrr::list_rbind(names_to = "blast_db_path")

  # TODO: these amplicons have the forward primer, but not the reverse. It needs
  # both or neither.

  checkmate::assert_names(
    names(result),
    permutation.of = c(
      "index",
      "blast_ordinal_id",
      # TODO: need to be more specific with this accession name (version? subject?)
      "accession",
      "sequence_hash_value",
      "sequence",
      "blast_db_path"
    )
  )

  result
}

# blastdbcmd_entry_batch_path <- file.path(output_directory, "entry_batch.txt")
# checkmate::assert_path_for_output(blastdbcmd_entry_batch_path)

# NOTE: It is *possible* that you have duplicate accession-range pairs, if the
# user provided BLAST DBs that have some overlapping sequences. This means that
# you could in theory have exact duplicates in the resulting amplicon
# files/data. (TODO may need to address this -- not sure if it's really a big
# deal, because when you blast those amplicons, you would get "duplicate" hits
# too...so potentially as long as we dedup the amplicon blast results, we can
# ignore this). The question is, do I care to attach the range to the resulting
# IDs? I think we will just leave them off for now. (this was the note in the old function)
.prepare_blastdbcmd_entry_batch_data <- function(amplicon_coordinates) {
  # Empty data frames shouldn't be passed in here, they should be handled before
  # you get to this function
  checkmate::assert_data_frame(
    amplicon_coordinates,
    min.rows = 1
  )
  checkmate::assert_names(
    names(amplicon_coordinates),
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

  # TODO: this shouldn't always be the reverse stop
  amplicon_coordinates |>
    # # Take only the columns that we need to construct the entry_batch file
    # dplyr::select("accession", "forward_start", "reverse_stop") |>
    dplyr::mutate(
      far_left_coordinate = pmin(
        .data$forward_start,
        .data$forward_stop,
        .data$reverse_start,
        .data$reverse_stop,
        na.rm = TRUE
      ),
      far_right_coordinate = pmax(
        .data$forward_start,
        .data$forward_stop,
        .data$reverse_start,
        .data$reverse_stop,
        na.rm = TRUE
      )
    ) |>
    dplyr::select("accession", "far_left_coordinate", "far_right_coordinate") |>

    # Original code uses forward_stop and reverse_stop, which I think is a bug. -- actually not sure if its a bug since they will change the reverse start/stop depending on orientation...but the primers should probably be included in the product so theirs is probably a bug regardless
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
    dplyr::select("accession", "range")
}

.write_to_tempfile <- function(data_frame, delim) {
  checkmate::assert_data_frame(data_frame)

  filename <- tempfile()

  data_frame |> readr::write_delim(file = filename, delim = delim)

  filename
}

.write_to_tempfasta <- function(.df, id_column, sequence_column) {
  checkmate::assert_data_frame(.df)
  checkmate::assert_names(
    names(.df),
    must.include = c(id_column, sequence_column)
  )

  filename <- tempfile()

  .df |>
    write_fasta(
      to_filename = filename,
      id_column = id_column,
      sequence_column = sequence_column
    )

  filename
}

.write_pieces_to_tempfiles <- function(.data_frame, chunks, delim) {
  print(".data_frame")
  print(.data_frame)
  print(dim(.data_frame))
  checkmate::assert_data_frame(.data_frame, min.rows = 1)
  checkmate::assert_count(chunks, positive = TRUE)
  checkmate::assert_string(delim, min.chars = 1)

  purrr::map_chr(
    .x = .chunk_data_frame(.data_frame, chunks),
    .f = function(df_chunk) {
      # If any of these fail, then the whole map will fail
      .write_to_tempfile(data_frame = df_chunk, delim = delim)
    }
  )
}

.write_pieces_to_tempfastas <- function(
  .df,
  chunks,
  id_column,
  sequence_column
) {
  checkmate::assert_data_frame(.df, min.rows = 1)
  checkmate::assert_count(chunks, positive = TRUE)
  checkmate::assert_data_frame(.df)
  checkmate::assert_names(
    names(.df),
    must.include = c(id_column, sequence_column)
  )

  purrr::map_chr(
    .x = .chunk_data_frame(.df, chunks),
    .f = function(df_chunk) {
      # If any of these fail, then the whole map will fail
      .write_to_tempfasta(
        .df = df_chunk,
        id_column = id_column,
        sequence_column = sequence_column
      )
    }
  )
}

.chunk_data_frame <- function(data_frame, chunks) {
  checkmate::assert_data_frame(data_frame, min.rows = 1)
  checkmate::assert_count(chunks, positive = TRUE)

  n <- nrow(data_frame)

  if (n < chunks) {
    chunks <- n
  }

  data_frame |>
    dplyr::mutate(
      chunk_id = ceiling(seq_len(n) * chunks / n)
    ) |>
    # TODO: dplyr docs say tihs may be removed in the future
    dplyr::group_split(.data$chunk_id) |>
    purrr::map(function(data_frame) {
      data_frame |> dplyr::select(-"chunk_id")
    })
}


# TODO: when running the blastdbcmd, you need to compare the accessions that
# actually made it into the stdout (which means they were found in the DB),
# compared to all that you tried. Because the original tracks the accessions
# that weren't found in the DB. (Note that they should all be found since we are
# pulling from the same DBs. The trick is that since we have multiple DBs you
# need to somehow aggregate this data across all DBs because we wouldn't expect
# every accession to be found in every DB, we only want to track the accessions
# that are found in none of the DBs)

# TODO: I need to make this function more general since I need to run blastdbcmd for other stuff too. (Could make an include sequences option or something)
.pull_sequences_from_blast_db <- function(
  blastdbcmd,
  entry_batch_path,
  blast_db_path
) {
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
      # This is blast ordinal ID (OID), accession, sequence hash value, sequence
      #
      # TODO: be careful about blast ordinal IDs. You could have a case where a
      # mostly identical hit has a different blast ordinal ID if it came from a
      # different DB. The e-values would likely not match though as the DB
      # wouldn't likely have the same exact stuff in it.
      "%o ~~~ %a ~~~ %h ~~~ %s"
    ),
    # TODO remove this
    echo_cmd = TRUE,
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
      print() |>
      purrr::map_lgl(.f = function(stderr_line) {
        print("stderr_line")
        print(stderr_line)
        # If the line is empty or it matches the "skipped" then it's okay
        a <- stderr_line == ""
        b <- stringr::str_detect(
          string = stderr_line,
          pattern = stringr::regex(
            "error.*blastdbcmd.*skipped",
            ignore_case = TRUE
          )
        )
        # TODO: inline these again
        a || b
      }) |>
      print() |>
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
      col_names = c(
        # You will need this for grouping sequences later
        "blast_ordinal_id",
        "accession",
        "sequence_hash_value",
        "sequence"
      ),
      col_types = "iccc",
      delim = " ~~~ "
    ) |>
    tibble::rowid_to_column("index")
}
