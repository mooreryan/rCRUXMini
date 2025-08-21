run_amplicon_blast <- function(
  blastn_executable,
  query_paths,
  db_paths,
  jobs = 1
) {
  checkmate::assert_character(query_paths, min.len = 1)
  lapply(query_paths, function(query_path) {
    checkmate::assert_file_exists(query_path)
  })

  checkmate::assert_character(db_paths, min.len = 1)
  checkmate::assert_int(jobs, lower = 1)
  checkmate::assert_character(db_paths)

  query_target_pairs <- expand.grid(QueryPath = query_paths, DbPath = db_paths)

  result <- parallel::mcmapply(
    query_target_pairs$QueryPath,
    query_target_pairs$DbPath,
    FUN = function(query_path, db_path) {
      # TODO
    },
    mc.preschedule = FALSE,
    mc.cores = jobs,
  )

  # TODO
}


# TODO: handle these.Zero rows error, actual blast command errors...

run_blast <- function(
  blast_executable,
  query_paths,
  db_paths,
  empty_blast_result_callback,
  blast_result_column_specification,
  job_failed_callback = function(query_path, db_path, exit_status) {},
  parse_failed_callback = function(query_path, db_path, error_condition) {},
  blast_executable_directory = NULL,
  # A list with extra args to pass to blast,
  # e.g., c("-blah", "whatever", "-yo", "234")
  extra_blast_arguments = NULL,
  jobs = 1
) {
  sys_which(blast_executable, blast_executable_directory)

  checkmate::assert_character(query_paths, min.len = 1)
  checkmate::assert_character(db_paths, min.len = 1)

  checkmate::assert_int(jobs, lower = 1)

  checkmate::assert_character(extra_blast_arguments, null.ok = TRUE)

  # "unique" checks that names are unique and that none are "" or NA
  checkmate::assert_list(blast_result_column_specification, names = "unique")

  checkmate::assert_function(
    empty_blast_result_callback,
    nargs = 0
  )
  checkmate::assert_function(
    job_failed_callback,
    args = c("query_path", "db_path", "exit_status")
  )

  checkmate::assert_function(
    parse_failed_callback,
    args = c("query_path", "db_path", "error_condition")
  )

  # This works fine even when each argument is length 1
  query_db_pairs <- expand.grid(query_path = query_paths, db_path = db_paths) |>
    checkmate::assert_data_frame(
      min.rows = 1,
      col.names = c("query_path", "db_path"),
      types = c("character", "character")
    )

  blast_results <- parallel::mcmapply(
    query_db_pairs$query_path,
    query_db_pairs$db_path,
    FUN = function(query_path, db_path) {
      process_result <- processx::run(
        command = blast_executable,
        # NOTE: append handles `extra_blast_arguments` being NULL
        args = append(
          c(
            "-db",
            db_path,
            "-query",
            query_path
          ),
          extra_blast_arguments
        ),
        error_on_status = FALSE,
        # These two args make mapply work like lapply so we get back a list.
        SIMPLIFY = FALSE,
        USE.NAMES = FALSE
      )

      if (rlang::is_na(process_result$status) || process_result$status != 0) {
        # The job failed. Either need to retry or do something else.
        job_failed_callback(
          query_path = query_path,
          db_path = db_path,
          exit_status = process_result$status
        )
        empty_blast_result_callback()
      } else {
        # Try to read the BLAST data
        rlang::try_fetch(
          expr = {
            blast_hit_data <- process_result$stdout |>
              readr::read_tsv(
                col_names = names(blast_result_column_specification),
                col_types = blast_result_column_specification
              )
            readr::stop_for_problems(blast_hit_data)
            # NOTE: since we specify col_names and col_types, if there are no
            # rows, this will be the good default value.
            blast_hit_data
          },
          # If there are any problems reading the data, report the error and
          # return
          # the default value.
          error = function(error_condition) {
            parse_failed_callback(
              query_path = query_path,
              db_path = db_path,
              error_condition = error_condition
            )
            empty_blast_result_callback()
          }
        )
      }
    },
    mc.preschedule = FALSE,
    mc.cores = jobs,
  )

  checkmate::assert_list(blast_results, min.len = 1)
  sapply(blast_results, function(blast_result) {
    checkmate::check_data_frame(blast_result)
  })

  dplyr::bind_rows(blast_results)
}
