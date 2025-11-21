# TODO: would be nice to include the blast_db_path in the output here.
#
#' Wrapper for \code{SnailBLAST::crawl} to run \code{blastn} that specifies
#' standard rCRUX error handlers.
#'
run_blastn <- function(
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
      # TODO: logging and also I think this gets truncated? (It does, regular message doesn't)
      warning(paste(
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
      warning(paste(
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
