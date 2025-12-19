# # This is blast ordinal ID (OID), accession, sequence hash value, taxid
# "%o ~~~ %a ~~~ %h ~~~ %T"
#

#' Run blastdbcmd
#'
#' @param blastdbcmd todo
#' @param entry_batch_path todo
#' @param blast_db_path todo
#' @param outfmt todo
#'
#' @returns An rcrux "result", either ok like \code{list(blast_db_path = blast_db_path, stdout = result$stdout)}, or an
#'   error, with a list containing the keys message, processx_object, condition.
#'   The last two may be NULL.
#'
run_blastdbcmd <- function(
  blastdbcmd,
  entry_batch_path,
  blast_db_path,
  outfmt
) {
  rlang::try_fetch(
    {
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
          outfmt
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
        #

        is_blank <- function(s) {
          trimws(s) == ""
        }

        is_skipped <- function(s) {
          stringr::str_detect(
            string = s,
            pattern = stringr::regex(
              "error.*blastdbcmd.*skipped",
              ignore_case = TRUE
            )
          )
        }

        error_message_is_okay <- function(stderr_line) {
          is_blank(stderr_line) || is_skipped(stderr_line)
        }

        stderr_lines <- stringr::str_split_1(result$stderr, pattern = "\r?\n")

        all_error_messaegs_are_okay <- purrr::map_lgl(
          .x = stderr_lines,
          .f = error_message_is_okay
        ) |>
          all()

        if (!all_error_messaegs_are_okay) {
          # Then at least some of the errors were ones that we can't recover
          # from, so fail.
          #
          # TODO: document the return type
          return(result__error(list(
            message = "There was some error when running blastdbcmd",
            processx_object = result,
            condition = NULL
          )))
        }
      }

      # If we are all good, return the stdout. It should be parsed by the
      # caller.
      result__ok(list(blast_db_path = blast_db_path, stdout = result$stdout))
    },
    # TODO: I'm actually not 100% sure why this try_fetch is needed. If any of
    # the things that are left in the try part give an error, then it is a
    # bug/programmer error, not blastdbcmd failing....unless something truly
    # weird happens in processx potentially.
    error = function(condition) {
      # TODO: we could recover by providing a default value here...
      print("something bad happend in run_blastdbcmd")
      for (name in names(condition)) {
        print(paste0("condition$", name, ": ", condition[[name]]))
      }

      result__error(list(
        message = "There was some very strange and unexpected error when running blastdbcmd",
        processx_object = NULL,
        condition = condition
      ))
    }
  )
}
