parse_primer_blast2 <- function(
  parse_primer_blast_command,
  primer_blast_tsvs,
  output_tsv,
  maximum_mismatches,
  minimum_length,
  maximum_length,
  num_threads
) {
  checkmate::assert_character(primer_blast_tsvs, min.len = 1)
  checkmate::assert_string(output_tsv)
  checkmate::assert_int(maximum_mismatches)
  checkmate::assert_int(minimum_length)
  checkmate::assert_int(maximum_length)
  checkmate::assert_int(num_threads)

  command_result <- processx::run(
    command = parse_primer_blast_command,
    args = c(
      "--out-file",
      output_tsv,
      "--threads",
      toString(num_threads),
      "--maximum-mismatches",
      toString(maximum_mismatches),
      "--minimum-length",
      toString(minimum_length),
      "--maximum-length",
      toString(maximum_length),
      primer_blast_tsvs
    ),
    echo_cmd = FALSE,
    error_on_status = FALSE,
    # Putting "current" inherits vars from parent env too
    env = c("current", RUST_BACKTRACE = "full")
  )

  if (is.na(command_result$status) || command_result$status != 0) {
    abort_rcrux_mini_error(
      "There was some error when running the parse primer blast script.",
      processx_object = command_result
    )
  }
}
