describe("the command line interface", {
  # This is a basic smoke test for the CLI interface for rCRUXMini
  it("works", {
    config_data <- create_config_for_test()
    on.exit(config_data$on_exit(), add = TRUE)

    script_path <- testthat::test_path("../../inst/cli/rCRUXMini.R")
    repo_path <- testthat::test_path("../..")

    run_result <- processx::run(
      "Rscript",
      args = c("--vanilla", script_path, config_data$config_file, repo_path),
      error_on_status = FALSE,
      echo_cmd = TRUE
    )
    print("status")
    cat(run_result$status)
    print("stdout")
    cat(run_result$stdout)
    print("stderr")
    cat(run_result$stderr)

    output_files <- list.files(
      path = config_data$output_directory,
      recursive = TRUE,
      full.names = FALSE
    ) |>
      purrr::discard(function(file) {
        stringr::str_ends(file, pattern = ".rds")
      }) |>
      purrr::map_chr(redact_tmp_filename)

    output_files |> expect_snapshot()
  })
})
